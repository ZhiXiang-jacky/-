# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 15:41:05 2022

@author: ZhiXiang
"""

# -*- coding: utf-8 -*-
"""
0606NSGAII+TOPSIS傳統排序 (機台0.1、總延誤0.6、總換模時間0.3) 挑最小總延誤(交配下降0.02)
將這個程式檔重整成
已經將資料處理成
1.NSGAII +TOPSIS 權重設定為W =[0.1,0.6,0.3]  OK
1-1. 傳統GA分兩段 Part1 和 Part 2 
2.把Total Completion Time 重新調整  半OK(先暫用)
3.修改連續次數停止條件 0601
"""

'''
1874筆新資料
關鍵變數名稱如下
'''


# 1.在計算每筆延誤天數：Delayed_days_df
# 2.甘特圖每個任務幾點開始幾點結束，所用的機台及物產品料號：GanttChart
# 3.紀錄每代的第一條前緣線即透過TOPSIS計算出來的結果：TOPSIS_Good_table


def clear_all():
    """Clears all the variables from the workspace of the spyder application."""
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue

        del globals()[var]


if __name__ == "__main__":
    clear_all()

import re
import csv
import pandas as pd
import numpy as np
import math
import time
import random
import copy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# import Draw_GanttChart as Draw
import plotly.express as px
from collections import Counter
from plotly.offline import download_plotlyjs, init_notebook_mode, plot
from plotly.figure_factory import create_gantt
# from operator import itemgetter
import datetime
import plotly as py
from itertools import compress  # 創三維空間大量資料回傳List 回傳index很有用
import xlrd
# from numba import jit
import TOPSISs as tp;  # TOPSIS 演算法(已將套件改成適合此案例方式)
# import plotly.io as pio
import Data_Normalization as DN;  # 嘉謙paper正歸化目標資料(已將套件改成適合此案例方式)
import multiprocessing
import pdb

nb_threads = 5

Screening_of_site = "OP3"  # 將BOM表"工序說明"篩選OP3
Machine_model_1 = 'PM2VA001'  # 篩選機台銑床型號此OP3篩選出'PM2VA001'
Machine_model_2 = 'PM9VA001'  # 篩選機台銑床型號此OP3篩選出'PM9VA001'

Machine_model_Number = []  # 為了去算有幾種機型
Machine_model_Number.append(Machine_model_1)
Machine_model_Number.append(Machine_model_2)

Number_Of_Job = 1874  # 0606改成最新版1874筆資料
Number_Of_Machine = 381  # 平行機台雖然最多是381 但用限制式 把有些訂單原本可以丟到PM9VA001的機台都重新選機至其他台
Number_Of_Chromosome = 50
Number_Of_Obj = 3

Chrossover_Rate = 0.8
Mutation_Rate = 0.2
'''learning_rate   幫助交配範圍大一點 (建議0.2~0.25) 
  【主要看訂單數決定 Ex 訂單數10個 若組成基因長度 會是 10*2 + 2種機型數(開始、結束)】 共 24個基因 大概交配一半以上就算多
  若訂單為3716筆，learning_rate 用0.25
  若訂單為10筆，learning_rate 用0.1

'''
learning_rate = 0.25

set_Timer = 70 * 60  # 分鐘*秒數
num_iterations = 9999999999999  # 迭代數
Continuous_termination = 500  # 連續Obj不變幾代就要停止
'''為了計算目標值 TSTP = Total completion time * Penalty_value_weight_TCT + Sum_Tardiness * Penalty_value_weight + Setup_times * Penalty_value_weight_Setup'''
Penalty_value_weight_TCT = 1  # Total_completion_time 權重設為1
Penalty_value_weight = 1  # Tardiness  權重設為1000
Penalty_value_weight_Setup = 1  # Setup_times  權重設為100

# [#1.5] Plot [幾個機台畫成一張圖]
plot_Machinenumber = 20;

'''NSGAII參數'''
# [#1.4] NSGA_II 參數 Rules variables
weight1 = 1;  # Target_Machineusing_KPI_Weight
weight2 = 1;  # Target_DeliveryRate_KPI_Weight
weight3 = 1;  # Target_ChangeLine_KPI_Weight
# weight4 = 1; # Target_WIP_KPI_Weight

Weights = np.array([weight1, weight2, weight3]);
Weights = Weights / Weights.sum();  # 讓權重總和為 1

'''TOPSIS參數'''
# [#1.4] TPOSIS 參數 Rules variables
TOPSISweight1 = 0.1;  # Target_Machineusing_KPI_Weight
TOPSISweight2 = 0.6;  # Target_DeliveryRate_KPI_Weight
TOPSISweight3 = 0.3;  # Target_ChangeLine_KPI_Weight
# # weight4 = 1; # Target_WIP_KPI_Weight

TOPSIS_Weights = np.array([TOPSISweight1, TOPSISweight2, TOPSISweight3]);
TOPSIS_Weights = TOPSIS_Weights / TOPSIS_Weights.sum();  # 讓權重總和為 1

criterias = ["-", "-", "-"];  # 正反指標 (看指標是望大還是望小特性)

# ===================================學校實驗室==========================================

'''讀取BOM表'''
BOM = pd.read_csv("0501_未整理BOM.csv",encoding= 'utf-8',dtype ={'orderBomKey':str,'工序代碼':str,
                                                    '工序說明':str,
                                                    'seqNo':int,
                                                    '加工時間':np.float64,
                                                    '設備/資源':str,
                                                    })
'''將BOM表作工序篩選Screening_of_site = OP3'''
BOM_filt = BOM["工序說明"] == Screening_of_site
BOM_DF = BOM.loc[BOM_filt]
BOM_DF.reset_index(drop=True, inplace=True)  # 重新更新Index_0424

# '''讀取BOM表'''
# BOM = pd.read_csv("5_01僅平行機台BOM.csv",dtype ={'orderBomKey':str,'工序代碼':str,
#                                                     '工序說明':str,
#                                                     'seqNo':int,
#                                                     '加工時間':np.float64,
#                                                     '設備/資源':str,
#                                                     })
'''將BOM表作工序篩選Screening_of_site = OP3'''
# =============================================================================
# # BOM_filt = BOM["工序說明"] == Screening_of_site
# # BOM_DF = BOM.loc[BOM_filt]
# # BOM_DF.reset_index(drop=True, inplace=True) #重新更新Index_0424
# =============================================================================

'''新增機台數量表'''
# Machine_Quantity_Table = pd.read_csv("其他.csv",encoding= 'ANSI',dtype ={'設備/資源':str,'數量':int})


# '''讀取Order檔'''

# Order = pd.read_csv("Order.csv",dtype ={'訂單編號':str,
#                                                     '產品料號':str,
#                                                     '目標數量':int,
#                                                     '最早開始時間':str,
#                                                     '交期':str,
#                                                     '優先度':int,
#                                                     'orderBomKey':str
#                                                     })

# ==============================04/24(不用看了)===============================================
# Order = pd.read_csv("0424剔除資料Order.csv", encoding= 'ANSI',dtype ={'訂單編號':str,
#                                                     '產品料號':str,
#                                                     '目標數量':int,
#                                                     '最早開始時間':str,
#                                                     '交期':str,
#                                                     '優先度':int,
#                                                     'orderBomKey':s.0tr
#                                                     })


# ======================0519剔除逾期Order檔1874筆=======================================================

Order = pd.read_csv("0519剔除逾期Order檔.csv",encoding= 'ANSI',dtype ={'訂單編號':str,
                                                    '產品料號':str,
                                                    '目標數量':int,
                                                    '最早開始時間':str,
                                                    '交期':str,
                                                    '優先度':int,
                                                    'orderBomKey':str
                                                    })

# Order = pd.read_csv("January_Order.csv",dtype ={'訂單編號':str,
#                                                     '產品料號':str,
#                                                     '目標數量':int,
#                                                     '最早開始時間':str,
#                                                     '交期':str,
#                                                     '優先度':int,
#                                                     'orderBomKey':str
#                                                     })
# Order = pd.read_csv("March_Order.csv",encoding= 'ANSI',dtype ={'訂單編號':str,
#                                                     '產品料號':str,
#                                                     '目標數量':int,
#                                                     '最早開始時間':str,
#                                                     '交期':str,
#                                                     '優先度':int,
#                                                     'orderBomKey':str
#                                                     })


# np.random.seed(10)  #測試

# Order_RandomRow = np.random.choice(a=3716, size=Number_Of_Job, replace=False, p=None)

# Order = Order.loc[Order_RandomRow]
'''---------------------新增機型、數量、Job對應加工時間-----------------------------------'''
# Order = pd.merge(Order,BOM_DF.loc[:,["產品料號","設備/資源"]],how = 'inner',on = "產品料號")
# Order = pd.merge(Order,Machine_Quantity_Table.loc[:,["設備/資源","機台數量"]],how = 'inner',on = "設備/資源")
# Order = pd.merge(Order,BOM_DF.loc[:,["產品料號","加工時間"]],how = 'inner',on = "產品料號")

'''0424資料預處理'''


def cut_text(text):
    textArr = re.findall(r'\w{8}', text)
    return textArr


dict_canRunTool = {}
for i in range(len(BOM_DF)):
    '''使用成'''
    dict_canRunTool[BOM_DF.iloc[i]['orderBomKey']] = BOM_DF.iloc[i]["設備/資源"]

for index, BomKey in enumerate(dict_canRunTool.keys()):
    '''(new)產品料號對應機台型號有多少機台數'''
    '''使用成字典BomKey 對應 機台，把String 改成List，方便以後選機'''
    dict_canRunTool[BomKey] = cut_text(dict_canRunTool[BomKey])

dict_BomKey_Processing = {}
for i in range(len(BOM_DF)):
    '''使用成字典ProcessTime'''
    dict_BomKey_Processing[BOM_DF.iloc[i]['orderBomKey']] = np.round(BOM_DF.iloc[i]["加工時間"], 3)

# =============================================================================
# Order_pt = []
# for i in range(len(Order)):
#     '''處理Order對應的ProcessTime 成List'''
#     Order_pt.append(dict_BomKey_Processing[Order.at[i,"orderBomKey"]])
# Order_pt = pd.DataFrame(Order_pt)  #將List轉換成DataFrame
# =============================================================================

Order = pd.merge(Order, BOM_DF.loc[:, ["orderBomKey", "設備/資源"]], how='inner', on="orderBomKey")  # 新增對應的可用機台Order表
Order = pd.merge(Order, BOM_DF.loc[:, ["orderBomKey", "加工時間"]], how='inner', on="orderBomKey")  # 新增加工時間至Order表

Order["總加工時間"] = round(Order['目標數量'] * Order['加工時間'] * 60, 3)  # 轉成秒數

'''讀取資源群組對應表(OP3所有機台)'''
# DF_All_Machines_Name = pd.read_csv("資源群組對應表.csv",encoding= 'ANSI')

start = time.time()
'''讀取Mould(換線資訊表)資料'''
Mould =pd.read_csv("0504最新版換線資訊表v2.csv",dtype ={'設備/資源編號':str,
                                                    '當前產品料號':str,
                                                    '目標產品料號':str,
                                                    '換線時間':int})

'''篩選機台銑床型號 PM2VA001或 PM2VD001'''
# Mould_filt = (Mould["設備/資源編號"] == Machine_model_1) | (Mould["設備/資源編號"] == Machine_model_2)
# Mould_DF = Mould.loc[Mould_filt]

'''預處理OP3所使用機型與機台數'''
# Mc_model_Quantity = Machine_Quantity_Table[(Machine_Quantity_Table["設備/資源"] == Machine_model_1) |(Machine_Quantity_Table["設備/資源"] == Machine_model_2)]
# Mc_model_Quantity.reset_index(inplace=True, drop=True)  #Index 從0開始
# Mc_model_Quantity = dict(np.array(Mc_model_Quantity)) #Dict 形式

'''機台數限制5%'''
# =============================================================================
# Mumber_of_machines_limit = (Mc_model_Quantity[Machine_model_1] + Mc_model_Quantity[Machine_model_2]) * 0.9
# =============================================================================

# =============================================================================
# # =============================================================================
# # '''改成array快一點點'''
# # # Mould_DF = Mould.loc[Mould_filt].values   #改成array


'''字典第一代(方法一)'''
Mould["當前_目標"] = Mould.apply(lambda x: x['當前產品料號'] + x['目標產品料號'], axis=1)
Mould_dict = Mould[['當前_目標', '換線時間']]
Mould_dict.index = Mould['當前_目標']
Mould_dict = Mould_dict.drop('當前_目標', axis=1)
Mould_dict = Mould_dict["換線時間"].to_dict()

# =============================================================================
# '''(new)產品料號對應機台型號有多少機台數'''
# Order['產料_機型'] = Order.apply(lambda x: x['產品料號'] + x['設備/資源'] , axis = 1)
# Order_dict = Order[['產料_機型','機台數量']]
# Order_dict.index = Order['產料_機型']
# Order_dict = Order_dict.drop('產料_機型',axis = 1)
# Order_dict = Order_dict['機台數量'].to_dict()
# =============================================================================

end = time.time()
'''我的電腦估計花49秒左右，嘉千電腦大概23秒挖糙'''
print('模具時間', end - start)

Machine_Start= pd.read_csv("5_01平行機台開始時間.csv")
'''資料預處理-訂單名稱'''
OrderName = []
for i in range(len(Order)):
    str_Name = str(Order.iloc[i]['訂單編號'])
    OrderName.append(str_Name)
'''資料預處理-產品料號'''
Product_Part_Number = []
for i in range(len(Order)):
    str_Name = str(Order.iloc[i]['產品料號'])
    Product_Part_Number.append(str_Name)

'''===新增的訂單對應機台型號預處理(或許會用到，結果沒用到之後可刪)==='''
Order_Corres_Machine_Model = []
for i in range(len(Order)):
    str_Name = str(Order.iloc[i]['設備/資源'])
    Order_Corres_Machine_Model.append(str_Name)

'''資料預處理-ArrivalTime'''
Arrival_Time = []
for i in range(len(Order)):
    str_time = str(Order.iloc[i]['最早開始時間'])
    timeStamp = int(time.mktime(time.strptime(str_time, "%Y/%m/%d %H:%M")))
    Arrival_Time.append(timeStamp)
'''資料預處理-RecoverTime'''
RecoverTime = []
for i in range(len(Machine_Start)):
    str_time = str(Machine_Start.iloc[i]['機台開始時間'])
    timeStamp = int(time.mktime(time.strptime(str_time, "%Y/%m/%d %H:%M")))
    RecoverTime.append(timeStamp)

'''資料預處理-ProcessTime'''
Process_Time = []
for i in range(len(Order)):
    '''已轉換成秒數'''
    timeStamp = Order.iloc[i]['總加工時間']
    Process_Time.append(timeStamp)

'''資料預處理-交期_time'''
Due_Date_times = []
for i in range(len(Order)):
    str_time = str(Order.iloc[i]['交期'])
    timeStamp = int(time.mktime(time.strptime(str_time, "%Y/%m/%d %H:%M")))
    Due_Date_times.append(timeStamp)


def read_excel_data():
    '''甘特圖106種產品料號顏色整理'''
    filename = '顏色產品料號.xlsx'
    data = xlrd.open_workbook(filename)
    table = data.sheet_by_name('工作表1')
    row_num = table.nrows  # 行数
    # col_num = table.ncols  # 列数
    datas = dict([])  # 这步也要转字典类型
    for i in range(row_num):
        colorss = dict([table.row_values(i)])  # 这一步就要给它转字典类型，不然update没法使用
        datas.update(colorss)
    # print(datas)
    return datas


if __name__ == "__main__":
    colorss = read_excel_data()


# %%交配突變
def single_point_crossover(A, B, X):
    A_new = np.append(A[:X], B[X:])
    B_new = np.append(B[:X], A[X:])
    # print(A_new)
    # print(B_new)
    return A_new, B_new


def multi_point_crossover(A, B, X):
    for i in X:
        A, B = single_point_crossover(A, B, i)
    return A, B


def Crossover(C, D):
    '''交配，採兩點交配!!!'''
    Crossover1 = C
    Crossover2 = D
    '''例如 有10個Job ,最多可能交配點為  [1,9]，不會到10,否則就變成單點交配'''
    # 這裡到最後去思考會不會可用機台數基因最後一個結束時間不會變更!!
    X = np.random.choice(range(1, Number_Of_Job * 2), size=2, replace=False)
    '''下面說明當X任兩點若小於Number_Of_Job * (1+learning_rate) 範圍 就不做交配!!，learning_rate在上面可以做變更'''
    while (abs(X[0] - X[1]) < round((Number_Of_Job) * (learning_rate))):  # 0602修改不加1的learning_rate 因為基因數沒有到很長
        X = np.random.choice(range(1, Number_Of_Job * 2), size=2, replace=False)
    return Crossover1, Crossover2, X


def Mutation(index):
    '''採2點突變法'''
    '''突變採任兩點取出不放回做突變'''
    temp_position1, temp_position2 = np.random.choice(range(1, Number_Of_Job * 2), size=2, replace=False)
    temp_RandomNumber1, temp_RandomNumber2 = [random.random() for i in range(2)]
    # print("第幾格被換 %d ,換成哪個一個亂數%f"%(temp_position,temp_RandomNumber))

    Temp_Mutation[index][temp_position1] = temp_RandomNumber1
    Temp_Mutation[index][temp_position2] = temp_RandomNumber2

    # print('新的',Temp_Mutation[index])

    return Temp_Mutation


def Mc_Mutation(Mc_Gene):
    '''採2點突變法'''
    '''突變採任兩點取出不放回做突變'''
    # temp_position1,temp_position2 =np.random.choice(range(0,4),size = 2,replace = False)
    # temp_RandomNumber1,temp_RandomNumber2 = [random.random() for i in range(2)]
    # Mc_Gene[temp_position1] = temp_RandomNumber1
    # Mc_Gene[temp_position2] = temp_RandomNumber2

    temp_position1 = np.random.choice(range(0, 4), size=1, replace=False)
    temp_RandomNumber1 = [random.random() for i in range(1)]
    # print("第幾格被換 %d ,換成哪個一個亂數%f"%(temp_position,temp_RandomNumber))

    Mc_Gene[temp_position1] = temp_RandomNumber1

    # print('新的',Temp_Mutation[index])

    return Mc_Gene


def List_to_Np(List):
    '''將temp_Mc__Corres_To_Job_End改成矩陣形式'''
    n = len(max(List, key=len))
    # print(n)

    # Make the lists equal in length
    lst_2 = [x + [0] * (n - len(x)) for x in List]
    # print(lst_2)
    matrix = np.array(lst_2) / (60 * 60 * 24)  # 已換成天數單位
    return matrix


def ArrivalTime(Arrivaltime_stamp):
    '''將染色體Part2 改成使用ArrivalTime排序'''
    timeString = Arrivaltime_stamp  # 時間格式為字串
    struct_time = time.strptime(timeString, "%Y/%m/%d %H:%M")  # 轉成時間元組
    time_stamp = int(time.mktime(struct_time))  # 轉成時間戳
    return time_stamp


def Machine_Corres_To_Job_ArrivalTime_Pt_PPN(temp_Machine_Corresponds_To_Job_sorted):
    '''直接用temp_Machine_Corresponds_To_Job_sorted 的Job index 去抓取 預處理的ArrivalTime 、PT 各自時間'''
    '''處理成相對位置!!'''
    temp_ArrivalTime = []
    temp_Mc_Corres_To_Job_Pt = []
    temp_Mc_Corres_To_Job_Product_Part_Number = []

    # 總共要跑幾次機台
    for Num_Mc in range(len(temp_Machine_Corresponds_To_Job_sorted)):
        temp = []
        temp_TotalTime = []
        temp_Product_PN = []
        # 每個機台內跑幾個Job數
        for McToJobLengh in range(len(temp_Machine_Corresponds_To_Job_sorted[Num_Mc])):

            time_stamp = Arrival_Time[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh]] - RecoverTime[
                Num_Mc]  # 轉成時間戳(有減機台相對位置)

            temp_TotalTime.append(Process_Time[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][
                McToJobLengh]])  # temp_TotalTime增加每個Job總加工時間

            temp_Product_PN.append(Product_Part_Number[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][
                McToJobLengh]])  # temp_Product_PN增加每個Job的產品料號
            if time_stamp < 0:
                '''若time_stamp <= 0 ，轉換為0 ，代表貨到達此刻可以加工'''
                time_stamp = 0

            temp.append(time_stamp)

        temp_ArrivalTime.append(temp)
        temp_Mc_Corres_To_Job_Pt.append(temp_TotalTime)
        temp_Mc_Corres_To_Job_Product_Part_Number.append(temp_Product_PN)

    return temp_ArrivalTime, temp_Mc_Corres_To_Job_Pt, temp_Mc_Corres_To_Job_Product_Part_Number


def Setup_Machine_Corres_To_Job_PPN(temp_Machine_Corresponds_To_Job_sorted):
    '''0526  直接用Zip_OrderNumber_Stime_Endtime 的Job index 去抓取 預處理的各自Job_Product_Part_Number時間'''
    '''換模專用'''
    temp_Mc_Corres_To_Job_Product_Part_Number = []

    # 總共要跑幾次機台
    for Num_Mc in range(len(temp_Machine_Corresponds_To_Job_sorted)):
        temp_Product_PN = []
        # 每個機台內跑幾個Job數
        for McToJobLengh in range(len(temp_Machine_Corresponds_To_Job_sorted[Num_Mc])):
            temp_Product_PN.append(Product_Part_Number[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][
                McToJobLengh]])  # temp_Product_PN增加每個Job的產品料號

        temp_Mc_Corres_To_Job_Product_Part_Number.append(temp_Product_PN)

    return temp_Mc_Corres_To_Job_Product_Part_Number


def Thebest_Machine_Corres_To_Job_ArrivalTime_Pt_PPN(temp_Machine_Corresponds_To_Job_sorted):
    '''處理成相對位置!!'''
    temp_Mc_Corres_To_Job_Index = []
    temp_ArrivalTime = []
    temp_Mc_Corres_To_Job_Pt = []
    temp_Mc_Corres_To_Job_Product_Part_Number = []

    # 總共要跑幾次機台
    for Num_Mc in range(len(temp_Machine_Corresponds_To_Job_sorted) - 1):

        temp_Job_index = []
        temp_Arr = []
        temp_TotalTime = []
        temp_Product_PN = []
        # 每個機台內跑幾個Job數
        for McToJobLengh in range(len(temp_Machine_Corresponds_To_Job_sorted[Num_Mc])):
            # print("機台%s, 第 %s 個Job" % (Num_Mc,McToJobLengh))
            if temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][0] == 'Stop':
                temp_Arr.append(temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][1] - RecoverTime[Num_Mc])
            else:
                temp_Job_index.append(temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][0])  # 新增Job的Index

                temp_TotalTime.append(
                    Process_Time[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][0]])  # 新增總加工時間

                temp_Product_PN.append(
                    Product_Part_Number[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][0]])  # 新增產品料號

                time_stamp = Arrival_Time[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][0]] - \
                             RecoverTime[Num_Mc]  # 到達時間轉成時間戳(有減機台相對位置)

                if time_stamp < 0:
                    '''若time_stamp <= 0 ，轉換為0 ，代表貨到達此刻可以加工'''
                    time_stamp = 0

                temp_Arr.append(time_stamp)

        temp_Mc_Corres_To_Job_Index.append(temp_Job_index)
        # print(temp_Mc_Corres_To_Job_Index)
        temp_ArrivalTime.append(temp_Arr)
        # print(temp_ArrivalTime)
        temp_Mc_Corres_To_Job_Pt.append(temp_TotalTime)
        # print(temp_Mc_Corres_To_Job_Pt)
        temp_Mc_Corres_To_Job_Product_Part_Number.append(temp_Product_PN)
        # print(temp_Mc_Corres_To_Job_Product_Part_Number)

    return temp_Mc_Corres_To_Job_Index, temp_ArrivalTime, temp_Mc_Corres_To_Job_Pt, temp_Mc_Corres_To_Job_Product_Part_Number


# %%'''計算非空機台，意旨真的機台數量'''
def Culate_Num_Mc_Actual_Open(Chromosom_Loop):
    '''計算非空機台'''
    Num_Machine_Actual_Open = []
    for i in range(len(Chromosom_Loop)):
        temp_Machine_Empty = []
        temp_Number_Total_Empty = 0
        temp_Open_number_mc = len(Chromosom_Loop[i])
        for j in range(len(Chromosom_Loop[i])):
            if not Chromosom_Loop[i][j]:
                # print(f'{i}條染色體{j}List is empty')
                temp_Machine_Empty.append(j)
                temp_Number_Total_Empty += 1
        temp_Open_number_mc = temp_Open_number_mc - temp_Number_Total_Empty
        Num_Machine_Actual_Open.append(temp_Open_number_mc)
    return Num_Machine_Actual_Open


# %%將訂單編號、開始結束綁再一起 [(J1,St1,Et1),(J2,St2,Et2)]
def OrderNumber_Stime_Endtime_Zip(Name, St, End, eachMaxmakespan):
    '''將訂單編號、開始結束綁再一起 [(J1,St1,Et1),(J2,St2,Et2)]'''
    temp_OrderNumber_Stime_Endtime = []
    for Num_Mc in range(len(End)):
        temp = list(zip(Name[Num_Mc], St[Num_Mc], End[Num_Mc]))
        temp_OrderNumber_Stime_Endtime.append(temp)
    temp_OrderNumber_Stime_Endtime.append(eachMaxmakespan)
    return temp_OrderNumber_Stime_Endtime


# %%將(訂單編號、St、Endt)改成[[訂單編號、St、Endt]]
def Zip_OrderNumber_Stime_Endtime_converList(Zip_converList, eachMaxmakespan=0):
    '''將(訂單編號、St、Endt)改成[[訂單編號、St、Endt]]'''
    temp_Zip_to_List = []
    for Num_Mc in range(Number_Of_Machine):
        # print(Num_Mc)
        temp = []
        for Num_Job in Zip_converList[Num_Mc]:
            temp.append(list(Num_Job))
        temp_Zip_to_List.append(temp)
    temp_Zip_to_List.append(eachMaxmakespan)

    return temp_Zip_to_List


# %%將訂單編號、開始結束綁再一起 [(J1,ArrivalT1,Pt1),(J2,ArrivalT2,Pt2)]

def Thebest_OrderNumber_St_Pt_Zip(Name, ArrivalT, Pt):
    '''此Function跟上面OrderNumber_St_Pt_Zip差異是 這個函數使用在最好的那條解碼身上!!'''
    temp_OrderNumber_Stime_Endtime = []
    for Num_Mc in range(len(Pt)):
        temp = list(zip(Name[Num_Mc], ArrivalT[Num_Mc], Pt[Num_Mc]))
        temp_OrderNumber_Stime_Endtime.append(temp)

    return temp_OrderNumber_Stime_Endtime


# %%將(訂單編號、ArrivalT、Endt)改成[[訂單編號、ArrivalT、Endt]]

def Thebest_Zip_OrderNumber_Stime_Endtime_converList(Zip_converList):
    temp_Zip_to_List = []
    for Num_Mc in range(Number_Of_Machine):
        # print(Num_Mc)
        temp = []
        for Num_Job in Zip_converList[Num_Mc]:
            temp.append(list(Num_Job))
        temp_Zip_to_List.append(temp)

    # print(temp_Zip_to_List)
    return temp_Zip_to_List


def Thebest_List_to_Np(List):
    '''將temp_Mc__Corres_To_Job_End改成矩陣形式，為了方便後續Total_completion_time，將每個Job的結束時間做加總'''
    temp = []
    for i in range(len(List) - 1):
        try:
            a = pd.DataFrame(List[i])
            temp_sum = a[2].sum()
            temp.append(temp_sum)
        # print(a[2].sum())
        except:
            temp.append(0)
        temp_Total_completion_time = np.sum(temp)

    return temp_Total_completion_time / (60 * 60 * 24)  # 已換成天數單位


def get_between_days(start_sec, end_sec):
    '''獲得兩個日期之間的天數(前插法時會考慮)'''
    work_days = round(((end_sec - start_sec) / (24 * 60 * 60)), 3)  # 已換成天數單位
    # work_hours = int((end_sec - start_sec)/(60*60))
    # work_minutes = int((end_sec - start_sec)/(60))
    return work_days


def convergence_Fitness(nice_makespan, k):
    '''Elitist_Chromosomes_max_makespan的nice_Objective 呼叫近來當y座標,X座標是Iteration'''

    ypt = nice_makespan

    xpt = list(range(k))

    plt.title("TSTP")
    plt.xlabel("Iteration")
    plt.ylabel("Objective")
    # plt.scatter(xpt, ypt,s = 50, color = 'y')
    plt.plot(xpt, ypt, linestyle='solid', color='y')

    plt.show()
    return


def convergence_Fitness_subplot(Nb_Mc, Tardiness, TCT, k):
    # Import necessary libraries
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams['font.sans-serif'] = ['Microsoft JhengHei']
    plt.rcParams['axes.unicode_minus'] = False

    # Change the figure size
    plt.figure(figsize=[9, 7])

    # Preparing the data to subplots
    x = list(range(k))
    y1 = Nb_Mc.values
    y2 = Tardiness.values
    y3 = TCT.values
    # y1 = df_excel["機台數"].values
    # y2 = df_excel["總延誤時間"].values
    # y3 = df_excel["總完工時間"].values

    # Plot the subplots
    # Plot 1
    '''機台數'''
    plt.subplot(2, 2, 1).title.set_text("機台數")
    # plt.xlabel("品牌名稱")
    # plt.ylabel("文章數量")

    plt.subplot(2, 2, 1)
    plt.plot(x, y1, 'g', linewidth=2)

    # Plot 2
    '''總延誤時間'''
    plt.subplot(2, 2, 2).title.set_text("總延誤時間")
    # plt.xlabel("品牌名稱")
    # plt.ylabel("文章數量")
    plt.subplot(2, 2, 2)
    # plt.scatter(x, y2, color='k', linewidth=2)
    # plt.scatter(x, y2, color='k', linewidth=2)
    plt.plot(x, y2, 'r', linewidth=2)
    # Plot 3
    '''總完工時間'''
    plt.subplot(2, 2, 4).title.set_text("總換模時間")
    # plt.xlabel("品牌名稱")
    # plt.ylabel("文章數量")
    plt.subplot(2, 2, 4)
    plt.plot(x, y3, '-.y', linewidth=3)

    plt.show()
    return


def Number_Of_Consecutive_Iterations_Occurs(a):
    b = a.iloc[-1]  # b的值为要求的数字
    t = 0
    w = 1
    for k, v in enumerate(a):
        if k > 0:
            if v == b and a[k - 1] == b:
                t += 1
                if w < t:
                    w = t
            else:
                t = 1

    print("%s連續出現%s次" % (np.round(b, 4), t))
    return t


# =============================================================================
# def Order_Delay_Cost():
#     Over_Due_List = []
#     count = 0
#     for Num_Mc in range(len(Zip_OrderNumber_Stime_Endtime)-1):
#         temp = []
#
#         for McToJobLengh in range(len(Zip_OrderNumber_Stime_Endtime[Num_Mc])):
#             '''把交期設為相對時間'''
#             str_time = Due_Date_times[Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][0]]
#             Due_timeStamp = str_time - RecoverTime[Num_Mc]
#             if Due_timeStamp < 0 :
#                 Due_timeStamp = 0
#             if Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][2] > Due_timeStamp:
#                 count +=1
#
#             # print(Num_Mc,timeStamp)
#             temp.append(Due_timeStamp)
#         Over_Due_List.append(temp)
#     return
# =============================================================================
# %%
'''選擇機制'''
from pymoo.util.dominator import Dominator
import numpy as np
import pandas as pd
import copy


# 2. use crowd distance to sort the solution,so each soluton have two new definition
# caculate crowd distance
# filter_out_duplicates 過濾掉重複項

def fast_non_dominated_sort(F, **kwargs):
    """
    Parameters
    ----------
    F: numpy.ndarray
        objective values for each individual.
    strategy: str
        search strategy, can be "sequential" or "binary".
    Returns
    -------
        fronts: list
            Indices of the individuals in each front.
    References
    ----------
    X. Zhang, Y. Tian, R. Cheng, and Y. Jin,
    An efficient approach to nondominated sorting for evolutionary multiobjective optimization,
    IEEE Transactions on Evolutionary Computation, 2015, 19(2): 201-213.
    copy in pymoo packages by original NSGA-II author
    """

    # M = Dominator.calc_domination_matrix(F)
    M = Dominator.calc_domination_matrix(F.values)  # 要再查一下再做什麼

    # calculate the dominance matrix
    n = M.shape[0]

    fronts = []

    if n == 0:
        return fronts

    # final rank that will be returned
    n_ranked = 0
    ranked = np.zeros(n, dtype=int)  # 每一條染色體給ranked

    # for each individual a list of all individuals that are dominated by this one  #每個個體可以凌掠的解集合
    is_dominating = [[] for _ in range(n)]  # 每一個p 的  Sp空集合

    # storage for the number of solutions dominated this one
    n_dominated = np.zeros(n)  # 每個個體被別人凌掠的個數  就是np

    current_front = []  # 目前的前緣線的解 (裡面會放個體)

    for i in range(n):
        '''i是不同的個體'''
        for j in range(i + 1, n):
            '''倆倆去比較  如果rel  ==1 就要增加進去凌掠的解集合，這裡之後花時間研究'''
            rel = M[i, j]
            if rel == 1:
                is_dominating[i].append(j)
                n_dominated[j] += 1
            elif rel == -1:
                is_dominating[j].append(i)
                n_dominated[i] += 1

        if n_dominated[i] == 0:  # 若都沒有被別人凌掠就是代表第一名
            current_front.append(i)
            ranked[i] = 1.0  # 單純標記
            n_ranked += 1  # 有幾次front  = 0 的解

    # append the first front to the current front
    fronts.append(current_front)

    # while not all solutions are assigned to a pareto front
    while n_ranked < n:

        next_front = []  # 新的front集合用空集合

        # for each individual in the current front
        for i in current_front:

            # all solutions that are dominated by this individuals
            for j in is_dominating[i]:
                n_dominated[j] -= 1
                if n_dominated[j] == 0:
                    next_front.append(j)
                    ranked[j] = 1.0  # 單純標記
                    n_ranked += 1

        fronts.append(next_front)
        current_front = next_front

    return fronts


##2. use crowd distance to sort the solution,so each soluton have two new definition
## caculate crowd distance
## filter_out_duplicates 過濾掉重複項
def cdist(A, B, **kwargs):
    import scipy
    return scipy.spatial.distance.cdist(A.astype(float), B.astype(float), **kwargs)  # 計算兩個輸入集合中每對之間的距離


def find_duplicates(X, epsilon=1e-16):  # (duplicates = 重複)
    # calculate the distance matrix from each point to another
    D = cdist(X, X)

    # set the diagonal to infinity  #np.triu_indices(len(X)) 可用成上三角形(包含自己跟自己 例如 A-> A 距離 0 )
    D[np.triu_indices(len(X))] = np.inf

    # set as duplicate if a point is really close to this one
    is_duplicate = np.any(D <= epsilon, axis=1)

    return is_duplicate


def calc_crowding_distance(F, Weights, filter_out_duplicates=True):
    '''跟網頁的一樣'''
    n_points, n_obj = F.shape  # (表示目標值的形狀 例如  幾條染色體(每個點代表一個個體) * 目標式數量)

    if n_points <= 2:
        '''np.inf 意指 +∞，是沒有確切的數值的,類型為浮點型，只有兩個點，這兩個點都會是+∞ 因為要保留下來'''
        return np.full(n_points, np.inf)

    else:

        if filter_out_duplicates:
            # filter out solutions which are duplicates - duplicates get a zero finally
            is_unique = np.where(np.logical_not(find_duplicates(F, epsilon=1e-32)))[0]
        else:
            # set every point to be unique without checking it
            is_unique = np.arange(n_points)

        # index the unique points of the array
        # _F = F[is_unique]
        _F = F.values[is_unique].astype(float)

        # sort each column and get index  #每一個目標去由小排到大，並回傳index
        I = np.argsort(_F, axis=0, kind='mergesort')

        # sort the objective space values for the whole matrix  #這個學起來!! 好用
        _F = _F[I, np.arange(n_obj)]  # 可以問看看概念

        # calculate the distance from each point to the last and next    #要再想一下為什麼要用這個!!
        dist = np.row_stack([_F, np.full(n_obj, np.inf)]) - np.row_stack([np.full(n_obj, -np.inf), _F])

        # calculate the norm for each objective - set to NaN if all values are equal
        norm = np.max(_F, axis=0) - np.min(_F, axis=0)  # '''f1max - f1min 跟 f2max -f2min'''
        norm[norm == 0] = np.nan

        # prepare the distance to last and next vectors
        dist_to_last, dist_to_next = dist, np.copy(dist)
        dist_to_last, dist_to_next = dist_to_last[:-1] / norm, dist_to_next[1:] / norm

        # _________________20220216 新增 權重於 NSGAII_____________________________
        for i in range(n_obj):
            dist_to_last[:, i] = dist_to_last[:, i] * Weights[i]
            dist_to_next[:, i] = dist_to_next[:, i] * Weights[i]

        # if we divide by zero because all values in one columns are equal replace by none
        dist_to_last[np.isnan(dist_to_last)] = 0.0
        dist_to_next[np.isnan(dist_to_next)] = 0.0

        # sum up the distance to next and last and norm by objectives - also reorder from sorted list
        J = np.argsort(I, axis=0)  # 不太懂幹嘛在argsort 一次(轉回來??)
        _cd = np.sum(dist_to_last[J, np.arange(n_obj)] + dist_to_next[J, np.arange(n_obj)], axis=1) / n_obj

        # save the final vector which sets the crowding distance for duplicates to zero to be eliminated
        crowding = np.zeros(n_points)
        crowding[is_unique] = _cd

    # crowding[np.isinf(crowding)] = 1e+14
    return crowding


# 3. Cacaluate Paroto_Optimial Front
def ValuesCal(Total_value, pop_size):
    temp_Paroto_Optimial = fast_non_dominated_sort(Total_value);  # 計算Front線
    front = np.zeros(pop_size * 2, dtype=int);  # 再檢查一下pop_size 是否要*2
    for i in range(len(temp_Paroto_Optimial)):  # 將front(前緣線)把值填進去
        for j in temp_Paroto_Optimial[i]:
            front[j] = i + 1;
            '''故意加1，用意是可以方便讓我們看前緣線是1開始 不是從0開始!!!'''
    fronts = copy.copy(front)
    Total_value["Front_value"] = pd.DataFrame(front);  # 會對應front值
    crowding_of_front = np.zeros(pop_size * 2);

    for k, fronts in enumerate(temp_Paroto_Optimial):  # 先按照 0,1,2,... 前緣線，同一條前緣線就去計算crowding_distance
        # calculate the crowding distance of the front
        crowding_of_front[fronts] = calc_crowding_distance(Total_value.iloc[fronts, :-1],
                                                           Weights)  # 不要包含Front_value值 (注意式iloc 還是 loc)
    # =============================================================================
    #     crowding_of_front[np.isinf(crowding_of_front)] = -1      #遇到inf 就直接等於-1 就是代表要直接保留住 (嘗試改成99999999999)
    # =============================================================================
    ##3.1 Compare there front and distenace values
    Total_value["Crowding_Distance_Value"] = pd.DataFrame(crowding_of_front)
    return Total_value


# 4 Cacaluate Rank
def RankforNSGAII(Total_value):
    cols = ['Front_value', 'Crowding_Distance_Value'];  # 設定要取行名
    # tups = Total_value[cols].sort_values(cols,ascending = [True,False]).apply(tuple,1) #排名 (由大到小)
    Total_valu = Total_value[cols].sort_values(cols, ascending=[True, False])
    # f,i = pd.factorize(tups)
    # factorized = pd.Series(f+1,tups.index).rank(ascending = False , method = 'min')
    return Total_value


# Total_value.assign(Rank = factorized)
# =============================================================================
# def Forward_Insertion_Method_Improved(Zip_OrderNumber_Stime_Endtime, Improve_Interation, temp_Select_Fitness_Df1):
#     # %%改善前插法
#     '''處理子群體(Offstring)'''
#     '''要做檢查!!'''
#     '''將Zip_OrderNumber_Stime_Endtime 的 每個Job 訂單index ArrivalTime、PT、產品料號拆開，為了後續前插法'''
#     best_temp_Mc_Corres_To_Job_Index, best_temp_Mc_Corres_To_Job_ArrivalTime, best_temp_Mc_Corres_To_Job_Pt, best_temp_Mc_PPNumber = Thebest_Machine_Corres_To_Job_ArrivalTime_Pt_PPN(
#         Zip_OrderNumber_Stime_Endtime)
#
#     '''轉換成[(訂單Index,到達時間,加工時間)]'''
#     re_best = Thebest_OrderNumber_St_Pt_Zip(Name=best_temp_Mc_Corres_To_Job_Index,
#                                             ArrivalT=best_temp_Mc_Corres_To_Job_ArrivalTime,
#                                             Pt=best_temp_Mc_Corres_To_Job_Pt)
#
#     re_best_Zip = Thebest_Zip_OrderNumber_Stime_Endtime_converList(re_best)
#
#     deep_Df_re_best_Zip = copy.deepcopy(re_best_Zip)
#
#     for re_arrange_Mc in range(len(re_best)):
#         '''N台機台'''
#         # Timestamp_Limit = 1598889600 #2020/9/1
#         Machine_struct_time = time.strptime(Machine_Start.at[re_arrange_Mc, "機台開始時間"], "%Y/%m/%d %H:%M")  # 轉成時間元組
#         '''下面Machine_time_stamp是每台機台的起始時間，變動會隨著停機開始時間csv檔'''
#         ''''目前都從0開始 只能先int(time.mktime(Machine_struct_time)) - int(time.mktime(Machine_struct_time))'''
#         Machine_time_stamp_First = int(time.mktime(Machine_struct_time)) - int(time.mktime(Machine_struct_time))
#
#         '''從第二個間閣開始要尋找插入的ArrivalTime 的相對位置'''
#         Machine_time_stamp_Second = int(time.mktime(Machine_struct_time))
#
#         for re_arrange_Job in range(len(re_best[re_arrange_Mc]) + 1):
#             '''間格˙'''
#             # print("第%s機台，整理 %s Job" %(re_arrange_Mc,re_arrange_Job))
#             for inspection_Job in range(re_arrange_Job + 1, len(re_best[re_arrange_Mc])):
#                 # print(re_arrange_Job,inspection_Job)
#
#                 if re_arrange_Job == 0:
#                     '''跑第0個間格'''
#
#                     insert_Job_Pt_day = re_best_Zip[re_arrange_Mc][inspection_Job][2] / (24 * 60 * 60)  # 加工時間換算成天數
#                     # print(insert_Job_Pt_day)
#                     insert_Job_Pt_second = re_best_Zip[re_arrange_Mc][inspection_Job][2]  # 加工時間換算成秒數
#                     # print(insert_Job_Pt_second)
#
#                     Setup_OrderName_PPN = Product_Part_Number[re_best_Zip[re_arrange_Mc][inspection_Job][0]]  # 目前產品料號
#
#                     Setup_RightOrderName_PPN = Product_Part_Number[
#                         Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]]  # 下一個的產品料號
#
#                     if Setup_OrderName_PPN != Setup_RightOrderName_PPN:
#                         '''判斷目前產品料號與右邊(下一個)產品料號是否不同，若不同則需加上換模時間(秒數)'''
#                         SetupRight_PPN_Time = Mould_dict[Setup_OrderName_PPN + Setup_RightOrderName_PPN] * 60
#                     else:
#                         SetupRight_PPN_Time = 0
#
#                     insert_JobPt_and_setup = (insert_Job_Pt_second + SetupRight_PPN_Time) / (24 * 60 * 60)
#                     insert_Job_Arrival = re_best_Zip[re_arrange_Mc][inspection_Job][1]
#
#                     '''最晚開始時間，若超過最晚開始一定會延遲(關鍵)'''
#                     Latest_start_time = Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][
#                                             1] - SetupRight_PPN_Time - insert_Job_Pt_second
#                     if Latest_start_time < 0:
#                         '''若最晚開始時間 < 0 代表可立即開始'''
#                         Latest_start_time = 0
#                     '''判斷get_between_days(Machine_time_stamp_First , Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][1]) 此間格天數 是否 大於 要插入的Job 加工天數 而且要符合 插入的 ArrivalTime 必須要 <= 最晚開始時間， 這樣才可做前插法'''
#                     if get_between_days(Machine_time_stamp_First,
#                                         Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][1]) >= round(
#                             insert_JobPt_and_setup, 3) and (insert_Job_Arrival <= Latest_start_time):
#                         insert_Job_Start = max(insert_Job_Arrival, Machine_time_stamp_First)
#                         insert_Job_End = insert_Job_Start + insert_Job_Pt_second
#                         Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].pop(inspection_Job)
#                         # Machine_Corresponds_To_Job_Module[0][re_arrange_Mc].pop(inspection_Job)
#
#                         Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].insert(re_arrange_Job, [
#                             re_best_Zip[re_arrange_Mc][inspection_Job][0], insert_Job_Start, insert_Job_End])
#
#                         # 下面五列程式碼是更正前插法的順序
#                         temp = copy.deepcopy(re_best_Zip[re_arrange_Mc][inspection_Job])
#                         re_best_Zip[re_arrange_Mc].pop(inspection_Job)
#                         re_best_Zip[re_arrange_Mc].insert(re_arrange_Job, temp)
#
#                         best_temp_Mc_PPNumber[re_arrange_Mc].pop(inspection_Job)
#                         best_temp_Mc_PPNumber[re_arrange_Mc].insert(re_arrange_Job, Order[
#                             Order.index == Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]][
#                             "產品料號"].values[0])
#
#                         break
#                     else:
#                         continue
#                 else:
#                     '''跑第一個間格'''
#                     insert_Job_Pt_day = re_best_Zip[re_arrange_Mc][inspection_Job][2] / (24 * 60 * 60)
#                     # print(insert_Job_Pt_day)
#                     insert_Job_Pt_second = re_best_Zip[re_arrange_Mc][inspection_Job][2]
#                     # print(insert_Job_Pt_second)
#                     '''目前料號 與 上一個料號名稱'''
#                     Setup_OrderName_PPN = Product_Part_Number[re_best_Zip[re_arrange_Mc][inspection_Job][0]]
#                     Setup_LeftOrderName_PPN = Product_Part_Number[
#                         Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job - 1][0]]
#
#                     '''計算目前料號 與 上一個料號換模時間'''
#                     if Setup_OrderName_PPN != Setup_LeftOrderName_PPN:
#                         SetupLeft_PPN_Time = Mould_dict[Setup_OrderName_PPN + Setup_LeftOrderName_PPN] * 60
#                     else:
#                         SetupLeft_PPN_Time = 0
#
#                     '''目前料號 與 下一個料號名稱'''
#                     Setup_RightOrderName_PPN = Product_Part_Number[
#                         Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]]
#
#                     '''計算目前料號 與 下一個料號換模時間'''
#                     if Setup_OrderName_PPN != Setup_RightOrderName_PPN:
#                         SetupRight_PPN_Time = Mould_dict[Setup_OrderName_PPN + Setup_RightOrderName_PPN] * 60
#                     else:
#                         SetupRight_PPN_Time = 0
#
#                     insert_JobPt_and_setup = (SetupLeft_PPN_Time + insert_Job_Pt_second + SetupRight_PPN_Time) / (
#                                 24 * 60 * 60)
#
#                     insert_Job_Arrival = Arrival_Time[
#                                              re_best_Zip[re_arrange_Mc][inspection_Job][0]] - Machine_time_stamp_Second
#                     if insert_Job_Arrival < 0:
#                         insert_Job_Arrival = 0
#
#                     Latest_start_time = Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][
#                                             1] - SetupRight_PPN_Time - insert_Job_Pt_second
#                     if Latest_start_time < 0:
#                         Latest_start_time = 0
#                     if get_between_days(Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job - 1][2],
#                                         Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][1]) >= round(
#                             insert_JobPt_and_setup, 3) and (insert_Job_Arrival <= Latest_start_time):
#                         '''判斷要前插的空間是否夠插入，且插入的ArrvalTime必須<=最晚開始時間，否則無法插入'''
#                         insert_Job_Start = max(insert_Job_Arrival,
#                                                Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job - 1][
#                                                    2] + SetupLeft_PPN_Time)
#                         insert_Job_End = insert_Job_Start + insert_Job_Pt_second
#
#                         Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].pop(inspection_Job)
#
#                         Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].insert(re_arrange_Job, [
#                             re_best_Zip[re_arrange_Mc][inspection_Job][0], insert_Job_Start, insert_Job_End])
#
#                         temp = copy.deepcopy(re_best_Zip[re_arrange_Mc][inspection_Job])
#                         re_best_Zip[re_arrange_Mc].pop(inspection_Job)
#                         re_best_Zip[re_arrange_Mc].insert(re_arrange_Job, temp)
#
#                         best_temp_Mc_PPNumber[re_arrange_Mc].pop(inspection_Job)
#                         best_temp_Mc_PPNumber[re_arrange_Mc].insert(re_arrange_Job, Order[
#                             Order.index == Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]][
#                             "產品料號"].values[0])
#
#                         break
#                     else:
#                         continue
#
#     '''要增加修改'''
#     '''修正開始時間結束時間'''
#
#     Setup_count = 0
#     Setup_times = 0  # 把這個取消掉 Sum_Tardiness 就跑出來  超怪的!!
#     for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime) - 1):
#         Deep_re_best_Zip_DF = pd.DataFrame(deep_Df_re_best_Zip[Nb_Of_Machine_index], columns=["訂單編號", "最早到達時間", "加工時間"])
#         for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
#             if Nb_Of_Job_index == 0:
#                 continue
#             else:
#
#                 if best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index - 1] != \
#                         best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]:
#                     Setup = Mould_dict[best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index - 1] +
#                                        best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]] * 60
#                     Setup_count += 1
#                     Setup_times += Setup
#                 else:
#                     Setup = 0
#
#                 Pt = Process_Time[re_best_Zip[Nb_Of_Machine_index][Nb_Of_Job_index][0]]
#                 Job_ArrivalTime = Deep_re_best_Zip_DF[
#                     Deep_re_best_Zip_DF["訂單編號"] == re_best_Zip[Nb_Of_Machine_index][Nb_Of_Job_index][0]][
#                     "最早到達時間"].values[0]
#                 Job_Start = max(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index - 1][2] + Setup,
#                                 Job_ArrivalTime)
#                 Job_End = Job_Start + Pt
#                 '''從第二個Job開始的開始時間'''
#                 Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][1] = Job_Start
#                 Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2] = Job_End
#                 # print(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][1] ,Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2] )
#
#     '''0520前插法後的新的 Total_completion_time(Flow time) 前插法後的ArrivalTime需再確認'''
#     # =============================================================================
#     # Total_completion_time = Thebest_List_to_Np(Zip_OrderNumber_Stime_Endtime)
#     # =============================================================================
#     Each_Total_EndTime = Thebest_List_to_Np(Zip_OrderNumber_Stime_Endtime)  # ˊ轉換成矩陣模式，方便做下一步np.sum的運算
#     Each_temp_Mc_Co_To_Job_Arrival_matrix = List_to_Np(best_temp_Mc_Corres_To_Job_ArrivalTime)
#     Each_Total_ArrivalTime = np.sum(Each_temp_Mc_Co_To_Job_Arrival_matrix)
#     Total_completion_time = Each_Total_EndTime - Each_Total_ArrivalTime
#
#     '''處理目標，each_Total_completion_Tardiness_Time = Total_completion_time + Tardiness * w'''
#     Sum_Tardiness = 0
#     Sum_Tardiness_cost = 0
#     for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime) - 1):
#         Machine_struct_time = time.strptime(Machine_Start.at[Nb_Of_Machine_index, "機台開始時間"], "%Y/%m/%d %H:%M")  # 轉成時間元組
#         '''下面Machine_time_stamp是每台機台的起始時間，變動會隨著停機開始時間csv檔'''
#         Machine_time_stamp = int(time.mktime(Machine_struct_time))
#
#         for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
#             Completion_time = Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2]
#
#             # Due_Date_time = int(time.mktime(time.strptime(Order.at[Order[Order["訂單編號"] ==  Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0]].index[0],'交期'], "%Y/%m/%d %H:%M")))-Machine_time_stamp
#             Due_Date_time = int(time.mktime(time.strptime(Order.at[Order[Order.index == Zip_OrderNumber_Stime_Endtime[
#                 Nb_Of_Machine_index][Nb_Of_Job_index][0]].index[0], '交期'], "%Y/%m/%d %H:%M"))) - Machine_time_stamp
#
#             if Completion_time - Due_Date_time > 0:
#                 Sum_Tardiness += (Completion_time - Due_Date_time)
#                 '''計算此單延誤損失成本'''
#                 # Sum_Tardiness_cost += Order.at[Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0],"單張延誤成本"]
#
#     Sum_Tardiness = Sum_Tardiness / (60 * 60 * 24)
#     Setup_times = Setup_times / (60 * 60 * 24)
#
#     # =============================================================================
#     #     '''0531下面目標式用Setup_times'''
#     #     each_Single_Target = Machine_index * Penalty_value_weight_Machine + Sum_Tardiness * Penalty_value_weight_TT + Setup_times * Penalty_value_weight_Setup
#     # =============================================================================
#     temp_Select_Fitness_Df1[2], temp_Select_Fitness_Df1[3] = Sum_Tardiness, Setup_times
#
#     # Zip_OrderNumber_Stime_Endtime[-1] =  np.array([Improve_Interation,Machine_index,Sum_Tardiness,Setup_times,each_Single_Target])
#
#     return Zip_OrderNumber_Stime_Endtime, temp_Select_Fitness_Df1
# =============================================================================
def Forward_Insertion_Method_Improved(q1, q2, arg1, Zip_OrderNumber_Stime_Endtime, temp_Select_Fitness_Df1):
    # %%改善前插法
    '''處理子群體(Offstring)'''
    '''要做檢查!!'''
    '''將Zip_OrderNumber_Stime_Endtime 的 每個Job 訂單index ArrivalTime、PT、產品料號拆開，為了後續前插法'''
    # pdb.set_trace()
    # for i in range(arg1, len(temp_Select_Fitness_Df1), nb_threads):
    #     print(i)
    # pdb.set_trace()
    # pdb.set_trace()
    best_temp_Mc_Corres_To_Job_Index, best_temp_Mc_Corres_To_Job_ArrivalTime, best_temp_Mc_Corres_To_Job_Pt, best_temp_Mc_PPNumber = Thebest_Machine_Corres_To_Job_ArrivalTime_Pt_PPN(
        Zip_OrderNumber_Stime_Endtime)

    '''轉換成[(訂單Index,到達時間,加工時間)]'''
    re_best = Thebest_OrderNumber_St_Pt_Zip(Name=best_temp_Mc_Corres_To_Job_Index,
                                            ArrivalT=best_temp_Mc_Corres_To_Job_ArrivalTime,
                                            Pt=best_temp_Mc_Corres_To_Job_Pt)

    re_best_Zip = Thebest_Zip_OrderNumber_Stime_Endtime_converList(re_best)

    deep_Df_re_best_Zip = copy.deepcopy(re_best_Zip)
    # pdb.set_trace()
    for re_arrange_Mc in range(len(re_best)):
        '''N台機台'''
        # Timestamp_Limit = 1598889600 #2020/9/1
        Machine_struct_time = time.strptime(Machine_Start.at[re_arrange_Mc, "機台開始時間"], "%Y/%m/%d %H:%M")  # 轉成時間元組
        '''下面Machine_time_stamp是每台機台的起始時間，變動會隨著停機開始時間csv檔'''
        ''''目前都從0開始 只能先int(time.mktime(Machine_struct_time)) - int(time.mktime(Machine_struct_time))'''
        Machine_time_stamp_First = int(time.mktime(Machine_struct_time)) - int(time.mktime(Machine_struct_time))

        '''從第二個間閣開始要尋找插入的ArrivalTime 的相對位置'''
        Machine_time_stamp_Second = int(time.mktime(Machine_struct_time))

        for re_arrange_Job in range(len(re_best[re_arrange_Mc]) + 1):
            '''間格˙'''
            # print("第%s機台，整理 %s Job" %(re_arrange_Mc,re_arrange_Job))
            for inspection_Job in range(re_arrange_Job + 1, len(re_best[re_arrange_Mc])):
                # print(re_arrange_Job,inspection_Job)

                if re_arrange_Job == 0:
                    '''跑第0個間格'''

                    insert_Job_Pt_day = re_best_Zip[re_arrange_Mc][inspection_Job][2] / (24 * 60 * 60)  # 加工時間換算成天數
                    # print(insert_Job_Pt_day)
                    insert_Job_Pt_second = re_best_Zip[re_arrange_Mc][inspection_Job][2]  # 加工時間換算成秒數
                    # print(insert_Job_Pt_second)

                    Setup_OrderName_PPN = Product_Part_Number[
                        re_best_Zip[re_arrange_Mc][inspection_Job][0]]  # 目前產品料號

                    Setup_RightOrderName_PPN = Product_Part_Number[
                        Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]]  # 下一個的產品料號

                    if Setup_OrderName_PPN != Setup_RightOrderName_PPN:
                        '''判斷目前產品料號與右邊(下一個)產品料號是否不同，若不同則需加上換模時間(秒數)'''
                        SetupRight_PPN_Time = Mould_dict[Setup_OrderName_PPN + Setup_RightOrderName_PPN] * 60
                    else:
                        SetupRight_PPN_Time = 0

                    insert_JobPt_and_setup = (insert_Job_Pt_second + SetupRight_PPN_Time) / (24 * 60 * 60)
                    insert_Job_Arrival = re_best_Zip[re_arrange_Mc][inspection_Job][1]

                    '''最晚開始時間，若超過最晚開始一定會延遲(關鍵)'''
                    Latest_start_time = Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][
                                            1] - SetupRight_PPN_Time - insert_Job_Pt_second
                    if Latest_start_time < 0:
                        '''若最晚開始時間 < 0 代表可立即開始'''
                        Latest_start_time = 0
                    '''判斷get_between_days(Machine_time_stamp_First , Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job][1]) 此間格天數 是否 大於 要插入的Job 加工天數 而且要符合 插入的 ArrivalTime 必須要 <= 最晚開始時間， 這樣才可做前插法'''
                    if get_between_days(Machine_time_stamp_First,
                                        Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][1]) >= round(
                        insert_JobPt_and_setup, 3) and (insert_Job_Arrival <= Latest_start_time):
                        insert_Job_Start = max(insert_Job_Arrival, Machine_time_stamp_First)
                        insert_Job_End = insert_Job_Start + insert_Job_Pt_second
                        Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].pop(inspection_Job)
                        # Machine_Corresponds_To_Job_Module[0][re_arrange_Mc].pop(inspection_Job)

                        Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].insert(re_arrange_Job, [
                            re_best_Zip[re_arrange_Mc][inspection_Job][0], insert_Job_Start, insert_Job_End])

                        # 下面五列程式碼是更正前插法的順序
                        temp = copy.deepcopy(re_best_Zip[re_arrange_Mc][inspection_Job])
                        re_best_Zip[re_arrange_Mc].pop(inspection_Job)
                        re_best_Zip[re_arrange_Mc].insert(re_arrange_Job, temp)

                        best_temp_Mc_PPNumber[re_arrange_Mc].pop(inspection_Job)
                        best_temp_Mc_PPNumber[re_arrange_Mc].insert(re_arrange_Job, Order[
                            Order.index == Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]][
                            "產品料號"].values[0])

                        break
                    else:
                        continue
                else:
                    '''跑第一個間格'''
                    insert_Job_Pt_day = re_best_Zip[re_arrange_Mc][inspection_Job][2] / (24 * 60 * 60)
                    # print(insert_Job_Pt_day)
                    insert_Job_Pt_second = re_best_Zip[re_arrange_Mc][inspection_Job][2]
                    # print(insert_Job_Pt_second)
                    '''目前料號 與 上一個料號名稱'''
                    Setup_OrderName_PPN = Product_Part_Number[re_best_Zip[re_arrange_Mc][inspection_Job][0]]
                    Setup_LeftOrderName_PPN = Product_Part_Number[
                        Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job - 1][0]]

                    '''計算目前料號 與 上一個料號換模時間'''
                    if Setup_OrderName_PPN != Setup_LeftOrderName_PPN:
                        SetupLeft_PPN_Time = Mould_dict[Setup_OrderName_PPN + Setup_LeftOrderName_PPN] * 60
                    else:
                        SetupLeft_PPN_Time = 0

                    '''目前料號 與 下一個料號名稱'''
                    Setup_RightOrderName_PPN = Product_Part_Number[
                        Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]]

                    '''計算目前料號 與 下一個料號換模時間'''
                    if Setup_OrderName_PPN != Setup_RightOrderName_PPN:
                        SetupRight_PPN_Time = Mould_dict[Setup_OrderName_PPN + Setup_RightOrderName_PPN] * 60
                    else:
                        SetupRight_PPN_Time = 0

                    insert_JobPt_and_setup = (SetupLeft_PPN_Time + insert_Job_Pt_second + SetupRight_PPN_Time) / (
                            24 * 60 * 60)

                    insert_Job_Arrival = Arrival_Time[re_best_Zip[re_arrange_Mc][inspection_Job][
                        0]] - Machine_time_stamp_Second
                    if insert_Job_Arrival < 0:
                        insert_Job_Arrival = 0

                    Latest_start_time = Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][
                                            1] - SetupRight_PPN_Time - insert_Job_Pt_second
                    if Latest_start_time < 0:
                        Latest_start_time = 0
                    if get_between_days(Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job - 1][2],
                                        Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][1]) >= round(
                        insert_JobPt_and_setup, 3) and (insert_Job_Arrival <= Latest_start_time):
                        '''判斷要前插的空間是否夠插入，且插入的ArrvalTime必須<=最晚開始時間，否則無法插入'''
                        insert_Job_Start = max(insert_Job_Arrival,
                                               Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job - 1][
                                                   2] + SetupLeft_PPN_Time)
                        insert_Job_End = insert_Job_Start + insert_Job_Pt_second

                        Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].pop(inspection_Job)

                        Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].insert(re_arrange_Job, [
                            re_best_Zip[re_arrange_Mc][inspection_Job][0], insert_Job_Start, insert_Job_End])

                        temp = copy.deepcopy(re_best_Zip[re_arrange_Mc][inspection_Job])
                        re_best_Zip[re_arrange_Mc].pop(inspection_Job)
                        re_best_Zip[re_arrange_Mc].insert(re_arrange_Job, temp)

                        best_temp_Mc_PPNumber[re_arrange_Mc].pop(inspection_Job)
                        best_temp_Mc_PPNumber[re_arrange_Mc].insert(re_arrange_Job, Order[
                            Order.index == Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]][
                            "產品料號"].values[0])

                        break
                    else:
                        continue

    '''要增加修改'''
    '''修正開始時間結束時間'''

    Setup_count = 0
    Setup_times = 0  # 把這個取消掉 Sum_Tardiness 就跑出來  超怪的!!
    for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime) - 1):
        Deep_re_best_Zip_DF = pd.DataFrame(deep_Df_re_best_Zip[Nb_Of_Machine_index],
                                           columns=["訂單編號", "最早到達時間", "加工時間"])
        for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
            if Nb_Of_Job_index == 0:
                continue
            else:

                if best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index - 1] != \
                        best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]:
                    Setup = Mould_dict[best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index - 1] +
                                       best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]] * 60
                    Setup_count += 1
                    Setup_times += Setup
                else:
                    Setup = 0

                Pt = Process_Time[re_best_Zip[Nb_Of_Machine_index][Nb_Of_Job_index][0]]
                Job_ArrivalTime = Deep_re_best_Zip_DF[
                    Deep_re_best_Zip_DF["訂單編號"] == re_best_Zip[Nb_Of_Machine_index][Nb_Of_Job_index][0]][
                    "最早到達時間"].values[0]
                Job_Start = max(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index - 1][2] + Setup,
                                Job_ArrivalTime)
                Job_End = Job_Start + Pt
                '''從第二個Job開始的開始時間'''
                Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][1] = Job_Start
                Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2] = Job_End
                # print(Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index][Nb_Of_Job_index][1] ,Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2] )

    '''0520前插法後的新的 Total_completion_time(Flow time) 前插法後的ArrivalTime需再確認'''
    # =============================================================================
    # Total_completion_time = Thebest_List_to_Np(Zip_OrderNumber_Stime_Endtime)
    # =============================================================================
    # pdb.set_trace()
    Each_Total_EndTime = Thebest_List_to_Np(Zip_OrderNumber_Stime_Endtime)  # ˊ轉換成矩陣模式，方便做下一步np.sum的運算
    Each_temp_Mc_Co_To_Job_Arrival_matrix = List_to_Np(best_temp_Mc_Corres_To_Job_ArrivalTime)
    Each_Total_ArrivalTime = np.sum(Each_temp_Mc_Co_To_Job_Arrival_matrix)
    Total_completion_time = Each_Total_EndTime - Each_Total_ArrivalTime

    '''處理目標，each_Total_completion_Tardiness_Time = Total_completion_time + Tardiness * w'''
    Sum_Tardiness = 0
    Sum_Tardiness_cost = 0
    for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime) - 1):
        Machine_struct_time = time.strptime(Machine_Start.at[Nb_Of_Machine_index, "機台開始時間"],
                                            "%Y/%m/%d %H:%M")  # 轉成時間元組
        '''下面Machine_time_stamp是每台機台的起始時間，變動會隨著停機開始時間csv檔'''
        Machine_time_stamp = int(time.mktime(Machine_struct_time))

        for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
            Completion_time = Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2]

            # Due_Date_time = int(time.mktime(time.strptime(Order.at[Order[Order["訂單編號"] ==  Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0]].index[0],'交期'], "%Y/%m/%d %H:%M")))-Machine_time_stamp
            Due_Date_time = int(time.mktime(time.strptime(Order.at[Order[Order.index ==
                                                                         Zip_OrderNumber_Stime_Endtime[
                                                                             Nb_Of_Machine_index][Nb_Of_Job_index][
                                                                             0]].index[0], '交期'],
                                                          "%Y/%m/%d %H:%M"))) - Machine_time_stamp

            if Completion_time - Due_Date_time > 0:
                Sum_Tardiness += (Completion_time - Due_Date_time)
                '''計算此單延誤損失成本'''
                # Sum_Tardiness_cost += Order.at[Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0],"單張延誤成本"]
    # pdb.set_trace()

    Sum_Tardiness = Sum_Tardiness / (60 * 60 * 24)
    Setup_times = Setup_times / (60 * 60 * 24)

    # =============================================================================
    #     '''0531下面目標式用Setup_times'''
    #     each_Single_Target = Machine_index * Penalty_value_weight_Machine + Sum_Tardiness * Penalty_value_weight_TT + Setup_times * Penalty_value_weight_Setup
    # =============================================================================
    temp_Select_Fitness_Df1[2], temp_Select_Fitness_Df1[3] = Sum_Tardiness, Setup_times

    # Zip_OrderNumber_Stime_Endtime[-1] =  np.array([Improve_Interation,Machine_index,Sum_Tardiness,Setup_times,each_Single_Target])
    Zip_OrderNumber_Stime_Endtime[-1] = temp_Select_Fitness_Df1.tolist()
    # q.put(Zip_OrderNumber_Stime_Endtime)
    # print(temp_Select_Fitness_Df1)
    q1.put(temp_Select_Fitness_Df1)
    q2.put(Zip_OrderNumber_Stime_Endtime)
    # print(q1)
    # print(q2)
    # return Zip_OrderNumber_Stime_Endtime, temp_Select_Fitness_Df1


# %%
# ==============一個Job總共有幾個機台可以用=============#
if __name__ == '__main__':
    Wip_canMachines = []
    The_Number_Of_Mc = []
    '''Wip_canMachines 為 每個Job可用那些機台(要依照BOM表去選機有2種型號)'''

    for Order_index in range(len(Order)):
        '''按照Order表的 Order_index順序去搜尋'''
        temp = len(dict_canRunTool[Order.at[Order_index, "orderBomKey"]])
        Wip_canMachines.append(temp)
        The_Number_Of_Mc.append(temp)
    '''每個訂單所可用的機台'''

    # The_Number_Of_Mc =[]
    # '''The_Number_Of_Mc 為  每個Job可用的機台數'''

    # for i in range(len(Wip_canMachines)):
    #     temp = 0
    #     for j in range(len(Wip_canMachines[i])):
    #         # print(j)
    #         if True:
    #             temp +=1
    #     The_Number_Of_Mc.append(temp)

    # %%所有機台名稱
    All_Machines_Name = list(map(str, Machine_Start["設備/資源編號"]))

    # Each_Iteration_Objective1, Each_Iteration_Objective2, Each_Iteration_Objective3 = [],[],[] #儲存每一代最好的Objective
    Each_Iteration_Objective = pd.DataFrame()
    # Each_Iteration_Objective = []
    t1_start = time.time();  # 計算終止條件的開始時間

    Start_Iterations = 0  # 從第0代開始計算
    new_gen_arr = 0  # 初始化
    for k in range(num_iterations):

        if k != 0 and k % 20 == 0:
            '''此限制式是若Chrossover_Rate 低於0.3 則固定為 0.3 ; 反之Mutation_Rate 若高於0.7則固定於0.7'''
            if Chrossover_Rate <= 0.3 and Mutation_Rate >= 0.7:
                Chrossover_Rate = 0.3
                Mutation_Rate = 0.7
            else:
                '''每20代Chrossover_Rate 就 * 1 - 0.05; Mutation_Rate  1 - Chrossover_Rate '''
                Chrossover_Rate *= (1 - 0.02)
                Mutation_Rate = (1 - Chrossover_Rate)

        if k == 0:
            # %%開始建立染色體條數(Ex 把所有chromosome 染色體數*(基因Job*2(包含選機、排序)) 都變成 0 )
            # 若遇到初始染色體是基數就自然多+1條，方便計算交配突變 (例如:設定母體51條，會與52條結果一樣)
            if (Number_Of_Chromosome % 2) == 1:
                Number_Of_Chromosome += 1
                chromosome = np.zeros(shape=(Number_Of_Chromosome * 2, Number_Of_Job * 2))
                # 建立Shuffle打亂機制Index(按照原創染色體數目Index去打亂)
                Temp_Index = np.arange(Number_Of_Chromosome)

                Temp_Shuffle = np.random.shuffle(Temp_Index)
            else:
                chromosome = np.zeros(shape=(Number_Of_Chromosome * 2, Number_Of_Job * 2))

                # 建立Shuffle打亂機制Index(按照原創染色體數目Index去打亂)
                Temp_Index = np.arange(Number_Of_Chromosome)

                Temp_Shuffle = np.random.shuffle(Temp_Index)

            # %% 直接建立多條染色體亂數(Number_Of_Chromosome)(Ex 剛開始 50條母體染色體亂數)
            for i in range(Number_Of_Chromosome):
                temp = np.array(np.random.rand(1, Number_Of_Job * 2))  # 建立一個空的DataFrame，然後下面跑Number_Of_Chromosome次!!
                chromosome[i] = temp
        else:
            # 重新創立母體chromosome大小
            chromosome = np.zeros(shape=(Number_Of_Chromosome * 2, Number_Of_Job * 2))

            for i in range(Number_Of_Chromosome):
                for j in range(len(chromosome[i])):
                    '''從第二代開始把上一代的選擇Number_Of_Chromosome 條 的基因染色體丟回 chromosome內'''
                    chromosome[i][j] = new_gen_arr[i][j]

        # %%整理基因第三部分(因為不希望機台數太小 希望在原來機台數的上下5~10%)

        # st = time.time()
        # for i in range(len(chromosome)):
        #     for  j in range(len(Machine_model_Number)):
        #         # print("第 %s 條染色體, 第 %s 基因 " % (i,j))
        #         a = chromosome[i][Number_Of_Job*2:][j*2:(j+1)*2]
        #         # print('舊的基因',a)
        #         while abs(np.diff(a)) < 0.005:
        #             '''a[:]重新給予新的可用機台基因'''
        #             a[:] = np.array(np.random.rand(len(Machine_model_Number)))
        #             # print("新的資料",a)
        # endt = time.time()
        # # print("整理基因第三部分",endt - st)

        # %%=====================交配突變=====================#
        Temp_Chrossover = np.zeros(shape=(np.rint((Number_Of_Chromosome * Chrossover_Rate)).astype(int),
                                          Number_Of_Job * 2))  # 建立Temp_Chrossover空的array(Number_Of_Chromosome * Chrossover_Rate)
        Temp_Chrossover_A_B = np.zeros(shape=(np.rint((Number_Of_Chromosome * Chrossover_Rate)).astype(int),
                                              Number_Of_Job * 2))  # Temp_Chrossover_A_B空的array(Number_Of_Chromosome * Chrossover_Rate) (跟上面差異為交配部分是放在Temp_Chrossover_A_B)
        # Temp_Mutation = np.zeros(shape=(np.rint((Number_Of_Chromosome * Mutation_Rate)).astype(int),Number_Of_Job*2))  #建立Temp_Mutation空的array(Number_Of_Chromosome * Mutation_Rate)
        Temp_Mutation = np.zeros(shape=(np.rint(Number_Of_Chromosome - len(Temp_Chrossover)).astype(int),
                                        Number_Of_Job * 2))  # 建立Temp_Mutation空的array(Number_Of_Chromosome * Mutation_Rate)
        # %%未交配先將 Temp_Index 打亂過的 去回推原來染色體 丟進Temp_Chrossover ，共有len(Temp_Index)*Chrossover_Rate條
        for i in range(np.rint((len(Temp_Index) * Chrossover_Rate)).astype(int)):
            '''共有len(Temp_Index)*Chrossover_Rate 條染色體需要做交配(小心可能在141代交配0.5587 突變 0.4413遇到有一列沒辦法去分配亂數)'''
            temp = chromosome[Temp_Index[i]]
            Temp_Chrossover[i] = temp

        # %% '''執行交配'''
        if len(Temp_Chrossover_A_B) % 2 == 0:
            '''判斷Temp_Chrossover_A_B格子數是否為『偶數』，判斷完再丟入Function ：Crossover及　multi_point_crossover'''
            for i in range(0, int((len(Temp_Index) * Chrossover_Rate)), 2):
                C, D = Temp_Chrossover[i], Temp_Chrossover[i + 1]
                # print(C,D)
                '''上面C,D 和下面 C1,D1 都是一樣染色體只是 丟進Crossover(C,D)為了要求出 X 分割哪兩點去做交配!!'''
                C1, D1, X = Crossover(C, D)
                '''X求出來若是 array([2,16])，丟進下面multi_point_crossover 會變成指交換 [2,16) 不包含 16index，因為index從0開始!!!'''
                A, B = multi_point_crossover(C1, D1, X)

                Temp_Chrossover_A_B[i], Temp_Chrossover_A_B[i + 1] = A, B
        else:
            '''判斷Temp_Chrossover_A_B格子數是否為『奇數』，判斷完再丟入Function ：Crossover及　multi_point_crossover'''
            temp = 0
            for i in range(0, int((len(Temp_Index) * Chrossover_Rate)) - 1, 2):
                # print(i)
                C, D = Temp_Chrossover[i], Temp_Chrossover[i + 1]
                # print(C,D)
                '''上面C,D 和下面 C1,D1 都是一樣染色體只是 丟進Crossover(C,D)為了要求出 X 分割哪兩點去做交配!!'''
                C1, D1, X = Crossover(C, D)
                '''X求出來若是 array([2,16])，丟進下面multi_point_crossover 會變成指交換 [2,16) 不包含 16index，因為index從0開始!!!'''
                A, B = multi_point_crossover(C1, D1, X)

                Temp_Chrossover_A_B[i], Temp_Chrossover_A_B[i + 1] = A, B
            '''下面是因為奇數少一條染色體所以任抓兩條去交配並新增一條回來˙!!'''
            r1, r2 = np.random.choice(range(Number_Of_Chromosome), size=2, replace=False)
            C, D = chromosome[r1], chromosome[r2]
            C1, D1, X = Crossover(C, D)
            A, B = multi_point_crossover(C1, D1, X)

            Temp_Chrossover_A_B[-1] = A

        chromosome[Number_Of_Chromosome:int(
            round(Number_Of_Chromosome * (1 + Chrossover_Rate), 0))] = Temp_Chrossover_A_B  # 子代交配出來的(要放在50~90空格中)
        # %%整理基因第三部分(最後可以取消)

        # st = time.time()
        # for i in range(len(chromosome)):
        #     for  j in range(len(Machine_model_Number)):
        #         # print("第 %s 條染色體, 第 %s 基因 " % (i,j))
        #         a = chromosome[i][Number_Of_Job*2:][j*2:(j+1)*2]
        #         # print('舊的基因',a)
        #         while abs(np.diff(a)) < 0.05:
        #             '''a[:]重新給予新的可用機台基因'''
        #             a[:] = np.array(np.random.rand(len(Machine_model_Number)))
        #             # print("新的資料",a)
        # endt = time.time()
        # # print("整理基因第三部分",endt - st)
        # %%突變
        '''(注意)到最後突變數會變'''
        '''小心這裡要用np.rint () .astype(int)，以免少一條染色體'''
        '''任選10條染色體index'''
        temp = np.random.choice(range(int(Number_Of_Chromosome * (1 + Chrossover_Rate))),
                                size=np.rint(Number_Of_Chromosome * Mutation_Rate).astype(int),
                                replace=False)  # EX：89條(母體及交配)隨機從裡面挑10條
        if len(temp) == 0:
            '''固定Mutatation 不改突變'''
            temp = np.random.choice(range(int(Number_Of_Chromosome * (1 + Chrossover_Rate))),
                                    size=np.rint(10).astype(int), replace=False)  # EX：89條(母體及交配)隨機從裡面挑10條
            for index, values in enumerate(temp):
                Temp_Mutation[index] = chromosome[values]
        else:
            for index, values in enumerate(temp):
                # print(index,values)
                '''先借放要準備突變的染色體在Temp_Mutation'''
                Temp_Mutation[index] = chromosome[values]

            for i in range(len(Temp_Mutation)):
                # print("原本",Temp_Mutation[i])
                '''去Mutation的 Function進行突變'''
                Temp_Mutation = Mutation(i)

        chromosome[int(round(Number_Of_Chromosome * (1 + Chrossover_Rate), 0)):] = Temp_Mutation

        # %%解碼
        #####---------正式解碼----------####
        # =============================================================================
        # '''選機，將機子分成幾等份 (選擇輪盤第幾塊)'''
        # if k == 0:
        #     '''first Interation'''
        #     Slice_of_pie = []  #Slice_of_pie 儲存每個Job被分配到第幾塊大餅
        #     Storage_Opt_Machine = []
        #     # temp_Storage_Opt_Machine = {} # 小小佔存器
        #     for i in range(len(chromosome)):
        #         temp_list = []
        #         temp_Storage_Opt_Machine = {} # 小小佔存器
        #         st = time.time()
        #         '''第幾塊餅，例如第二塊餅，但他的index是1'''
        #         # temp_number_mc = sum([len(x) for x in temp_Storage_Opt_Machine.values()])

        #         '''修正機台數量，從染色體修正'''

        #                 temp = 0
        #                 temp1 = 2
        #                 for index,value in enumerate(Mc_model_Quantity):
        #                     if value == Machine_model_1:
        #                         Available_Machine_model_Gene_St,Available_Machine_model_Gene_End = chromosome[i][Number_Of_Job*2+temp:Number_Of_Job*2+temp1]
        #
        #                         Num_Mc_Start = math.floor( Mc_model_Quantity[value]* Available_Machine_model_Gene_St)            #機台數 * 可用機台數 開始
        #                         Num_Mc_End   = int(np.round(Mc_model_Quantity[value] * Available_Machine_model_Gene_End ))  #機台數 * 可用機台數 結束
        #
        #                         diff_Num_Mc1 = Num_Mc_End - Num_Mc_Start #機台數結束 - 機台數開始 = 可用機台範圍
        #                         if diff_Num_Mc1 < 0:
        #                             '''交換 Start 跟 End'''
        #                             temp = Num_Mc_End
        #                             Num_Mc_End = Num_Mc_Start
        #                             Num_Mc_Start = temp
        #                             diff_Num_Mc1 = Num_Mc_End - Num_Mc_Start
        #                         elif diff_Num_Mc1 == 0:
        #                             '''機台開始與結束一樣，ex 開始58台 結束也58台，那就範圍 開始58 結束 59即可[58:59]'''
        #                             Num_Mc_End += 1
        #                         temp = temp1
        #                         temp1*=2
        #                     else:
        #                         Available_Machine_model_Gene_St,Available_Machine_model_Gene_End = chromosome[i][Number_Of_Job*2+temp:Number_Of_Job*2+temp1]
        #                         '''若某訂單的機型是Machine_model_2，就去增加他的機台型號數量'''
        #                         Num_Mc_Start = math.floor(Mc_model_Quantity[value] * Available_Machine_model_Gene_St)
        #                         Num_Mc_End   = int(np.round(Mc_model_Quantity[value] * Available_Machine_model_Gene_End))
        #                         diff_Num_Mc2 = Num_Mc_End - Num_Mc_Start
        #                         if diff_Num_Mc2 < 0:
        #                             '''交換 Start 跟 End'''
        #                             temp = Num_Mc_End
        #                             Num_Mc_End = Num_Mc_Start
        #                             Num_Mc_Start = temp
        #                             diff_Num_Mc2 = Num_Mc_End - Num_Mc_Start
        #                         elif diff_Num_Mc2 == 0:
        #                             '''機台開始與結束一樣，ex 開始58台 結束也58台，那就範圍 開始58 結束 59即可[58:59]'''
        #                             Num_Mc_End += 1
        #                 Total_Num_Mc = diff_Num_Mc1+diff_Num_Mc2
        #
        #                 while Total_Num_Mc < Mumber_of_machines_limit :
        #                     '''如果跑出來的機台數少於當時設定的Mumber_of_machines_limit  就要重新設定目前的染色體基因'''
        # # =============================================================================
        # #                     # for j in range(len(Mc_model_Quantity)):
        # #                     #     a = chromosome[i][Number_Of_Job*2:][j*2:(j+1)*2]
        # #                     #     # print('舊的基因',a)
        # #                     #     '''a[:]重新給予新的可用機台基因'''
        # #                     #     a[:] = np.array(np.random.rand(len(Machine_model_Number)))
        # # =============================================================================
        #                     '''新的雙點機台基因突變'''
        #                     a = chromosome[i][Number_Of_Job*2:]
        #                     # print('舊的基因',a)
        #                     '''a[:]重新給予新的可用機台基因'''
        #                     # temp_position1,temp_position2 =np.random.choice(range(0,4),size = 2,replace = False)
        #                     # temp_RandomNumber1,temp_RandomNumber2 = [random.random() for i in range(2)]
        #                     a[:] = Mc_Mutation(a)
        #
        #
        #                     temp = 0
        #                     temp1 = 2
        #                     for index,value in enumerate(Mc_model_Quantity):
        #                         if value == Machine_model_1:
        #                             Available_Machine_model_Gene_St,Available_Machine_model_Gene_End = chromosome[i][Number_Of_Job*2+temp:Number_Of_Job*2+temp1]
        #
        #                             Num_Mc_Start = math.floor( Mc_model_Quantity[value]* Available_Machine_model_Gene_St)            #機台數 * 可用機台數 開始
        #                             Num_Mc_End   = int(np.round(Mc_model_Quantity[value] * Available_Machine_model_Gene_End ))  #機台數 * 可用機台數 結束
        #
        #                             diff_Num_Mc1 = Num_Mc_End - Num_Mc_Start #機台數結束 - 機台數開始 = 可用機台範圍
        #                             if diff_Num_Mc1 < 0:
        #                                 '''交換 Start 跟 End'''
        #                                 temp = Num_Mc_End
        #                                 Num_Mc_End = Num_Mc_Start
        #                                 Num_Mc_Start = temp
        #                                 diff_Num_Mc1 = Num_Mc_End - Num_Mc_Start
        #                             elif diff_Num_Mc1 == 0:
        #                                 '''機台開始與結束一樣，ex 開始58台 結束也58台，那就範圍 開始58 結束 59即可[58:59]'''
        #                                 Num_Mc_End += 1
        #                             temp = temp1
        #                             temp1*=2
        #                             # temp_Storage_Opt_Machine[Machine_model_1] = All_Machines_Name[Num_Mc_Start:Num_Mc_End]
        #                             temp_Storage_Opt_Machine[Machine_model_1] = All_Machines_Name[Num_Mc_Start:Num_Mc_End]
        #                         else:
        #                             Available_Machine_model_Gene_St,Available_Machine_model_Gene_End = chromosome[i][Number_Of_Job*2+temp:Number_Of_Job*2+temp1]
        #                             '''若某訂單的機型是Machine_model_2，就去增加他的機台型號數量'''
        #                             Num_Mc_Start = math.floor(Mc_model_Quantity[value] * Available_Machine_model_Gene_St)
        #                             Num_Mc_End   = int(np.round(Mc_model_Quantity[value] * Available_Machine_model_Gene_End))
        #                             diff_Num_Mc2 = Num_Mc_End - Num_Mc_Start
        #                             if diff_Num_Mc2 < 0:
        #                                 '''交換 Start 跟 End'''
        #                                 temp = Num_Mc_End
        #                                 Num_Mc_End = Num_Mc_Start
        #                                 Num_Mc_Start = temp
        #                                 diff_Num_Mc2 = Num_Mc_End - Num_Mc_Start
        #                             elif diff_Num_Mc2 == 0:
        #                                 '''機台開始與結束一樣，ex 開始58台 結束也58台，那就範圍 開始58 結束 59即可[58:59]'''
        #                                 Num_Mc_End += 1
        #
        #                             temp_Storage_Opt_Machine[Machine_model_2] = All_Machines_Name[Mc_model_Quantity[Machine_model_1] + Num_Mc_Start -1 : Mc_model_Quantity[Machine_model_1] + Num_Mc_End]
        #
        #                     Total_Num_Mc = diff_Num_Mc1+diff_Num_Mc2   #重新計算總共多少台機台
        #
        #                 endt = time.time()
        #                 # print("第%s條染色體，計算Part3運算速度：%s" % (i,endt -st))
        #
        #
        #                 st = time.time()
        #                 length_temp_strge_opt_machine_model1 = (1/len(temp_Storage_Opt_Machine[Machine_model_1]))
        #                 length_temp_strge_opt_machine_model2 = (1/len(temp_Storage_Opt_Machine[Machine_model_2]))
        #                 for Order_index in range(len(chromosome[i][:Number_Of_Job])):
        #                     '''判斷是哪個機型!!(切忌少用iloc，at 速度更快!!)'''
        #                     if Order.at[Order_index, '設備/資源'] == Machine_model_1:
        #                     # if Order.iloc[Order_index]['設備/資源'] == Machine_model_1:
        #                         temp  = math.ceil(chromosome[i][Order_index] / length_temp_strge_opt_machine_model1)
        #                         # print(temp)
        #
        #                     else:
        #                         temp  = math.ceil(chromosome[i][Order_index] / length_temp_strge_opt_machine_model2) #主要計算選機部分式第幾塊大餅
        #                         # print(temp)
        #                     temp_list.append(temp)
        #                 # print(temp_list)
        #                 endt = time.time()
        #                 # print("第%s條染色體，第幾%s個基因解碼判斷速度：%s" % (i,Order_index,endt -st))
        #                 Storage_Opt_Machine.append(temp_Storage_Opt_Machine)
        #                 Slice_of_pie.append(temp_list)
        #
        #             #新增表格
        #             sttake = time.time()
        #             temp_Select_Fitness_Df1 = np.zeros((2*Number_Of_Chromosome,Number_Of_Obj+1))  #+1是因為要包含Index
        #
        #             for index,value in enumerate(Storage_Opt_Machine):
        #                 # print(index,sum([len(x) for x in Storage_Opt_Machine[index].values()]))
        #                 temp_Select_Fitness_Df1[index] = index,sum([len(x) for x in Storage_Opt_Machine[index].values()]),0,0
        #
        #             # temp_Select_Fitness_Df1 = pd.DataFrame(temp_Select_Fitness_Df1,columns = ['Index',"機台數"])
        #             endtake = time.time()
        #             # print("新增表格",endtake-sttake,"秒")
        #             #%%選機，利用Slice_of_pie 之index-1 回傳Wip_canMachines 正確機台，Slice_of_pie[i]-1是因為index 從0開始 所以要-1
        #             Select_Machine = []
        #             st = time.time()
        #             for i in range(len(Slice_of_pie)):
        #                 #正確選機Index
        #                 # print(i)
        #                 temp_list = []
        #                 for Order_index in range(len(Slice_of_pie[i])):
        #                     if Order.at[Order_index, '設備/資源'] == Machine_model_1:
        #                     # if Order.iloc[Order_index]['設備/資源'] == Machine_model_1:
        #                         temp_list.append(Storage_Opt_Machine[i][Machine_model_1][Slice_of_pie[i][Order_index]-1])
        #                         # print(temp_list)
        #                     else:
        #                         temp_list.append(Storage_Opt_Machine[i][Machine_model_2][Slice_of_pie[i][Order_index]-1])
        #                         # print(temp_list)
        #                 Select_Machine.append(temp_list)
        #             # print(Select_Machine)
        #             endt = time.time()
        #             # print("第幾%s個基因選機解碼判斷速度：%s" % (Order_index,endt -st))
        #
        #         else:
        #
        #             '''other Interation'''
        #             Slice_of_pie = []
        #             Storage_Opt_Machine = []
        #             # temp_Storage_Opt_Machine = {} # 小小佔存器
        #             for i in range(Number_Of_Chromosome,len(chromosome)):
        #                 temp_list = []
        # # =============================================================================
        # #                 temp_Storage_Opt_Machine = {} # 小小佔存器
        # # =============================================================================
        #                 temp = 0
        #                 temp1 = 2
        #                 for index,value in enumerate(Mc_model_Quantity):
        #                     if value == Machine_model_1:
        #                         Available_Machine_model_Gene_St,Available_Machine_model_Gene_End = chromosome[i][Number_Of_Job*2+temp:Number_Of_Job*2+temp1]
        #
        #                         Num_Mc_Start = math.floor( Mc_model_Quantity[value]* Available_Machine_model_Gene_St)            #機台數 * 可用機台數 開始
        #                         Num_Mc_End   = int(np.round(Mc_model_Quantity[value] * Available_Machine_model_Gene_End ))  #機台數 * 可用機台數 結束
        #
        #                         diff_Num_Mc1 = Num_Mc_End - Num_Mc_Start #機台數結束 - 機台數開始 = 可用機台範圍
        #                         if diff_Num_Mc1 < 0:
        #                             '''交換 Start 跟 End'''
        #                             temp = Num_Mc_End
        #                             Num_Mc_End = Num_Mc_Start
        #                             Num_Mc_Start = temp
        #                             diff_Num_Mc1 = Num_Mc_End - Num_Mc_Start
        #                         elif diff_Num_Mc1 == 0:
        #                             '''機台開始與結束一樣，ex 開始58台 結束也58台，那就範圍 開始58 結束 59即可[58:59]'''
        #                             Num_Mc_End += 1
        #                         temp = temp1
        #                         temp1*=2
        #                     else:
        #                         Available_Machine_model_Gene_St,Available_Machine_model_Gene_End = chromosome[i][Number_Of_Job*2+temp:Number_Of_Job*2+temp1]
        #                         '''若某訂單的機型是Machine_model_2，就去增加他的機台型號數量'''
        #                         Num_Mc_Start = math.floor(Mc_model_Quantity[value] * Available_Machine_model_Gene_St)
        #                         Num_Mc_End   = int(np.round(Mc_model_Quantity[value] * Available_Machine_model_Gene_End))
        #                         diff_Num_Mc2 = Num_Mc_End - Num_Mc_Start
        #                         if diff_Num_Mc2 < 0:
        #                             '''交換 Start 跟 End'''
        #                             temp = Num_Mc_End
        #                             Num_Mc_End = Num_Mc_Start
        #                             Num_Mc_Start = temp
        #                             diff_Num_Mc2 = Num_Mc_End - Num_Mc_Start
        #                         elif diff_Num_Mc2 == 0:
        #                             '''機台開始與結束一樣，ex 開始58台 結束也58台，那就範圍 開始58 結束 59即可[58:59]'''
        #                             Num_Mc_End += 1
        #                 Total_Num_Mc = diff_Num_Mc1+diff_Num_Mc2
        #
        #                 while Total_Num_Mc < Mumber_of_machines_limit :
        #                     '''如果跑出來的機台數少於當時設定的Mumber_of_machines_limit  就要重新設定目前的染色體基因'''
        # # =============================================================================
        # #                     # for j in range(len(Mc_model_Quantity)):
        # #                     #     a = chromosome[i][Number_Of_Job*2:][j*2:(j+1)*2]
        # #                     #     # print('舊的基因',a)
        # #                     #     '''a[:]重新給予新的可用機台基因'''
        # #                     #     a[:] = np.array(np.random.rand(len(Machine_model_Number)))
        # # =============================================================================
        #                     '''a[:]重新給予新的可用機台基因'''
        #                     a = chromosome[i][Number_Of_Job*2:]
        #                     # print('舊的基因',a)
        #                     a[:] = Mc_Mutation(a)
        #
        #                     temp = 0
        #                     temp1 = 2
        #                     temp_Storage_Opt_Machine = {} # 小小佔存器
        #                     for index,value in enumerate(Mc_model_Quantity):
        #                         if value == Machine_model_1:
        #                             Available_Machine_model_Gene_St,Available_Machine_model_Gene_End = chromosome[i][Number_Of_Job*2+temp:Number_Of_Job*2+temp1]
        #
        #                             Num_Mc_Start = math.floor( Mc_model_Quantity[value]* Available_Machine_model_Gene_St)            #機台數 * 可用機台數 開始
        #                             Num_Mc_End   = int(np.round(Mc_model_Quantity[value] * Available_Machine_model_Gene_End ))  #機台數 * 可用機台數 結束
        #
        #                             diff_Num_Mc1 = Num_Mc_End - Num_Mc_Start #機台數結束 - 機台數開始 = 可用機台範圍
        #                             if diff_Num_Mc1 < 0:
        #                                 '''交換 Start 跟 End'''
        #                                 temp = Num_Mc_End
        #                                 Num_Mc_End = Num_Mc_Start
        #                                 Num_Mc_Start = temp
        #                                 diff_Num_Mc1 = Num_Mc_End - Num_Mc_Start
        #                             elif diff_Num_Mc1 == 0:
        #                                 '''機台開始與結束一樣，ex 開始58台 結束也58台，那就範圍 開始58 結束 59即可[58:59]'''
        #                                 Num_Mc_End += 1
        #                             temp = temp1
        #                             temp1*=2
        #                             # temp_Storage_Opt_Machine[Machine_model_1] = All_Machines_Name[Num_Mc_Start:Num_Mc_End]
        #                             temp_Storage_Opt_Machine[Machine_model_1] = All_Machines_Name[Num_Mc_Start:Num_Mc_End]
        #                         else:
        #                             Available_Machine_model_Gene_St,Available_Machine_model_Gene_End = chromosome[i][Number_Of_Job*2+temp:Number_Of_Job*2+temp1]
        #                             '''若某訂單的機型是Machine_model_2，就去增加他的機台型號數量'''
        #                             Num_Mc_Start = math.floor(Mc_model_Quantity[value] * Available_Machine_model_Gene_St)
        #                             Num_Mc_End   = int(np.round(Mc_model_Quantity[value] * Available_Machine_model_Gene_End))
        #                             diff_Num_Mc2 = Num_Mc_End - Num_Mc_Start
        #                             if diff_Num_Mc2 < 0:
        #                                 '''交換 Start 跟 End'''
        #                                 temp = Num_Mc_End
        #                                 Num_Mc_End = Num_Mc_Start
        #                                 Num_Mc_Start = temp
        #                                 diff_Num_Mc2 = Num_Mc_End - Num_Mc_Start
        #                             elif diff_Num_Mc2 == 0:
        #                                 '''機台開始與結束一樣，ex 開始58台 結束也58台，那就範圍 開始58 結束 59即可[58:59]'''
        #                                 Num_Mc_End += 1
        #
        #                             temp_Storage_Opt_Machine[Machine_model_2] = All_Machines_Name[Mc_model_Quantity[Machine_model_1] + Num_Mc_Start -1 : Mc_model_Quantity[Machine_model_1] + Num_Mc_End]
        #
        #                     Total_Num_Mc = diff_Num_Mc1+diff_Num_Mc2   #重新計算總共多少台機台
        # # =============================================================================
        #
        # # =============================================================================
        #                 st = time.time()
        #                 length_temp_strge_opt_machine_model1 = (1/len(temp_Storage_Opt_Machine[Machine_model_1]))
        #                 length_temp_strge_opt_machine_model2 = (1/len(temp_Storage_Opt_Machine[Machine_model_2]))
        #                 for Order_index in range(len(chromosome[i][:Number_Of_Job])):
        #                     '''判斷是哪個機型!!(切忌少用iloc，at 速度更快!!)'''
        #                     if Order.at[Order_index, '設備/資源'] == Machine_model_1:
        #                     # if Order.iloc[Order_index]['設備/資源'] == Machine_model_1:
        #                         temp  = math.ceil(chromosome[i][Order_index] / length_temp_strge_opt_machine_model1)
        #                         # print(temp)
        #
        #                     else:
        #                         temp  = math.ceil(chromosome[i][Order_index] / length_temp_strge_opt_machine_model2) #主要計算選機部分式第幾塊大餅
        #                         # print(temp)
        #                     temp_list.append(temp)
        #                 # print(temp_list)
        #                 endt = time.time()
        #                 # print("第%s條染色體，第幾%s個基因解碼判斷速度：%s" % (i,Order_index,endt -st))
        #                 Storage_Opt_Machine.append(temp_Storage_Opt_Machine)
        #                 Slice_of_pie.append(temp_list)
        #
        #             #%%選機，利用Slice_of_pie 之index-1 回傳Wip_canMachines 正確機台，Slice_of_pie[i]-1是因為index 從0開始 所以要-1
        #             Select_Machine = []
        #             st = time.time()
        #             for i in range(Number_Of_Chromosome):
        #                 # print(i)
        #                 temp_list = []
        #                 for Order_index in range(len(Slice_of_pie[i])):
        #                     if Order.at[Order_index, '設備/資源'] == Machine_model_1:
        #                     # if Order.iloc[Order_index]['設備/資源'] == Machine_model_1:
        #                         temp_list.append(Storage_Opt_Machine[i][Machine_model_1][Slice_of_pie[i][Order_index]-1])
        #                         # print(temp_list)
        #                     else:
        #                         temp_list.append(Storage_Opt_Machine[i][Machine_model_2][Slice_of_pie[i][Order_index]-1])
        #                         # print(temp_list)
        #
        #                 Select_Machine.append(temp_list)
        #             # print(Select_Machine)
        #             endt = time.time()
        #             # print("第幾%s個基因選機解碼判斷速度：%s" % (Order_index,endt -st))
        #             '''新增表格'''
        #             sttake = time.time()
        #             temp_Select_Fitness_Df1 = np.zeros((Number_Of_Chromosome,Number_Of_Obj+1))  #+1是因為要包含Index
        #
        #             for index,value in enumerate(Storage_Opt_Machine):
        #                 '''只要解剩下Number_Of_Chromosome 次的 目標式解碼'''
        #                 # print(index,sum([len(x) for x in Storage_Opt_Machine[index].values()]))
        #                 temp_Select_Fitness_Df1[index] = index,sum([len(x) for x in Storage_Opt_Machine[index].values()]),0,0
        #
        #             # temp_Select_Fitness_Df1 = pd.DataFrame(temp_Select_Fitness_Df1,columns = ['Index',"機台數"])
        #             endtake = time.time()
        #             print("新增表格",endtake-sttake,"秒")

        # %%解碼
        #####---------正式解碼----------####

        '''選機，將機子分成幾等份 (選擇輪盤第幾塊)'''
        if k == 0:
            '''first Interation'''
            Slice_of_pie = []  # Slice_of_pie 儲存每個Job被分配到第幾塊大餅
            # st = time.time()
            for i in range(len(chromosome)):
                temp_list = []
                '''第幾塊餅，例如第二塊餅，但他的index是1'''
                st = time.time()
                for j in range(len(chromosome[i][:Number_Of_Job])):
                    temp = math.ceil(chromosome[i][j] / (1 / The_Number_Of_Mc[j]))  # 主要計算選機部分式第幾塊大餅
                    temp_list.append(temp)
                # print(temp_list)
                Slice_of_pie.append(temp_list)
                endt = time.time()
            # print("解碼判斷速度：",endt -st)
            # %%選機，利用Slice_of_pie 之index-1 回傳Wip_canMachines 正確機台，Slice_of_pie[i]-1是因為index 從0開始 所以要-1
            Select_Machine = []
            for i in range(len(Slice_of_pie)):
                # print(i)
                '''正確選機Index'''
                temp_list = []
                for j in range(len(Slice_of_pie[i])):
                    # print(Order.at[j,"orderBomKey"])
                    check_Machine = dict_canRunTool[Order.at[j, "orderBomKey"]][Slice_of_pie[i][j] - 1]
                    if check_Machine == "PM9VA001":
                        while check_Machine == "PM9VA001":
                            '''renew_random重新更新亂數，尋找非PM9VA001的機台'''
                            renew_random = np.random.random()
                            temp = math.ceil(renew_random / (1 / The_Number_Of_Mc[j]))
                            Slice_of_pie[i][j] = temp
                            check_Machine = dict_canRunTool[Order.at[j, "orderBomKey"]][Slice_of_pie[i][j] - 1]

                    temp_list.append(check_Machine)
                Select_Machine.append(temp_list)
            # print(Select_Machine)

        else:
            '''other Interation'''
            '''這裡要注意從第二代只解碼子代'''
            Slice_of_pie = []
            for i in range(Number_Of_Chromosome, len(chromosome)):
                temp_list = []

                for j in range(len(chromosome[i][:Number_Of_Job])):
                    temp = math.ceil(chromosome[i][j] / (1 / The_Number_Of_Mc[j]))  # 主要計算選機部分式第幾塊大餅
                    # print(temp)
                    temp_list.append(temp)

                Slice_of_pie.append(temp_list)

            # %%選機，利用Slice_of_pie 之index-1 回傳Wip_canMachines 正確機台，Slice_of_pie[i]-1是因為index 從0開始 所以要-1
            Select_Machine = []
            for i in range(Number_Of_Chromosome):
                # print(i)
                '''上面Slice_of_pie是呼叫用哪個index，而這裡是將index改成˙正確機台型號'''
                temp_list = []
                for j in range(len(Slice_of_pie[i])):
                    # print(j)
                    check_Machine = dict_canRunTool[Order.at[j, "orderBomKey"]][Slice_of_pie[i][j] - 1]
                    if check_Machine == "PM9VA001":
                        while check_Machine == "PM9VA001":
                            '''renew_random重新更新亂數，尋找非PM9VA001的機台'''
                            renew_random = np.random.random()
                            temp = math.ceil(renew_random / (1 / The_Number_Of_Mc[j]))
                            Slice_of_pie[i][j] = temp
                            check_Machine = dict_canRunTool[Order.at[j, "orderBomKey"]][Slice_of_pie[i][j] - 1]

                    temp_list.append(check_Machine)
                Select_Machine.append(temp_list)
            # print(Select_Machine)
        # %%創三維空間 創每部機台對應到的Job
        # Ex:最外層用Chromosome_Machine_Corresponds_To_Job
        #    中間層用len(chromosome)100條，然後最內層為 每條染色體50台機台

        if k == 0:
            '''必須修正606列-617計算時間(已修正新的)'''
            # ============================================================================
            '''加快速度605~639列'''
            st = time.time()
            Chromosome_Machine_Corresponds_To_Job_dict = []

            for Num_Chromosome in range(len(chromosome)):
                '''這邊可以檢查每條染色體都用都少台機台'''
                machines = {}
                for i in range(len(Select_Machine[Num_Chromosome])):
                    # print(not Select_Machine[Num_Chromosome][i] in machines)
                    if not Select_Machine[Num_Chromosome][i] in machines:
                        machines[Select_Machine[Num_Chromosome][i]] = []
                    machines[Select_Machine[Num_Chromosome][i]].append(i)
                    # print(machines[Select_Machine[Num_Chromosome][i]])
                    # print(machines)

                Chromosome_Machine_Corresponds_To_Job_dict.append(machines)
            endt = time.time()

            # print("Machine_Corresponds_To_Job 把Job 改成 後基因100個(未排序)",endt - st)

            st = time.time()
            '''為了將dict改成三維List'''
            Chromosome_Machine_Corresponds_To_Job = []
            for Num_Chromosome in range(len(chromosome)):
                # print(Num_Chromosome)
                Machine_Corresponds_To_Job = []
                for index, value in enumerate(All_Machines_Name):

                    '''(注意)下面 i 只是為了方便讓電腦看Job index從0~3715訂單'''
                    if value not in Chromosome_Machine_Corresponds_To_Job_dict[Num_Chromosome].keys():
                        '''若某個機台都沒有被分配到Job，需保留一個空的機台空間，否則後續會出現錯誤'''
                        Machine_Corresponds_To_Job.append([])
                    else:
                        Machine_Corresponds_To_Job.append(
                            Chromosome_Machine_Corresponds_To_Job_dict[Num_Chromosome][value])

                Chromosome_Machine_Corresponds_To_Job.append(Machine_Corresponds_To_Job)
            endt = time.time()
            # print("Machine_Corresponds_To_Job 把Job 改成 後基因100個(未排序)",endt - st)
            # %%計算非空機機台至Storage_Opt_Machine(04/26新增)
            Storage_Opt_Machine = Culate_Num_Mc_Actual_Open(Chromosome_Machine_Corresponds_To_Job)
            temp_Select_Fitness_Df1 = np.zeros(
                (2 * Number_Of_Chromosome, Number_Of_Obj + 1))  # +1是因為要包含Index，(0531修正) +1 是因為整合權重至單目標 0608 減掉目標 最後用合併

            for index, value in enumerate(Storage_Opt_Machine):
                # print(index,sum([len(x) for x in Storage_Opt_Machine[index].values()]))
                temp_Select_Fitness_Df1[index] = index, Storage_Opt_Machine[
                    index], 0, 0  # 0531多一個0代表要存放整合單一目標 0608減掉目標 最後用合併
            # %%要把 Machine_Corresponds_To_Job 把Job 改成 後基因100個(未排序)
            st = time.time()
            Machine_Corresponds_To_Job_part2 = []  # part2 是指染色體的第二部分 排序基因!!!
            # 會跑幾條染色體 EX 100次(100條)
            for Num_Chromosome in range(len(Chromosome_Machine_Corresponds_To_Job)):
                temp_part2 = []

                # EX下面跑50次(因50個機台)
                for Nb_Of_Machine_index in range(len(Chromosome_Machine_Corresponds_To_Job[Num_Chromosome])):
                    temp = []
                    # 每台機台內有幾個Job
                    for Len_Machine_index in range(
                            len(Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index])):
                        # print(Len_Machine_index)
                        ##-----------------不能動到Chromosome_Machine_Corresponds_To_Job，因為他要回去找染色體index---------##
                        '''下面分為傳統GA排序 及 Arrivaltime'''
                        temp.append(chromosome[Num_Chromosome][
                                        Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index][
                                            Len_Machine_index] + Number_Of_Job])
                        # =========================================part2改成ArrivalTime====================================
                        # Arrivaltime_stamp = Order.at[Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index][Len_Machine_index],"最早開始時間"]
                        # convert_Arrivaltemp = ArrivalTime(Arrivaltime_stamp)
                        '''[Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index][Len_Machine_index]+Number_Of_Job] 不用-1 因為Machine_Corresponds_To_Job沒有index +1'''
                        # temp.append(convert_Arrivaltemp)
                        # =============================================================================
                        # print(temp)
                    temp_part2.append(temp)
                Machine_Corresponds_To_Job_part2.append(temp_part2)
            endt = time.time()
            # print("已排序",endt - st)

        else:

            '''必須修正734列-742計算時間(已修正新的)'''
            '''因跑第二代染色體都只要算一半的數量就好'''
            st = time.time()
            Chromosome_Machine_Corresponds_To_Job_dict = []

            for Num_Chromosome in range(len(Select_Machine)):
                # print(Num_Chromosome)
                Machine_Corresponds_To_Job = []
                machines = {}
                for i in range(len(Select_Machine[Num_Chromosome])):
                    if not Select_Machine[Num_Chromosome][i] in machines:
                        machines[Select_Machine[Num_Chromosome][i]] = []
                    machines[Select_Machine[Num_Chromosome][i]].append(i)

                Chromosome_Machine_Corresponds_To_Job_dict.append(machines)
            endt = time.time()
            # print("Machine_Corresponds_To_Job 把Job 改成 後基因100個(未排序)",endt - st)

            # %%為了將dict改成三維List
            st = time.time()
            '''為了將dict改成三維List'''
            Chromosome_Machine_Corresponds_To_Job = []
            for Num_Chromosome in range(len(Select_Machine)):
                # print(Num_Chromosome)
                Machine_Corresponds_To_Job = []
                for index, value in enumerate(All_Machines_Name):

                    '''(注意)下面 i 只是為了方便讓電腦看Job index從0~3715訂單'''
                    if value not in Chromosome_Machine_Corresponds_To_Job_dict[Num_Chromosome].keys():
                        '''若某個機台都沒有被分配到Job，需保留一個空的機台空間，否則後續會出現錯誤'''
                        Machine_Corresponds_To_Job.append([])
                    else:
                        Machine_Corresponds_To_Job.append(
                            Chromosome_Machine_Corresponds_To_Job_dict[Num_Chromosome][value])

                Chromosome_Machine_Corresponds_To_Job.append(Machine_Corresponds_To_Job)
            endt = time.time()
            # print("Machine_Corresponds_To_Job 把Job 改成 後基因100個(未排序)",endt - st)

            # %%計算非空機機台至Storage_Opt_Machine(04/26新增)
            Storage_Opt_Machine = Culate_Num_Mc_Actual_Open(Chromosome_Machine_Corresponds_To_Job)
            temp_Select_Fitness_Df1 = np.zeros((Number_Of_Chromosome,
                                                Number_Of_Obj + 1))  # +1是因為要包含Index ，(0531修正) +1 是因為整合權重至單目標 #解碼只解一半 不需要*2   0608 減掉目標 最後用合併

            for index, value in enumerate(Storage_Opt_Machine):
                # print(index,sum([len(x) for x in Storage_Opt_Machine[index].values()]))
                temp_Select_Fitness_Df1[index] = index, Storage_Opt_Machine[
                    index], 0, 0  # 0531多一個0代表要存放整合單一目標 0608減掉目標 最後用合併
            # %%要把 Machine_Corresponds_To_Job 把Job 改成 後基因100個(未排序)
            Machine_Corresponds_To_Job_part2 = []  # part2 是指染色體的第二部分 排序基因!!!
            # 會跑幾條染色體 EX 100次(100條)
            for Num_Chromosome in range(len(Chromosome_Machine_Corresponds_To_Job)):
                # print(Num_Chromosome)
                temp_part2 = []

                # EX下面跑50次(因50個機台)
                for Nb_Of_Machine_index in range(len(Chromosome_Machine_Corresponds_To_Job[Num_Chromosome])):
                    temp = []
                    # 每台機台內有幾個Job
                    for Len_Machine_index in range(
                            len(Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index])):
                        '''第二次Interation之後都要從Num_Chromosome+Number_Of_Chromosome(子代開始計算)'''
                        temp.append(chromosome[Num_Chromosome + Number_Of_Chromosome][
                                        Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index][
                                            Len_Machine_index] + Number_Of_Job])
                    # =============================================================================
                    #                         Arrivaltime_stamp = Order.at[Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index][Len_Machine_index],"最早開始時間"]
                    #                         convert_Arrivaltemp = ArrivalTime(Arrivaltime_stamp)
                    #                         '''[Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index][Len_Machine_index]+Number_Of_Job] 不用-1 因為Machine_Corresponds_To_Job沒有index +1'''
                    #                         temp.append(convert_Arrivaltemp)
                    # =============================================================================
                    temp_part2.append(temp)
                Machine_Corresponds_To_Job_part2.append(temp_part2)
        # %% 另一種方法使用DataFrame Concat合併，要排序!! Machine_Corresponds_To_Job_sort 做排序
        # 順便做開始、過程、結束時間
        # ---*****可嘗試用成def*******---#
        Machine_Corresponds_To_Job_sorted = []
        Machine_Corresponds_To_Job_Pt = []
        Machine_Corresponds_To_Job_End = []
        Machine_Corresponds_To_Job_Start = []
        Machine_Corresponds_To_Job_Module = []
        OrderIndex_Converse_OrderNumber = []
        List_OrderNumber_Stime_Endtime = []
        st = time.time()
        # Totol_completion_Time = []
        if k == 0:
            # Maxmakespan = []
            Total_completion_Tardiness_Times = []  # 這是每個訂單的完工時間加總
        else:
            # Maxmakespan = new_Interation_Makespan
            # =============================================================================
            #             Total_completion_Tardiness_Times = new_Interation_Fitness #新一代這是每個訂單的完工時間加總
            # =============================================================================
            pass
        if k == 0:  # 解碼比較久需要用成表格
            temp_Machine_Corresponds_To_Job_sorted = [[] for i in range(len(chromosome))]

            for Num_Chromosome in range(len(chromosome)):  # 因為第一代母代及子代都未解碼，range是100條染色體都需解碼
                temp_Machine = [[] for i in range(Number_Of_Machine)]
                for Nb_Of_Machine_index in range(len(Chromosome_Machine_Corresponds_To_Job[Num_Chromosome])):
                    '''(注意)下面這個為了要用降冪所以先轉成array(將染色體後半段的基因排序拿來做排序)(速度處理得比較快)'''
                    array_Prob = np.array(Machine_Corresponds_To_Job_part2[Num_Chromosome][Nb_Of_Machine_index])
                    # =============================================================================
                    #                     ind = np.argsort(-array_Prob) #np.argsort(-array_Prob)可以回傳 array_Prob由大到小的Index值
                    # =============================================================================
                    ind = np.argsort(
                        array_Prob)  # np.argsort(array_Prob)可以回傳 array_Prob由小到大的Index值 (05/03 修改ArrivalTime排序)
                    temp_Machine[Nb_Of_Machine_index] = ind  # 把ind array得形式加入至 temp_Machine[Nb_Of_Machine_index]
                temp_Machine_Corresponds_To_Job_sorted[Num_Chromosome] = temp_Machine

                temp_chromosome = []
                for Nb_Of_Machine_index in range(len(Chromosome_Machine_Corresponds_To_Job[Num_Chromosome])):
                    temp_Machine = []
                    for Len_Machine_index in range(
                            len(temp_Machine_Corresponds_To_Job_sorted[Num_Chromosome][Nb_Of_Machine_index])):
                        temp_Machine.append(Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index][
                                                temp_Machine_Corresponds_To_Job_sorted[Num_Chromosome][
                                                    Nb_Of_Machine_index][Len_Machine_index]])
                        # print(temp_Machine)
                    temp_chromosome.append(temp_Machine)
                Machine_Corresponds_To_Job_sorted.append(temp_chromosome)  # 存100條排序過後的訂單index
                # 計算每個機台的Job之Arrival、ProcessTime時間(相對時間)、每部機台之Job產品料號
                temp_Mc_Corres_To_Job_ArrivalTime, temp_Mc_Corres_To_Job_Pt, temp_Mc_PPNumber = Machine_Corres_To_Job_ArrivalTime_Pt_PPN(
                    Machine_Corresponds_To_Job_sorted[Num_Chromosome])

                # 開始時間計算
                temp_Mc_Corres_To_Job_Start = []
                temp_Mc_Corres_To_Job_End = []
                st = time.time()
                Setup_count = 0  # 記錄換幾次模
                Setup_times = 0  # 記錄整體換模時間
                for Nb_Of_Machine_index in range(len(temp_Mc_Corres_To_Job_ArrivalTime)):
                    temp_Job_Start = []
                    temp_Job_End = []
                    for Nb_Of_Job_index in range(len(temp_Mc_Corres_To_Job_ArrivalTime[Nb_Of_Machine_index])):
                        # print(Nb_Of_Machine_index,Nb_Of_Job_index )
                        if Nb_Of_Job_index == 0:
                            '''每個機台的第一筆Job皆不用考慮換模時間'''
                            Job_Start = temp_Mc_Corres_To_Job_ArrivalTime[Nb_Of_Machine_index][Nb_Of_Job_index]
                            Job_End = Job_Start + temp_Mc_Corres_To_Job_Pt[Nb_Of_Machine_index][Nb_Of_Job_index]
                            temp_Job_Start.append(Job_Start)
                            temp_Job_End.append(Job_End)
                        else:
                            '''每個機台的第二筆Job後皆考慮換模時間，以下的換模時間有很多種方法，最好的方法是"用字典第一代"方式搜尋時間最快'''

                            '''用字典第一代'''
                            if temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index - 1] != \
                                    temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]:
                                Setup = Mould_dict[temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index - 1] +
                                                   temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]] * 60
                                Setup_count += 1
                                Setup_times += Setup
                            else:
                                Setup = 0  # 換模時間0

                            '''考慮模具 所以Job_Start要  Max(前一個Job結束時間+換模時間,此Job的ArrivalTime) 取大的當作目前Job的開始時間'''
                            Job_Start = max(temp_Job_End[Nb_Of_Job_index - 1],
                                            temp_Job_End[Nb_Of_Job_index - 1] + Setup,
                                            temp_Mc_Corres_To_Job_ArrivalTime[Nb_Of_Machine_index][
                                                Nb_Of_Job_index])  # 加模具

                            Job_End = Job_Start + temp_Mc_Corres_To_Job_Pt[Nb_Of_Machine_index][
                                Nb_Of_Job_index]  # Job結束時間 = Job開始 + ProcessTime
                            temp_Job_Start.append(Job_Start)
                            temp_Job_End.append(Job_End)
                    temp_Mc_Corres_To_Job_Start.append(temp_Job_Start)
                    temp_Mc_Corres_To_Job_End.append(temp_Job_End)

                '''(重要  要判斷交期  利用每個Job的完工時間與交期作依據 ，Completion_time 和 Due_Date_time!!!!!!!)'''
                Each_temp_Mc_Co_To_Job_End_matrix = List_to_Np(temp_Mc_Corres_To_Job_End)  # ˊ轉換成矩陣模式，方便做下一步np.sum的運算
                Each_temp_Mc_Co_To_Job_Arrival_matrix = List_to_Np(temp_Mc_Corres_To_Job_ArrivalTime)
                Each_Total_EndTime = np.sum(
                    Each_temp_Mc_Co_To_Job_End_matrix)  # 為了計算Total_completion_time，將每個Job的結束時間做加總
                Each_Total_ArrivalTime = np.sum(Each_temp_Mc_Co_To_Job_Arrival_matrix)
                Total_completion_time = Each_Total_EndTime - Each_Total_ArrivalTime
                '''處理目標，each_Total_completion_Tardiness_Time = Total_completion_time + Tardiness * w'''
                Sum_Tardiness = 0
                Sum_Tardiness_cost = 0

                st1 = time.time()
                for Nb_Of_Machine_index in range(len(temp_Mc_Corres_To_Job_ArrivalTime)):

                    for Nb_Of_Job_index in range(len(temp_Mc_Corres_To_Job_ArrivalTime[Nb_Of_Machine_index])):
                        Completion_time = temp_Mc_Corres_To_Job_End[Nb_Of_Machine_index][Nb_Of_Job_index]  # 每個Job的完工時間

                        Due_Date_time = Due_Date_times[
                                            Machine_Corresponds_To_Job_sorted[Num_Chromosome][Nb_Of_Machine_index][
                                                Nb_Of_Job_index]] - RecoverTime[
                                            Nb_Of_Machine_index]  # 每個Job的交期時間(-RecoverTime是為了轉成相對時間)
                        if Completion_time - Due_Date_time > 0:
                            '''若完成時間超過交期則會 >0，因此需加總在Sum_Tardiness'''
                            Sum_Tardiness += (Completion_time - Due_Date_time)
                            # '''計算此單延誤損失成本'''
                            # Sum_Tardiness_cost += Order.at[Machine_Corresponds_To_Job_sorted[Num_Chromosome][Nb_Of_Machine_index][Nb_Of_Job_index],"單張延誤成本"]
                endt1 = time.time()
                # print(f'測試1條原來速度{endt1-st1}')

                Sum_Tardiness = Sum_Tardiness / (60 * 60 * 24)  # 原本是秒數，已轉換為天數
                Setup_times = Setup_times / (60 * 60 * 24)  # 原本是秒數，將Setup_times轉為天數
                # f'第{k}條染色體，延誤總成本{Sum_Tardiness_cost}'#總延誤成本取代print

                # temp_Select_Fitness_Df1 = pd.DataFrame(temp_Select_Fitness_Df1,columns = ['Index',"機台數"])

                # ===================================0503 總延誤、總完工目標==========================================
                #                 temp_Select_Fitness_Df1[Num_Chromosome][2] , temp_Select_Fitness_Df1[Num_Chromosome][3]= Sum_Tardiness,Total_completion_time
                # =============================================================================
                '''總延誤時間、總換模時間'''
                temp_Select_Fitness_Df1[Num_Chromosome][2], temp_Select_Fitness_Df1[Num_Chromosome][
                    3] = Sum_Tardiness, Setup_times

                # print(Setup_times/(60*24))  #換成天數
                '''下面目標式用Setup_count'''
                # each_Total_completion_Tardiness_Time =  Total_completion_time * Penalty_value_weight_TCT + Sum_Tardiness * Penalty_value_weight + Setup_count * Penalty_value_weight_Setup
                # '''下面目標式用Setup_times'''
                # each_Total_completion_Tardiness_Time =  Total_completion_time * Penalty_value_weight_TCT + Sum_Tardiness * Penalty_value_weight + Setup_times * Penalty_value_weight_Setup

                '''下面目標僅有Tardiness'''
                # each_Total_completion_Tardiness_Time =   Total_completion_time * Penalty_value_weight_TCT

                # Total_completion_Tardiness_Times.append(each_Total_completion_Tardiness_Time)

                Machine_Corresponds_To_Job_Start.append(temp_Mc_Corres_To_Job_Start)
                Machine_Corresponds_To_Job_End.append(temp_Mc_Corres_To_Job_End)
                Machine_Corresponds_To_Job_Pt.append(temp_Mc_Corres_To_Job_Pt)

                # 將訂單編號、Job開始時間、結束時間綁在一起[(訂單Index,開始時間,結束時間)]
                # temp_OrderNumber_Stime_Endtime = OrderNumber_Stime_Endtime_Zip(Name = Machine_Corresponds_To_Job_sorted[Num_Chromosome], St = temp_Mc_Corres_To_Job_Start, End = temp_Mc_Corres_To_Job_End, eachMaxmakespan=Total_completion_Tardiness_Times[Num_Chromosome])
                temp_OrderNumber_Stime_Endtime = OrderNumber_Stime_Endtime_Zip(
                    Name=Machine_Corresponds_To_Job_sorted[Num_Chromosome], St=temp_Mc_Corres_To_Job_Start,
                    End=temp_Mc_Corres_To_Job_End, eachMaxmakespan=0)

                # 下面這個將原本Zip改成 List [[訂單Index,開始時間,結束時間]] 以方便畫圖時把ST、End時間做更改  因tuple 無法更改!!
                # temp_OrderNumber_Stime_Endtime_list = Zip_OrderNumber_Stime_Endtime_converList(temp_OrderNumber_Stime_Endtime, eachMaxmakespan = Total_completion_Tardiness_Times[Num_Chromosome])
                temp_OrderNumber_Stime_Endtime_list = Zip_OrderNumber_Stime_Endtime_converList(
                    temp_OrderNumber_Stime_Endtime, eachMaxmakespan=0)

                List_OrderNumber_Stime_Endtime.append(temp_OrderNumber_Stime_Endtime_list)

            # Each_Iteration_Objective1.append(min(Total_completion_Tardiness_Times)) #儲存每一代最好(小)的Objective
            '''第一代前插法'''
            #                for Improve_Interation in range(len(chromosome)):
            #                    st = time.time()
            #                    List_OrderNumber_Stime_Endtime[Improve_Interation], temp_Select_Fitness_Df1[Improve_Interation] = Forward_Insertion_Method_Improved(List_OrderNumber_Stime_Endtime[Improve_Interation],Improve_Interation,temp_Select_Fitness_Df1[Improve_Interation])
            #                    List_OrderNumber_Stime_Endtime[Improve_Interation][-1] = temp_Select_Fitness_Df1[Improve_Interation]    #只會放[Index, 機台數, 總延誤時間, 總完工時間]
            #                    endt = time.time()
            #                    print("第%s代，第%s條染色體，前插改善時間 %s"%(k,Improve_Interation,(endt-st)))

            '''Multiprocessing第一代前插法'''
            t = []
            processing_st = time.time()
            # print("original_D1", temp_Select_Fitness_Df1)
            m = multiprocessing.Manager()
            q1 = m.Queue()  # thread可放入process同樣的queue中
            q2 = m.Queue()
            for arg1 in range(nb_threads):
                '''thread有幾個'''
                for j in range(arg1, len(temp_Select_Fitness_Df1), nb_threads):
                    '''某個threading 要處理幾條染色體'''
                    t.append(multiprocessing.Process(target=Forward_Insertion_Method_Improved, args=(
                        q1, q2, arg1, List_OrderNumber_Stime_Endtime[j], temp_Select_Fitness_Df1[j])))
                # t.append(threading.Thread(target = Forward_Insertion_Method_Improved, args=(i, List_OrderNumber_Stime_Endtime,temp_Select_Fitness_Df1)))
            # print('>>>>>>>>>>>>>>>>>>>>>>>>> thread_create')

            for _t in t:
                _t.start()
            # pdb.set_trace()

            # print('<<<<<<<<<<<<<<<<<<<<<<<<< thread_join')
            for _t in t:
                _t.join()
            # ent = time.time()
            # print("前插法平行",ent-st)

            temp_Select_Fitness_Df1 = np.zeros_like(temp_Select_Fitness_Df1)  # 重新創建temp_Select_Fitness_Df1 空的array

            '''#重新創建List_OrderNumber_Stime_Endtime 空的，為了方便對應原來的Index (花超多時間解決這個問題)'''
            List_OrderNumber_Stime_Endtime = [[] for i in range(len(List_OrderNumber_Stime_Endtime))]

            # print("List length", len(List_OrderNumber_Stime_Endtime))
            for arg1 in range(nb_threads):
                '''thread有幾個'''
                for j in range(arg1, len(temp_Select_Fitness_Df1), nb_threads):
                    '''某個threading 要處理幾條染色體'''
                    small_temp = q1.get()  # temp_Select_Fitness_Df1 get 回來     (一條一條回來) get只能用一次要不然會不見

                    small_temp_List = q2.get()  # List_OrderNumber_Stime_Endtime get 回來(一條一條回來) get只能用一次要不然會不見

                    # print(small_temp_List[-1][0])  # 讀取List 最後index
                    '''下列就是為了找List KPI [index, 機台數 , 總延誤時間  , 總換模時間]， index 放到正確位子'''
                    List_OrderNumber_Stime_Endtime[int(small_temp_List[-1][0])].insert(int(small_temp_List[-1][0]),
                                                                                       small_temp_List)

                    temp_Select_Fitness_Df1[int(small_temp[0])] = small_temp
            processing_endt = time.time()
            print("平行前插法時間", processing_endt - processing_st)
            # breakpoint()
            # for i in range(len(List_OrderNumber_Stime_Endtime)):
            #     print(List_OrderNumber_Stime_Endtime[i][-1][-1])    #p List_OrderNumber_Stime_EndNumber_Stime_Endtime[i][-1][-1] or List_OrderNumber_Stime_Endtime[i][0][-1]  OK

        else:
            '''第二代之後要做的!! 將上一代選擇解碼過 Number_Of_Chromosome 條 選擇放進List_OrderNumber_Stime_Endtime'''
            for i in range(len(new_decode)):
                List_OrderNumber_Stime_Endtime.append(new_decode[i])

            for Num_Chromosome in range(Number_Of_Chromosome):  # 要改range改成50條染色體(為了加快速度解子代的染色體)
                temp_Machine = [[] for i in range(Number_Of_Machine)]
                for Nb_Of_Machine_index in range(len(Chromosome_Machine_Corresponds_To_Job[Num_Chromosome])):
                    '''(注意)下面這個為了要用降冪所以先轉成array(將染色體後半段的基因排序拿來做排序)(速度處理得比較快)'''
                    array_Prob = np.array(Machine_Corresponds_To_Job_part2[Num_Chromosome][Nb_Of_Machine_index])
                    # =============================================================================
                    #                     ind = np.argsort(-array_Prob) #np.argsort(-array_Prob)可以回傳 array_Prob由大到小的Index值
                    # =============================================================================
                    ind = np.argsort(
                        array_Prob)  # np.argsort(array_Prob)可以回傳 array_Prob由小到大的Index值 (05/03 修改ArrivalTime排序)
                    temp_Machine[Nb_Of_Machine_index] = ind  # 把ind array得形式加入至 temp_Machine[Nb_Of_Machine_index]
                temp_Machine_Corresponds_To_Job_sorted[Num_Chromosome] = temp_Machine

                temp_chromosome = []
                for Nb_Of_Machine_index in range(len(Chromosome_Machine_Corresponds_To_Job[Num_Chromosome])):
                    temp_Machine = []
                    for Len_Machine_index in range(
                            len(temp_Machine_Corresponds_To_Job_sorted[Num_Chromosome][Nb_Of_Machine_index])):
                        '''下面是將Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index] 每個機台內的Job按照temp_Machine_Corresponds_To_Job_sorted排序Index去做排序'''
                        temp_Machine.append(Chromosome_Machine_Corresponds_To_Job[Num_Chromosome][Nb_Of_Machine_index][
                                                temp_Machine_Corresponds_To_Job_sorted[Num_Chromosome][
                                                    Nb_Of_Machine_index][Len_Machine_index]])

                    temp_chromosome.append(temp_Machine)
                Machine_Corresponds_To_Job_sorted.append(temp_chromosome)  # 存100條排序過後的訂單index
                # 計算每個機台的Job之Arrival、ProcessTime時間(相對時間)、每部機台之Job產品料號
                '''Machine_Corres_To_Job_ArrivalTime_Pt_PPN   Function處理每個機台的Job之Arrival、PT時間(相對時間)'''
                temp_Mc_Corres_To_Job_ArrivalTime, temp_Mc_Corres_To_Job_Pt, temp_Mc_PPNumber = Machine_Corres_To_Job_ArrivalTime_Pt_PPN(
                    Machine_Corresponds_To_Job_sorted[Num_Chromosome])

                # 開始、結束時間計算
                temp_Mc_Corres_To_Job_Start = []
                temp_Mc_Corres_To_Job_End = []
                Setup_count = 0  # 記錄換幾次模
                Setup_times = 0  # 紀錄整體換模時間
                for Nb_Of_Machine_index in range(len(temp_Mc_Corres_To_Job_ArrivalTime)):
                    temp_Job_Start = []
                    temp_Job_End = []
                    for Nb_Of_Job_index in range(len(temp_Mc_Corres_To_Job_ArrivalTime[Nb_Of_Machine_index])):
                        # print(Nb_Of_Machine_index,Nb_Of_Job_index )
                        if Nb_Of_Job_index == 0:
                            '''每個機台的第一筆Job皆不用考慮換模時間'''
                            Job_Start = temp_Mc_Corres_To_Job_ArrivalTime[Nb_Of_Machine_index][Nb_Of_Job_index]
                            Job_End = Job_Start + temp_Mc_Corres_To_Job_Pt[Nb_Of_Machine_index][Nb_Of_Job_index]
                            temp_Job_Start.append(Job_Start)
                            temp_Job_End.append(Job_End)
                        else:
                            '''每個機台的第二筆Job後皆考慮換模時間，以下的換模時間有很多種方法，最好的方法是"用字典第一代"方式搜尋時間最快'''

                            '''要考慮模具'''

                            '''用字典第一代'''
                            if temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index - 1] != \
                                    temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]:
                                Setup = Mould_dict[temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index - 1] +
                                                   temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]] * 60
                                Setup_count += 1
                                Setup_times += Setup
                            else:
                                Setup = 0

                            '''考慮模具 所以Job_Start要  Max(前一個Job結束時間+換模時間,此Job的ArrivalTime) 取大的當作目前Job的開始時間'''
                            Job_Start = max(temp_Job_End[Nb_Of_Job_index - 1],
                                            temp_Job_End[Nb_Of_Job_index - 1] + Setup,
                                            temp_Mc_Corres_To_Job_ArrivalTime[Nb_Of_Machine_index][
                                                Nb_Of_Job_index])  # 加模具

                            Job_End = Job_Start + temp_Mc_Corres_To_Job_Pt[Nb_Of_Machine_index][
                                Nb_Of_Job_index]  # Job結束時間 = Job開始 + ProcessTime
                            temp_Job_Start.append(Job_Start)
                            temp_Job_End.append(Job_End)
                    temp_Mc_Corres_To_Job_Start.append(temp_Job_Start)
                    temp_Mc_Corres_To_Job_End.append(temp_Job_End)

                '''(重要  要判斷交期  利用每個Job的完工時間與交期作依據 ，Completion_time 和 Due_Date_time!!!!!!!)'''
                # =============================================================================
                #                 Each_temp_Mc_Co_To_Job_End_matrix = List_to_Np(temp_Mc_Corres_To_Job_End) #ˊ轉換成矩陣模式，方便做下一步np.sum的運算
                #                 Total_completion_time = np.sum(Each_temp_Mc_Co_To_Job_End_matrix) #為了計算Total_completion_time，將每個Job的結束時間做加總
                # =============================================================================
                Each_temp_Mc_Co_To_Job_End_matrix = List_to_Np(temp_Mc_Corres_To_Job_End)  # ˊ轉換成矩陣模式，方便做下一步np.sum的運算
                Each_temp_Mc_Co_To_Job_Arrival_matrix = List_to_Np(temp_Mc_Corres_To_Job_ArrivalTime)
                Each_Total_EndTime = np.sum(
                    Each_temp_Mc_Co_To_Job_End_matrix)  # 為了計算Total_completion_time，將每個Job的結束時間做加總
                Each_Total_ArrivalTime = np.sum(Each_temp_Mc_Co_To_Job_Arrival_matrix)
                Total_completion_time = Each_Total_EndTime - Each_Total_ArrivalTime

                '''處理目標，each_Total_completion_Tardiness_Time = Total_completion_time + Tardiness * w'''
                Sum_Tardiness = 0
                Sum_Tardiness_cost = 0
                for Nb_Of_Machine_index in range(len(temp_Mc_Corres_To_Job_ArrivalTime)):

                    for Nb_Of_Job_index in range(len(temp_Mc_Corres_To_Job_ArrivalTime[Nb_Of_Machine_index])):
                        Completion_time = temp_Mc_Corres_To_Job_End[Nb_Of_Machine_index][Nb_Of_Job_index]

                        Due_Date_time = Due_Date_times[
                                            Machine_Corresponds_To_Job_sorted[Num_Chromosome][Nb_Of_Machine_index][
                                                Nb_Of_Job_index]] - RecoverTime[
                                            Nb_Of_Machine_index]  # 每個Job的交期時間(-RecoverTime是為了轉成相對時間)
                        if Completion_time - Due_Date_time > 0:
                            '''若完成時間超過交期則會 >0，因此需加總在Sum_Tardiness'''
                            Sum_Tardiness += Completion_time - Due_Date_time
                            # '''計算此單延誤損失成本'''
                            # Sum_Tardiness_cost += Order.at[Machine_Corresponds_To_Job_sorted[Num_Chromosome][Nb_Of_Machine_index][Nb_Of_Job_index],"單張延誤成本"]

                Sum_Tardiness = Sum_Tardiness / (60 * 60 * 24)  # 將Sum_Tardiness轉換為天數
                Setup_times = Setup_times / (60 * 60 * 24)  # 將Setup_times轉為天數
                # print(Sum_Tardiness)

                # ===================================0503 總延誤、總完工目標==========================================
                #                 temp_Select_Fitness_Df1[Num_Chromosome][2] , temp_Select_Fitness_Df1[Num_Chromosome][3]= Sum_Tardiness,Total_completion_time
                # =============================================================================
                '''總延誤時間、總換模時間05/03修正'''
                temp_Select_Fitness_Df1[Num_Chromosome][2], temp_Select_Fitness_Df1[Num_Chromosome][
                    3] = Sum_Tardiness, Setup_times

                # print(Setup_times/(60*24))  #換成天數
                '''下面目標式用Setup_count'''
                # each_Total_completion_Tardiness_Time =  Total_completion_time * Penalty_value_weight_TCT + Sum_Tardiness * Penalty_value_weight + Setup_count * Penalty_value_weight_Setup
                '''下面目標式用Setup_times'''
                # each_Total_completion_Tardiness_Time =  Total_completion_time * Penalty_value_weight_TCT + Sum_Tardiness * Penalty_value_weight + Setup_times * Penalty_value_weight_Setup

                '''下面目標僅有Tardiness'''
                # each_Total_completion_Tardiness_Time =   Total_completion_time * Penalty_value_weight_TCT

                # 全部染色體的Total_completion_Tardiness_Times
                # Total_completion_Tardiness_Times.append(each_Total_completion_Tardiness_Time)

                Machine_Corresponds_To_Job_Start.append(temp_Mc_Corres_To_Job_Start)
                Machine_Corresponds_To_Job_End.append(temp_Mc_Corres_To_Job_End)
                Machine_Corresponds_To_Job_Pt.append(temp_Mc_Corres_To_Job_Pt)

                # 將訂單編號、Job開始時間、結束時間綁在一起[(訂單Index,開始時間,結束時間)]
                temp_OrderNumber_Stime_Endtime = OrderNumber_Stime_Endtime_Zip(
                    Name=Machine_Corresponds_To_Job_sorted[Num_Chromosome], St=temp_Mc_Corres_To_Job_Start,
                    End=temp_Mc_Corres_To_Job_End, eachMaxmakespan=temp_Select_Fitness_Df1[Num_Chromosome])

                # 下面這個將原本Zip改成 List [[訂單Index,開始時間,結束時間]]以方便畫圖時把ST、End時間做更改 因tuple 無法更改!!
                temp_OrderNumber_Stime_Endtime_list = Zip_OrderNumber_Stime_Endtime_converList(
                    temp_OrderNumber_Stime_Endtime, eachMaxmakespan=temp_Select_Fitness_Df1[Num_Chromosome])

                List_OrderNumber_Stime_Endtime.append(temp_OrderNumber_Stime_Endtime_list)

            # Each_Iteration_Objective.append(min(Total_completion_Tardiness_Times))
            #                for Improve_Interation in range(Number_Of_Chromosome):
            #                    st = time.time()
            #                    List_OrderNumber_Stime_Endtime[Number_Of_Chromosome + Improve_Interation], temp_Select_Fitness_Df1[Improve_Interation] = Forward_Insertion_Method_Improved(List_OrderNumber_Stime_Endtime[Number_Of_Chromosome + Improve_Interation],Improve_Interation,temp_Select_Fitness_Df1[Improve_Interation])
            #                    List_OrderNumber_Stime_Endtime[Number_Of_Chromosome + Improve_Interation][-1] = temp_Select_Fitness_Df1[Improve_Interation]
            #                    endt = time.time()
            '''Multiprocessing第二代前插法'''
            t = []
            processing_st = time.time()
            # print("original_D1", temp_Select_Fitness_Df1)
            m = multiprocessing.Manager()
            q1 = m.Queue()  # thread可放入process同樣的queue中
            q2 = m.Queue()
            Interation_later = []
            '''examine code1998 - 2052'''
            '''0621 OKOK'''
            for arg1 in range(nb_threads):
                '''thread有幾個'''
                for j in range(arg1, len(temp_Select_Fitness_Df1), nb_threads):
                    '''某個threading 要處理幾條染色體'''
                    # print("thread:%s, List:%s , temp_D1:%s" % (arg1,j+len(temp_Select_Fitness_Df1),j ))

                    t.append(multiprocessing.Process(target=Forward_Insertion_Method_Improved, args=(
                        q1, q2, arg1, List_OrderNumber_Stime_Endtime[j + len(temp_Select_Fitness_Df1)],
                        temp_Select_Fitness_Df1[j])))
                # t.append(threading.Thread(target = Forward_Insertion_Method_Improved, args=(i, List_OrderNumber_Stime_Endtime,temp_Select_Fitness_Df1)))

            # print('>>>>>>>>>>>>>>>>>>>>>>>>> thread_create')

            for _t in t:
                _t.start()
            # pdb.set_trace()

            # print('<<<<<<<<<<<<<<<<<<<<<<<<< thread_join')
            for _t in t:
                _t.join()
            # ent = time.time()
            # print("前插法平行",ent-st)

            temp_Select_Fitness_Df1 = np.zeros_like(temp_Select_Fitness_Df1)  # 重新創建temp_Select_Fitness_Df1 空的array

            '''#重新創建List_OrderNumber_Stime_Endtime 空的，為了方便對應原來的Index (花超多時間解決這個問題)'''
            temp_List_OrderNumber_Stime_Endtime = [[] for i in range(int(len(List_OrderNumber_Stime_Endtime) / 2))]

            for arg1 in range(nb_threads):
                '''thread有幾個'''
                for j in range(arg1, len(temp_Select_Fitness_Df1), nb_threads):
                    # print("thread:%s, List:%s , temp_D1:%s" % (arg1,j+len(temp_Select_Fitness_Df1),j ))
                    '''某個threading 要處理幾條染色體'''
                    small_temp = q1.get()  # temp_Select_Fitness_Df1 get 回來     (一條一條回來) get只能用一次要不然會不見

                    small_temp_List = q2.get()  # List_OrderNumber_Stime_Endtime get 回來(一條一條回來) get只能用一次要不然會不見
                    # print("thread:%s, List:%s , temp_D1:%s" % (arg1, j + len(temp_Select_Fitness_Df1), j))
                    # print('proess Name',small_temp_List[-1][0])  # 讀取List 最後index
                    # breakpoint()
                    '''下列就是為了找List KPI [index, 機台數 , 總延誤時間  , 總換模時間]， index 放到正確位子'''

                    temp_List_OrderNumber_Stime_Endtime[int(small_temp_List[-1][0])].insert(int(small_temp_List[-1][0]),
                                                                                            small_temp_List)

                    temp_Select_Fitness_Df1[int(small_temp[0])] = small_temp

                    List_OrderNumber_Stime_Endtime[int(small_temp[0]) + len(temp_Select_Fitness_Df1)] = \
                    temp_List_OrderNumber_Stime_Endtime[int(small_temp_List[-1][0])][0]

            processing_endt = time.time()
            print("平行前插法時間", processing_endt - processing_st)
        endt = time.time()
        print('計算開始時間結束時間', endt - st)
        # Fitness_result = DN.Normalization(temp_Select_Fitness_Df1[:,1:],Penalty_Weights).calc()     #計算fitness
        # temp_Select_Fitness_Df1  = np.concatenate((temp_Select_Fitness_Df1, Fitness_result), axis=1)  #做結合
        # if k == 0:
        #     '''要將Fitness放進去'''
        #     for Num_Chromosome in range(len(chromosome)):  #因為第一代母代及子代都未解碼，range是100條染色體都需解碼
        #         '''0608下面目標僅有Fitness'''
        #         List_OrderNumber_Stime_Endtime[Num_Chromosome][-1] = temp_Select_Fitness_Df1[Num_Chromosome]
        # else:
        #     for Num_Chromosome in range(int(len(chromosome) /2) ):  #因為第一代母代及子代都未解碼，range是100條染色體都需解碼
        #         '''0608下面目標僅有Fitness'''
        #         List_OrderNumber_Stime_Endtime[Number_Of_Chromosome + Num_Chromosome][-1] = temp_Select_Fitness_Df1[Num_Chromosome]

        # %%選擇NSGAII
        new_decode = []  # 要放新一代菁英跟 輪盤的解碼
        new_moudule = []  # 要放解碼後的模具
        # new_Interation_NSGA_Decoding =
        st = time.time()
        # 下面Fitness 就是TSTP !!!! 要把原來好的上一代挑出來的 Fitness 存起來
        # temp_Select_Fitness_Df = pd.DataFrame([ [index,value]  for index,value in enumerate(Total_completion_Tardiness_Times)],columns = ['Index',"Fitness"])
        '''下面新的'''
        # =============================================================================
        #         if k ==0:
        #             temp_Select_Fitness_Df1 = pd.DataFrame(temp_Select_Fitness_Df1,columns = ['Index',"機台數","總延誤時間","總完工時間"])
        #             Fix_temp_Select_Fitness_Df1 = copy.copy(temp_Select_Fitness_Df1)
        #         else:
        #             temp_Select_Fitness_Df1 = pd.DataFrame(temp_Select_Fitness_Df1,columns = ['Index',"機台數","總延誤時間","總完工時間"])
        #             Fix_temp_Select_Fitness_Df1 = Fix_temp_Select_Fitness_Df1.append(temp_Select_Fitness_Df1)
        #             Fix_temp_Select_Fitness_Df1.reset_index(inplace=True, drop=True)
        #             temp_Select_Fitness_Df1 = copy.copy(Fix_temp_Select_Fitness_Df1)
        # =============================================================================
        temp_Select_Fitness_Df1 = pd.DataFrame(temp_Select_Fitness_Df1, columns=['Index', "機台數", "總延誤時間", "總換模時間"])

        if k == 0:
            # Fix_temp_Select_Fitness_Df1 = copy.copy(temp_Select_Fitness_Df1)
            temp_Select_Fitness_Df1 = ValuesCal(temp_Select_Fitness_Df1.loc[:, ["機台數", "總延誤時間", "總換模時間"]],
                                                Number_Of_Chromosome)
            Fix_temp_Select_Fitness_Df1 = copy.copy(temp_Select_Fitness_Df1)
        else:
            '''母代(上一次選擇出來的)加上這次解碼過後的子代'''
            Fix_temp_Select_Fitness_Df1 = Fix_temp_Select_Fitness_Df1.append(temp_Select_Fitness_Df1)
            Fix_temp_Select_Fitness_Df1.reset_index(inplace=True, drop=True)

            temp_Select_Fitness_Df1 = copy.copy(Fix_temp_Select_Fitness_Df1)
            '''temp_Select_Fitness_Df1會有較新的front 跟 cw 出現(下面temp_Select_Fitness_Df1 跟Fix_temp_Select_Fitness_Df1 很像一樣動作 )'''
            temp_Select_Fitness_Df1 = ValuesCal(temp_Select_Fitness_Df1.loc[:, ["機台數", "總延誤時間", "總換模時間"]],
                                                Number_Of_Chromosome)
            Fix_temp_Select_Fitness_Df1 = ValuesCal(temp_Select_Fitness_Df1.loc[:, ["機台數", "總延誤時間", "總換模時間"]],
                                                    Number_Of_Chromosome)

        # test_NSGA_TEST = RankforNSGAII(temp_Select_Fitness_Df1)
        # =============================================================================
        #         if k ==0:
        #             temp_Select_Fitness_Df1.insert(0,"Index",[i for i in range(Number_Of_Chromosome*2)],True)
        #             Fix_temp_Select_Fitness_Df1.insert(0,"Index",[i for i in range(Number_Of_Chromosome*2)],True)
        #         else:
        #             temp_Select_Fitness_Df1.insert(0,"Index",[i for i in range(Number_Of_Chromosome*2)],True)
        #             Fix_temp_Select_Fitness_Df1.insert(0,"Index",[i for i in range(Number_Of_Chromosome*2)],True)
        # =============================================================================
        temp_Select_Fitness_Df1.insert(0, "Index", [i for i in range(Number_Of_Chromosome * 2)], True)
        Fix_temp_Select_Fitness_Df1.insert(0, "Index", [i for i in range(Number_Of_Chromosome * 2)], True)
        '''將test_NSGA_TEST  Front_value由小排到大'''
        # temp_Select_Fitness_Df1.sort_values("Front_value",ascending = True,inplace = True)
        N = 0
        Min_whole_Front = min(temp_Select_Fitness_Df1["Front_value"])
        Select_Fitness_Elitism_new = np.zeros(shape=(Number_Of_Chromosome, Number_Of_Job * 2))

        NSGAII_Select_population_Index = []  # 儲存前緣線及擁擠距離到下一代的Index
        while N + len(temp_Select_Fitness_Df1[temp_Select_Fitness_Df1[
                                                  "Front_value"] == Min_whole_Front]) <= Number_Of_Chromosome:  # 0524 將< 改成 <=
            '''把Front前面都加入至N 裡面若某個Front的個數超過Number_Of_Chromosome，就要用CW比較'''
            Select_population_Index = \
            temp_Select_Fitness_Df1[temp_Select_Fitness_Df1["Front_value"] == Min_whole_Front]["Index"].tolist()

            '''把菁英制挑選出來的前 int(Number_Of_Chromosome*0.2 條的染色體 存入temp_Select_Fitness_Elitism_new'''
            for i, v in enumerate(Select_population_Index):
                Select_Fitness_Elitism_new[N + i] = chromosome[v]
                NSGAII_Select_population_Index.append(v)

            # new_Interation_Fitness.append(temp_Select_Fitness_Df1[temp_Select_Fitness_Df1["Front_value"])

            N += len(temp_Select_Fitness_Df1[temp_Select_Fitness_Df1["Front_value"] == Min_whole_Front])
            Min_whole_Front += 1

            # print("測試")
        '''把前N筆資料拿掉(因為前N筆資料式菁英制的)，後100-N筆資料是輪盤法(CW法)'''

        # new_Interation_Fitness = temp_Select_Fitness_Df1.iloc[NSGAII_Select_population_Index]["總完工時間"].tolist()
        # 下面temp用意是 為了 drop 所找的index list
        temp = temp_Select_Fitness_Df1.iloc[NSGAII_Select_population_Index].index[:].tolist()
        # temp_Select_Fitness_Df1 = temp_Select_Fitness_Df1.drop(temp_Select_Fitness_Df1.iloc[temp_Select_Fitness_Df1["Front_value"]].index)
        temp_Select_Fitness_Df1 = temp_Select_Fitness_Df1.drop(temp, 0)

        # '''排出第幾名Rank，同一個Front_value，越大的的CW Rank越大; 反之越小的CW Rank越小'''
        sort_obj = ['Front_value', 'Crowding_Distance_Value'];  # 設定要取行名
        # tups = Total_value[cols].sort_values(cols,ascending = [True,False]).apply(tuple,1) #排名 (由大到小)
        temp_Select_Fitness_Df1 = temp_Select_Fitness_Df1.sort_values(sort_obj, ascending=[True, False])
        # temp_Select_Fitness_Df1 = RankforNSGAII(temp_Select_Fitness_Df1)

        for i in range(Number_Of_Chromosome - N):
            '''剩下的要跑CW法'''
            Min_whole_Front = min(temp_Select_Fitness_Df1["Front_value"])  # 更新Min_whole_Front
            cw = temp_Select_Fitness_Df1.iloc[0]["Crowding_Distance_Value"]
            same_Ft_cw = np.logical_and(temp_Select_Fitness_Df1["Front_value"] == Min_whole_Front,
                                        temp_Select_Fitness_Df1["Crowding_Distance_Value"] == cw)
            # if True  in same_Ft_cw:
            temp = temp_Select_Fitness_Df1[same_Ft_cw]["Index"].values[:].tolist()
            cw_choice = int(np.random.choice(range(len(temp)), size=1, replace=False))  # 0524要修改這串
            Select_Fitness_Elitism_new[i + N] = chromosome[temp[cw_choice]]
            NSGAII_Select_population_Index.append(temp[cw_choice])

            temp_Select_Fitness_Df1 = temp_Select_Fitness_Df1.drop(temp[cw_choice])

        # '''把前10筆資料拿掉(因為前10筆資料式菁英制的)，後N-10筆資料是輪盤法'''
        # temp_Select_Fitness_No_Elitism = temp_Select_Fitness_Df1.drop(temp_Select_Fitness_Df1.index[:int(Number_Of_Chromosome*0.2)])
        # # print(temp_Select_Fitness_No_Elitism)

        # '''排出第幾名Rank，同一個Front_value，越大的的CW Rank越大; 反之越小的CW Rank越小'''
        # temp_Select_Fitness_No_Elitism = RankforNSGAII(temp_Select_Fitness_No_Elitism)

        # temp_Select_Fitness_No_Elitism.sort_values("Rank",ascending = False,inplace = True)

        # '''把全部名次加總起來分之Rank，可以得到相對應的機率'''
        # temp_Select_Fitness_No_Elitism["ratio"]  =  temp_Select_Fitness_No_Elitism["Rank"] /np.sum(temp_Select_Fitness_No_Elitism["Rank"])
        # # print(temp_Select_Fitness_No_Elitism)

        # '''製作累加的Fitness'''
        # temp_Select_Fitness_No_Elitism["Accumulate"] = np.add.accumulate(temp_Select_Fitness_No_Elitism["ratio"])
        # # print(temp_Select_Fitness_No_Elitism)

        # '''隨機亂數製造出Number_Of_Chromosome*0.8條亂數，為了要丟出飛鏢'''
        # selection_rand=[np.random.rand() for i in range(int(Number_Of_Chromosome*0.8))]  #選出Number_Of_Chromosome*0.8 條 隨機機率 為了要用輪盤法挑選
        # # print(selection_rand)

        # '''此次出現的Select_population_Index 是為了存取菁英制前 [:int(Number_Of_Chromosome*0.2)] 條'''
        # Select_population_Index = list(map(int,temp_Select_Fitness_Df1.iloc[:int(Number_Of_Chromosome*0.2)]["Index"]))  #菁英制選了哪幾個，輪盤選了哪即條母體染色體
        # # print(Select_population_Index)

        # '''保留前10條染色體的Fitness(菁英制)，目前僅有Tardiness'''
        # new_Interation_Fitness = list(map(float,temp_Select_Fitness_Df1.iloc[:int(Number_Of_Chromosome*0.2)]["總延誤時間"]))  #菁英制選了哪幾個Fitness，輪盤選了哪即條母體染色體Fitness

        # '''創造輪盤法要存的染色體 空 array'''
        # temp_Select_Fitness_RWS_new = np.zeros(shape=(int(Number_Of_Chromosome*0.8),Number_Of_Job*2 ))

        # '''把菁英制挑選出來的前 int(Number_Of_Chromosome*0.2 條的染色體 存入temp_Select_Fitness_Elitism_new'''
        # temp_Select_Fitness_Elitism_new = chromosome[temp_Select_Fitness_Df1.iloc[:int(Number_Of_Chromosome*0.2)]["Index"]]

        # for i in range(int(Number_Of_Chromosome*0.8)):
        #     '''利用輪盤法要產生Number_Of_Chromosome*0.8條'''
        #     '''等級式輪盤法'''
        #     if selection_rand[i] <= temp_Select_Fitness_No_Elitism["Accumulate"].iloc[0]:

        #         temp_Select_Fitness_RWS_new[i] = chromosome[temp_Select_Fitness_No_Elitism["Accumulate"].index[0]]
        #         # print(temp_Select_Fitness_RWS_new )

        #         Select_population_Index.append(temp_Select_Fitness_No_Elitism["Accumulate"].index[0])
        #         new_Interation_Fitness.append(temp_Select_Fitness_No_Elitism["總延誤時間"].iloc[0])
        #         # print(new_Interation_Fitness)
        #     else:

        #         for j in range(0,len(temp_Select_Fitness_No_Elitism)):

        #             if selection_rand[i] > temp_Select_Fitness_No_Elitism["Accumulate"].iloc[j] and selection_rand[i]<= temp_Select_Fitness_No_Elitism["Accumulate"].iloc[j+1]:
        #                 # print(selection_rand[i] > temp_Select_Fitness_No_Elitism["Accumulate"].iloc[j] and selection_rand[i]<= temp_Select_Fitness_No_Elitism["Accumulate"].iloc[j+1])
        #                 temp_Select_Fitness_RWS_new[i] = chromosome[temp_Select_Fitness_No_Elitism["Accumulate"].index[j+1]]  #+1原因是因為要找 a < x <= b 要找到b 的值回傳index
        #                 # print(temp_Select_Fitness_RWS_new)
        #                 '''把輪盤法的Select_population_Index加進去'''
        #                 Select_population_Index.append(temp_Select_Fitness_No_Elitism["Accumulate"].index[j+1]) #把 菁英制挑出染色體的Index 跟 輪盤法挑出染色體 Index存下來
        #                 new_Interation_Fitness.append(temp_Select_Fitness_No_Elitism["總延誤時間"].iloc[j+1])
        #                 # print(new_Interation_Fitness)
        #                 break
        # %%
        # '''將TSTP由小排到大'''
        # temp_Select_Fitness_Df.sort_values("Fitness",ascending = True,inplace = True)
        # # print(temp_Select_Fitness_Df)

        # '''把前10筆資料拿掉(因為前10筆資料式菁英制的)，後N-10筆資料是輪盤法'''
        # temp_Select_Fitness_No_Elitism = temp_Select_Fitness_Df.drop(temp_Select_Fitness_Df.index[:int(Number_Of_Chromosome*0.2)])
        # # print(temp_Select_Fitness_No_Elitism)

        # '''排出第幾名Rank，越小的Fitness Rank越大; 反之越大的Fitness Rank越小'''
        # temp_Select_Fitness_No_Elitism["Rank"] = temp_Select_Fitness_No_Elitism["Fitness"].rank(ascending = False , method = 'min')
        # # print(temp_Select_Fitness_No_Elitism)

        # '''把全部名次加總起來分之Rank，可以得到相對應的機率'''
        # temp_Select_Fitness_No_Elitism["ratio"]  =  temp_Select_Fitness_No_Elitism["Rank"] /np.sum(temp_Select_Fitness_No_Elitism["Rank"])
        # # print(temp_Select_Fitness_No_Elitism)

        # '''製作累加的Fitness'''
        # temp_Select_Fitness_No_Elitism["Accumulate"] = np.add.accumulate(temp_Select_Fitness_No_Elitism["ratio"])
        # # print(temp_Select_Fitness_No_Elitism)

        # '''隨機亂數製造出Number_Of_Chromosome*0.8條亂數，為了要丟出飛鏢'''
        # selection_rand=[np.random.rand() for i in range(int(Number_Of_Chromosome*0.8))]  #選出Number_Of_Chromosome*0.8 條 隨機機率 為了要用輪盤法挑選
        # # print(selection_rand)

        # '''此次出現的Select_population_Index 是為了存取菁英制前 [:int(Number_Of_Chromosome*0.2)] 條'''
        # Select_population_Index = list(map(int,temp_Select_Fitness_Df.iloc[:int(Number_Of_Chromosome*0.2)]["Index"]))  #菁英制選了哪幾個，輪盤選了哪即條母體染色體
        # # print(Select_population_Index)

        # '''保留前10條染色體的Fitness(菁英制)'''
        # new_Interation_Fitness = list(map(float,temp_Select_Fitness_Df.iloc[:int(Number_Of_Chromosome*0.2)]["Fitness"]))  #菁英制選了哪幾個Fitness，輪盤選了哪即條母體染色體Fitness

        # '''創造輪盤法要存的染色體 空 array'''
        # temp_Select_Fitness_RWS_new = np.zeros(shape=(int(Number_Of_Chromosome*0.8),Number_Of_Job*2 ))

        # '''把菁英制挑選出來的前 int(Number_Of_Chromosome*0.2 條的染色體 存入temp_Select_Fitness_Elitism_new'''
        # temp_Select_Fitness_Elitism_new = chromosome[temp_Select_Fitness_Df.iloc[:int(Number_Of_Chromosome*0.2)]["Index"]]

        # for i in range(int(Number_Of_Chromosome*0.8)):
        #     '''利用輪盤法要產生Number_Of_Chromosome*0.8條'''
        #     '''等級式輪盤法'''
        #     if selection_rand[i] <= temp_Select_Fitness_No_Elitism["Accumulate"].iloc[0]:

        #         temp_Select_Fitness_RWS_new[i] = chromosome[temp_Select_Fitness_No_Elitism["Accumulate"].index[0]]
        #         # print(temp_Select_Fitness_RWS_new )

        #         Select_population_Index.append(temp_Select_Fitness_No_Elitism["Accumulate"].index[0])
        #         new_Interation_Fitness.append(temp_Select_Fitness_No_Elitism["Fitness"].iloc[0])
        #         # print(new_Interation_Fitness)
        #     else:

        #         for j in range(0,len(temp_Select_Fitness_No_Elitism)):

        #             if selection_rand[i] > temp_Select_Fitness_No_Elitism["Accumulate"].iloc[j] and selection_rand[i]<= temp_Select_Fitness_No_Elitism["Accumulate"].iloc[j+1]:
        #                 # print(selection_rand[i] > temp_Select_Fitness_No_Elitism["Accumulate"].iloc[j] and selection_rand[i]<= temp_Select_Fitness_No_Elitism["Accumulate"].iloc[j+1])
        #                 temp_Select_Fitness_RWS_new[i] = chromosome[temp_Select_Fitness_No_Elitism["Accumulate"].index[j+1]]  #+1原因是因為要找 a < x <= b 要找到b 的值回傳index
        #                 # print(temp_Select_Fitness_RWS_new)
        #                 '''把輪盤法的Select_population_Index加進去'''
        #                 Select_population_Index.append(temp_Select_Fitness_No_Elitism["Accumulate"].index[j+1]) #把 菁英制挑出染色體的Index 跟 輪盤法挑出染色體 Index存下來
        #                 new_Interation_Fitness.append(temp_Select_Fitness_No_Elitism["Fitness"].iloc[j+1])
        #                 # print(new_Interation_Fitness)
        #                 break

        # 創造新的染色體(將菁英制temp_Select_Fitness_Elitism_new 與 輪盤法 temp_Select_Fitness_RWS_new 合併再一起)
        # new_gen_arr = np.concatenate((temp_Select_Fitness_Elitism_new, temp_Select_Fitness_RWS_new))
        '''下面意旨  要前緣線最小又要總延誤時間最小，當作此次的最小折線圖/
        也就是蒐集目標式 numpy  ["Index","機台數","總延誤時間","總完工時間"]/
        到最後用 Index 當作最好的去做前插法!!'''

        '''#######################累加每代最好值進 Each_Iteration_ObjectiveI _使用總延誤時間##################################'''
        # temp_Each_Iteration_Objective = Fix_temp_Select_Fitness_Df1[np.logical_and(Fix_temp_Select_Fitness_Df1["Front_value"]==1,Fix_temp_Select_Fitness_Df1["總延誤時間"] ==min(Fix_temp_Select_Fitness_Df1["總延誤時間"]))]
        # if len(temp_Each_Iteration_Objective) !=1:
        #     '''(檢查)如果temp_Each_Iteration_Objective篩選出來不止一條最小的解(兩組一樣只是染色體不同條)，就隨機挑一條'''
        #     temp_rand = np.random.choice(range(len(temp_Each_Iteration_Objective)))  #EX：2條最好解隨機從裡面挑1條
        #     temp_Each_Iteration_Objective = temp_Each_Iteration_Objective.iloc[temp_rand].to_frame().T

        # Each_Iteration_Objective = Each_Iteration_Objective.append(temp_Each_Iteration_Objective)

        '''#######################累加每代最好值進 Each_Iteration_ObjectiveI _使用總延誤時間##################################'''
        TOPSIS_result = tp.Topsis(
            Fix_temp_Select_Fitness_Df1.iloc[np.where(Fix_temp_Select_Fitness_Df1["Front_value"] == 1)].iloc[:, 1:4],
            TOPSIS_Weights, criterias).calc()
        '''要將Fix表 跟 TOPSIS算好的Rank 做結合'''
        Fix_table = Fix_temp_Select_Fitness_Df1.iloc[np.where(Fix_temp_Select_Fitness_Df1["Front_value"] == 1)]
        Fix_table.reset_index(drop=True, inplace=True)
        TOPSIS_Good_table = pd.concat([Fix_table, TOPSIS_result], axis=1)  # 做結合

        if len(TOPSIS_Good_table[TOPSIS_Good_table['Rank'] == max(TOPSIS_Good_table['Rank'])]) != 1:
            '''不等於 1 代表 有多重解(想成周老師的 兩個0.9方案)，就要選最重要的目標 例如 總延誤時間'''
            temp_Each_Iteration_Objective = TOPSIS_Good_table.iloc[
                TOPSIS_Good_table[TOPSIS_Good_table['Rank'] == max(TOPSIS_Good_table['Rank'])][
                    '總延誤時間'].idxmin()].to_frame().T  # TOPSIS最大Rank，並f2 最後要改成 此問題最關鍵的目標  ex 總延誤時間最小
        else:
            temp_Each_Iteration_Objective = TOPSIS_Good_table[
                TOPSIS_Good_table['Rank'] == max(TOPSIS_Good_table['Rank'])]

        Each_Iteration_Objective = Each_Iteration_Objective.append(temp_Each_Iteration_Objective)

        '''#######################累加每代最好值進 Each_Iteration_Objective II _使用CW法最大##################################'''
        # temp_Each_Iteration_Objective = Fix_temp_Select_Fitness_Df1[Fix_temp_Select_Fitness_Df1["Front_value"]==0].sort_values(sort_obj,ascending = [True,False])
        # Each_Iteration_Objective = Each_Iteration_Objective.append(temp_Each_Iteration_Objective.iloc[0])
        # Each_Iteration_Objective.append(np.array(Fix_temp_Select_Fitness_Df1[np.logical_and(Fix_temp_Select_Fitness_Df1["Front_value"]==0,Fix_temp_Select_Fitness_Df1["總完工時間"] ==min(Fix_temp_Select_Fitness_Df1["總完工時間"]))]).tolist())
        new_gen_arr = Select_Fitness_Elitism_new
        for index, value in enumerate(NSGAII_Select_population_Index):  # 會儲存處理到下一次的
            # print(index,value)
            '''包含new_gen_arr已包括(菁英及輪盤挑選的染色體)'''
            new_decode.append(List_OrderNumber_Stime_Endtime[value])  # new_decode 儲存新的解碼

        '''下一代母體(50條)'''
        Fix_temp_Select_Fitness_Df1 = Fix_temp_Select_Fitness_Df1.iloc[
            NSGAII_Select_population_Index]  # 選50條到下一代選擇一開始的50條

        Start_Iterations += 1
        # print("第 %s 迭代 , 機台數為 %s 總延誤時間%s  總完工時間%s" %(Start_Iterations,Each_Iteration_Objective[-1][0][:-2][1],Each_Iteration_Objective[-1][0][:-2][2],Each_Iteration_Objective[-1][0][:-2][3]))
        print("第 %s 迭代 " % (Start_Iterations))
        print(
            f'機台數：{np.round(Each_Iteration_Objective["機台數"][-1:].values[0], 4)} \n總延誤時間：{np.round(Each_Iteration_Objective["總延誤時間"][-1:].values[0], 4)}\n總換模時間：{np.round(Each_Iteration_Objective["總換模時間"][-1:].values[0], 4)}')

        t1_e = time.time()
        t1 = (t1_e - t1_start)
        # count_TCT = Each_Iteration_Objective["總完工時間"].tolist()
        # =============================================================================
        #         count_TCT = Each_Iteration_Objective["Fitness"].tolist() #0531改成連續出現Fitness
        #         test_Maxcount = Counter(count_TCT)  #Counter 可以判斷出Fitness同樣出現幾次 處理成Counter({Fitness,出現同樣次數})
        #         test_Maxcount = test_Maxcount.most_common()  #若出現同時最多次會被記錄下來 處理成 [(Fitness,出現同樣次數)]
        #         test_Maxcount.sort(key=lambda tup: tup[0])  #要排序才能抓出最好的那次解看有沒有重複出現!!
        # =============================================================================
        Each_Iteration_Objective.reset_index(drop=True, inplace=True)
        test_Maxcount = Number_Of_Consecutive_Iterations_Occurs(Each_Iteration_Objective["總延誤時間"])  # 0601新增連續方法
        # print(test_Maxcount)
        endt = time.time()
        # print('選擇',endt - st)
        '''終止條件為 若 < set_Timer  以及 連續出現Fitness <= Continuous_termination 次數就必須終止'''
        if (t1 < set_Timer) and test_Maxcount <= Continuous_termination:
            print(t1)
        else:
            print(t1)
            break

        # if(t1 < set_Timer):
        #     print(t1)
        # else:
        #     break

    '''將Zip_OrderNumber_Stime_Endtime 將每條染色體的Fitness排序由小到大排'''
    '''取用Each_Iteration_Objective 最後一個[-1]的Index 去回推List_OrderNumber_Stime_Endtime[拿這個當作前插法]'''
    # Zip_OrderNumber_Stime_Endtime = sorted(List_OrderNumber_Stime_Endtime, key=lambda e:e[-1])
    Zip_OrderNumber_Stime_Endtime = List_OrderNumber_Stime_Endtime[int(Each_Iteration_Objective[-1:]["Index"])]
    '''收斂圖'''
    # convergence_Fitness(Each_Iteration_Objective[:k+1],k+1)
    # =============================================================================
    # # convergence_Fitness(count_TCT[:k+1],k+1)
    # =============================================================================
    '''新的多樣化收斂圖'''
    convergence_Fitness_subplot(Each_Iteration_Objective["機台數"], Each_Iteration_Objective["總延誤時間"],
                                Each_Iteration_Objective["總換模時間"], k + 1)
    '''deep_Zip_OrderNumber_Stime_Endtime 前插法前的紀錄(最好的那條解碼)'''
    deep_Zip_OrderNumber_Stime_Endtime = copy.deepcopy(Zip_OrderNumber_Stime_Endtime)
    '''---------------------------------------------------------------'''
    '''0526計算Zip_OrderNumber_Stime_Endtime 的換模開始及結束(為了畫圖需要)'''
    # 開始、結束時間計算

    Record_Job_Index = []
    for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime) - 1):
        temp_Job_Index = []
        temp_Setup_Start = []
        temp_Setup_End = []
        for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
            # print(Nb_Of_Machine_index,Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0])
            temp_Job_Index.append(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0])
        Record_Job_Index.append(temp_Job_Index)

    Record_Setup_St_Endt_Index = []
    Record_Job_PPN = Setup_Machine_Corres_To_Job_PPN(Record_Job_Index)  # Zip處理成Record_Job_PPN

    temp_Record_Setup_St_Endt_Index = []
    temp_Mc_Corres_To_Setup_Index = []
    temp_Mc_Corres_To_Setup_Start = []
    temp_Mc_Corres_To_Setup_End = []
    Setup_count = 0  # 記錄換幾次模
    Setup_times = 0  # 紀錄整體換模時間
    for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime) - 1):
        temp_Setup_Index = []
        temp_Setup_Start = []
        temp_Setup_End = []
        for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
            # print(Nb_Of_Machine_index,Nb_Of_Job_index )
            if Nb_Of_Job_index == 0:
                '''每個機台的第一筆Job皆不用考慮換模時間'''
                Setup_Start = Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2]
                Setup_End = Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2]
                temp_Setup_Index.append("Setup")
                temp_Setup_Start.append(Setup_Start)
                temp_Setup_End.append(Setup_End)
            else:
                '''每個機台的第二筆Job後皆考慮換模時間，以下的換模時間有很多種方法，最好的方法是"用字典第一代"方式搜尋時間最快'''

                '''要考慮模具'''

                '''用字典第一代'''
                if Record_Job_PPN[Nb_Of_Machine_index][Nb_Of_Job_index - 1] != Record_Job_PPN[Nb_Of_Machine_index][
                    Nb_Of_Job_index]:
                    Setup = Mould_dict[Record_Job_PPN[Nb_Of_Machine_index][Nb_Of_Job_index - 1] +
                                       Record_Job_PPN[Nb_Of_Machine_index][Nb_Of_Job_index]] * 60
                    Setup_count += 1
                    Setup_times += Setup
                else:
                    Setup = 0

                '''考慮模具 所以Job_Start要  Max(前一個Job結束時間+換模時間,此Job的ArrivalTime) 取大的當作目前Job的開始時間'''
                Setup_End = Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][
                    1]  # Setup結束時間 = 目前Job開始時間
                Setup_Start = Setup_End - Setup  # 加模具

                temp_Setup_Index.append("Setup")
                temp_Setup_Start.append(Setup_Start)
                temp_Setup_End.append(Setup_End)

        temp_Mc_Corres_To_Setup_Index.append(temp_Setup_Index)
        temp_Mc_Corres_To_Setup_Start.append(temp_Setup_Start)
        temp_Mc_Corres_To_Setup_End.append(temp_Setup_End)
        # Record_Setup_St_Endt_Index.append(temp_Job_Start)

    # 將訂單編號、Job開始時間、結束時間綁在一起[(訂單Index,開始時間,結束時間)]
    Record_Setup_St_Endt_Index = OrderNumber_Stime_Endtime_Zip(Name=temp_Mc_Corres_To_Setup_Index,
                                                               St=temp_Mc_Corres_To_Setup_Start,
                                                               End=temp_Mc_Corres_To_Setup_End, eachMaxmakespan=0)

    # 下面這個將原本Zip改成 List [[訂單Index,開始時間,結束時間]]以方便畫圖時把ST、End時間做更改 因tuple 無法更改!!
    Record_Setup_St_Endt_Index = Zip_OrderNumber_Stime_Endtime_converList(Record_Setup_St_Endt_Index)

    '''--------------------------------------------------------------'''


    def TakeSt_minute_Func(Nice_Chromosome):
        # ---將 每個Job對應到的開始時間讀取出來到下一個stTime_out_Func---#
        temp_List = []  # 紀錄有50台機台
        for Num_Mc in range(Number_Of_Machine):

            temp = []  # 用temp主要是因為有50台機台想要紀錄每台機台內的Job數
            for Num_Job in range(len(Nice_Chromosome[Num_Mc])):
                # print('第 %s Job :開始時間 %s' % (Num_Job,Nice_Chromosome[Num_Mc][Num_Job][1]))
                TakeSt_minute = converge_Time(stsec=Nice_Chromosome[Num_Mc][Num_Job][1] + RecoverTime[Num_Mc])
                # print(TakeSt_minute)    #換算好的時間戳%Y-%m-%d %H:%M:%S
                temp.append(TakeSt_minute)
            temp_List.append(temp)

        return temp_List  # 回傳至draw_machine函數裡的st


    def converge_Time(stsec):
        # ---將converge_Time 把時間戳，轉回"%Y-%m-%d %H:%M:%S"---#
        Time_stamp = stsec  # 設定timeStamp
        struct_time = time.localtime(Time_stamp)  # 轉成時間元組
        timeString = time.strftime("%Y-%m-%d %H:%M:%S", struct_time)  # 轉成字串
        # print(timeString)   #將時間戳轉換成%Y-%m-%d %H:%M:%S
        return timeString


    # ------------------------------------下面是結束時間的換算----------------------------------------------#

    def TakeEnd_minute_Func(Nice_Chromosome):
        # ---將 每個Job對應到的結束時間讀取出來到下一個EndTime_out_Func---#
        temp_List = []  # 紀錄有50台機台
        for Num_Mc in range(Number_Of_Machine):

            temp = []  # 用temp主要是因為有50台機台想要紀錄每台機台內的Job數
            for Num_Job in range(len(Nice_Chromosome[Num_Mc])):
                # print('第 %s Job :結束時間 %s' % (Num_Job,Nice_Chromosome[Num_Mc][Num_Job][2]))
                TakeEnd_minute = converge_Time_Func(Endsec=Nice_Chromosome[Num_Mc][Num_Job][2] + RecoverTime[Num_Mc])
                # print(TakeEnd_minute)  #換算好的時間戳%Y-%m-%d %H:%M:%S
                temp.append(TakeEnd_minute)
            temp_List.append(temp)

        return temp_List  # 回傳至draw_machine函數裡的endt


    def converge_Time_Func(Endsec):
        # ---將converge_Time 把時間戳，轉回"%Y-%m-%d %H:%M:%S"---#
        Time_stamp = Endsec  # 設定timeStamp
        struct_time = time.localtime(Time_stamp)  # 轉成時間元組
        timeString = time.strftime("%Y-%m-%d %H:%M:%S", struct_time)  # 轉成字串
        # print(timeString)  #將時間戳轉換成%Y-%m-%d %H:%M:%S
        return timeString


    # =============================================================================

    df = []


    def draw_machine(index, value):
        '''把Zip_OrderNumber_Stime_Endtime排序過後的第一條拿來畫圖!!'''

        endt = TakeEnd_minute_Func(Zip_OrderNumber_Stime_Endtime)
        # print(endt)
        st = TakeSt_minute_Func(Zip_OrderNumber_Stime_Endtime)
        # print(st)
        each_Job = [j for j in Zip_OrderNumber_Stime_Endtime[index]]

        for i in range(len(each_Job)):
            if each_Job[i][0] == 'setup':
                df.append(dict(Task='%s' % 'setup', Start='%s' % st[index][i],
                               Finish='%s' % endt[index][i], Resource=value, item="換模", Mould="換模"))
            # elif each_Job[i][0] == 'Stop':
            #     df.append(dict(Task='%s'% 'Stop', Start='%s' % st[index][i],
            #     Finish='%s'  % endt[index][i],Resource= value,item = "停機",Mould = "停機"))
            else:
                df.append(dict(Task='%s' % OrderName[each_Job[i][0]], Start='%s' % st[index][i],
                               Finish='%s' % endt[index][i], Resource=value, item=Product_Part_Number[each_Job[i][0]]))

        return df


    for index, value in enumerate(All_Machines_Name):
        # print(index,value)
        draw_machine(index, value)

    df_setup = []


    def draw_Setup_machine(index, value):
        '''把Zip_OrderNumber_Stime_Endtime排序過後的第一條拿來畫圖!!'''

        endt = TakeEnd_minute_Func(Record_Setup_St_Endt_Index)
        # print(endt)
        st = TakeSt_minute_Func(Record_Setup_St_Endt_Index)
        # print(st)
        each_Job = [j for j in Record_Setup_St_Endt_Index[index]]

        for i in range(len(each_Job)):
            df_setup.append(dict(Task='Setup', Start='%s' % st[index][i],
                                 Finish='%s' % endt[index][i], Resource=value, item='換模'))

        return df_setup


    for index, value in enumerate(All_Machines_Name):
        # print(index,value)
        draw_Setup_machine(index, value)
    # ------------------------------------下面是開始時間的換算----------------------------------------------#

    # =============================================================================
    # #呈現圖表(0526前的程式碼)
    # fig = px.timeline(df + df_setup, x_start="Start", x_end="Finish",y="Resource", hover_name = "item" ,hover_data=["item"],color="item",color_discrete_map = colorss,text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    # # fig = px.timeline(df, x_start="Start", x_end="Finish",y="Resource",color = "Task", hover_name = "Task",text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    # fig.update_traces(textposition='inside',marker_line_color='rgb(8,48,107)')
    # # fig.update_layout(yaxis={'categoryorder':'category ascending'})  #將Y軸由小排到大
    # #將Y軸由大排到小
    #
    # fig.update_layout(yaxis={'categoryorder':'category descending'}, title={   # 设置整个标题的名称和位置
    #         "text":"甘特圖",
    #         "y":0.96,  # y轴数值
    #         "x":0.5,  # x轴数值
    #         "xanchor":"center",  # x、y轴相对位置
    #         "yanchor":"top"
    #     })
    # =============================================================================
    # 呈現圖表(0526後新的，按照All_Machines_Name排甘特圖)

    fig = px.timeline(df + df_setup, x_start="Start", x_end="Finish", y="Resource", hover_name="item",
                      hover_data=["item"], color="item", color_discrete_map=colorss, text="Task",
                      color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    # fig = px.timeline(df, x_start="Start", x_end="Finish",y="Resource",color = "Task", hover_name = "Task",text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    fig.update_traces(textposition='inside', marker_line_color='rgb(8,48,107)')
    # fig.update_layout(yaxis={'categoryorder':'category ascending'})  #將Y軸由小排到大
    # 將Y軸由大排到小
    fig.update_yaxes(categoryorder='array', categoryarray=All_Machines_Name)
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(title={  # 设置整个标题的名称和位置
        "text": "甘特圖",
        "y": 0.96,  # y轴数值
        "x": 0.5,  # x轴数值
        "xanchor": "center",  # x、y轴相对位置
        "yanchor": "top"
    })

    # fig.show()

    # plot(fig)

    '''下載出甘特圖及分解'''
    GanttChart = pd.DataFrame(df + df_setup)
    GanttChart = GanttChart.sort_values(by=['Resource'])
    GanttChart.reset_index(drop=True, inplace=True)
    GanttChart["Group"] = 0  # 創一欄 Group 欄位 裡面都放0
    All_Machine_uniqe = GanttChart["Resource"]
    All_Machine_uniqe = All_Machine_uniqe.drop_duplicates().tolist()
    # print(GanttChart)
    Machine_group_temp = 0
    for Machine_group in range(len(GanttChart)):
        '''可嘗試把 All_Machines_Name改成有畫圖的機台˙就好 可能不需要全跑'''
        if Machine_group == 0:
            GanttChart["Group"][Machine_group] = Machine_group_temp
        elif GanttChart["Resource"][Machine_group - 1] != GanttChart["Resource"][Machine_group]:
            Machine_group_temp += 1
            GanttChart["Group"][Machine_group] = Machine_group_temp
        else:
            GanttChart["Group"][Machine_group] = Machine_group_temp

        # GanttChart["Group"][GanttChart["Resource"] == All_Machines_Name[Machine_group]] = Machine_group
    Machinecount = 0;

    prev_temp = -1
    lower = -1
    upper = 0

    for i in range(int(np.ceil(len(GanttChart["Resource"].unique()) / plot_Machinenumber))):

        temp = plot_Machinenumber * (i + 1) - 1;
        # print(temp)
        if temp >= int(len(GanttChart["Group"].unique())):
            Machinecount = np.hstack((Machinecount, (len(GanttChart)) + 1))  # 0526有改len(部分)因為影響到plot_endpoint
            # print(Machinecount)
            break;
        temp1 = GanttChart[GanttChart.Group == temp].index.max()
        # print(temp1)
        Machinecount = np.hstack((Machinecount, temp1 + 1))
        # print(Machinecount)

    plot_startpoint = Machinecount[:-1]
    plot_endpoint = Machinecount[1:]

    for i in range(len(plot_startpoint)):
        # fig = px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i],],x_start= 'Start',x_end = 'Finish', y='Resource', hover_name = "item" ,hover_data=["item", "Mould"],color="item")
        # fig = px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i],], x_start="Start", x_end="Finish",y="Resource", hover_name = "item" ,hover_data=["item", "Mould"],color="item",color_discrete_map = colors,text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
        # 呈現圖表
        fig = px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i], ], x_start="Start", x_end="Finish",
                          y="Resource", hover_name="item", hover_data=["item"], color="item",
                          color_discrete_map=colorss, text="Task", color_discrete_sequence=px.colors.qualitative.Plotly,
                          title="甘特圖")
        # fig =px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i],], x_start="Start", x_end="Finish",y="Resource", hover_name = "item" ,hover_data=["item"],color="item",color_discrete_map = colorss,text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
        fig.update_traces(textposition='inside', marker_line_color='rgb(8,48,107)')

        # fig.update_layout(yaxis={'categoryorder':'category ascending'})  #將Y軸由小排到大
        # 將Y軸由大排到小
        # =============================================================================
        #
        #     fig.update_layout(yaxis={'categoryorder':'category descending'}, title={   # 设置整个标题的名称和位置
        #             "text":"甘特圖",
        #             "y":0.96,  # y轴数值
        #             "x":0.5,  # x轴数值
        #             "xanchor":"center",  # x、y轴相对位置
        #             "yanchor":"top"
        #         })
        #     py.offline.plot(fig, filename='C:/Users/ZhiXiang/Downloads/未前插'+ i.__str__()+'.html')    #下載html的路徑  可以自行更改路徑
        # =============================================================================
        fig.update_yaxes(categoryorder='array', categoryarray=All_Machine_uniqe)
        fig.update_yaxes(autorange="reversed")
        fig.update_layout(title={  # 设置整个标题的名称和位置
            "text": "甘特圖",
            "y": 0.96,  # y轴数值
            "x": 0.5,  # x轴数值
            "xanchor": "center",  # x、y轴相对位置
            "yanchor": "top"
        })
        # py.offline.plot(fig, filename='C:/Users/ZhiXiang/Downloads/法二傳統GA未前插'+ i.__str__()+'.html')    #下載html的路徑  可以自行更改路徑
    #
    # plot(fig)

    '''要做檢查!!'''
    '''將Zip_OrderNumber_Stime_Endtime 的 每個Job 訂單index ArrivalTime、PT、產品料號拆開，為了後續前插法'''
    best_temp_Mc_Corres_To_Job_Index, best_temp_Mc_Corres_To_Job_ArrivalTime, best_temp_Mc_Corres_To_Job_Pt, best_temp_Mc_PPNumber = Thebest_Machine_Corres_To_Job_ArrivalTime_Pt_PPN(
        Zip_OrderNumber_Stime_Endtime)

    '''0520計算前插法前的 Total_completion_time(Flow Time)'''
    Each_Total_EndTime = Thebest_List_to_Np(Zip_OrderNumber_Stime_Endtime)  # ˊ轉換成矩陣模式，方便做下一步np.sum的運算
    Each_temp_Mc_Co_To_Job_Arrival_matrix = List_to_Np(best_temp_Mc_Corres_To_Job_ArrivalTime)
    Each_Total_ArrivalTime = np.sum(Each_temp_Mc_Co_To_Job_Arrival_matrix)
    Total_completion_time = Each_Total_EndTime - Each_Total_ArrivalTime

    # =============================================================================
    # Total_completion_time = Thebest_List_to_Np(Zip_OrderNumber_Stime_Endtime)
    # =============================================================================
    '''處理目標，each_Total_completion_Tardiness_Time = Total_completion_time + Tardiness * w'''
    Sum_Tardiness = 0
    Setup_times = 0
    Setup_count = 0
    for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime) - 1):
        Machine_struct_time = time.strptime(Machine_Start.at[Nb_Of_Machine_index, "機台開始時間"], "%Y/%m/%d %H:%M")  # 轉成時間元組
        '''下面Machine_time_stamp是每台機台的起始時間，變動會隨著停機開始時間csv檔'''
        Machine_time_stamp = int(time.mktime(Machine_struct_time))

        for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
            Completion_time = Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2]

            Due_Date_time = int(time.mktime(time.strptime(Order.at[Order[Order.index == Zip_OrderNumber_Stime_Endtime[
                Nb_Of_Machine_index][Nb_Of_Job_index][0]].index[0], '交期'], "%Y/%m/%d %H:%M"))) - Machine_time_stamp

            if Completion_time - Due_Date_time > 0:
                Sum_Tardiness += (Completion_time - Due_Date_time)

            if Nb_Of_Job_index == 0:
                '''每個機台的第一筆Job皆不用考慮換模時間'''
                continue
            else:
                '''每個機台的第二筆Job後皆考慮換模時間，以下的換模時間有很多種方法，最好的方法是"用字典第一代"方式搜尋時間最快'''

                '''考慮模具'''
                if best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index - 1] != \
                        best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]:
                    Setup = Mould_dict[best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index - 1] +
                                       best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]] * 60
                    Setup_count += 1
                    Setup_times += Setup
                else:
                    Setup = 0

    Sum_Tardiness = Sum_Tardiness / (60 * 60 * 24)
    Setup_times = Setup_times / (60 * 60 * 24)
    print("=============================================================")
    print("未修正換模次數", Setup_count)
    print("未修正Obj", Zip_OrderNumber_Stime_Endtime[-1][1:])

    print("未修正Sum_Tardiness*W", Sum_Tardiness * Penalty_value_weight)
    # print("未修正Total_completion_time",Total_completion_time)

    print("未修正Setup_times*W", Setup_times * Penalty_value_weight_Setup)
    print("未修正Setup_count*W", Setup_count)
    print("總完工時間", Total_completion_time)
    '''計算前插完前閒置時間'''
    Sum_Idle_time = 0
    for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime) - 1):
        try:
            Machine_last_time = np.array(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][-1][-1])
        except:
            Machine_last_time = np.array(0)
        # print('機台最後時間',Machine_last_time)
        # print('每個機台總加工時間',np.sum(np.array(best_temp_Mc_Corres_To_Job_Pt[Nb_Of_Machine_index ])))
        Sum_Idle_time += (Machine_last_time - np.sum(np.array(best_temp_Mc_Corres_To_Job_Pt[Nb_Of_Machine_index]))) / (
                    60 * 60 * 24)
    each_Mc_Idle_time = Sum_Idle_time / Number_Of_Machine
    print('總閒置時間為:', Sum_Idle_time)

    print('平均一台閒置時間為:', each_Mc_Idle_time)
    Total_running_time_of_the_machine = 0
    for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime) - 1):
        try:
            Total_running_time_of_the_machine += np.array(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][-1][-1])
        except:
            Total_running_time_of_the_machine += np.array(0)
    Percentage_of_total_idle_time = Sum_Idle_time / (Total_running_time_of_the_machine / (60 * 60 * 24))
    print('總閒置時間比例=閒置時間(天)/總開機時間(天)', Percentage_of_total_idle_time)

    Over_Due_List = []
    count = 0
    for Num_Mc in range(len(Zip_OrderNumber_Stime_Endtime) - 1):
        temp = []

        for McToJobLengh in range(len(Zip_OrderNumber_Stime_Endtime[Num_Mc])):
            '''把交期設為相對時間'''
            str_time = Due_Date_times[Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][0]]
            Due_timeStamp = str_time - RecoverTime[Num_Mc]
            if Due_timeStamp < 0:
                Due_timeStamp = 0
            if Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][2] > Due_timeStamp:
                count += 1

            # print(Num_Mc,timeStamp)
            temp.append(Due_timeStamp)
        Over_Due_List.append(temp)
    print('在%s筆訂單內有 %s延誤' % (Number_Of_Job, count))
    print('延誤筆數/總訂單筆數 比率 %s' % (count / Number_Of_Job))
    # Old_OverDue_percent = count/Number_Of_Job

    # '''轉換成[(訂單Index,到達時間,加工時間)]'''
    # re_best =  Thebest_OrderNumber_St_Pt_Zip(Name =best_temp_Mc_Corres_To_Job_Index,ArrivalT = best_temp_Mc_Corres_To_Job_ArrivalTime,Pt = best_temp_Mc_Corres_To_Job_Pt)
    #
    # re_best_Zip = Thebest_Zip_OrderNumber_Stime_Endtime_converList(re_best)
    #
    # deep_Df_re_best_Zip = copy.deepcopy(re_best_Zip)
    #
    # for re_arrange_Mc in range(len(re_best)):
    #     '''N台機台'''
    #     # Timestamp_Limit = 1598889600 #2020/9/1
    #     Machine_struct_time  = time.strptime(Machine_Start.at[re_arrange_Mc,"機台開始時間"], "%Y/%m/%d %H:%M") # 轉成時間元組
    #     '''下面Machine_time_stamp是每台機台的起始時間，變動會隨著停機開始時間csv檔'''
    #     ''''目前都從0開始 只能先int(time.mktime(Machine_struct_time)) - int(time.mktime(Machine_struct_time))'''
    #     Machine_time_stamp_First   =  int(time.mktime(Machine_struct_time)) - int(time.mktime(Machine_struct_time))
    #
    #     '''從第二個間閣開始要尋找插入的ArrivalTime 的相對位置'''
    #     Machine_time_stamp_Second  =  int(time.mktime(Machine_struct_time))
    #
    #     for re_arrange_Job in range(len(re_best[re_arrange_Mc])+1):
    #         '''間格˙'''
    #         # print("第%s機台，整理 %s Job" %(re_arrange_Mc,re_arrange_Job))
    #         for inspection_Job in range(re_arrange_Job +1,len(re_best[re_arrange_Mc])):
    #             # print(re_arrange_Job,inspection_Job)
    #
    #             if re_arrange_Job == 0:
    #                 '''跑第0個間格'''
    #
    #                 insert_Job_Pt_day = re_best_Zip[re_arrange_Mc][inspection_Job][2]/(24*60*60)  #加工時間換算成天數
    #                 # print(insert_Job_Pt_day)
    #                 insert_Job_Pt_second = re_best_Zip[re_arrange_Mc][inspection_Job ][2]  #加工時間換算成秒數
    #                 # print(insert_Job_Pt_second)
    #
    #                 Setup_OrderName_PPN = Product_Part_Number[re_best_Zip[re_arrange_Mc][inspection_Job][0]] #目前產品料號
    #
    #                 Setup_RightOrderName_PPN = Product_Part_Number[Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]] #下一個的產品料號
    #
    #                 if Setup_OrderName_PPN != Setup_RightOrderName_PPN:
    #                     '''判斷目前產品料號與右邊(下一個)產品料號是否不同，若不同則需加上換模時間(秒數)'''
    #                     SetupRight_PPN_Time = Mould_dict[ Setup_OrderName_PPN + Setup_RightOrderName_PPN]*60
    #                 else:
    #                     SetupRight_PPN_Time = 0
    #
    #                 insert_JobPt_and_setup = (insert_Job_Pt_second  + SetupRight_PPN_Time)/(24*60*60)
    #                 insert_Job_Arrival = re_best_Zip[re_arrange_Mc][inspection_Job][1]
    #
    #                 '''最晚開始時間，若超過最晚開始一定會延遲(關鍵)'''
    #                 Latest_start_time = Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][1] - SetupRight_PPN_Time  - insert_Job_Pt_second
    #                 if Latest_start_time  < 0:
    #                     '''若最晚開始時間 < 0 代表可立即開始'''
    #                     Latest_start_time = 0
    #                 '''判斷get_between_days(Machine_time_stamp_First , Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][1]) 此間格天數 是否 大於 要插入的Job 加工天數 而且要符合 插入的 ArrivalTime 必須要 <= 最晚開始時間， 這樣才可做前插法'''
    #                 if get_between_days(Machine_time_stamp_First , Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][1])>= round(insert_JobPt_and_setup,3) and (insert_Job_Arrival <=Latest_start_time):
    #                     insert_Job_Start = max(insert_Job_Arrival,Machine_time_stamp_First )
    #                     insert_Job_End   = insert_Job_Start  +insert_Job_Pt_second
    #                     Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].pop(inspection_Job)
    #                     # Machine_Corresponds_To_Job_Module[0][re_arrange_Mc].pop(inspection_Job)
    #
    #                     Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].insert(re_arrange_Job ,[re_best_Zip[re_arrange_Mc][inspection_Job][0],insert_Job_Start,insert_Job_End])
    #
    #                     #下面五列程式碼是更正前插法的順序
    #                     temp = copy.deepcopy(re_best_Zip[re_arrange_Mc][inspection_Job])
    #                     re_best_Zip[re_arrange_Mc].pop(inspection_Job)
    #                     re_best_Zip[re_arrange_Mc].insert(re_arrange_Job,temp)
    #
    #                     best_temp_Mc_PPNumber[re_arrange_Mc].pop(inspection_Job)
    #                     best_temp_Mc_PPNumber[re_arrange_Mc].insert(re_arrange_Job,Order[Order.index == Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]]["產品料號"].values[0])
    #
    #                     break
    #                 else:
    #                     continue
    #             else:
    #                 '''跑第一個間格'''
    #                 insert_Job_Pt_day = re_best_Zip[re_arrange_Mc][inspection_Job][2]/(24*60*60)
    #                 # print(insert_Job_Pt_day)
    #                 insert_Job_Pt_second = re_best_Zip[re_arrange_Mc][inspection_Job ][2]
    #                 # print(insert_Job_Pt_second)
    #                 '''目前料號 與 上一個料號名稱'''
    #                 Setup_OrderName_PPN = Product_Part_Number[re_best_Zip[re_arrange_Mc][inspection_Job][0]]
    #                 Setup_LeftOrderName_PPN = Product_Part_Number[Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job - 1][0]]
    #
    #                 '''計算目前料號 與 上一個料號換模時間'''
    #                 if Setup_OrderName_PPN != Setup_LeftOrderName_PPN:
    #                     SetupLeft_PPN_Time = Mould_dict[ Setup_OrderName_PPN + Setup_LeftOrderName_PPN]*60
    #                 else:
    #                     SetupLeft_PPN_Time = 0
    #
    #                 '''目前料號 與 下一個料號名稱'''
    #                 Setup_RightOrderName_PPN = Product_Part_Number[Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]]
    #
    #                 '''計算目前料號 與 下一個料號換模時間'''
    #                 if Setup_OrderName_PPN != Setup_RightOrderName_PPN:
    #                     SetupRight_PPN_Time = Mould_dict[ Setup_OrderName_PPN + Setup_RightOrderName_PPN]*60
    #                 else:
    #                     SetupRight_PPN_Time = 0
    #
    #                 insert_JobPt_and_setup = (SetupLeft_PPN_Time + insert_Job_Pt_second  + SetupRight_PPN_Time)/(24*60*60)
    #
    #                 insert_Job_Arrival = Arrival_Time[re_best_Zip[re_arrange_Mc][inspection_Job][0]] -  Machine_time_stamp_Second
    #                 if insert_Job_Arrival   < 0:
    #                     insert_Job_Arrival  = 0
    #
    #                 Latest_start_time = Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][1] - SetupRight_PPN_Time  - insert_Job_Pt_second
    #                 if Latest_start_time  < 0:
    #                     Latest_start_time = 0
    #                 if get_between_days(Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job-1][2] ,Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][1])>= round(insert_JobPt_and_setup,3) and (insert_Job_Arrival <=Latest_start_time):
    #                     '''判斷要前插的空間是否夠插入，且插入的ArrvalTime必須<=最晚開始時間，否則無法插入'''
    #                     insert_Job_Start = max(insert_Job_Arrival,Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job - 1][2] + SetupLeft_PPN_Time)
    #                     insert_Job_End   = insert_Job_Start  +insert_Job_Pt_second
    #
    #                     Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].pop(inspection_Job)
    #
    #                     Zip_OrderNumber_Stime_Endtime[re_arrange_Mc].insert(re_arrange_Job ,[re_best_Zip[re_arrange_Mc][inspection_Job][0],insert_Job_Start,insert_Job_End])
    #
    #                     temp = copy.deepcopy(re_best_Zip[re_arrange_Mc][inspection_Job])
    #                     re_best_Zip[re_arrange_Mc].pop(inspection_Job)
    #                     re_best_Zip[re_arrange_Mc].insert(re_arrange_Job,temp)
    #
    #                     best_temp_Mc_PPNumber[re_arrange_Mc].pop(inspection_Job)
    #                     best_temp_Mc_PPNumber[re_arrange_Mc].insert(re_arrange_Job,Order[Order.index == Zip_OrderNumber_Stime_Endtime[re_arrange_Mc][re_arrange_Job][0]]["產品料號"].values[0])
    #
    #                     break
    #                 else:
    #                     continue
    #
    # '''要增加修改'''
    # '''修正開始時間結束時間'''
    #
    # Setup_count = 0
    # Setup_times = 0  #把這個取消掉 Sum_Tardiness 就跑出來  超怪的!!
    # for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime)-1):
    #     Deep_re_best_Zip_DF = pd.DataFrame(deep_Df_re_best_Zip[Nb_Of_Machine_index],columns =["訂單編號","最早到達時間","加工時間"])
    #     for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
    #         if Nb_Of_Job_index == 0:
    #             continue
    #         else:
    #
    #             if best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index-1] != best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]:
    #                 Setup = Mould_dict[best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index-1]+best_temp_Mc_PPNumber[Nb_Of_Machine_index][Nb_Of_Job_index]]*60
    #                 Setup_count += 1
    #                 Setup_times += Setup
    #             else:
    #                 Setup = 0
    #
    #             Pt = Process_Time[re_best_Zip[Nb_Of_Machine_index][Nb_Of_Job_index][0]]
    #             Job_ArrivalTime = Deep_re_best_Zip_DF[Deep_re_best_Zip_DF["訂單編號" ]== re_best_Zip[Nb_Of_Machine_index][Nb_Of_Job_index][0]]["最早到達時間"].values[0]
    #             Job_Start= max(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index-1][2]+Setup , Job_ArrivalTime)
    #             Job_End  = Job_Start  + Pt
    #             '''從第二個Job開始的開始時間'''
    #             Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][1] =  Job_Start
    #             Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2] = Job_End
    #             # print(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][1] ,Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2] )
    #
    #
    # '''0520前插法後的新的 Total_completion_time(Flow time) 前插法後的ArrivalTime需再確認'''
    # # =============================================================================
    # # Total_completion_time = Thebest_List_to_Np(Zip_OrderNumber_Stime_Endtime)
    # # =============================================================================
    # Each_Total_EndTime = Thebest_List_to_Np(Zip_OrderNumber_Stime_Endtime)  #ˊ轉換成矩陣模式，方便做下一步np.sum的運算
    # Each_temp_Mc_Co_To_Job_Arrival_matrix = List_to_Np(best_temp_Mc_Corres_To_Job_ArrivalTime)
    # Each_Total_ArrivalTime = np.sum(Each_temp_Mc_Co_To_Job_Arrival_matrix)
    # Total_completion_time =  Each_Total_EndTime - Each_Total_ArrivalTime
    #
    # '''處理目標，each_Total_completion_Tardiness_Time = Total_completion_time + Tardiness * w'''
    # Sum_Tardiness = 0
    # Sum_Tardiness_cost = 0
    # for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime)-1):
    #     Machine_struct_time  = time.strptime(Machine_Start.at[Nb_Of_Machine_index,"機台開始時間"], "%Y/%m/%d %H:%M") # 轉成時間元組
    #     '''下面Machine_time_stamp是每台機台的起始時間，變動會隨著停機開始時間csv檔'''
    #     Machine_time_stamp   =  int(time.mktime(Machine_struct_time))
    #
    #     for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
    #         Completion_time = Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2]
    #
    #         # Due_Date_time = int(time.mktime(time.strptime(Order.at[Order[Order["訂單編號"] ==  Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0]].index[0],'交期'], "%Y/%m/%d %H:%M")))-Machine_time_stamp
    #         Due_Date_time = int(time.mktime(time.strptime(Order.at[Order[Order.index ==  Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0]].index[0],'交期'], "%Y/%m/%d %H:%M")))-Machine_time_stamp
    #
    #         if Completion_time - Due_Date_time > 0:
    #             Sum_Tardiness += (Completion_time - Due_Date_time)
    #             '''計算此單延誤損失成本'''
    #             # Sum_Tardiness_cost += Order.at[Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0],"單張延誤成本"]
    #
    #
    # Sum_Tardiness = Sum_Tardiness/(60*60*24)
    # Setup_times = Setup_times/(60*60*24)
    # # =============================================================================
    # # # each_Total_completion_Tardiness_Time = Total_completion_time * Penalty_value_weight_TCT + Sum_Tardiness * Penalty_value_weight + Setup_times * Penalty_value_weight_Setup
    # # =============================================================================
    # # each_Total_completion_Tardiness_Time = Machine_Fixcost + Sum_Tardiness_cost
    #
    # '''要找時間修改'''
    # # =============================================================================
    # # Each_Iteration_Objective.append(each_Total_completion_Tardiness_Time )
    # # =============================================================================
    # Zip_OrderNumber_Stime_Endtime[-1] =  np.array([0,np.round(Each_Iteration_Objective["機台數"][-1:].values[0],4),Sum_Tardiness,Setup_times])
    #
    # '''0606改善後前插法排程指標及KPI'''
    # improve = copy.deepcopy(Each_Iteration_Objective.iloc[-1])
    # improve["總延誤時間"],improve["總換模時間"] = Zip_OrderNumber_Stime_Endtime[-1][2], Zip_OrderNumber_Stime_Endtime[-1][3]
    #
    # Each_Iteration_Objective = Each_Iteration_Objective.append(improve.to_frame().T)
    #
    # '''收斂圖'''
    #
    # # =============================================================================
    # # # convergence_Fitness(count_TCT[:k+1],k+1)
    # # =============================================================================
    # '''新的多樣化收斂圖'''
    # convergence_Fitness_subplot(Each_Iteration_Objective["機台數"],Each_Iteration_Objective["總延誤時間"],Each_Iteration_Objective["總換模時間"],k+1+1)
    #
    #
    #
    #
    # print("=============================================================")
    # print("最新的換模次數" ,Setup_count)
    # print("已修正Obj" ,Zip_OrderNumber_Stime_Endtime[-1][1:])
    #
    # print("已修正Sum_Tardiness*W" ,Sum_Tardiness * Penalty_value_weight)
    # # print("已修正Total_completion_time",Total_completion_time)
    # print("已修正Setup_times * Penalty_value_weight_Setup" ,Setup_times * Penalty_value_weight_Setup)
    # print("已修正Setup_count * Penalty_value_weight_Setup" ,Setup_count)
    # print("總完工時間",Total_completion_time)
    # # print(f"最新的換模次數：{Setup_count}")
    #
    # # print(f'已修正Total_completion_time：{Total_completion_time}')
    # # print(f'已修正Sum_Tardiness*W：{Sum_Tardiness * Penalty_value_weight}')
    # # # print("未修正Setup_count*W" ,Setup_count * Penalty_value_weight_Setup)
    # # print(f'已修正Setup_times * Penalty_value_weight_Setup：{Setup_times * Penalty_value_weight_Setup}')
    #
    #
    # '''計算前插完後閒置時間'''
    # Sum_Idle_time  = 0
    # for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime)-1):
    #     try:
    #         Machine_last_time = np.array(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][-1][-1])
    #     except:
    #         Machine_last_time = np.array(0)
    #     # print('機台最後時間',Machine_last_time)
    #     # print('每個機台總加工時間',np.sum(np.array(best_temp_Mc_Corres_To_Job_Pt[Nb_Of_Machine_index ])))
    #     Sum_Idle_time += (Machine_last_time-np.sum(np.array(best_temp_Mc_Corres_To_Job_Pt[Nb_Of_Machine_index ])))/(60*60*24)
    # each_Mc_Idle_time = Sum_Idle_time/Number_Of_Machine
    # print('總閒置時間為:',Sum_Idle_time)
    #
    # print('平均一台閒置時間為:',each_Mc_Idle_time)
    # Total_running_time_of_the_machine = 0
    # for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime)-1):
    #     try:
    #         Total_running_time_of_the_machine += np.array(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][-1][-1])
    #     except:
    #         Total_running_time_of_the_machine +=np.array(0)
    # Percentage_of_total_idle_time = Sum_Idle_time/(Total_running_time_of_the_machine/(60*60*24))
    # print('總閒置時間比例=閒置時間(天)/總開機時間(天)',Percentage_of_total_idle_time)
    #
    # Over_Due_List = []
    # count = 0   #累加延誤筆數
    # Add_up_Delayed_days = 0  #累加延誤天數
    # # Delayed_days_df = pd.DataFrame(columns= ['訂單名稱','延誤天數'])
    # '''在計算每筆延誤天數'''
    # Delayed_days_df = [] #為了能創造出二維List儲存每筆延誤天數，在後續轉成csv(新)
    # for Num_Mc in range(len(Zip_OrderNumber_Stime_Endtime)-1):
    #     temp = []
    #
    #     for McToJobLengh in range(len(Zip_OrderNumber_Stime_Endtime[Num_Mc])):
    #         '''把交期設為相對時間'''
    #         temp_Delayed_days_list = []  #每個訂單會存['訂單名稱','延誤天數']
    #         Each_Job_Delayed_days = 0  #個別訂單重新歸零計算是否超過交期
    #
    #         str_time = Due_Date_times[Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][0]]
    #         Due_timeStamp = str_time - RecoverTime[Num_Mc]
    #         if Due_timeStamp < 0 :
    #             Due_timeStamp = 0
    #         if Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][2] > Due_timeStamp:
    #             Add_up_Delayed_days += abs(Due_timeStamp - Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][2])/(60*60*24)  #總延誤天數
    #
    #             Each_Job_Delayed_days += abs(Due_timeStamp - Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][2])/(60*60*24)  #個別訂單延誤天數
    #
    #             temp_Delayed_days_list.append(Order[Order.index == Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][0]]['訂單編號'].values[0])
    #             temp_Delayed_days_list.append(Each_Job_Delayed_days)
    #             # print(Order[Order.index == Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][0]]['訂單編號'].values[0] +' 延誤 ' + str(Each_Job_Delayed_days)+'天' )
    #             Delayed_days_df.append(temp_Delayed_days_list) #加入到二維List裡面
    #             count +=1
    #         else:
    #             # print(Order[Order.index == Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][0]]['訂單編號'].values[0] +' 延誤 ' + str(Each_Job_Delayed_days)+'天' )
    #             temp_Delayed_days_list.append(Order[Order.index == Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][0]]['訂單編號'].values[0])
    #             temp_Delayed_days_list.append(Each_Job_Delayed_days)
    #             # print(Order[Order.index == Zip_OrderNumber_Stime_Endtime[Num_Mc][McToJobLengh][0]]['訂單編號'].values[0] +' 延誤 ' + str(Each_Job_Delayed_days)+'天' )
    #             Delayed_days_df.append(temp_Delayed_days_list)
    #         # print(Num_Mc,timeStamp)
    #         temp.append(Due_timeStamp)
    #     Over_Due_List.append(temp)
    # print('在%s筆訂單內有 %s延誤'%(Number_Of_Job,count))
    # print('延誤筆數/總訂單筆數 比率 %s'%(count/Number_Of_Job))
    # print('總延誤天數為 %s' % Add_up_Delayed_days)
    #
    # '''0604計算Zip_OrderNumber_Stime_Endtime 的換模開始及結束(為了畫圖需要)'''
    # #開始、結束時間計算
    #
    # Record_Job_Index = []
    # for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime)-1):
    #     temp_Job_Index =[]
    #     temp_Setup_Start = []
    #     temp_Setup_End   = []
    #     for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
    #         # print(Nb_Of_Machine_index,Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0])
    #         temp_Job_Index.append(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][0])
    #     Record_Job_Index.append(temp_Job_Index)
    #
    # Record_Setup_St_Endt_Index = []
    # Record_Job_PPN = Setup_Machine_Corres_To_Job_PPN(Record_Job_Index)  #Zip處理成Record_Job_PPN
    #
    # temp_Record_Setup_St_Endt_Index = []
    # temp_Mc_Corres_To_Setup_Index =[]
    # temp_Mc_Corres_To_Setup_Start =[]
    # temp_Mc_Corres_To_Setup_End = []
    # Setup_count = 0  #記錄換幾次模
    # Setup_times = 0  #紀錄整體換模時間
    # for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime)-1):
    #     temp_Setup_Index = []
    #     temp_Setup_Start = []
    #     temp_Setup_End   = []
    #     for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index])):
    #         # print(Nb_Of_Machine_index,Nb_Of_Job_index )
    #         if Nb_Of_Job_index == 0:
    #             '''每個機台的第一筆Job皆不用考慮換模時間'''
    #             Setup_Start = Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2]
    #             Setup_End   = Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][2]
    #             temp_Setup_Index.append("Setup")
    #             temp_Setup_Start.append(Setup_Start)
    #             temp_Setup_End.append(Setup_End)
    #         else:
    #             '''每個機台的第二筆Job後皆考慮換模時間，以下的換模時間有很多種方法，最好的方法是"用字典第一代"方式搜尋時間最快'''
    #
    #             '''要考慮模具'''
    #
    #             '''用字典第一代'''
    #             if Record_Job_PPN[Nb_Of_Machine_index][Nb_Of_Job_index-1] != Record_Job_PPN[Nb_Of_Machine_index][Nb_Of_Job_index]:
    #                 Setup = Mould_dict[Record_Job_PPN[Nb_Of_Machine_index][Nb_Of_Job_index-1]+Record_Job_PPN[Nb_Of_Machine_index][Nb_Of_Job_index]]*60
    #                 Setup_count += 1
    #                 Setup_times += Setup
    #             else:
    #                 Setup = 0
    #
    #
    #             '''考慮模具 所以Job_Start要  Max(前一個Job結束時間+換模時間,此Job的ArrivalTime) 取大的當作目前Job的開始時間'''
    #             Setup_End  = Zip_OrderNumber_Stime_Endtime[Nb_Of_Machine_index][Nb_Of_Job_index][1] #Setup結束時間 = 目前Job開始時間
    #             Setup_Start= Setup_End - Setup #加模具
    #
    #             temp_Setup_Index.append("Setup")
    #             temp_Setup_Start.append(Setup_Start)
    #             temp_Setup_End.append(Setup_End)
    #
    #     temp_Mc_Corres_To_Setup_Index.append(temp_Setup_Index)
    #     temp_Mc_Corres_To_Setup_Start.append(temp_Setup_Start)
    #     temp_Mc_Corres_To_Setup_End.append(temp_Setup_End)
    #     # Record_Setup_St_Endt_Index.append(temp_Job_Start)
    #
    # #將訂單編號、Job開始時間、結束時間綁在一起[(訂單Index,開始時間,結束時間)]
    # Record_Setup_St_Endt_Index = OrderNumber_Stime_Endtime_Zip(Name = temp_Mc_Corres_To_Setup_Index, St = temp_Mc_Corres_To_Setup_Start, End = temp_Mc_Corres_To_Setup_End, eachMaxmakespan=0)
    #
    # #下面這個將原本Zip改成 List [[訂單Index,開始時間,結束時間]]以方便畫圖時把ST、End時間做更改 因tuple 無法更改!!
    # Record_Setup_St_Endt_Index = Zip_OrderNumber_Stime_Endtime_converList(Record_Setup_St_Endt_Index)
    #
    #
    #
    # # # 開啟輸出的 CSV 檔案
    # # with open('真實資料無ArrivalTime_output.csv', 'w', newline='') as csvfile:
    # #   # 建立 CSV 檔寫入器
    # #   writer = csv.writer(csvfile)
    # #
    # #   # 寫入一列資料
    # #   writer.writerow(['訂單編號', '延誤天數'])
    # #
    # #   # 寫入另外幾列資料
    # #   for i in range(len(Delayed_days_df)):
    # #       writer.writerow(Delayed_days_df[i])
    # #
    # # # New_OverDue_percent = count/Number_Of_Job
    # #
    # # # OverDue_percent = abs((New_OverDue_percent - Old_OverDue_percent)/int(Old_OverDue_percent))
    # # # print("進步多少：%s" % OverDue_percent )
    # #
    # #
    # df = []
    # def draw_machine(index,value):
    #     '''把Zip_OrderNumber_Stime_Endtime排序過後的第一條拿來畫圖!!'''
    #
    #     endt = TakeEnd_minute_Func(Zip_OrderNumber_Stime_Endtime)
    #     # print(endt)
    #     st   = TakeSt_minute_Func(Zip_OrderNumber_Stime_Endtime)
    #     # print(st)
    #     each_Job  = [j for j in Zip_OrderNumber_Stime_Endtime[index]]
    #
    #
    #
    #     for i in range(len(each_Job)):
    #         if each_Job[i][0] == 'setup':
    #             df.append(dict(Task='%s'% 'setup', Start='%s' % st[index][i],
    #             Finish='%s'  % endt[index][i],Resource= value,item = "換模",Mould = "換模"))
    #         elif each_Job[i][0] == 'Stop':
    #             df.append(dict(Task='%s'% 'Stop', Start='%s' % st[index][i],
    #             Finish='%s'  % endt[index][i],Resource= value,item = "停機",Mould = "停機"))
    #         else:
    #             df.append(dict(Task='%s'% OrderName[each_Job[i][0]], Start='%s' % st[index][i],
    #             Finish='%s'  % endt[index][i],Resource= value,item = Product_Part_Number[each_Job[i][0]]))
    #
    #
    #     return df
    #
    # df_setup = []
    # def draw_Setup_machine(index,value):
    #     '''把Zip_OrderNumber_Stime_Endtime排序過後的第一條拿來畫圖!!'''
    #
    #     endt = TakeEnd_minute_Func(Record_Setup_St_Endt_Index)
    #     # print(endt)
    #     st   = TakeSt_minute_Func(Record_Setup_St_Endt_Index)
    #     # print(st)
    #     each_Job  = [j for j in Record_Setup_St_Endt_Index[index]]
    #
    #     for i in range(len(each_Job)):
    #
    #         df_setup.append(dict(Task='Setup', Start='%s' % st[index][i],
    #         Finish='%s'  % endt[index][i],Resource= value,item = '換模'))
    #
    #
    #     return df_setup
    #
    # for index,value in enumerate(All_Machines_Name):
    #     # print(index,value)
    #     draw_Setup_machine(index,value)
    # #------------------------------------下面是開始時間的換算----------------------------------------------#
    #
    # def TakeSt_minute_Func(Nice_Chromosome):
    #     #---將 每個Job對應到的開始時間讀取出來到下一個stTime_out_Func---#
    #     temp_List  =[] #紀錄有50台機台
    #     for Num_Mc in range(Number_Of_Machine):
    #
    #         temp = [] #用temp主要是因為有50台機台想要紀錄每台機台內的Job數
    #         for Num_Job in range(len(Nice_Chromosome[Num_Mc])):
    #             # print('第 %s Job :開始時間 %s' % (Num_Job,Nice_Chromosome[Num_Mc][Num_Job][1]))
    #             TakeSt_minute = converge_Time(stsec = Nice_Chromosome[Num_Mc][Num_Job][1] + RecoverTime[Num_Mc])
    #             # print(TakeSt_minute)    #換算好的時間戳%Y-%m-%d %H:%M:%S
    #             temp.append(TakeSt_minute)
    #         temp_List.append(temp)
    #
    #     # print(temp_List)
    #     return temp_List #回傳至draw_machine函數裡的st
    #
    #
    # def converge_Time(stsec):
    #     #---將converge_Time 把時間戳，轉回"%Y-%m-%d %H:%M:%S"---#
    #         Time_stamp = stsec # 設定timeStamp
    #         struct_time = time.localtime(Time_stamp) # 轉成時間元組
    #         timeString = time.strftime("%Y-%m-%d %H:%M:%S", struct_time) # 轉成字串
    #         # print(timeString)   #將時間戳轉換成%Y-%m-%d %H:%M:%S
    #         return timeString
    # #------------------------------------下面是結束時間的換算----------------------------------------------#
    #
    # def TakeEnd_minute_Func(Nice_Chromosome):
    #     #---將 每個Job對應到的結束時間讀取出來到下一個EndTime_out_Func---#
    #     temp_List  =[] #紀錄有50台機台
    #     for Num_Mc in range(Number_Of_Machine):
    #         # print('機台 %s' % Num_Mc)
    #
    #         temp = []   #用temp主要是因為有50台機台想要紀錄每台機台內的Job數
    #         for Num_Job in range(len(Nice_Chromosome[Num_Mc])):
    #             # print('第 %s Job :結束時間 %s' % (Num_Job,Nice_Chromosome[Num_Mc][Num_Job][2]))
    #             TakeEnd_minute = converge_Time_Func(Endsec = Nice_Chromosome[Num_Mc][Num_Job][2] + RecoverTime[Num_Mc])
    #             # print(TakeEnd_minute)  #換算好的時間戳%Y-%m-%d %H:%M:%S
    #             temp.append(TakeEnd_minute)
    #         temp_List.append(temp)
    #
    #     # print(temp_List)
    #     return temp_List   #回傳至draw_machine函數裡的endt
    #
    #
    #
    # def converge_Time_Func(Endsec):
    #     #---將converge_Time 把時間戳，轉回"%Y-%m-%d %H:%M:%S"---#
    #     Time_stamp = Endsec # 設定timeStamp
    #     struct_time = time.localtime(Time_stamp) # 轉成時間元組
    #     timeString = time.strftime("%Y-%m-%d %H:%M:%S", struct_time) # 轉成字串
    #     # print(timeString)  #將時間戳轉換成%Y-%m-%d %H:%M:%S
    #     return timeString
    # # =============================================================================
    #
    #
    #
    #
    # for index,value in enumerate(All_Machines_Name):
    #     draw_machine(index,value)
    #
    #
    # #呈現圖表
    # fig = px.timeline(df+df_setup, x_start="Start", x_end="Finish",y="Resource", hover_name = "item" ,hover_data=["item"],color="item",color_discrete_map = colorss,text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    # # fig = px.timeline(df, x_start="Start", x_end="Finish",y="Resource",color = "Task", hover_name = "Task",text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    # fig.update_traces(textposition='inside',marker_line_color='rgb(8,48,107)')
    # # fig.update_layout(yaxis={'categoryorder':'category ascending'})  #將Y軸由小排到大
    # #將Y軸由大排到小
    # fig.update_yaxes(categoryorder='array', categoryarray= All_Machine_uniqe)
    # fig.update_yaxes(autorange="reversed")
    # fig.update_layout( title={   # 设置整个标题的名称和位置
    #         "text":"甘特圖",
    #         "y":0.96,  # y轴数值
    #         "x":0.5,  # x轴数值
    #         "xanchor":"center",  # x、y轴相对位置
    #         "yanchor":"top"
    #     })
    #
    #
    #
    # # fig.show()
    #
    # # plot(fig)
    #
    # '''下載出甘特圖及分解'''
    # GanttChart = pd.DataFrame(df+df_setup)
    # GanttChart = GanttChart.sort_values(by=['Resource'])
    # GanttChart.reset_index(drop=True, inplace=True)
    # GanttChart["Group"] = 0  # 創一欄 Group 欄位 裡面都放0
    # All_Machine_uniqe = GanttChart["Resource"]
    # All_Machine_uniqe = All_Machine_uniqe.drop_duplicates().tolist()
    # # print(GanttChart)
    # Machine_group_temp = 0
    # for Machine_group in range(len(GanttChart)):
    #     '''可嘗試把 All_Machines_Name改成有畫圖的機台˙就好 可能不需要全跑'''
    #     if Machine_group == 0:
    #         GanttChart["Group"][Machine_group] = Machine_group_temp
    #     elif GanttChart["Resource"][Machine_group-1] != GanttChart["Resource"][Machine_group]:
    #         Machine_group_temp +=1
    #         GanttChart["Group"][Machine_group] = Machine_group_temp
    #     else:
    #         GanttChart["Group"][Machine_group] = Machine_group_temp
    #
    #     # GanttChart["Group"][GanttChart["Resource"] == All_Machines_Name[Machine_group]] = Machine_group
    # Machinecount = 0;
    #
    # prev_temp = -1
    # lower = -1
    # upper = 0
    #
    # for i in range(int(np.ceil(len(GanttChart["Resource"].unique())/plot_Machinenumber))):
    #
    #     temp = plot_Machinenumber * (i+1) -1;
    #     # print(temp)
    #     if temp >= int(len(GanttChart["Group"].unique())):
    #         Machinecount = np.hstack((Machinecount,(len(GanttChart))+1))  #0526有改len(部分)因為影響到plot_endpoint
    #         # print(Machinecount)
    #         break;
    #     temp1 = GanttChart[GanttChart.Group ==temp].index.max()
    #     # print(temp1)
    #     Machinecount = np.hstack((Machinecount,temp1+1))
    #     # print(Machinecount)
    #
    # plot_startpoint = Machinecount[:-1]
    # plot_endpoint = Machinecount[1:]
    #
    # for i in range(len(plot_startpoint)):
    #     # fig = px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i],],x_start= 'Start',x_end = 'Finish', y='Resource', hover_name = "item" ,hover_data=["item", "Mould"],color="item")
    #     # fig = px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i],], x_start="Start", x_end="Finish",y="Resource", hover_name = "item" ,hover_data=["item", "Mould"],color="item",color_discrete_map = colors,text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    #     #呈現圖表
    #     fig = px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i],], x_start="Start", x_end="Finish",y="Resource", hover_name = "item" ,hover_data=["item"],color="item",color_discrete_map = colorss ,text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    #     # fig =px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i],], x_start="Start", x_end="Finish",y="Resource", hover_name = "item" ,hover_data=["item"],color="item",color_discrete_map = colorss,text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    #     fig.update_traces(textposition='inside',marker_line_color='rgb(8,48,107)')
    #
    #     # fig.update_layout(yaxis={'categoryorder':'category ascending'})  #將Y軸由小排到大
    #     #將Y軸由大排到小
    # # =============================================================================
    # #
    # #     fig.update_layout(yaxis={'categoryorder':'category descending'}, title={   # 设置整个标题的名称和位置
    # #             "text":"甘特圖",
    # #             "y":0.96,  # y轴数值
    # #             "x":0.5,  # x轴数值
    # #             "xanchor":"center",  # x、y轴相对位置
    # #             "yanchor":"top"
    # #         })
    # #     py.offline.plot(fig, filename='C:/Users/ZhiXiang/Downloads/未前插'+ i.__str__()+'.html')    #下載html的路徑  可以自行更改路徑
    # # =============================================================================
    #     fig.update_yaxes(categoryorder='array', categoryarray= All_Machine_uniqe)
    #     fig.update_yaxes(autorange="reversed")
    #     fig.update_layout( title={   # 设置整个标题的名称和位置
    #             "text":"甘特圖",
    #             "y":0.96,  # y轴数值
    #             "x":0.5,  # x轴数值
    #             "xanchor":"center",  # x、y轴相对位置
    #             "yanchor":"top"
    #         })
    #     # py.offline.plot(fig, filename='C:/Users/ZhiXiang/Downloads/法二傳統GA已前插'+ i.__str__()+'.html')    #下載html的路徑  可以自行更改路徑
    #
    #
    #
    #
    #
    #
    #
    #
    #
    # # =============================================================================
    # # '''0531以前舊的下載出甘特圖及分解'''
    # # GanttChart = pd.DataFrame(df)
    # # GanttChart["Group"] = 0  # 創一欄 Group 欄位 裡面都放0
    # # # print(GanttChart)
    # # Machine_group_temp = 0
    # # for Machine_group in range(len(GanttChart)):
    # #     '''可嘗試把 All_Machines_Name改成有畫圖的機台˙就好 可能不需要全跑'''
    # #     if Machine_group == 0:
    # #         GanttChart["Group"][Machine_group] = Machine_group_temp
    # #     elif GanttChart["Resource"][Machine_group-1] != GanttChart["Resource"][Machine_group]:
    # #         Machine_group_temp +=1
    # #         GanttChart["Group"][Machine_group] = Machine_group_temp
    # #     else:
    # #         GanttChart["Group"][Machine_group] = Machine_group_temp
    # #     # GanttChart["Group"][GanttChart["Resource"] == All_Machines_Name[Machine_group]] = Machine_group
    # # Machinecount = 0;
    # #
    # # for i in range(int(np.ceil(len(GanttChart["Resource"].unique())/plot_Machinenumber))):
    # #
    # #     temp = plot_Machinenumber * (i+1) -1;
    # #     # print(temp)
    # #     if temp >= int(len(GanttChart["Group"].unique())):
    # #         Machinecount = np.hstack((Machinecount,(len(Order))+1))
    # #         # print(Machinecount)
    # #         break;
    # #     temp1 = GanttChart[GanttChart.Group ==temp].index.max()
    # #     # print(temp1)
    # #     Machinecount = np.hstack((Machinecount,temp1+1))
    # #     # print(Machinecount)
    # #
    # # plot_startpoint = Machinecount[:-1]
    # # plot_endpoint = Machinecount[1:]
    # #
    # # for i in range(len(plot_startpoint)):
    # #     # fig = px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i],],x_start= 'Start',x_end = 'Finish', y='Resource', hover_name = "item" ,hover_data=["item", "Mould"],color="item")
    # #     # fig = px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i],], x_start="Start", x_end="Finish",y="Resource", hover_name = "item" ,hover_data=["item", "Mould"],color="item",color_discrete_map = colors,text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    # #     #呈現圖表
    # #     fig = px.timeline(GanttChart.iloc[plot_startpoint[i]:plot_endpoint[i],], x_start="Start", x_end="Finish",y="Resource", hover_name = "item" ,hover_data=["item"],color="item",color_discrete_map = colorss ,text = "Task",color_discrete_sequence=px.colors.qualitative.Plotly, title="甘特圖")
    # #     fig.update_traces(textposition='inside',marker_line_color='rgb(8,48,107)')
    # #     # fig.update_layout(yaxis={'categoryorder':'category ascending'})  #將Y軸由小排到大
    # #     #將Y軸由大排到小
    # #     fig.update_yaxes(categoryorder='array', categoryarray= All_Machine_uniqe)
    # #     fig.update_yaxes(autorange="reversed")
    # #     fig.update_layout(title={   # 设置整个标题的名称和位置
    # #             "text":"甘特圖",
    # #             "y":0.96,  # y轴数值
    # #             "x":0.5,  # x轴数值
    # #             "xanchor":"center",  # x、y轴相对位置
    # #             "yanchor":"top"
    # #         })
    # #
    # #     py.offline.plot(fig, filename='C:/Users/ZhiXiang/Downloads/新方法已前插'+ i.__str__()+'.html')    #下載html的路徑  可以自行更改路徑
    # #
    # # # plot(fig)
    # # =============================================================================
    # #
    # # =============================================================================
    # '''印出最好的機台'''
    # # =============================================================================
    # # best_Opt_Machine =[]
    # # opt_machine = pd.DataFrame(GanttChart["Resource"].unique())
    # # for i in range(len(opt_machine)):
    # #     # print(i)
    # #     temp_best_Storage_Opt_Machine ={}
    # #     for index,value in enumerate(Mc_model_Quantity):
    # #         '''每一條染色體一開始計算可用機台有哪幾台!!'''
    # #         if value == Machine_model_1:
    # #             temp_best_Storage_Opt_Machine[Machine_model_1] = opt_machine[opt_machine[0].str.contains(Machine_model_1)][0].values.tolist()
    # #         else:
    # #             temp_best_Storage_Opt_Machine[Machine_model_2] = opt_machine[opt_machine[0].str.contains(Machine_model_2)][0].values.tolist()
    # # best_Opt_Machine.append(temp_best_Storage_Opt_Machine)
    # #
    # # with open('C:/Users/ZhiXiang/Downloads/最佳化機台.csv', 'w') as f:
    # #     writer = csv.writer(f)
    # #     for i in range(len(best_Opt_Machine)):
    # #         for k, v in best_Opt_Machine[i].items():
    # #            writer.writerow([k, v])
    # #            writer.writerow([k,len(v)])
    # # =============================================================================
    #
