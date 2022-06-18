# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 16:26:25 2022

@author: TiffanyTang
"""
import re
import csv
import pandas as pd
import numpy as np
import math
import time
import random

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# import Draw_GanttChart as Draw
import plotly.express as px
from collections import Counter
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
from plotly.figure_factory import create_gantt
# from operator import itemgetter
import datetime
import plotly as py
from itertools import compress  #創三維空間大量資料回傳List 回傳index很有用
import xlrd
from numba import jit
import pickle
import Data_Normalization as DN;         #嘉謙paper正歸化目標資料(已將套件改成適合此案例方式)
import threading
import pdb
import pickle
import copy
with open('test.pickle_List', 'rb') as f:
    List_OrderNumber_Stime_Endtime = pickle.load(f)
    
Deecopy_List_OrderNumber_Stime_Endtime = copy.deepcopy(List_OrderNumber_Stime_Endtime)
with open('test.pickle_temp_D1', 'rb') as f:
    temp_Select_Fitness_Df1 = pickle.load(f)

Deecopy_Select_Fitness_Df1 = copy.deepcopy(temp_Select_Fitness_Df1)



nb_threads = 10
# =============================================================================
# # from test_Xiang import Forward_Insertion_Method_Improved
# =============================================================================
# import TOPSISs as tp;         #TOPSIS 演算法(已將套件改成適合此案例方式)
# import plotly.io as pio

Screening_of_site = "OP3"   #將BOM表"工序說明"篩選OP3
Machine_model_1 = 'PM2VA001'  #篩選機台銑床型號此OP3篩選出'PM2VA001'
Machine_model_2 = 'PM9VA001'  #篩選機台銑床型號此OP3篩選出'PM9VA001'

Machine_model_Number= []  #為了去算有幾種機型
Machine_model_Number.append(Machine_model_1)
Machine_model_Number.append(Machine_model_2)

Number_Of_Job = 1874  #0519改成最新版1874筆資料   362
Number_Of_Machine = 381    #平行機台雖然最多是381 但用限制式 把有些訂單原本可以丟到PM9VA001的機台都重新選機至其他台
Number_Of_Chromosome = 10
Number_Of_Obj = 3

Chrossover_Rate = 0.8
Mutation_Rate   = 0.2
'''learning_rate   幫助交配範圍大一點 (建議0.2~0.25) 
  【主要看訂單數決定 Ex 訂單數10個 若組成基因長度 會是 10種機型數(開始、結束)】 共 10個基因 大概交配一半以上就算多
  若訂單為3716筆，learning_rate 用0.25
  若訂單為10筆，learning_rate 用0.1

'''
learning_rate   = 0.25          

set_Timer =  90 * 60 #分鐘*秒數
num_iterations = 1 #迭代數
Continuous_termination =500 #連續Obj不變幾代就要停止



'''為了計算目標值 TSTP = Total completion time * Penalty_value_weight_TCT + Sum_Tardiness * Penalty_value_weight + Setup_times * Penalty_value_weight_Setup'''
Penalty_value_weight_Machine = 0.1    #Total_completion_time 權重設為1
Penalty_value_weight_TT = 0.6      #Tardiness  權重設為1000
Penalty_value_weight_Setup = 0.3 #Setup_times  權重設為100

Penalty_Weights = np.array([Penalty_value_weight_Machine,Penalty_value_weight_TT,Penalty_value_weight_Setup]);
Penalty_Weights = Penalty_Weights/Penalty_Weights.sum();  #讓權重總和為 1 

#[#1.5] Plot [幾個機台畫成一張圖]
plot_Machinenumber = 20;

'''NSGAII參數'''
#[#1.4] NSGA_II 參數 Rules variables
weight1 = 1; # Target_Machineusing_KPI_Weight
weight2 = 1; # Target_DeliveryRate_KPI_Weight
weight3 = 1; # Target_ChangeLine_KPI_Weight
# weight4 = 1; # Target_WIP_KPI_Weight

Weights = np.array([weight1,weight2,weight3]);
Weights = Weights/Weights.sum();  #讓權重總和為 1 

# =============================================================================
# '''TOPSIS參數'''
# #[#1.4] TPOSIS 參數 Rules variables
# TOPSISweight1 = 0.1; # Target_Machineusing_KPI_Weight
# TOPSISweight2 = 0.6; # Target_DeliveryRate_KPI_Weight  #0605改0.3
# TOPSISweight3 = 0.3; # Target_ChangeLine_KPI_Weight    #0605改0.6
# # # weight4 = 1; # Target_WIP_KPI_Weight
# 
# TOPSIS_Weights = np.array([TOPSISweight1,TOPSISweight2,TOPSISweight3]);
# TOPSIS_Weights = TOPSIS_Weights/TOPSIS_Weights.sum();  #讓權重總和為 1 
# 
# criterias = ["-"  , "-" , "-" ];  #正反指標 (看指標是望大還是望小特性)
# =============================================================================




# ===================================學校實驗室==========================================

'''讀取BOM表'''
BOM = pd.read_csv("0501_未整理BOM.csv",encoding= 'ANSI',dtype ={'orderBomKey':str,'工序代碼':str,
                                                    '工序說明':str,
                                                    'seqNo':int,
                                                    '加工時間':np.float64,
                                                    '設備/資源':str,
                                                    })
'''將BOM表作工序篩選Screening_of_site = OP3'''
BOM_filt = BOM["工序說明"] == Screening_of_site
BOM_DF = BOM.loc[BOM_filt]
BOM_DF.reset_index(drop=True, inplace=True) #重新更新Index_0424


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
    dict_canRunTool[BOM_DF.iloc[i]['orderBomKey']] =BOM_DF.iloc[i]["設備/資源"]

for index,BomKey in enumerate(dict_canRunTool.keys()):
    '''(new)產品料號對應機台型號有多少機台數'''
    '''使用成字典BomKey 對應 機台，把String 改成List，方便以後選機'''
    dict_canRunTool[BomKey] = cut_text(dict_canRunTool[BomKey])

dict_BomKey_Processing = {} 
for i in range(len(BOM_DF)):
    '''使用成字典ProcessTime'''
    dict_BomKey_Processing[BOM_DF.iloc[i]['orderBomKey']] = np.round(BOM_DF.iloc[i]["加工時間"],3)

# =============================================================================
# Order_pt = []
# for i in range(len(Order)):
#     '''處理Order對應的ProcessTime 成List'''
#     Order_pt.append(dict_BomKey_Processing[Order.at[i,"orderBomKey"]])
# Order_pt = pd.DataFrame(Order_pt)  #將List轉換成DataFrame
# =============================================================================

Order = pd.merge(Order,BOM_DF.loc[:,["orderBomKey","設備/資源"]],how = 'inner',on = "orderBomKey") #新增對應的可用機台Order表
Order = pd.merge(Order,BOM_DF.loc[:,["orderBomKey","加工時間"]],how = 'inner',on = "orderBomKey")  #新增加工時間至Order表

Order["總加工時間"] = round(Order['目標數量']*Order['加工時間']*60,3)   #轉成秒數


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
Mould["當前_目標"] = Mould.apply(lambda x: x['當前產品料號']+x['目標產品料號'],axis = 1)
Mould_dict = Mould[['當前_目標','換線時間']]
Mould_dict.index = Mould['當前_目標']
Mould_dict = Mould_dict.drop('當前_目標',axis =1)
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
print('模具時間',end-start)

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
    timeStamp=int(time.mktime(time.strptime(str_time,"%Y/%m/%d %H:%M")))
    Arrival_Time.append(timeStamp)
'''資料預處理-RecoverTime'''
RecoverTime = []
for i in range(len(Machine_Start)):
    str_time = str(Machine_Start.iloc[i]['機台開始時間'])
    timeStamp=int(time.mktime(time.strptime(str_time,"%Y/%m/%d %H:%M")))
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
    timeStamp=int(time.mktime(time.strptime(str_time,"%Y/%m/%d %H:%M")))
    Due_Date_times.append(timeStamp)

def read_excel_data():
    '''甘特圖106種產品料號顏色整理'''
    filename = '顏色產品料號.xlsx'
    data = xlrd.open_workbook(filename)
    table = data.sheet_by_name('工作表1')
    row_num = table.nrows  # 行数
    # col_num = table.ncols  # 列数
    datas = dict([]) # 这步也要转字典类型
    for i in range(row_num):
        colorss = dict([table.row_values(i)]) # 这一步就要给它转字典类型，不然update没法使用
        datas.update(colorss)
    # print(datas)
    return datas
if __name__ == "__main__":
    colorss = read_excel_data()

#%%交配突變
def single_point_crossover(A,B,X):
    A_new = np.append(A[:X],B[X:])
    B_new = np.append(B[:X],A[X:])
    # print(A_new)
    # print(B_new)
    return A_new, B_new


def multi_point_crossover(A,B,X):
    for i in X :
        A ,B = single_point_crossover(A,B,i)
    return A, B


def Crossover (C,D):
    '''交配，採兩點交配!!!'''
    Crossover1 = C
    Crossover2 = D
    '''例如 有10個Job ,最多可能交配點為  [1,9]，不會到10,否則就變成單點交配'''
    #這裡到最後去思考會不會可用機台數基因最後一個結束時間不會變更!!
    X = np.random.choice(range(1,Number_Of_Job*2),size = 2,replace = False)
    '''下面說明當X任兩點若小於Number_Of_Job * (1+learning_rate) 範圍 就不做交配!!，learning_rate在上面可以做變更'''
    while ( abs(X[0]-X[1]) < round((Number_Of_Job)* (learning_rate))) : #0602修改不加1的learning_rate 因為基因數沒有到很長
        X = np.random.choice(range(1,Number_Of_Job*2 ),size = 2,replace = False)
    return Crossover1,Crossover2,X

# def Mutation(index):
#     '''採2點突變法'''
#     '''突變採任兩點取出不放回做突變'''
#     temp_position1,temp_position2 =np.random.choice(range(1,Number_Of_Job*2),size = 2,replace = False)
#     temp_RandomNumber1,temp_RandomNumber2 = [random.random() for i in range(2)]
#     # print("第幾格被換 %d ,換成哪個一個亂數%f"%(temp_position,temp_RandomNumber))

#     Temp_Mutation[index][temp_position1] = temp_RandomNumber1
#     Temp_Mutation[index][temp_position2] = temp_RandomNumber2
    
#     # print('新的',Temp_Mutation[index])

#     return Temp_Mutation


def Mc_Mutation(Mc_Gene):
    '''採2點突變法'''
    '''突變採任兩點取出不放回做突變'''
    # temp_position1,temp_position2 =np.random.choice(range(0,4),size = 2,replace = False)
    # temp_RandomNumber1,temp_RandomNumber2 = [random.random() for i in range(2)]
    # Mc_Gene[temp_position1] = temp_RandomNumber1
    # Mc_Gene[temp_position2] = temp_RandomNumber2
    
    temp_position1 =np.random.choice(range(0,4),size = 1,replace = False)
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
    lst_2 = [x + [0]*(n-len(x)) for x in List] 
    # print(lst_2)
    matrix = np.array(lst_2)/(60*60*24) #已換成天數單位
    return matrix

def ArrivalTime(Arrivaltime_stamp):
    '''將染色體Part2 改成使用ArrivalTime排序'''
    timeString = Arrivaltime_stamp# 時間格式為字串
    struct_time = time.strptime(timeString, "%Y/%m/%d %H:%M") # 轉成時間元組
    time_stamp = int(time.mktime(struct_time)) # 轉成時間戳
    return time_stamp

def Machine_Corres_To_Job_ArrivalTime_Pt_PPN(temp_Machine_Corresponds_To_Job_sorted):
    '''直接用temp_Machine_Corresponds_To_Job_sorted 的Job index 去抓取 預處理的ArrivalTime 、PT 各自時間'''
    '''處理成相對位置!!'''
    temp_ArrivalTime= []
    temp_Mc_Corres_To_Job_Pt = []
    temp_Mc_Corres_To_Job_Product_Part_Number = []

    #總共要跑幾次機台
    for Num_Mc in range(len(temp_Machine_Corresponds_To_Job_sorted)):
        temp = []
        temp_TotalTime = []
        temp_Product_PN = []
        #每個機台內跑幾個Job數
        for  McToJobLengh in range(len(temp_Machine_Corresponds_To_Job_sorted[Num_Mc])):
            
            time_stamp  = Arrival_Time[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh]] - RecoverTime[Num_Mc]# 轉成時間戳(有減機台相對位置)

            temp_TotalTime.append(Process_Time[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh]])  #temp_TotalTime增加每個Job總加工時間
            
            temp_Product_PN.append(Product_Part_Number[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh]])  #temp_Product_PN增加每個Job的產品料號
            if time_stamp < 0 :
                '''若time_stamp <= 0 ，轉換為0 ，代表貨到達此刻可以加工'''
                time_stamp = 0

            temp.append(time_stamp)
            
        temp_ArrivalTime.append(temp)
        temp_Mc_Corres_To_Job_Pt.append(temp_TotalTime)
        temp_Mc_Corres_To_Job_Product_Part_Number.append(temp_Product_PN)

    return temp_ArrivalTime,temp_Mc_Corres_To_Job_Pt,temp_Mc_Corres_To_Job_Product_Part_Number

def Setup_Machine_Corres_To_Job_PPN(temp_Machine_Corresponds_To_Job_sorted):
    '''0526  直接用Zip_OrderNumber_Stime_Endtime 的Job index 去抓取 預處理的各自Job_Product_Part_Number時間'''
    '''換模專用'''
    temp_Mc_Corres_To_Job_Product_Part_Number = []

    #總共要跑幾次機台
    for Num_Mc in range(len(temp_Machine_Corresponds_To_Job_sorted)):
        temp_Product_PN = []
        #每個機台內跑幾個Job數
        for  McToJobLengh in range(len(temp_Machine_Corresponds_To_Job_sorted[Num_Mc])):
            
            temp_Product_PN.append(Product_Part_Number[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh]])  #temp_Product_PN增加每個Job的產品料號

        temp_Mc_Corres_To_Job_Product_Part_Number.append(temp_Product_PN)

    return temp_Mc_Corres_To_Job_Product_Part_Number

def Thebest_Machine_Corres_To_Job_ArrivalTime_Pt_PPN(temp_Machine_Corresponds_To_Job_sorted):
    '''處理成相對位置!!'''
    temp_Mc_Corres_To_Job_Index =[]
    temp_ArrivalTime= []
    temp_Mc_Corres_To_Job_Pt = []
    temp_Mc_Corres_To_Job_Product_Part_Number = []
    
    #總共要跑幾次機台
    for Num_Mc in range(len(temp_Machine_Corresponds_To_Job_sorted)-1):
        
        temp_Job_index = []
        temp_Arr = []
        temp_TotalTime = []
        temp_Product_PN  = []
        #每個機台內跑幾個Job數
        for  McToJobLengh in range(len(temp_Machine_Corresponds_To_Job_sorted[Num_Mc])):
            # print("機台%s, 第 %s 個Job" % (Num_Mc,McToJobLengh))
            if temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][0] == 'Stop':
                temp_Arr.append(temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][1] -  RecoverTime[Num_Mc])
            else:
                temp_Job_index.append(temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][0]) #新增Job的Index
                
                temp_TotalTime.append(Process_Time[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][0]]) #新增總加工時間
                
                temp_Product_PN.append(Product_Part_Number[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][0]]) #新增產品料號
                
                time_stamp  = Arrival_Time[temp_Machine_Corresponds_To_Job_sorted[Num_Mc][McToJobLengh][0]] - RecoverTime[Num_Mc]# 到達時間轉成時間戳(有減機台相對位置)

                if time_stamp < 0 :
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

    return temp_Mc_Corres_To_Job_Index,temp_ArrivalTime,temp_Mc_Corres_To_Job_Pt,temp_Mc_Corres_To_Job_Product_Part_Number

#%%'''計算非空機台，意旨真的機台數量'''
def Culate_Num_Mc_Actual_Open(Chromosom_Loop):
    '''計算非空機台'''
    Num_Machine_Actual_Open =[]
    for i in range(len(Chromosom_Loop)):
        temp_Machine_Empty = []
        temp_Number_Total_Empty = 0
        temp_Open_number_mc =len(Chromosom_Loop[i])
        for j in range(len(Chromosom_Loop[i])):
            if not Chromosom_Loop[i][j]:
                # print(f'{i}條染色體{j}List is empty')
                temp_Machine_Empty.append(j)
                temp_Number_Total_Empty +=1
        temp_Open_number_mc = temp_Open_number_mc - temp_Number_Total_Empty
        Num_Machine_Actual_Open.append(temp_Open_number_mc)
    return Num_Machine_Actual_Open

#%%將訂單編號、開始結束綁再一起 [(J1,St1,Et1),(J2,St2,Et2)]
def OrderNumber_Stime_Endtime_Zip(Name,St,End,eachMaxmakespan):
    '''將訂單編號、開始結束綁再一起 [(J1,St1,Et1),(J2,St2,Et2)]'''
    temp_OrderNumber_Stime_Endtime = []
    for Num_Mc in range(len(End)):
        temp = list(zip(Name[Num_Mc],St[Num_Mc],End[Num_Mc]))
        temp_OrderNumber_Stime_Endtime.append(temp)
    temp_OrderNumber_Stime_Endtime.append(eachMaxmakespan)
    return temp_OrderNumber_Stime_Endtime

#%%將(訂單編號、St、Endt)改成[[訂單編號、St、Endt]]
def Zip_OrderNumber_Stime_Endtime_converList(Zip_converList,eachMaxmakespan = 0):
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

#%%將訂單編號、開始結束綁再一起 [(J1,ArrivalT1,Pt1),(J2,ArrivalT2,Pt2)]

def Thebest_OrderNumber_St_Pt_Zip(Name,ArrivalT,Pt):
    '''此Function跟上面OrderNumber_St_Pt_Zip差異是 這個函數使用在最好的那條解碼身上!!'''
    temp_OrderNumber_Stime_Endtime = []
    for Num_Mc in range(len(Pt)):
        temp = list(zip(Name[Num_Mc],ArrivalT[Num_Mc],Pt[Num_Mc]))
        temp_OrderNumber_Stime_Endtime.append(temp)
    
    return temp_OrderNumber_Stime_Endtime



#%%將(訂單編號、ArrivalT、Endt)改成[[訂單編號、ArrivalT、Endt]]

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
    for i in range(len(List)-1):
        try:
            a = pd.DataFrame(List[i])
            temp_sum = a[2].sum()
            temp.append(temp_sum)
        # print(a[2].sum())
        except:
            temp.append(0)
        temp_Total_completion_time = np.sum(temp)
        
    return temp_Total_completion_time /(60*60*24)   #已換成天數單位

def get_between_days(start_sec, end_sec):
    '''獲得兩個日期之間的天數(前插法時會考慮)'''
    work_days = round(((end_sec - start_sec)/(24*60*60)),3)   #已換成天數單位
    # work_hours = int((end_sec - start_sec)/(60*60))
    # work_minutes = int((end_sec - start_sec)/(60))
    return work_days

def convergence_Fitness(nice_makespan,k):
    '''Elitist_Chromosomes_max_makespan的nice_Objective 呼叫近來當y座標,X座標是Iteration'''

    ypt = nice_makespan

    xpt = list(range(k))

    plt.title("Fitness")
    plt.xlabel("Iteration")
    plt.ylabel("Objective")
    # plt.scatter(xpt, ypt,s = 50, color = 'y')
    plt.plot(xpt, ypt,linestyle='solid', color = 'y')

    plt.show()
    return

def convergence_Fitness_subplot(Nb_Mc,Tardiness,TCT,k):
    # Import necessary libraries
    import matplotlib.pyplot as plt
    import numpy as np
    
    plt.rcParams['font.sans-serif'] = ['Microsoft JhengHei'] 
    plt.rcParams['axes.unicode_minus'] = False
    
    #Change the figure size
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
    plt.subplot(2,2, 4).title.set_text("總換模時間")
    # plt.xlabel("品牌名稱") 
    # plt.ylabel("文章數量")
    plt.subplot(2,2, 4)
    plt.plot(x, y3, '-.y', linewidth=3)
    
    
    
    plt.show()
    return
def Number_Of_Consecutive_Iterations_Occurs(a):
    b = a.iloc[-1]			# b的值为要求的数字
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
        
    print("%s連續出現%s次"%(np.round(b,4),t))
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
#%%
'''選擇機制'''
from pymoo.util.dominator import Dominator
import numpy as np
import pandas as pd
import copy

#2. use crowd distance to sort the solution,so each soluton have two new definition 
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

    
    #M = Dominator.calc_domination_matrix(F)
    M = Dominator.calc_domination_matrix(F.values)   #要再查一下再做什麼

    # calculate the dominance matrix
    n = M.shape[0]

    fronts = []

    if n == 0:
        return fronts

    # final rank that will be returned
    n_ranked = 0
    ranked = np.zeros(n, dtype=int)   #每一條染色體給ranked

    # for each individual a list of all individuals that are dominated by this one  #每個個體可以凌掠的解集合
    is_dominating = [[] for _ in range(n)]  #每一個p 的  Sp空集合

    # storage for the number of solutions dominated this one
    n_dominated = np.zeros(n)  #每個個體被別人凌掠的個數  就是np

    current_front = []  #目前的前緣線的解 (裡面會放個體)

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

        if n_dominated[i] == 0:  #若都沒有被別人凌掠就是代表第一名
            current_front.append(i)
            ranked[i] = 1.0  #單純標記
            n_ranked += 1    #有幾次front  = 0 的解

    # append the first front to the current front
    fronts.append(current_front)

    # while not all solutions are assigned to a pareto front
    while n_ranked < n:

        next_front = []  #新的front集合用空集合

        # for each individual in the current front
        for i in current_front:

            # all solutions that are dominated by this individuals
            for j in is_dominating[i]:
                n_dominated[j] -= 1
                if n_dominated[j] == 0:
                    next_front.append(j)
                    ranked[j] = 1.0  #單純標記
                    n_ranked += 1

        fronts.append(next_front)
        current_front = next_front

    return fronts

##2. use crowd distance to sort the solution,so each soluton have two new definition 
## caculate crowd distance
## filter_out_duplicates 過濾掉重複項
def cdist(A, B, **kwargs):
    import scipy
    return scipy.spatial.distance.cdist(A.astype(float), B.astype(float), **kwargs)  #計算兩個輸入集合中每對之間的距離

def find_duplicates(X, epsilon=1e-16): #(duplicates = 重複)
    # calculate the distance matrix from each point to another
    D = cdist(X, X)

    # set the diagonal to infinity  #np.triu_indices(len(X)) 可用成上三角形(包含自己跟自己 例如 A-> A 距離 0 )
    D[np.triu_indices(len(X))] = np.inf  

    # set as duplicate if a point is really close to this one
    is_duplicate = np.any(D <= epsilon, axis=1)

    return is_duplicate

def calc_crowding_distance(F,Weights ,filter_out_duplicates=True):
    '''跟網頁的一樣'''
    n_points, n_obj = F.shape  #(表示目標值的形狀 例如  幾條染色體(每個點代表一個個體) * 目標式數量)

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
        #_F = F[is_unique]
        _F = F.values[is_unique].astype(float)
        
        
        # sort each column and get index  #每一個目標去由小排到大，並回傳index
        I = np.argsort(_F, axis=0,kind='mergesort')
        
        

        # sort the objective space values for the whole matrix  #這個學起來!! 好用
        _F = _F[I, np.arange(n_obj)]     #可以問看看概念

        # calculate the distance from each point to the last and next    #要再想一下為什麼要用這個!!
        dist = np.row_stack([_F, np.full(n_obj, np.inf)]) - np.row_stack([np.full(n_obj, -np.inf), _F])

        # calculate the norm for each objective - set to NaN if all values are equal
        norm = np.max(_F, axis=0) - np.min(_F, axis=0)  #'''f1max - f1min 跟 f2max -f2min'''
        norm[norm == 0] = np.nan

        # prepare the distance to last and next vectors
        dist_to_last, dist_to_next = dist, np.copy(dist)
        dist_to_last, dist_to_next = dist_to_last[:-1] / norm, dist_to_next[1:] / norm
        
        
        # _________________20220216 新增 權重於 NSGAII_____________________________        
        for i in range(n_obj):
            dist_to_last[:,i] = dist_to_last[:,i]*Weights[i]
            dist_to_next[:,i] = dist_to_next[:,i]*Weights[i]

        # if we divide by zero because all values in one columns are equal replace by none
        dist_to_last[np.isnan(dist_to_last)] = 0.0
        dist_to_next[np.isnan(dist_to_next)] = 0.0

        # sum up the distance to next and last and norm by objectives - also reorder from sorted list
        J = np.argsort(I, axis=0)   #不太懂幹嘛在argsort 一次(轉回來??)
        _cd = np.sum(dist_to_last[J, np.arange(n_obj)] + dist_to_next[J, np.arange(n_obj)], axis=1) / n_obj

        # save the final vector which sets the crowding distance for duplicates to zero to be eliminated
        crowding = np.zeros(n_points)
        crowding[is_unique] = _cd

    # crowding[np.isinf(crowding)] = 1e+14
    return crowding

        
#3. Cacaluate Paroto_Optimial Front
def ValuesCal(Total_value,pop_size):
    temp_Paroto_Optimial = fast_non_dominated_sort(Total_value);  #計算Front線
    front = np.zeros(pop_size*2,dtype=int);      #再檢查一下pop_size 是否要*2
    for i in range(len(temp_Paroto_Optimial)):  #將front(前緣線)把值填進去
        for j in temp_Paroto_Optimial[i] :
            front[j] = i+1;         '''故意加1，用意是可以方便讓我們看前緣線是1開始 不是從0開始!!!'''
    fronts = copy.copy(front)
    Total_value["Front_value"] = pd.DataFrame(front);     #會對應front值
    crowding_of_front = np.zeros(pop_size*2);
    
    for k, fronts in enumerate(temp_Paroto_Optimial):  #先按照 0,1,2,... 前緣線，同一條前緣線就去計算crowding_distance
    # calculate the crowding distance of the front
        crowding_of_front[fronts] = calc_crowding_distance(Total_value.iloc[fronts, :-1],Weights)  #不要包含Front_value值 (注意式iloc 還是 loc)
# =============================================================================
#     crowding_of_front[np.isinf(crowding_of_front)] = -1      #遇到inf 就直接等於-1 就是代表要直接保留住 (嘗試改成99999999999)
# =============================================================================
##3.1 Compare there front and distenace values  
    Total_value["Crowding_Distance_Value"]= pd.DataFrame(crowding_of_front)    
    return Total_value

#4 Cacaluate Rank 
def RankforNSGAII(Total_value):
    cols = ['Front_value','Crowding_Distance_Value']; #設定要取行名
    # tups = Total_value[cols].sort_values(cols,ascending = [True,False]).apply(tuple,1) #排名 (由大到小)
    Total_valu = Total_value[cols].sort_values(cols,ascending = [True,False])
    # f,i = pd.factorize(tups)
    # factorized = pd.Series(f+1,tups.index).rank(ascending = False , method = 'min')
    return Total_value
    
def Forward_Insertion_Method_Improved(arg1, Zip_OrderNumber_Stime_Endtime, temp_Select_Fitness_Df1):
    # %%改善前插法
    '''處理子群體(Offstring)'''
    '''要做檢查!!'''
    '''將Zip_OrderNumber_Stime_Endtime 的 每個Job 訂單index ArrivalTime、PT、產品料號拆開，為了後續前插法'''
    # pdb.set_trace()
    for i in range(arg1, len(temp_Select_Fitness_Df1), nb_threads):
        print(i)
        # pdb.set_trace()
        best_temp_Mc_Corres_To_Job_Index, best_temp_Mc_Corres_To_Job_ArrivalTime, best_temp_Mc_Corres_To_Job_Pt, best_temp_Mc_PPNumber = Thebest_Machine_Corres_To_Job_ArrivalTime_Pt_PPN(
            Zip_OrderNumber_Stime_Endtime[i])

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
                            Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job][0]]  # 下一個的產品料號

                        if Setup_OrderName_PPN != Setup_RightOrderName_PPN:
                            '''判斷目前產品料號與右邊(下一個)產品料號是否不同，若不同則需加上換模時間(秒數)'''
                            SetupRight_PPN_Time = Mould_dict[Setup_OrderName_PPN + Setup_RightOrderName_PPN] * 60
                        else:
                            SetupRight_PPN_Time = 0

                        insert_JobPt_and_setup = (insert_Job_Pt_second + SetupRight_PPN_Time) / (24 * 60 * 60)
                        insert_Job_Arrival = re_best_Zip[re_arrange_Mc][inspection_Job][1]

                        '''最晚開始時間，若超過最晚開始一定會延遲(關鍵)'''
                        Latest_start_time = Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job][
                                                1] - SetupRight_PPN_Time - insert_Job_Pt_second
                        if Latest_start_time < 0:
                            '''若最晚開始時間 < 0 代表可立即開始'''
                            Latest_start_time = 0
                        '''判斷get_between_days(Machine_time_stamp_First , Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job][1]) 此間格天數 是否 大於 要插入的Job 加工天數 而且要符合 插入的 ArrivalTime 必須要 <= 最晚開始時間， 這樣才可做前插法'''
                        if get_between_days(Machine_time_stamp_First,
                                            Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job][1]) >= round(
                                insert_JobPt_and_setup, 3) and (insert_Job_Arrival <= Latest_start_time):
                            insert_Job_Start = max(insert_Job_Arrival, Machine_time_stamp_First)
                            insert_Job_End = insert_Job_Start + insert_Job_Pt_second
                            Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc].pop(inspection_Job)
                            # Machine_Corresponds_To_Job_Module[0][re_arrange_Mc].pop(inspection_Job)

                            Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc].insert(re_arrange_Job, [
                                re_best_Zip[re_arrange_Mc][inspection_Job][0], insert_Job_Start, insert_Job_End])

                            # 下面五列程式碼是更正前插法的順序
                            temp = copy.deepcopy(re_best_Zip[re_arrange_Mc][inspection_Job])
                            re_best_Zip[re_arrange_Mc].pop(inspection_Job)
                            re_best_Zip[re_arrange_Mc].insert(re_arrange_Job, temp)

                            best_temp_Mc_PPNumber[re_arrange_Mc].pop(inspection_Job)
                            best_temp_Mc_PPNumber[re_arrange_Mc].insert(re_arrange_Job, Order[
                                Order.index == Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job][0]][
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
                            Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job - 1][0]]

                        '''計算目前料號 與 上一個料號換模時間'''
                        if Setup_OrderName_PPN != Setup_LeftOrderName_PPN:
                            SetupLeft_PPN_Time = Mould_dict[Setup_OrderName_PPN + Setup_LeftOrderName_PPN] * 60
                        else:
                            SetupLeft_PPN_Time = 0

                        '''目前料號 與 下一個料號名稱'''
                        Setup_RightOrderName_PPN = Product_Part_Number[
                            Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job][0]]

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

                        Latest_start_time = Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job][
                                                1] - SetupRight_PPN_Time - insert_Job_Pt_second
                        if Latest_start_time < 0:
                            Latest_start_time = 0
                        if get_between_days(Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job - 1][2],
                                            Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job][1]) >= round(
                                insert_JobPt_and_setup, 3) and (insert_Job_Arrival <= Latest_start_time):
                            '''判斷要前插的空間是否夠插入，且插入的ArrvalTime必須<=最晚開始時間，否則無法插入'''
                            insert_Job_Start = max(insert_Job_Arrival,
                                                   Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job - 1][
                                                       2] + SetupLeft_PPN_Time)
                            insert_Job_End = insert_Job_Start + insert_Job_Pt_second

                            Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc].pop(inspection_Job)

                            Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc].insert(re_arrange_Job, [
                                re_best_Zip[re_arrange_Mc][inspection_Job][0], insert_Job_Start, insert_Job_End])

                            temp = copy.deepcopy(re_best_Zip[re_arrange_Mc][inspection_Job])
                            re_best_Zip[re_arrange_Mc].pop(inspection_Job)
                            re_best_Zip[re_arrange_Mc].insert(re_arrange_Job, temp)

                            best_temp_Mc_PPNumber[re_arrange_Mc].pop(inspection_Job)
                            best_temp_Mc_PPNumber[re_arrange_Mc].insert(re_arrange_Job, Order[
                                Order.index == Zip_OrderNumber_Stime_Endtime[i][re_arrange_Mc][re_arrange_Job][0]][
                                "產品料號"].values[0])

                            break
                        else:
                            continue

        '''要增加修改'''
        '''修正開始時間結束時間'''

        Setup_count = 0
        Setup_times = 0  # 把這個取消掉 Sum_Tardiness 就跑出來  超怪的!!
        for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime[i]) - 1):
            Deep_re_best_Zip_DF = pd.DataFrame(deep_Df_re_best_Zip[Nb_Of_Machine_index],
                                               columns=["訂單編號", "最早到達時間", "加工時間"])
            for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index])):
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
                    Job_Start = max(Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index][Nb_Of_Job_index - 1][2] + Setup,
                                    Job_ArrivalTime)
                    Job_End = Job_Start + Pt
                    '''從第二個Job開始的開始時間'''
                    Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index][Nb_Of_Job_index][1] = Job_Start
                    Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index][Nb_Of_Job_index][2] = Job_End
                    # print(Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index][Nb_Of_Job_index][1] ,Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index][Nb_Of_Job_index][2] )

        '''0520前插法後的新的 Total_completion_time(Flow time) 前插法後的ArrivalTime需再確認'''
        # =============================================================================
        # Total_completion_time = Thebest_List_to_Np(Zip_OrderNumber_Stime_Endtime[i])
        # =============================================================================
        # pdb.set_trace()
        Each_Total_EndTime = Thebest_List_to_Np(Zip_OrderNumber_Stime_Endtime[i])  # ˊ轉換成矩陣模式，方便做下一步np.sum的運算
        Each_temp_Mc_Co_To_Job_Arrival_matrix = List_to_Np(best_temp_Mc_Corres_To_Job_ArrivalTime)
        Each_Total_ArrivalTime = np.sum(Each_temp_Mc_Co_To_Job_Arrival_matrix)
        Total_completion_time = Each_Total_EndTime - Each_Total_ArrivalTime

        '''處理目標，each_Total_completion_Tardiness_Time = Total_completion_time + Tardiness * w'''
        Sum_Tardiness = 0
        Sum_Tardiness_cost = 0
        for Nb_Of_Machine_index in range(len(Zip_OrderNumber_Stime_Endtime[i]) - 1):
            Machine_struct_time = time.strptime(Machine_Start.at[Nb_Of_Machine_index, "機台開始時間"],
                                                "%Y/%m/%d %H:%M")  # 轉成時間元組
            '''下面Machine_time_stamp是每台機台的起始時間，變動會隨著停機開始時間csv檔'''
            Machine_time_stamp = int(time.mktime(Machine_struct_time))

            for Nb_Of_Job_index in range(len(Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index])):
                Completion_time = Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index][Nb_Of_Job_index][2]

                # Due_Date_time = int(time.mktime(time.strptime(Order.at[Order[Order["訂單編號"] ==  Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index][Nb_Of_Job_index][0]].index[0],'交期'], "%Y/%m/%d %H:%M")))-Machine_time_stamp
                Due_Date_time = int(time.mktime(time.strptime(Order.at[Order[Order.index ==
                                                                             Zip_OrderNumber_Stime_Endtime[i][
                                                                                 Nb_Of_Machine_index][Nb_Of_Job_index][
                                                                                 0]].index[0], '交期'],
                                                              "%Y/%m/%d %H:%M"))) - Machine_time_stamp

                if Completion_time - Due_Date_time > 0:
                    Sum_Tardiness += (Completion_time - Due_Date_time)
                    '''計算此單延誤損失成本'''
                    # Sum_Tardiness_cost += Order.at[Zip_OrderNumber_Stime_Endtime[i][Nb_Of_Machine_index][Nb_Of_Job_index][0],"單張延誤成本"]

        Sum_Tardiness = Sum_Tardiness / (60 * 60 * 24)
        Setup_times = Setup_times / (60 * 60 * 24)

        # =============================================================================
        #     '''0531下面目標式用Setup_times'''
        #     each_Single_Target = Machine_index * Penalty_value_weight_Machine + Sum_Tardiness * Penalty_value_weight_TT + Setup_times * Penalty_value_weight_Setup
        # =============================================================================
        temp_Select_Fitness_Df1[i][2], temp_Select_Fitness_Df1[i][3] = Sum_Tardiness, Setup_times

        # Zip_OrderNumber_Stime_Endtime[i][-1] =  np.array([Improve_Interation,Machine_index,Sum_Tardiness,Setup_times,each_Single_Target])
        Zip_OrderNumber_Stime_Endtime[i][-1] = temp_Select_Fitness_Df1[i].tolist()
    return Zip_OrderNumber_Stime_Endtime, temp_Select_Fitness_Df1

if __name__ =='__main__':
    t= []
    for i in range(nb_threads):
        t.append(threading.Thread(target = Forward_Insertion_Method_Improved, args=(i, List_OrderNumber_Stime_Endtime,temp_Select_Fitness_Df1)))
    
    print('>>>>>>>>>>>>>>>>>>>>>>>>> thread_create')
    st = time.time()
    for _t in t:
        _t.start()
    # pdb.set_trace()
    time.sleep(1)
    print('<<<<<<<<<<<<<<<<<<<<<<<<< thread_join')
    for _t in t:
        _t.join()
    ent = time.time()
    print("前插法平行",ent-st)
    k = 0
    Fitness_result = DN.Normalization(temp_Select_Fitness_Df1[:,1:],Penalty_Weights).calc()     #計算fitness
    temp_Select_Fitness_Df1  = np.concatenate((temp_Select_Fitness_Df1, Fitness_result), axis=1)  #做結合
    
    new_decode  = []  #要放新一代菁英跟 輪盤的解碼
    new_moudule = []  #要放解碼後的模具
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
    if k ==0:
        temp_Select_Fitness_Df1 = pd.DataFrame(temp_Select_Fitness_Df1,columns = ['Index',"機台數","總延誤時間","總換模時間","Fitness"])
        Fix_temp_Select_Fitness_Df1 = 0
    # c = Forward_Insertion_Method_Improved(List_OrderNumber_Stime_Endtime,temp_Select_Fitness_Df1)
    # w = threading.Thread(target=Forward_Insertion_Method_Improved)
