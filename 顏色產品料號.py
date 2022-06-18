# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 16:45:01 2021

@author: ZhiXiang
"""

import xlrd

def read_excel_data():
    filename = '數字.xlsx'
    data = xlrd.open_workbook(filename)
    table = data.sheet_by_name('工作表1')
    row_num = table.nrows  # 行数
    # col_num = table.ncols  # 列数
    datas = dict([]) # 这步也要转字典类型
    for i in range(row_num):
        colorss = dict([table.row_values(i)]) # 这一步就要给它转字典类型，不然update没法使用
        datas.update(colorss)
    print(datas)
    return datas
if __name__ == "__main__":
    colorss = read_excel_data()