# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 18:05:27 2022

@author: ZhiXiang
"""

import numpy as np
import pandas as pd
import warnings
from scipy.stats import rankdata

#tp.Topsis(TOPSIS_value, Weights,criterias).calc()
class Normalization():
    evaluation_matrix = np.array([])  # Matrix
    weighted_normalized = np.array([])  # Weight matrix
    normalized_decision = np.array([])  # Normalisation matrix
    M = 0  # Number of rows
    N = 0  # Number of columns

    '''
	Create an evaluation matrix consisting of m alternatives and n criteria,
	with the intersection of each alternative and criteria given as {\displaystyle x_{ij}}x_{ij},
	we therefore have a matrix {\displaystyle (x_{ij})_{m\times n}}(x_{{ij}})_{{m\times n}}.
	'''

    def __init__(self, evaluation_matrix, weight_matrix):
        # M×N matrix
        self.evaluation_matrix = np.array(evaluation_matrix, dtype="float")   #建構評分矩陣

        # M alternatives (options)
        self.row_size = len(self.evaluation_matrix)

        # N attributes/criteria
        self.column_size = len(self.evaluation_matrix[0])

        # N size weight matrix
        self.weight_matrix = np.array(weight_matrix, dtype="float")

    '''
	# Step 1
	Determine the worst alternative {\displaystyle (A_{w})}(A_{w}) and the best alternative {\displaystyle (A_{b})}(A_{b}):
	'''
    '''將 Wi / max(Fi) 計算出來'''
    def step_1(self):
        self.worst_alternatives = np.zeros(self.column_size)
        self.best_alternatives = np.zeros(self.column_size)
        for i in range(self.column_size):

            self.best_alternatives[i] = self.weight_matrix[i] /max(self.evaluation_matrix[:, i])   #理想解乃是望大的屬性上達最大值
        
        
    '''
	# Step 2
	The matrix {\displaystyle (x_{ij})_{m\times n}}(x_{{ij}})_{{m\times n}} is then normalised to form the matrix
	'''
    '''每個原始目標 乘上權重'''
    def step_2(self):
        # normalized scores
        self.normalized_decision = np.copy(self.evaluation_matrix)
        self.sum_Fitness = np.zeros((self.row_size,1))
        for i in range(self.row_size):
            '''每一個屬性平方加總(計算邏輯要想一下很讚)'''
            for j in range(self.column_size):
                '''已經將Wi/max(Fi) * Fi 算成self.normalized_decision'''
                self.normalized_decision[i][j] *= self.best_alternatives[j]
            
            for j in range(self.column_size):
                '''將 self.normalized_decision 每個方案的目標正歸化後加起來'''
                self.sum_Fitness[i] +=  self.normalized_decision[i][j]
            
            # print(self.sum_Fitness)
                
                # sqrd_sum[j] += self.evaluation_matrix[i, j]**2
        # for i in range(self.row_size):
        #     for j in range(self.column_size):
        #         self.normalized_decision[i,
        #                                  j] = self.evaluation_matrix[i, j]/(sqrd_sum[j]**0.5)
        
        self.normalized_decision = pd.DataFrame(self.normalized_decision).fillna(0);
        # self.sum_Fitness= pd.DataFrame(self.sum_Fitness).fillna(0);
        

    def calc(self):

        self.step_1()

        self.step_2()
     
        # self.sum_Fitness.columns = ['Fitness']
        # self.sum_Fitness.fillna(0) # 若結果一直收斂，計算的直將會變成 nan 因此計算不出距離需為 0 20210605  
        #columns = ['Score','Rank']
        # self.sum_Fitness['Fitness'] = self.sum_Fitness['Fitness'].astype(float)
        return self.sum_Fitness
