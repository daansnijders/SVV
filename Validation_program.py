# -*- coding: utf-8 -*-
import scipy
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

node = []
element_set = []
n_set = []
j = 0

with open("./data/F100-19.inp") as g:
    for oneline in g:
        if oneline.startswith('*') == False:
            if j < 3244:
                #print("First list node", j)
                oneline = oneline.strip('\n')
                line_cont = [float(y) for y in oneline.split(',')]
                node.append(line_cont)
                
            elif j > 3243 and j < (3243+3251):
                #print("2nd list point ",j)
                oneline = oneline.strip('\n')
                line_cont = [float(y) for y in oneline.split(',')]
                element_set.append(line_cont)
             
            elif j > 7020 and j < (7038):
                print (oneline)
                #oneline = oneline.strip('\n')
                line_cont = [float(y) for y in oneline.split(',')]
                n_set.append(line_cont)                
        j = j + 1


#
#with open("./data/F100-19.inp") as g:
#    for oneline in g:
#        for j in range(0, 3235):
#            if oneline.startswith('*') == False:
#                oneline = oneline.strip('\n')
#                line_cont = [float(y) for y in oneline.split(',')]
#                node.append(line_cont)
#        for j in range(3236, 6485):
#            if oneline.startswith('*') == False:
#                oneline = oneline.strip('\n')
#                line_cont = [float(y) for y in oneline.split(',')]
#                element_set.append(line_cont)
#        #for j in range(7009,7031):
#           # if oneline.startswith('*') == False:
#               # print("SOMEHTIN")
#            
            