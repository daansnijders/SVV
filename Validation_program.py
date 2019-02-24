# -*- coding: utf-8 -*-
import scipy
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

node = []
j = 0
with open("./data/F100-19.inp") as g:
    for oneline in g:
        if oneline.startswith('*') == False:
            if j < 3244:
                j =+ 1
                line_cont = [str(y) for y in oneline.split(',')]
                node.append(line_cont)
          
            
        else:
            j =+ 1
