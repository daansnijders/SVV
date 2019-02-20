# -*- coding: utf-8 -*-
import math
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import Definitions
import numpy as np

ha = 0.161 #height of thickest part of airfoil
Ca = 0.505 #dist from trailing to leading edge
n = 11 #number of stringers
list_length = n+3 #number of idealised booms INCLUDING two added AND Ghost nodepos
tskin = 0.0011
t_stiff = 0.0012
h_stiff = 0.013
w_stiff = 0.017

spacing, Cr, alpharad = Definitions.boom_spacing(ha, Ca, n)
#nodepos = Definitions.boom_location(spacing, Cr, alpharad, list_length)

alphadeg = math.degrees(alpharad)

nodepos = Definitions.boom_location(spacing, Cr, alpharad, list_length, ha)
print (nodepos[0:14][2])

plt.plot(nodepos[0:14][2], nodepos[0:14][1])
plt.show()

