# -*- coding: utf-8 -*-
import math
import Definitions
import numpy as np


ha = 0.161 #height of thickest part of airfoil
Ca = 0.505 #dist from trailing to leading edge
n = 11 #number of stringers
list_length = n+3 #number of idealised booms INCLUDING two added AND Ghost nodepos
tskin = 0.0011
tspar = 0.0024
t_stiff = 0.0012
h_stiff = 0.013
w_stiff = 0.017



spacing, Cr, alpharad, Ct = Definitions.boom_spacing(ha, Ca, n)

alphadeg = math.degrees(alpharad)

nodepos, arc, dist = Definitions.boom_location(spacing, Cr, alpharad, list_length, ha)

area_stiff = Definitions.area_stiff(t_stiff, h_stiff, w_stiff)




ycg, zcg = Definitions.centroid_nonidealized(tskin, ha, Ca, Ct, tspar, nodepos, area_stiff)

#print (boom_area_incl_skin)




My = 0.1
Mz = 0.1


"""Sybren calculate centroid of non idealized structure"""
""""Sybren calculate MOI 30 deg and 0 deg"""


boom_area = Definitions.boom_area_updater(tskin, spacing, Mz, My, Izz, Iyy, area_stiff, zcg, nodepos, dist, arc, tspar, ha)
Ixx, Iyy, Izz = Definitions.boom_inertia(list_length, nodepos, boom_area)


#Moment of inertia assumed as matrix with [[Iyy, Izy];[Izy, Izz]] for beam bending
I = np.array([[np.ones(n*4+1)*4.38537e-04, np.ones(n*4+1)*6.475e-04],[np.ones(n*4+1)*6.475e-04, np.ones(n*4+1)*6.0755e-04]])


###### Beam Bending ###################################################################################33

q = -4530
ndis = 100

l1 = 0.153
l2 = 1.281
l3 = 2.681
l4 = 2.771
xa = 0.28

Ca = 0.547
ha = 0.225

#Theta2 is actuator 2 and theta is actuator 1
theta = np.ones(n*6+1)*26/180*np.pi
theta1 = theta[2*n]
theta2 = theta[4*n]
zsc = 0


d1 = 0.01103
d2 = 0
d3 = 0.01642


E = 73.1e+09
#moment of inertia assumed as matrix with [[Iyy, Izy];[Izy, Izz]]
I = np.array([[np.ones(n*6+1)*8.26559e-05, np.ones(n*6+1)*-3.42075e-05],[np.ones(n*6+1)*-3.42075e-05, np.ones(n*6+1)*2.92041e-05]])

P2 = 91700

spread = 0.001

v2, u2, xt, r1, r2, r3, My, Mz, rz1, rz2, rz3, P1 = Definitions.bendingconvergence(q,ndis,l1,l2,l3,l4,E,I,d1,d2,d3,P2,xa,Ca,ha,theta2,theta1,spread)












