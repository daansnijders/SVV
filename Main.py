# -*- coding: utf-8 -*-
import math
import Definitions


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

q = 3860
n = 200

l1 = 0.125
l2 = 0.498
l3 = 1.494
l4 = 1.611
xa = 0.245

Ca = 0.245
theta2 = 30*np.pi/180
ha = 0.161

d1 = 0.389
d2 = 0
d3 = 1.245

E = 72000000


rz3 = 50000
P2 = 49200
r3 = 150000

spread = 0.0001

v2, u2, xt, r1, r2, r3, My, Mz, rz1, rz2, rz3, P1 = bendingconvergence(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,spread)












