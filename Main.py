# -*- coding: utf-8 -*-
import math
import Definitions
import numpy as np


##### Variable Definitions #####

q = -3860 #load distribution, + upwards
P2 = 49200 #Actuator II force in negative direction

E = 73.1e+09 #E-modulus

#Convergence Variables

ndis = 100 #No. of Sections per section discretized
spread = 0.001 #dx for Jacobian Convergence


#Aileron Geometry

l1 = 0.125
l2 = 0.498
l3 = 1.494
l4 = 1.611
xa = 0.245

Ca = 0.505
ha = 0.161


#Beam Deflection and Mx Distribution Variables


inittwist = 0 #Twist of rib C in degrees (counterclockwise upwards)

theta = np.ones(ndis*6+1)*inittwist/180*np.pi #theta[4*n] is actuator 2 and theta[2*n] is actuator 1


#Y deflections of hinges

d1 = 0.389
d2 = 0
d3 = 1.245


#Boom Area Variables

n = 11 #No. of stringers
list_length = n+3 #No. of idealised booms INCLUDING two added AND Ghost nodepos
tskin = 0.0011
tspar = 0.0024
t_stiff = 0.0012
h_stiff = 0.013
w_stiff = 0.017

zsc = 0 #Shear Center Location (Required but left at 0)


#####General Code######

#Initial Centroid

spacing, Cr, alpharad, Ct = Definitions.boom_spacing(ha, Ca, n)

nodepos, arc, dist = Definitions.boom_location(spacing, Cr, alpharad, list_length, ha)

area_stiff = Definitions.area_stiff(t_stiff, h_stiff, w_stiff)

ycg, zcg = Definitions.centroid_nonidealized(tskin, ha, Ca, Ct, tspar, nodepos, area_stiff)


#Initial Moment of Inertia

I = np.array([[np.ones(ndis*6+1)*9.434e-05, np.ones(ndis*6+1)*0],[np.ones(ndis*6+1)*0, np.ones(ndis*6+1)*1.252e-05]])
Izz = I[1][1]




#Initial reaction forces

r1, rz1, rx2, r2, rz2, r3, rz3, P1 = Definitions.ReactionForces(theta[4*n],P2,-q,Ca,ha,E,Izz[3*n],l1,l2,l3,l4,xa,d1,d3)
P1 = -P1



##spacing, Cr, alpharad, Ct = Definitions.boom_spacing(ha, Ca, n)
##
##alphadeg = math.degrees(alpharad)
##
##nodepos, arc, dist = Definitions.boom_location(spacing, Cr, alpharad, list_length, ha)
##
##area_stiff = Definitions.area_stiff(t_stiff, h_stiff, w_stiff)
##
##
##
##
##ycg, zcg = Definitions.centroid_nonidealized(tskin, ha, Ca, Ct, tspar, nodepos, area_stiff)
##
###print (boom_area_incl_skin)
##
##
##
##
##My = 0.1
##Mz = 0.1
##
##
##"""Sybren calculate centroid of non idealized structure"""
##""""Sybren calculate MOI 30 deg and 0 deg"""
##
##
##boom_area = Definitions.boom_area_updater(tskin, spacing, Mz, My, Izz, Iyy, area_stiff, zcg, nodepos, dist, arc, tspar, ha)
##Ixx, Iyy, Izz = Definitions.boom_inertia(list_length, nodepos, boom_area)
##
##
###Moment of inertia assumed as matrix with [[Iyy, Izy];[Izy, Izz]] for beam bending
##
##
##
######## Beam Bending ###################################################################################33
##
##
##
##
##
##
##
###moment of inertia assumed as matrix with [[Iyy, Izy];[Izy, Izz]]
##I = np.array([[np.ones(n*6+1)*8.26559e-05, np.ones(n*6+1)*-3.42075e-05],[np.ones(n*6+1)*-3.42075e-05, np.ones(n*6+1)*2.92041e-05]])
##
##
##
##v2, u2, xt, r1, r2, r3, My, Mz, rz1, rz2, rz3, P1 = Definitions.bendingconvergence(q,ndis,l1,l2,l3,l4,E,I,d1,d2,d3,P2,xa,Ca,ha,theta2,theta1,spread)












