# -*- coding: utf-8 -*-
import math
import Definitions
import numpy as np
import matplotlib.pyplot as plt
import sys
import time

################################# Variable Definition ########################################################################

q = -3860 #load distribution, + upwards
P2 = 49200 #Actuator II force in negative direction

E = 73.1e+09 #E-modulus
G = 28e+09

###Convergence Variables ###

ndis = 100 #No. of Sections per section discretized
spread = 0.001 #dx for Jacobian Convergence


###Aileron Geometry ###

l1 = 0.125
l2 = 0.498
l3 = 1.494
l4 = 1.611
xa = 0.245

Ca = 0.505
ha = 0.161


###Beam Deflection and Mx Distribution Variables ###


inittwist = 0 #Twist of rib C in degrees (counterclockwise upwards)

theta = np.ones(ndis*6+1)*inittwist/180*np.pi #theta[4*n] is actuator 2 and theta[2*n] is actuator 1


###Y deflections of hinges ###

d1 = 0.00389
d2 = 0
d3 = 0.01245


###Boom Area Variables ###

n = 11 #No. of stringers
list_length = n+3 #No. of idealised booms INCLUDING two added AND Ghost nodepos
tskin = 0.0011
tspar = 0.0024
t_stiff = 0.0012
h_stiff = 0.013
w_stiff = 0.017

A1=0.0101791529 # Area of cell 1 
A2=0.03417225 # Area of cell 2

zsc = 0 #Shear Center Location (Required but left at 0)


##################################### General Code #############################################################################



#Initial Centroid - Working

spacing, Cr, alpharad, Ct = Definitions.boom_spacing(ha, Ca, n)

nodepos, arc, dist = Definitions.boom_location(spacing, Cr, alpharad, list_length, ha)

area_stiff = Definitions.area_stiff(t_stiff, h_stiff, w_stiff)

ycg, zcg = Definitions.centroid_nonidealized(tskin, ha, Ca, Ct, tspar, nodepos, area_stiff)


#Plot Nodepos

Definitions.scatter(nodepos)

#Initial Moment of Inertia - Working

I = Definitions.ExactMOIdiscretisation(q,ndis,l1,l2,l3,l4,tskin,tspar,t_stiff,w_stiff,h_stiff,zcg,n,spacing,nodepos,xa,Ca,ha,theta,zsc)

zcentroidglobal, ycentroidglobal = Definitions.centroidglobal(zcg,theta)

#Initial reaction forces



##r1, rz1, rx2, r2, rz2, r3, rz3, P1 = Definitions.ReactionForces(theta[4*n],P2,-q,Ca,ha,E,I[1][1][3*ndis],l1,l2,l3,l4,xa,d1,d3)
##P1 = -P1

v2, u2, xt, r1, r2, r3, Vy, Vz, My, Mz, rz1, rz2, rz3, P1 = Definitions.bendingconvergence(q,ndis,l1,l2,l3,l4,E,I,d1,d2,d3,P2,xa,Ca,ha,theta,spread)





#Initial Twist Calculation

Mx,xt = Definitions.torque(q,ndis,l1,l2,l3,l4,P1,P2,xa,Ca,ha,theta,zsc)

theta, rate_twist_lst,xt = Definitions.overalltwist(Mx,A1,A2,arc,Cr,ha,xa,G,tskin,l1,l2,l3,l4,ndis,inittwist)


##### Iteration ######



iteration = 0


while iteration < 2:


    iteration += 1
    print('Iteration no. ' + str(iteration)+'\n')

    #Beam Deflection Convergence

    v2, u2, xt, r1, r2, r3, Vy, Vz, My, Mz, rz1, rz2, rz3, P1 = Definitions.bendingconvergence(q,ndis,l1,l2,l3,l4,E,I,d1,d2,d3,P2,xa,Ca,ha,theta,spread)

    #Initial Moment of Inertia - Working

    I = Definitions.ExactMOIdiscretisation(q,ndis,l1,l2,l3,l4,tskin,tspar,t_stiff,w_stiff,h_stiff,zcg,n,spacing,nodepos,xa,Ca,ha,theta,zsc)
    zcentroidglobal, ycentroidglobal = Definitions.centroidglobal(zcg,theta)

    #Initial Twist Calculation

    Mx,xt = Definitions.torque(q,ndis,l1,l2,l3,l4,P1,P2,xa,Ca,ha,theta,zsc)

    theta, rate_twist_lst,xt = Definitions.overalltwist(-Mx,A1,A2,arc,Cr,ha,xa,G,tskin,l1,l2,l3,l4,ndis,inittwist)

   
    print('\n'+'Ry1 = ' , float(r1[0]) ,' Ry2 = ', float(r2[0]) , ' Ry3 = ', float(r3[0]) , '\n\r'+' Rz1 = ', float(rz1[0]) , ' Rz2 = ', float(rz2[0]) ,' Rz3 = ',float(rz3[0]), '\r')



#Discretised Boom Area Twist rate and Qrib

boom_area, twist_rate, qrib_1, qrib_2 = Definitions.ratetwistandshearflowdiscretisation(tskin, tspar, spacing, l1,l2,l3,l4,xa, Mz, My, Mx, Vy, Vz, I, area_stiff, zcg, nodepos, dist, arc, Ca, ha, G, theta, alpharad, ndis)













