# -*- coding: utf-8 -*-
import math
import Definitions
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
from mpl_toolkits.mplot3d import Axes3D

################################# Variable Definition ########################################################################

q = -3860 #load distribution, + upwards
P2 = 49200 #Actuator II force in negative direction

E = 73.1e+09 #E-modulus
G = 28e+09

###Convergence Variables ###

ndis = 100 #No. of Sections per section discretized
spread = 0.00001 #dx for Jacobian Convergence
tol = 0.0005
c = 0.01


###Aileron Geometry ###

l1 = 0.125
l2 = 0.498
l3 = 1.494
l4 = 1.611
xa = 0.245

Ca = 0.505
ha = 0.161


###Beam Deflection and Mx Distribution Variables ###


inittwist = 30 #Twist of rib C in degrees (counterclockwise upwards)

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

#Initial Moment of Inertia - Working

I, Ilocal = Definitions.ExactMOIdiscretisation(q,ndis,l1,l2,l3,l4,tskin,tspar,t_stiff,w_stiff,h_stiff,zcg,n,spacing,nodepos,xa,Ca,ha,theta,zsc)


#Initial reaction forces

v2, u2, xt, r1, r2, r3, Vy, Vz, My, Mz, rz1, rz2, rz3, P1 = Definitions.bendingconvergence(q,ndis,l1,l2,l3,l4,E,I,d1,d2,d3,P2,xa,Ca,ha,theta,spread)

#Initial Twist Calculation

Mx,xt = Definitions.torque(q,ndis,l1,l2,l3,l4,P1,P2,xa,Ca,ha,theta,zsc)

theta, rate_twist_lst,xt = Definitions.overalltwist(Mx,A1,A2,arc,Cr,ha,xa,G,tskin,l1,l2,l3,l4,ndis,inittwist)


##### Iteration ######



iteration = 0

error = 1

while error > tol:

    iteration += 1
    print('Iteration no. ' + str(iteration)+'\n')

    #Geometrical Update

    if iteration >1:
        boom_areaold = boom_area
    
    boom_area, twist_rate, qrib_1, qrib_2 = Definitions.ratetwistandshearflowdiscretisation(tskin, tspar, spacing, l1,l2,l3,l4,xa, Mz, My, Mx, Vy, Vz, Ilocal, area_stiff, zcg, nodepos, dist, arc, Ca, ha, G, theta, alpharad, ndis)

    if iteration >1:
        boom_area = boom_areaold + (boom_area-boom_areaold)*c

    #Initial Moment of Inertia - Working

    I, Ilocal = Definitions.idealisedMOIdiscretisation(ndis,l1,l2,l3,l4,xa,list_length, nodepos, boom_area, theta)
    
    #Beam Deflection Convergence

    r1old, r2old, r3old, rz1old, rz2old, rz3old, P1old = r1, r2,r3,rz1,rz2,rz3,P1

    v2, u2, xt, r1, r2, r3, Vy, Vz, My, Mz, rz1, rz2, rz3, P1 = Definitions.bendingconvergence(q,ndis,l1,l2,l3,l4,E,I,d1,d2,d3,P2,xa,Ca,ha,theta,spread)

    #Initial Twist Calculation

    Mx,xt = Definitions.torque(q,ndis,l1,l2,l3,l4,P1,P2,xa,Ca,ha,theta,zsc)
    
    theta, rate_twist_lst,xt = Definitions.overalltwist(-Mx,A1,A2,arc,Cr,ha,xa,G,tskin,l1,l2,l3,l4,ndis,inittwist)

    #Print Results

    error = np.max([(abs(r1-r1old)/r1old),abs((r2-r2old)/r2old), abs((r3-r3old)/r3old),abs((rz1-rz1old)/rz1old), abs((rz2-rz2old)/rz2old), abs((rz3-rz3old)/rz3old), abs((P1-P1old)/P1old)])
       
    print('\n'+'Ry1 = ' , float(r1) ,' Ry2 = ', float(r2) , ' Ry3 = ', float(r3) , '\r'+'\n'+' Rz1 = ', float(rz1) , ' Rz2 = ', float(rz2) ,' Rz3 = ',float(rz3), '\n', 'P1 = ', P1, '\n', 'Error = ', float(error))


#Update q-rib1 and qrib2 and find nodepositions

boom_area, twist_rate, qrib_1, qrib_2 = Definitions.ratetwistandshearflowdiscretisation(tskin, tspar, spacing, l1,l2,l3,l4,xa, Mz, My, Mx, Vy, Vz, Ilocal, area_stiff, zcg, nodepos, dist, arc, Ca, ha, G, theta, alpharad, ndis)
nodepos2, rot = Definitions.offset(zcg, theta, nodepos, v2, u2, xt)


#Print Forces

print('\n'+'Ry1 = ' , float(r1[0]) ,' Ry2 = ', float(r2[0]) , ' Ry3 = ', float(r3[0]) , '\r'+'\n'+' Rz1 = ', float(rz1[0]) , ' Rz2 = ', float(rz2[0]) ,' Rz3 = ',float(rz3[0]), '\n', 'P1 = ', P1, 'P2 = ',P2)


#Plot node positions

plt.subplot(2,2,1)
plt.plot(nodepos2[0][:,0],nodepos2[0][:,1])
plt.subplot(2,2,2)
plt.plot(nodepos2[0][:,0], nodepos2[0][:,2])

plt.subplot(2,2,3)
plt.plot(nodepos2[6][:,0],nodepos2[6][:,1])
plt.subplot(2,2,4)
plt.plot(nodepos2[6][:,0], nodepos2[6][:,2])

plt.show()



















