# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 17:15:09 2019

@author: Daan
"""

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



##r1, rz1, rx2, r2, rz2, r3, rz3, P1 = Definitions.ReactionForces(theta[4*n],P2,-q,Ca,ha,E,I[1][1][3*ndis],l1,l2,l3,l4,xa,d1,d3)
##P1 = -P1

v2, u2, xt, r1, r2, r3, Vy, Vz, My, Mz, rz1, rz2, rz3, P1 = Definitions.bendingconvergence(q,ndis,l1,l2,l3,l4,E,I,d1,d2,d3,P2,xa,Ca,ha,theta,spread)





#Initial Twist Calculation

Mx,xt = Definitions.torque(q,ndis,l1,l2,l3,l4,P1,P2,xa,Ca,ha,theta,zsc)


theta, rate_twist_lst,xt = Definitions.overalltwist(Mx,A1,A2,arc,Cr,ha,xa,G,tskin,l1,l2,l3,l4,ndis,inittwist)


##### Iteration ######



iteration = 0

while iteration < 20:

    iteration += 1
    print('Iteration no. ' + str(iteration)+'\n')

    #Geometrical Update
    
    boom_area, twist_rate, qrib_1, qrib_2 = Definitions.ratetwistandshearflowdiscretisation(tskin, tspar, spacing, l1,l2,l3,l4,xa, Mz, My, Mx, Vy, Vz, Ilocal, area_stiff, zcg, nodepos, dist, arc, Ca, ha, G, theta, alpharad, ndis)

    #Initial Moment of Inertia - Working

    I, Ilocal = Definitions.idealisedMOIdiscretisation(ndis,l1,l2,l3,l4,xa,list_length, nodepos, boom_area, theta)

    #Beam Deflection Convergence

    v2, u2, xt, r1, r2, r3, Vy, Vz, My, Mz, rz1, rz2, rz3, P1 = Definitions.bendingconvergence(q,ndis,l1,l2,l3,l4,E,I,d1,d2,d3,P2,xa,Ca,ha,theta,spread)

    #Initial Twist Calculation

    Mx,xt = Definitions.torque(q,ndis,l1,l2,l3,l4,P1,P2,xa,Ca,ha,theta,zsc)

    theta, xt =  Definitions.overalltwist2(twist_rate,xa,G,l1,l2,l3,l4,ndis,inittwist)
##    theta, rate_twist_lst,xt = Definitions.overalltwist(-Mx,A1,A2,arc,Cr,ha,xa,G,tskin,l1,l2,l3,l4,ndis,inittwist)

   
    print('\n'+'Ry1 = ' , float(r1[0]) ,' Ry2 = ', float(r2[0]) , ' Ry3 = ', float(r3[0]) , '\r'+'\n'+' Rz1 = ', float(rz1[0]) , ' Rz2 = ', float(rz2[0]) ,' Rz3 = ',float(rz3[0]), '\n', 'P1 = ', P1)


#Update q-rib1 and qrib2 and find nodepositions

boom_area, twist_rate, qrib_1, qrib_2 = Definitions.ratetwistandshearflowdiscretisation(tskin, tspar, spacing, l1,l2,l3,l4,xa, Mz, My, Mx, Vy, Vz, Ilocal, area_stiff, zcg, nodepos, dist, arc, Ca, ha, G, theta, alpharad, ndis)
nodepos2 = Definitions.offset(zcg, theta, nodepos, v2, u2, xt)


#Print Forces

print('\n'+'Ry1 = ' , float(r1[0]) ,' Ry2 = ', float(r2[0]) , ' Ry3 = ', float(r3[0]) , '\r'+'\n'+' Rz1 = ', float(rz1[0]) , ' Rz2 = ', float(rz2[0]) ,' Rz3 = ',float(rz3[0]), '\n', 'P1 = ', P1, 'P2 = ',P2)

import numpy as np
import math
tring_booms = [1,2,3,4,12,13,8,9,10,11]
circ_booms = [12,5,6,7,13]
alpharad=0.187
spacing=0.1015480896
list_length=14
ha=0.161
Cr=0.0805
#tring_qsum=[17,27,37,47,57,67,77,87,97,107]
#circ_qsum=[12,13,14,15,16]
tskin=0.0011
tsk=tskin
tspar=0.0024
Ca=0.505
ha=0.161
G=28000000000
theta=(math.pi)/6.0
Izz0=4.7189e-6
Iyy0=4.72552e-5

Izz=3.6621e-5
Iyy=1.35357e-5
Mx=1664.25

#Mz=-1.14161513e+04
#My=3.74664695e+03

Vz=49200  #value to be checked 
Vy=3860   # value to be checked 

stiff_area=3.5999999999999994e-05
zcg=-0.09993502900006332
dist=0.43206 # WHAT VALUES IS THIS ?
arc=0.2528982086



n=100
#T = np.ones(6*n+1)*3.76 # list containing the torque values at each spanwise location of the aileron 
A1=0.0101791529 # Area of cell 1 
A2=0.03417225 # Area of cell 2
arc=0.2528982086
l=0.43206
ha=0.161
G=28000000000
t=0.0011
#n=400.
l1=0.125
l2=0.498    
l3=1.494
l4=1.611
xa=0.245

q=3860.
#P1=53585

q=-3860
P1=57561.755

P2=49200
Ca=0.505
#theta=np.ones(6*n+1)*30*(np.pi/180)
zsc=0.















def boom_location(spacing, Cr, alpharad, list_length, ha):
    alphadeg = math.degrees(alpharad)
    nodepos = [[] for _ in range(list_length)]
    #appending the n nodes to location list
    nodepos[0] = [0, 0., -Cr]
    nodepos[1] = [0, nodepos[0][1] + 0.5 * spacing * math.sin(alpharad), nodepos[0][2] + 0.5 * spacing * math.cos(alpharad)]
    nodepos[2] = [0, nodepos[1][1] + spacing * math.sin(alpharad), nodepos[1][2] + spacing * math.cos(alpharad)]
    nodepos[3] = [0, nodepos[2][1] + spacing * math.sin(alpharad), nodepos[2][2] + spacing * math.cos(alpharad)]    
    nodepos[4] = [0, nodepos[3][1] + spacing * math.sin(alpharad), nodepos[3][2] + spacing * math.cos(alpharad)]
    
    
    s_y = ha/2. - nodepos[4][1]
    s_z = -nodepos[4][2]
    dist = math.sqrt(s_y**2 + s_z**2)
    arc = spacing - dist
    theta = 2*arc/ha

    #boog
    nodepos[5] = [0, ha/2*math.cos(theta), ha/2*math.sin(theta)]
    nodepos[6] = [0, 0, ha/2]
    nodepos[7] = [0, -nodepos[5][1], nodepos[5][2]]
    
    check = spacing + arc
    check_two = 1/4*2*math.pi*ha/2
    #print (check)
    #print (check_two)
    
    #negative of upper side
    nodepos[8] = [0, -nodepos[4][1], nodepos[4][2]]
    nodepos[9] = [0, -nodepos[3][1], nodepos[3][2]]
    nodepos[10] = [0, -nodepos[2][1], nodepos[2][2]]
    nodepos[11] = [0, -nodepos[1][1], nodepos[1][2]]
    
    #added booms for spar
    nodepos[12] = [0, ha/2, 0]
    nodepos[13] = [0, -ha/2, 0]

    
    #print (nodepos)
    
    return nodepos, arc, dist

nodepos=boom_location(spacing, Cr, alpharad, list_length, ha)[0]

def boom_area_updater(tsk, spacing, Mz, My, Izz, Iyy, stiff_area, zcg, nodepos, dist, arc, tspar, ha, theta):
    
    My = My*math.cos(theta)+Mz*math.sin(theta)
    Mz = Mz*math.cos(theta)-My*math.sin(theta)
    
    a= -1
    b= (Mz*Iyy)/(My*Izz)
    c= zcg
    
    d=[0]
    for i in range(13):
        x=0.
        x=abs(a*nodepos[i+1][2]+b*nodepos[i+1][1]+c)/(math.sqrt(a**2+b**2))
        if nodepos[i+1][1]+(c+a*nodepos[i+1][2])/b < 0:
            x=x*(-1)
        else:
            x=x
        d.append(x)
    
    boom_area=[0]
    boom_area.append(stiff_area+(tsk*spacing/6)*(4+(d[11]+d[2])/d[1]))
    for i in range(2):
        boom_area.append(stiff_area+(tsk*spacing/6)*(4+(d[i+1]+d[i+3])/d[i+2]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(2+d[3]/d[4])+(tsk*dist/6)*(2+d[12]/d[4]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(2+d[6]/d[5])+(tsk*arc/6)*(2+d[12]/d[5]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(4+(d[5]+d[7])/d[6]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(2+d[6]/d[7])+(tsk*arc/6)*(2+d[13]/d[7]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(2+d[9]/d[8])+(tsk*dist/6)*(2+d[13]/d[8]))
    for i in range(2):
        boom_area.append(stiff_area+(tsk*spacing/6)*(4+(d[i+8]+d[i+10])/d[i+9]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(4+(d[10]+d[1])/d[11]))
    boom_area.append((tsk*dist/6)*(2+d[4]/d[12])+(tsk*arc/6)*(2+d[5]/d[12])+(tspar*ha/6)*(2+d[13]/d[12]))
    boom_area.append((tsk*dist/6)*(2+d[8]/d[13])+(tsk*arc/6)*(2+d[7]/d[13])+(tspar*ha/6)*(2+d[12]/d[13]))
    return boom_area


def shear_flow_finder(boom_area_inclskin, Izz, Iyy, theta, node_pos, Ca, ha, Mx, Vy, Vz, G, tskin, tspar):
    #decompose forces to be local
    Vy = Vy*math.cos(theta)+Vz*math.sin(theta)
    Vz = Vz*math.cos(theta)-Vy*math.sin(theta)
    #counter-clockwise movement
    baes = boom_area_inclskin
    bxyz = node_pos
    #define the order of movement for the triangular and the circular section
    tring_booms = [12,13,8,9,10,11,1,2,3,4]
    circ_booms = [13,12,5,6,7]
    #circle section is I and triangular section is II
    #find areas of sections if not found
    s=math.sqrt((ha/2)**2+(Ca-ha/2)**2)
    peri=(ha/2)*math.pi
    AI=0.5*math.pi*(ha/2)**2
    AII=0.5*(Ca-ha/2)*(ha/2)
    #define distances between boom and next one and associated thicknesses
    spac = 0.1015545
    edge = 0.0766246
    tring_dist = [ha, edge, spac, spac, spac, spac, spac, spac, spac, edge]
    circ_dist = [ha,spac-edge, spac, spac, spac-edge]
    tring_thicc = [tspar, tsk, tsk, tsk, tsk, tsk, tsk, tsk, tsk, tsk]
    circ_thicc = [tspar, tsk, tsk, tsk, tsk]
    #find base shear flow for each cell
    tring_q = [0]
    circ_q = [0]
    for i in tring_booms:
        tring_q.append((-Vy/Izz)*baes[i]*bxyz[i][1]+(-Vz/Iyy)*baes[i]*bxyz[i][2]+tring_q[-1])
    for i in circ_booms:
        circ_q.append((-Vy/Izz)*baes[i]*bxyz[i][1]+(-Vz/Iyy)*baes[i]*bxyz[i][2]+circ_q[-1])
    
    #find force produced by each boom due to shear flows
    tring_fz=[]
    tring_fy=[]
    circ_fz= []
    circ_fy= []
    for i in range (len(tring_booms)):
        if i == len(tring_booms)-1:
            tring_fz.append(tring_q[i+1]*(bxyz[tring_booms[0]][2]-bxyz[tring_booms[i]][2]))
        else:
            tring_fz.append(tring_q[i+1]*(bxyz[tring_booms[i+1]][2]-bxyz[tring_booms[i]][2]))

    for i in range (len(tring_booms)):
        if i == len(tring_booms)-1:
            tring_fy.append(tring_q[i+1]*(bxyz[tring_booms[0]][1]-bxyz[tring_booms[i]][1]))
        else:
            tring_fy.append(tring_q[i+1]*(bxyz[tring_booms[i+1]][1]-bxyz[tring_booms[i]][1]))

    for i in range (len(circ_booms)):
        if i == len(circ_booms)-1:
            circ_fz.append(circ_q[i+1]*(bxyz[circ_booms[0]][2]-bxyz[circ_booms[i]][2]))
        else:
            circ_fz.append(circ_q[i+1]*(bxyz[circ_booms[i+1]][2]-bxyz[circ_booms[i]][2]))

    for i in range (len(circ_booms)):
        if i == len(circ_booms)-1:
            circ_fy.append(circ_q[i+1]*(bxyz[circ_booms[0]][1]-bxyz[circ_booms[i]][1]))
        else:
            circ_fy.append(circ_q[i+1]*(bxyz[circ_booms[i+1]][1]-bxyz[circ_booms[i]][1]))

    
    #find moment due to force
    #counter-clockwise positive
    moments=0
    for i in range(len(tring_fz)):
        moments += tring_fz[i]*bxyz[tring_booms[i]][1]

    for i in range(len(tring_fy)):
        moments += tring_fy[i]*(-1)*bxyz[tring_booms[i]][2]

    for i in range(len(circ_fz)):
        moments += circ_fz[i]*bxyz[tring_booms[i]][1]          

    for i in range(len(circ_fy)):
        moments += circ_fy[i]*(-1)*bxyz[tring_booms[i]][2] 
    print (moments)
    #find line integral of (qbi*ds)/(t*G)
    tring_li=0
    circ_li=0
    for i in range(len(tring_dist)):
        tring_li += (tring_q[i+1]*tring_dist[i])/(tring_thicc[i]*G)

    for i in range(len(circ_dist)):
        circ_li += (circ_q[i+1]*circ_dist[i])/(circ_thicc[i]*G)
    #set up matrix
    A=np.matrix([[1/(2*AI)*(peri/(tsk*G)+ha/(tspar*G)), -ha/(2*AI*tspar*G), -1], [-ha/(2*AII*tspar*G), 1/(2*AII)*(2*s/(tsk*G)+ha/(tspar*G)), -1], [2*AI, 2*AII, 0]])
    b=np.matrix([[1/(-2*AI)*circ_li], [1/(-2*AII)*tring_li], [Mx-moments]])
    #solve matrix for redundant shear flows and rate of twist
    x = np.linalg.solve(A,b)
    qs0I = x.item(0)
    qs0II = x.item(1)
    print (qs0I,qs0II, 1234)
    twist_rate = x.item(2)
    #define order of output shear flows
    circoo=[2,3,4,5,1]
    tringoo=[7,8, 9,10,1,2,3,4,5,6]
    #find total shear flows
    circ_qt=[]
    tring_qt=[]    
    for i in circoo:
        circ_qt.append(circ_q[i]+qs0I)
    for i in tringoo:
        tring_qt.append(tring_q[i]+qs0II)
    
    return twist_rate, circ_qt, tring_qt

boom_area = boom_area_updater(tsk, spacing, Mz, My, Izz, Iyy, stiff_area, zcg, nodepos, dist, arc, tspar, ha, theta)
twist_rate, circ_qt, tring_qt =shear_flow_finder(boom_area, Izz, Iyy, theta, nodepos, Ca, ha, Mx, Vy, Vz, G, tskin, tspar)


sigma_v = Definitions.von_mises_stress (nodepos2, Ilocal, ndis, Mx, My, Mz,circ_qt,tring_qt, tskin, tspar)
