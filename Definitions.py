# -*- coding: utf-8 -*-
import math
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np

"""EXAMPLE AREAS, TO BE DEFINED / PROGRAMED / FILLED IN"""
B = [[0.0], [0.1], [0.2], [0.3], [0.4], [0.5], [0.6], [0.5], [0.4], [0.3], [0.2], [0.1], [0.12], [0.13]] #This is now a 2D list but could be transfered to 1D?

def beambending(ray,rby,rcy,rdy,q,dx,moi,centroid):
    print('test')


def centroid(nodepos, boom_area, list_length):  #TO BE CHECKED
    frac1 = 0.
    frac2 = 0.
    for i in range (list_length):    
       frac1 = frac1 + nodepos[i][2]*boom_area[i]
       frac2 = frac2 + boom_area[i]
    
    centroid = frac1/frac2
    
    return centroid

def centroid_nonidealized(tskin, ha, Ca, Ct, tspar, nodepos, area_stiff):
    ycg = 0
    
    #moment LE
    A_LE = tskin*0.5*math.pi*ha
    Qyy_LE = A_LE * (2 * (ha/2) / math.pi)
    
    #moment TE
    A_TE = 2 * tskin * Ct
    Qyy_TE = A_TE * -(Ca - ha/2) / 2
    
    #moment spar
    A_sp = tspar * ha
    Qyy_sp = A_sp * 0
    
    Qyy_stiff = 0
    
    for i in range (11):
        Qyy_stiff =+ area_stiff * nodepos[i+1][2]
        
    A_stiff = 11 * area_stiff
    Qyy_sum = Qyy_LE + Qyy_TE + Qyy_sp + Qyy_stiff
    A_sum = A_LE + A_TE + A_sp + A_stiff
    
    zcg = Qyy_sum / A_sum
    
    return ycg, zcg
  
#ha = 0.161 m, Ca = 0.505 and n = 11
def boom_spacing(ha, Ca, n):
    Cr = Ca-ha/2.
    Ct = math.sqrt(Cr**2+(ha/2)**2)
    Cb = Ct
    Cc = math.pi*ha/2.
    Circ = Ct+Cb+Cc
    spacing = Circ/n
    alpharad = math.atan2((ha/2), Cr)
    
    return spacing, Cr, alpharad, Ct

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

def area_stiff(t_stiff, h_stiff, w_stiff):
    area_stiff = t_stiff * h_stiff + w_stiff * t_stiff
    return area_stiff

def boom_area_inclskin(tskin, tspar, spacing, nodepos, area_stiff, dist, arc, ha):
    boom_area = 14*[0]
    boom_area[1] = area_stiff + tskin*spacing/6*(2+nodepos[11][1]/nodepos[1][1]) + tskin*spacing/6*(2+nodepos[2][1]/nodepos[1][1])
    boom_area[2] = area_stiff + tskin*spacing/6*(2+nodepos[1][1]/nodepos[2][1]) + tskin*spacing/6*(2+nodepos[3][1]/nodepos[2][1])
    boom_area[3] = area_stiff + tskin*spacing/6*(2+nodepos[2][1]/nodepos[3][1]) + tskin*spacing/6*(2+nodepos[4][1]/nodepos[3][1])
    boom_area[4] = area_stiff + tskin*spacing/6*(2+nodepos[3][1]/nodepos[4][1]) + tskin*dist/6*(2+nodepos[12][1]/nodepos[4][1])
    boom_area[5] = area_stiff + tskin*arc/6*(2+nodepos[12][1]/nodepos[5][1]) + tskin*spacing/6*(2+nodepos[6][1]/nodepos[5][1])
    boom_area[6] = area_stiff 
    boom_area[7] = area_stiff + tskin*spacing/6*(2+nodepos[6][1]/nodepos[7][1]) + tskin*arc/6*(2+nodepos[13][1]/nodepos[7][1])
    boom_area[8] = area_stiff + tskin*dist/6*(2+nodepos[13][1]/nodepos[8][1]) + tskin*spacing/6*(2+nodepos[9][1]/nodepos[8][1])
    boom_area[9] = area_stiff + tskin*spacing/6*(2+nodepos[8][1]/nodepos[9][1]) + tskin*spacing/6*(2+nodepos[10][1]/nodepos[9][1])
    boom_area[10] = area_stiff + tskin*spacing/6*(2+nodepos[9][1]/nodepos[10][1]) + tskin*spacing/6*(2+nodepos[11][1]/nodepos[10][1])
    boom_area[11] = area_stiff + tskin*spacing/6*(2+nodepos[10][1]/nodepos[11][1]) + tskin*spacing/6*(2+nodepos[1][1]/nodepos[11][1])
    boom_area[12] = tskin*dist/6*(2+nodepos[4][1]/nodepos[12][1]) + tskin*arc/6*(2+nodepos[5][1]/nodepos[12][1]) + tspar*ha/6*(2+nodepos[13][1]/nodepos[12][1])
    boom_area[13] = tskin*spacing/6*(2+nodepos[11][1]/nodepos[13][1]) + tskin*spacing/6*(2+nodepos[2][1]/nodepos[13][1]) + tspar*ha/6*(2+nodepos[12][1]/nodepos[13][1])

    return boom_area
    
def boom_area_exclskin(area_stiff, nodepos, tspar, ha):
    boom_area = 14*[0]
    boom_area[1] = area_stiff 
    boom_area[2] = area_stiff 
    boom_area[3] = area_stiff 
    boom_area[4] = area_stiff 
    boom_area[5] = area_stiff 
    boom_area[6] = area_stiff 
    boom_area[7] = area_stiff 
    boom_area[8] = area_stiff 
    boom_area[9] = area_stiff 
    boom_area[10] = area_stiff
    boom_area[11] = area_stiff 
    boom_area[12] = tspar*ha/6*(2+nodepos[13][1]/nodepos[12][1])
    boom_area[13] = tspar*ha/6*(2+nodepos[12][1]/nodepos[13][1])

    return boom_area
    

def boom_inertia(list_length, nodepos): #TO BE CHECKED
    
    #nodepos = boom_location() #getting positions from previous function

    Ixx = [0 for _ in range(list_length)] #Defining lists of Ixx Iyy and Izz
    Iyy = [] #Etc.
    Izz = [] #Etc.
    
    for i in range(list_length):
        Iyy.append( B[i][0] * (nodepos[i][2]) ** 2) # Iyy = Boom Area * Z distance squared
        Izz.append( B[i][0] * (nodepos[i][1]) ** 2) # Izz = Boom Area * Y distance squared
        
    #print()
    #print("Ixx is: ",Ixx)
    #print()
    #print("Iyy is: ",Iyy)
    #print("Izz is: ",Izz)
    #print()    
    #print()
    
    Ixx_final = 0.
    Iyy_final = 0.
    Izz_final = 0.
    for i in range (14):
        Ixx_final = Ixx_final + Ixx[i]
        Iyy_final = Iyy_final + Iyy[i]
        Izz_final = Izz_final + Izz[i]
    
    return Ixx_final, Iyy_final, Izz_final

def scatter(nodepos):
    y = []
    z = []
    for i in range (14):
        y.append(nodepos[i][1])
        z.append(nodepos[i][2])
    
    plt.ylim((-0.2,0.2))
    plt.xlim((-0.5,0.1))
    plt.scatter(z,y)
    plt.show()
    
def find_shear_center(boom_area_excluding_skin,Izz,node_pos,ha):
    #finding shear center
    #counter-clockwise movement
    baes = boom_area_excluding_skin
    bxyz = node_pos
    #define the order of booms for circular and triangular cells
    tring_booms = [1,2,3,4,12,13,8,9,10,11]
    circ_booms = [12,5,6,7,13]
    #define distances between boom and the next one 
    spac = 0.1015545
    edge = 0.0766246
    tring_dist = [spac, spac, spac, edge, ha, edge, spac, spac, spac, spac]
    circ_dist = [spac-edge, spac, spac, spac-edge, ha]
    #find base shear flow for each cell
    tring_q = [0]
    circ_q = [0]
    for i in tring_booms:
        tring_q.append(-(1/Izz)*baes[i]*bxyz[i][1]+tring_q[-1])
    for j in circ_booms:
        circ_q.append(-(1/Izz)*baes[j]*bxyz[j][1]+circ_q[-1])
    #find redundant shear flow
    tring_qr = 0
    circ_qr = 0
    for i in range(len(tring_dist)):
        tring_qr += tring_q[i+1]*tring_dist[i]
    tring_qr=tring_qr/(sum(tring_dist))

    for j in range(len(circ_dist)):
        circ_qr += circ_q[j+1]*circ_dist[i]
    circ_qr=circ_qr/(sum(circ_dist))
    #find total shear flow
    tring_qt = [x+tring_qr for x in tring_q]
    circ_qt = [x+circ_qr for x in circ_q]
    #find force produced by each boom
    tring_fz=[]
    tring_fy=[]
    circ_fz= []
    circ_fy= []
    for i in range (len(tring_booms)):
        if i == len(tring_booms)-1:
            tring_fz.append(tring_qt[i+1]*(bxyz[tring_booms[0]][2]-bxyz[tring_booms[i]][2]))
        else:
            tring_fz.append(tring_qt[i+1]*(bxyz[tring_booms[i+1]][2]-bxyz[tring_booms[i]][2]))

    for i in range (len(tring_booms)):
        if i == len(tring_booms)-1:
            tring_fy.append(tring_qt[i+1]*(bxyz[tring_booms[0]][1]-bxyz[tring_booms[i]][1]))
        else:
            tring_fy.append(tring_qt[i+1]*(bxyz[tring_booms[i+1]][1]-bxyz[tring_booms[i]][1]))

    for i in range (len(circ_booms)):
        if i == len(circ_booms)-1:
            circ_fz.append(circ_qt[i+1]*(bxyz[circ_booms[0]][2]-bxyz[circ_booms[i]][2]))
        else:
            circ_fz.append(circ_qt[i+1]*(bxyz[circ_booms[i+1]][2]-bxyz[circ_booms[i]][2]))

    for i in range (len(circ_booms)):
        if i == len(circ_booms)-1:
            circ_fy.append(circ_qt[i+1]*(bxyz[circ_booms[0]][1]-bxyz[circ_booms[i]][1]))
        else:
            circ_fy.append(circ_qt[i+1]*(bxyz[circ_booms[i+1]][1]-bxyz[circ_booms[i]][1]))

    #find sum of moments about center of coordinate system
    #clockwise positive
    moments=0
    for i in range(len(tring_fz)):
        moments += tring_fz[i]*(-1)*bxyz[tring_booms[i]][1]

    for i in range(len(tring_fy)):
        moments += tring_fy[i]*bxyz[tring_booms[i]][2]

    for i in range(len(circ_fz)):
        moments += circ_fz[i]*(-1)*bxyz[tring_booms[i]][1]          

    for i in range(len(circ_fy)):
        moments += circ_fy[i]*bxyz[tring_booms[i]][2] 



    #value will equal z distance of shear center
    sc_position= [0,0,moments]
    return sc_position
                           

def shear_flow_shear(boom_area_incl_skin, node_pos, Vy, Vz):
    #finding shear center
    #counter-clockwise movement
    baes = boom_area_incl_skin
    bxyz = node_pos
    #define the order of booms for circular and triangular cells
    tring_booms = [1,2,3,4,12,13,8,9,10,11]
    circ_booms = [12,5,6,7,13]
    #define distances between boom and the next one 
    spac = 0.1015545
    edge = 0.0766246
    tring_dist = [spac, spac, spac, edge, ha, edge, spac, spac, spac, spac]
    circ_dist = [spac-edge, spac, spac, spac-edge, ha]
    #find base shear flow for each cell
    tring_q = [0]
    circ_q = [0]
    for i in tring_booms:
        tring_q.append(-(Vy/Izz)*baes[i-1]*bxyz[i-1][1]+tring_q[-1]-(Vz/Iyy)*baes[i-1]*bxyz[i-1][1]+tring_q[-1])
    for j in circ_booms:
        circ_q.append(-(Vy/Izz)*baes[j-1]*bxyz[j-1][1]+circ_q[-1]-(Vz/Iyy)*baes[j-1]*bxyz[j-1][1]+circ_q[-1]) 
    
     #find redundant shear flow
    tring_qr = 0
    circ_qr = 0
    for i in range(len(tring_dist)):
        tring_qr += tring_q[i+1]*tring_dist[i]
    tring_qr=tring_qr/(sum(tring_dist))

    for j in range(len(circ_dist)):
        circ_qr += circ_q[j+1]*circ_dist[i]
    circ_qr=circ_qr/(sum(circ_dist))
    #find total shear flow
    tring_qt = [x+tring_qr for x in tring_q]
    circ_qt = [x+circ_qr for x in circ_q]

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
