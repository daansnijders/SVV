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
    
  
#ha = 0.161 m, Ca = 0.505 and n = 11
def boom_spacing(ha, Ca, n):
    Cr = Ca-ha/2.
    Ct = math.sqrt(Cr**2+(ha/2)**2)
    Cb = Ct
    Cc = math.pi*ha/2.
    Circ = Ct+Cb+Cc
    spacing = Circ/n
    alpharad = math.atan2((ha/2), Cr)
    
    return spacing, Cr, alpharad

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
    boom_area[6] = area_stiff + tskin*spacing/6*(2+nodepos[5][1]/nodepos[6][1]) + tskin*spacing/6*(2+nodepos[7][1]/nodepos[6][1])
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


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    