# -*- coding: utf-8 -*-
import math

ha = 0.161 #height of thickest part of airfoil
Ca = 0.505 #dist from trailing to leading edge
n = 11 #number of stringers
list_length = n+3 #number of idealised booms INCLUDING two added AND Ghost nodeloc


"""EXAMPLE AREAS, TO BE DEFINED / PROGRAMED / FILLED IN"""
B = [[0], [0.1], [0.2], [0.3], [0.4], [0.5], [0.6], [0.5], [0.4], [0.3], [0.2], [0.1], [0.12], [0.13]]
"""END OF EXAMPLE AREAS"""

def beambending(ray,rby,rcy,rdy,q,dx,moi,centroid):

    print('test')


def centroid():
    
    print ('Test')
    
  
#ha = 0.161 m, Ca = 0.505 and n = 11
def boom_spacing(ha, Ca, n):
    Cr = Ca-ha/2.
    Ct = math.sqrt(Cr**2+(ha/2)**2)
    Cb = Ct
    Cc = math.pi*ha
    Circ = Ct+Cb+Cc
    spacing = Circ/n
    alpharad = math.atan2((ha/2), Cr)
    
    return spacing, Cr, alpharad

def boom_location():
    spacing, Cr, alpharad = boom_spacing(ha, Ca, n)
    alphadeg = math.degrees(alpharad)
    nodeloc = [[] for _ in range(list_length)]
    #appending the n nodes to location list
    nodeloc[0] = [0, 0, -Cr]
    nodeloc[1] = [0, 0.5 * spacing * math.sin(alphadeg), 0.5 * spacing * math.cos(alphadeg)]
    nodeloc[2] = [0, nodeloc[1][1] + spacing * math.sin(alphadeg), nodeloc[1][2] + spacing * math.cos(alphadeg)]
    nodeloc[3] = [0, nodeloc[2][1] + spacing * math.sin(alphadeg), nodeloc[2][2] + spacing * math.cos(alphadeg)]
    nodeloc[3] = [0, nodeloc[3][1] + spacing * math.sin(alphadeg), nodeloc[3][2] + spacing * math.cos(alphadeg)]    
    
    print (nodeloc)
    
    
    
    
    
    
    
    
def boom_inertia():

    Ixx = [[] for _ in range(list_length)] #Moments of inertia over the x axis
    Iyy = [[] for _ in range(list_length)] #Moments of inertia over the y axis
    Izz = [[] for _ in range(list_length)] #Etc.
    
    
    """ As its idealized as booms, only the Steiner terms are aken into account:
        eg:     Ixx = A * dy ^ 2            """ 
    
    Ixx = Ixx.append(0)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    