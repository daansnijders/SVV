# -*- coding: utf-8 -*-
import math

ha = 0.161
Ca = 0.505
n = 11

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

def boom_locations():
    spacing, Cr, alpharad = boom_spacing(ha, Ca, n)
    alphadeg = math.degrees(alpharad)
    list_length = n+3
    node = [[] for _ in range(list_length)]
    #appending the n nodes to location list
    node[0] = [0, 0, -Cr]
    node[1] = [0, 0.5 * spacing * math.sin(alphadeg), 0.5 * spacing * math.cos(alphadeg)]
    node[2] = [0, node[1][1] + spacing * math.sin(alphadeg), node[1][2] + spacing * math.cos(alphadeg)]
    node[3] = [0, node[2][1] + spacing * math.sin(alphadeg), node[2][2] + spacing * math.cos(alphadeg)]
    node[3] = [0, node[3][1] + spacing * math.sin(alphadeg), node[3][2] + spacing * math.cos(alphadeg)]    
    
    
    
    print (node[0:4])