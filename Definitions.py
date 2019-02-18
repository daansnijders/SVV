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
   
    
    return spacing, Cr

def boom_placement():
    spacing, Cr = boom_spacing(ha, Ca, n)
    node = list()
    node = node.append([0, 0, -Cr])
    print (node[0])