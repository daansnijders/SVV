# -*- coding: utf-8 -*-
import math

def beambending(ray,rby,rcy,rdy,q,dx,moi,centroid):

    print('test')


def centroid():
    
    print ('Test')
    
    
def Boom_spacing(ha, Ca, n):
    Cr = Ca-ha/2.
    Ct = math.sqrt(Cr**2+(ha/2)**2)
    Cb = Ct
    Cc = math.pi*ha
    Circ = Ct+Cb+Cc
    Spacing = Circ/n
    print ('The spacing is', Spacing)
    
    return Spacing