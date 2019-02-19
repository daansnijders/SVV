# -*- coding: utf-8 -*-
import math

ha = 0.161 #height of thickest part of airfoil
Ca = 0.505 #dist from trailing to leading edge
n = 11 #number of stringers
list_length = n+3 #number of idealised booms INCLUDING two added AND Ghost nodepos

"""EXAMPLE AREAS, TO BE DEFINED / PROGRAMED / FILLED IN"""
B = [[0.0], [0.1], [0.2], [0.3], [0.4], [0.5], [0.6], [0.5], [0.4], [0.3], [0.2], [0.1], [0.12], [0.13]] #This is now a 2D list but could be transfered to 1D?

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
    nodepos = [[] for _ in range(list_length)]
    #appending the n nodes to location list
    nodepos[0] = [0, 0., -Cr]
    nodepos[1] = [0, 0.5 * spacing * math.sin(alphadeg), 0.5 * spacing * math.cos(alphadeg)]
    nodepos[2] = [0, nodepos[1][1] + spacing * math.sin(alphadeg), nodepos[1][2] + spacing * math.cos(alphadeg)]
    nodepos[3] = [0, nodepos[2][1] + spacing * math.sin(alphadeg), nodepos[2][2] + spacing * math.cos(alphadeg)]
    nodepos[3] = [0, nodepos[3][1] + spacing * math.sin(alphadeg), nodepos[3][2] + spacing * math.cos(alphadeg)]    
    """THE FOLLOWING ONES ARE NOT YET CORRECT BUT REQUIRE FURTHER CODING"""
    nodepos[4] = [0, nodepos[1][1] + spacing * math.sin(alphadeg), nodepos[1][2] + spacing * math.cos(alphadeg)]
    nodepos[5] = [0, nodepos[2][1] + spacing * math.sin(alphadeg), nodepos[2][2] + spacing * math.cos(alphadeg)]
    nodepos[6] = [0, nodepos[3][1] + spacing * math.sin(alphadeg), nodepos[3][2] + spacing * math.cos(alphadeg)]
    nodepos[7] = [0, nodepos[1][1] + spacing * math.sin(alphadeg), nodepos[1][2] + spacing * math.cos(alphadeg)]
    nodepos[8] = [0, nodepos[2][1] + spacing * math.sin(alphadeg), nodepos[2][2] + spacing * math.cos(alphadeg)]
    nodepos[9] = [0, nodepos[3][1] + spacing * math.sin(alphadeg), nodepos[3][2] + spacing * math.cos(alphadeg)]
    nodepos[10] = [0, nodepos[1][1] + spacing * math.sin(alphadeg), nodepos[1][2] + spacing * math.cos(alphadeg)]
    nodepos[11] = [0, nodepos[2][1] + spacing * math.sin(alphadeg), nodepos[2][2] + spacing * math.cos(alphadeg)]
    nodepos[12] = [0, nodepos[3][1] + spacing * math.sin(alphadeg), nodepos[3][2] + spacing * math.cos(alphadeg)]
    nodepos[13] = [0, nodepos[1][1] + spacing * math.sin(alphadeg), nodepos[1][2] + spacing * math.cos(alphadeg)]

    
    print (nodepos)
    
    return nodepos

    
def boom_inertia(): #TO BE CHECKED
    
    nodepos = boom_location() #getting positions from previous function

    Ixx = [0 for _ in range(list_length)] #Defining lists of Ixx Iyy and Izz
    Iyy = [] #Etc.
    Izz = [] #Etc.
    
    for i in range(list_length):
        Iyy.append( B[i][0] * (nodepos[i][2]) ** 2) # Iyy = Boom Area * Z distance squared
        Izz.append( B[i][0] * (nodepos[i][1]) ** 2) # Izz = Boom Area * Y distance squared
        
    
    print()
    print("Ixx is: ",Ixx)
    print()
    print("Iyy is: ",Iyy)
    print("Izz is: ",Izz)
    print()    
    print()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    