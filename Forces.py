# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:35:45 2019

@author: Sybren
"""
import numpy as np
from math import *
import Definitions

def ReactionForces(theta,P,q,Ca,ha,E,Izz,x1,x2,x3,xa,span,d1,d3): 
    """ Reaction Forces in x direction """
    R2x = 0.                                    #sum of forces in x
    
    """ Reaction Forces in y direction """
    
    #Matrix of equations for solving reaction forces in y
    eq1 = [1.,1.,1.,0.,0.]                      #Sum of forces in y
    eq2 = [-(x2 - x1),0.,(x3-x2),0.,0.]         #sum of moments around hinge 2
    eq3 = [0.,0.,0.,x1,1.]                      #deflection of hinge 1
    eq4 = [(-(x2-x1)**3.)/6., 0.,0.,x2,1.]      #deflection of hinge 2
    eq5 = [(-(x3-x1)**3.)/6., (-(x3-x2)**3.)/6., 0., x3, 1.]    #deflection of hinge 3
    ans1 = [span*q]
    ans2 = [(span/2. - x2)*span*q]
    ans3 = [E*Izz*d1 - q/24.* x1**4.]
    ans4 = [-q/24.*(x2**4.)]
    ans5 = [E*Izz*d3 - q/24.* x3**4.]
    
    A = np.array([eq1, eq2, eq3, eq4, eq5])
    b = np.array([ans1,ans2,ans3,ans4,ans5])
    # x = [R1y, R2y, R3y, A, B], where A and B are the integration constants from the bending equation
    x = np.linalg.solve(A, b)
    
    #check solution
    check1 = np.allclose(np.dot(A, x), b)
    
    #Get reaction forces from x
    R1y = float(x[0])
    R2y = float(x[1])
    R3y = float(x[2])
    
    """ Reaction Forces in z direction """
    A1 = -span*q/(tan(theta)) - P                #force in actuator 1 (sum of moments around hinge)
    
    eq1 = [1.,1.,1.,0.,0.]                      #sum of forces in z
    eq2 = [-(x2 - x1),0.,(x3-x2),0.,0.]         #sum of moments round hinge 2
    eq3 = [0.,0.,0.,x1,1.]                      #zero deflection of hinge 1
    eq4 = [((x2-x1)**3.)/6., 0.,0.,x2,1.]       #zero deflection of hinge 2
    eq5 = [((x3-x1)**3.)/6., ((x3-x2)**3.)/6., 0., x3, 1.]    #zero deflection of hinge 3
    ans1 = [-(A1+P)]
    ans2 = [(A1*(xa/2.) - P*(xa/2.))]
    ans3 = [0.]
    ans4 = [-A1/6.*((xa/2.)**3.)]
    ans5 = [-A1/6.*(x3-x2+(xa/2.))**3. -P/6.*(x3-x2-(xa/2.))**3.]
    
    A = np.array([eq1, eq2, eq3, eq4, eq5])
    b = np.array([ans1,ans2,ans3,ans4,ans5])
    # y = [R1z, R2z, R3z, A, B], where A and B are the integration constants from the bending equation
    y = np.linalg.solve(A, b)
    
    #check solution
    check2 = np.allclose(np.dot(A, y), b)
    
    #Get reaction forces from y
    R1z = float(y[0])
    R2z = float(y[1])
    R3z = float(y[2])
    
    return R1y, R1z, R2x, R2y, R2z, R3y, R3z, A1
    
def ExactMOI(theta,Ca,ha,tskin,Ct,tspar,nodepos,area_stiff)

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
    
def Torsion(theta,P,q,Ca,ha,A1):
    move = 0.25*Ca - ha/2.
    movez = move*cos(theta)
    movey = move*sin(theta)
    
    M_add_P = P*movey
    M_add_q = q*movex

    


theta = (30./180.*pi)
P = 49200.
q = 3860.
Ca = 0.505
ha = 0.161
x1 = 0.125
x2 = 0.498
x3 = 1.494
xa = 0.245
span = 1.611
d1 = 0.00389
d3 = 0.01245
E = 73084430000         #N/m2
Izz = 0.0013            #m4
test = ReactionForces(theta,P,q,Ca,ha,E,Izz,x1,x2,x3,xa,span,d1,d3)