# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:35:45 2019

@author: Sybren
"""
import numpy as np
from math import *

def ReactionForces(theta,P,q,Ca,ha,E,Ixx,Iyy,x1,x2,x3,xa,span,d1,d3): 
    """ Reaction Forces in x direction """
    Rx2 = 0.                                    #sum of forces in x
    
    """ Reaction Forces in y direction """
    
    #Matrix of equations for solving reaction forces in y
    eq1 = [1.,1.,1.,0.,0.]                      #Sum of forces in y
    eq2 = [-(x2 - x1),0.,(x3-x2),0.,0.]         #sum of moments around hinge 2
    eq3 = [0.,0.,0.,x1,1.]                      #deflection of hinge 1
    eq4 = [(-(x2-x1)**3.)/6., 0.,0.,x2,1.]      #deflection of hinge 2
    eq5 = [(-(x3-x1)**3.)/6., (-(x3-x2)**3.)/6., 0., x3, 1.]    #deflection of hinge 3
    ans1 = [span*q]
    ans2 = [(span/2. - x2)*span*q]
    ans3 = [E*Iyy*d1 - q/24.* x1**4.]
    ans4 = [-q/24.*(x2**4.)]
    ans5 = [E*Ixx*d3 - q/24.* x3**4.]
    
    A = np.array(eq1, eq2, eq3, eq4, eq5)
    b = np.array(ans1,ans2,ans3,ans4,ans5)
    # x = [R1y, R2y, R3y, A, B], where A and B are the integration constants from the bending equation
    x = np.linalg.solve(A, b)
    
    #check solution
    check = np.allclose(np.dot(A, x), b)
    
    #Get reaction forces from x
    Ry1 = x[0]
    Ry2 = x[1]
    Ry3 = x[2]
    
    """ Reaction Forces in z direction """
    A1 = -span*q/tan(theta) - P                #force in actuator 1 (sum of moments around hinge)
    
    eq1 = [1.,1.,1.,0.,0.]                      #sum of forces in z
    eq2 = [-(x2 - x1),0.,(x3-x2),0.,0.]         #sum of moments round hinge 2
    eq3 = [0.,0.,0.,x1,1.]                      #zero deflection of hinge 1
    eq4 = [((x2-x1)**3.)/6., 0.,0.,x2,1.]      #zero deflection of hinge 2
    eq5 = [((x3-x1)**3.)/6., ((x3-x2)**3.)/6., 0., x3, 1.]    #zero deflection of hinge 3
    ans1 = [-(A1+P)]
    ans2 = [(A1*(xa/2.) - P*(xa/2.))]
    ans3 = [0.]
    ans4 = [-A1/6.*((xa/2.)**3.)]
    ans5 = [-A1/6.*(x3-x2+(xa/2.))**3. -P/6.*(x3-x2-(xa/2.))**3.]
    
    A = np.array(eq1, eq2, eq3, eq4, eq5)
    b = np.array(ans1,ans2,ans3,ans4,ans5)
    # y = [R1z, R2z, R3z, A, B], where A and B are the integration constants from the bending equation
    y = np.linalg.solve(A, b)
    
    #check solution
    check = np.allclose(np.dot(A, y), b)
    
    #Get reaction forces from y
    Rz1 = y[0]
    Rz2 = y[1]
    Rz3 = y[2]
    
    return R2x, R1y, R2y, R3y, R1z, R2z, R3z
    
    
    
def Torsion(theta,P,q,Ca,ha):
    move = 0.25*Ca - ha/2.
    movez = move*cos(theta)
    movey = move*sin(theta)
    
    M_add_P = P*movey
    M_add_q = q*movex
    

