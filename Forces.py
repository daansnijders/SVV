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
    
    return R2x, R1y, R2y, R3y, R1z, R2z, R3z
    
        
def ExactMOI(theta,Ca,ha,t_sk,t_sp,t_st,w_st,h_st,zcg,no_st,spacing)
    
    """ Stringer MOI """
    y_bar = (t_st*w*t_st/2. + (h_st-t_st)*t*((h_st-t_st)/2. + t_st)) / (w_st*t_st + (h_st-t_st)*t_st)
    A = (w_st*t_st + (h_st-t_st)*t_st)
    Izz_st = 1./12. *w_st*t_st**3. + (w_st*t_st)*(y_bar - t_st/2.) + 1./12.*(h_st - t_st)**3.*t_st + (h_st-t_st)*t_st * ((h_st-t_st)/2. - y_bar)
    Iyy_st = 1./12 *(h_st-t_st)*t_st**3. + 1./12. * w_st**3. *t
    
    """ Half arc MOI """
    Izz_arc = 1./2. * pi * ha**3. * t_sk
    Iyy_arc = 1./2. * pi * ha**3. * t_sk
    
    """ Spar MOI """
    Izz_sp = 1./12. * ha**3. * t_sp
    Iyy_sp = 1./12. * t_sp**3. * ha
    
    """ Beam MOI """
    a = sqrt((Ca - ha/2.)**2. + (ha/2.)**2.)
    angle = arctan(ha/(2.*(Ca - ha/2.)))
    Izz_beam = t_sk * a**3. *(sin(angle))**2. / 12.
    Iyy_beam = t_sk * a**3. *(cos(angle))**2. / 12.
    
    """ Overall MOI """
    Izz = 0.
    Iyy = 0.
    loc = []
    #For loop for all stringers
    for i in range(len(no_st)):
        if i <= (no_st - 1.)/2. - 2):
            z_loc = Ca - ha/2. - zcg - (i+0.5)*spacing/cos(theta) 
            y_loc = (i+0.5)*spacing/sin(theta)
            Izz_new = (Iyy_st + Izz_st)/2. + (Izz_st - Iyy_st)/2. * cos(2.*angle)
            Iyy_new = (Iyy_st + Izz_st)/2. - (Izz_st - Iyy_st)/2. * cos(2.*angle)
        elif i == (no_st - 1.)/2. - 1):
            z_loc = Ca - ha/2. - zcg - (i+0.5)*spacing/cos(theta) 
            y_loc = (i-0.5)*spacing/sin(theta)
            Izz_new = (Iyy_st + Izz_st)/2. + (Izz_st - Iyy_st)/2. * cos(-2.*angle)
            Iyy_new = (Iyy_st + Izz_st)/2. - (Izz_st - Iyy_st)/2. * cos(-2.*angle)
        elif i == (no_st - 1.)/2. :
            z_loc = - zg - ha/2.
            y_loc = 0.
            Izz_new = Ixx_st
            Iyy_new = Izz_st
        elif i == (no_st - 1.)/2. + 1):
            z_loc = Ca - ha/2. - zcg - (no_st - 0.5 - i)*spacing/cos(theta) 
            y_loc = (no_st - 1.5 - i)*spacing/sin(theta)
            Izz_new = (Iyy_st + Izz_st)/2. + (Izz_st - Iyy_st)/2. * cos(2.*angle)
            Iyy_new = (Iyy_st + Izz_st)/2. - (Izz_st - Iyy_st)/2. * cos(2.*angle)
        else:
            z_loc = Ca - ha/2. - zcg - (no_st - 0.5 - i)*spacing/cos(theta) 
            y_loc = (no_st - 0.5 - i)*spacing/sin(theta)
            Izz_new = (Iyy_st + Izz_st)/2. + (Izz_st - Iyy_st)/2. * cos(-2.*angle)
            Iyy_new = (Iyy_st + Izz_st)/2. - (Izz_st - Iyy_st)/2. * cos(-2.*angle)
        
        Izz += Izz_new + A * (z_loc)**2.
        Iyy += Iyy_new + A * (y_loc)**2.
        loc.append([z_loc, y_loc])
    
    #Add half arc
    Izz += Izz_arc
    Iyy += Iyy_arc + (pi*((ha/2.)**2 - (ha/2. - t_sk)**2.)/2.)* (2.*ha/2./pi + zcg)**2.
    
    #Add spar
    Izz += Izz_sp
    Iyy += Iyy_sp + (t_sp * ha)*zcg**2.
    
    #Add beams
    Izz += (Izz_beam + (Ca-ha/2.)*cos(angle)*t_sk * (((Ca-ha/2.)/2.)*tan(angle))**2.)*2.
    Iyy += (Iyy_beam + (Ca-ha/2.)*cos(angle)*t_sk *(((Ca-ha/2.)/2.)-zcg)**2.)
    
    Izz_0 = Izz
    Iyy_0 = Iyy        
    
    Iyy_theta = (Izz_0 + Iyy_0)/2. + (Izz_0 - Iyy_0)/2. * cos(2.*theta)
    Izz_theta = (Izz_0 + Iyy_0)/2. - (Izz_0 - Iyy_0)/2. * cos(2.*theta)
    Izy_theta = (Izz_0 - Iyy_0)/2. *sin(2.*theta)
    
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