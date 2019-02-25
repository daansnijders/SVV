# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:35:45 2019

@author: Sybren
"""
import numpy as np
import math
import Definitions
import matplotlib.pyplot as plt

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
    A1 = -span*q*(0.25*Ca - ha/2.)/(ha/2. * (math.cos(theta) - math.sin(theta))) - P  #force in actuator 1 (sum of moments around hinge)
    
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
    
        
def ExactMOI(theta,Ca,ha,t_sk,t_sp,t_st,w_st,h_st,zcg,n,spacing,nodepos):
    
    """ Stringer MOI """
    y_bar = (t_st*w_st*t_st/2. + (h_st-t_st)*t_st*((h_st-t_st)/2. + t_st)) / (w_st*t_st + (h_st-t_st)*t_st)
    A = (w_st*t_st + (h_st-t_st)*t_st)
    Izz_st = 1./12. *w_st*t_st**3. + (w_st*t_st)*(y_bar - t_st/2.) + 1./12.*(h_st - t_st)**3.*t_st + (h_st-t_st)*t_st * ((h_st-t_st)/2. - y_bar)
    Iyy_st = 1./12 *(h_st-t_st)*t_st**3. + 1./12. * w_st**3. *t_st
    
    """ Half arc MOI """
    Izz_arc = 1./2. * math.pi * ha**3. * t_sk
    Iyy_arc = 1./2. * math.pi * ha**3. * t_sk
    
    """ Spar MOI """
    Izz_sp = 1./12. * ha**3. * t_sp
    Iyy_sp = 1./12. * t_sp**3. * ha
    
    """ Beam MOI """
    a = math.sqrt((Ca - ha/2.)**2. + (ha/2.)**2.)
    angle = math.atan(ha/(2.*(Ca - ha/2.)))
    Izz_beam = t_sk * a**3. *(math.sin(angle))**2. / 12.
    Iyy_beam = t_sk * a**3. *(math.cos(angle))**2. / 12.
    
    """ Overall MOI """
    Izz = 0.
    Iyy = 0.
    locy = []
    locz = []
    #For loop for all stringers
    for i in range(n):
        if i <= ((n - 1.)/2. - 2):
            z_loc = Ca - ha/2. - abs(zcg) - (i+0.5)*spacing*math.cos(angle) 
            y_loc = (i+0.5)*spacing*math.sin(angle)
            Izz_new = (Iyy_st + Izz_st)/2. + (Izz_st - Iyy_st)/2. * math.cos(2.*angle)
            Iyy_new = (Iyy_st + Izz_st)/2. - (Izz_st - Iyy_st)/2. * math.cos(2.*angle)
        elif i == ((n - 1.)/2. - 1):
            z_loc = -abs(nodepos[i+1][2]) - abs(zcg) 
            y_loc = abs(nodepos[i+1][1])
            Izz_new = (Iyy_st + Izz_st)/2. + (Izz_st - Iyy_st)/2. * math.cos(-2.*angle)
            Iyy_new = (Iyy_st + Izz_st)/2. - (Izz_st - Iyy_st)/2. * math.cos(-2.*angle)
        elif i == (n - 1.)/2. :
            z_loc = - abs(zcg) - ha/2.
            y_loc = 0.
            Izz_new = Iyy_st
            Iyy_new = Izz_st
        elif i == ((n - 1.)/2. + 1):
            z_loc = -abs(nodepos[i+1][2]) - abs(zcg)
            y_loc = nodepos[i+1][1]
            Izz_new = (Iyy_st + Izz_st)/2. + (Izz_st - Iyy_st)/2. * math.cos(2.*angle)
            Iyy_new = (Iyy_st + Izz_st)/2. - (Izz_st - Iyy_st)/2. * math.cos(2.*angle)
        else:
            z_loc = Ca - ha/2. - abs(zcg) - (n - 0.5 - i)*spacing*math.cos(angle) 
            y_loc = -(n - 0.5 - i)*spacing*math.sin(angle)
            Izz_new = (Iyy_st + Izz_st)/2. + (Izz_st - Iyy_st)/2. * math.cos(-2.*angle)
            Iyy_new = (Iyy_st + Izz_st)/2. - (Izz_st - Iyy_st)/2. * math.cos(-2.*angle)
        
        Izz += Izz_new + A * (z_loc)**2.
        Iyy += Iyy_new + A * (y_loc)**2.
        locy.append(y_loc)
        locz.append(z_loc)
        
    plt.plot(locz, locy, "bo")
    plt.grid(True)
    plt.show
    
    #Add half arc
    Izz += Izz_arc
    Iyy += Iyy_arc + (math.pi*((ha/2.)**2 - (ha/2. - t_sk)**2.)/2.)* (2.*ha/2./math.pi + zcg)**2.
    
    #Add spar
    Izz += Izz_sp
    Iyy += Iyy_sp + (t_sp * ha)*zcg**2.
    
    #Add beams
    Izz += (Izz_beam + (Ca-ha/2.)*math.cos(angle)*t_sk * (((Ca-ha/2.)/2.)*math.tan(angle))**2.)*2.
    Iyy += (Iyy_beam + (Ca-ha/2.)*math.cos(angle)*t_sk *(((Ca-ha/2.)/2.)-zcg)**2.)
    
    Izz_0 = Izz
    Iyy_0 = Iyy        
    
    Iyy_theta = (Izz_0 + Iyy_0)/2. + (Izz_0 - Iyy_0)/2. * math.cos(2.*theta)
    Izz_theta = (Izz_0 + Iyy_0)/2. - (Izz_0 - Iyy_0)/2. * math.cos(2.*theta)
    Izy_theta = (Izz_0 - Iyy_0)/2. *math.sin(2.*theta)
    
    return Iyy_0, Izz_0, Iyy_theta, Izz_theta, Izy_theta
    
def InternalMoment(x,q,A1,P,x1,x2,x3,xa,R1y,R2y,R3y,R1z,R2z,R3z):
    # X starts from the left most point on the aileron, before hinge 1
    Mz = q/2. * x**2.
    My = 0.
    
    if x > x1 :
        Mz += -R1y * (x - x1)
        My +=  R1z * (x - x1)
        
    if x > (x2 - xa/2.) : 
        My +=  A1 * (x - (x2 - xa/2.))
    
    if x > x2:
        Mz += -R2y * (x - x2)
        My +=  R2z * (x - x2)
        
    if x > (x2 + xa/2.) : 
        My +=  P * (x - (x2 + xa/2.))
        
    if x > x3:
        Mz += -R3y * (x - x3)
        My +=  R3z * (x - x3)
    
    return Mz, My

theta = (30./180.*math.pi)
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
t_sk = 1.1 * 10**(-3)
t_sp = 2.4 * 10**(-3)
t_st = 1.2 * 10**(-3)
w_st = 0.017
h_st = 0.013
A = (w_st*t_st + (h_st-t_st)*t_st)
n = 11
list_length = 14
spacing, Cr, alpharad, Ct = Definitions.boom_spacing(ha, Ca, n)
nodepos, arc, dist = Definitions.boom_location(spacing, Cr, alpharad, list_length, ha)
ycg, zcg = Definitions.centroid_nonidealized(t_sk, ha, Ca, Ct, t_sp, nodepos, A)
Iyy_0, Izz_0, Iyy_theta, Izz_theta, Izy_theta = ExactMOI(theta,Ca,ha,t_sk,t_sp,t_st,w_st,h_st,zcg,n,spacing,nodepos)
test = ReactionForces(theta,P,q,Ca,ha,E,Izz_theta,x1,x2,x3,xa,span,d1,d3)


