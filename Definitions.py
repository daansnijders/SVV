# -*- coding: utf-8 -*-
import math
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np


def deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1):
    #theta2 is the angle of section 2 where hinge 2 is
    
    # calculation of reaction forces by static equalibrium
    
    r2 = (-q*(l4)*(l4/2-l1)-r3*(l3-l1))/(l2-l1)
    r1 = -q*l4-r2-r3

    #Theta2 is actuator 2 and theta is actuator 1

    P1 = (P2*(-ha/2*np.sin(theta2)+ha/2*np.cos(theta2))-q*l4*(0.25*Ca-ha/2)*np.cos(theta2))/(-ha/2*np.sin(theta1)+ha/2*np.cos(theta1))
    rz2 = (-P1*(l2-l1-xa/2)+P2*(l2-l1+xa/2)-rz3*(l3-l1))/(l2-l1)
    rz1 = P2-P1-rz3-rz2

    #Moment, curvature, slope and deflection array assignment

    Vy = np.array([0])
    Vz = np.array([0])

    
    Mz = np.array([0])
    My = np.array([0])
    
    vdouble = np.array([0])
    vsingle = np.array([0])
    v = np.array([0])

    udouble = np.array([0])
    usingle = np.array([0])
    u = np.array([0])

    xt = np.array([0])


    i = 1

    for x in np.linspace(0,l1,n+1)[1:]:
        xt = np.append(xt, x)
        dx = l1/(n)

        Vy = np.append(Vy, q*x)
        Vz = np.append(Vz, 0)
        
        Mz = np.append(Mz, -q*(x**2)/2)
        My = np.append(My, 0)

        
        vdouble = np.append(vdouble,-(Mz[i-1]*I[0][0][i-1]-My[i-1]*I[0][1][i-1])/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))
        vsingle = np.append(vsingle,(-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**3)/6)*I[0][0][i-1]-(0)*I[0][1][i-1]))
        v = np.append(v, (-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**4)/24)*I[0][0][i-1]-(0)*I[0][1][i-1]))

        udouble = np.append(udouble, -(-Mz[i-1]*I[0][1][i]+My[i-1]*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        usingle = np.append(usingle,-(-(-q*(x**3)/6)*I[0][1][i]+0*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        u = np.append(u, -(-(-q*(x**4)/24)*I[0][1][i]+0*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))

        i = i+1

        
    for x in np.linspace(l1,l2-xa/2,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l2-xa/2-l1)/(n)

        Vy = np.append(Vy, q*x+r1)
        Vz = np.append(Vz, rz1)
        
        Mz = np.append(Mz, -q*(x**2)/2-r1*(x-l1))
        My = np.append(My, rz1*(x-l1))


        vdouble = np.append(vdouble,-(Mz[i-1]*I[0][0][i-1]-My[i-1]*I[0][1][i-1])/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))
        vsingle = np.append(vsingle,(-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**3)/6-r1*(x-l1)**2/2)*I[0][0][i-1]-
                                                    (rz1*(x-l1)**2/2)*I[0][1][i-1]))
        v = np.append(v, (-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**4)/24-r1*(x-l1)**3/6)*I[0][0][i-1]-
                                                    (rz1*(x-l1)**3/6)*I[0][1][i-1]))

        udouble = np.append(udouble, -(-Mz[i-1]*I[0][1][i]+My[i-1]*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        usingle = np.append(usingle,-(-(-q*(x**3)/6-r1*(x-l1)**2/2)*I[0][1][i]+(rz1*(x-l1)**2/2)*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        u = np.append(u, -(-(-q*(x**4)/24-r1*(x-l1)**3/6)*I[0][1][i]+(rz1*(x-l1)**3/6)*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))

        i = i+1

    for x in np.linspace(l2-xa/2,l2,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l2-l2+xa/2)/(n)

        Vy = np.append(Vy, q*x+r1)
        Vz = np.append(Vz, rz1+P1)
        
        Mz = np.append(Mz, -q*(x**2)/2-r1*(x-l1))
        My = np.append(My, rz1*(x-l1)+P1*(x-l2+xa/2))


        vdouble = np.append(vdouble,-(Mz[i-1]*I[0][0][i-1]-My[i-1]*I[0][1][i-1])/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))
        vsingle = np.append(vsingle,(-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**3)/6-r1*(x-l1)**2/2)*I[0][0][i-1]-
                                                    (rz1*(x-l1)**2/2+P1*(x-l2+xa/2)**2/2)*I[0][1][i-1]))
        v = np.append(v, (-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**4)/24-r1*(x-l1)**3/6)*I[0][0][i-1]-
                                                    (rz1*(x-l1)**3/6+P1*(x-l2+xa/2)**3/6)*I[0][1][i-1]))

        udouble = np.append(udouble, -(-Mz[i-1]*I[0][1][i]+My[i-1]*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        usingle = np.append(usingle,-(-(-q*(x**3)/6-r1*(x-l1)**2/2)*I[0][1][i]+(rz1*(x-l1)**2/2+P1*(x-l2+xa/2)**2/2)*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        u = np.append(u, -(-(-q*(x**4)/24-r1*(x-l1)**3/6)*I[0][1][i]+(rz1*(x-l1)**3/6+P1*(x-l2+xa/2)**3/6)*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))

        i = i+1

        
    for x in np.linspace(l2,l2+xa/2,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l2+xa/2-l2)/(n)

        Vy = np.append(Vy, q*x+r1+r2)
        Vz = np.append(Vz, rz1+P1+rz2)
        
        Mz = np.append(Mz, -q*(x**2)/2-r1*(x-l1)-r2*(x-l2))
        My = np.append(My, rz1*(x-l1)+P1*(x-l2+xa/2)+rz2*(x-l2))

        vdouble = np.append(vdouble,-(Mz[i-1]*I[0][0][i-1]-My[i-1]*I[0][1][i-1])/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))
        vsingle = np.append(vsingle,(-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**3)/6-r1*(x-l1)**2/2-r2*(x-l2)**2/2)*I[0][0][i-1]-
                                                    (rz1*(x-l1)**2/2+P1*(x-l2+xa/2)**2/2+
                                                     rz2*(x-l2)**2/2)*I[0][1][i-1]))
        v = np.append(v, (-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**4)/24-r1*(x-l1)**3/6-r2*(x-l2)**3/6)*I[0][0][i-1]-
                                                    (rz1*(x-l1)**3/6+P1*(x-l2+xa/2)**3/6+
                                                     rz2*(x-l2)**3/6)*I[0][1][i-1]))



        udouble = np.append(udouble, -(-Mz[i-1]*I[0][1][i]+My[i-1]*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        usingle = np.append(usingle,-(-(-q*(x**3)/6-r1*(x-l1)**2/2-r2*(x-l2)**2/2)*I[0][1][i]+(rz1*(x-l1)**2/2+P1*(x-l2+xa/2)**2/2+
                                                     rz2*(x-l2)**2/2)*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        u = np.append(u, -(-(-q*(x**4)/24-r1*(x-l1)**3/6-r2*(x-l2)**3/6)*I[0][1][i]+(rz1*(x-l1)**3/6+P1*(x-l2+xa/2)**3/6+
                                                     rz2*(x-l2)**3/6)*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        
        i = i+1

    for x in np.linspace(l2+xa/2,l3,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l3-l2-xa/2)/(n)

        Vy = np.append(Vy, q*x+r1+r2)
        Vz = np.append(Vz, rz1+P1+rz2-P2)
        
        Mz = np.append(Mz, -q*(x**2)/2-r1*(x-l1)-r2*(x-l2))
        My = np.append(My, rz1*(x-l1)+P1*(x-l2+xa/2)+rz2*(x-l2)-P2*(x-l2-xa/2))

        vdouble = np.append(vdouble,-(Mz[i-1]*I[0][0][i-1]-My[i-1]*I[0][1][i-1])/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))
        vsingle = np.append(vsingle,(-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**3)/6-r1*(x-l1)**2/2-r2*(x-l2)**2/2)*I[0][0][i-1]-
                                                    (rz1*(x-l1)**2/2+P1*(x-l2+xa/2)**2/2+
                                                     rz2*(x-l2)**2/2-P2*(x-l2-xa/2)**2/2)*I[0][1][i-1]))
        v = np.append(v, (-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**4)/24-r1*(x-l1)**3/6-r2*(x-l2)**3/6)*I[0][0][i-1]-
                                                    (rz1*(x-l1)**3/6+P1*(x-l2+xa/2)**3/6+
                                                     rz2*(x-l2)**3/6-P2*(x-l2-xa/2)**3/6)*I[0][1][i-1]))



        udouble = np.append(udouble, -(-Mz[i-1]*I[0][1][i]+My[i-1]*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        usingle = np.append(usingle,-(-(-q*(x**3)/6-r1*(x-l1)**2/2-r2*(x-l2)**2/2)*I[0][1][i]+(rz1*(x-l1)**2/2+P1*(x-l2+xa/2)**2/2+
                                                     rz2*(x-l2)**2/2-P2*(x-l2-xa/2)**2/2)*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        u = np.append(u, -(-(-q*(x**4)/24-r1*(x-l1)**3/6-r2*(x-l2)**3/6)*I[0][1][i]+(rz1*(x-l1)**3/6+P1*(x-l2+xa/2)**3/6+
                                                     rz2*(x-l2)**3/6-P2*(x-l2-xa/2)**3/6)*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        
        i = i+1
    for x in np.linspace(l3,l4,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l4-l3)/(n)

        Vy = np.append(Vy, q*x+r1+r2+r3)
        Vz = np.append(Vz, rz1+P1+rz2-P2+rz3)
        
        Mz = np.append(Mz, -q*(x**2)/2-r1*(x-l1)-r2*(x-l2)-r3*(x-l3))
        My = np.append(My, rz1*(x-l1)+P1*(np.sign((x-l2+xa/2)) == 1)*(x-l2+xa/2)+rz2*(x-l2)-P2*(np.sign((x-l2-xa/2)) == 1)*(x-l2-xa/2)+rz3*(x-l3))
                                                
        
        vdouble = np.append(vdouble,-(Mz[i-1]*I[0][0][i-1]-My[i-1]*I[0][1][i-1])/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))
        vsingle = np.append(vsingle,(-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**3)/6-r1*(x-l1)**2/2-r2*(x-l2)**2/2-r3*(x-l3)**2/2)*I[0][0][i-1]-
                                                    (rz1*(x-l1)**2/2+P1*(np.sign((x-l2+xa/2)) == 1)*(x-l2+xa/2)**2/2+
                                                     rz2*(x-l2)**2/2-P2*(np.sign((x-l2-xa/2)) == 1)*(x-l2-xa/2)**2/2+rz3*(x-l3)**2/2)*I[0][1][i-1]))
        v = np.append(v, (-1/(E*(I[0][0][i-1]*I[1][1][i-1]-I[0][1][i-1]**2)))*((-q*(x**4)/24-r1*(x-l1)**3/6-r2*(x-l2)**3/6-r3*(x-l3)**3/6)*I[0][0][i-1]-
                                                    (rz1*(x-l1)**3/6+P1*(np.sign((x-l2+xa/2)) == 1)*(x-l2+xa/2)**3/6+
                                                     rz2*(x-l2)**3/6-P2*(np.sign((x-l2-xa/2)) == 1)*(x-l2-xa/2)**3/6+rz3*(x-l3)**3/6)*I[0][1][i-1]))

        udouble = np.append(udouble, -(-Mz[i-1]*I[0][1][i]+My[i-1]*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        usingle = np.append(usingle,-(-(-q*(x**3)/6-r1*(x-l1)**2/2-r2*(x-l2)**2/2-r3*(x-l3)**2/2)*I[0][1][i]+(rz1*(x-l1)**2/2+P1*(np.sign((x-l2+xa/2)) == 1)*(x-l2+xa/2)**2/2+
                                                     rz2*(x-l2)**2/2-P2*(np.sign((x-l2-xa/2)) == 1)*(x-l2-xa/2)**2/2+rz3*(x-l3)**2/2)*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        u = np.append(u, -(-(-q*(x**4)/24-r1*(x-l1)**3/6-r2*(x-l2)**3/6-r3*(x-l3)**3/6)*I[0][1][i]+(rz1*(x-l1)**3/6+P1*(np.sign((x-l2+xa/2)) == 1)*(x-l2+xa/2)**3/6+
                                                     rz2*(x-l2)**3/6-P2*(np.sign((x-l2-xa/2)) == 1)*(x-l2-xa/2)**3/6+rz3*(x-l3)**3/6)*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        
        i = i+1

    
    deltav = ((v[3*n]-d2)*l1-(v[n]-d1)*l2)/(l1-l2)
    deltavsingle = ((v[n]-d1)-deltav)/l1

    v2 = -deltavsingle*xt-deltav+v

    deltau = ((u[3*n]-0)*l1-(u[n]-0)*l2)/(l1-l2)
    deltausingle = ((u[n]-0)-deltau)/l1
    u2 = -deltausingle*xt-deltau+u
            
    return [v2, u2, xt, r1, r2, Vy, Vz, My, Mz, rz1, rz2, P1]



def bendingconvergence(q,n,l1,l2,l3,l4,E,I,d1,d2,d3,P2,xa,Ca,ha,theta2,theta1,spread):

    #spread refers to the dx term in the jacobian matrix

    r3 = 0
    rz3 = 0
    
    def dfdx(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1,spread):
        
        dvdr3 = (deflection(q,n,r3+spread/2,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)[0][5*n]-deflection(q,n,r3-spread/2,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)[0][5*n])/spread
        dvdrz3 = (deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3+spread/2,P2,xa,Ca,ha,theta2,theta1)[0][5*n]-deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3-spread/2,P2,xa,Ca,ha,theta2,theta1)[0][5*n])/spread    
        dudr3 = (deflection(q,n,r3+spread/2,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)[1][5*n]-deflection(q,n,r3-spread/2,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)[1][5*n])/spread
        dudrz3 = (deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3+spread/2,P2,xa,Ca,ha,theta2,theta1)[1][5*n]-deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3-spread/2,P2,xa,Ca,ha,theta2,theta1)[1][5*n])/spread
        
        return np.matrix([[dvdr3, dvdrz3],[dudr3,dudrz3]])

    v,u = deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)[0][5*n], deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)[1][5*n]

    while round(v, 12) != d3 or round(u, 12) != 0:
    
        invjacobian = np.linalg.inv(dfdx(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1,spread))

        dx = invjacobian * np.matrix([[deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)[0][5*n]-d3],[deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)[1][5*n]]])


        r3,rz3 = r3-dx[0], rz3-dx[1]

        v,u = deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)[0][5*n], deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)[1][5*n]

        print('V-Error = ',abs(v-d3),' and U-Error = ', u)

    
    
    v2, u2, xt, r1, r2, Vy, Vz, My, Mz, rz1, rz2, P1 = deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2,theta1)
    
    return v2, u2, xt, r1, r2, r3, Vy, Vz, My, Mz, rz1, rz2, rz3, P1


def torque(q,n,l1,l2,l3,l4,P1,P2,xa,Ca,ha,theta,zsc):
    
    Mx = np.array([0])

    xt = np.array([0])


    i = 1

    for x in np.linspace(0,l1,n+1)[1:]:
        xt = np.append(xt, x)
        dx = l1/(n)

        Mx = np.append(Mx, (q*x)*(-zsc+0.25*Ca-ha/2)*np.cos(theta[i-1]))


        i = i+1
    for x in np.linspace(l1,l2-xa/2,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l2-l1)/(n)
        
       
        Mx = np.append(Mx, (q*x)*(-zsc+0.25*Ca-ha/2)*np.cos(theta[i-1]))

        i = i+1

    for x in np.linspace(l2-xa/2,l2,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l2-l1)/(n)
        
       
        Mx = np.append(Mx, (q*x)*(-zsc+0.25*Ca-ha/2)*np.cos(theta[i-1])+P1*(-ha/2*np.sin(theta[2*n])+ha/2*np.cos(theta[2*n])))

        i = i+1
    for x in np.linspace(l2,l2+xa/2,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l3-l2)/(n)
        
        Mx = np.append(Mx, (q*x)*(-zsc+0.25*Ca-ha/2)*np.cos(theta[i-1])+P1*(-zsc*np.sin(theta[2*n])+(-ha/2*np.sin(theta[2*n])+ha/2*np.cos(theta[2*n]))))
        
        i = i+1
    for x in np.linspace(l2+xa/2,l3,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l3-l2)/(n)
        
        Mx = np.append(Mx, (q*x)*(-zsc+0.25*Ca-ha/2)*np.cos(theta[i-1])+P1*(-zsc*np.sin(theta[2*n])+(-ha/2*np.sin(theta[2*n])+ha/2*np.cos(theta[2*n])))
                       - P2*(-zsc*np.sin(theta[4*n])+(-ha/2*np.sin(theta[4*n])+ha/2*np.cos(theta[4*n]))))
        
        i = i+1
    for x in np.linspace(l3,l4,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l4-l3)/(n)

        Mx = np.append(Mx, (q*x)*(-zsc+0.25*Ca-ha/2)*np.cos(theta[i-1])+P1*(-zsc*np.sin(theta[2*n])+(-ha/2*np.sin(theta[2*n])+ha/2*np.cos(theta[2*n])))
                       - P2*(-zsc*np.sin(theta[4*n])+(-ha/2*np.sin(theta[4*n])+ha/2*np.cos(theta[4*n]))))
       

        i = i+1

    return Mx, xt


def centroid(nodepos, boom_area, list_length):  #TO BE CHECKED
    ycg = 0
    frac1 = 0.
    frac2 = 0.
    for i in range (list_length):    
       frac1 = frac1 + nodepos[i][2]*boom_area[i]
       frac2 = frac2 + boom_area[i]
    
    zcg = frac1/frac2
    
    return ycg, zcg

def centroid_nonidealized(tskin, ha, Ca, Ct, tspar, nodepos, area_stiff):
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
  
#ha = 0.161 m, Ca = 0.505 and n = 11
def boom_spacing(ha, Ca, n):
    Cr = Ca-ha/2.
    Ct = math.sqrt(Cr**2+(ha/2)**2)
    Cb = Ct
    Cc = math.pi*ha/2.
    Circ = Ct+Cb+Cc
    spacing = Circ/n
    alpharad = math.atan2((ha/2), Cr)
    
    return spacing, Cr, alpharad, Ct

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
#
#def boom_area_inclskin(tskin, tspar, spacing, nodepos, area_stiff, dist, arc, ha):
#    boom_area = 14*[0]
#    boom_area[1] = area_stiff + tskin*spacing/6*(2+nodepos[11][1]/nodepos[1][1]) + tskin*spacing/6*(2+nodepos[2][1]/nodepos[1][1])
#    boom_area[2] = area_stiff + tskin*spacing/6*(2+nodepos[1][1]/nodepos[2][1]) + tskin*spacing/6*(2+nodepos[3][1]/nodepos[2][1])
#    boom_area[3] = area_stiff + tskin*spacing/6*(2+nodepos[2][1]/nodepos[3][1]) + tskin*spacing/6*(2+nodepos[4][1]/nodepos[3][1])
#    boom_area[4] = area_stiff + tskin*spacing/6*(2+nodepos[3][1]/nodepos[4][1]) + tskin*dist/6*(2+nodepos[12][1]/nodepos[4][1])
#    boom_area[5] = area_stiff + tskin*arc/6*(2+nodepos[12][1]/nodepos[5][1]) + tskin*spacing/6*(2+nodepos[6][1]/nodepos[5][1])
#    boom_area[6] = area_stiff 
#    boom_area[7] = area_stiff + tskin*spacing/6*(2+nodepos[6][1]/nodepos[7][1]) + tskin*arc/6*(2+nodepos[13][1]/nodepos[7][1])
#    boom_area[8] = area_stiff + tskin*dist/6*(2+nodepos[13][1]/nodepos[8][1]) + tskin*spacing/6*(2+nodepos[9][1]/nodepos[8][1])
#    boom_area[9] = area_stiff + tskin*spacing/6*(2+nodepos[8][1]/nodepos[9][1]) + tskin*spacing/6*(2+nodepos[10][1]/nodepos[9][1])
#    boom_area[10] = area_stiff + tskin*spacing/6*(2+nodepos[9][1]/nodepos[10][1]) + tskin*spacing/6*(2+nodepos[11][1]/nodepos[10][1])
#    boom_area[11] = area_stiff + tskin*spacing/6*(2+nodepos[10][1]/nodepos[11][1]) + tskin*spacing/6*(2+nodepos[1][1]/nodepos[11][1])
#    boom_area[12] = tskin*dist/6*(2+nodepos[4][1]/nodepos[12][1]) + tskin*arc/6*(2+nodepos[5][1]/nodepos[12][1]) + tspar*ha/6*(2+nodepos[13][1]/nodepos[12][1])
#    boom_area[13] = tskin*spacing/6*(2+nodepos[11][1]/nodepos[13][1]) + tskin*spacing/6*(2+nodepos[2][1]/nodepos[13][1]) + tspar*ha/6*(2+nodepos[12][1]/nodepos[13][1])
#
#    return boom_area
#    
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
    

def boom_inertia(list_length, nodepos, B): #TO BE CHECKED
    
    #nodepos = boom_location() #getting positions from previous function

    Ixx = [0 for _ in range(list_length)] #Defining lists of Ixx Iyy and Izz
    Iyy = [] #Etc.
    Izz = [] #Etc.
    
    for i in range(list_length):
        Iyy.append( B[i] * (nodepos[i][2]) ** 2) # Iyy = Boom Area * Z distance squared
        Izz.append( B[i] * (nodepos[i][1]) ** 2) # Izz = Boom Area * Y distance squared
        
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
    
def find_shear_center(boom_area_excluding_skin,Izz,node_pos,ha):
    #finding shear center
    #counter-clockwise movement
    baes = boom_area_excluding_skin
    bxyz = node_pos
    #define the order of booms for circular and triangular cells
    tring_booms = [1,2,3,4,12,13,8,9,10,11]
    circ_booms = [12,5,6,7,13]
    #define distances between boom and the next one 
    spac = 0.1015545
    edge = 0.0766246
    tring_dist = [spac, spac, spac, edge, ha, edge, spac, spac, spac, spac]
    circ_dist = [spac-edge, spac, spac, spac-edge, ha]
    #find base shear flow for each cell
    tring_q = [0]
    circ_q = [0]
    for i in tring_booms:
        tring_q.append(-(1/Izz)*baes[i]*bxyz[i][1]+tring_q[-1])
    for j in circ_booms:
        circ_q.append(-(1/Izz)*baes[j]*bxyz[j][1]+circ_q[-1])
    #find redundant shear flow
    tring_qr = 0
    circ_qr = 0
    for i in range(len(tring_dist)):
        tring_qr += tring_q[i+1]*tring_dist[i]
    tring_qr=-tring_qr/(sum(tring_dist))

    for j in range(len(circ_dist)):
        circ_qr += circ_q[j+1]*circ_dist[i]
    circ_qr=-circ_qr/(sum(circ_dist))
    #find total shear flow
    tring_qt = [x+tring_qr for x in tring_q]
    circ_qt = [x+circ_qr for x in circ_q]
    #find force produced by each boom
    tring_fz=[]
    tring_fy=[]
    circ_fz= []
    circ_fy= []
    for i in range (len(tring_booms)):
        if i == len(tring_booms)-1:
            tring_fz.append(tring_qt[i+1]*(bxyz[tring_booms[0]][2]-bxyz[tring_booms[i]][2]))
        else:
            tring_fz.append(tring_qt[i+1]*(bxyz[tring_booms[i+1]][2]-bxyz[tring_booms[i]][2]))

    for i in range (len(tring_booms)):
        if i == len(tring_booms)-1:
            tring_fy.append(tring_qt[i+1]*(bxyz[tring_booms[0]][1]-bxyz[tring_booms[i]][1]))
        else:
            tring_fy.append(tring_qt[i+1]*(bxyz[tring_booms[i+1]][1]-bxyz[tring_booms[i]][1]))

    for i in range (len(circ_booms)):
        if i == len(circ_booms)-1:
            circ_fz.append(circ_qt[i+1]*(bxyz[circ_booms[0]][2]-bxyz[circ_booms[i]][2]))
        else:
            circ_fz.append(circ_qt[i+1]*(bxyz[circ_booms[i+1]][2]-bxyz[circ_booms[i]][2]))

    for i in range (len(circ_booms)):
        if i == len(circ_booms)-1:
            circ_fy.append(circ_qt[i+1]*(bxyz[circ_booms[0]][1]-bxyz[circ_booms[i]][1]))
        else:
            circ_fy.append(circ_qt[i+1]*(bxyz[circ_booms[i+1]][1]-bxyz[circ_booms[i]][1]))

    #find sum of moments about center of coordinate system
    #clockwise positive
    moments=0
    for i in range(len(tring_fz)):
        moments += tring_fz[i]*(-1)*bxyz[tring_booms[i]][1]

    for i in range(len(tring_fy)):
        moments += tring_fy[i]*bxyz[tring_booms[i]][2]

    for i in range(len(circ_fz)):
        moments += circ_fz[i]*(-1)*bxyz[tring_booms[i]][1]          

    for i in range(len(circ_fy)):
        moments += circ_fy[i]*bxyz[tring_booms[i]][2] 



    #value will equal z distance of shear center
    sc_position= [0,0,moments]
    return sc_position, tring_booms, circ_booms
                           

def shear_flow_shear(boom_area_incl_skin, node_pos, Vy, Vz,ha,Izz,Iyy):
    #finding shear center
    #counter-clockwise movement
    baes = boom_area_incl_skin
    bxyz = node_pos
    #define the order of booms for circular and triangular cells
    tring_booms = [1,2,3,4,12,13,8,9,10,11]
    circ_booms = [12,5,6,7,13]
    #define distances between boom and the next one 
    spac = 0.1015545
    edge = 0.0766246
    tring_dist = [spac, spac, spac, edge, ha, edge, spac, spac, spac, spac]
    circ_dist = [spac-edge, spac, spac, spac-edge, ha]
    #find base shear flow for each cell
    tring_q = [0]
    circ_q = [0]
    for i in tring_booms:
        tring_q.append(-(Vy/Izz)*baes[i-1]*bxyz[i-1][1]-(Vz/Iyy)*baes[i-1]*bxyz[i-1][1]+tring_q[-1])
    for j in circ_booms:
        circ_q.append(-(Vy/Izz)*baes[j-1]*bxyz[j-1][1]-(Vz/Iyy)*baes[j-1]*bxyz[j-1][1]+circ_q[-1]) 
    
     #find redundant shear flow
    tring_qr = 0
    circ_qr = 0
    for i in range(len(tring_dist)):
        tring_qr += tring_q[i+1]*tring_dist[i]

    tring_qr=-tring_qr/(sum(tring_dist))

    for j in range(len(circ_dist)):
        circ_qr += circ_q[j+1]*circ_dist[i]
    circ_qr=-circ_qr/(sum(circ_dist))

    tring_qr= -tring_qr/(sum(tring_dist))

    for j in range(len(circ_dist)):
        circ_qr += circ_q[j+1]*circ_dist[i]
    circ_qr= -circ_qr/(sum(circ_dist))

    #find total shear flow
    tring_qt = [x+tring_qr for x in tring_q]
    circ_qt = [x+circ_qr for x in circ_q]

    
    return tring_qt, circ_qt

def shear_flow_torsion(T,A1,A2,arc,l,ha,G,t):
    
    # T= resultant torque applied to the cross section
    # A= area cell
    # arc= lenght of the leading edge semicircle
    # ha= diameter of the leading edge semi circle
    # l= lenght of the triangular section (from tip of the triangle to intersection skin-spar)
    # G=shear modulus
    # t= skin thickness 
    

    A=np.matrix([[0,2*A1,2*A2],[-1,(arc+ha)/(2*A1*G*t),-ha/(2*A1*G*t)],[-1,-ha/(2*A2*G*t),(2*l+ha)/(2*A2*G*t)]])
    b=np.matrix([[T],[0],[0]])
    x = np.linalg.solve(A,b)
    rate_twist=x.item(0)
    q1=x.item(1) # shear flow due to torsion in cell 1
    q2=x.item(2) # shear flow due to torsion in cell 2

    
    return rate_twist,q1,q2

def shear_flow_total(tring_qt,circ_qt,q1,q2):
    
    tring_qsum=[] # list containing all the total shearflows along the triangular section 
    # given by the contribution of the shear flow due to shear and the shearflow due to torsion
    
    circ_qsum=[]  # list containing all the total shearflows along the circular section 
    # given by the contribution of the shear flow due to shear and the shearflow due to torsion
    
    for i in tring_qt:
        
        tring_qsum_element=i+q1
        tring_qsum.append(tring_qsum_element)
        
    for j in circ_qt:
        
        circ_qsum_element=j+q2
        circ_qsum.append(circ_qsum_element)
        
        return tring_qsum, circ_qsum

def shear_flow_rib(tring_qsum,circ_qsum,nodepos,ha,circ_booms,tring_booms,alpharad):
    
    lst_tri=[]  # list containing the wing skin shear flows for the triangular cell excluding the shear flows along the spar
    lst_circ=[] # list containing the wing skin shear flows for the circular cell excluding the shear flows along the spar
    
    
    for i in range(len(tring_qsum)):
        
        lst_tri.append(tring_qsum(i))
        
        if i==5:
            
            lst_tri.revove(tring_qsum(i))
    
    for j in range(len(circ_qsum)):
        
        lst_circ.append(circ_qsum(i))
        
        if i==5:
            lst_circ.remove(circ_qsum(i))
    
    # rib shear flow circular cell
    
    Sy1=0.
    r=0.
    for i in lst_circ:
        
        r=r+1
        Sy1=Sy1+i*(nodepos[circ_booms(r)][1]-nodepos[circ_booms(r-1)][1]) # total vertical shear force acting due to the wing skin shear flow in the circular cell
    
    qrib_1 = Sy1/ha
    
    # rib shear flow triangular cell
    
    Sy2=0.
    r=0.
    
    for j in lst_tri:
        
        r=r+1
        Sy2=Sy2+i*(nodepos[tring_booms(r)][1]-nodepos[tring_booms(r-1)][1]) # total vertical shear force acting due to the wing skin shear flow in the triangular cell
        
    lst_tri_mom=[] # list containing the wing skin triangular shear flow which do create a moment about boom 13
    
    i=0.
    while i<=3:
        lst_tri_mom.append(lst_tri(i))
        i+=1
    
    
    mom_sum=0. # total sum of the moments genarated by the wing skin shear flows about boom n 13
    mom_sumy=0. # sum of the moments generated by the wing skin shear forces acting in the y direction about boom 13
    mom_sumz=0. # sum of the moments generated by the wing skin shear forces acting in the z direction about boom 13
    r=0.
    
    for i in lst_tri_mom:
        
        r=r+1
        mom_sumy=mom_sumy+i*(nodepos[r][1]-nodepos[r-1][1])*(nodepos[r][0]-nodepos[r-1][0])
        mom_sumz=mom_sumz+i*(nodepos[r][0]-nodepos[r-1][0])*(ha/2+(nodepos[r][1]-nodepos[r-1][1]))
        mom_sum=mom_sum+mom_sumz+mom_sumy
        
    Pz=mom_sum/ha # flange force in the z direction
    P=Pz/np.cos(alpharad)
    Py=P*np.sin(alpharad)
    
    #Shear force carried by the web
    
    Sw=Sy2-2*Py
    qrib_2=Sw/ha # shear flow carried by web 2
    
    return qrib_1,qrib_2
        
    
        
    
    
def boom_area_updater(tsk, spacing, Mz, My, Izz, Iyy, stiff_area, zcg, nodepos, dist, arc, tspar, ha):
    a= -1
    b= -(Mz*Iyy)/(My*Izz)
    c= zcg
    
    d=[0]
    for i in range(13):
        x=0.
        x=abs(a*nodepos[i+1][2]+b*nodepos[i+1][1]+c)/(math.sqrt(a**2+b**2))
        if nodepos[i+1][1]+(c+a*nodepos[i+1][2])/b < 0:
            x=x*(-1)
        else:
            x=x
        d.append(x)
    
    boom_area=[0]
    boom_area.append(stiff_area+(tsk*spacing/6)*(4+(d[11]+d[2])/d[1]))
    for i in range(2):
        boom_area.append(stiff_area+(tsk*spacing/6)*(4+(d[i+1]+d[i+3])/d[i+2]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(2+d[3]/d[4])+(tsk*dist/6)*(2+d[12]/d[4]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(2+d[6]/d[5])+(tsk*arc/6)*(2+d[12]/d[5]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(4+(d[5]+d[7])/d[6]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(2+d[6]/d[7])+(tsk*arc/6)*(2+d[13]/d[7]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(2+d[9]/d[8])+(tsk*dist/6)*(2+d[13]/d[8]))
    for i in range(2):
        boom_area.append(stiff_area+(tsk*spacing/6)*(4+(d[i+8]+d[i+10])/d[i+9]))
    boom_area.append(stiff_area+(tsk*spacing/6)*(4+(d[10]+d[1])/d[11]))
    boom_area.append((tsk*dist/6)*(2+d[4]/d[12])+(tsk*arc/6)*(2+d[5]/d[12])+(tspar*ha/6)*(2+d[13]/d[12]))
    boom_area.append((tsk*dist/6)*(2+d[8]/d[13])+(tsk*arc/6)*(2+d[7]/d[13])+(tspar*ha/6)*(2+d[12]/d[13]))
    return boom_area
            
        
    
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
    y_bar = (h_st**2. * t_st + t_st**2. * (w_st - t_st)) / (2. * (w_st*t_st + (h_st - t_st)*t_st))
    A = (w_st*t_st + (h_st-t_st)*t_st)
    Izz_st = 1./12. *w_st*t_st**3. + (w_st*t_st)*(y_bar - t_st/2.)**2. + 1./12.*(h_st - t_st)**3.*t_st + (h_st-t_st)*t_st * ((h_st-t_st)/2. - y_bar)**2.
    Iyy_st = 1./12 *(h_st-t_st)*t_st**3. + 1./12. * w_st**3. *t_st
    
    """ Half arc MOI """
    Izz_arc = 1./2. * math.pi * (ha/2.)**3. * t_sk
    Iyy_arc = 1./2. * math.pi * (ha/2.)**3. * t_sk
    
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
        
        Izz += Izz_new + A * (y_loc)**2.
        Iyy += Iyy_new + A * (z_loc)**2.
        locy.append(y_loc)
        locz.append(z_loc)
        
#    plt.plot(locz, locy, "bo")
#    plt.grid(True)
#    plt.show
    
    #Add half arc
    Izz += Izz_arc
    Iyy += Iyy_arc + (math.pi*((ha/2.)**2 - (ha/2. - t_sk)**2.)/2.)* (2.*ha/2./math.pi + abs(zcg))**2.
    
    #Add spar
    Izz += Izz_sp
    Iyy += Iyy_sp + (t_sp * ha)*(abs(zcg))**2.
    
    #Add beams
    Izz += (Izz_beam + (Ca-ha/2.)/math.cos(angle)*t_sk * (((Ca-ha/2.)/2.)*math.tan(angle))**2.)*2.
    Iyy += (Iyy_beam + (Ca-ha/2.)/math.cos(angle)*t_sk *(((Ca-ha/2.)/2.)-abs(zcg))**2.)*2.
    
    Izz_0 = Izz
    Iyy_0 = Iyy        
    
    Iyy_theta = (Izz_0 + Iyy_0)/2. + (Izz_0 - Iyy_0)/2. * math.cos(2.*theta)
    Izz_theta = (Izz_0 + Iyy_0)/2. - (Izz_0 - Iyy_0)/2. * math.cos(2.*theta)
    Izy_theta = (Izz_0 - Iyy_0)/2. *math.sin(2.*theta)
    
    return Iyy_0, Izz_0, Iyy_theta, Izz_theta, Izy_theta

def ExactMOIdiscretisation(q,ndis,l1,l2,l3,l4,t_sk,t_sp,t_st,w_st,h_st,zcg,n,spacing,nodepos,xa,Ca,ha,theta,zsc):


    I = np.array([[],[]],[[],[]])

    xt = np.array([0])


    i = 1

    for x in np.linspace(0,l1,ndis+1)[1:]:
        xt = np.append(xt, x)
        dx = l1/(n)


        Iyy, Izz, Izy = ExactMOI(theta[i-1],Ca,ha,t_sk,t_sp,t_st,w_st,h_st,zcg,n,spacing,nodepos)[2:]
        I = np.append(I, [[Iyy,Izy],[Izy, Izz]])
        


        i = i+1
    for x in np.linspace(l1,l2-xa/2,ndis+1)[1:]:
        xt = np.append(xt, x)
        dx = (l2-l1)/(n)
        
        Iyy, Izz, Izy = ExactMOI(theta[i-1],Ca,ha,t_sk,t_sp,t_st,w_st,h_st,zcg,n,spacing,nodepos)[2:]
        I = np.append(I, [[Iyy,Izy],[Izy, Izz]])

        i = i+1

    for x in np.linspace(l2-xa/2,l2,ndis+1)[1:]:
        xt = np.append(xt, x)
        dx = (l2-l1)/(n)
        
       
        Iyy, Izz, Izy = ExactMOI(theta[i-1],Ca,ha,t_sk,t_sp,t_st,w_st,h_st,zcg,n,spacing,nodepos)[2:]
        I = np.append(I, [[Iyy,Izy],[Izy, Izz]])

        i = i+1
    for x in np.linspace(l2,l2+xa/2,ndis+1)[1:]:
        xt = np.append(xt, x)
        dx = (l3-l2)/(n)
        
        Iyy, Izz, Izy = ExactMOI(theta[i-1],Ca,ha,t_sk,t_sp,t_st,w_st,h_st,zcg,n,spacing,nodepos)[2:]
        I = np.append(I, [[Iyy,Izy],[Izy, Izz]])
        
        i = i+1
    for x in np.linspace(l2+xa/2,l3,ndis+1)[1:]:
        xt = np.append(xt, x)
        dx = (l3-l2)/(n)
        
        Iyy, Izz, Izy = ExactMOI(theta[i-1],Ca,ha,t_sk,t_sp,t_st,w_st,h_st,zcg,n,spacing,nodepos)[2:]
        I = np.append(I, [[Iyy,Izy],[Izy, Izz]])
        
        i = i+1
    for x in np.linspace(l3,l4,ndis+1)[1:]:
        xt = np.append(xt, x)
        dx = (l4-l3)/(n)

        Iyy, Izz, Izy = ExactMOI(theta[i-1],Ca,ha,t_sk,t_sp,t_st,w_st,h_st,zcg,n,spacing,nodepos)[2:]
        I = np.append(I, [[Iyy,Izy],[Izy, Izz]])
       

        i = i+1

        return I

def shear_flow_finder(boom_area_inclskin, Izz, Iyy, theta, node_pos, Ca, ha, Mx, Vy, Vz, G, tsk, tspar):
    #decompose forces to be local
    Vy = Vy*math.cos(theta)+Vz*math.sin(theta)
    Vz = Vz*math.cos(theta)-Vy*math.sin(theta)
    #counter-clockwise movement
    baes = boom_area_inclskin
    bxyz = node_pos
    #define the order of movement for the triangular and the circular section
    tring_booms = [12,13,8,9,10,11,1,2,3,4]
    circ_booms = [13,12,5,6,7]
    #circle section is I and triangular section is II
    #find areas of sections if not found
    s=math.sqrt((ha/2)**2+(Ca-ha/2)**2)
    peri=(ha/2)*math.pi
    AI=0.5*math.pi*(ha/2)**2
    AII=0.5*(Ca-ha/2)*(ha/2)
    #define distances between boom and next one and associated thicknesses
    spac = 0.1015545
    edge = 0.0766246
    tring_dist = [ha, edge, spac, spac, spac, spac, spac, spac, spac, edge]
    circ_dist = [ha,spac-edge, spac, spac, spac-edge]
    tring_thicc = [tspar, tsk, tsk, tsk, tsk, tsk, tsk, tsk, tsk, tsk]
    circ_thicc = [tspar, tsk, tsk, tsk, tsk]
    #find base shear flow for each cell
    tring_q = [0]
    circ_q = [0]
    for i in tring_booms:
        tring_q.append((-Vy/Izz)*baes[i]*bxyz[i][1]+(-Vz/Iyy)*baes[i]*bxyz[i][2]+tring_q[-1])
    for j in circ_booms:
        circ_q.append((-Vy/Izz)*baes[j]*bxyz[j][1]+(-Vz/Iyy)*baes[i]*bxyz[i][2]+circ_q[-1])
    
    #find force produced by each boom due to shear flows
    tring_fz=[]
    tring_fy=[]
    circ_fz= []
    circ_fy= []
    for i in range (len(tring_booms)):
        if i == len(tring_booms)-1:
            tring_fz.append(tring_q[i+1]*(bxyz[tring_booms[0]][2]-bxyz[tring_booms[i]][2]))
        else:
            tring_fz.append(tring_q[i+1]*(bxyz[tring_booms[i+1]][2]-bxyz[tring_booms[i]][2]))

    for i in range (len(tring_booms)):
        if i == len(tring_booms)-1:
            tring_fy.append(tring_q[i+1]*(bxyz[tring_booms[0]][1]-bxyz[tring_booms[i]][1]))
        else:
            tring_fy.append(tring_q[i+1]*(bxyz[tring_booms[i+1]][1]-bxyz[tring_booms[i]][1]))

    for i in range (len(circ_booms)):
        if i == len(circ_booms)-1:
            circ_fz.append(circ_q[i+1]*(bxyz[circ_booms[0]][2]-bxyz[circ_booms[i]][2]))
        else:
            circ_fz.append(circ_q[i+1]*(bxyz[circ_booms[i+1]][2]-bxyz[circ_booms[i]][2]))

    for i in range (len(circ_booms)):
        if i == len(circ_booms)-1:
            circ_fy.append(circ_q[i+1]*(bxyz[circ_booms[0]][1]-bxyz[circ_booms[i]][1]))
        else:
            circ_fy.append(circ_q[i+1]*(bxyz[circ_booms[i+1]][1]-bxyz[circ_booms[i]][1]))

    
    #find moment due to force
    #counter-clockwise positive
    moments=0
    for i in range(len(tring_fz)):
        moments += tring_fz[i]*bxyz[tring_booms[i]][1]

    for i in range(len(tring_fy)):
        moments += tring_fy[i]*(-1)*bxyz[tring_booms[i]][2]

    for i in range(len(circ_fz)):
        moments += circ_fz[i]*bxyz[tring_booms[i]][1]          

    for i in range(len(circ_fy)):
        moments += circ_fy[i]*(-1)*bxyz[tring_booms[i]][2] 
    
    #find line integral of (qbi*ds)/(t*G)
    tring_li=0
    circ_li=0
    for i in range(len(tring_dist)):
        tring_li += (tring_q[i+1]*tring_dist[i])/(tring_thicc[i]*G)

    for j in range(len(circ_dist)):
        circ_li += (circ_q[j+1]*circ_dist[i])/(circ_thicc[i]*G)
    #set up matrix
    A=np.matrix([[1/(2*AI)*(peri/(tsk*G)+ha/(tspar*G)), -ha/(2*AI*tspar*G), -1], [-ha/(2*AII*tspar*G), 1/(2*AII)*(2*s/(tsk*G)+ha/(tspar*G)), -1], [2*AI, 2*AII, 0]])
    b=np.matrix([[1/(-2*AI)*circ_li], [1/(-2*AII)*tring_li], [Mx-moments]])
    #solve matrix for redundant shear flows and rate of twist
    x = np.linalg.solve(A,b)
    qs0I = x.item(0)
    qs0II = x.item(1)
    twist_rate = x.item(2)
    #define order of output shear flows
    circoo=[2,3,4,5,1]
    tringoo=[7,8, 9,10,1,2,3,4,5,6]
    #find total shear flows
    circ_qt=[]
    tring_qt=[]    
    for i in circoo:
        circ_qt.append(circ_q[i]+qs0I)
    for i in tringoo:
        tring_qt.append(tring_q[i]+qs0II)
    
    return twist_rate, circ_qt, tring_qt


        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
