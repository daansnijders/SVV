import numpy as np
import matplotlib.pyplot as plt

n=10000
#T = np.ones(6*n+1)*3.76 # list containing the torque values at each spanwise location of the aileron 
A1=0.0101791529 # Area of cell 1 
A2=0.03417225 # Area of cell 2
arc=0.2528982086
l=0.43206
ha=0.161
G=28000000000
t=0.0011
#n=400.
l1=0.125
l2=0.498    
l3=1.494
l4=1.611
xa=0.245
q=3860
P1=53585
P2=49200
Ca=0.505
theta=np.ones(6*n+1)*30*(np.pi/180)
zsc=0.

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

T,xt=torque(q,n,l1,l2,l3,l4,P1,P2,xa,Ca,ha,theta,zsc) #torque vector!


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

def overalltwist(T,A1,A2,arc,l,ha,G,t,l1,l2,l3,l4,n):

    xt=np.array([0])
    theta = np.array([0])
    rate_twist_lst= np.array([])
   
    i = 1
    for x in np.linspace(0,l1,n+1)[1:]:
        xt=np.append(xt,x)
        dx = l1/(n)
        rate_twist=shear_flow_torsion(T[i-1],A1,A2,arc,l,ha,G,t)[0]
        rate_twist_lst = np.append(rate_twist_lst,rate_twist)
        theta_elem=theta[-1]+rate_twist_lst[i-1]*dx
        theta = np.append(theta,theta_elem)
        i = i+1
    
           
    for x in np.linspace(l1,l2-xa/2,n+1)[1:]:
        xt=np.append(xt,x)
        dx = (l2-l1)/(n)
        rate_twist=shear_flow_torsion(T[i-1],A1,A2,arc,l,ha,G,t)[0]
        rate_twist_lst = np.append(rate_twist_lst,rate_twist)
        theta_elem=theta[-1]+rate_twist_lst[i-1]*dx
        theta = np.append(theta,theta_elem)
        i = i+1
    
   # print(len(theta))
   # print(theta[-1])
    #fuck=theta_elem
    a=0
    for r in theta:
        Theta=r-theta_elem
        theta[a]=Theta
        a+=1
    #print(len(theta))
    #print(theta[-1])
       
    for x in np.linspace(l2-xa/2,l2,n+1)[1:]:
        xt=np.append(xt,x)
        dx = (l2-l1)/(n)
        rate_twist=shear_flow_torsion(T[i-1],A1,A2,arc,l,ha,G,t)[0]
        rate_twist_lst = np.append(rate_twist_lst,rate_twist)
        theta_elem=theta[-1]+rate_twist_lst[i-1]*dx
        theta = np.append(theta,theta_elem)
        i = i+1
    
        
        
    
    for x in np.linspace(l2,l2+xa/2,n+1)[1:]:
        xt=np.append(xt,x)
        dx = (l3-l2)/(n)
        rate_twist=shear_flow_torsion(T[i-1],A1,A2,arc,l,ha,G,t)[0]
        rate_twist_lst = np.append(rate_twist_lst,rate_twist)
        theta_elem=theta[-1]+rate_twist_lst[i-1]*dx
        theta = np.append(theta,theta_elem)
        i = i+1
    
    for x in np.linspace(l2+xa/2,l3,n+1)[1:]:
        xt=np.append(xt,x)
        dx = (l3-l2)/(n)
        rate_twist=shear_flow_torsion(T[i-1],A1,A2,arc,l,ha,G,t)[0]
        rate_twist_lst = np.append(rate_twist_lst,rate_twist)
        theta_elem=theta[-1]+rate_twist_lst[i-1]*dx
        theta = np.append(theta,theta_elem)
        i = i+1
    
        
        
    for x in np.linspace(l3,l4,n+1)[1:]:
        xt=np.append(xt,x)
        dx = (l4-l3)/(n)
        rate_twist=shear_flow_torsion(T[i-1],A1,A2,arc,l,ha,G,t)[0]
        rate_twist_lst = np.append(rate_twist_lst,rate_twist)
        theta_elem=theta[-1]+rate_twist_lst[i-1]*dx
        theta = np.append(theta,theta_elem)
        i = i+1
        
        #for i in range(len(theta)):
            #theta[i]=theta[i]-fuck

    return theta, rate_twist_lst,xt
xt=overalltwist(l1,l2,l3,l4,n)[2]
theta= overalltwist(l1,l2,l3,l4,n)[0]

#plt.axis([0.0,0.5,-0.02,0,0025])
#plt.xlim(0,0.37)
#plt.ylim(-0.001,0.0005)
plt.plot(xt,theta)
plt.show()

