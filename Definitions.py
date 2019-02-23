# -*- coding: utf-8 -*-
import math
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np


def deflection(q,n,r3,l1,l2,l3,l4,E,I,d1,d2,d3,rz3,P2,xa,Ca,ha,theta2):
    #theta2 is the angle of section 2 where hinge 2 is
    
    # calculation of reaction forces by static equalibrium
    r2 = (q*(l4)*(1.611/2-l1)-r3*(l3-l1))/(l2-l1)
    r1 = q*l4-r2-r3

    P1 = (r2*(0.25*Ca-ha/2)*np.cos(theta2))/(-ha/2*np.sin(theta2)-ha/2*np.cos(theta2))+P2
    rz2 = (-P1*(l2-l1-xa/2)+P2*(l2-l1+xa/2)-rz3*(l3-l1))/(l2-l1)
    rz1 = P2-P1-rz3-rz2

    #Moment, curvature, slope and deflection calculations
    
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

        Mz = np.append(Mz, q*(x**2)/2)
        My = np.append(My, 0)

        vdouble = np.append(vdouble,-(Mz[i]*I[0][0][i]-My[i]*I[0][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        vsingle = np.append(vsingle,vsingle[i-1]+vdouble[i]*dx)
        v = np.append(v, v[i-1] + vsingle[i-1]*dx + vdouble[i]*(dx**2)/2)

        udouble = np.append(udouble, -(-Mz[i]*I[0][1][i]-My[i]*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        usingle = np.append(usingle,usingle[i-1]+udouble[i]*dx)
        u = np.append(u, u[i-1] + usingle[i-1]*dx + udouble[i]*(dx**2)/2)

        i = i+1
    for x in np.linspace(l1,l2,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l2-l1)/(n)
        
        Mz = np.append(Mz, q*(x**2)/2-r1*(x-l1))
        My = np.append(My, rz1*(x-l1)+P1*(np.sign((x-l2+xa/2)) == 1)*(x-l2+xa/2))

        vdouble = np.append(vdouble,-(Mz[i]*I[0][0][i]-My[i]*I[0][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        vsingle = np.append(vsingle,vsingle[i-1]+vdouble[i]*dx)
        v = np.append(v, v[i-1] + vsingle[i-1]*dx + vdouble[i]*(dx**2)/2)

        udouble = np.append(udouble, -(-Mz[i]*I[0][1][i]-My[i]*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        usingle = np.append(usingle,usingle[i-1]+udouble[i]*dx)
        u = np.append(u, u[i-1] + usingle[i-1]*dx + udouble[i]*(dx**2)/2)

        i = i+1
    for x in np.linspace(l2,l3,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l3-l2)/(n)
        
        Mz = np.append(Mz, q*(x**2)/2-r1*(x-l1)-r2*(x-l2))
        My = np.append(My, rz1*(x-l1)+P1*(np.sign((x-l2+xa/2)) == 1)*(x-l2+xa/2)+rz2*(x-l2)-P2*(np.sign((x-l2-xa/2)) == 1)*(x-l2-xa/2))

        
        vdouble = np.append(vdouble,-(Mz[i]*I[0][0][i]-My[i]*I[0][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        vsingle = np.append(vsingle,vsingle[i-1]+vdouble[i]*dx)
        v = np.append(v, v[i-1] + vsingle[i-1]*dx + vdouble[i]*(dx**2)/2)

        udouble = np.append(udouble, -(-Mz[i]*I[0][1][i]-My[i]*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        usingle = np.append(usingle,usingle[i-1]+udouble[i]*dx)
        u = np.append(u, u[i-1] + usingle[i-1]*dx + udouble[i]*(dx**2)/2)
        
        i = i+1
    for x in np.linspace(l3,l4,n+1)[1:]:
        xt = np.append(xt, x)
        dx = (l4-l3)/(n)
        
        Mz = np.append(Mz, q*(x**2)/2-r1*(x-l1)-r2*(x-l2)-r3*(x-l3))
        My = np.append(My, 0)
        
        vdouble = np.append(vdouble,-(Mz[i]*I[0][0][i]-My[i]*I[0][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        vsingle = np.append(vsingle,vsingle[i-1]+vdouble[i]*dx)
        v = np.append(v, v[i-1] + vsingle[i-1]*dx + vdouble[i]*(dx**2)/2)

        udouble = np.append(udouble, -(-Mz[i]*I[0][1][i]-My[i]*I[1][1][i])/(E*(I[0][0][i]*I[1][1][i]-I[0][1][i]**2)))
        usingle = np.append(usingle,usingle[i-1]+udouble[i]*dx)
        u = np.append(u, u[i-1] + usingle[i-1]*dx + udouble[i]*(dx**2)/2)
        
        i = i+1
    
    deltav = ((v[2*n]-d2)*l1-(v[n]-d1)*l2)/(l1-l2)
    deltavsingle = ((v[n]-d1)-deltav)/l1

    v2 = -deltavsingle*xt-deltav+v

    deltau = ((u[2*n]-0)*l1-(u[n]-0)*l2)/(l1-l2)
    deltausingle = ((u[n]-0)-deltau)/l1
    u2 = -deltausingle*xt-deltau+u
            
    return [v2, u2, xt, r1, r2, My, Mz, rz1, rz2]

def bendingconvergence(q,n,r3a,r3b,l1,l2,l3,l4,E,I,d1,d2,d3,rz3a,rz3b,P2,xa,Ca,ha,theta2):
    #r3a and r3b relate to the initial reaction forces used to iterate bisection method
    #n is for the number of sections per section
    #I is an array that must follow the discretisation of xt, v2, thus len(I)=n*4+1

    #v2 is the array of vertical positions for a all points, where xt is the x position of the point

    v2r3a = deflection(q,n,r3a,l1,l2,l3,l4,E,I,d1,d2,d3,rz3a,P2,xa,Ca,ha,theta2)[0]

    v2r3b = deflection(q,n,r3b,l1,l2,l3,l4,E,I,d1,d2,d3,rz3b,P2,xa,Ca,ha,theta2)[0]

    u2rz3a = deflection(q,n,r3a,l1,l2,l3,l4,E,I,d1,d2,d3,rz3a,P2,xa,Ca,ha,theta2)[1]

    u2rz3b = deflection(q,n,r3b,l1,l2,l3,l4,E,I,d1,d2,d3,rz3b,P2,xa,Ca,ha,theta2)[1]  

    if v2r3a[3*n] > d3 or v2r3b[3*n] < d3:
        print('Please select a more apropriate y reaction force')
        
    if u2rz3a[3*n] > 0 or u2rz3b[3*n] < 0:
        print('Please select a more apropriate z reaction force')
        
    else:

        r3mid = (r3b+r3a)/2

        rz3mid = (rz3b+rz3a)/2

        v2, u2 = deflection(q,n,r3mid,l1,l2,l3,l4,E,I,d1,d2,d3,rz3mid,P2,xa,Ca,ha,theta2)[0:2]

        
        while round(v2[3*n],12) != d3 or round(u2[3*n],12) != 0:
            
            v2, u2, xt, r1, r2, My, Mz, rz1, rz2 = deflection(q,n,r3mid,l1,l2,l3,l4,E,I,d1,d2,d3,rz3mid,P2,xa,Ca,ha,theta2)

            if v2[3*n] < d3:
                r3a = r3mid
                r3b = r3b
            elif v2[3*n] > d3:
                r3a = r3a
                r3b = r3mid
            if u2[3*n] < 0:
                rz3a = rz3mid
                rz3b = rz3b
            elif u2[3*n] > 0:
                rz3a = rz3a
                rz3b = rz3mid
                
            r3mid = (r3b+r3a)/2
            rz3mid = (rz3b+rz3a)/2
            
    r3 = r3mid
    rz3 = rz3mid
    
    return v2, u2, xt, r1, r2, r3, rz1, rz2, rz3



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
            
        
    
    



        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
