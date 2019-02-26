import numpy as np
import matplotlib.pyplot as plt

n=400
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
q=3860.
P1=53585
P2=49200
Ca=0.505
inittwist = 30
theta=np.ones(6*n+1)*inittwist*(np.pi/180)
zsc=0.


tring_booms = [1,2,3,4,12,13,8,9,10,11]
circ_booms = [12,5,6,7,13]
alpharad=0.187
spacing=0.1015480896
list_length=14
ha=0.161
Cr=0.0805
#tring_qsum=[17,27,37,47,57,67,77,87,97,107]
#circ_qsum=[12,13,14,15,16]
tskin=0.0011
tsk=tskin
tspar=0.0024
Ca=0.505
ha=0.161
G=28000000000
theta=30
Izz0=4.7189e-6
Iyy0=4.72552e-5

Izz=3.6621e-5
Iyy=1.35357e-5
Mx=1664.25

Mz=-1.14161513e+04
My=3.74664695e+03

Vz=49200  #value to be checked 
Vy=3860   # value to be checked 

stiff_area=3.5999999999999994e-05
zcg=-0.09993502900006332
dist=0.43206 # WHAT VALUES IS THIS ?
arc=0.2528982086

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


# =============================================================================
# def shear_flow_torsion(T,A1,A2,arc,l,ha,G,t):
#     
#     # T= resultant torque applied to the cross section
#     # A= area cell
#     # arc= lenght of the leading edge semicircle
#     # ha= diameter of the leading edge semi circle
#     # l= lenght of the triangular section (from tip of the triangle to intersection skin-spar)
#     # G=shear modulus
#     # t= skin thickness 
#     
# 
#     A=np.matrix([[0,2*A1,2*A2],[-1,(arc+ha)/(2*A1*G*t),-ha/(2*A1*G*t)],[-1,-ha/(2*A2*G*t),(2*l+ha)/(2*A2*G*t)]])
#     b=np.matrix([[T],[0],[0]])
#     x = np.linalg.solve(A,b)
#     rate_twist=x.item(0)
#     q1=x.item(1) # shear flow due to torsion in cell 1
#     q2=x.item(2) # shear flow due to torsion in cell 2
#    
# 
#     
#     return rate_twist,q1,q2
# 
# =============================================================================

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

nodepos=boom_location(spacing, Cr, alpharad, list_length, ha)[0]

def boom_area_updater(tsk, spacing, Mz, My, Izz, Iyy, stiff_area, zcg, nodepos, dist, arc, tspar, ha, theta):
    
    My = My*math.cos(theta)+Mz*math.sin(theta)
    Mz = Mz*math.cos(theta)-My*math.sin(theta)
    
    a= -1
    b= (Mz*Iyy)/(My*Izz)
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


def shear_flow_finder(boom_area_inclskin, Izz, Iyy, theta, node_pos, Ca, ha, Mx, Vy, Vz, G, tskin, tspar):
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
    for i in circ_booms:
        circ_q.append((-Vy/Izz)*baes[i]*bxyz[i][1]+(-Vz/Iyy)*baes[i]*bxyz[i][2]+circ_q[-1])
    
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
    print (moments)
    #find line integral of (qbi*ds)/(t*G)
    tring_li=0
    circ_li=0
    for i in range(len(tring_dist)):
        tring_li += (tring_q[i+1]*tring_dist[i])/(tring_thicc[i]*G)

    for i in range(len(circ_dist)):
        circ_li += (circ_q[i+1]*circ_dist[i])/(circ_thicc[i]*G)
    #set up matrix
    A=np.matrix([[1/(2*AI)*(peri/(tsk*G)+ha/(tspar*G)), -ha/(2*AI*tspar*G), -1], [-ha/(2*AII*tspar*G), 1/(2*AII)*(2*s/(tsk*G)+ha/(tspar*G)), -1], [2*AI, 2*AII, 0]])
    b=np.matrix([[1/(-2*AI)*circ_li], [1/(-2*AII)*tring_li], [Mx-moments]])
    #solve matrix for redundant shear flows and rate of twist
    x = np.linalg.solve(A,b)
    qs0I = x.item(0)
    qs0II = x.item(1)
    print (qs0I,qs0II, 1234)
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

boom_area = boom_area_updater(tsk, spacing, Mz, My, Izz, Iyy, stiff_area, zcg, nodepos, dist, arc, tspar, ha, theta)
twist_rate, circ_qt, tring_qt =shear_flow_finder(boom_area, Izz, Iyy, theta, nodepos, Ca, ha, Mx, Vy, Vz, G, tskin, tspar)

























def overalltwist(T,A1,A2,arc,l,ha,xa,G,t,l1,l2,l3,l4,n,inittwist):

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
    a=0
    for r in theta:
        Theta=r-theta_elem
        theta[a]=Theta
        a+=1
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
        

    #theta = theta + (inittwist*np.pi/180-theta[2*n])

    return theta, rate_twist_lst,xt


xt=overalltwist(T,A1,A2,arc,l,ha,xa,G,t,l1,l2,l3,l4,n,inittwist)[2]
theta= overalltwist(T,A1,A2,arc,l,ha,xa,G,t,l1,l2,l3,l4,n,inittwist)[0]

#plt.axis([0.0,0.5,-0.02,0,0025])
#plt.xlim(0,0.37)
#plt.ylim(-0.001,0.0005)
plt.plot(xt,theta)
plt.show()

