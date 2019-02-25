# -*- coding: utf-8 -*-
import scipy
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

node = []
element_set = []
node_set = []
j = 0

with open("./data/F100-19.inp") as g:
    for oneline in g:
        if oneline.startswith('*') == False:
            if j < 3244:    #Putting all the nodes in a list
                oneline = oneline.strip('\n')
                line_cont = [float(y) for y in oneline.split(',')]
                node.append(line_cont)
                
            elif j > 3243 and j < (3243+3251):
                #putting all the elements in a list
                oneline = oneline.strip('\n')
                line_cont = [float(y) for y in oneline.split(',')]
                element_set.append(line_cont)
             
            elif j > 7020 and j < (7041):
                #Assigning which nodes belong to Ribs
                oneline = oneline.strip(',\n')
                line_cont = [float(y) for y in oneline.split(',')]
                node_set.append(line_cont)                
        j = j + 1
        
#Dividing which of the selected nodes belong to rib A, B, C, or D (still in lists)
noderib_a = node_set[:4]
noderib_b = node_set[4:8]
noderib_c = node_set[8:12]
noderib_d = node_set[12:16]

#Convert the lists into one array (per rib) 
noderib_a = np.array([j for i in noderib_a for j in i]) - 1
noderib_b = np.array([j for i in noderib_b for j in i]) - 1
noderib_c = np.array([j for i in noderib_c for j in i]) - 1
noderib_d = np.array([j for i in noderib_d for j in i]) - 1

#Make empty lists of the elements to later link them to the nodes
noderib_a_element = []
noderib_b_element = []
noderib_c_element = []
noderib_d_element = []

for q in noderib_a:
    noderib_a_element.append(element_set[int(q)-1])
for q in noderib_b:
    noderib_b_element.append(element_set[int(q)-1])
for q in noderib_c:
    noderib_c_element.append(element_set[int(q)-1])
for q in noderib_d:
    noderib_d_element.append(element_set[int(q)-1])

noderib_a_element = np.array(noderib_a_element)
noderib_b_element = np.array(noderib_b_element)
noderib_c_element = np.array(noderib_c_element)
noderib_d_element = np.array(noderib_d_element)

von_misses_stress = []
with open ("./data/F100_SLC1.rpt") as g_i:
    l = (li.rstrip().split() for li in g_i)
    l = list(li for li in l if li)
    
    i = 0
    for li in l:
        if li[0].isdigit() == True and i < 13045:
            von_misses_stress.append(li)
        i =+ 1
        
    von_misses_stress = np.array(von_misses_stress)
    von_misses_stress = von_misses_stress.astype(float)    
    von_misses_stress = np.column_stack((von_misses_stress[:,0],((von_misses_stress[:,2]+von_misses_stress[:,3])/2)))
    von_misses_stress = von_misses_stress[np.argsort(von_misses_stress[:,0])]

#------------------------------------

k = 0
von_misses_stress_element = []
for i in range (0, 13045):
    if von_misses_stress[k][0] == von_misses_stress[k+3][0] and k < 13045:
        von_misses_stress_element.append((von_misses_stress[k][1]+von_misses_stress[k+1][1]+von_misses_stress[k+2][1]+von_misses_stress[k+3][1])/4)
        k =+ 4
von_misses_stress_element = np.array(von_misses_stress_element)
von_misses_stress_element = np.column_stack((np.arange(1, len(von_misses_stress_element)+1), von_misses_stress_element))


z = []
y = []
x = []
j = 0

for j in range(len(node)):
    z.append([node[j][0]]+[node[j][2]])
    y.append([node[j][0]]+[node[j][2]])
    x.append(node[j][:2])
    j =+ 1

node = np.array(node)
x,y,z = node[:,1], node[:,2], node[:,3]

display = plt.figure()
axis = plt.axes(projection='3d')
axis.scatter3D(x,-z,y, c=z, cmap='Greens');
axis.set_xlim3d(-200,2000)
axis.set_ylim3d(-500,500)
axis.set_zlim3d(-200,200)
plt.show()


#"""Open the F100_ULC1_rpt file """
deflection = []
l = []

with open("./data/F100_ULC1.rpt") as f_in:
    l = (line.rstrip().split() for line in f_in)
    l = list(line for line in l if line)

i = 0
for line in l:
    if l[i][0].isdigit() == True:
        deflection.append(l[i])
    i = i + 1
        
deflection = np.array(deflection).astype(np.float)
deflection = deflection[:,[0,5,6,7,8]]
#making an array with the node numbers and the deflectionections which correspons to it.
deflection = deflection[0:3254] #selecting the oldest dataset.


#now we need to sort the nodes for which the z and y coordinates are the same , with a varying x - coordinate
location_y = []
with open("./data/F100_UR1.rpt") as f_in:
    l = (li.rstrip().split() for li in f_in)
    l = list(li for li in l if li)
i = 0
for li in l:
    if l[i][0].isdigit() == True:
        location_y.append(l[i])
    i = i + 1
    
location_y = np.array(location_y).astype(np.float)
location_y = location_y[:,[0,5,6,7,8]]
location_y = location_y[0:3254]

#we can do this by checking at which indices y,z are the max or min valua
i = 0
yn = np.argwhere(y == 0)
LEzi = np.argwhere(z == max(z))
TEzi = np.argwhere(z == min(z))
TEzn = np.take(node[:,0],TEzi)
LEzn = np.take(node[:,0],LEzi)
#now that we know the indices where the trailing edge and leading edge are, we can extract the indices of x
#we now know all the nodes of the trailing and the leading edge.
#now for a 2d scatter plot we need to take the x values and the deflectionection in y values
dypos_trailingedge = np.take(deflection[:,3],TEzi)
dypos_leadingedge = np.take(deflection[:,3],LEzi)

ypos_trailingedge = np.take(location_y[:,3],TEzi)
ypos_leadingedge = np.take(location_y[:,3],LEzi)
xpos_trailingedge = np.take(node[:,1],TEzi)
xpos_leadingedge = np.take(node[:,1],LEzi)
ypos_trailingedge = np.take(node[:,2],TEzi)
ypos_leadingedge = np.take(node[:,2],LEzi)
ypos_trailingedge = (dypos_trailingedge - ypos_trailingedge)
ypos_leadingedge = (dypos_leadingedge - ypos_leadingedge)
dist = np.column_stack((xpos_trailingedge ,dypos_trailingedge))
dist = dist[np.argsort(dist[:,0])]
    
"""PLOTTING THE """

slopes_val = []
steps_val = []
for i in range(dist[:,0].size):
    if i+1 < dist[:,1].size :
        dx = dist[i+1,0] - dist[i,0]
        slope = (dist[i+1,1] - dist[i,1])
        slopes_val.append(slope)
        steps_val.append(dist[i,0])

#the code underneath is only used for plotting the slope graph. Otherwise the other lines forthe plotting is used.
#plt.plot(steps_num , slopes_num , 'ro')
#plt.plot(steps_val , slopes_val , 'bo')
#plt.plot(steps_val , slopes_val , 'blue ')


#to plot the deflection over the x-position
plt.plot(dist[:,0],dist[:,1], 'bo',label = 'Validation Data' )
#plt.plot(numx ,numdefl , 'ro',label = 'Numerical data')
plt.plot(dist[:,0],dist[:,1], 'k' )
#plt.plot(numx ,numdefl , 'k')

plt.xlabel("X-position along span (mm)")
plt.ylabel("Y-position (mm)")
plt.grid('on')
axes = plt.gca()

axes.set_xlim([-150,1800])
axes.set_ylim([150,600])
#plt.scatter(x_LE , dy_LE , alpha=0.5)
plt.legend(['Validation data','Numerical data'])
plt.show()




