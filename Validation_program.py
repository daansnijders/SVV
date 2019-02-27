# -*- coding: utf-8 -*-
import scipy
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import Definitions

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
        
#Dividing which of the selected node-lists/lines belong to rib A, B, C, or D (still in lists)
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

#Cross reference the elements
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


"""VON MISSES CALCULATIONS"""

von_misses_stress = []

#Open and extract the von mise stresses per element from the following document:
with open ("./data/F100_SLC1.rpt") as g2:
    l = (li.rstrip().split() for li in g2)
    l = list(li for li in l if li)
    
    #Append the von misses stresses to the empty list
    i = 0
    for li in l:
        if li[0].isdigit() == True and i < 13045:
            von_misses_stress.append(li)
        i =+ 1
    
        
    von_misses_stress = np.array(von_misses_stress)
    von_misses_stress = von_misses_stress.astype(float)    
    #Now take the average of the last two columns: vm_stres on the outside and inside
    von_misses_stress = np.column_stack((von_misses_stress[:,0],((von_misses_stress[:,2]+von_misses_stress[:,3])/2)))
    von_misses_stress = von_misses_stress[np.argsort(von_misses_stress[:,0])]

k = 0
von_misses_stress_element = []

#Taking the average of the four elements surrounding a node (after checking if they match)
for i in range (0, 13045):
    if von_misses_stress[k][0] == von_misses_stress[k+3][0] and k < 13045:
        von_misses_stress_element.append((von_misses_stress[k][1]+von_misses_stress[k+1][1]+von_misses_stress[k+2][1]+von_misses_stress[k+3][1])/4)
        k =+ 4
        
#Two arryas with stresses and the associated element nr        
von_misses_stress_element = np.array(von_misses_stress_element)
von_misses_stress_element = np.column_stack((np.arange(1, len(von_misses_stress_element)+1), von_misses_stress_element))

"""PLOTTING THE AIRFOIL 3D:"""
#We have now used the node list to find the elements and their associated vm stresses,
#But will now set it again to make a plot 
x = [] 
z = []
y = []


j = 0

#Appending x, y, and z to the selected nodes
for j in range(len(node)):
    x.append([node[j][0]]+[node[j][1]])
    z.append([node[j][0]]+[node[j][3]])
    y.append([node[j][0]]+[node[j][2]])
    
    j += 1

#Appending the coordinates to the nodes lists.
node = np.array(node)
x,y,z = node[:,1], node[:,2], node[:,3]

display = plt.figure()
plt.title('Selected nodes and their resp. x,y,z.')
axis = plt.axes(projection='3d')
axis.scatter3D(x,-z,y, c=z, cmap='Blues');
axis.set_xlim3d(-200,2000)
axis.set_ylim3d(-500,500)
axis.set_zlim3d(-200,200)
plt.show()


"""CALC DEFLECTIONS"""
deflection = []
l = []

with open("./data/F100_ULC1.rpt") as g3:
    l = (regel.rstrip().split() for regel in g3)
    l = list(regel for regel in l if regel)

j = 0
#append the whole line if it starts with a digit
for regel in l:
    if l[j][0].isdigit() == True:
        deflection.append(l[j])
    j += 1
       
      
deflection = np.array(deflection)
deflection = deflection.astype(np.float)
#since we do not care about the RF:
deflection = deflection[:,[0,5,6,7,8]]
#Only selecting the datarows we need:
deflection = deflection[0:3254] 

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

location_y = np.array(location_y)
location_y = location_y.astype(np.float)
location_y = location_y[:,[0,5,6,7,8]]
location_y = location_y[0:3254]

#Strip down to nodes laying on y = 0.
y_n = np.argwhere(y == 0)
#Find out which nodes lay on the leading and trailing edge:
LEzi = np.argwhere(z == max(z))
TEzi = np.argwhere(z == min(z))
TEzn = np.take(node[:,0],TEzi)
LEzn = np.take(node[:,0],LEzi)


dypos_trailingedge = np.take(deflection[:,3],TEzi)
dypos_leadingedge = np.take(deflection[:,3],LEzi)
ypos_trailingedge = np.take(location_y[:,3],TEzi)
ypos_leadingedge = np.take(location_y[:,3],LEzi)
xpos_trailingedge = np.take(node[:,1],TEzi)
xpos_leadingedge = np.take(node[:,1],LEzi)

ypos_trailingedge = np.take(node[:,2],TEzi)
ypos_leadingedge = np.take(node[:,2],LEzi)
#Taking into account the 
ypos_trailingedge = (dypos_trailingedge - ypos_trailingedge)
ypos_leadingedge = (dypos_leadingedge - ypos_leadingedge)

dist = np.column_stack((xpos_trailingedge , ypos_trailingedge))
dist = dist[np.argsort(dist[:,0])]
distdy = np.column_stack((xpos_trailingedge , dypos_trailingedge))
distdy = distdy[np.argsort(distdy[:,0])]


"""PLOTTING SLOPE"""
#Finding (and later plotting) slope
validation_slope = []
validation_step = []
for j in range(dist[:,0].size):
    if j+1 < dist[:,1].size :
        s = dist[j+1,0] - dist[j,0]
        pitch = (dist[j+1,1] - dist[j,1])
        validation_slope.append(pitch)
        validation_step.append(dist[j,0])

#To plot the slope:
#plt.plot(validation_step, validation_slope , 'ro')
#plt.plot(validation_step, validation_slope , 'red')



"""Importing nummerical Data from TE and LE """

nodepos2, rot = Definitions.offset(zcg, theta, nodepos, v2, u2, xt)

#plt.plot(steps_num , slopes_num , 'ro') - NUM
#plt.plot(numx ,numdefl , 'ro',label = 'Numerical data') - NUM
#plt.plot(numx ,numdefl , 'k') - NUM



"""PLOTING DEFLECTIONS""" 
#Deflection of y (?) along x:
plt.plot(dist[:,0],dist[:,1], 'r',label = 'Validation Data' )
#plt.plot(distdy[:,0],distdy[:,1], 'g' )

#plt.plot(nodepos2[0][:,0],nodepos2[0][:,1])

plt.xlabel("X along span b (mm)")
plt.ylabel("Y-Deflection of Validated (red line) and Numerical (Blue)???")
plt.grid('on')
axes = plt.gca()

axes.set_xlim([-100,1700])
axes.set_ylim([210,250])
#plt.scatter(xpos_leadingedge , ypos_leadingedge , alpha=0.5)
#plt.legend(['Data validation','Data numerical model'])
plt.show()







     



"""MORE ON VONS MISE"""   
element_set = np.array(element_set)
element_hinges = element_set[0][1:],element_set[2][1:],element_set[4][1:]
element_hinges = np.array([i for sl in element_hinges for i in sl]) - 1
#Calculate the stresses on the hinges'elements 
stress_hinge_1 = (von_misses_stress_element[[int(element_hinges[0])],1]+von_misses_stress_element[[int(element_hinges[1])],1] + von_misses_stress_element[[int(element_hinges[2])],1] + von_misses_stress_element[[int(element_hinges[3])],1])/4
stress_hinge_2 = (von_misses_stress_element[[int(element_hinges[4])],1]+von_misses_stress_element[[int(element_hinges[5])],1] + von_misses_stress_element[[int(element_hinges[6])],1] +von_misses_stress_element[[int(element_hinges[7])],1])/4
stress_hinge_3 = (von_misses_stress_element[[int(element_hinges[8])],1]+von_misses_stress_element[[int(element_hinges[9])],1] + von_misses_stress_element[[int(element_hinges[10])],1] + von_misses_stress_element[[int(element_hinges[11])],1])/4
oppervlakte = ((max(x)-min(x))/(len(TEzn)-1))**2





"""INSET REACTION FORCES IN Y-DIR AT HINGES"""
yh1 = 1.
yh2 = 1.
yh3 = 1.

"""INSET MAX STRESSES AT RIBS """
ribas = 1.
ribbs = 1.
ribcs = 1
ribds = 1.

noderib_a_element = np.array([i for sl in noderib_a_element for i in sl]) - 1
noderib_a_element = noderib_a_element.astype(np.int)
stress_rib_a = np.take(von_misses_stress_element[:,1],noderib_a_element)

noderib_b_element = np.array([i for sl in noderib_b_element for i in sl]) - 1
noderib_b_element = noderib_b_element.astype(np.int)
stress_rib_b = np.take(von_misses_stress_element[:,1],noderib_b_element)

noderib_c_element = np.array([i for sl in noderib_c_element for i in sl]) - 1
noderib_c_element = noderib_c_element.astype(np.int)
stress_rib_c = np.take(von_misses_stress_element[:,1],noderib_c_element)

noderib_d_element = np.array([i for sl in noderib_d_element for i in sl]) - 1
noderib_d_element = noderib_d_element.astype(np.int)
stress_rib_d = np.take(von_misses_stress_element[:,1],noderib_d_element)



#print (stress_hinge_1*1000,stress_hinge_2*1000,stress_hinge_3*1000)
#print (stress_hinge_1*(4*oppervlakte),stress_hinge_2*(4*oppervlakte),stress_hinge_3*(4*oppervlakte))
#print ("ABS",(stress_hinge_1*(4*oppervlakte) -yh1), stress_hinge_2*(4*oppervlakte) -yh2, stress_hinge_3*(4*oppervlakte)- yh3)
#print ("Abs" ,(((stress_hinge_1*(4*oppervlakte) -yh1)/yh1))*100, ((stress_hinge_2*(4*oppervlakte) -yh2)/yh2)*100, ((stress_hinge_3*(4*oppervlakte)- yh3)/yh3)*100)
#print (max(stress_rib_a)*1000)
#print (max(stress_rib_b)*1000)
#print (max(stress_rib_c)*1000)
#print (max(stress_rib_d)*1000)
print (((max(stress_rib_a)*1000) - ribas) ,((max(stress_rib_b)*1000) -ribbs) ,((max(stress_rib_c)*1000)-ribcs) ,((max(stress_rib_d)*1000) -ribds))


