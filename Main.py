# -*- coding: utf-8 -*-
import math
import Definitions


ha = 0.161 #height of thickest part of airfoil
Ca = 0.505 #dist from trailing to leading edge
n = 11 #number of stringers
list_length = n+3 #number of idealised booms INCLUDING two added AND Ghost nodepos
tskin = 0.0011
tspar = 0.0024
t_stiff = 0.0012
h_stiff = 0.013
w_stiff = 0.017



spacing, Cr, alpharad, Ct = Definitions.boom_spacing(ha, Ca, n)

alphadeg = math.degrees(alpharad)

nodepos, arc, dist = Definitions.boom_location(spacing, Cr, alpharad, list_length, ha)

area_stiff = Definitions.area_stiff(t_stiff, h_stiff, w_stiff)




ycg, zcg = Definitions.centroid_nonidealized(tskin, ha, Ca, Ct, tspar, nodepos, area_stiff)

#print (boom_area_incl_skin)

Ixx, Iyy, Izz = Definitions.boom_inertia(list_length, nodepos, boom_area)


My = 0.1
Mz = 0.1


boom_area = Definitions.boom_area_updater(tskin, spacing, Mz, My, Izz, Iyy, area_stiff, zcg, nodepos, dist, arc, tspar, ha)

