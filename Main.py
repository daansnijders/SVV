# -*- coding: utf-8 -*-
import math
import Definitions

ha = 0.161 #height of thickest part of airfoil
Ca = 0.505 #dist from trailing to leading edge
n = 11 #number of stringers
list_length = n+3 #number of idealised booms INCLUDING two added AND Ghost nodepos


spacing, Cr, alpharad = Definitions.boom_spacing(ha, Ca, n)
nodepos = Definitions.boom_location(spacing, Cr, alpharad, list_length)