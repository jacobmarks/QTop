 #
 # QTop
 #
 # Copyright (c) 2016 Jacob Marks (jacob.marks@yale.edu)
 #
 # This file is part of QTop.
 #
 # QTop is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.

import sys
sys.path.append('../')
from src import common, color_codes, error_models, visualization
sys.path.append('decoders/')
from gcc import *
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


##################   Testing ##################


L, d, p = 13, 2, 0.04

code = color_codes.Color_6_6_6(L,d)

model = error_models.CodeCapacity()
code = code.CodeCycle(model, p)
copy = code.Primal.copy()

visualization.PlotPlaquette(code, "Color Code Plaquettes",1)

decoder = GCC_decoder()
code = decoder(code)

if code.hasLogicalError():
	print "ERROR"
	print [node for node in copy.nodes() if copy.node[node]['charge']['Z']!= 0]
else:
	print "GOOD JOB!"

visualization.PlotPlaquette(code, "Logical Error", 2)
plt.show()





# datas = [(7.0, 3.464), (14.5, 6.062), (8.5, 4.33), (6.5, 9.526), (6.5, 6.062), (9.5, 11.258)]
# for dat in datas:
# 	code.Primal.node[dat]['charge']['Z'] = 1

# d = (2.5, 2.598)
# code.Primal.node[d]['charge']['Z'] = 2
# code.Primal.node[(4.0, 3.464)]['charge']['Z'] = 5
# code.Primal.node[(5.0, 3.464)]['charge']['Z'] = 5

# for m in code.Stabilizers['blue']:
# 	if len(code.Stabilizers['blue'][m]['data']) == 6:
# 		d0, d1, d2 = code.Stabilizers['blue'][m]['order'][0], code.Stabilizers['blue'][m]['order'][1], code.Stabilizers['blue'][m]['order'][2]
# 		code.Primal.node[d0]['charge']['Z'] = 5
# 		code.Primal.node[d1]['charge']['Z'] = 1
# 		code.Primal.node[d2]['charge']['Z'] = 3
# 		break

# for m in code.Stabilizers['red']:
# 	if len(code.Stabilizers['red'][m]['data']) == 6:
# 		d0, d1, d5 = code.Stabilizers['red'][m]['order'][0], code.Stabilizers['red'][m]['order'][1], code.Stabilizers['red'][m]['order'][5]
# 		code.Primal.node[d0]['charge']['Z'] = 2
# 		code.Primal.node[d1]['charge']['Z'] = 4
# 		code.Primal.node[d5]['charge']['Z'] = 3
# 		break
	
# for m in code.Stabilizers['green']:
# 	for d1 in code.Plaquette(m,'green'):
# 		if d1 in code.Primal.nodes():
# 			code.Primal.node[d1]['charge']['Z'] = 2
# 			break
# 	break

