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
from gcc3 import *
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


##################   Testing ##################


L, d, p = 13, 2, 0.04

code = color_codes.Color_6_6_6(L,d)


# for node in code.Primal.nodes():
# 	if node in [(7.0, 10.392), (11.5, 11.258), (9.5, 7.794), (14.5, 2.598), (8.5, 6.062), (17.5, 2.598), (12.5, 7.794), (9.5, 4.33), (14.0, 6.928), (5.5, 4.33), (11.0, 6.928)]:
# 		code.Primal.node[node]['charge']['Z'] = 1


# for node in code.Primal.nodes():
# 	if node in [(7.0, 3.464), (13.0, 8.66), (3.5, 2.598), (14.5, 7.794), (7.0, 10.392), (11.5, 9.526), (11.5, 4.33), (16.0, 3.464), (16.0, 5.196)]:
# 		code.Primal.node[node]['charge']['Z'] = 1


# for node in code.Primal.nodes():
# 	if node in [(7.0, 6.928), (11.5, 7.794), (11.0, 8.66), (8.0, 5.196), (11.0, 12.124), (8.5, 12.99), (9.5, 6.062)]:
# 		code.Primal.node[node]['charge']['Z'] = 1

# code.Primal.node[(11.5, 6.062)]['charge']['Z'] = 1
# (2.5, 2.598)
# code.Primal.node[(3.5, 2.598)]['charge']['Z'] = 1
model = error_models.CodeCapacity()
ERR = False
# visualization.PlotPlaquette(code, "Color Code Plaquettes",1)
# visualization.PlotPlaquette(code, "Logical Error", 2)

code = code.CodeCycle(model, p)


copy = code.Primal.copy()

visualization.PlotPlaquette(code, "Color Code Plaquettes",1)

decoder = GCC_decoder()
code = decoder(code)

ERR = code.hasLogicalError()

visualization.PlotPlaquette(code, "After Error", 2)

if code.hasLogicalError():
	print "ERROR"
	print [node for node in copy.nodes() if copy.node[node]['charge']['Z']!= 0]

plt.show()

# while not ERR:
# 	fig = plt.figure(1)
# 	fig.clear()
# 	fig = plt.figure(2)
# 	fig.clear()
# 	code = code.CodeCycle(model, p)


# 	copy = code.Primal.copy()

# 	visualization.PlotPlaquette(code, "Color Code Plaquettes",1)

# 	decoder = GCC_decoder()
# 	code = decoder(code)

# 	ERR = code.hasLogicalError()
	
# 	visualization.PlotPlaquette(code, "Logical Error", 2)

# print [node for node in copy.nodes() if copy.node[node]['charge']['Z']!= 0]
# plt.show()


