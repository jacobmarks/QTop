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
from src import color_codes, error_models, visualization
sys.path.append('decoders/')
from dsp_test import *
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


##################   Testing ##################


L, d, p = 5, 2, 0.04

# d = (2.5, 2.598)
# code.Primal.node[(4.0, 3.464)]['charge']['Z'] = 5
# code.Primal.node[(5.0, 3.464)]['charge']['Z'] = 5

code = color_codes.Color_6_6_6(L,d)
# print code.Primal.nodes()
# d = (5.5, 2.598)
# code.Primal.node[d]['charge']['Z'] = 1

model = error_models.CodeCapacity()
code = code.CodeCycle(model, p)
visualization.PlotPlaquette(code, "Before Decoding", 1)

decoder = DSP_decoder()
code = decoder(code)
if code.hasLogicalError():
	print "ERROR"
else:
	print "GOOD JOB!"
# visualization.PlotPrimal(code, "Bound Data", 2)
visualization.PlotPlaquette(code, "After Decoding", 3)
plt.show()




