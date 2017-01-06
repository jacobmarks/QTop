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
from dsp import *
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt


##################   Testing ##################


L, d, p = 13, 2, 0.08

code = color_codes.Color_6_6_6(L,d)
model = error_models.PhaseError()
code = code.CodeCycle(model, p)
visualization.PlotPlaquette(code, "Before Decoding", 1)

decoder = DSP_decoder()
code = decoder(code)
if code.hasLogicalError():
	print "ERROR"
else:
	print "GOOD JOB!"
visualization.PlotPlaquette(code, "After Decoding", 3)
plt.show()

