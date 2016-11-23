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

from math import *
import networkx as nx
import matplotlib.pyplot as plt
from surface_codes import *
from color_codes import *
from error_models import *
from decoders import *
from visualization import *


##################   Testing ##################


L, d, p = 13, 7, 0.0

code = Color_6_6_6(L,d)

for m in code.Stabilizers['red']:
	for data in code.Plaquette(m,'red'):
		code.Primal.node[data]['charge']['Z'] = 4
		break
	break



model = CodeCapacity()
code = code.CodeCycle(model, p)




decoder = GSP_decoder()
code = decoder(code)
PlotPlaquette(code, "Before Decoding", 1)
plt.show()