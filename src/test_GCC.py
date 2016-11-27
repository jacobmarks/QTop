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
from color_codes import *
from error_models import *
from decoders import gcc
from visualization import *


##################   Testing ##################


L, d, p = 13, 7, 0.0

code = Color_6_6_6(L,d)


for m in code.Stabilizers['red']:
	for d1 in code.Plaquette(m,'red'):
		code.Primal.node[d1]['charge']['Z'] = 5
		break
	break
	# for d2 in code.Plaquette(m,'red'):
	# 	if d2 != d1:
	# 		code.Primal.node[d2]['charge']['Z'] = 2
	# 		break
	# break



model = CodeCapacity()
code = code.CodeCycle(model, p)


PlotPlaquette(code, "Before Decoding", 1)

decoder = gcc.GCC_decoder()
code = decoder(code)
PlotPlaquette(code, "After Decoding", 2)
plt.show()