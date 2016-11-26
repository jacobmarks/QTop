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
from decoders import dsp
from visualization import *


##################   Testing ##################


L, d, p = 13, 2, 0.04

code = Color_6_6_6(L,d)
model = CodeCapacity()
code = code.CodeCycle(model, p)
# PlotPlaquette(code, "Before Decoding", 1)

decoder = dsp.DSP_decoder()
code = decoder(code)
PlotPrimal(code, "Bound Data", 1)
# PlotPlaquette(code, "After Decoding", 2)
plt.show()