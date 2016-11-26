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

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from surface_codes import *
from decoders import rg
from error_models import *
from visualization import *

L, d, p = 7, 7, .1

code = SurfaceCode(L, d)
model = CodeCapacity()
code = code.CodeCycle(model, p)
PlotPlaquette(code, "Before Decoding", 1)

decoder = rg.HDRG_decoder()
code = decoder(code)
PlotPlaquette(code, "After Decoding", 2)
plt.show()

