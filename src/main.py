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
sys.path.insert(0, 'Decoders/')
from common import *
from surface_codes import *
from error_models import *
from decoders import *
from visualization import *
from threshold import *
from simulation import *



################## Surface Code Simulation ##################

L = 8
d = 2
p = .1
model = CodeCapacity()
decoder = HDRG_decoder()
sim = simulation(2, model, decoder)

sim(L,p)
# run(sim, L_array, p_array, num_trials)


