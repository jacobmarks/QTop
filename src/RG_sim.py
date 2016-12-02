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

from common import *
from surface_codes import *
from error_models import *
from decoders import rg
# from visualization import *
# from threshold import *
from simulation import *

################## Surface Code Simulation ##################


model = CodeCapacity()
decoder = rg.HDRG_decoder()
L_vals = [3,5,9,11,17]
p_vals = np.logspace(-1.8,-.6,20)
num_trials = 10000
dims = [2,3,5,13,25]

for d in dims:
	sim = simulation(d, 'Surface Code', [model, 'Code Capacity'], [decoder, 'RG'])
	run(sim, L_vals, p_vals, num_trials)

