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
from color_codes import *
from error_models import *
from decoders import gcc
from simulation import *

################## Surface Code Simulation ##################

model = CodeCapacity()
decoder = gcc.GCC_decoder()
L_vals = [9,13,17,21]
p_vals = np.logspace(-2.5,-1.5,10)
num_trials = 10000
d = 25
sim = simulation(d, '6-6-6 Color Code', [model, 'Code Capacity'], [decoder, 'GCC'])
run(sim, L_vals, p_vals, num_trials)

