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
import numpy as np
import sys
sys.path.append('../')
from src import common, simulation, error_models
sys.path.append('../src')
from decoders import mwpm
################## Surface Code Simulation ##################

path_to = str(sys.argv[1])
model = error_models.CodeCapacity()
decoder = mwpm.MWPM_decoder()
sim = simulation.simulation(2, 'Surface Code', [model, 'Code Capacity'], [decoder, 'MWPM'], path_to)
L_vals = [3,5]
p_vals = np.logspace(-1.8,-.6,20)
num_trials = 10000
simulation.run(sim, L_vals, p_vals, num_trials)

