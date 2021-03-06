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
model = error_models.Depolarizing()
decoder = mwpm.MWPM_decoder()
sim = simulation.simulation(2, 'Surface Code', [model, 'Depolarizing Channel'], [decoder, 'MWPM'], path_to)
L_vals = [3,5,7,9,11,13]
p_vals = np.linspace(0.08,0.15,15)
num_trials = 30000
simulation.run(sim, L_vals, p_vals, num_trials)

