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
import pickle
import sys
import time
import random
import sets
from common import *
from threshold import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

name = str(sys.argv[1])

with open('../data_runs/' + name + '.pickle', 'rb') as handle:
  d = pickle.load(handle)

succ, phys, L = d['p_succ'], d['p_phys'], d['L']
num_depths = len(set(L))
num_probs = len(L)/num_depths
phys_probs = phys[0:num_probs-1]

for i in range(num_depths):
	depth = L[i*num_probs]
	success_probs = succ[i*num_probs:(i+1)*num_probs - 1]
	error_probs = [1 - p for p in success_probs]
	plt.plot(phys_probs, error_probs, label=str(depth))


X = [phys,L]

param_bounds=([0,-np.inf,-np.inf,-np.inf,-np.inf,0,-np.inf],[1,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])
params, _ = curve_fit(form, X, succ, bounds = param_bounds, max_nfev = 10e7)
# max_nfev = 10e7
threshold = params[0]
title = "threshold = " + str(threshold)
plt.title(str(title))
plt.xlabel("Physical Error Rate")
plt.ylabel("Logical Error Rate")
plt.legend(loc='upper left')
plt.savefig('../plots/' + name + '.png')
plt.show()