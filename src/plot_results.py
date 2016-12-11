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

name, path_from, path_to = str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3])

with open(path_from + name + '.pickle', 'rb') as handle:
  d = pickle.load(handle)

succ, phys, L = d['p_succ'], d['p_phys'], d['L']
num_depths = len(set(L))
num_probs = len(L)/num_depths
phys_probs = phys[0:num_probs-1]

for i in range(num_depths):
	depth = L[i*num_probs]
	success_probs = succ[i*num_probs:(i+1)*num_probs - 1]
	error_probs = [1 - p for p in success_probs]
	# plt.semilogy(phys_probs, error_probs, label=str(depth))
	# plt.loglog(phys_probs, error_probs, label=str(depth))
	plt.plot(phys_probs, error_probs, label=str(depth))


X = [phys,L]

param_bounds=([0,-np.inf,-np.inf,-np.inf,-np.inf,0,-np.inf],[1,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])
params, pcov = curve_fit(form, X, succ, bounds = param_bounds, max_nfev = 10e7)

error = []
for i in range(len(params)):
	try:
		error.append(np.absolute(pcov[i][i])**0.5)
	except:
		error.append( 0.00 )

p_err = np.array(error)
threshold, threshold_uncert = params[0], p_err[0]
code, decoder, model, trials, dim = d['code_type'], d['decoder_type'], d['model_type'], d['trials'], d['dimension']
# title = "threshold = " + str(round(threshold, 3)) + "$\pm$" + str(round(threshold_uncert, 3))
print "threshold = " + str(round(threshold, 3)) + "$\pm$" + str(round(threshold_uncert, 3))
# plt.axvline(x=threshold, linewidth=2, color='k', ls = 'dashed', label='threshold')
title = "d = " + str(dim) + " "+ str(code) + " under " + str(model) + " Error Model"
plt.title(str(title))
plt.xlabel("Physical Error Rate")
plt.ylabel("Logical Error Rate")
plt.legend(loc='upper left', title = "Code Depths")
plt.savefig(path_to + name + '.png')
plt.show()