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
from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

# performs approximate threshold by fitting logical success probability 
# as a function of physical error rate for each code depths. 

# inputs: 
# N: Num trials per data point
# data: 2d dict with p_succ for L, p combinations


# Fitting form used here is from
# C. Wang, J. Harrington, and J. Preskill, 2003

class scaling(object):

	def __init__(self, N = False, data = False):
		self.num_trials = N
		self.values = data['values']
		self.p = data['p']
		self.L = data['L']
		self.threshold = False


	def fit(self):
		p_array, L_array, p_succ = [], [], []
		for L in self.L:
			for p in self.p:
				L_array.append(L)
				p_array.append(p)
				p_succ.append(self.values[L][p])

		X = [p_array, L_array]
		print X, p_succ
		self.params, _ = curve_fit(self.form, X, p_succ, maxfev=int(10e7))
		self.threshold = self.params[1]
		print self.threshold
		return self

	def get_threshold(self):
		if not self.threshold:
			self = self.fit()
		return self.threshold

	def form(self, X, p_th, v):

		return (X[0] - p_th) * np.power(X[1], float(1)/v)
	def uncertainty(self, p_succ):
		N = self.num_trials
		return np.sqrt(float(p_succ*(1 - p_succ))/N )

	def plot(self, LOG = True, save_name = 'threshold_plot.png'):
		if not self.threshold:
			self = self.fit()

		fig = plt.figure()
		if LOG:
			ax = plt.subplot()
			ax.set_xscale("log", nonposx='clip')
			ax.set_yscale("log", nonposy='clip')

		for L in self.L:
			p_array, sim_array = [], []
			# if self.num_trials != False:
			# 	sigma = []
			for p in self.p:
				p_array.append(p)
				sim_prob = self.values[L][p]
				sim_array.append(sim_prob)
				plt.scatter(p, 1 - sim_prob)
				# if self.num_trials != False:
					# sigma.append(self.uncertainty(sim_prob))
			fit_array = self.form([p_array, [L]], *self.params)
			success_array = 1 - fit_array
			plt.plot(p_array, success_array, label= "fitted data for L = " + str(L))
			# if self.num_trials != False:
			# 	plt.errorbar(p_array, fit_array, yerr = sigma)

		
	
		plt.axvline(x=self.threshold,color='k',ls='dashed', label = "threshold = " + str(self.threshold))
		plt.title("Physical Error Rate vs. Logical Error Rate")
		plt.xlabel("Physical Error Rate")
		plt.ylabel("Logical Error Rate")
		plt.legend(loc='upper left', prop={'size':9})
		plt.savefig(save_name)
		plt.show()




	def __call__(self, N = False, data = False):
		# if N == False:
		# 	if self.num_trials == False:
		# 		print 'Error: give number of trials as input' 
		# 		return
		# else:
		self.num_trials = N

		if data == False:
			if self.num_trials == False:
				print 'Error: give data array as input' 
				return
		else:
			self.values, self.L, self.p = data['values'], data['L'], data['p']

		self.plot()
		
