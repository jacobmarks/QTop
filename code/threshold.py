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


# make class for fitting
# so that users can make their own fits if they want to


class scaling(object):

	def __init__(self, N = False, data = False):
		self.num_trials = N
		self.data = data
		self.threshold = False


	def fit(self):
		p_array, L_array, p_succ = [], [], []
		for L in self.data:
			for p in self.data[L]:
				L_array.append(L)
				p_array.append(p)
				p_succ.append(self.data[L][p])

		X = [p_array, L_array]
		self.params, _ = curve_fit(self.form, X, p_succ, maxfev=int(10e7))
		print self.params
		# self.threshold = self.params[4]
		self.threshold = self.params[4]
		return self

	def get_threshold(self):
		if not self.threshold:
			self = self.fit()

		return self.threshold

	# def form(self, X, a, b):
	def form(self, X, A, B, C, D, p_th, u, v):

		return A + B * (X[0] - p_th) * np.power(X[1], float(1)/v) + C * (X[0] - p_th) * np.power(X[1], float(1)/v) ** 2 + D * np.power(X[1], - float(1)/u)
		# return X[0] + a * X[1] + b
	def uncertainty(self, p_succ):
		N = self.num_trials
		return sqrt(float(p_succ*(1 - p_succ))/N )

	def plot(self, log_plot = False, save_name = 'threshold_plot.png'):
		if not self.threshold:
			self = self.fit()

		fig = plt.figure()
		for L in self.data:
			p_array, sim_array, sigma = [], [], []
			for p in self.data[L]:
				p_array.append(p)
				sim_prob = self.data[L][p]
				sim_array.append(sim_prob)
				plt.scatter(p, sim_prob)
				# sigma.append(self.uncertainty(sim_prob))
			fit_array = self.form([p_array, [L]], *self.params)
			# fit_array = self.form([p_array, [L]], self.params[0], self.params[1],self.params[2],self.params[3], self.params[4], self.params[5], self.params[6])
			plt.plot(p_array, fit_array, label= "fitted data for L = " + str(L))

		plt.axvline(x=self.threshold,color='k',ls='dashed', label = "threshold")
		plt.title("Physical Error Rate vs. Logical Success Rate")
		plt.xlabel("Physical Error Rate")
		plt.ylabel("Logical Success Rate")
		plt.legend(loc='upper right', prop={'size':9})
		plt.savefig(save_name)
		plt.show()




	def __call__(self, N = False, data = False):
		if N == False:
			if self.num_trials == False:
				print 'Error: give number of trials as input' 
				return
		else:
			self.num_trials = N

		if data == False:
			if self.num_trials == False:
				print 'Error: give data array as input' 
				return
		else:
			self.data = data

		self.plot()
		

