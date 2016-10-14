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

 # from common import *
# from geometry import *
from error_models import *
from kitaev import *
from color_codes import *
from threshold import *
from geometry import *
from common import *
from decoders import *

def GenerateCode(code_type, geometry, depth, dimension):

	if code_type == 'kitaev':
		code = KTC(depth, dimension)
	if code_type == '4_8_8':
		code = Toric_4_8_8(depth, dimension)
	if code_type == '6_6_6':
		code = Toric_6_6_6(depth, dimension)
	if code_type == '4_6_12':
		code = Toric_4_6_12(depth, dimension)
	return code



class simulation:

	def __init__(self, code_type, geometry, dimension, model, decoder):
		self.model = model
		self.decoder = decoder
		self.code_type = code_type
		self.geometry = geometry
		self.dimension = dimension

	def __call__(self, L, p):

		# re-initialize code (i.e. set everything to ground state)
		code = GenerateCode(self.code_type, self.geometry, L, self.dimension)
		# Apply quantum error channel from model
		code = code.CodeCycle(self.model, p)

		# Find syndrome through stabilizer measurements
		syn = syndrome(code)

		# Classical decoding
		code = self.decoder(code, syn)

		# Check for logical errors
		return Assessment(code)


def run(simulation, L_array, p_array, num_trials):

	data = {'values':{}, 'L':L_array, 'p':p_array}

	for L in L_array:
		data['values'][L] = {}
		for p in p_array:
			successes = 0

			for i in range(num_trials):
				print 'L:', L, 'p:', p, 'trial:', i
				if simulation(L, p):
					successes += 1
			data['values'][L][p] = float(successes)/num_trials

	fitting = scaling(num_trials, data)
	threshold = fitting.get_threshold()
	print threshold
	# fitting()






			