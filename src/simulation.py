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


class simulation:

	def __init__(self, code, model, decoder):
		self.code = code
		self.model = model
		self.decoder = decoder

	def __call__(self):

		# re-initialize code (i.e. set everything to ground state)

		# Apply quantum error channel from model
		self.code = self.code.CodeCycle(self.model)

		# Find syndrome through stabilizer measurements
		syndrome = syndrome(self.code)

		# Classical decoding
		self.code = Decode(self.code, syndrome, self.decoder)

		# Check for logical errors
		return Assessment(self.code)


def run(simulation, L_array, p_array, num_trials):

	data = {'values':{}, 'L':L_array, 'p':p_array}

	for L in L_array:
		data['values'][L] = {}
		for p in p_array:

			successes = 0

			for i in range(num_trials):
				if simulation():
					successes += 1
			data['values'][L][p] = float(successes)/num_trials

	fitting = scaling(num_trials, data)
	threshold = fitting.get_threshold()
	print threshold
	fitting()






			