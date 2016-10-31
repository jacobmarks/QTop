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
from surface_codes import *
# from color_codes import *
from threshold import *
from common import *
from decoders import *


class simulation:

	def __init__(self, dimension, model, decoder):
		self.model = model
		self.decoder = decoder
		self.dimension = dimension

	def __call__(self, L, p):

		code = SurfaceCode(L, self.dimension)
		code = code.CodeCycle(self.model, p)
		Decode(code, self.decoder)
		print code.Assessment()
		return code.Assessment()


# def run(simulation, L_array, p_array, num_trials):

# 	data = {'values':{}, 'L':L_array, 'p':p_array}

# 	for L in L_array:
# 		data['values'][L] = {}
# 		for p in p_array:
# 			successes = 0

# 			for i in range(num_trials):
# 				print 'L:', L, 'p:', p, 'trial:', i
# 				if simulation(L, p):
# 					successes += 1
# 			data['values'][L][p] = float(successes)/num_trials

# 	fitting = scaling(num_trials, data)
# 	threshold = fitting.get_threshold()
# 	print threshold
# 	# fitting()






			