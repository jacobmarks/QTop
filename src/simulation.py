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
from error_models import *
from surface_codes import *
from color_codes import *
from threshold import *
from decoders import decoders


class simulation:

	def __init__(self, dimension, model, decoder):
		self.model = model
		self.decoder = decoder
		self.dimension = dimension

	def __call__(self, L, p):

		code = SurfaceCode(L, self.dimension)
		code = code.CodeCycle(self.model, p)
		decoders.Decode(code, self.decoder)
		return code.Assessment()


def run(sim, L_vals, p_vals, num_trials):
	L_array, p_phys_array, p_error_array = [], [], []
	for L in L_vals:
		pL_error_vals = []
		for p in p_vals:
			L_array.append(L)
			p_phys_array.append(p)


			errors = 0

			for t in range(num_trials):
				print 'L', L, 'p', p, 'trial ', t+1, '/', num_trials
				if not sim(L,p):
					errors += 1
			p_error = float(errors)/num_trials
			pL_error_vals.append(p_error)

		p_error_array += pL_error_vals
		plt.plot(p_vals, pL_error_vals, label=str(L))


	p_success_array = [1 - p for p in p_error_array]

	X = [p_phys_array,L_array]
	param_bounds=([0,-np.inf,-np.inf,-np.inf,-np.inf,0,-np.inf],[1,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])


	params, _ = curve_fit(form, X, p_success_array, bounds = param_bounds)
	threshold = params[0]
	print threshold
	print params
	title = "threshold = " + str(threshold)
	plt.title(str(title))
	plt.xlabel("Physical Error Rate")
	plt.ylabel("Logical Error Rate")
	plt.legend(loc='upper left')
	plt.savefig('plot.png')
	plt.show()







			