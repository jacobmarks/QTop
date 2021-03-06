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
import time
import random
import signal
from common import *
from error_models import *
from surface_codes import *
from color_codes import *
from decoders import decoders

class TimeoutException(Exception):
	pass

def timeout_handler(signum, frame):
	raise TimeoutException


class simulation:

	def __init__(self, dimension, code_type, model, decoder, path_to):
		self.model = model
		[self.decoder, self.decoder_type] = decoder
		[self.model, self.model_type] = model
		self.dimension = dimension
		self.code_type = code_type
		self.path = path_to

	def __call__(self, L, p):
		signal.signal(signal.SIGALRM, timeout_handler)
		signal.alarm(10)
		if self.code_type == 'Surface Code':
			code = SurfaceCode(L, self.dimension)
		if self.code_type == '6-6-6 Color Code':
			code = Color_6_6_6(L, self.dimension)
		code = code.CodeCycle(self.model, p)
		try:
			decoders.Decode(code, self.decoder)
			return code.Assessment()
		except TimeoutException:
			return self(L, p)
		except ValueError:
			return self(L, p)


def run(sim, L_vals, p_vals, num_trials):
	L_array, p_phys_array, p_error_array = [], [], []
	for L in L_vals:
		pL_error_vals = []
		ps = len(p_vals)
		for i in range(ps):
			p = p_vals[i]
			L_array.append(L)
			p_phys_array.append(p)
			errors = 0
			
			for t in range(num_trials):
				print 'L', L, 'p', i+1, '/', ps, 'trial ', t+1, '/', num_trials
				if not sim(L,p):
					errors += 1
			p_error = float(errors)/num_trials
			pL_error_vals.append(p_error)

		p_error_array += pL_error_vals


	p_success_array = [1 - p for p in p_error_array]

	d = {'p_succ': p_success_array, 'p_phys': p_phys_array, 'L':L_array}
	d['trials'] = num_trials
	d['dimension'] = sim.dimension
	d['code_type'], d['model_type'], d['decoder_type'] = sim.code_type, sim.model_type, sim.decoder_type

	current_date_time = time.strftime("%c")
	current_date_time = current_date_time.replace("  ", "_")
	current_date_time = current_date_time.replace(" ", "_")
	file_name = sim.path + current_date_time + ".pickle"
	with open(file_name, 'wb') as handle:
		pickle.dump(d, handle)
	print str(current_date_time)





			