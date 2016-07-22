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
 
from random import random, randint
from math import *

from common import *

# Triple [X,Y,Z] of Pauli Error Rates for each gate
# Fourier generalizes Hadamard to qudits
# SUM generalizes CNOT to qudits

def pauli_X(data, dim):
	data.charge['X'] = (data.charge['X'] + randint(1,dim))%dim
	return data

def pauli_Z(data, dim):
	data.charge['Z'] = (data.charge['Z'] + randint(1,dim))%dim
	return data 

def pauli_Y(data, dim):
	data.charge['X'] = (data.charge['X'] + randint(1,dim))%dim
	data.charge['Z'] = (data.charge['Z'] + randint(1,dim))%dim
	return data



# bit flip and phase flip channel
def BP_Channel(data, dim, error_rates):
	p_error = sum(error_rates)
	if random() < p_error:
		err_type = random()
		if err_type < float(error_rates[0])/p_error:
			data = pauli_X(data, dim)
		elif err_type > float(error_rates[0] + error_rates[1])/p_error:
			data = pauli_Y(data, dim)
		else:
			data = pauli_Z(data, dim)
	return data

class ErrorModel:

	def __init__(self, **kwargs):
		prop_defaults = {'initialize': [0,0,0],
		'identity': [0,0,0],
		'fourier': [0,0,0],
		'measure': [0,0,0],
		'sum': {'target':[0,0,0], 'control': [0,0,0]},
		}
		
		for (prop, default) in prop_defaults.iteritems():
			setattr(self, prop, kwargs.get(prop, default))

	def Initialize(self, code, type):
		dimension = code.dimension
		for measure in code.stabilizers[type]:
			code.stabilizers[type][measure] = BP_Channel(measure, dimension, self.initialize)

	def Identity(self, code):
		dimension = code.dimension
		for data in code.primal.nodes():
			code.primal.node[data] = BP_Channel(data, dimension, self.identity)

	def Fourier(self, code, type):
		dimension = code.dimension
		for measure in code.stabilizers[type]:
			code.stabilizers[type][measure] = BP_Channel(measure, dimension, self.fourier)

	def Measure(self, code, type):
		dimension = code.dimension
		for measure in code.stabilizers[type]:
			code.stabilizers[type][measure] = BP_Channel(measure, dimension, self.measure)

	def Sum(self, code, order, type, charge_type):
		return code
		dimension = code.dimension
		for target in code.stabilizers[type]:
			for count in range(code.types[type]['num_sides']):
				if count not in code.stabilizers[type][target].order:
					continue
				control = code.stabilizers[type][target].order[count]
				print control.position, target
				if count%2 ==0:
					sign = 1
					# Add control to target
				else:
					sign = -1
					# subtract control from target

				# control_charge = control.charge[charge_type]
				# print code.stabilizers[type][target].center
				# target_charge = code.stabilizers[type][target].center.charge[charge_type]
				# print control_charge, target_charge, sign
				# code.stabilizers[type][target].center.charge[charge_type] = (target_charge + sign * control_charge)%dimension
				# print code.stabilizers[type][target].center.charge[charge_type]
				# new_control = BP_Channel(control, dimension, self.sum['control'])
				# code.stabilizers[type][target].data[count].charge = new_control.charge
				# code.stabilizers[type][target].order[count].charge = code.stabilizers[type][target].data[count].charge
				# code.stabilizers[type][target].center = BP_Channel(target, dimension, self.sum['target'])

class CodeCapacity(ErrorModel):
	def __init__(self, data_error_rate):
		ErrorModel.__init__(self)
		p = float(data_error_rate)
		self.identity = [p/3,p/3,p/3]

class Phenomenological(ErrorModel):
	def __init__(self, physical_error_rate):
		ErrorModel.__init__(self)
		p = float(physical_error_rate)
		self.identity = [p/3,p/3,p/3]
		self.measure = [p, p, p]

class CircuitLevel(ErrorModel):
	def __init__(self, error_rate):
		ErrorModel.__init__(self)
		p = float(error_rate)
		self.initialize = [p/3,p/3,p/3]
		self.identity = [p/3,p/3,p/3]
		self.fourier = [p/3,p/3,p/3]
		self.measure = [p, p, p]
		self.sum = {'target':[p/4,p/4,p/4], 'control': [p/4,p/4,p/4]}
