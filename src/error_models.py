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

def pauli_X(qubit, dim):
	qubit.charge['X'] = (qubit.charge['X'] + randint(1,dim - 1))%dim
	return qubit

def pauli_Z(qubit, dim):
	# print 'Z', data.position
	qubit.charge['Z'] = (qubit.charge['Z'] + randint(1,dim - 1))%dim
	return qubit 

def pauli_Y(qubit, dim):
	# print 'Y', data.position
	qubit.charge['X'] = (qubit.charge['X'] + randint(1,dim - 1))%dim
	qubit.charge['Z'] = (qubit.charge['Z'] + randint(1,dim - 1))%dim
	return qubit



# bit flip and phase flip channel
def BP_Channel(qubit, dim, error_rates):
	p_error = sum(error_rates)
	if random() < p_error:
		err_type = random()
		if err_type < float(error_rates[0])/p_error:
			qubit = pauli_X(qubit, dim)
		elif err_type > float(error_rates[0] + error_rates[1])/p_error:
			qubit = pauli_Y(qubit, dim)
		else:
			qubit = pauli_Z(qubit, dim)
	return qubit

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
		for measure_position in code.stabilizers[type]:
			measure_qubit = code.syndromes[measure_position]
			code.syndromes[measure_position] = BP_Channel(measure_qubit, dimension, self.initialize)
		return code

	def Identity(self, code):
		dimension = code.dimension
		for data_position in code.data:
			data_qubit = code.data[data_position]
			code.data[data_position] = BP_Channel(data_qubit, dimension, self.identity)
		return code

	def Fourier(self, code, type):
		dimension = code.dimension
		for measure_position in code.stabilizers[type]:
			measure_qubit = code.syndromes[measure_position]
			code.syndromes[measure_position] = BP_Channel(measure_qubit, dimension, self.fourier)
		return code

	def Measure(self, code, type):
		dimension = code.dimension
		for measure_position in code.stabilizers[type]:
			measure_qubit = code.syndromes[measure_position]
			code.syndromes[measure_position] = BP_Channel(measure_qubit, dimension, self.measure)
		return code

	def Sum(self, code, count, type, charge_type):
		if count %2 == 0:
			sign = 1 # Add control to target
		else:
			sign = -1 # subtract control from target

		dimension = code.dimension

		for target_position in code.stabilizers[type]:
			stabilizer = code.stabilizers[type][target_position]
			target_qubit = code.syndromes[target_position]
			if count in stabilizer.order:
				control_position = stabilizer.order[count]
				control_qubit = code.data[control_position]
				control_charge = control_qubit.charge[charge_type]
				target_charge = target_qubit.charge[charge_type]

				code.syndromes[target_position].charge[charge_type] = (target_charge + sign * control_charge)%dimension
				code.data[control_position] = BP_Channel(control_qubit, dimension, self.sum['control'])
			code.syndromes[target_position] = BP_Channel(target_qubit, dimension, self.sum['target'])
		return code


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
