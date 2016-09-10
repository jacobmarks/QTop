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
	qubit.charge['Z'] = (qubit.charge['Z'] + randint(1,dim - 1))%dim
	return qubit 

def pauli_Y(qubit, dim):
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

# Function handle: for perfect gates, also performs gate with 0 errors
PerfectGate = lambda p: [0,0,0]

# For imperfect Gates
def error_probs(n):
	return lambda p: [float(p)/n, float(p)/n, float(p)/n]

class ErrorModel:

	def __init__(self, **kwargs):

		prop_defaults = {}

		for gate in ['initialize','identity','fourier','measure']:
			prop_defaults[gate] = PerfectGate
		prop_defaults['sum'] = {'target':PerfectGate, 'control':PerfectGate}

		
		for (prop, default) in prop_defaults.iteritems():
			setattr(self, prop, kwargs.get(prop, default))

	def Initialize(self, code, type, p):
		dimension = code.dimension
		for measure_position in code.stabilizers[type]:
			code.syndrome[measure_position].charge = Charge()
			measure_qubit = code.syndrome[measure_position]
			code.syndrome[measure_position] = BP_Channel(measure_qubit, dimension, self.initialize(p))
		return code

	def Identity(self, code, p):
		dimension = code.dimension
		for data_position in code.data:
			data_qubit = code.data[data_position]
			code.data[data_position] = BP_Channel(data_qubit, dimension, self.identity(p))
		return code

	def Fourier(self, code, type, p):
		dimension = code.dimension
		for measure_position in code.stabilizers[type]:
			measure_qubit = code.syndrome[measure_position]
			code.syndrome[measure_position] = BP_Channel(measure_qubit, dimension, self.fourier(p))
		return code

	def Measure(self, code, type, p):
		dimension = code.dimension
		for measure_position in code.stabilizers[type]:
			measure_qubit = code.syndrome[measure_position]
			code.syndrome[measure_position] = BP_Channel(measure_qubit, dimension, self.measure(p))
		return code

	def Sum(self, code, count, type, charge_type, p):
		num_sides = code.types[type]['num_sides']
		if count in range(num_sides/2):
			sign = 1 # Add control to target
		else:
			sign = -1 # subtract control from target

		dimension = code.dimension

		for target_position in code.stabilizers[type]:
			stabilizer = code.stabilizers[type][target_position]
			target_qubit = code.syndrome[target_position]
			if count in stabilizer.order:
				control_position = stabilizer.order[count]
				control_qubit = code.data[control_position]
				control_charge = control_qubit.charge[charge_type]
				target_charge = target_qubit.charge[charge_type]

				code.syndrome[target_position].charge[charge_type] = (target_charge + sign * control_charge)%dimension
				code.data[control_position] = BP_Channel(control_qubit, dimension, self.sum['control'](p))
			code.syndrome[target_position] = BP_Channel(target_qubit, dimension, self.sum['target'](p))
		return code


class CodeCapacity(ErrorModel):
	def __init__(self):
		ErrorModel.__init__(self)
		self.identity = error_probs(3)

class Phenomenological(ErrorModel):
	def __init__(self):
		ErrorModel.__init__(self)
		self.identity = error_probs(3)
		self.measure = error_probs(1)

class CircuitLevel(ErrorModel):
	def __init__(self):
		ErrorModel.__init__(self)
		self.initialize = error_probs(3)
		self.identity = error_probs(3)
		self.fourier = error_probs(3)
		self.measure = error_probs(1)
		self.sum = {'target':error_probs(4), 'control': error_probs(4)}




