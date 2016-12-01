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

def pauli_X(charge, dim):
	charge['X'] = (charge['X'] + randint(1,dim - 1))%dim
	return charge

def pauli_Z(charge, dim):
	charge['Z'] = (charge['Z'] + randint(1,dim - 1))%dim
	return charge 

def pauli_Y(charge, dim):
	charge['X'] = (charge['X'] + randint(1,dim - 1))%dim
	charge['Z'] = (charge['Z'] + randint(1,dim - 1))%dim
	return charge



# bit flip and phase flip channel
def BP_Channel(charge, dim, error_rates):
	p_error = sum(error_rates)
	if random() < p_error:
		err_type = random()
		if err_type < float(error_rates[0])/p_error:
			charge = pauli_X(charge, dim)
		elif err_type > float(error_rates[0] + error_rates[1])/p_error:
			charge = pauli_Y(charge, dim)
		else:
			charge = pauli_Z(charge, dim)
	return charge

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
		dim = code.dimension
		for measure_qubit in code.Stabilizers[type]:
			if measure_qubit not in code.External[type]:
				charge = code.Stabilizers[type][measure_qubit]['charge']
				code.Stabilizers[type][measure_qubit]['charge'] = BP_Channel(charge, dim, self.initialize(p))
		return code

	def Identity(self, code, p):
		dim = code.dimension
		for qubit in code.Primal.nodes():
			charge = code.Primal.node[qubit]['charge']
			code.Primal.node[qubit]['charge'] = BP_Channel(charge, dim, self.identity(p))
		return code

	def Fourier(self, code, type, p):
		dim = code.dimension
		for measure_qubit in code.Stabilizers[type]:
			if measure_qubit not in code.External[type]:
				charge = code.Stabilizers[type][measure_qubit]['charge']
				code.Stabilizers[type][measure_qubit]['charge'] = BP_Channel(charge, dim, self.initialize(p))
		return code

	def Measure(self, code, type, p):
		dim = code.dimension
		for measure_qubit in code.Stabilizers[type]:
			if measure_qubit not in code.External[type]:
				charge = code.Stabilizers[type][measure_qubit]['charge']
				code.Stabilizers[type][measure_qubit]['charge'] = BP_Channel(charge, dim, self.measure(p))
		return code

	def Sum(self, code, count, num_sides, type, charge_type, p):
		dim = code.dimension

		# sign = Sign(count, num_sides)
		sign = code.Sign(count, num_sides)


		for measure_qubit in code.Stabilizers[type]:
			if measure_qubit not in code.External[type]:
				measure_charge = code.Stabilizers[type][measure_qubit]['charge']
				if count in code.Stabilizers[type][measure_qubit]['order']:
					data_qubit = code.Stabilizers[type][measure_qubit]['order'][count]
					data_charge = code.Primal.node[data_qubit]['charge']
					code.Stabilizers[type][measure_qubit]['charge'][charge_type] = (measure_charge[charge_type] + sign * data_charge[charge_type])%dim
					code.Primal.node[data_qubit]['charge'] = BP_Channel(data_charge, dim, self.sum['control'](p))
				code.Stabilizers[type][measure_qubit]['charge'] = BP_Channel(measure_charge, dim, self.sum['target'](p))
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




