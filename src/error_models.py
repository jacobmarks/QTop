 #
 # QTop
 #
 # Copyright (c) 2017 Jacob Marks (jacob.marks@yale.edu)
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



# Triple [X, Z, XZ] of Pauli Error Rates for each gate
# Fourier generalizes Hadamard to qudits
# SUM generalizes CNOT to qudits

def PauliX(charge, dim):
	charge['X'] = (charge['X'] + randint(1, dim - 1))%dim
	return charge

def PauliZ(charge, dim):
	charge['Z'] = (charge['Z'] + randint(1, dim - 1))%dim
	return charge 

# PauliXZ is correlated Pauli X and Z errors
def PauliXZ(charge, dim):
	r = randint(1,dim - 1)
	charge['X'] = (charge['X'] + r)%dim
	charge['Z'] = (charge['Z'] + r)%dim
	return charge


def Identity_Channel(charge):
	return charge

def PhaseFlip_Channel(charge, dim, p):
	if random() < p:
		charge = PauliZ(charge, dim)
	return charge

def BitFlip_Channel(charge, dim, p):
	if random() < p:
		charge = PauliX(charge, dim)
	return charge

def Correlated_Channel(charge, dim, p):
	if random() < p:
		charge = PauliXZ(charge, dim)
	return charge

# Independently Compose Bit Flip and Phase Flip Channels
def BP_Channel(charge, dim, p):
	return PhaseFlip_Channel(BitFlip_Channel(charge, dim, p), dim, p)

def Depolarizing_Channel(charge, dim, p):
	if random() < p:
		charge = Depolarizing_Helper(charge, dim)
	return charge

def Depolarizing_Helper(charge, dim):
	r1, r2 = randint(0,dim - 1), randint(0,dim - 1)
	if r1 + r2 == 0:
		return Depolarizing_Helper(charge, dim)

	c = charge
	c['X'] = (charge['X'] + r1)%dim
	c['Z'] = (charge['Z'] + r2)%dim
	return c


# Error Channel
def Error_Channel(charge, dim, p, type):
	if type == "perfect gate":
		return Identity_Channel(charge)
	if type == "bp":
		return BP_Channel(charge, dim, p)
	if type == "dp":
		return Depolarizing_Channel(charge, dim, p)
	if type == "bit flip":
		return BitFlip_Channel(charge, dim, p)
	if type == "phase flip":
		return PhaseFlip_Channel(charge, dim, p)
	if type == "correlated":
		return Correlated_Channel(charge, dim, p)
	else:
		return Identity_Channel(charge)



class ErrorModel:

	def __init__(self, **kwargs):

		prop_defaults = {}

		for gate in ['initialize','identity','fourier','measure']:
			prop_defaults[gate] = "perfect gate"
		prop_defaults['sum'] = {'target':"perfect gate", 'control':"perfect gate"}

		
		for (prop, default) in prop_defaults.iteritems():
			setattr(self, prop, kwargs.get(prop, default))

	def Initialize(self, code, type, p):
		dim = code.dimension
		for measure_qubit in code.Stabilizers[type]:
			if measure_qubit not in code.External[type]:
				charge = code.Stabilizers[type][measure_qubit]['charge']
				code.Stabilizers[type][measure_qubit]['charge'] = Error_Channel(charge, dim, p, self.initialize)
		return code

	def Identity(self, code, p):
		dim = code.dimension
		for qubit in code.Primal.nodes():
			charge = code.Primal.node[qubit]['charge']
			code.Primal.node[qubit]['charge'] = Error_Channel(charge, dim, p, self.identity)
		return code

	def Fourier(self, code, type, p):
		dim = code.dimension
		for measure_qubit in code.Stabilizers[type]:
			if measure_qubit not in code.External[type]:
				charge = code.Stabilizers[type][measure_qubit]['charge']
				code.Stabilizers[type][measure_qubit]['charge'] = Error_Channel(charge, dim, p, self.initialize)
		return code

	def Measure(self, code, type, p):
		dim = code.dimension
		for measure_qubit in code.Stabilizers[type]:
			if measure_qubit not in code.External[type]:
				charge = code.Stabilizers[type][measure_qubit]['charge']
				code.Stabilizers[type][measure_qubit]['charge'] = Error_Channel(charge, dim, p, self.measure)
		return code

	def Sum(self, code, count, num_sides, type, charge_type, control_p, target_p):
		dim = code.dimension
		sign = code.Sign(count, num_sides)


		for measure_qubit in code.Stabilizers[type]:
			if measure_qubit not in code.External[type]:
				measure_charge = code.Stabilizers[type][measure_qubit]['charge']
				if count in code.Stabilizers[type][measure_qubit]['order']:
					data_qubit = code.Stabilizers[type][measure_qubit]['order'][count]
					data_charge = code.Primal.node[data_qubit]['charge']
					code.Stabilizers[type][measure_qubit]['charge'][charge_type] = (measure_charge[charge_type] + sign * data_charge[charge_type])%dim
					code.Primal.node[data_qubit]['charge'] = Error_Channel(data_charge, dim, control_p, self.sum['control'])
				code.Stabilizers[type][measure_qubit]['charge'] = Error_Channel(measure_charge, dim, target_p, self.sum['target'])
		return code

class CodeCapacity(ErrorModel):
	def __init__(self):
		ErrorModel.__init__(self)
		self.identity = "bp"

class Depolarizing(ErrorModel):
	def __init__(self):
		ErrorModel.__init__(self)
		self.identity = "dp"

class Phenomenological(ErrorModel):
	def __init__(self):
		ErrorModel.__init__(self)
		self.identity = "bp"
		self.measure = "bp"
		self.sum['target'] = "bp"
		self.sum['control'] = "bp"

class CircuitLevel(ErrorModel):
	def __init__(self):
		ErrorModel.__init__(self)
		self.identity = "bp"
		self.measure = "bp"
		self.fourier = "bp"
		self.initialize = "bp"
		self.sum['target'] = "bp"
		self.sum['control'] = "bp"

class PhaseFlip(ErrorModel):
	def __init__(self):
		ErrorModel.__init__(self)
		self.identity = "phase flip"

class BitFlip(ErrorModel):
	def __init__(self):
		ErrorModel.__init__(self)
		self.identity = "bit flip"

