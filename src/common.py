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
 
import networkx as nx
from random import random, randint
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def Charge(X_charge = 0, Z_charge = 0):
	return {'X':X_charge,'Z':Z_charge}

class Qubit:
	def __init__(self, position, charge = False, type = 'data'):
		self.position = position
		self.type = type
		if charge == False:
			self.charge = Charge()
		else:
			self.charge = charge
		

class Stabilizer:
	def __init__(self, type, data, order = False):
		self.data = data
		self.type = type
		if order == False:
			self.order = Order(data)
		else:
			self.order = order

def Order(data_list):
	count, order = 0, {}
	for data_position in data_list:
		order[count] = data_position
		count += 1
	return order

class Code:
	'Common base class for all types of topological codes'

	def __init__(self, depth, dimension = 2):
		self.depth = depth
		self.dimension = dimension
		self.data = {}
		self.syndromes = {}
		self.memberships = {}
		self.generateProperties()
		self.generateStabilizers()
		self.generatePrimal()
		self.generateDual()

	def generateStabilizerTypes(self):
		self.stabilizers = {}
		for type in self.types:
			self.stabilizers[type] = {}

	def generateProperties(self):
		self.generateTypes()
		self.generateColors()
		self.generateStabilizerTypes()

	def generateStabilizerData(self, measure_position, scale, num_sides, angle = 0):
		data = []
		for k in range(num_sides):
			x = scale * float(cos(2*pi*k/num_sides))
			y = scale * float(sin(2*pi*k/num_sides))
			if angle != 0:
				x_prime = float(cos(angle))*x - float(sin(angle))*y
				y_prime = float(sin(angle))*x + float(cos(angle))*y
				x, y = x_prime, y_prime
			position = (round(x + measure_position[0], 3) ,round(y  + measure_position[1], 3))
			data.append(position)
		return data

	def generatePrimalEdges(self, stabilizer):
		stabilizer_type = stabilizer.type
		num_sides = self.types[stabilizer_type]['num_sides']
		for k in range(num_sides):
			if k in stabilizer.order:
				vertex1 = stabilizer.order[k]
				while (k+1)%num_sides not in stabilizer.order:
					k += 1
				vertex2 = stabilizer.order[(k+1)%num_sides]
				self.primal.add_edge(*(vertex1,vertex2), color = 'black')


	def generatePrimal(self):
		self.primal = nx.Graph()
		for type in self.types:
			for measure_position in self.stabilizers[type]:
				stabilizer = self.stabilizers[type][measure_position]
				for count in stabilizer.order:
					data_position = stabilizer.order[count]
					
					if data_position not in self.memberships:
						self.memberships[data_position] = {}
					
					if type not in self.memberships[data_position]:
						self.memberships[data_position][type] = []
					
					self.memberships[data_position][type].append(measure_position)
					if data_position not in self.primal.nodes():
						self.primal.add_node(data_position)	
						self.data[data_position] = Qubit(data_position, Charge(), type = 'data')
						
				self.generatePrimalEdges(stabilizer)


	def distance(self, qubit1, qubit2, lattice_type = 'dual'):
		if lattice_type == 'dual':
			lattice = self.dual.copy()
		else:
			lattice = self.shrunk[lattice_type].copy()
			
		if qubit1 in lattice.nodes() and  qubit2 in lattice.nodes():
			return nx.shortest_path_length(self.dual, qubit1, qubit2)
		elif qubit1 in lattice.nodes() and qubit2 not in lattice.nodes():
			type = self.syndromes[qubit2].type
			qubit2 = self.external[type][0]
			return nx.shortest_path_length(lattice, qubit1, qubit2) + 1
		elif qubit1 not in lattice.nodes() and qubit2 in lattice.nodes():
			type = self.syndromes[qubit1].type
			qubit1 = self.external[type][0]
			return nx.shortest_path_length(lattice, qubit1, qubit2) + 1
 		else:
 			type1 = self.syndromes[qubit1].type
			qubit1 = self.external[type1][0]
 			type2 = self.syndromes[qubit2].type
			qubit2 = self.external[type2][0]
 			return nx.shortest_path_length(lattice, qubit1, qubit2) + 2

	# def plot(self, plot_number, title, charge_types = ['X','Z'], lattice = 'primal'):
	# 	if lattice == 'primal':
	# 		plot_primal(self, plot_number, title, charge_types)
	# 	elif lattice == 'dual':
	# 		plot_dual(self, plot_number, title, charge_types)
	# 	else:
	# 		# input is shrunk lattice type
	# 		plot_shrunk(self, lattice, plot_number, title, charge_types)
