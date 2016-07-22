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
	def __init__(self, position, charge, type = 'data'):
		self.position = position
		self.charge = Charge()
		self.type = type

class Stabilizer:
	def __init__(self, center, data, order):
		self.center = center
		self.data = data
		self.order = order

def Order(data_list):
	count, order = 0, {}
	for qubit in data_list:
		order[count] = qubit
		count += 1
	return order

class Code:
	'Common base class for all types of topological codes'

	def __init__(self, depth, dimension = 2):
		self.depth = depth
		self.dimension = dimension
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

	def generateStabilizerData(self, center, scale, num_sides, angle = 0):
		data = []
		for k in range(num_sides):
			x = scale * float(cos(2*pi*k/num_sides))
			y = scale * float(sin(2*pi*k/num_sides))
			if angle != 0:
				x_prime = float(cos(angle))*x - float(sin(angle))*y
				y_prime = float(sin(angle))*x + float(cos(angle))*y
				x, y = x_prime, y_prime
			position = (round(x + center.position[0], 3) ,round(y  + center.position[1], 3))
			data.append(Qubit(position, charge = Charge()))
		return data

	def generatePrimalEdges(self, stabilizer):
		stabilizer_type = stabilizer.center.type
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
			for center in self.stabilizers[type]:
				stabilizer = self.stabilizers[type][center]
				for order in stabilizer.order:
					data = stabilizer.order[order]
					
					if data.position not in self.memberships:
						self.memberships[data.position] = {}
					
					if type not in self.memberships[data.position]:
						self.memberships[data.position][type] = []
					
					self.memberships[data.position][type].append(center)
					if data not in self.primal.nodes():
						self.primal.add_node(data)	
						
				self.generatePrimalEdges(stabilizer)


	def distance(self, qubit1, qubit2):
		if qubit1 in self.dual.nodes() and  qubit2 in self.dual.nodes():
			return nx.shortest_path_length(self.dual, qubit1, qubit2)
		elif qubit1 in self.dual.nodes() and qubit2 not in self.dual.nodes():
			qubit2 = self.external[qubit2.type][qubit2]
			return nx.shortest_path_length(self.dual, qubit1, qubit2) + 1
		elif qubit1 not in self.dual.nodes() and qubit2 in self.dual.nodes():
			qubit1 = self.external[qubit1.type][qubit1]
			return nx.shortest_path_length(self.dual, qubit1, qubit2) + 1
 		else:
 			qubit1 = self.external[qubit1.type][qubit1]
 			qubit2 = self.external[qubit2.type][qubit2]
 			return nx.shortest_path_length(self.dual, qubit1, qubit2) + 2

	def plot(self, charge_type, plot_number, title, lattice = 'primal'):
		if lattice == 'primal':
			plot_primal(self, charge_type, plot_number, title)
		elif lattice == 'dual':
			plot_dual(self, charge_type, plot_number, title)
		else:
			# input is shrunk lattice type
			plot_shrunk(self, lattice, charge_type, plot_number, title)
