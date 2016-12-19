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

from math import *
import networkx as nx
# import matplotlib.pyplot as plt
from common import *
from error_models import *


class ColorCode(Code):
	'Common base class for all types of color codes'

	def __init__(self, depth, dimension = 2):
		self.code = 'color'
		Code.__init__(self, depth, dimension)

	def generateColors(self):
		self.colors = {'red':'red', 'blue':'blue', 'green':'green', 'data':'black'}
		
	def complementaryType(self, types):
		for type in self.types:
			if type not in types:
				complement = type
		return complement

	def complementaryTypes(self, type):
		for t1 in self.types:
			if t1 != type:
				t2 = self.complementaryType([type, t1])
				return t1, t2

	def Sign(self, count, num_sides = 6):
		if count % 2 == 0:
			return 1
		else:
			return -1

	def Plaquette(self, measure_qubit, type):
		checks = {}

		sides = self.types[type]['sides']
		angle = self.types[type]['angle']
		scale = self.types[type]['scale']

		for k in range(sides):
			x = scale * float(cos(2*pi*k/sides))
			y = scale * float(sin(2*pi*k/sides))
			if angle != 0:
				x_prime = float(cos(angle))*x - float(sin(angle))*y
				y_prime = float(sin(angle))*x + float(cos(angle))*y
				x, y = x_prime, y_prime
			data_qubit = (round(x + measure_qubit[0], 3) ,round(y  + measure_qubit[1], 3))
			checks[data_qubit] = k
		return checks

	def generateDual(self):

		for type1 in self.Stabilizers:
			for m1 in self.Stabilizers[type1]:
				for count in self.Stabilizers[type1][m1]['order']:
					data = self.Stabilizers[type1][m1]['order'][count]
					for type2 in self.Primal.node[data]['measures']:
						if type2 != type1:
							for m2 in self.Primal.node[data]['measures'][type2]:
								edge_type = self.complementaryType([type1, type2])
								self.Dual[edge_type].add_edge(*(m1, m2), type = edge_type)


	######## Code Cycle #######

	def CodeCycle(self, model, p):

	# length of code cycle is 2*(max number of sides + 2)
	# 'X' and 'Z' checks each require one step for initialization,
	# one step for each data qubit check operation, and a final measurement
		max_sides = max([self.types[type]['sides'] for type in self.types])


		# Initialization
		for type in self.types:
			self = model.Initialize(self, type, p)
		
		self = model.Identity(self, p)

		for charge_type in ['X','Z']:
			# Stabilizer Check Operations
			for count in range(max_sides):
				for type in self.types:
					sides = self.types[type]['sides']
					self = model.Sum(self, count, sides, type, charge_type, p)

			# Measurement
			for type in self.types:
				self = model.Measure(self, type, p)

		return self

	def generatePrimalEdges(self):
		for type in self.types:
			num_sides = self.types[type]['sides']
			for m in self.Stabilizers[type]:
				stabilizer = self.Stabilizers[type][m]
				for k in range(num_sides):
					if k in stabilizer['order']:
						vertex1 = stabilizer['order'][k]
						while (k+1)%num_sides not in stabilizer['order']:
							k += 1
						vertex2 = stabilizer['order'][(k+1)%num_sides]
						self.Primal.add_edge(*(vertex1,vertex2), color = 'black')

	def PrimalBound(self, count, type, measures):
		# need to make list of measures not contain externals
		dict1 = self.Stabilizers[type][measures[type]]['order']
		dict2 = self.Stabilizers[type][measures[type]]['data']
		if count in dict1:
			node = dict1[count]
			self.Stabilizers[type][measures[type]]['data'] = removekey(dict2, node)
			self.Stabilizers[type][measures[type]]['order'] = removekey(dict1, count)
			self.Primal.remove_node(node)
		return self 

	def hasLogicalError(self):
		# check if correction commutes or anticommutes with logical operator on boundary
		d = self.dimension
		# for ct in ['X','Z']:
		for ct in ['Z']:
			Sum = 0
			for ext in self.External['red']:
				for data in self.Stabilizers['red'][ext]['data']:
					count = self.Stabilizers['red'][ext]['data'][data]
					sign = self.Sign(count)
					c = self.Primal.node[data]['charge'][ct]
					Sum += sign * c
		
			if Sum % d != 0:
				return True
		return False


#############   6 - 6 - 6 Color Codes #############

class Color_6_6_6(ColorCode):

	def __init__(self, depth, dimension = 2):
		self.code = '6_6_6'
		ColorCode.__init__(self, depth, dimension)

	def generateCode(self):
		depth = self.depth

		self.types = {}
		for type in ['red', 'blue', 'green']:
			self.types[type] = {'sides':6,'scale':1, 'angle':0}
		

		for type in self.types:
			self.Dual[type] = nx.Graph()
			self.Stabilizers[type] = {}
			self.External[type] = []

		N = int(float(depth)/2)
		for i in range(N + 1):
			for j in range(N + 1 - i):
				mG = (round(3 * (i + float(j)/2), 3), round(float(3* sqrt(3) * j)/2, 3))
				mB = (round(mG[0] + float(3)/2, 3), round(mG[1] + float(sqrt(3))/2, 3))
				mR = (round(mG[0], 3), round(float(mG[1] + sqrt(3)), 3))
				measures = {'red':mR, 'blue':mB, 'green':mG}

				for m in measures:
					if j == 0:
						if m != 'red' or i == 0:
							continue
						self.External['red'].append(measures[m])
					if j == N - i:
						if m == 'blue':
							self.External['blue'].append(measures[m])

					if i == 0:
						if m == 'red':
							continue
						if m == 'green':
							if j == 0:
								continue
							self.External['green'].append(measures[m])

					self.Stabilizers[m][measures[m]] = {'data': {}, 'charge': Charge(), 'order':{}, 'sides':6}
					P = self.Plaquette(measures[m], m)
					for data in P:
						count = P[data]
						self.Stabilizers[m][measures[m]]['order'][count] = data
						self.Stabilizers[m][measures[m]]['data'][data] = count
						if data not in self.Primal.nodes():
							self.Primal.add_node(data, charge = Charge())
							self.Primal.node[data]['measures'] = {'red':[], 'blue':[], 'green':[]}
						self.Primal.node[data]['measures'][m].append(measures[m])

		self.generatePrimalEdges()

		for data in self.Primal.nodes():
			if len([type for type in self.types if self.Primal.node[data]['measures'][type] != []]) == 1:
				t = [type for type in self.types if self.Primal.node[data]['measures'][type] != []][0]
				m = self.Primal.node[data]['measures'][t][0]
				count = self.Stabilizers[t][m]['data'][data]
				del self.Stabilizers[t][m]['order'][count]
				del self.Stabilizers[t][m]['data'][data]
				self.Primal.remove_node(data)
				continue

		for data in self.Primal.nodes():
			for type in self.External:
				p0, p1 = self.External[type][0:2]
				if data not in self.Primal.nodes():
					break

				if colinear(p0, p1, data):
					for t in self.Primal.node[data]['measures']:
						if self.Primal.node[data]['measures'][t] == []:
							continue
						m = self.Primal.node[data]['measures'][t][0]
						count = self.Stabilizers[t][m]['data'][data]
						del self.Stabilizers[t][m]['order'][count]
						del self.Stabilizers[t][m]['data'][data]
					self.Primal.remove_node(data)
					break



