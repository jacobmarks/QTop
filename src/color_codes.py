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
import matplotlib.pyplot as plt
from common import *
from error_models import *
from visualization import *


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

		# If stabilizer has "external" then add to dual, remove from stabilizers?
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

		for charge_type in ['X','Z']:

			# Initialization
			for type in self.types:
				self = model.Initialize(self, type, p)
			for type in self.types:
				self = model.Identity(self, p)

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
		d = self.dimension
		for charge_type in ['X','Z']:
			E_FLAG = True
			for type in self.types:
				Sum = 0
				for node in self.Boundary[type]:
					Sum += self.Primal.node[node]['charge'][charge_type]
				if Sum % d == 0:
					E_FLAG = False
			if E_FLAG:
				return True
		return False


#############   6 - 6 - 6 Color Codes #############

class Color_6_6_6(ColorCode):

	def __init__(self, depth, dimension = 2):
		self.code = '6_6_6'
		ColorCode.__init__(self, depth, dimension)

	def generateCode(self):
		self.types = {}
		self.types['red'] = {'sides':6,'scale':1, 'angle':0}
		self.types['blue'] = {'sides':6,'scale':1, 'angle':0}
		self.types['green'] = {'sides':6,'scale':1, 'angle':0}
		depth = self.depth

		for type in self.types:
			self.Dual[type] = nx.Graph()
			self.Stabilizers[type] = {}
			self.Boundary[type] = []
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
						continue
					if j == N - i:
						if m == 'blue':
							self.External['blue'].append(measures[m])
							continue

					if i == 0:
						if m == 'red':
							continue
						if m == 'green':
							if j == 0:
								continue
							self.External['green'].append(measures[m])
							continue

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

					if j == 1 and m == 'green':
						for count in [4,5]:
							self = self.PrimalBound(count, m, measures)
					if i == 0 and m == 'blue':
						for count in [2,3]:
							self = self.PrimalBound(count, m, measures)
					if m == 'red' and i == N - j:
						for count in [0,1]:
							self = self.PrimalBound(count, m, measures)

		self.generatePrimalEdges()



#############   4 - 8 - 8 Color Codes #############



class Color_4_8_8(ColorCode):

	def __init__(self, depth, dimension = 2):
		self.code = '4_8_8'
		ColorCode.__init__(self, depth, dimension)

	def generateCode(self):
		self.types = {}
		self.types['red'] = {'sides':4, 'scale':float(1)/sqrt(2), 'angle':0}
		self.types['blue'] = {'sides':8, 'scale':float(1)/(2*sin(float(pi)/8)), 'angle':float(pi)/8}
		self.types['green'] = {'sides':8,'scale':float(1)/(2*sin(float(pi)/8)), 'angle':float(pi)/8}
		depth = self.depth

		for type in self.types:
			self.Dual[type] = nx.Graph()
			self.Stabilizers[type] = {}
			self.Boundary[type] = {}
			self.External[type] = []

		N = int(float(depth -1)/2)
		x = 1 + sqrt(2)

		for j in range(N + 1):
			for i in range(N + 1 - j):
				mG = (round(x * ( 2 * i + j), 3), round(x * j,3))
				mB = (round(mG[0] + x, 3), mG[1])
				mR1 = (round(mG[0] - float(x)/2, 3), round(float(mG[1] + float(x)/2), 3))
				mR2 = (round(mG[0] + float(x)/2, 3), round(float(mG[1] + float(x)/2), 3))
				measures = {'red1':mR1, 'red2':mR2, 'blue':mB, 'green':mG}

				for m in measures:
					if i == 0:
						if m == 'red1' or m == 'red2':
							continue
						elif type == 'green':
							if j == 0:
								continue
							self.External['red'].append(mG)
							self.External['blue'].append(mG)
							continue

					if i == N - j:
						if m == 'red2':
							 continue
						if m == 'blue' and j != 0:
							self.External['green'].append(mB)
							self.External['red'].append(mB)
							continue

					if j == 0:
						if (m == 'blue') or (m == 'green') or i == 0:
							continue
						self.External['blue'].append(measures[m])
						self.External['green'].append(measures[m])
						continue

					if m == 'red1' or m == 'red2':
						type = 'red'
					else:
						type = m

					sides = self.types[type]['sides']
					self.Stabilizers[type][measures[m]] = {'data': {}, 'charge': Charge(), 'order':{}, 'sides':sides}
					P = self.Plaquette(measures[m], type)
					for data in P:
						count = P[data]
						self.Stabilizers[type][measures[m]]['order'][count] = data
						self.Stabilizers[type][measures[m]]['data'][data] = count
						if data not in self.Primal.nodes():
							self.Primal.add_node(data, charge = Charge())
							self.Primal.node[data]['measures'] = {'red':[], 'blue':[], 'green':[]}
						self.Primal.node[data]['measures'][type].append(measures[m])


					if (j == 1 and type in ['blue', 'green']):
						for count in [5,6]:
							self = self.PrimalBound(count, m, measures)
					if type == 'blue' and i == 0:
						if j == 1:
							for count in [2,3,5,6]:
								self = self.PrimalBound(count, m, measures)
						else:
							for count in [2,3]:
								self = self.PrimalBound(count, m, measures)
					if type == 'green' and i == N - j:
						for count in [0,1]:
							self = self.PrimalBound(count, m, measures)

		self.generatePrimalEdges()






