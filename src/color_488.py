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
from color_codes import *
from common import *
from error_models import *
from visualization import *




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