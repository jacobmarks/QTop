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
from mpl_toolkits.mplot3d import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from common import *
from geometry import *

# Planar and toric versions of 
# Alexander Kitaev's quantum double model




class KitaevCode(Code):
	'Common base class for all surface codes'

	def __init__(self, depth, dimension = 2):
		Code.__init__(self, depth, dimension)

	def generateTypes(self):
		self.types = {'X':{'num_sides':4},'Z':{'num_sides':4} }
		self.code = 'kitaev'

	def generateColors(self):
		self.colors = {'X':'red', 'Z':'blue','data':'black'}

	def CodeCycle(self, model):

		# find length of code cycle:
		max_num_sides = 0
		for type in self.types:
			if self.types[type]['num_sides'] > max_num_sides:
				max_num_sides = self.types[type]['num_sides']


		# Step 1:
		model.Identity(self)
		model.Initialize(self, 'X')

		# Step 2:
		model.Initialize(self, 'Z')
		model.Fourier(self, 'X')

		# # Steps 3-6:
		for order in range(max_num_sides):
			for type in self.types:
				model.Sum(self, order, type, 'Z')

		# Step 7:
		model.Fourier(self, 'X')
		model.Measure(self, 'Z')

		# Step 8:
		model.Measure(self, 'X')

		return self

	def generateDual(self):
		self.dual = nx.Graph()
		self.shrunk = {}
		for type in self.types:
			self.shrunk[type] = nx.Graph()

		for type1 in self.stabilizers:
			for measure1 in self.stabilizers[type1]:
				stabilizer1 = self.stabilizers[type1][measure1]
				center1 = stabilizer1.center
				for count in stabilizer1.order:
					data = stabilizer1.order[count]
					for type2 in self.memberships[data.position]:
						if type2 == type1:
							for measure2 in self.memberships[data.position][type2]:
								if measure2 != measure1:
									stabilizer2 = self.stabilizers[type2][measure2]
									center2 = stabilizer2.center
									edge_type = type1
									self.shrunk[edge_type].add_edge(*(center1, center2), type = edge_type)
									self.dual.add_edge(*(center1, center2), type = edge_type)


# Kitaev Surface Code
class KSC(PlanarCode, KitaevCode):
	'Planar Subclass of Kitaev Code'

	def __init__(self, depth, dimension = 2):
		PlanarCode.__init__(self, dimension)
		KitaevCode.__init__(self, depth, dimension)


	def generateStabilizers(self):
		self.boundary_data = {}
		for type in self.types:
			self.external[type] = {}
			self.boundary_measures[type] = {}
			self.boundary_data[type] = {'upper':[],'lower':[]}


		depth = self.depth
		for i in range(depth + 1):
			for j in range(depth):
				plaquettes = []
				X_position = (round(2 * i, 3), round(2 * j, 3))
				Z_position = (round(2 * j + 1, 3), round(2 * i - 1, 3))
				
				plaquettes.append({'type':'X','position':X_position})
				plaquettes.append({'type':'Z','position':Z_position})
				for measure in plaquettes:
					type, position = measure['type'], measure['position']
					center = Qubit(position, Charge(), type)
					num_sides = self.types[type]['num_sides']
					data = self.generateStabilizerData(center, 1, num_sides, 0)
				
					if i == 0 or i == depth:
						if i == 0:
							boundary_type, sign = 'lower', 1
						elif i == depth:
							boundary_type, sign = 'upper', -1

						if type == 'X':
							measure = (round(position[0] + sign * 2, 3), round(position[1], 3))
							bound_data = (round(position[0] + sign * 1, 3), round(position[1], 3))
						else:
							measure = (round(position[0], 3), round(position[1] + sign * 2, 3))
							bound_data = (round(position[0], 3), round(position[1] + sign * 1, 3))


						self.external[type][position] = {'measure':measure, 'data':bound_data}
						self.boundary_data[type][boundary_type].append(bound_data)
						self.boundary_measures[type][measure] = {'external':position,'data':bound_data}	
						continue

					if j == 0 or j == depth - 1:
						order = {}
						for count in range(num_sides):
							qubit = data[count]
							if type == 'X':
								if (j == 0 and count == 3) or (j == depth - 1 and count == 1):
									continue
							elif type == 'Z':
								if (j == 0 and count == 2) or (j == depth - 1 and count == 0):
									continue

							order[count] = qubit

					else:
						order = Order(data)

					self.stabilizers[type][center.position] = Stabilizer(center,data, order)
				


# Kitaev Toric Code
class KTC(ToricCode, KitaevCode):
	'Toric Subclass of Kitaev Code'

	def __init__(self, depth, dimension = 2):
		ToricCode.__init__(self, depth, dimension)
		KitaevCode.__init__(self, depth, dimension)

		
	def generateStabilizers(self):
		depth = self.depth
		length = 2 * depth
		self.length1, self.length2 = length, False

		for i in range(depth):
			for j in range(depth):
				X_position = (2 * i, 2 * j)
				Z_position = (2 * j + 1, 2 * i - 1)
				pos_types = {'X':X_position, 'Z':Z_position}


				for type in pos_types:
					position = pos_types[type]
					center = Qubit(position, Charge(), type)
					stabilizer_data = self.generateStabilizerData(center, 1, 4, 0)
					order = Order(stabilizer_data)
					toric_center = PlanarToToric(center, length)
					self.stabilizers[type][toric_center.position] = Stabilizer(toric_center,stabilizer_data, order)
					# Save planar coordinates 
					self.stabilizers[type][toric_center.position].planar_coords = position
				
