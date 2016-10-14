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
import sys

# Planar and toric versions of 
# Alexander Kitaev's quantum double model




class KitaevCode(Code):
	'Common base class for all surface codes'

	def __init__(self, depth, dimension = 2):
		self.code = 'kitaev'
		Code.__init__(self, depth, dimension)


	def generateTypes(self):
		self.types = {'X':{'num_sides':4},'Z':{'num_sides':4} }

	def generateColors(self):
		self.colors = {'X':'red', 'Z':'blue','data':'black'}

	def CodeCycle(self, model, p = 0):

		# find length of code cycle:
		num_sides = 4


		# Step 1:
		self = model.Identity(self, p)
		self = model.Initialize(self, 'X', p)

		# # Step 2:
		self = model.Initialize(self, 'Z', p)
		self = model.Fourier(self, 'X', p)

		# # # Steps 3-6:
		for count in range(num_sides):
			for type in self.types:
				charge_type = type
				self = model.Sum(self, count, type, charge_type, p)

		# # Step 7:
		self = model.Fourier(self, 'X', p)
		self = model.Measure(self, 'Z', p)

		# # Step 8:
		self = model.Measure(self, 'X', p)

		return self

	def generateDual(self):
		self.dual = nx.Graph()
		self.shrunk = {}
		for type in self.types:
			self.shrunk[type] = nx.Graph()

		for type1 in self.stabilizers:
			for measure1_position in self.stabilizers[type1]:
				stabilizer1 = self.stabilizers[type1][measure1_position]
				measure1_qubit = self.syndrome[measure1_position]
				self.syndrome[measure1_position] = Qubit(measure1_position, Charge(), type = type1)
				for count in stabilizer1.order:
					data_position = stabilizer1.order[count]
					for type2 in self.memberships[data_position]:
						if type2 == type1:
							for measure2_position in self.memberships[data_position][type2]:
								if measure2_position != measure1_position:
									edge_type = type1
									self.shrunk[edge_type].add_edge(*(measure1_position, measure2_position), type = edge_type)
									self.dual.add_edge(*(measure1_position, measure2_position), type = edge_type)



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
					planar_measure_position = pos_types[type]
					measure_qubit = Qubit(planar_measure_position, Charge(), type)
					stabilizer_data = self.generateStabilizerData(planar_measure_position, 1, 4, 0)
					order = Order(stabilizer_data)
					toric_measure_position = ToricCoordinates(planar_measure_position, length)
					self.stabilizers[type][toric_measure_position] = Stabilizer(type,stabilizer_data, order)
					# Save planar coordinates 
					self.stabilizers[type][toric_measure_position].planar_coords = planar_measure_position
					self.syndrome[toric_measure_position] = Qubit(toric_measure_position, charge = Charge(), type = type)
