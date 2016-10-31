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
from error_models import *

# Planar and toric versions of 
# color codes proposed by Bombin and Delgado

# Including all 3 3-regular tilings of the plane
# 6-6-6, 4-8-8, and 4-6-12

class ColorCode(Code):
	'Common base class for all types of color codes'

	def __init__(self, depth, dimension = 2):
		self.code = 'color'
		Code.__init__(self, depth, dimension)


	def generateTypes(self):
		self.types = {'red', 'blue', 'green'}

	def generateColors(self):
		self.colors = {'red':'red', 'blue':'blue', 'green':'green', 'data':'black'}
		
	def complementaryType(self, types):
		for type in self.types:
			if type not in types:
				complement = type
		return complement

	def generateDual(self):
		self.dual = nx.Graph()
		self.shrunk = {}
		for type in self.types:
			self.shrunk[type] = nx.Graph()

		for type1 in self.stabilizers:
			for measure1_position in self.stabilizers[type1]:
				stabilizer1 = self.stabilizers[type1][measure1_position]
				measure1_qubit = self.syndromes[measure1_position]

				for count in stabilizer1.order:
					data = stabilizer1.order[count]
					for type2 in self.memberships[data]:
						if type2 != type1:
							for measure2_position in self.memberships[data][type2]:
								stabilizer2 = self.stabilizers[type2][measure2_position]
								measure2_qubit = self.syndromes[measure2_position]
								edge_type = self.complementaryType([type1, type2])
								self.shrunk[edge_type].add_edge(*(measure1_position, measure2_position), type = edge_type)
								self.dual.add_edge(*(measure1_position, measure2_position), type = edge_type)


	######## Code Cycles #######

	# default code cycle is 1-qubit code cycle
	# but many other possible code cycles are allowed
	# including interleaved cycles


	def OneQubitCodeCycle(self, model, p):

		# length of code cycle is 2*(max number of sides + 2)
		# 'X' and 'Z' checks each require one step for initialization,
		# one step for each data qubit check operation, and a final measurement

		max_num_sides = 0
		for type in self.types:
			if self.types[type]['num_sides'] > max_num_sides:
				max_num_sides = self.types[type]['num_sides']

		for charge_type in ['X','Z']:

			# Initialization
			for type in self.stabilizers:
				self = model.Initialize(self, type, p)
			for type in self.stabilizers:
				self = model.Identity(self, p)

			# Stabilizer Check Operations
			for order in range(max_num_sides):
				for type in self.types:
					self = model.Sum(self, order, type, charge_type, p)

			# Measurement
			for type in self.types:
				self = model.Measure(self, type, p)

		return self


class Color_4_8_8(ColorCode):

	def __init__(self, depth, dimension = 2):
		self.code = '4_8_8'
		ColorCode.__init__(self, depth, dimension)

	def generateTypes(self):
		self.types = {}
		self.types['red'] = {'num_sides':4, 'scale':float(1)/sqrt(2)}
		self.types['blue'] = {'num_sides':8, 'scale':float(1)/(2*sin(float(pi)/8))}
		self.types['green'] = {'num_sides':8,'scale':float(1)/(2*sin(float(pi)/8))}
						

class Toric_4_8_8(ToricCode, Color_4_8_8):

	def __init__(self, depth, dimension = 2):
		ToricCode.__init__(self, depth, dimension)
		Color_4_8_8.__init__(self, depth, dimension)
		

	def generateStabilizers(self):

		depth = self.depth
		N = int(float(depth)/4)
		# length = (2 + sqrt(2)) * N
		x_dist = round(2 * (float(1)/2 + float(1)/sqrt(2)), 3)
		length = N * x_dist
		
		self.length1, self.length2 = length, False

		for i in range(N):
			for j in range(N):
				plaquettes = []
				green_position = ( round(x_dist * ( 2 * i + j), 3), round(x_dist * j,3))
				plaquettes.append({'type':'green','angle':float(pi)/8,'position':green_position})
				plaquettes.append({'type':'blue','angle':float(pi)/8,'position':(round(green_position[0] + x_dist, 3), round(green_position[1], 3))})
				plaquettes.append({'type':'red','angle':0,'position':( round(green_position[0] - float(x_dist)/2, 3), round(green_position[1] + float(x_dist)/2, 3))})
				plaquettes.append({'type':'red','angle':0,'position':(round(green_position[0] + float(x_dist)/2, 3), round(green_position[1] + float(x_dist)/2, 3))})
				for face in plaquettes:
					type, angle, position = face['type'], face['angle'], face['position']
					num_sides, scale = self.types[type]['num_sides'], self.types[type]['scale']
					measure = Qubit(position, Charge(), self.colors[type])
					stabilizer_data = self.generateStabilizerData(position, scale, num_sides, angle)
					order = Order(stabilizer_data)
					toric_position = ToricCoordinates(position, length)
					self.stabilizers[type][toric_position] = Stabilizer(type,stabilizer_data, order)
					self.stabilizers[type][toric_position].planar_coords = position
					self.syndromes[toric_position] = Qubit(toric_position, charge = Charge(), type = type)


class Color_6_6_6(ColorCode):

	def __init__(self, depth, dimension = 2):
		self.code = '6_6_6'
		ColorCode.__init__(self, depth, dimension)

	def generateTypes(self):
		self.types = {'red':{'num_sides':6,'scale':1}, 'blue':{'num_sides':6,'scale':1}, 'green':{'num_sides':6,'scale':1}}



class Toric_6_6_6(ToricCode, Color_6_6_6):

	def __init__(self, depth, dimension = 2):
		ToricCode.__init__(self, depth, dimension)
		Color_6_6_6.__init__(self, depth, dimension)
		

	def generateStabilizers(self):

		depth = self.depth
		N = int(float(depth)/2)
		length1, length2 = float(3 * N),3 * sqrt(3) * float(N)/2
		self.length1, self.length2 = length1, length2
		for i in range(N):
			for j in range(N):
				plaquettes = []
				green_position = (round(3 * (i + float(j)/2), 3), round(float(3* sqrt(3) * j)/2, 3))
				plaquettes.append({'type':'green','angle':0,'position':green_position})
				plaquettes.append({'type':'blue','angle':0,'position':(round(green_position[0] + float(3)/2, 3), round(green_position[1] + float(sqrt(3))/2, 3))})
				plaquettes.append({'type':'red','angle':0,'position':(round(green_position[0], 3), round(green_position[1] + sqrt(3), 3))})
				for face in plaquettes:
					type, angle, position = face['type'], face['angle'], face['position']
					num_sides, scale = self.types[type]['num_sides'], self.types[type]['scale']
					measure = Qubit(position, Charge(), self.colors[type])
					stabilizer_data = self.generateStabilizerData(position, scale, num_sides, angle)
					order = Order(stabilizer_data)
					toric_position = ToricCoordinates(position, length1, length2)
					self.stabilizers[type][toric_position] = Stabilizer(type,stabilizer_data, order)
					self.stabilizers[type][toric_position].planar_coords = position
					self.syndromes[toric_position] = Qubit(toric_position, charge = Charge(), type = type)


