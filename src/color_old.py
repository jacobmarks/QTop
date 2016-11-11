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

# color codes proposed by Bombin and Delgado

# Including all 3 3-regular tilings of the plane
# 6-6-6, 4-8-8, and 4-6-12

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

	def generateCode(self):
		self.types = {'red', 'blue', 'green'}

		for type in self.types:
            self.Dual[type] = nx.Graph()
            self.Stabilizers[type] = {}
            self.Boundary[type] = {}
            self.External[type] = {}

	def generateDual(self):
		self.dual = nx.Graph()
		self.shrunk = {}
		for type in self.types:
			self.shrunk[type] = nx.Graph()

		for type1 in self.Stabilizers:
			for measure1_position in self.Stabilizers[type1]:
				Stabilizer1 = self.Stabilizers[type1][measure1_position]
				measure1_qubit = self.syndromes[measure1_position]

				for count in Stabilizer1.order:
					data = Stabilizer1.order[count]
					for type2 in self.memberships[data]:
						if type2 != type1:
							for measure2_position in self.memberships[data][type2]:
								Stabilizer2 = self.Stabilizers[type2][measure2_position]
								measure2_qubit = self.syndromes[measure2_position]
								edge_type = self.complementaryType([type1, type2])
								self.shrunk[edge_type].add_edge(*(measure1_position, measure2_position), type = edge_type)
								self.dual.add_edge(*(measure1_position, measure2_position), type = edge_type)


	######## Code Cycles #######

	# default code cycle is 1-qubit code cycle
	# but many other possible code cycles are allowed
	# including interleaved cycles


	# def OneQubitCodeCycle(self, model, p):

	# 	# length of code cycle is 2*(max number of sides + 2)
	# 	# 'X' and 'Z' checks each require one step for initialization,
	# 	# one step for each data qubit check operation, and a final measurement

	# 	max_num_sides = 0
	# 	for type in self.types:
	# 		if self.types[type]['num_sides'] > max_num_sides:
	# 			max_num_sides = self.types[type]['num_sides']

	# 	for charge_type in ['X','Z']:

	# 		# Initialization
	# 		for type in self.Stabilizers:
	# 			self = model.Initialize(self, type, p)
	# 		for type in self.Stabilizers:
	# 			self = model.Identity(self, p)

	# 		# Stabilizer Check Operations
	# 		for order in range(max_num_sides):
	# 			for type in self.types:
	# 				self = model.Sum(self, order, type, charge_type, p)

	# 		# Measurement
	# 		for type in self.types:
	# 			self = model.Measure(self, type, p)

	# 	return self


class Color_4_8_8(ColorCode):

	def __init__(self, depth, dimension = 2):
		self.code = '4_8_8'
		ColorCode.__init__(self, depth, dimension)

	def generateTypes(self):
		self.types = {}
		self.types['red'] = {'num_sides':4, 'scale':float(1)/sqrt(2)}
		self.types['blue'] = {'num_sides':8, 'scale':float(1)/(2*sin(float(pi)/8))}
		self.types['green'] = {'num_sides':8,'scale':float(1)/(2*sin(float(pi)/8))}
						

class Triangular_4_8_8(Color_4_8_8):

	def __init__(self, depth, dimension = 2):
		Color_4_8_8.__init__(self, depth, dimension)


	def generateCode(self):
		for type in self.types:
			self.External[type] = {}
			self.Boundary[type] = {}

		depth = self.depth
		N = int(float(depth -1)/2)
		x_dist = round(1 + sqrt(2)), 3)

		for j in range(N + 1):
			for i in range(N + 1 - j):
				plaquettes = []
				green_position = ( round(x_dist * ( 2 * i + j), 3), round(x_dist * j,3))
				plaquettes.append({'type':'green','angle':float(pi)/8,'position':green_position})
				plaquettes.append({'type':'blue','angle':float(pi)/8,'position':(round(green_position[0] + x_dist, 3), round(green_position[1], 3))})
				plaquettes.append({'type':'red','angle':0,'position':( round(green_position[0] - float(x_dist)/2, 3), round(green_position[1] + float(x_dist)/2, 3))})
				plaquettes.append({'type':'red','angle':0,'position':(round(green_position[0] + float(x_dist)/2, 3), round(green_position[1] + float(x_dist)/2, 3))})
				for measure_qubit in plaquettes:
					type, angle, position = measure_qubit['type'], measure_qubit['angle'], measure_qubit['position']
					self.syndromes[position] = Qubit(position, charge = Charge(), type = type)
					num_sides, scale = self.types[type]['num_sides'], self.types[type]['scale']
					
					if i == 0:
						if type == 'red':
							continue
						elif type == 'green':
							if j == 0:
								continue
							red_partners = [( round(position[0] + float(x_dist)/2, 3), round(position[1] + float(x_dist)/2, 3))]
							blue_partners = []
							blue_partners.append(( round(position[0] + x_dist, 3), round(position[1], 3)))
							blue_partners.append(( round(position[0], 3), round(position[1] + x_dist, 3)))
							self.External['green'][position] = {'red':red_partners, 'blue':blue_partners}
							continue

					if i == N - j:
						if type == 'red' and position == (round(green_position[0] + float(x_dist)/2, 3), round(green_position[1] + float(x_dist)/2, 3)):
							 continue
						if type == 'blue' and j != 0:
							red_partners = [( round(position[0] - float(x_dist)/2, 3), round(position[1] - float(x_dist)/2, 3))]
							green_partners = []
							green_partners.append(( round(position[0] - x_dist, 3), round(position[1], 3)))
							green_partners.append(( round(position[0], 3), round(position[1] - x_dist, 3)))
							self.External['blue'][position] = {'red':red_partners, 'green':green_partners}
							continue

					if j == 0:
						if type != 'red' or i == 0:
							continue
						else:
							if position == ( round(position[0] - float(x_dist)/2, 3), round(position[1] + float(x_dist)/2, 3)):
								green_partners = [( round(position[0] - float(x_dist)/2, 3), round(position[1] + float(x_dist)/2, 3))]
								blue_partners = [( round(position[0] + float(x_dist)/2, 3), round(position[1] + float(x_dist)/2, 3))]
							else:
								green_partners = [( round(position[0] + float(x_dist)/2, 3), round(position[1] + float(x_dist)/2, 3))]
								blue_partners = [( round(position[0] - float(x_dist)/2, 3), round(position[1] + float(x_dist)/2, 3))]
							self.External['red'][position] = {'green':green_partners, 'blue':blue_partners}
							continue
					

					measure = Qubit(position, Charge(), type)
					data = self.generateStabilizerData(position, scale, num_sides, angle)


					if (j == 1 and type in ['blue', 'green']) or (i == 0 and type == 'blue') or (type == 'green' and i == N - j):
						order = {}
						
						if (j == 1 and type in ['blue', 'green']):
							for count in range(num_sides):
								if count in [5,6]:
									continue
								if type == 'blue' and i == 0 and count in [2,3]:
									continue
								if type == 'green' and i == N - j and count in [0,1]:
									continue
								order[count] = data[count]

						if i == 0 and type == 'blue':
							self.boundary_syndromes['blue'][position] = []
							self.boundary_syndromes['blue'][position].append(( round(position[0] - float(x_dist)/2, 3), round(position[1], 3)))
							self.boundary_syndromes['blue'][position].append(( round(position[0], 3), round(position[1] + float(x_dist)/2, 3)))
							red_position = ( round(position[0] - float(x_dist)/4, 3), round(position[1] - float(x_dist)/4, 3))
							self.boundary_syndromes['red'][position] = red_position
							self.boundary_data['green'].append((round(red_position[0], 3), round(red_position[1] + float(1)/sqrt(2)), 3))
							if j != 0:
								self.boundary_data['green'].append((round(red_position[0] - float(1)/sqrt(2), 3), round(red_position[1]), 3))
							for count in range(num_sides):
								if count in [2,3]:
									continue
								if j == 1 and count in [5,6]:
									continue
								order[count] = data[count]

						if type == 'green' and i == N - j:
							self.boundary_syndromes['green'][position] = []
							self.boundary_syndromes['green'][position].append(( round(position[0] + float(x_dist)/2, 3), round(position[1], 3)))
							self.boundary_syndromes['green'][position].append(( round(position[0], 3), round(position[1] + float(x_dist)/2, 3)))
							red_position = ( round(position[0] + float(x_dist)/4, 3), round(position[1] - float(x_dist)/4, 3))
							self.boundary_syndromes['red'][position] = red_position
							self.boundary_data['blue'].append((round(red_position[0], 3), round(red_position[1] + float(1)/sqrt(2), 3)))
							if j != 0:
								self.boundary_data['blue'].append((round(red_position[0] + float(1)/sqrt(2), 3), round(red_position[1], 3)))
							for count in range(num_sides):
								if count in[0,1]:
									continue
								if j == 1 and count in [5,6]:
									continue
								order[count] = data[count]

							

				
					else:
						order = Order(data)
					self.Stabilizers[type][position] = Stabilizer(type, data, order)
						
						

class Color_6_6_6(ColorCode):

	def __init__(self, depth, dimension = 2):
		self.code = '6_6_6'
		ColorCode.__init__(self, depth, dimension)

	def generateTypes(self):
		self.types = {'red':{'num_sides':6,'scale':1}, 'blue':{'num_sides':6,'scale':1}, 'green':{'num_sides':6,'scale':1}}


class Triangular_6_6_6(Color_6_6_6):

	def __init__(self, depth, dimension = 2):
		Color_6_6_6.__init__(self, depth, dimension)

	def generateStabilizers(self):
		for type in self.types:
			self.External[type] = {}
			self.boundary_syndromes[type] = {}
			self.boundary_data[type] = []

		depth = self.depth
		N = int(float(depth)/2)
		for i in range(N + 1):
			for j in range(N + 1 - i):
				plaquettes = []
				green_position = (round(3 * (i + float(j)/2), 3), round(float(3* sqrt(3) * j)/2, 3))
				plaquettes.append({'type':'green','angle':0,'position':green_position})
				plaquettes.append({'type':'blue','angle':0,'position':(round(green_position[0] + float(3)/2, 3), round(green_position[1] + float(sqrt(3))/2, 3))})
				plaquettes.append({'type':'red','angle':0,'position':(round(green_position[0], 3), round(float(green_position[1] + sqrt(3)), 3))})
				for face in plaquettes:
					type, angle, position = face['type'], face['angle'], face['position']
					self.syndromes[position] = Qubit(position, charge = Charge(), type = type)
					num_sides, scale = self.types[type]['num_sides'], self.types[type]['scale']
					
					if j == 0:
						if type != 'red' or i == 0:
							continue
						green_partners = []
						green_partners.append(( round(position[0] + float(3)/2, 3), round(position[1] + float(sqrt(3))/2, 3)))
						green_partners.append(( round(position[0] - float(3)/2, 3), round(position[1] + float(sqrt(3))/2, 3)))
						blue_partners = [( round(position[0], 3), round(position[1] + sqrt(3), 3))]
						self.External['red'][position] = {'green':green_partners, 'blue':blue_partners}
						continue
					if j == N - i:
						if type == 'blue':
							green_partners = [( round(position[0] - float(3)/2, 3), round(position[1] - float(sqrt(3))/2, 3))]
							red_partners = []
							red_partners.append(( round(position[0] + float(3)/2, 3), round(position[1] + float(sqrt(3))/2, 3)))
							red_partners.append(( round(position[0] - float(3)/2, 3), round(position[1] + float(sqrt(3))/2, 3)))
							self.External['blue'][position] = {'green':green_partners, 'red':red_partners}
							continue

					if i == 0:
						if type == 'red':
							continue
						if type == 'green':
							red_partners = [( round(position[0] + float(3)/2, 3), round(position[1] - float(sqrt(3))/2, 3))]
							blue_partners = []
							blue_partners.append(( round(position[0] + float(3)/2, 3), round(position[1] + float(sqrt(3))/2, 3)))
							blue_partners.append(( round(position[0], 3), round(position[1] + sqrt(3), 3)))
							self.External['green'][position] = {'blue':blue_partners, 'red':red_partners}
							continue

					measure = Qubit(position, Charge(), type)
					data = self.generateStabilizerData(position, scale, num_sides, angle)
					

					if (j == 1 and type == 'green') or (i == 0 and type == 'blue') or (i == N - j and type == 'red'):
						order = {}
						for count in range(num_sides):
							if (j == 1 and type == 'green'):
								if count in [4,5]:
									continue
							if i == 0 and type == 'blue':
								if count in [2,3]:
									continue
							if type == 'red' and i == N - j:
								if count in[0,1]:
									continue

							order[count] = data[count]
					else:
						order = Order(data)
					
					self.Stabilizers[type][position] = Stabilizer(type, data, order)
		

