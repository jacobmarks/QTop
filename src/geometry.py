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


# planar and toric geometries for simulating 
# topological codes


############ Planar Code class and methods ############


class PlanarCode(Code):

	def __init__(self, dimension = 2):
		self.geometry = 'planar'
		self.boundary_syndrome = {}
		self.boundary_data = {}
		self.external = {}

	def hasLogicalError(self, charge_type):
		dimension = self.dimension
		netCharge = 0
		if self.code == 'color':
			for type in self.boundary_data:
				NetCharge = 0
				for bound_position in self.boundary_data[type]:
					NetCharge += self.data[bound_position].charge[charge_type]

				if netCharge% dimension != 0:
					return True
			return False

		elif self.code == 'kitaev':
			for type in self.boundary_data:
				for boundary_type in self.boundary_data[type]:
					NetCharge = 0
					for bound_position in self.boundary_data[type][boundary_type]:
						NetCharge += self.data[bound_position].charge[charge_type]

					if netCharge% dimension != 0:
						return True
			return False


	def plot_primal(self, plot_number, title, charge_types = ['X','Z']):
		Primal = plt.figure(plot_number)
		d = self.dimension

		for charge_type in charge_types:
			for type in self.types:
				for measure_position in self.stabilizers[type]:
					measure_qubit = self.syndrome[measure_position]
					charge = measure_qubit.charge[charge_type]
					if charge != 0:
						plt.scatter(*measure_position,marker="*", color = self.colors[type], s=200*float(charge)/(d-1))
					else:
						plt.scatter(*measure_position, color = self.colors[type])


			for type in self.types:
				for ext in self.external[type]:
					plt.scatter(*ext, marker="D", color = self.colors[type])

			for position in self.primal.nodes():
				charge = self.data[position].charge[charge_type]
				# if charge != 0:
				# 	plt.scatter(*position,marker="*", color= self.colors['data'], s=200*float(charge)/(d-1))
				# else:
				plt.scatter(*position, color = self.colors['data'])


		for edge in self.primal.edges():
			plt.plot((edge[0][0],edge[1][0]), (edge[0][1],edge[1][1]), color = self.colors['data'])

		plt.title(str(title))
		return plt.figure(plot_number)

	def plot_dual(self, plot_number, title, charge_types = ['X','Z']):
		Dual = plt.figure(plot_number)
		d = self.dimension

		for charge_type in ['X', 'Z']:

			for type in self.types:
				for measure_position in self.stabilizers[type]:
					measure_qubit = self.syndrome[measure_position]

					charge = measure_qubit.charge[charge_type]
					if charge != 0:
						plt.scatter(*measure_position,marker="*", color = self.colors[type], s=200*float(charge)/(d-1))

		for type in self.types:
			for ext in self.external[type]:
				plt.scatter(*ext, marker="D", color = self.colors[type])

		for edge in self.dual.edges():
			type = self.dual.get_edge_data(*edge)['type']
			[x0,y0] = edge[0][0], edge[0][1]
			[x1,y1] = edge[1][0], edge[1][1]
			plt.plot(*((x0, x1), (y0, y1)), color = self.colors[type])

		plt.title(str(title))
		return plt.figure(plot_number)

	def plot_shrunk(self, shrunk_type, plot_number, title, charge_types = ['X','Z']):
		Dual = plt.figure(plot_number)
		d = self.dimension

		for type in self.stabilizers:
			if type == shrunk_type and self.code == 'kitaev':
				for ext in self.external[type]:
					plt.scatter(*ext, marker="D", color = self.colors[type])
				continue
			for measure_position in self.stabilizers[type]:
				measure_qubit = self.syndrome[measure_position]
				for charge_type in charge_types:
					charge = measure_qubit.charge[charge_type]
					if charge != 0:
						plt.scatter(*measure_position,marker="*", color = self.colors[type], s=200*float(charge)/(d-1))

			if type != shrunk_type and self.code == 'color':
				for ext in self.external[type]:
					plt.scatter(*ext, marker="D", color = self.colors[type])
			

		for edge in self.shrunk[shrunk_type].edges():
			[x0,y0] = edge[0][0], edge[0][1]
			[x1,y1] = edge[1][0], edge[1][1]
			plt.plot(*((x0, x1), (y0, y1)), color = self.colors[shrunk_type])

		plt.title(str(title))
		return plt.figure(plot_number)


############ Toric Code class and methods ############

class ToricCode:

	def __init__(self, depth, dimension = 2):
		self.geometry = 'toric'
		self.length1, self.length2 = False, False

	def generateStabilizerData(self, measure_position, scale, num_sides, angle = 0):
		planar_positions = Code.generateStabilizerData(self, measure_position, scale, num_sides, angle)
		toric_positions = []
		for position in planar_positions:
			toric_positions.append(ToricCoordinates(position, self.length1, self.length2))
		return toric_positions

	

	def hasLogicalError(self, charge_type):
		# Make graph from subgraph of dual using only stabilizers with some, but not all
		# of their member data qubits with non-zero charge
		# look for self-loops with non-trivial crossing number
		# if one exists, we have a logical error

		# store list of positions of nonzero data
		charged_data = []
		for position in self.data:
			data = self.data[position]
			if data.charge[charge_type] != 0:
				charged_data.append(position)

		subgraph_nodes = []
		HasTrivialData, HasNontrivialData = False, False
		for type in self.stabilizers:
			for measure_position in self.stabilizers[type]:
				stabilizer = self.stabilizers[type][measure_position]
			
			for data in stabilizer.data:
				if data in charged_data:
					HasNontrivialData = True
				else:
					HasTrivialData = True

			if HasTrivialData and HasNontrivialData:
				subgraph_nodes.append(measure_position)



		Homologies = self.dual.subgraph(subgraph_nodes)
		# Find all cycles, then decide if each is trivial
		for cycle in nx.cycle_basis(Homologies):
			print "cycle: ", cycle
			if self.notTrivialCycle(cycle):
				return True

		return False

	def notTrivialCycle(self, cycle):
		# If odd number of crossing points for x or y, 
		# then we have a non-trival homology
		# i.e. a logical operator
		x_crosses, z_crosses = 0, 0
		dim = self.dimension
		length1 = self.length1
		if self.length2 != False:
			length2 = self.length2
		else:
			length2 = length1

		cycle_length = len(cycle)
		for link in range(cycle_length):
			node1 = cycle[link]
			node2 = cycle[(link+ 1)%cycle_length]
			type1, type2 = node1.type, node2.type

			pos1 = self.stabilizers[type1][node1.position].planar_coords
			pos2 = self.stabilizers[type1][node2.position].planar_coords

			x_dist, z_dist = abs(pos2[0] - pos1[0]), abs(pos2[1] - pos1[1])
			if (pos1[0] < x_dist and (length1 - pos2[0]) < x_dist):
				x_crosses += 1
			elif (pos2[0] < x_dist and (length1 - pos1[0]) < x_dist):
				x_crosses -= 1

			if (pos1[1] < z_dist and (length2 - pos2[1]) < z_dist):
				z_crosses += 1
			elif (pos2[1] < z_dist and (length2 - pos1[1]) < z_dist):
				z_crosses -= 1

		if x_crosses%dim != 0 or z_crosses%dim != 0:
			return True
		else:
			return False


	def plot_primal(self, plot_number, title, charge_types = ['X','Z']):
		Primal = plt.figure(plot_number)
		ax = Primal.add_subplot(111, projection='3d')
		d = self.dimension

		for type in self.stabilizers:
			for target_position in self.stabilizers[type]:
				stabilizer = self.stabilizers[type][target_position]
				target_qubit = self.syndrome[target_position]
				color = self.colors[target_qubit.type]
				[x,y,z] = target_position[0], target_position[1], target_position[2]
				for charge_type in charge_types:
					charge = target_qubit.charge[charge_type]
					if charge != 0:
						ax.scatter(x,y,z,marker="*",color=color,s=200*float(charge)/(d-1))


		for data_position in self.data:
			data_qubit = self.data[data_position]
			color = self.colors[data_qubit.type]
			[x,y,z] = data_position[0], data_position[1], data_position[2]
			for charge_type in charge_types:
				charge = data_qubit.charge[charge_type]
				if charge != 0:
					ax.scatter(x,y,z,marker="*",color=color,s=200*float(charge)/(d-1))

		for edge in self.primal.edges():
			[x0,y0,z0] = edge[0][0], edge[0][1], edge[0][2]
			[x1,y1,z1] = edge[1][0], edge[1][1], edge[1][2]
			ax.plot([x0,x1], [y0,y1], [z0,z1], color = 'black')

		plt.title(str(title))
		return plt.figure(plot_number)


	def plot_dual(self, plot_number, title, charge_types = ['X','Z']):

		Dual = plt.figure(plot_number)
		ax = Dual.add_subplot(111, projection='3d')
		d = self.dimension

		for edge in self.dual.edges():
			type = self.dual.get_edge_data(*edge)['type']
			[x0,y0,z0] = edge[0][0], edge[0][1], edge[0][2]
			[x1,y1,z1] = edge[1][0], edge[1][1], edge[1][2]
			ax.plot([x0,x1], [y0,y1], [z0,z1], color = self.colors[type])
		
		for type in self.stabilizers:
			for target_position in self.stabilizers[type]:
				target_qubit = self.syndrome[target_position]
				color = self.colors[target_qubit.type]
				[x,y,z] = target_position[0], target_position[1], target_position[2]
				for charge_type in ['X', 'Z']:
					charge = target_qubit.charge[charge_type]
					if charge != 0:
						ax.scatter(x,y,z,marker="*",color=color,s=200*float(charge)/(d-1))

		plt.title(str(title))
		return plt.figure(plot_number)

	def plot_shrunk(self, shrunk_type, plot_number, title, charge_types = ['X','Z']):
		# plot all edges of shrunk type

		Dual = plt.figure(plot_number)
		ax = Dual.add_subplot(111, projection='3d')
		d = self.dimension

		for edge in self.shrunk[shrunk_type].edges():
			[x0,y0,z0] = edge[0][0], edge[0][1], edge[0][2]
			[x1,y1,z1] = edge[1][0], edge[1][1], edge[1][2]
			ax.plot([x0,x1], [y0,y1], [z0,z1], color = self.colors[shrunk_type])
		
		for type in self.stabilizers:
			if type == shrunk_type:
				continue

			for target_position in self.stabilizers[type]:
				target_qubit = self.syndrome[target_position]
				color = self.colors[target_qubit.type]
				[x,y,z] = target_position[0], target_position[1], target_position[2]
				for charge_type in ['X', 'Z']:
					charge = target_qubit.charge[charge_type]
					if charge != 0:
						ax.scatter(x,y,z,marker="*",color=color,s=200*float(charge)/(d-1))

		plt.title(str(title))
		return plt.figure(plot_number)



def ToricCoordinates(planar_position,length1, length2 = False):
	if length2 == False:
		length2 = length1
	u = 2 * pi * float(planar_position[0])/length1
	v = 2 * pi * float(planar_position[1])/length2
	x = round(float(cos(u)) * (3 + float(cos(v))), 3)
	y = round(float(sin(u)) * (3 + float(cos(v))), 3)
	z = round(float(sin(v)), 3)
	return (x, y, z)
















