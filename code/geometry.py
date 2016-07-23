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
		self.boundary_measures = {}
		self.boundary_data = {}
		self.external = {}

	def hasLogicalError(self, charge_type):
		dimension = self.dimension
		netCharge = 0
		if self.code == 'color':
			for type in self.boundary_data:
				for qubit in self.primal.nodes():
					if qubit.position in self.boundary_data[type]:
						charge = qubit.charge[charge_type]
						netCharge += charge

				if netCharge% dimension != 0:
					return True
			return False

		elif self.code == 'kitaev':
			for type in self.boundary_data:
				for boundary_type in self.boundary_data[type]:
					for qubit in self.primal.nodes():
						if qubit.position in self.boundary_data[type][boundary_type]:
							charge = qubit.charge[charge_type]
							netCharge += charge

					if netCharge% dimension != 0:
						return True
			return False


	def plot_primal(self, charge_type, plot_number, title):
		Primal = plt.figure(plot_number)
		d = self.dimension

		for type in self.types:
			for target in self.stabilizers[type]:
				center = self.stabilizers[type][target].center
				position = center.position
				charge = center.charge[charge_type]
				if charge != 0:
					plt.scatter(*position,marker="*", color = self.colors[type], s=200*float(charge)/(d-1))
				else:
					plt.scatter(*position, color = self.colors[type])


		for type in self.types:
			for ext in self.external[type]:
				plt.scatter(*ext, marker="D", color = self.colors[type])

		for data in self.primal.nodes():
			type = data.type
			position = data.position
			charge = data.charge[charge_type]
			if charge != 0:
				plt.scatter(*position,marker="*", color= self.colors[type], s=200*float(charge)/(d-1))
			else:
				plt.scatter(*position, color = self.colors[type])


		for edge in self.primal.edges():
			plt.plot((edge[0].position[0],edge[1].position[0]), (edge[0].position[1],edge[1].position[1]), color = 'black')

		plt.title(str(title))
		return plt.figure(plot_number)

	def plot_dual(self, charge_type, plot_number, title):
		Dual = plt.figure(plot_number)
		d = self.dimension

		for type in self.types:
			for target in self.stabilizers[type]:
				center = self.stabilizers[type][target].center
				position = center.position

				charge = center.charge[charge_type]
				if charge != 0:
					plt.scatter(*position,marker="*", color = self.colors[type], s=200*float(charge)/(d-1))

		for type in self.types:
			for ext in self.external[type]:
				plt.scatter(*ext, marker="D", color = self.colors[type])

		for edge in self.dual.edges():
			type = self.dual.get_edge_data(*edge)['type']
			[x0,y0] = edge[0].position[0], edge[0].position[1]
			[x1,y1] = edge[1].position[0], edge[1].position[1]
			plt.plot(*((x0, x1), (y0, y1)), color = self.colors[type])

		plt.title(str(title))
		return plt.figure(plot_number)

	def plot_shrunk(self, shrunk_type, charge_type, plot_number, title):
		Dual = plt.figure(plot_number)
		d = self.dimension

		for type in self.stabilizers:
			if type == shrunk_type:
				continue
			for target in self.stabilizers[type]:
				center = self.stabilizers[type][target].center
				position = center.position
				charge = center.charge[charge_type]
				if charge != 0:
					plt.scatter(*position,marker="*", color = self.colors[type], s=200*float(charge)/(d-1))

			for ext in self.external[type]:
				plt.scatter(*ext, marker="D", color = self.colors[type])

		for edge in self.shrunk[shrunk_type].edges():
			[x0,y0] = edge[0].position[0], edge[0].position[1]
			[x1,y1] = edge[1].position[0], edge[1].position[1]
			plt.plot(*((x0, x1), (y0, y1)), color = self.colors[shrunk_type])

		plt.title(str(title))
		return plt.figure(plot_number)


############ Toric Code class and methods ############

class ToricCode:

	def __init__(self, depth, dimension = 2):
		self.geometry = 'toric'

	def generateStabilizerData(self, center, scale, num_sides, angle = 0):
		planar_data = Code.generateStabilizerData(self, center, scale, num_sides, angle)
		toric_data = []
		for qubit in planar_data:
			toric_data.append(PlanarToToric(qubit, self.length1, self.length2))
		return toric_data

	

	def hasLogicalError(self, charge_type):
		# Make graph from subgraph of dual using only stabilizers with some, but not all
		# of their member data qubits qith non-zero charge
		# look for self-loops with non-trivial winding number
		# if one exists, we have a logical error

		# store list of positions of nonzero data
		charged_data = []
		for data in self.primal.nodes():
			if data.charge[charge_type] != 0:
				charged_data.append(data.position)

		subgraph_nodes = []
		for measure in self.dual:
			HasTrivialData, HasNontrivialData = False, False
			for type in self.stabilizers:
				if measure.position in self.stabilizers[type]:
					stabilizer = self.stabilizers[type][measure.position]
			
			for data in stabilizer.data:
				if data.position in charged_data:
					HasNontrivialData = True
				else:
					HasTrivialData = True

			if HasTrivialData and HasNontrivialData:
				subgraph_nodes.append(measure)



		Homologies = self.dual.subgraph(subgraph_nodes)
		# Find all cycles, then decide if each is trivial
		for cycle in nx.cycle_basis(Homologies):
			if notTrivialCycle(self, cycle):
				return True

		return False


	def plot_primal(self, charge_type, plot_number, title):
		Primal = plt.figure(plot_number)
		ax = Primal.add_subplot(111, projection='3d')
		d = self.dimension

		for type in self.stabilizers:
			for target in self.stabilizers[type]:
				stabilizer = self.stabilizers[type][target]
				center = stabilizer.center
				color = self.colors[center.type]
				[x,y,z] = center.position[0], center.position[1], center.position[2]
				charge = center.charge[charge_type]
				if charge != 0:
					ax.scatter(x,y,z,marker="*",color=color,s=200*float(charge)/(d-1))


		for data in self.primal.nodes():
			color = self.colors[data.type]
			[x,y,z] = data.position[0], data.position[1], data.position[2]
			charge = data.charge[charge_type]
			if charge != 0:
				ax.scatter(x,y,z,marker="*",color=color,s=200*float(charge)/(d-1))

		for edge in self.primal.edges():
			[x0,y0,z0] = edge[0].position[0], edge[0].position[1], edge[0].position[2]
			[x1,y1,z1] = edge[1].position[0], edge[1].position[1], edge[1].position[2]
			ax.plot([x0,x1], [y0,y1], [z0,z1], color = 'black')

		plt.title(str(title))
		return plt.figure(plot_number)


	def plot_dual(self, charge_type, plot_number, title):

		Dual = plt.figure(plot_number)
		ax = Dual.add_subplot(111, projection='3d')
		d = self.dimension

		for edge in self.dual.edges():
			type = self.dual.get_edge_data(*edge)['type']
			[x0,y0,z0] = edge[0].position[0], edge[0].position[1], edge[0].position[2]
			[x1,y1,z1] = edge[1].position[0], edge[1].position[1], edge[1].position[2]
			ax.plot([x0,x1], [y0,y1], [z0,z1], color = self.colors[type])
		
		for type in self.stabilizers:
			for target in self.stabilizers[type]:
				stabilizer = self.stabilizers[type][target]
				center = stabilizer.center
				color = self.colors[center.type]
				[x,y,z] = center.position[0], center.position[1], center.position[2]
				charge = center.charge[charge_type]
				if charge != 0:
					ax.scatter(x,y,z,marker="*",color=color,s=200*float(charge)/(d-1))

		plt.title(str(title))
		return plt.figure(plot_number)

	def plot_shrunk(self, shrunk_type, charge_type, plot_number, title):
		# plot all edges of shrunk type

		Dual = plt.figure(plot_number)
		ax = Dual.add_subplot(111, projection='3d')
		d = self.dimension

		for edge in self.shrunk[shrunk_type].edges():
			[x0,y0,z0] = edge[0].position[0], edge[0].position[1], edge[0].position[2]
			[x1,y1,z1] = edge[1].position[0], edge[1].position[1], edge[1].position[2]
			ax.plot([x0,x1], [y0,y1], [z0,z1], color = self.colors[shrunk_type])
		
		for type in self.stabilizers:
			if type == shrunk_type:
				continue

			for target in self.stabilizers[type]:
				stabilizer = self.stabilizers[type][target]
				center = stabilizer.center
				color = self.colors[center.type]
				[x,y,z] = center.position[0], center.position[1], center.position[2]
				charge = center.charge[charge_type]
				if charge != 0:
					ax.scatter(x,y,z,marker="*",color=color,s=200*float(charge)/(d-1))

		plt.title(str(title))
		return plt.figure(plot_number)




def ToricCoordinates(qubit,length1, length2 = False):
	if length2 == False:
		length2 = length1
	u = 2 * pi * float(qubit.position[0])/length1
	v = 2 * pi * float(qubit.position[1])/length2
	x = round(float(cos(u)) * (3 + float(cos(v))), 3)
	y = round(float(sin(u)) * (3 + float(cos(v))), 3)
	z = round(float(sin(v)), 3)
	return (x, y, z)

def PlanarToToric(qubit, length1, length2 = False):
	qubit.position = ToricCoordinates(qubit, length1, length2)
	return qubit


def notTrivialCycle(self, cycle):
	# If odd number of crossing points for x or y, 
	# then we have a non-trival homology
	# i.e. a logical operator
	x_crosses, y_crosses = 0, 0

	length1 = self.length1
	if self.length2 != False:
		length2 = self.length2
	else:
		length2 == self.length1

	cycle_length = len(cycle)
	for link in cycle_length:
		node1 = cycle[link]
		node2 = cycle[link%cycle_length]
		type1, type2 = self.dual.node[node1]['type'], self.dual.node[node2]['type']

		pos1 = self.stabilizers[type1][node1.position].planar_coords
		pos2 = self.stabilizers[type1][node2.position].planar_coords

		x_dist, z_dist = abs(pos2[0] - pos1[0]]), abs(pos2[1] - pos1[1]])
		if (pos1[0] < x_dist and (length1 - pos2[0]) < x_dist):
			x_crosses += 1
		elif (pos2[0] < x_dist and (length1 - pos1[0]) < x_dist):
			x_crosses += 1

		if (pos1[1] < y_dist and (length2 - pos2[1]) < y_dist):
			y_crosses += 1
		elif (pos2[1] < y_dist and (length2 - pos1[1]) < y_dist):
			y_crosses += 1

	if x_crosses%2 == 1 or y_crosses%2 == 1:
		return True
	else:
		return False













