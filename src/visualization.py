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

import matplotlib.pyplot as plt
from common import *

def PlotPrimal(code, title, plot_number = 1):
	dim = code.dimension
	Primal = plt.figure(plot_number)
	for node in code.Primal.nodes():
		charge = code.Primal.node[node]['charge']
		# if charge['X'] != 0 and charge['Z'] != 0:
		# 	plt.scatter(*node,marker="*",color='orange',s=200*float(charge['Z'])/(dim-1))
		if charge['Z'] != 0:
			plt.scatter(*node,marker="*",color='blue',s=200*float(charge['Z'])/(dim-1))
		# elif charge['X'] != 0:
		# 	plt.scatter(*node,marker="*",color='red',s=200*float(charge['X'])/(dim-1))
		else:
			plt.scatter(*node,color='black')
	for edge in code.Primal.edges():
		plt.plot((edge[0][0],edge[1][0]),(edge[0][1],edge[1][1]),color='black')
	
	for type in code.External:
		for node in code.External[type]:
			plt.scatter(*node,color=code.colors[type],marker='D')


	plt.title(str(title))

	return plt.figure(plot_number)
	

def PlotDual(code, title, plot_number = 2):
	dim = code.dimension
	Dual = plt.figure(plot_number)
	for type in code.types:
		for qubit in code.Stabilizers[type]:
			charge = code.Stabilizers[type][qubit]['charge']
			if charge['X'] != 0 and charge['Z'] != 0:
				plt.scatter(*qubit,marker="*",color='orange',s=200*float(charge['Z'])/(dim-1))
			elif charge['Z'] != 0:
				plt.scatter(*qubit,marker="*",color='blue',s=200*float(charge['Z'])/(dim-1))
			elif charge['X'] != 0:
				plt.scatter(*qubit,marker="*",color='red',s=200*float(charge['X'])/(dim-1))
	
		for edge in code.Dual[type].edges():
			plt.plot((edge[0][0],edge[1][0]),(edge[0][1],edge[1][1]),color=code.colors[type])
	plt.title(str(title))

	return plt.figure(plot_number)

def PlotPlaquette(code, title, plot_number = 3):
	dim = code.dimension
	Plaquette = plt.figure(plot_number)
	for type in code.types:
		for qubit in code.Stabilizers[type]:
			charge = code.Stabilizers[type][qubit]['charge']
			# if charge['X'] != 0 and charge['Z'] != 0:
			# 	plt.scatter(*qubit,marker="*",color=code.colors[type],s=200*float(charge['Z'])/(dim-1))
			if charge['Z'] != 0:
				plt.scatter(*qubit,marker="*",color=code.colors[type],s=400*float(charge['Z'])/(dim-1))
			# elif charge['X'] != 0:
			# 	plt.scatter(*qubit,marker="*",color=code.colors[type],s=200*float(charge['X'])/(dim-1))
			else:
				plt.scatter(*qubit,color=code.colors[type])

		for edge in code.Dual[type].edges():
			if any(edge[0] in code.External[t] for t in code.types) or any(edge[1] in code.External[t] for t in code.types):
				plt.plot((edge[0][0],edge[1][0]),(edge[0][1],edge[1][1]),color=code.colors[type], linestyle = '--')
			else:
				plt.plot((edge[0][0],edge[1][0]),(edge[0][1],edge[1][1]),color=code.colors[type])
	
	for node in code.Primal.nodes():
		charge = code.Primal.node[node]['charge']
		# if charge['X'] != 0 and charge['Z'] != 0:
		# 	plt.scatter(*node,marker="*",color='orange',s=200*float(charge['Z'])/(dim-1))
		if charge['Z'] != 0:
			plt.scatter(*node,marker="*",color='orange',s=400*float(charge['Z'])/(dim-1))
		# elif charge['X'] != 0:
		# 	plt.scatter(*node,marker="*",color='red',s=200*float(charge['X'])/(dim-1))
		# else:
		plt.scatter(*node,color='black')
	
	for type in code.types:
		for node in code.External[type]:
			plt.scatter(*node,color=code.colors[type],marker='D')


	plt.title(str(title))

	return plt.figure(plot_number)


def PlotShrunk(code, type, title, plot_number = 4):
	dim = code.dimension
	Shrunk = plt.figure(plot_number)

	for t in code.types:
		if t != type:
			for qubit in code.Stabilizers[t]:
				charge = code.Stabilizers[t][qubit]['charge']
				if charge['Z'] != 0:
					plt.scatter(*qubit,marker="*",color=code.colors[t],s=200*float(charge['Z'])/(dim-1))
	
			for qubit in code.External[t]:
				plt.scatter(*qubit,color=code.colors[t],marker='D')

	for edge in code.Dual[type].edges():
			plt.plot((edge[0][0],edge[1][0]),(edge[0][1],edge[1][1]),color=code.colors[type])
	plt.title(str(title))

	return plt.figure(plot_number)





