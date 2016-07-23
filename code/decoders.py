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

from common import *
from math import *
import networkx as nx


# common base class for decoding functions
# as well as decoding based on 
# Edmunds' Blossom algorithm and
# Hard Decision Renormalization Group Methods

def syndrome(code, type, charge_type = 'Z'):
    Syndrome = nx.Graph()

    # Find all non-trivial check operators
    for measure in code.stabilizers[type]:
    	center = code.stabilizers[type][measure].center
    	if center.charge[charge_type] != 0:
    		Syndrome.add_node(center)

    return Syndrome


## Call with Decoder(code, syndrome)
class Decoder:

	def __init__(self):
		self.code = code

	def __call__(self, code, syndrome):
		correction = []
		return correction



# Min Weight Perfect Matching based in Edmunds' Blossom algorithm
# only works for qubit systems

def MWPM(code, syndrome, type, charge_type):
	for check1 in syndrome.nodes():
		for check2 in syndrome.nodes():
			if check1 != check2:
				weight = - distance(check1, check2, type)
				syndrome.add_edge(*(check1, check2), weight=weight)

	if code.geometry != 'toric':
		syndrome = addExternalNodes(code, syndrome, type)

	temp_matching = nx.max_weight_matching(syndrome, maxcardinality=True)
	# remove repeats
	matching = {}
	for node, neighbor in temp_matching:
		if neighbor not in matching:
			if node in self.dual.nodes() or neighbor in self.dual.nodes():
				matching[node] = neighbor

	return matching


def addExternalNodes(code, syndrome, type):
	external = nx.Graph()

	for node in syndrome.nodes():
		external_node = code.associatedExternal(node, type)
		external.add_node(external_node)
		weight = - distance(node, external_node, type)
		syndrome.add_edge(*(node, external_node))

	# Ensure even number of elements in syndrome
    # so min weight matching can proceed successfully

    if len(syndrome.nodes()) % 2 != 0:
    	removed_node = external.nodes()[0]
    	min_weight_edge = syndrome.edges(removed_node)[0]
    	min_weight = syndrome.get_edge_data(*min_weight_edge)['weight']
    	for node in external.nodes():
    		edge = syndrome.edges(node)[0]
    		weight = syndrome.get_edge_data(*edge)['weight']
    		if weight < min_weight:
    			removed_node = node
    			min_weight = weight

    	external.remove_node(removed_node)
    	syndrome.removed_node(removed_node)

    for ext1 in external.nodes():
    	for ext2 in external.nodes():
    		if ext1 != ext2:
    			syndrome.add_edge(*(ext1, ext2), weight = 0)

    return syndrome









