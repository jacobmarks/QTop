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

# from common import *
from decoders import *
from matplotlib import path
from math import floor
# from common import *
import sys
sys.path.append('../../')
from src import common
# from common import *
import networkx as nx
import numpy as np



 ############# Minimum Weight Perfect Matching ###############

# Given a network of weighted edges, finds the 
# minimum weight perfect matching of the nodes.

class MWPM_decoder(surface_decoder):

    def __call__(self, code):
        return surface_decoder.__call__(self, code)

    def algorithm(self):
        return MWPM()


class MWPM(matching_algorithm):

    def __init__(self):
        pass

    def __call__(self, code, Syndrome, type, charge_type):
        matching = MinWeightMatching(code, Syndrome, type, charge_type)
        code = Recovery(code, matching, type, charge_type)
        return code



def MinWeightMatching(code, Syndrome, type, charge_type):
    dim = code.dimension

    # Fully connect check operators
    for check1 in Syndrome.nodes():
        for check2 in Syndrome.nodes():
            if check1 != check2:
                # weight = - code.distance(type, check1, check2)
                weight = - common.euclidean_dist(check1, check2)
                Syndrome.add_edge(*(check1, check2), weight=weight)

    # Generate Boundary Graph
    External_Graph = nx.Graph()

    for node in Syndrome.nodes():
        charge = Syndrome.node[node]['charge']
        external_node = AssociatedExternal(node, code.Dual[type], code.External[type])
        External_Graph.add_node(external_node, charge=(-charge) % dim)
        # weight = -code.distance(type, node, external_node)
        weight = - common.euclidean_dist(node, external_node)
        Syndrome.add_edge(*(node, external_node), weight=weight)

    # Ensure even number of elements in Syndrome
    # so min weight matching can proceed successfully
    if len(Syndrome.nodes()) % 2 != 0:
        removed_external = External_Graph.nodes()[0]
        edge = Syndrome.edges(removed_external)[0]
        min_weight = Syndrome.get_edge_data(*edge)['weight']
        for external_node in External_Graph.nodes():
            edge = Syndrome.edges(external_node)[0]
            weight = Syndrome.get_edge_data(*edge)['weight']
            if weight < min_weight:
                removed_external = external_node
                min_weight = weight

        External_Graph.remove_node(removed_external)
        Syndrome.remove_node(removed_external)

    # Connect External Nodes
    for ext1 in External_Graph:
        for ext2 in External_Graph:
            if ext1 != ext2:
                Syndrome.add_edge(*(ext1, ext2), weight=0)

    TempMatching = nx.max_weight_matching(Syndrome, maxcardinality=True)

    Matching = {}
    # each edge appears twice in TempMatching
    # Should only appear once in Matching
    for node, neighbor in TempMatching.items():
        if neighbor not in Matching:
            if node in code.Dual[type].nodes() or neighbor in code.Dual[type].nodes():
                Matching[node] = neighbor

    return Matching


# Recovery Operations
# Generate recovery chains to correct for errors during code cycle. 

def Recovery(code, Matching, type, charge_type):
    for terminal1, terminal2 in Matching.items():
        if terminal1 in code.Dual[type].nodes() and terminal2 in code.Dual[type].nodes():
            code = InternalRecovery(code, terminal1, terminal2, type, charge_type)
        else:
            code = BoundaryRecovery(code, terminal1, terminal2, type, charge_type)

    return code

def InternalRecovery(code, terminal1, terminal2, type, charge_type):
    measure_chain = nx.shortest_path(code.Dual[type], terminal1, terminal2)
    chain_length = nx.shortest_path_length(code.Dual[type], terminal1, terminal2)
    for link in range(chain_length):
        vertex1 = measure_chain[link]
        vertex2 = measure_chain[link + 1]
        for data_qubit in code.Stabilizers[type][vertex1]['data']:
            if data_qubit in code.Stabilizers[type][vertex2]['data']:
                prior_charge = code.Primal.node[data_qubit]['charge'][charge_type]
                code.Primal.node[data_qubit]['charge'][charge_type] = (prior_charge + 1) % 2

    return code

def BoundaryRecovery(code, terminal1, terminal2, type, charge_type):
    if terminal1 not in code.Stabilizers[type]:
        data_near_boundary = code.External[type][terminal1]['data']
        prior_charge = code.Primal.node[data_near_boundary]['charge'][charge_type]
        code.Primal.node[data_near_boundary]['charge'][charge_type] = (prior_charge + 1) % 2

        terminal1 = code.External[type][terminal1]['measure']
        InternalRecovery(code, terminal1, terminal2, type, charge_type)
    else:
        data_near_boundary = code.External[type][terminal2]['data']
        prior_charge = code.Primal.node[data_near_boundary]['charge'][charge_type]
        code.Primal.node[data_near_boundary]['charge'][charge_type] = (prior_charge + 1) % 2
        terminal2 = code.External[type][terminal2]['measure']
        InternalRecovery(code, terminal1, terminal2, type, charge_type)

    return code


def AssociatedExternal(node, Dual, External):
    associate = External.iterkeys().next()
    min_dist = nx.shortest_path_length(Dual, node, External[associate]['measure']) + 1

    for candidate in External:
        distance = nx.shortest_path_length(Dual, node, External[candidate]['measure']) + 1
        if distance < min_dist:
            min_dist = distance
            associate = candidate
    return associate