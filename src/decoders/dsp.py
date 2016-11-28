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

from decoders import *
from matplotlib import path
from math import floor
import sys
sys.path.append('../../')
from src import common
import networkx as nx
import numpy as np



def DSP_Matching(Syndrome, External, dim):

    # Fully connect check operators
    for check1 in Syndrome.nodes():
        for check2 in Syndrome.nodes():
            if check1 != check2:
                weight = - common.euclidean_dist(check1, check2)
                Syndrome.add_edge(*(check1, check2), weight=weight)

    # Generate Boundary Graph
    External_Graph = nx.Graph()

    for node in Syndrome.nodes():
        charge = Syndrome.node[node]['charge']
        external_node = DSP_AssociatedExternal(node, External)
        External_Graph.add_node(external_node, charge=(-charge) % dim)
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
        if neighbor not in Matching and node not in Matching:
            if node not in External or neighbor not in External:
                Matching[node] = neighbor

    return Matching


def DSP_AssociatedExternal(int_node, external_nodes):
    associate = external_nodes[0]
    min_dist = common.euclidean_dist(int_node, associate)

    for candidate in external_nodes:
        distance = common.euclidean_dist(int_node, candidate)
        if distance < min_dist:
            min_dist = distance
            associate = candidate
    return associate

def DSP_AssociatedInternal(ext_node, internal_nodes):
    associate = internal_nodes[0]
    min_dist = common.euclidean_dist(ext_node, associate)

    for candidate in internal_nodes:
        distance = common.euclidean_dist(ext_node, candidate)
        if distance < min_dist:
            min_dist = distance
            associate = candidate
    return associate

def DSP_Path(DualGraph, terminal1, terminal2):
    if terminal1 in DualGraph.nodes() and terminal2 in DualGraph.nodes():
        return nx.shortest_path(DualGraph, terminal1, terminal2)
    elif terminal1 not in DualGraph.nodes():
        ext = terminal1
        terminal1 = DSP_AssociatedInternal(ext, DualGraph.nodes())
        return DSP_Path(DualGraph, terminal2, terminal1) + [ext]
    else:
        return DSP_Path(DualGraph, terminal2, terminal1)

class DSP(matching_algorithm):

    def __init__(self):
        pass

    def __call__(self, code, charge_type):
        errors = {}
        for type in code.types:
            errors[type] = code.Syndrome(type, charge_type)
            

        shrunk_errors, matches = {}, {}
        loops_graph = nx.Graph()
        for t1 in code.types:
            comps = [errors[t2] for t2 in code.types if t2 != t1]


            shrunk_errors[t1] = nx.union(comps[0],comps[1])
            matches[t1] = DSP_Matching(shrunk_errors[t1], code.External[t1], 2)

            for start in matches[t1]:
                end = matches[t1][start]
                chain = DSP_Path(code.Dual[t1], start, end)
                links = len(chain) -1

                for i in range(links):
                    node1, node2 = chain[i], chain[i+1]
                    edge = (node1, node2)
                    if edge in loops_graph.edges():
                        loops_graph.remove_edge(*edge)
                    else:
                        loops_graph.add_edge(*edge)
        Exts = code.External['red']+code.External['blue']+code.External['green']
        code, loops_graph = correctLoops(code, loops_graph, charge_type)
        for node1 in loops_graph.nodes():
            if node1 in Exts:
                for node2 in loops_graph.nodes():
                    if node2 in Exts and node2 != node1:
                        if nx.has_path(loops_graph,node1,node2):
                            loops_graph.add_edge(*(node1, node2))
        code, loops_graph = correctLoops(code, loops_graph, charge_type)
        return code

class DSP_decoder(decoder):

    def __call__(self, code):
        matching = self.algorithm()
        for charge_type in ['X','Z']:
            code = matching(code, charge_type)

        code = reset_measures(code)
        return code

    def algorithm(self):
        return DSP()

def correctLoops(code, loops_graph, charge_type):
    while nx.cycle_basis(loops_graph) != []:
            cycle = nx.cycle_basis(loops_graph)[0]
            loop = path.Path(cycle)
            for data in code.Primal.nodes():
                if loop.contains_points([data]) == [True]:
                    charge = code.Primal.node[data]['charge']
                    code.Primal.node[data]['charge'][charge_type] = (charge[charge_type] + 1)%2

            l = len(cycle)
            for i in range(l):
                n1, n2 = cycle[i], cycle[(i+1)%l]
                loops_graph.remove_edge(*(n1,n2))

    return code, loops_graph