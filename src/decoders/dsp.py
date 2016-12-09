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
import itertools

class DSP_decoder(decoder):

    def __call__(self, code):
        matching = self.algorithm()
        # for charge_type in ['X','Z']:
        for charge_type in ['Z']:
            code = matching(code, charge_type)

        code = reset_measures(code)
        return code

    def algorithm(self):
        return DSP()

def DSP_Matching(Syndrome, External, dim):

    # Fully connect check operators
    for check1 in Syndrome.nodes():
        for check2 in Syndrome.nodes():
            if check1 != check2:
                weight = - common.euclidean_dist(check1, check2) +.1
                Syndrome.add_edge(*(check1, check2), weight=weight)

    # Generate Boundary Graph
    External_Graph = nx.Graph()

    for m in Syndrome.nodes():
        ext = DSP_AssociatedExternal(m, External)
        External_Graph.add_node(ext)
        weight = - common.euclidean_dist(m, ext)
        Syndrome.add_edge(*(m, ext), weight=weight)

    # Ensure even number of elements in Syndrome
    # so min weight matching can proceed successfully
    if len(Syndrome.nodes()) % 2 != 0:
        far_ext = min(External_Graph.nodes(), key = lambda x:Syndrome.get_edge_data(*Syndrome.edges(x)[0])['weight'])
        External_Graph.remove_node(far_ext)
        Syndrome.remove_node(far_ext)

    # Connect External Nodes
    edges = itertools.combinations(External_Graph,2)
    Syndrome.add_edges_from(edges, weight = 0)

    TempMatching = nx.max_weight_matching(Syndrome, maxcardinality=True)

    Matching = {}
    # each edge appears twice in TempMatching
    # Should only appear once in Matching
    for node, neighbor in TempMatching.items():
        if neighbor not in Matching and node not in Matching:
            if node not in External or neighbor not in External:
                Matching[node] = neighbor

    print Matching
    return Matching

class DSP(matching_algorithm):

    def __init__(self):
        pass

    def __call__(self, code, charge_type):
        errors = {}
        for type in code.types:
            errors[type] = code.Syndrome(type, charge_type)

        shrunk_errs, shrunk_exts, matches = {}, {}, {}
        loops_graph = nx.Graph()
        for t1 in code.types:
            [t2, t3] = code.complementaryTypes(t1)
            shrunk_errs[t1] = nx.union(errors[t2], errors[t3])
            shrunk_exts[t1] = code.External[t2] + code.External[t3]
            print t1, shrunk_exts[t1]
            matches[t1] = DSP_Matching(shrunk_errs[t1], shrunk_exts[t1], 2)

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
        
        print hasConnectedBoundaries(code, loops_graph, Exts)
        print "EDGES", loops_graph.edges()
        # print nx.has_path(loops_graph,(4.5, 2.598), (6.0, 1.732))


        code, loops_graph = correctLoops(code, loops_graph, charge_type)
        # print "EDGES", loops_graph.edges()
        while hasConnectedBoundaries(code, loops_graph, Exts):
            print "EDGES", loops_graph.edges()
            # print "GOTCHA"
            node, ext1, ext2 = connectedBoundaries(loops_graph, Exts)
            code, loops_graph = makeBoundLoop(code, loops_graph, node, ext1, ext2)
            code, loops_graph = correctLoops(code, loops_graph, charge_type)
            print "REMOVED", loops_graph.edges()
            # print "UN-CONNECTEDS"
            print node, ext1, ext2
            print loops_graph.edges()
        return code

def correctLoops(code, loops_graph, charge_type):
    while nx.cycle_basis(loops_graph) != []:
            cycle = nx.cycle_basis(loops_graph)[0]
            print "CYCLE", cycle
            loop = path.Path(cycle)
            for data in code.Primal.nodes():
                if loop.contains_points([data]) == [True]:
                    charge = code.Primal.node[data]['charge'][charge_type]
                    code.Primal.node[data]['charge'][charge_type] = (charge + 1)%2

            l = len(cycle)
            for i in range(l):
                n1, n2 = cycle[i], cycle[(i+1)%l]
                loops_graph.remove_edge(*(n1,n2))

    return code, loops_graph

def DSP_Path(DualGraph, terminal1, terminal2):
    return nx.shortest_path(DualGraph, terminal1, terminal2)

def DSP_AssociatedExternal(int_node, external_nodes):
    return min(external_nodes, key = lambda x:common.euclidean_dist(int_node, x))

def hasConnectedBoundaries(code, loops_graph, Exts):
    for node in loops_graph.nodes():
        for ext1 in Exts:
            for ext2 in Exts:
                if ext1 in loops_graph.nodes() and ext2 in loops_graph.nodes() and node not in Exts:
                    if nx.has_path(loops_graph,node,ext1) and nx.has_path(loops_graph,node,ext2) and ext1 != ext2:
                        print "LOOP"
                        print nx.shortest_path(loops_graph, node, ext1)
                        print nx.shortest_path(loops_graph, node, ext2)
                        return True
                        # for t in code.types:
                        #     if ext1 in code.Stabilizers[t]:
                        #         t1 = t
                        #     if ext2 in code.Stabilizers[t]:
                        #         t2 = t
                        #     if node in code.Stabilizers[t]:
                        #         t_n = t
                        # if t1 != t2 and t2!=t_n and t_n!=t1:
                        #     return True
                        # # elif t1 == t2 and t2 == t_n and t_n == t1:
                        # elif t1 == t2:
                        #     return True
    return False

def connectedBoundaries(loops_graph, Exts):
    for node in loops_graph.nodes():
        for ext1 in Exts:
            for ext2 in Exts:
                if ext1 in loops_graph.nodes() and ext2 in loops_graph.nodes() and node not in Exts:
                    if nx.has_path(loops_graph,node,ext1) and nx.has_path(loops_graph,node,ext2) and ext1 != ext2:
                        return node, ext1, ext2

def makeBoundLoop(code, loops_graph, int, e1, e2):
    for t1 in code.External:
        if e1 in code.External[t1]:
            break
    for t2 in code.External:
        if e2 in code.External[t2]:
            break
    if t1 == t2:
        code, loops_graph = makeSingleBoundLoop(code, loops_graph, int, e1, e2)
    else:
        code, loops_graph = makeDoubleBoundLoop(code, loops_graph, e1, e2, t1, t2)
    return code, loops_graph

def makeSingleBoundLoop(code, loops_graph, int, e1, e2):
    loops_graph.add_edge(*(e1, e2))
    return code, loops_graph

def makeDoubleBoundLoop(code, loops_graph, e1, e2, t1, t2):
    t3 = code.complementaryType([t1,t2])
    for ext1 in code.External[t1]:
        if any(ext2 in code.External[t2] for ext2 in code.Dual[t3].neighbors(ext1)):
            break
    for ext2 in code.External[t2]:
        if ext2 in code.Dual[t3].neighbors(ext1):
            break
    loops_graph.add_edge(*(e1,ext1))
    loops_graph.add_edge(*(e2,ext2))
    loops_graph.add_edge(*(ext1,ext2))

    return code, loops_graph

