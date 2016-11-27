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
from decoders import *
from matplotlib import path
from math import floor


def GSP_Path(DualGraph, start, end):
    terminal1, terminal2 = start[0], end[0]
    if terminal1 in DualGraph.nodes() and terminal2 in DualGraph.nodes():
        return nx.shortest_path(DualGraph, terminal1, terminal2)
    elif terminal1 not in DualGraph.nodes():
        ext = terminal1
        terminal1 = DSP_AssociatedInternal(ext, DualGraph.nodes())
        return DSP_Path(DualGraph, terminal2, terminal1) + [ext]
    else:
        return DSP_Path(DualGraph, terminal2, terminal1)


################ GCC ###################

class GCC_decoder(decoder):

    def __call__(self, code):
        matching = self.algorithm()
        # really want ['X', Z'] but this is easier for testing
        for charge_type in ['Z']:
            code = matching(code, charge_type)

        # code = reset_measures(code)

        return code

    def algorithm(self):
        return GCC()


class GCC(matching_algorithm):

    def __init__(self):
        pass

    def __call__(self, code, charge_type):
        l,d = code.depth, code.dimension
        s = {}
        for type in code.types:
            s[type] = code.Syndrome(type, charge_type)

        unclustered_graph = nx.union(s['green'], nx.union(s['red'], s['blue']))

        for edge in code.Primal.edges():
            break
        scale = 2*euclidean_dist(edge[0], edge[1])

        loop_graph = nx.Graph()
        # for iter in range(int(float(l/3))):
        for iter in range(4):
            clusters = GCC_Partition(unclustered_graph, (iter+2)*scale)
            for cluster in clusters:
            	code, unclustered_graph = GCC_Annihilate(cluster, code, unclustered_graph, charge_type)
  
        return code

def GCC_Partition(UnclusteredGraph, scale):
# Make edges on Unclustered graph
# between all nodes separated by distance 'scale'
	for node1 in UnclusteredGraph.nodes():
		for node2 in UnclusteredGraph.nodes():
			if node1 != node2:
				dist = euclidean_dist(node1, node2)
				if dist <= scale:
					UnclusteredGraph.add_edge(*(node1, node2), weight=dist)
	Clusters = []
	subgraphs = nx.connected_component_subgraphs(UnclusteredGraph)
	for i, sg in enumerate(subgraphs):
		Clusters.append(sg.nodes(data=True))

	return Clusters

def GCC_Annihilate(cluster, code, unclustered_graph, ct):
	color_clusters = {}
	for type in code.types:
		color_clusters[type] = []

	for node in cluster:
		type = node[1]['type']
		color_clusters[type].append(node)

	for type in code.types:
		color_clusters[type], unclustered_graph, code = GCC_One_Color_Simplify(color_clusters[type], unclustered_graph, code, type, ct)
	
	color_clusters[type], unclustered_graph, code = GCC_Two_Color_Simplify(color_clusters, unclustered_graph, code, ct)

	return code, unclustered_graph


def GCC_One_Color_Simplify(cc, uc, code, t, ct):
	d = code.dimension
	while len(cc) > 1 :
		start, end = cc[0], cc[1]
		k1, k2 = start[1]['charge'], end[1]['charge']
		cc.remove(start)
		uc.remove_node(start[0])
		code.Stabilizers[t][start[0]]['charge'][ct] = 0
		
		cc.remove(end)
		if (k1 - k2)%d == 0:
			uc.remove_node(end[0])
			code.Stabilizers[t][end[0]]['charge'][ct] = 0
		else:
			code.Stabilizers[t][end[0]]['charge'][ct] = (k1 - k2)%d
			uc.node[end[0]]['charge'] = (k1 - k2)%d
			new_end = (end[0],{'charge':(k1 - k2)%d,'type':t})
			cc.append(new_end)

		code = GCC_One_Color_Transport(start, end, code, t, ct)
	

	return cc, uc, code


def GCC_One_Color_Transport(start, end, code, t, ct):
	k = start[1]['charge']
	d = code.dimension
	t1, t2 = code.complementaryTypes('t')
	dual1 = nx.shortest_path(code.Dual[t1], start[0], end[0])
	num_loops = (len(dual1)-1)/2
	for i in range(num_loops):
		start, end = dual1[2*i +2], dual1[2*i]
		dual2 = nx.shortest_path(code.Dual[t2], start, end)
		triangle1 = dual1[2*i:(2*i+2)] + dual2[1:]
		triangle2 =  dual2[:2] + dual1[2*i+1:(2*i+3)]
		loop1 = path.Path(triangle1)
		for data in code.Primal.nodes():
			if loop1.contains_points([data]) == [True]:
				c = code.Primal.node[data]['charge'][ct]
				code.Primal.node[data]['charge'][ct] = (c+k)%d

		sign = code.transportSign(triangle2[1], triangle2[2])
		loop2 = path.Path(triangle2)
		for data in code.Primal.nodes():
			if loop2.contains_points([data]) == [True]:
				c = code.Primal.node[data]['charge'][ct]
				code.Primal.node[data]['charge'][ct] = (c+ sign * k)%d
	return code


def GCC_Two_Color_Simplify(cc, uc, code, ct):
	for t in code.types:
		if cc[t] == []:
			return cc, uc, code

	d = code.dimension
	ms = {}
	for t in code.types:
		ms[t] = cc[t][0]

	triangle, uc, code, ct = GCC_Connect(ms, uc, code, ct)

	uc, code = GCC_Two_Color_Transport(triangle, uc, code, ct)


	return cc, uc, code



def GCC_Connect(ms, uc, code, ct):

	t1, t2, t3 = 'red', 'blue', 'green'
	m1, m2, m3 = ms[t1], ms[t2], ms[t3]
	path1 = nx.shortest_path(code.Dual[t1], m2[0], m3[0])
	path2 = nx.shortest_path(code.Dual[t2], m1[0], m3[0])


	if len(path1) > 2:

		
		start1, end1 = m1, (path2[-2], {'type':t2, 'charge':0})
		code = GCC_One_Color_Transport(start1, end1, code, t1, ct)

		k1 = start1[1]['charge']
		uc.remove_node(start1[0])
		code.Stabilizers[t1][start1[0]]['charge'][ct] = 0

		uc.add_node(end1[0], type = t1, charge = k1)
		code.Stabilizers[t1][end1[0]]['charge'][ct] = k1
	else:
		end1 = m1

	if len(path2) > 2:

		start2, end2 = m2, (path1[-2], {'type':t1, 'charge':0})
		code = GCC_One_Color_Transport(start2, end2, code, t2, ct)
		
		k2 = start2[1]['charge']
		uc.remove_node(start2[0])
		code.Stabilizers[t2][start2[0]]['charge'][ct] = 0

		uc.add_node(end2[0], type = t2, charge = k2)
		code.Stabilizers[t2][end1[0]]['charge'][ct] = k2
	else:
		end2 = m2

	ms[t1], ms[t2] = end1, end2


	return ms, uc, code, ct



def GCC_Two_Color_Transport(triangle, uc, code, ct):
	k = triangle['red'][1]['charge']

	r, g, b = triangle['red'][0], triangle['green'][0], triangle['blue'][0]
	d = code.dimension
	sides = code.types['red']['sides']

	cycle = [r, g, b, r]
	loop = path.Path(cycle)
	for data in code.Primal.nodes():
		if loop.contains_points([data]) == [True]:
			c = code.Primal.node[data]['charge'][ct]
			count = code.Stabilizers['red'][cycle[0]]['data'][data]
			sign = Sign(count, sides)
			code.Primal.node[data]['charge'][ct] = (c - sign * k)%d
	

	uc.remove_node(r)
	code.Stabilizers['red'][r]['charge'][ct] = 0

	for color in ['green', 'blue']:
		m = triangle[color][0]
		c = triangle[color][1]['charge']
		sign = code.transportSign(r, m)
		charge = (c - sign * k)%d

		code.Stabilizers[color][m]['charge'][ct] = charge
		if charge == 0:
			uc.remove_node(m)
		else:
			print 'else'
			uc.node[m]['charge'] = charge

	print uc


	return uc, code



