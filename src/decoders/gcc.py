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
        scale = 2*common.euclidean_dist(edge[0], edge[1])

        loop_graph = nx.Graph()
        # i = 2
        # while unclustered_graph.nodes() != []
        	# i += 1
        for iter in range(4):
            clusters = GCC_Partition(unclustered_graph, (iter+2)*scale)
            for cluster in clusters:
            	code, unclustered_graph = GCC_Annihilate(cluster, code, unclustered_graph, charge_type, (iter+2)*scale)
        return code

def GCC_Partition(UnclusteredGraph, scale):
# Make edges on Unclustered graph
# between all nodes separated by distance 'scale'
	for node1 in UnclusteredGraph.nodes():
		for node2 in UnclusteredGraph.nodes():
			if node1 != node2:
				dist = common.euclidean_dist(node1, node2)
				if dist <= scale:
					UnclusteredGraph.add_edge(*(node1, node2), weight=dist)
	Clusters = []
	subgraphs = nx.connected_component_subgraphs(UnclusteredGraph)
	for i, sg in enumerate(subgraphs):
		Clusters.append(sg.nodes(data=True))

	return Clusters

def GCC_Annihilate(cluster, code, unclustered_graph, ct, scale):
	color_clusters = {}
	for type in code.types:
		color_clusters[type] = []

	for node in cluster:
		type = node[1]['type']
		color_clusters[type].append(node)

	for type in code.types:
		color_clusters[type], unclustered_graph, code = GCC_One_Color_Simplify(color_clusters[type], unclustered_graph, code, type, ct)
	
	color_clusters, unclustered_graph, code = GCC_Two_Color_Simplify(color_clusters, unclustered_graph, code, ct)
	
	ts = [t for t in color_clusters if color_clusters[t] != []]
	if len(ts) != 0:
		ms = [color_clusters[t] for t in ts][0]
		for type in code.External:
			for ext in code.External[type]:
				print ms
				if all(common.euclidean_dist(ext, m[0]) < scale for m in ms):
					color_clusters, unclustered_graph, code = GCC_Boundary_Simplify(color_clusters, unclustered_graph, code, ct, scale)
					break
			break
	return code, unclustered_graph

def GCC_Boundary_Simplify(cc, uc, code, ct, scale):
	while len([t for t in cc if cc[t] != []]) != 0:

		excited_types = [t for t in cc if cc[t] != []]
		if len(excited_types) == 2:
			cc, uc, code = GCC_Boundary_Two_Color_Transport(cc, uc, code, ct, scale)
	return cc, uc, code

def GCC_Boundary_Two_Color_Transport(cc, uc, code, ct, scale):
	d = code.dimension
	[t1, t2] = [t for t in cc if cc[t] != []]
	node1, node2 = cc[t1][0], cc[t2][0]
	k1, k2 = node1[1]['charge'], node2[1]['charge']
	ext_t = code.complementaryType([t1, t2])
	for ext in code.External[ext_t]:
		if common.euclidean_dist(ext, node1[0]) < scale and common.euclidean_dist(ext, node2[0]) < scale:
			print "Boundary Neutral!!!!"
			cycle = [node1[0], node2[0], ext, node1[0]]
			loop = path.Path(cycle)
			for data in code.Primal.nodes():
				if loop.contains_points([data]) == [True] and data in code.Stabilizers[t1][node1[0]]['data']:
					count = code.Stabilizers[t1][node1[0]]['data'][data]
					sign = code.Sign(count)
					break
			for data in code.Primal.nodes():
				if loop.contains_points([data]) == [True]:
					c = code.Primal.node[data]['charge'][ct]
					code.Primal.node[data]['charge'][ct] = (c - sign * k1)%d
			
			uc.remove_node(node1[0])
			cc[t1].remove(node1)
			code.Stabilizers[t1][node1[0]]['charge'][ct] = 0
			charge = (k2 - k1)%d

			code.Stabilizers[t2][node2[0]]['charge'][ct] = charge
			if charge == 0:
				uc.remove_node(node2[0])
				cc[t2].remove(node2)
			else:
				uc.node[node2[0]]['charge'] = charge
				cc[t2][0]['charge'] = charge
			break
	return cc, uc, code


def GCC_One_Color_Simplify(cc, uc, code, t, ct):
	d = code.dimension
	while len(cc) > 1 :
		start, end = cc[0], cc[1]
		print start, end
		k1, k2 = start[1]['charge'], end[1]['charge']
		print k1, k2
		cc.remove(start)
		uc.remove_node(start[0])
		code.Stabilizers[t][start[0]]['charge'][ct] = 0

		cc, uc, code = GCC_One_Color_Transport(start, end, cc, uc, code, t, ct)
	

	return cc, uc, code


def GCC_One_Color_Transport(s, e, cc, uc, code, t, ct):
	k1, k2 = s[1]['charge'], e[1]['charge']
	d = code.dimension
	t1, t2 = code.complementaryTypes(t)
	dual1 = nx.shortest_path(code.Dual[t1], s[0], e[0])
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
				count = code.Stabilizers[t][s[0]]['data'][data]
				sign = code.Sign(count)
				code.Primal.node[data]['charge'][ct] = (c - sign * k1)%d
				break
		loop2 = path.Path(triangle2)
		for data in code.Primal.nodes():
			if loop2.contains_points([data]) == [True]:
				c = code.Primal.node[data]['charge'][ct]
				code.Primal.node[data]['charge'][ct] = (c - sign * k1)%d

		cc.remove(e)
		if (k2 - sign * k1)%d == 0:
			uc.remove_node(e[0])
			code.Stabilizers[t][e[0]]['charge'][ct] = 0
			print 'YEES'
		else:
			code.Stabilizers[t][e[0]]['charge'][ct] = (k2 - sign * k1)%d
			uc.node[e[0]]['charge'] = (k2 - sign * k1)%d
			new_end = (e[0],{'charge':(k2 - sign * k1)%d,'type':t})
			cc.append(new_end)

	return cc, uc, code


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
			sign = code.Sign(count, sides)
			code.Primal.node[data]['charge'][ct] = (c - sign * k)%d
	

	uc.remove_node(r)
	code.Stabilizers['red'][r]['charge'][ct] = 0

	for color in ['green', 'blue']:
		m = triangle[color][0]
		c = triangle[color][1]['charge']
		charge = (c - k)%d

		code.Stabilizers[color][m]['charge'][ct] = charge
		if charge == 0:
			uc.remove_node(m)
		else:
			uc.node[m]['charge'] = charge



	return uc, code



