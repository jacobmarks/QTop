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
        scale = 5*common.euclidean_dist(edge[0], edge[1])

        i = 2

        while unclustered_graph.nodes() != []:
        	print unclustered_graph.nodes(data = True)
        	clusters = GCC_Partition(unclustered_graph, i*scale)
        	for cluster in clusters:
	        	# sys.exit(0)
        		code, unclustered_graph = GCC_Annihilate(cluster, code, unclustered_graph, charge_type, i*scale)
        	i += 1

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
		color_clusters[type] = [node for node in cluster if node[1]['type'] == type]
		color_clusters[type], unclustered_graph, code = GCC_One_Color_Simplify(color_clusters[type], unclustered_graph, code, type, ct)		
		print color_clusters[type]


	color_clusters, unclustered_graph, code = GCC_Two_Color_Simplify(color_clusters, unclustered_graph, code, ct)
	color_clusters, unclustered_graph, code = GCC_Boundary_Simplify(color_clusters, unclustered_graph, code, ct, scale)
	
	return code, unclustered_graph



def GCC_One_Color_Simplify(cc, uc, code, t, ct):
	d = code.dimension
	while len(cc) > 1 :
		start, end = cc[0], cc[1]
		cc, uc, code = GCC_One_Color_Transport(start, end, cc, uc, code, t, ct)
	

	return cc, uc, code


def GCC_One_Color_Transport(s, e, cc, uc, code, t, ct):
	k = s[1]['charge']

	# need to get k1 each time
	# and only remove start if charge = 0 afterward
	d = code.dimension
	t1, t2 = code.complementaryTypes(t)
	dual1 = nx.shortest_path(code.Dual[t1], s[0], e[0])
	num_loops = (len(dual1)-1)/2
	for i in range(num_loops):
		start, end = dual1[2*i], dual1[2*i+2]
		k1, k2 = code.Stabilizers[t][start]['charge'][ct], code.Stabilizers[t][end]['charge'][ct]
		print k1, k2
		print start, end
		# sys.exit(0)

		if start in uc:
			cc.remove((start, {'charge':k1,'type':t}))
			uc.remove_node(start)
			code.Stabilizers[t][start]['charge'][ct] = 0


		dual2 = nx.shortest_path(code.Dual[t2], end, start)
		triangle1 = dual1[2*i:(2*i+2)] + dual2[1:]
		triangle2 =  dual2[:2] + dual1[2*i+1:(2*i+3)]
		loop1 = path.Path(triangle1)
		for data in code.Primal.nodes():
			if loop1.contains_points([data]) == [True]:
				c = code.Primal.node[data]['charge'][ct]
				count = code.Stabilizers[t][start]['data'][data]
				sign = code.Sign(count)
				code.Primal.node[data]['charge'][ct] = (c - sign * k)%d
				break
		loop2 = path.Path(triangle2)
		for data in code.Primal.nodes():
			if loop2.contains_points([data]) == [True]:
				c = code.Primal.node[data]['charge'][ct]
				code.Primal.node[data]['charge'][ct] = (c - sign * k)%d



		end_charge = (k2 + k)%d
		print "END CHARGE", end_charge
		# end_charge = (k2 + sign * k)%d
		code.Stabilizers[t][end]['charge'][ct] = end_charge

		if end in uc.nodes():
			print "BOWSER"
			for m in cc:
				if m[0] == end:
					print m, "MMMMMM"
					cc.remove(m)
			if end_charge == 0:

				uc.remove_node(end)
			else:
				print "SANTA BABY"
				uc.node[end]['charge'] = end_charge
				cc.append((end,{'charge':end_charge,'type':t}))


		else:
			print "HOLLA"
			uc.add_node(end, charge = end_charge, type = t)
			cc.append((end,{'charge':end_charge,'type':t}))	

		print cc, uc.nodes()

	return cc, uc, code


def GCC_Two_Color_Simplify(cc, uc, code, ct):
	print "ITS A SMALL WORLD AFTER ALL"
	for t in code.types:
		if cc[t] == []:
			return cc, uc, code

	d = code.dimension
	ms = {}
	for t in code.types:
		ms[t] = cc[t][0]
	print ms
	# sys.exit(0)
	triangle, uc, code, ct = GCC_Connect(ms, cc, uc, code, ct)
	print triangle, "TRIANGLE"
	uc, code = GCC_Two_Color_Transport(triangle, uc, code, ct)
	return cc, uc, code



def GCC_Connect(ms, cc, uc, code, ct):
	d = code.dimension

	t1, t2, t3 = 'red', 'blue', 'green'
	m1, m2, m3 = ms[t1], ms[t2], ms[t3]
	print ms


	m1_data = code.Stabilizers['red'][m1[0]]['data']
	# print m1_data
	# print code.Stabilizers[t2][m2[0]]['data']
	# sys.exit(0)

	if any(node in code.Stabilizers[t2][m2[0]]['data'] for node in m1_data):
		m2_new = m2
		print "HURRAH"
		# sys.exit(0)

	else:
		for m in code.Stabilizers[t2]:
			if not any(m in code.External[t] for t in code.External):
				if any(node in code.Stabilizers[t2][m]['data'] for node in m1_data):
					break

		c = code.Stabilizers[t2][m]['charge'][ct]
		m2_new = (m, {'charge':c, 'type':t2})

		if m not in uc.nodes():
			uc.add_node(m, charge = c, type = t2)
			cc[t2].append(m2_new)

		cc[t2], uc, code = GCC_One_Color_Transport(m2, m2_new, cc[t2], uc, code, t2, ct)

	m2_data = code.Stabilizers[t2][m2_new[0]]['data']

	if any(node in code.Stabilizers[t3][m3[0]]['data'] for node in m1_data) and any(node in code.Stabilizers[t3][m3[0]]['data'] for node in m2_data) and not any(m in code.External[t] for t in code.External):
		m3_new = m3

	else:
		print "Cherry Tree"
		for m in code.Stabilizers[t3]:
			if m not in code.External[t3]:
				if any(node in code.Stabilizers[t3][m]['data'] for node in m1_data):
					if any(node in code.Stabilizers[t3][m]['data'] for node in m2_data):
						break
		print m
		c = code.Stabilizers[t3][m]['charge'][ct]
		print c
		m3_new = (m, {'charge':c, 'type':t3})
		print m3_new

		if m not in uc.nodes():
			uc.add_node(m, charge = c, type = t3)
			cc[t3].append(m3_new)

		print cc[t3]
		cc[t3], uc, code = GCC_One_Color_Transport(m3, m3_new, cc[t3], uc, code, t3, ct)


	m3_new = cc[t3][0]
	ms[t2], ms[t3] = m2_new, m3_new
	print t2, t3
	print ms
	# sys.exit(0)

	return ms, uc, code, ct



def GCC_Two_Color_Transport(triangle, uc, code, ct):
	print triangle, "HEYYYA"
	# sys.exit(0)
	
	k = triangle['red'][1]['charge']
	print k
	r, g, b = triangle['red'][0], triangle['green'][0], triangle['blue'][0]
	print r, g, b
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
		print c, k, "THIS IS IT"
		# sys.exit(0)
		charge = (c - k)%d

		code.Stabilizers[color][m]['charge'][ct] = charge
		if charge == 0:
			uc.remove_node(m)
		else:
			uc.node[m]['charge'] = charge


	print uc

	return uc, code

def GCC_Boundary_Simplify(cc, uc, code, ct, scale):
	# ts = [t for t in cc if cc[t] != []]

	# if len(ts) == 0:
	# 	return cc, uc, code

	# # Need len 1 case - make this into l2 case and pass along recursively

	# if len(ts) == 2:
	# 	ms = [cc[t] for t in ts][0]
	# 	ext_t = code.complementaryType(ts)
	# 	print ext_t
	# 	for ext in code.External[ext_t]:
	# 		if all(common.euclidean_dist(ext, m[0]) < scale for m in ms):
	# 			cc, uc, code = GCC_Boundary_Two_Color_Transport(ext, cc, uc, code, ct, scale)
	# 			return cc, uc, code
		# also the case where we treat each of the two as a separate one_color_transport

	return cc, uc, code
		
	

def GCC_Boundary_Two_Color_Transport(ext, cc, uc, code, ct, scale):
	d = code.dimension
	[t1, t2] = [t for t in cc if cc[t] != []]
	node1, node2 = cc[t1][0], cc[t2][0]
	k1, k2 = node1[1]['charge'], node2[1]['charge']

	cycle = [node1[0], node2[0], ext, node1[0]]
	loop = path.Path(cycle)

	changed_data = []
	for data in code.Primal.nodes():
		if loop.contains_points([data]) == [True]:
			changed_data.append(data)
			if data in code.Stabilizers[t1][node1[0]]['data']:
				count = code.Stabilizers[t1][node1[0]]['data'][data]
				sign = code.Sign(count)
	# print changed_data
	for data in changed_data:
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
	

			
	cc, uc, code = GCC_Boundary_Simplify(cc, uc, code, ct, scale)
	return cc, uc, code

