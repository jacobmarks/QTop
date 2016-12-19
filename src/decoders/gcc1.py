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
sys.path.append('../')
sys.path.append('mcl/')
from mcl_clustering import networkx_mcl
sys.path.append('community/')
from community_louvain import *
sys.path.append('../../')
from src import common, visualization
import networkx as nx
import numpy as np

################ GCC ###################

class GCC_decoder(decoder):

    def __call__(self, code):
        matching = self.algorithm()
        # for ct in ['X','Z']:
        for ct in ['Z']:
            code = matching(code, ct)
        return code


    def algorithm(self):
        return GCC()


class GCC(matching_algorithm):

    def __init__(self):
        pass

    def __call__(self, code, ct):
        l,d = code.depth, code.dimension
        s = {}
        for type in code.types:
            s[type] = code.Syndrome(type, ct)

        uc = nx.union(s['green'], nx.union(s['red'], s['blue']))
        for edge in code.Primal.edges():
            break
        scale = 2*common.euclidean_dist(edge[0], edge[1])

        i = 2
        last = []
        while uc.nodes() != []:
        	if last != uc.nodes():
	        	visualization.PlotPlaquette(code, "Logical Error", i+2)
        	last = uc.nodes()
        	clusters = GCC_Partition(uc, i*scale)
        	for cluster in clusters:
        		# print cluster
        		code, uc = GCC_Annihilate(cluster, code, uc, ct, (i+1)*scale)
        	i += 1
        	
        return code

def center(cluster):
	X = round(np.mean([x for (x,y) in cluster]),3)
	Y = round(np.mean([y for (x,y) in cluster]),3)
	return (X,Y)

def closest_to_point(nodes,pt):
    return min(nodes, key = lambda x:common.euclidean_dist(pt, x))

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
	
	# print UnclusteredGraph.edges()
	# print UnclusteredGraph.nodes()
	# if len(UnclusteredGraph.nodes()) == 1:
	# 	return [[UnclusteredGraph.nodes(data=True)[0]]]
	# parts = best_partition(UnclusteredGraph)
	# print parts
	# for c in parts:
	# 	if not isinstance( parts[c], ( int, long ) ):
	# 		return []
	# # print UnclusteredGraph.nodes()
	# values = [parts.get(node) for node in UnclusteredGraph.nodes()]
	# for i in range(max(values)+1):
	# 	Clusters.append([])
	# for node in UnclusteredGraph.nodes(data=True):
	# 	if node[0] in parts:
	# 		i = parts[node[0]]
	# 	Clusters[i].append(node)

	# print Clusters
	# sys.exit(0)
	subgraphs = nx.connected_component_subgraphs(UnclusteredGraph)
	for i, sg in enumerate(subgraphs):
		Clusters.append(sg.nodes(data=True))

	# print UnclusteredGraph.nodes()
	# _, clust_dict = networkx_mcl(UnclusteredGraph, expand_factor = 2)
	# print clust_dict
	# # print UnclusteredGraph.nodes()
	# while clust_dict != {}:
	# 	for key in clust_dict:
	# 		c = clust_dict[key]
	# 		clust = [UnclusteredGraph.nodes(data=True)[i] for i in c]
	# 		Clusters.append(clust)
	# 		for val in c:
	# 			del clust_dict[val]
	# 		break
	# print Clusters
	# sys.exit(0)


	return Clusters

def GCC_Annihilate(cluster, code, uc, ct, scale):
	# print cluster
	# sys.exit(0)
	nodes = [node[0] for node in cluster]
	cntr = center(nodes)
	cc = {}
	for type in code.types:
		cc[type] = [node for node in cluster if node[1]['type'] == type]
		cc[type], uc, code = GCC_One_Color_Simplify(cc[type], uc, code, type, ct, cntr, scale)		

	cc, uc, code = GCC_Two_Color_Simplify(cc, uc, code, ct, cntr)
	cc, uc, code = GCC_Boundary_Simplify(cc, uc, code, ct, cntr, scale)
	
	return code, uc



def GCC_One_Color_Simplify(cc, uc, code, t, ct, cntr, scale):
	d = code.dimension

	# center for each color
	nodes = [node[0] for node in cc]
	cntr = center(nodes)


	mid_m = closest_to_point(code.Stabilizers[t],cntr)
	# see if each node closer to boundary than mid_m
	# for m in cc:
		# partner = closest_to_point(cc[0], node)
		# if any(common.euclidean_dist(ext, m[0]) < (scale+ .1) for ext in code.External[t]):
		# 	ext = closest_to_point(code.External[t],m[0])
		# 	cc, uc, code = GCC_One_Color_Transport(m, (ext, {'charge':0,'type':t}), cc, uc, code, t, ct)




	while any(m[0]!=mid_m for m in cc):
		for start in cc:
			if start[0] != mid_m:
				break
		if any(end[0]==mid_m for end in cc):
			for end in cc:
				if end[0]==mid_m:
					break
		else:
			end = (mid_m, {'charge':0,'type':t})

		cc, uc, code = GCC_One_Color_Transport(start, end, cc, uc, code, t, ct)
	return cc, uc, code


def GCC_One_Color_Transport(s, e, cc, uc, code, t, ct):
	end_sign = 1
	k, k_end = s[1]['charge'], e[1]['charge']

	d = code.dimension
	t1, t2 = code.complementaryTypes(t)
	dual1 = nx.shortest_path(code.Dual[t1], s[0], e[0])
	if s[0] in uc:
		cc.remove((s[0], {'charge':k,'type':t}))
		uc.remove_node(s[0])
		code.Stabilizers[t][s[0]]['charge'][ct] = 0

	num_loops = (len(dual1)-1)/2
	for i in range(num_loops):
		start, end = dual1[2*i], dual1[2*i+2]
		k1, k2 = code.Stabilizers[t][start]['charge'][ct], code.Stabilizers[t][end]['charge'][ct]
		if any(start in code.External[type] for type in code.External) and any(end in code.External[type] for type in code.External):
			mid = dual1[2*i+1]
			for start_type in code.External:
				if start in code.External[start_type]:
					break
			for end_type in code.External:
				if start in code.External[end_type]:
					break
			for mid_type in code.Stabilizers:
				if mid in code.Stabilizers[mid_type]:
					break
			for data in code.Stabilizers[start_type][start]['data']:
				if data in code.Stabilizers[mid_type][mid]['data']:
					c = code.Primal.node[data]['charge'][ct]
					count = code.Stabilizers[start_type][start]['data'][data]
					sign = code.Sign(count)
					code.Primal.node[data]['charge'][ct] = (c - sign * k)%d
			for data in code.Stabilizers[end_type][end]['data']:
				if data in code.Stabilizers[mid_type][mid]['data']:
					c = code.Primal.node[data]['charge'][ct]
					code.Primal.node[data]['charge'][ct] = (c - sign * k)%d
			continue


		dual2 = nx.shortest_path(code.Dual[t2], end, start)
		triangle1 = dual1[2*i:(2*i+2)] + dual2[1:]
		triangle2 =  dual2[:2] + dual1[2*i+1:(2*i+3)]
		loop1 = path.Path(triangle1)
		if any(loop1.contains_points([data]) == [True] for data in code.Primal.nodes()):
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

	end_charge = (k_end + end_sign * k)%d
	code.Stabilizers[t][e[0]]['charge'][ct] = end_charge

	if e[0] in uc.nodes():
		for m in cc:
			if m[0] == e[0]:
				cc.remove(m)
		if end_charge == 0:
			uc.remove_node(e[0])
		else:
			uc.node[e[0]]['charge'] = end_charge
			cc.append((e[0],{'charge':end_charge,'type':t}))
	else:
		uc.add_node(e[0], charge = end_charge, type = t)
		cc.append((e[0],{'charge':end_charge,'type':t}))	

	return cc, uc, code


def GCC_Two_Color_Simplify(cc, uc, code, ct, cntr):
	if any(cc[t] == [] for t in cc):

		return cc, uc, code

	d = code.dimension
	ms = {}
	for t in code.types:
		ms[t] = cc[t][0]
	triangle, uc, code, ct = GCC_Connect(ms, cc, uc, code, ct)
	uc, code = GCC_Two_Color_Transport(triangle, uc, code, ct)
	return cc, uc, code


def GCC_Connect(ms, cc, uc, code, ct):

	t1 = 'red'
	for t in ms:
		if ms[t][0] in code.External[t]:
			t1 = t
			break
	d = code.dimension

	[t2, t3] = code.complementaryTypes(t1)
	m1, m2, m3 = ms[t1], ms[t2], ms[t3]
	m1_data = code.Stabilizers[t1][m1[0]]['data']

	if m1[0] in code.Dual[t3].neighbors(m2[0]):
		m2_new = m2

	else:
		for m in code.Stabilizers[t2]:
			if m1[0] in code.Dual[t3].neighbors(m):
				break
		m2_new = (m, {'charge':0, 'type':t2})
		cc[t2], uc, code = GCC_One_Color_Transport(m2, m2_new, cc[t2], uc, code, t2, ct)
	
	m2_new = cc[t2][0]
	m2_data = code.Stabilizers[t2][m2_new[0]]['data']
	if m1[0] in code.Dual[t2].neighbors(m3[0]) and m2[0] in code.Dual[t1].neighbors(m3[0]):
		m3_new = m3

	else:
		for m in code.Stabilizers[t3]:
			if m not in code.External[t3]:
				if any(node in code.Stabilizers[t3][m]['data'] for node in m1_data):
					if any(node in code.Stabilizers[t3][m]['data'] for node in m2_data):
						break
		m3_new = (m, {'charge':0, 'type':t3})

		cc[t3], uc, code = GCC_One_Color_Transport(m3, m3_new, cc[t3], uc, code, t3, ct)

	m3_new = cc[t3][0]
	ms[t2], ms[t3] = m2_new, m3_new

	return ms, uc, code, ct


def GCC_Two_Color_Transport(triangle, uc, code, ct):
	for t1 in triangle:
		if triangle[t1][1]['charge'] != 0:
			k = triangle[t1][1]['charge']
			m0 = triangle[t1][0]
			break

	r, g, b = triangle['red'][0], triangle['green'][0], triangle['blue'][0]
	d = code.dimension
	sides = code.types[t1]['sides']

	cycle = [r, g, b, r]
	for data in code.Stabilizers['red'][r]['data']:
		if data in code.Stabilizers['blue'][b]['data'] and data in code.Stabilizers['green'][g]['data']:
			c = code.Primal.node[data]['charge'][ct]
			count = code.Stabilizers[t1][triangle[t1][0]]['data'][data]
			sign = code.Sign(count, sides)
			code.Primal.node[data]['charge'][ct] = (c - sign * k)%d
	
	uc.remove_node(m0)
	code.Stabilizers[t1][m0]['charge'][ct] = 0

	for color in code.complementaryTypes(t1):
		m = triangle[color][0]
		c = triangle[color][1]['charge']
		charge = (c - k)%d

		code.Stabilizers[color][m]['charge'][ct] = charge
		if charge == 0:
			if m not in uc.nodes():
				visualization.PlotPlaquette(code, "Logical Error", 2)
				plt.show()
			if m in uc.nodes():
				uc.remove_node(m)
		else:
			uc.node[m]['charge'] = charge

	return uc, code

def GCC_Boundary_Simplify(cc, uc, code, ct, cntr, scale):
	ints = cc['red'] + cc['blue'] + cc['green']

	if len(ints) == 2:
		# print "TWO COLS"
		cc, uc, code = GCC_Boundary_Two_Color_Simplify(ints, cc, uc, code, ct, cntr, scale)
	elif len(ints) == 1:
		t, m = ints[0][1]['type'], ints[0][0]
		# print "ONE YALL"
		cc, uc, code = GCC_Boundary_One_Color_Simplify(m, cc, uc, code, t, ct, cntr, scale)
	
	for node in uc.nodes():
		if any(node in code.External[t] for t in code.External):
			uc.remove_node(node)
	for t in code.External:
		for ext in code.External[t]:
			code.Stabilizers[t][ext]['charge'][ct] = 0
	for t in cc:
		for ext in cc[t]:
			if ext[0] in code.External[t]:
				cc[t].remove(ext)
	return cc, uc, code
		
def GCC_Boundary_One_Color_Simplify(m, cc, uc, code, t, ct, cntr, scale):
	[t1,t2] = code.complementaryTypes(t)
	c = cc[t][0][1]['charge']
	d = code.dimension

	if any(common.euclidean_dist(ext, m) < scale for ext in code.External[t]):
		for ext in code.External[t]:
			if common.euclidean_dist(ext, m) < scale:
				uc.add_node(ext, charge = d-c, type = t)
				cc[t].append((ext,{'charge':d-c,'type':t}))
				code.Stabilizers[t][ext]['charge'][ct] = d-c
				cc[t], uc, code = GCC_One_Color_Simplify(cc[t], uc, code, t, ct, cntr, scale)		
				# print "CONNECTION", m, ext
				break

	elif any(common.euclidean_dist(ext, m) < scale for ext in code.External[t1]) and any(common.euclidean_dist(ext, m) < scale for ext in code.External[t2]):
		if any(ext1 in code.External[t1] for ext1 in code.Dual[t2].neighbors(m)) and any(ext2 in code.External[t2] for ext2 in code.Dual[t1].neighbors(m)):
			m_new = m
		else:
			for m_new in code.Stabilizers[t]:
				if any(ext1 in code.External[t1] for ext1 in code.Dual[t2].neighbors(m_new)) and any(ext2 in code.External[t2] for ext2 in code.Dual[t1].neighbors(m_new)):
					break
			s, e = cc[t][0], (m_new,{'charge':0,'type':t})
			uc.add_node(m_new, charge = 0, type = t)
			cc[t].append(e)
			cc[t], uc, code = GCC_One_Color_Transport(s, e, cc[t], uc, code, t, ct)

		for ext1 in code.External[t1]:
			if ext1 in code.Dual[t2].neighbors(m_new):
				break
		for ext2 in code.External[t2]:
			if ext2 in code.Dual[t1].neighbors(m_new):
				break
		# print "CONNECTION", m, ext1, ext2
		uc.add_node(ext1, charge = c, type = t1)
		cc[t1].append((ext1,{'charge':c,'type':t1}))
		uc.add_node(ext2, charge = c, type = t2)
		cc[t2].append((ext2,{'charge':c,'type':t2}))

		cc, uc, code = GCC_Two_Color_Simplify(cc, uc, code, ct, cntr)
	return cc, uc, code 

def GCC_Boundary_Two_Color_Simplify(ints, cc, uc, code, ct, cntr, scale):
	# print "HOWDY ROWDY"
	# print ints
	[m0, m1] = ints
	c0, c1 = m0[1]['charge'], m1[1]['charge']
	t0, t1 = m0[1]['type'], m1[1]['type']
	if c0 == c1:
		# print "Equality!!!"
		t2 = code.complementaryType([t0,t1])
		if any((common.euclidean_dist(ext, m0[0]) < scale and common.euclidean_dist(ext, m1[0]) < scale) for ext in code.External[t2]):
			# print "External"
			if any( (ext in code.Dual[t1].neighbors(m0[0])  and ext in code.Dual[t0].neighbors(m1[0])) for ext in code.External[t2]):
				for ext in code.External[t2]:
					if ext in code.Dual[t1].neighbors(m0[0]) and ext in code.Dual[t0].neighbors(m1[0]):
						# print "CONNECTION", m0, m1, ext
						break
			else:
				for ext in code.External[t2]:
					if common.euclidean_dist(ext, m0[0]) < scale and common.euclidean_dist(ext, m1[0]) < scale:
						# print "SECOND RATE", m0, m1, ext
						break

			uc.add_node(ext, charge = c0, type = t2)
			cc[t2].append((ext,{'charge':c0,'type':t2}))
			code.Stabilizers[t2][ext]['charge'][ct] = c0
			cc, uc, code = GCC_Two_Color_Simplify(cc, uc, code, ct, cntr)
		elif any(common.euclidean_dist(ext, m0[0]) < scale for ext in code.External[t0]) and any(common.euclidean_dist(ext, m1[0]) < scale for ext in code.External[t1]):
			# print "YEEEEHHAAAAAA"
			cc, uc, code = GCC_Boundary_One_Color_Simplify(m0[0], cc, uc, code, t0, ct, cntr, scale)
			cc, uc, code = GCC_Boundary_One_Color_Simplify(m1[0], cc, uc, code, t1, ct, cntr, scale)

	elif any(common.euclidean_dist(ext, m0[0]) < scale for ext in code.External[t0]) and any(common.euclidean_dist(ext, m1[0]) < scale for ext in code.External[t1]):
		# print "YEEEEHHAAAAAA"
		cc, uc, code = GCC_Boundary_One_Color_Simplify(m0[0], cc, uc, code, t0, ct, cntr, scale)
		cc, uc, code = GCC_Boundary_One_Color_Simplify(m1[0], cc, uc, code, t1, ct, cntr, scale)
	
	return cc, uc, code







