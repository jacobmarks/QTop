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


# Example usage as follows:
# clustering = MinWeightMatch()
# decoder = MWPM_Decoder(clustering)
# Decode(code, syndrome, decoder)

# or for color codes, you could have
# clustering = HDRG_cluster()
# decoder = DSP(clustering)
# Decode(code, syndrome, decoder)




########### Syndrome generation and manipulation ############

def syndrome(code, type, charge_type = 'Z'):
    syndrome = nx.Graph()

    # Find all non-trivial check operators
    for measure in code.stabilizers[type]:
    	center = code.stabilizers[type][measure].center
    	if center.charge[charge_type] != 0:
    		syndrome.add_node(center)

    return syndrome

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




############ Decode function and base decoder classes ############


# make sure we are using the correct syndrome for each decoding

## Call with Decode(code, syndrome, decoder)
def Decode(code, syndrome, decoder):
	decoder(code, syndrome)


class decoder(object):

	HDRG_cluster = HDRG_cluster()

	def __init__(self, cluster = HDRG_cluster):
		self.cluster = cluster

	def __call__(self, code, syndrome):
		clusters = {}
		for charge_type in ['X','Z']:
			clusters[charge_type] = {}
			for type in code.types:
				clusters[charge_type][type] = []
				clusters[charge_type][type] = self.cluster(code, syndrome, type, charge_type)
		
		code = self.recover(code, clusters)
		return code


class surface_decoder(decoder):

	def __init__(self, cluster):
		self.cluster = cluster
		self.recover = surface_recovery()

	def __call__(self, code, syndrome):
		decoder.__call__(self, code, syndrome)

class MWPM_Decoder(surface_decoder):

	def __call__(self, code, syndrome):
		surface_decoder.__call__(self, code, syndrome)

	def cluster(self, code, syndrome, type, charge_type):
		MinWeightMatch(code, syndrome, type, charge_type)

class HDRG_Decoder(surface_decoder):

	def __call__(self, code, syndrome):
		surface_decoder.__call__(self, code, syndrome)

	def cluster(self, code, syndrome, type, charge_type):
		HDRG_cluster(code, syndrome, type, charge_type)


# Delfosse Surface Projection decoder for color codes #

class DSP(decoder):

	def __init__(self, cluster):
		self.cluster = cluster
		self.recover = fill_recovery()

	def __call__(self, code, syndrome):
		decoder.__call__(self, code, syndrome)


	

	



############### Base clustering classes ################
# Class of algorithms to break syndrome up into local clusters
# of non-trivial stabilizers

class cluster(object):

	def __init__(self):
		pass


# Min Weight Perfect Matching
# Pairs up nodes to minimize overall weight
# using Edmunds' Blossom algorithm

class MinWeightMatch(cluster):

	def __init__(self):
		pass

	def __call__(self, code, syndrome, type, charge_type):
		for check1 in syndrome.nodes():
			for check2 in syndrome.nodes():
				if check1 != check2:
					weight = - distance(check1, check2, type)
					syndrome.add_edge(*(check1, check2), weight=weight)

		if code.geometry != 'toric':
			syndrome = addExternalNodes(code, syndrome, type)

		temp_matching = nx.max_weight_matching(syndrome, maxcardinality=True)
		matching = []

		for node, neighbor in temp_matching:
			if [neighbor,node] not in matching:
				if node in self.dual.nodes() or neighbor in self.dual.nodes():
					matching.append([node, neighbor])

	return matching



 # Hard Decision Renormalization Group clustering 
 # iteratively increase distance scale, annihilating
 # neutral clusters at each level

class HDRG_cluster(cluster):

	def __init__(self):
		pass

	def __call__(self, code, syndrome, type, charge_type):
 
    	# Add in external elements
    	if code.geometry != 'toric':
    		for node in code.external[type]:
    			syndrome.add_node(node)
 
 
    	scale = 1
    	NeutralClusters = []
    	dim = code.dimension
    	# Copy Graph for modification
    	unclustered_graph = syndrome.copy()
 
    	while StillClustering(unclustered_graph):
        	clusters = partition(code, unclustered_graph, scale, type, charge_type)
        	for cluster in clusters:
        		if cluster_charge
            	NetCharge = cluster_charge(unclustered_graph, cluster, dim) 
            	if NetCharge == 'Neutral':
                	NeutralClusters.append(cluster)
                	fuse(code, unclustered_graph, cluster, dim, type, charge_type)
                	annihilate(unclustered_graph, cluster)
            	elif NetCharge == 'BoundaryNeutral':
                	NeutralClusters.append(cluster)
                	bound_fuse(code, unclustered_graph, cluster, dim, type, charge_type)
                	annihilate(unclustered_graph, cluster)
         
        	scale += 1      # increase distance scale


############### Base recovery classes ################
# Class of algorithms to apply recovery chains

class recovery(object):

	def __init__(self):
		pass


class surface_recovery(recovery):

	def __init__(self):
		pass

	def __call__(self, code, clusters, type, charge_type):
		for charge_type in clusters:
			for type in clusters[charge_type]:
				for cluster in clusters[charge_type][type]:
					fuse(code, cluster, type, charge_type)


# Make closed loops on dual lattice and fill in the interior
class fill_recovery(recovery):

	def __init__(self):
		pass

############# Functions used in clustering and recovery ############

def annihilate(unclustered_graph, cluster):
    for node in cluster:
        unclustered_graph.remove_node(node)
 
 
# Boolean function, True if still have nodes to cluster
def StillClustering(code, unclustered_graph):
    for node in unclustered_graph.nodes():
    	for type in code.types:
        	if node not in code.external[type]:
            	return True
    return False


# Choose fixed node for cluster fusion
def central_node(code, unclustered_graph, cluster):
	internal = False
    # Initialize Center
    for node in cluster:
    	for type in code.types:
        	if node not in code.external[type]:
            	heart = node
        		max_degree = unclustered_graph.degree(heart)
        		internal = True
            	break
        if internal:
        	break
 
    for node in cluster:
        if unclustered_graph.degree(node) > max_degree:
        	for type in code.types:
        		if node not in code.external[type]:
            		heart = node
            		max_degree = unclustered_graph.degree(heart)
 
    return heart



# Cluster is either
# A) Neutral: no boundary elements, and sum of charges == 0 mod d
# B) Boundary-Neutral: has a boundary element
# C) Charged: no boundary elements, and sum of charges != 0 mod d
def cluster_charge(code, unclustered_graph, cluster, charge_type):
	dim = code.dimension

    NetCharge = 0
    BOUND = False
    Charged = True
 
    for node in cluster:
        if node not in code.dual.nodes():
            BOUND = True
            for mate in cluster:
                if node in code.dual.nodes():
                    Charged = False
 
        NetCharge += node.charge[charge_type]
     
    if BOUND:
        if Charged:
            return 'Charged'
        elif NetCharge % d == 0:
            # Neutral supercedes BoundaryNeutral
            removed_nodes = []
            for node in cluster:
                 
                if node not in code.dual.nodes():
                    removed_nodes.append(node)
            for node in removed_nodes:
                cluster.remove(node)
            return 'Neutral'
        else:
            return 'BoundaryNeutral'
 
    # No Boundary Elements
    else:
        if NetCharge % d == 0:
            return 'Neutral'
        else:
            return 'Charged'


# partitions nodes into maximally disjoint clusters
# for a given distance 
def partition(code, unclustered_graph, scale, type):
    # Make edges on Unclustered graph
    # between nodes separated by distance 'scale'
    for node1 in unclustered_graph.nodes():
        for node2 in unclustered_graph.nodes():
            if node1 != node2:
                distance = code.distance(node1, node2, lattice_type = type)
                if distance <= scale:
                    unclustered_graph.add_edge(*(node1, node2), weight=distance)
    Clusters = []
    # Networkx connected components analysis
    subgraphs = nx.connected_component_subgraphs(unclustered_graph)
    for i, sg in enumerate(subgraphs):
        Clusters.append(sg.nodes(data=True))

    return Clusters



def recovery_chain(code, terminal1, terminal2, type, charge_type):
	if terminal1 not in code.shrunk[type]:
		internal1 = code.external[type][terminal1]['measure'][0]
		EXT1 = True
	else:
		internal1 = terminal1


	if terminal2 not in code.shrunk[type]:
		internal2 = code.external[type][terminal2]['measure'][0]
		EXT2 = True
	else:
		internal2 = terminal2


	chain = nx.shortest_path(code.shrunk[type],internal1, internal2)
	if EXT1:
		chain = [terminal1] + chain
	if EXT2:
		chain = chain + [terminal2]


def fuse(code, cluster, type, charge_type):
    # use syndrome transport to move excitations around lattice
    fusion_graph = code.shrunk[type].subgraph(cluster)
    center = central_node(fusion_graph, cluster)
    for mate in cluster:
        if mate != center:
            transport(code, center, mate, type, charge_type)










def bound_fuse(code, cluster, type, charge_type):
    # Find Excess charge and give all charge to one boundary element
    # then delete the other boundary elements
 	dim = code.dimension
    NetCharge = 0
    for node in cluster:
        NetCharge += node.charge[charge_type]
    NetCharge = NetCharge % dim
 
    IsNeutral = False
    removed_boundaries = []
 
    for node in cluster:
       if node not in code.shrunk[type]:
            if not IsNeutral:
                charged_boundary = node
                charged_boundary.charge[charge_type] = (- NetCharge)% dim
                IsNeutral = True
            else:
                removed_boundaries.append(node)    
             
    for boundary in removed_boundaries:
        cluster.remove(boundary)
 
 
    center = central_node(code, unclustered_graph, cluster)
    for mate in cluster:
        if mate != center:
            transport(code, center, mate, type, charge_type)

def transport(code, fixed_node, mobile_node, type, charge_type):
	dim = code.dimension
 	chain = recovery_chain(code, mobile_node, fixed_node, type, charge_type)
 	chain_length = len(chain) - 1

    first_link = chain[0]
    if first_link not in code.dual.nodes():
    	charge = first_link.charge[charge_type]
    	for shared in code.external[type][first_link]['data']:
    		if shared in code.stabilizers[type][chain[1].position].data:

 
    for link in range(chain_length):
        previous, next = chain[link], chain[link + 1]
        for shared in code.stabilizers[type][previous.position].data:
        	if shared in code.stabilizers[type][next.position].data:
        		for count in code.stabilizers[type][previous.position].order:
        			if data = code.stabilizers[type][previous.position].order[data]:
        				sign = pow(-1, order%2)

    	delta_charge = sign * charge
    	data_charge = shared.charge[charge_type]
    	code.primal.node[shared].charge[charge_type] = (data_charge - delta_charge)%dim
 

def bound_transport(Tiling, Dual, fixed_node, mobile_node, d, type):
    if 'boundary' in mobile_node[1]:
        charge = mobile_node[1]['charge']
        first_link = Tiling['Boundary'][type][mobile_node[0]]['measure']
        shared = Tiling['Boundary'][type][mobile_node[0]]['data']
        Sign = Tiling['Boundary'][type][mobile_node[0]]['sign']
        delta_charge = Sign * charge
        first_link_charge = Tiling['Primal'].node[shared]['charge'][type]
        Tiling['Primal'].node[shared]['charge'][type] = (first_link_charge - delta_charge)%d
         
        mobile_node = first_link
        chain = nx.shortest_path(Dual[type], mobile_node, fixed_node[0])
        chain_length = len(chain) - 1
        for link in range(chain_length):
            previous, next = chain[link], chain[link + 1]
            for shared in Tiling['Faces'][type][previous]['data']:
                if shared in Tiling['Faces'][type][next]['data']:
                    Sign = Tiling['Faces'][type][previous]['data'][shared]
                    delta_charge = Sign * charge
                    data_charge = Tiling['Primal'].node[shared]['charge'][type]
                    Tiling['Primal'].node[shared]['charge'][type] = (data_charge - delta_charge)%d
 
    else:
        transport(Tiling, Dual, fixed_node, mobile_node, d, type)



def bulk_recovery(code, terminal1, terminal2, type, charge_type):

    for link in range(chain_length):
        vertex1 = measure_chain[link]
        vertex2 = measure_chain[link + 1]
        for data_qubit in PlaquetteLattice.neighbors(vertex1):
            if data_qubit in PlaquetteLattice.neighbors(vertex2):
                prior_charge = PlaquetteLattice.node[data_qubit][charge][charge_type]
                PlaquetteLattice.node[data_qubit][charge][charge_type] = (prior_charge + 1) % 2


# return the correction chains on dual lattice corresponding
# to the specified shrunk lattice
def SurfaceProjection(code, syndrome, shrunk_type, type):
	shrunk_lattice = code.shrunk[shrunk_type]
	shrunk_syndrome = syndrome.copy()
	for stabilizer in syndrome:
		if syndrome.node[stabilizer]['type'] == shrunk_type:
			shrunk_syndrome.remove(stabilizer)

	return shrunk_syndrome

	# now we have a surface code and its corresponding syndrome 

 

 

 
 
 

 








