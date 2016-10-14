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
from error_models import *
from math import *
import networkx as nx
import sys



# common base class for decoding functions
# as well as decoding based on
# Edmunds' Blossom algorithm and
# Hard Decision Renormalization Group Methods



############ Decode function and base decoder classes ############


# make sure we are using the correct syndrome for each decoding

## Call with Decode(code, syndrome, decoder)
def Decode(code, syndrome, decoder):
    decoder(code, syndrome)


class decoder(object):

    def __init__(self, match):
        self.match = match()

    def __call__(self, code, syndrome):
        pairs = self.match(code, syndrome)
        code = self.recover(code, pairs)
        return code


class surface_decoder(decoder):

    def __init__(self, match):
        self.match = self.match()
        self.recover = surface_recovery()

    def __call__(self, code):
        for type in code.types:
            for charge_type in ['X', 'Z']:
                syndrome = self.syndrome(code, type, charge_type)
                pairs = self.match(code, syndrome, type, charge_type)
        
                code = self.recover(code, pairs, type, charge_type)

        code = reset_measures(code)
        return code

    def syndrome(self, code, type, charge_type):
        syndrome = nx.Graph()
        # Find all non-trivial check operators
        for measure_position in code.stabilizers[type]:
            measure_qubit = code.syndrome[measure_position]
            charge = measure_qubit.charge[charge_type]
            if charge != 0:
                syndrome.add_node(measure_position, charge = charge)

        return syndrome

class MWPM_Decoder(surface_decoder):

    def __init__(self):
        self.match = self.matching()
        self.recover = surface_recovery()

    def __call__(self, code):
        return surface_decoder.__call__(self, code)

    def matching(self):
        return MinWeightMatch()

class HDRG_Decoder(surface_decoder):

    def __init__(self):
        self.match = self.matching()
        self.recover = surface_recovery()

    def __call__(self, code):
        return surface_decoder.__call__(self, code)


    def matching(self):
        return HDRG_Match()


# Delfosse Surface Projection decoder for color codes #

class DSP(decoder):

    def __init__(self, match):
        self.match = match
        self.recover = fill_recovery()

    def __call__(self, code):
        for charge_type in ['X', 'Z']:
            for shrunk_type in code.types:
                syndrome = self.syndrome(code, shrunk_type, charge_type)

            pairs[shrunk_type] = self.match(code, syndrome, shrunk_type, charge_type)

            code = self.recover(code, pairs, charge_type)

        code = reset_measures(code)
        return code

    def syndrome(self, code, shrunk_type, charge_type):
        syndrome = nx.Graph()
        # Find all non-trivial check operators
        for measure_position in code.syndrome:
            measure_qubit = code.syndrome[measure_position]
            if measure_qubit.type != shrunk_type:
                charge = measure_qubit.charge[charge_type]
                if charge != 0:
                    syndrome.add_node(measure_position, charge = charge)

        return syndrome



############### Base clustering classes ################
# Class of algorithms to break syndrome up into local clusters
# of non-trivial stabilizers

class Match(object):

    def __init__(self):
        pass


# Min Weight Perfect Matching
# Pairs up nodes to minimize overall weight
# using Edmunds' Blossom algorithm

class MinWeightMatch(Match):

    def __init__(self):
        pass

    def __call__(self, code, syndrome, type, charge_type):
        matches = []

        for check1 in syndrome.nodes():
            for check2 in syndrome.nodes():
                if check1 != check2:
                    weight = - code.distance(check1, check2, type)
                    syndrome.add_edge(*(check1, check2), weight=weight)

        temp_matching = nx.max_weight_matching(syndrome, maxcardinality=True)
        for node in temp_matching:
            neighbor = temp_matching[node]
            if [neighbor,node] not in matches:
                if node in code.dual.nodes() or neighbor in code.dual.nodes():
                    matches.append([node, neighbor])

        return matches


 # Hard Decision Renormalization Group matching
 # iteratively increase distance scale, annihilating
 # neutral clusters at each level

class HDRG_Match(Match):

    def __init__(self):
        pass

    def __call__(self, code, syndrome, type, charge_type):
        matches = []

        scale = 1
        NeutralClusters = []
        dim = code.dimension
        # Copy Graph for modification
        unclustered_graph = syndrome.copy()

        while StillClustering(code, unclustered_graph):
            clusters = partition(code, unclustered_graph, scale, type, charge_type)
            for cluster in clusters:
                if fuse_cluster(code, cluster, charge_type):
                    matches = PairOffCluster(matches, code, cluster, type, charge_type)
                    annihilate(unclustered_graph, cluster)
            scale += 1
        return matches
        



############### Base recovery classes ################
# Class of algorithms to apply recovery chains

def reset_measures(code):
        perfect_gates = ErrorModel()
        code = code.CodeCycle(perfect_gates, 0)
        return code

class recovery(object):

    def __init__(self):
        pass

    def __call__(self, code, matches, type, charge_type):
        code = self.correct(code, matches, type, charge_type)
        return code

    


class surface_recovery(recovery):

    def __init__(self):
        pass

    def correct(self, code, pairs, type, charge_type):
        for pair in pairs:
            code = fuse(code, pair, type, charge_type)
        return code

    def __call__(self, code, pairs, type, charge_type):
        code = recovery.__call__(self, code, pairs, type, charge_type)
        return code



# Make closed loops on dual lattice and fill in the interior
class fill_recovery(recovery):

    def __init__(self):
        pass

    def __call__(self, code, pairs, charge_type):
        return code

############# Functions used in matching, clustering and recovery ############

def annihilate(unclustered_graph, cluster):
    for node in cluster:
        unclustered_graph.remove_node(node)


# Boolean function, True if still have nodes to cluster
def StillClustering(code, unclustered_graph):
    for node in unclustered_graph.nodes():
        if node in code.syndrome:
            return True
        else:
            return False

def PairOffCluster(matches, code, cluster, type, charge_type):
    
    dim = code.dimension
    unmatched_nodes = cluster

    internal = []


    for node in cluster:
        MATCHED = False
        if node in unmatched_nodes:
            node_charge = code.syndrome[node].charge[charge_type]
            for partner in cluster:
                if partner in unmatched_nodes:
                    partner_charge = code.syndrome[partner].charge[charge_type]
                    if (node_charge + partner_charge)%dim == 0 and node != partner:
                        matches[charge_type].append([node, partner])
                        unmatched_nodes.remove(node)
                        unmatched_nodes.remove(partner)
                        MATCHED = True

    return matches



# Cluster is either
# A) Neutral: no boundary elements, and sum of charges == 0 mod d
# B) Charged: no boundary elements, and sum of charges != 0 mod d
# if neutral, fuse and annihilate cluster
def fuse_cluster(code, cluster, charge_type):
    dim = code.dimension
    NetCharge = 0
    Charged = True


    for node in cluster:
        NetCharge += code.syndrome[type][charge_type][node].charge[charge_type]

    NetCharge = NetCharge%dim

    if NetCharge == 0:
        return True
    else:
        return False

# partitions nodes into maximally disjoint clusters
# for a given distance
def partition(code, unclustered_graph, scale, type, charge_type):
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


def fuse(code, pair, type, charge_type):
    start_node, end_node = pair[1], pair[0]
    recovery_chain = nx.shortest_path(code.shrunk[type],start_node,end_node)
    print recovery_chain
    dim = code.dimension
    chain_length = len(recovery_chain) - 1

    
    charge = code.syndrome[start_node].charge[charge_type]



    for link in range(chain_length):
        first_node, second_node = recovery_chain[link], recovery_chain[(link + 1)%(chain_length+1)]
        for data in code.stabilizers[type][first_node].data:
            if data in code.stabilizers[type][second_node].data:
                if link != 0:
                    if (previous_data, data) not in code.primal.edges():
                        shared_data = data
                else:
                    shared_data = data

        for count in code.stabilizers[type][first_node].order:
            if code.stabilizers[type][first_node].order[count] == shared_data:
                num_sides = code.types[type]['num_sides']
                if count in range(num_sides/2):
                    sign = -1 # Add control to target
                else:
                    sign = 1 # subtract control from target
        previous_data = shared_data
        previous_charge = code.data[shared_data].charge[charge_type]
        code.data[shared_data].charge[charge_type] = (previous_charge + sign * charge)%dim
    return code


def Assessment(code):
    for charge_type in ['X', 'Z']:
        if code.hasLogicalError(charge_type):
            return False
    return True
