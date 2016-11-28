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


class HDRG_decoder(surface_decoder):

    def __call__(self, code):
        return surface_decoder.__call__(self, code)

    def algorithm(self):
        return HDRG()


class HDRG(matching_algorithm):

    def __init__(self):
        pass

    def __call__(self, code, Syndrome, type, charge_type):
        return Renormalization(code, Syndrome, type, charge_type)


def Renormalization(code, Syndrome, type, charge_type):

    dim = code.dimension

    # Add in Boundary elements
    for node in code.External[type]:
        Syndrome.add_node(node, charge=0, external=True)


    scale = 1
    NeutralClusters = []
    # Copy Graph for modification
    UnclusteredGraph = Syndrome.copy()

    while StillClustering(UnclusteredGraph):
        clusters = Partition(UnclusteredGraph, code, type, scale)
        for cluster in clusters:
            NetCharge = ClusterCharge(UnclusteredGraph, cluster, dim) 
            if NetCharge == 'Neutral':
                NeutralClusters.append(cluster)
                Fuse(code, UnclusteredGraph, cluster, dim, type, charge_type)
                Annihilate(UnclusteredGraph, cluster)
            elif NetCharge == 'BoundaryNeutral':
                NeutralClusters.append(cluster)
                BoundFuse(code, UnclusteredGraph, cluster, dim, type, charge_type)
                Annihilate(UnclusteredGraph, cluster)
        
        scale += 1      # increase distance scale

    return code



def Annihilate(UnclusteredGraph,cluster):
    for node in cluster:
        UnclusteredGraph.remove_node(node[0])


# Boolean function, True if still have nodes to cluster
def StillClustering(UnclusteredGraph):
    for node in UnclusteredGraph.nodes():
        if 'external' not in UnclusteredGraph.node[node]:
            return True
    return False


# Cluster is either
# A) Neutral
# B) Boundary-Neutral
# C) Charged
def ClusterCharge(UnclusteredGraph, cluster, dim):
    NetCharge = 0
    BOUND = False
    Charged = True

    for node in cluster:
        if 'external' in node[1]:
            BOUND = True
            for mate in cluster:
                if 'external' not in mate[1]:
                    Charged = False

        NetCharge += node[1]['charge']
    
    if BOUND:
        if Charged:
            return 'Charged'
        elif NetCharge % dim == 0:
            # Neutral supercedes BoundaryNeutral
            removed_nodes = []
            for node in cluster:
                
                if 'external' in node[1]:
                    removed_nodes.append(node)
            for node in removed_nodes:
                cluster.remove(node)
            return 'Neutral'
        else:
            return 'BoundaryNeutral'

    # No Boundary Elements
    else:
        if NetCharge % dim == 0:
            return 'Neutral'
        else:
            return 'Charged'
# partitions nodes into maximally disjoint clusters
# for a given distance


def Partition(UnclusteredGraph, code, type, scale):
    # Make edges on Unclustered graph
    # between all nodes separated by distance 'scale'
    # on Dual Lattice
    for node1 in UnclusteredGraph.nodes():
        for node2 in UnclusteredGraph.nodes():
            if node1 != node2:
                d = code.distance(type, node1, node2)
                if d <= scale:
                    UnclusteredGraph.add_edge(*(node1, node2), weight=d)
    Clusters = []
    # Networkx connected components analysis
    subgraphs = nx.connected_component_subgraphs(UnclusteredGraph)
    for i, sg in enumerate(subgraphs):
        Clusters.append(sg.nodes(data=True))
            

    return Clusters


# Choose fixed node for cluster fusion
def CentralNode(UnclusteredGraph, cluster):
    # Initialize Center
    for node in cluster:
        if 'external' not in node[1]:
            center = node
            degree = UnclusteredGraph.degree(center[0])
            break

    for node in cluster:
        if UnclusteredGraph.degree(node[0]) > degree and 'external' not in node[1]:
            center = node
            degree = UnclusteredGraph.degree(center[0])

    return center


def Fuse(code, UnclusteredGraph, cluster, dim, type, charge_type):
    # use syndrome transport to move excitations around lattice
    center = CentralNode(UnclusteredGraph, cluster)
    for mate in cluster:
        if mate != center:
            Transport(code, center, mate, dim, type, charge_type)

def Transport(code, fixed_node, mobile_node, dim, type, charge_type):

    chain = nx.shortest_path(code.Dual[type], mobile_node[0], fixed_node[0])
    chain_length = len(chain) - 1
    first_link = chain[0]

    charge = code.Stabilizers[type][first_link]['charge'][charge_type]

    for link in range(chain_length):
        previous, next = chain[link], chain[link + 1]
        for shared in code.Stabilizers[type][previous]['data']:
            if shared in code.Stabilizers[type][next]['data']:
                num_sides = code.Stabilizers[type][next]['sides']
                count = code.Stabilizers[type][previous]['data'][shared]
                sign = common.Sign(count, num_sides)
                delta_charge = sign * charge
                data_charge = code.Primal.node[shared]['charge'][charge_type]
                code.Primal.node[shared]['charge'][charge_type] = (data_charge - delta_charge)%dim



def BoundFuse(code, UnclusteredGraph, cluster, dim, type, charge_type):
    # Find Excess charge and give all charge to one external element
    # then delete the other external elements

    NetCharge = 0
    for node in cluster:
        NetCharge += node[1]['charge']
    NetCharge = NetCharge % dim 

    IsNeutral = False
    removed_externals = []

    for node in cluster:
        if 'external' in node[1]:
            if not IsNeutral:
                charged_external = node
                charged_external[1]['charge'] = - NetCharge
                IsNeutral = True
            else:
                removed_externals.append(node)    
            
    for external in removed_externals:
        cluster.remove(external)


    center = CentralNode(UnclusteredGraph, cluster)
    for mate in cluster:
        if mate != center:
            BoundTransport(code, center, mate, dim, type, charge_type)

def BoundTransport(code, fixed_node, mobile_node, dim, type, charge_type):
    if 'external' in mobile_node[1]:
        charge = mobile_node[1]['charge']
        first_link = code.External[type][mobile_node[0]]['measure']
        shared = code.External[type][mobile_node[0]]['data']

        count = code.External[type][mobile_node[0]]['order']
        num_sides = code.External[type][mobile_node[0]]['sides']
        sign = common.Sign(count, num_sides)
        delta_charge = sign * charge
        first_link_charge = code.Primal.node[shared]['charge'][charge_type]
        code.Primal.node[shared]['charge'][charge_type] = (first_link_charge - delta_charge)%dim
        
        mobile_node = first_link
        chain = nx.shortest_path(code.Dual[type], mobile_node, fixed_node[0])
        chain_length = len(chain) - 1
        for link in range(chain_length):
            previous, next = chain[link], chain[link + 1]
            for shared in code.Stabilizers[type][previous]['data']:
                if shared in code.Stabilizers[type][next]['data']:

                    num_sides = code.Stabilizers[type][previous]['sides']
                    count = code.Stabilizers[type][previous]['data'][shared]
                    sign = common.Sign(count, num_sides)
                    delta_charge = sign * charge
                    data_charge = code.Primal.node[shared]['charge'][charge_type]
                    code.Primal.node[shared]['charge'][charge_type] = (data_charge - delta_charge)%dim

    else:
        Transport(code, fixed_node, mobile_node, dim, type, charge_type)
