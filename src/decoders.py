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
from matplotlib import path

############ Decode function and base decoder classes ############

## Call with Decode(code, decoder)
def Decode(code, decoder):
    decoder(code)


class decoder(object):

    def __call__(self, code):
        return code

class surface_decoder(decoder):

    def __init__(self):
        self.recover = self.algorithm()

    def __call__(self, code):
        for type in code.types:
            for charge_type in ['X', 'Z']:
                syndrome = code.Syndrome(type, charge_type)
                code = self.recover(code, syndrome, type, charge_type)

        code = reset_measures(code)
        return code

class MWPM_decoder(surface_decoder):

    def __call__(self, code):
        return surface_decoder.__call__(self, code)

    def algorithm(self):
        return MWPM()

class HDRG_decoder(surface_decoder):

    def __call__(self, code):
        return surface_decoder.__call__(self, code)

    def algorithm(self):
        return HDRG()

def reset_measures(code):
    for type in code.types:
        for measure_qubit in code.Stabilizers[type]:
            code.Stabilizers[type][measure_qubit]['charge'] = Charge()
    return code



class matching_algorithm(object):

    def __init__(self):
        pass

class MWPM(matching_algorithm):

    def __init__(self):
        pass

    def __call__(self, code, Syndrome, type, charge_type):
        matching = MinWeightMatching(code, Syndrome, type, charge_type)
        code = Recovery(code, matching, type, charge_type)
        return code

class HDRG(matching_algorithm):

    def __init__(self):
        pass

    def __call__(self, code, Syndrome, type, charge_type):
        return Renormalization(code, Syndrome, type, charge_type)

 ############# Minimum Weight Perfect Matching ###############

# Given a network of weighted edges, finds the 
# minimum weight perfect matching of the nodes.

def MinWeightMatching(code, Syndrome, type, charge_type):
    dim = code.dimension

    # Fully connect check operators
    for check1 in Syndrome.nodes():
        for check2 in Syndrome.nodes():
            if check1 != check2:
                # weight = - code.distance(type, check1, check2)
                weight = - euclidean_dist(check1, check2)
                Syndrome.add_edge(*(check1, check2), weight=weight)

    # Generate Boundary Graph
    External_Graph = nx.Graph()

    for node in Syndrome.nodes():
        charge = Syndrome.node[node]['charge']
        external_node = AssociatedExternal(node, code.Dual[type], code.External[type])
        External_Graph.add_node(external_node, charge=(-charge) % dim)
        # weight = -code.distance(type, node, external_node)
        weight = - euclidean_dist(node, external_node)
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



###################### Renormalization Decoder ######################


# Takes as input a fully connected graph
# containing all non-trivial syndrome
# measurement results from code cycle

# Iteratively increases distance
# finding maximal disjoint clusters and
# annihilating neutral clusters at each distance

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
                sign = Sign(count, num_sides)
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
        sign = Sign(count, num_sides)
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
                    sign = Sign(count, num_sides)
                    delta_charge = sign * charge
                    data_charge = code.Primal.node[shared]['charge'][charge_type]
                    code.Primal.node[shared]['charge'][charge_type] = (data_charge - delta_charge)%dim

    else:
        Transport(code, fixed_node, mobile_node, dim, type, charge_type)




def DSP_Matching(Syndrome, External, dim):

    # Fully connect check operators
    for check1 in Syndrome.nodes():
        for check2 in Syndrome.nodes():
            if check1 != check2:
                weight = - euclidean_dist(check1, check2)
                Syndrome.add_edge(*(check1, check2), weight=weight)

    # Generate Boundary Graph
    External_Graph = nx.Graph()

    for node in Syndrome.nodes():
        charge = Syndrome.node[node]['charge']
        external_node = DSP_AssociatedExternal(node, External)
        External_Graph.add_node(external_node, charge=(-charge) % dim)
        weight = - euclidean_dist(node, external_node)
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
    min_dist = euclidean_dist(int_node, associate)

    for candidate in external_nodes:
        distance = euclidean_dist(int_node, candidate)
        if distance < min_dist:
            min_dist = distance
            associate = candidate
    return associate

def DSP_AssociatedInternal(ext_node, internal_nodes):
    associate = internal_nodes[0]
    min_dist = euclidean_dist(ext_node, associate)

    for candidate in internal_nodes:
        distance = euclidean_dist(ext_node, candidate)
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


