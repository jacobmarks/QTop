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

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

#################### Base Classes #####################

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

def manhattan_dist(A,B):
    return abs(A[0]-B[0]) + abs(A[1]-B[1])

def euclidean_dist(A,B):
    return round(pow(abs(A[0]-B[0]),2) + pow(abs(A[1]-B[1]),2),3)

def Charge(X_charge = 0, Z_charge = 0):
    return {'X':X_charge,'Z':Z_charge}

def Sign(count, num_sides):
    if count in range(num_sides/2):
        sign = 1 # Add control to target
    else:
        sign = -1 # subtract control from target
    return sign

class Code:
    'Common base class for all topological codes'

    def __init__(self, depth, dimension = 2):
        self.depth = depth
        self.dimension = dimension
        self.Primal = nx.Graph()
        self.Dual, self.Stabilizers = {}, {}
        self.Boundary, self.External = {}, {}
        self.generateColors()
        self.generateCode()
        self.generateDual()


    ##### Syndrome Generation #####
    # Creates a Syndrome Graph given non-trivial
    # check-operator measurement results

    def Syndrome(self, type, charge_type):
        Syndrome = nx.Graph()
        # Find all non-trivial check operators
        for measure_qubit in self.Stabilizers[type]:
            charge = self.Stabilizers[type][measure_qubit]['charge'][charge_type]
            if charge != 0:
                Syndrome.add_node(measure_qubit, charge=charge)
        return Syndrome

    def PhysicalErrors(self):
        PhysicalError = False
        for type in self.types:
            for measure_qubit in self.Stabilizers[type]:
                for charge_type in ['X','Z']:
                    if self.Stabilizers[type][measure_qubit]['charge'][charge_type]!= 0:
                        PhysicalError = True
                
        return PhysicalError

    ### distance Metric ###

    def distance(self, type, node1, node2):
        if node1 in self.Dual[type].nodes() and node2 in self.Dual[type].nodes():
            return nx.shortest_path_length(self.Dual[type], node1, node2)
        elif node1 in self.Dual[type].nodes() and node2 not in self.Dual[type].nodes():
            node2 = self.External[type][node2]['measure']
            return nx.shortest_path_length(self.Dual[type], node1, node2) + 1
        elif node1 not in self.Dual[type].nodes() and node2 in self.Dual[type].nodes():
            node1 = self.External[type][node1]['measure']
            return nx.shortest_path_length(self.Dual[type], node1, node2) + 1
        else:
            node1 = self.External[type][node1]['measure']
            node2 = self.External[type][node2]['measure']
            return nx.shortest_path_length(self.Dual[type], node1, node2) + 2


    # Re-initializes Measurement qubits
    def Assessment(self):
        if self.hasLogicalError():
            return False
        else:
            return True








