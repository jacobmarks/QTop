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
from common import *
from error_models import *

################# Surface Code Base Class ##############
# Up Left Right Down

def check0(qubit):
    return (qubit[0], qubit[1] + 1)

def check1(qubit):
    return (qubit[0] - 1, qubit[1])

def check2(qubit):
    return (qubit[0] + 1, qubit[1])

def check3(qubit):
    return (qubit[0], qubit[1] - 1)


class SurfaceCode(Code):
    'Common base class for all surface codes'

    def __init__(self, depth, dimension = 2):
        self.code = 'surface'
        Code.__init__(self, depth, dimension)

    def Plaquette(self, qubit):
        return {check0(qubit):0, check1(qubit):1, check2(qubit):2, check3(qubit):3}

    def generateColors(self):
        self.colors = {'X':'red', 'Z':'blue','data':'black'}


    def generateCode(self):
        self.types = ['X','Z']
        
        for type in self.types:
            self.Dual[type] = nx.Graph()
            self.Stabilizers[type] = {}
            self.Boundary[type] = {}
            self.External[type] = {}

        # Horizontal data qubits
        depth = self.depth

        for i in range(depth):
            for j in range(depth + 1):
                h = (2 * i, 2 * j - 1)
                self.Primal.add_node(h, charge=Charge())

        # Vertical data qubits
        for i in range(depth - 1):
            for j in range(depth):
                v = (2 * i + 1, 2 * j)
                self.Primal.add_node(v, charge=Charge())

        for i in range(depth):
            for j in range(depth):

                # Measure-Z qubits
                mZ = (2 * i, 2 * j)
                self.Stabilizers['Z'][mZ] = {'data': {}, 'charge': Charge(), 'order':{}, 'sides':4}
                for data in self.Plaquette(mZ):
                    if data in self.Primal.nodes():
                        count = self.Plaquette(mZ)[data]
                        self.Stabilizers['Z'][mZ]['order'][count] = data
                        self.Stabilizers['Z'][mZ]['data'][data] = count

                # External MeasureZ qubits
                if j == 0:
                    boundary_data, lower_boundary = (2 * i, -1), (2 * i, -2)
                    self.Boundary['Z'][boundary_data] = 'lower'
                    self.External['Z'][lower_boundary] = {'data': boundary_data, 'measure': mZ, 'order':0, 'sides':4}
                if j == depth - 1:
                    boundary_data, upper_boundary = (2 * i, 2 * depth - 1), (2 * i, 2 * depth)
                    self.Boundary['Z'][boundary_data] = 'upper'
                    self.External['Z'][upper_boundary] = {'data': boundary_data, 'measure': mZ, 'order':3, 'sides':4}


        for i in range(depth - 1):
            for j in range(depth + 1):

                # MeasureX qubits
                mX = (2 * i + 1, 2 * j - 1)
                self.Stabilizers['X'][mX] = {'data': {}, 'charge':Charge(), 'order':{}, 'sides':4}
                for data in self.Plaquette(mX):
                    if data in self.Primal.nodes():
                        count = self.Plaquette(mX)[data]
                        self.Stabilizers['X'][mX]['order'][count] = data
                        self.Stabilizers['X'][mX]['data'][data] = count

                # External MeasureX qubits
                if i == 0:
                    boundary_data, left_boundary = (0, 2 * j - 1), (-1, 2 * j - 1)
                    self.Boundary['X'][boundary_data] = 'left'
                    self.External['X'][left_boundary] = {'data': boundary_data, 'measure': mX, 'order':2, 'sides':4}

                if i == (depth - 2):
                    boundary_data, right_boundary = (2 * (depth - 1), 2 * j - 1), (2 * (depth - 1) + 1, 2 * j - 1)
                    self.Boundary['X'][boundary_data] = 'right'
                    self.External['X'][right_boundary] = {'data': boundary_data, 'measure': mX, 'order':1, 'sides':4}


        # Add Primal Lattice Edges
        for qubit in self.Primal.nodes():
            partners = [(qubit[0] - 1, qubit[1] - 1), (qubit[0] - 1, qubit[1] + 1), (qubit[0] + 1, qubit[1] - 1),
                        (qubit[0] + 1, qubit[1] + 1), (qubit[0] - 2, qubit[1]), (qubit[0] + 2, qubit[1]),
                        (qubit[0], qubit[1] - 2), (qubit[0], qubit[1] + 2)]
            for partner in partners:
                if partner in self.Primal.nodes():
                    self.Primal.add_edge(*(qubit, partner))

    def generateDual(self):
        for type in self.types:
            for node in self.Stabilizers[type]:
                self.Dual[type].add_node(node, Charge())
                partners = [(node[0] - 2, node[1]), (node[0] + 2, node[1]),
                            (node[0], node[1] - 2), (node[0], node[1] + 2)]
                for partner in partners:
                    if partner in self.Stabilizers[type]:
                        self.Dual[type].add_edge(*(node, partner))


    def CodeCycle(self, model, p = 0):

        # find length of code cycle:
        num_sides = 4


        # Step 1:
        self = model.Identity(self, p)
        self = model.Initialize(self, 'X', p)

        # # Step 2:
        self = model.Initialize(self, 'Z', p)
        self = model.Fourier(self, 'X', p)

        # # # Steps 3-6:
        for count in range(num_sides):
            for type in self.types:
                charge_type = type
                self = model.Sum(self, count, num_sides, type, charge_type, p)

        # # Step 7:
        self = model.Fourier(self, 'X', p)
        self = model.Measure(self, 'Z', p)

        # # Step 8:
        self = model.Measure(self, 'X', p)

        return self

    def hasLogicalError(self):
        for type in self.types:
            for charge_type in ['X','Z']:
                LogicalLattice = self.Primal.copy()
                for node in self.Primal.nodes():
                    if self.Primal.node[node]['charge'][charge_type] == 0:
                        LogicalLattice.remove_node(node)
                for node1 in LogicalLattice.nodes():
                    for node2 in LogicalLattice.nodes():
                        if node1 in self.Boundary[type] and node2 in self.Boundary[type]:
                            if self.Boundary[type][node1] != self.Boundary[type][node2]:
                                start, end = node1, node2
                                if start in LogicalLattice.nodes() and end in LogicalLattice.nodes():
                                    if nx.has_path(LogicalLattice, start, end):
                                        return True

        return False











