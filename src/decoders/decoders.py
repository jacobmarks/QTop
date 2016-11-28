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

# from common import *
from matplotlib import path
from math import floor
import sys
sys.path.append('../')
sys.path.append('../../')

from src import common

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


def reset_measures(code):
    for type in code.types:
        for measure_qubit in code.Stabilizers[type]:
            code.Stabilizers[type][measure_qubit]['charge'] = common.Charge()
    return code



class matching_algorithm(object):

    def __init__(self):
        pass




