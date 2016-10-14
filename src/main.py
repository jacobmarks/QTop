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
from kitaev import *
from decoders import *
import numpy as np
from simulation import *
from geometry import *

def probs():
  return lambda p: [float(p), 0, float(p)/2]

errors = probs()

def main():
	code = KTC(11,3)
	model = ErrorModel(initialize = errors)

	# model = CodeCapacity()
	code = code.CodeCycle(model,.1)
	# matching = MinWeightMatch()
	# decoder = MWPM_Decoder()
	code.plot_primal(1, 'Kitaev Toric Code')
	# code = decoder(code)
	# code.plot_primal(2, 'Primal')
	# code.plot_dual('Z', 3, 'Dual')
	# code.plot_shrunk('Z', 4, 'Shrunk')
	# print Assessment(code)
	plt.show()


if __name__ == '__main__':
	main()


