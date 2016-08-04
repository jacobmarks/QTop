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
from geometry import *
from error_models import *
from kitaev import *
from color_codes import *
from decoders import *

def main():
	code = KTC(8, 2)
	
	
	model = CodeCapacity(.01)
	code = code.CodeCycle(model)
	syn = syndrome(code)

	MW_cluster = MinWeightMatch()
	matching = MW_cluster(code, syn)
	print matching
	# rec = surface_recovery()
	# code = rec(code, matching)

	decoder = MWPM_Decoder(MW_cluster)
	# code = decoder(code, syn)
	Decode(code, syn, decoder)
	code.plot_shrunk('X', 'Z', 1, 'Shrunk Z Lattice')
	code.plot_shrunk('Z', 'X', 2, 'Shrunk X Lattice')
	perfect_gates = ErrorModel()
	code = code.CodeCycle(perfect_gates)
	# code = Decode(code, syn, decoder)

	# code.plot_primal('Z', 1, 'Primal Z Lattice')
	# code.plot_primal('blue', 2, 'Primal Z Lattice')
	# code.plot_dual('X', 3, 'Dual X Lattice')
	code.plot_shrunk('X', 'Z', 3, 'Shrunk Z Lattice')
	code.plot_shrunk('Z', 'X', 4, 'Shrunk X Lattice')
	# code.plot_shrunk('green', 'Z', 5, 'Shrunk green Lattice')
	# code.plot_shrunk('blue', 'Z', 6, 'Shrunk blue Lattice')
	plt.show()

# print code.hasLogicalError('Z')


if __name__ == '__main__':
	main()
