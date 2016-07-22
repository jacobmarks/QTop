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
# from geometry import *
from error_models import *
from kitaev import *
from color_codes import *

def main():
	code = KSC(10, 5)
	model = CodeCapacity(.05)
	code.plot_primal('Z', 1, 'Primal Lattice')
	plt.show()

# print code.hasLogicalError('Z')
# code = code.CodeCycle(model)

if __name__ == '__main__':
	main()
