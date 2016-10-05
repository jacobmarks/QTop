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
	return lambda p: [float(p)/1, float(p)/2, float(p)/3]


def main():
	errors = probs()
	model = ErrorModel(initialize = errors)
# print model.initialize(.2)

# code = Toric_6_6_6(8,2)
	code = KTC(5,2)

# There isn't an X component
# print code.dual.edges()


	model = CodeCapacity()
	code = code.CodeCycle(model,.1)
	matching = MinWeightMatch()
	decoder = MWPM_Decoder()
	code.plot_primal(1, 'Before')
	code = decoder(code)
	code.plot_primal(2, 'After')
	print Assessment(code)
	plt.show()


if __name__ == '__main__':
	main()



# 

# p_array = np.linspace(0.15,0.4,100)
# L_array = [3,5]
# num_trials = 10

# sim = simulation('kitaev', 'toric', 2, model, decoder)
# run(sim, L_array, p_array, num_trials)

# p = 0.01
# # for data in code.data:
# # 	print code.data[data]
# for p in np.linspace(0.08,0.15,10)
# 	failures = 
# 	for trial in range(1000):
# 	# initialize
# 	code = code.CodeCycle(model,p)
# 	code = code.CodeCycle(model, p)
# 	

# 	matching = MinWeightMatch()
# 	matches = matching(code, syn)

# 	recover = surface_recovery()

# 	A = recover(code, matches)
# 	print Assessment(code)


# plt.show()



