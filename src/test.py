
from common import *
from error_models import *
from kitaev import *
from decoders import *
import numpy as np
from simulation import *
from geometry import *




def probs():
	return lambda p: [float(p)/1, float(p)/2, float(p)/3]

errors = probs()
# model = ErrorModel(initialize = errors)
# print model.initialize(.3)

# code = Triangular_6_6_6(15,5)
code = KSC(15,5)
model = CodeCapacity()
code = code.CodeCycle(model,.03)
matching = HDRG_Match()
decoder = MWPM_Decoder()
# decoder = DSP(matching)
decoder(code)
code.plot_dual(1, 'Kitaev Toric Code')
plt.show()
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















