
from common import *
from error_models import *
from kitaev import *


# charge = Charge(2,4)
# position = (5,2)
# type = 'X'
# data = generateStabilizerData(position, 1, 6, 0)
# print data

# qubit = Qubit(position, charge, type)
# print 'qubit attributes:'
# print 'position:', qubit.position
# print 'X charge:', qubit.charge['X']
# print 'Z charge:', qubit.charge['Z']
# print 'type:', qubit.type

# data = [(0,0),(1,1),(2,2),(3,3)]
# order = Order(data)

# stabilizer = Stabilizer(type, data, order)
# print 'stabilizer attributes:'
# print 'type:', stabilizer.type
# print 'data:', stabilizer.data
# print 'order:', stabilizer.order


# Error Model Works
# model = CodeCapacity(.3)
# print model.identity
# print model.initialize

# dim = 3
# error_rates = [1,0,0]
# Paulis and BP_Channel work
# qubit = pauli_X(qubit,3)
# qubit = BP_Channel(qubit, dim, error_rates)

# print 'X charge:', qubit.charge['X']
# print 'Z charge:', qubit.charge['Z'] 
def probs():
	return lambda p: [float(p)/1, float(p)/2, float(p)/3]

errors = probs()
# model = ErrorModel(initialize = errors)
# print model.initialize(.3)

Surface = KSC(11,3)
model = ErrorModel(identity = errors)
p = 0.03
# for data in Surface.data:
# 	print Surface.data[data]
Surface = Surface.CodeCycle(model, p)
Surface.plot_primal(1, 'Kitaev Surface Code')
plt.show()
# print Surface.primal.nodes()
# print Surface.data
# Surface = Surface.CodeCycle(model)
# for position in Surface.data:
# 	for type in Surface.data[position].charge:
# 		if Surface.data[position].charge[type] != 0:
# 			print position, type, Surface.data[position].charge[type]
# for type in Surface.types:
# 	for measure in Surface.stabilizers[type]:
# 		print Surface.stabilizers[type][measure].data

# print Surface.stabilizers
# Surface.plot_primal('Z', 1, 'Primal Z Lattice')
# Surface.plot_dual('X', 3, 'Dual X Lattice')
# Surface.plot_shrunk('X', 'X', 1, 'Shrunk X Lattice')
# plt.show()


# class TEST:

# 	def __init__(self, initialize):
# 		self.initialize = initialize

# 	def Initialize(self, p):
# 		print self.initialize(p)

# Probs = lambda p: [float(p)/n, float(p)/n, float(p)/n]
# def probs():
# 	return lambda p: [float(p)/1, float(p)/2, float(p)/3]

# f = probs()
# print f

# A = TEST(f)
# A.Initialize(3)



# PerfectGate = lambda p: 0
# print PerfectGate(5)













