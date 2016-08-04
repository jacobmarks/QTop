
from common import *
from error_models import *
from kitaev import *


charge = Charge(2,4)
position = (3,5)
type = 'red'

qubit = Qubit(position, charge, type)
# print 'qubit attributes:'
# print 'position:', qubit.position
# print 'X charge:', qubit.charge['X']
# print 'Z charge:', qubit.charge['Z']
# print 'type:', qubit.type

data = [(0,0),(1,1),(2,2),(3,3)]
order = Order(data)

stabilizer = Stabilizer(type, data, order)
# print 'stabilizer attributes:'
# print 'type:', stabilizer.type
# print 'data:', stabilizer.data
# print 'order:', stabilizer.order


# Error Model Works
model = CodeCapacity(.3)
# print model.identity
# print model.initialize

dim = 3
error_rates = [1,0,0]
# Paulis and BP_Channel work
# qubit = pauli_X(qubit,3)
# qubit = BP_Channel(qubit, dim, error_rates)

# print 'X charge:', qubit.charge['X']
# print 'Z charge:', qubit.charge['Z'] 
Surface = KTC(5)

# print Surface.primal.nodes()
# print Surface.data
Surface = Surface.CodeCycle(model)
for position in Surface.data:
	for type in Surface.data[position].charge:
		if Surface.data[position].charge[type] != 0:
			print position, type, Surface.data[position].charge[type]
# for type in Surface.types:
# 	for measure in Surface.stabilizers[type]:
# 		print Surface.stabilizers[type][measure].data

# print Surface.stabilizers
# Surface.plot_primal('Z', 1, 'Primal Z Lattice')
# Surface.plot_dual('X', 3, 'Dual X Lattice')
Surface.plot_shrunk('X', 'X', 1, 'Shrunk X Lattice')
plt.show()

