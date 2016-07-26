# Test decoder format

class decoder:

	def __init__(self, MWPM = 3):
		self.MWPM = MWPM


	def __call__(self, code):
		print code + 2 + self.MWPM

A = decoder()
def Decode(code, MWPM, decoder = A):
	decoder(code)

# Decode(4,2)

A = decoder()
A(3)

Decode(7,1,A)



# Make a decoder, with matching algorithm as optional parameter.
# RG as default

# then Decode(code, syndrome, decoder)