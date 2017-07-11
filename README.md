# QTop
## Version 0.2 - 11 July 2017

QTop is an open-source python module for simulation and visualization of 
topological quantum codes. QTop is object-oriented and easy to read,
facilitating the addition and assessment of new topological computing
architectures, as well as the development of novel decoders.

QTop allows for the simulation of topologies with arbitrary code depth,
qudit dimension, and error models. Currently, QTop features Kitaev quantum 
double models, as well as color codes in 3-regular planar tilings. In 
the future, I hope to add 3-dimensional color codes, gauge color codes, and 
more exotic topological systems.

Special thanks to Vlad Gheorghiu for his help with plotting, and Tomas Jochym-O'Connor for help implementing DSP and GCC decoders.


Copyright (c) 2017 Jacob Marks, jamarks AT stanford DOT edu.

---
QTop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

QTop is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with QTop.  If not, see <http://www.gnu.org/licenses/>.

---

# Brief Tutorial

## Basics
QTop's base classes are `Code`, `ErrorModel`, and `decoder`. The methods of `Code` are included in the 'common' module. In this first installment of the tutorial, you will learn how to create instances of topological quantum error correcting codes.

In QTop, `Code` is the base class for topological error correcting codes. From this, surface and color codes are derived as subclasses. In addition, this base class can easily be extended to 3D color codes, gauge color codes, and more exotic topological codes. `Code` contains all information about the data and measurement qudits, the dimension of each qudit, and the code depth, i.e. the smallest number of physical qubits that must present errors for there to be a logical error. In addition, `Code` stores the network connections on the primal and dual lattices.

Let's say we want to instantiate a 6-6-6 color code of code depth 11, and qudit dimension 3. The relevant subclass of `Code` is `Color666`, and this class resides in the `color_codes` module. 

```python
code = color_codes.Color666(11, 3)
```

We can view this code with the plotting methods

```python
visualization.PlotPlaquette(code, "Color Code")
plt.show()
```

which produces the following plot:

![alt text](/visualizations/plaquette.png)

Black dots represent data qubits. Measurement qubits, which come in types 'red', 'blue', and 'green', lie at the center of each hexagonal stabilizer. The edges of each color are the edges of the dual shrunk lattice of that color. Diamonds of a given color show external 'pseudo' measurement qudits, and their 'pseudo' connections to internal measurement qudits are pictured as dashed lines.

If instead, we wanted to view only the primal lattice (connections between data qudits) we could type:

```python
visualization.PlotPrimal(code, "Color Code")
plt.show()
```

![alt text](/visualizations/primal.png)

## Error Models

Once we have our code, we need to choose an error model under which to simulate the code's time evolution. In QTop, error models are represented by objects containing the Pauli X, Y and Z error probabilities associated with each gate. All error model objects are instances of the base class `ErrorModel`, or a subclass thereof. Inherent subclasses include `CodeCapacity`, `Phenomenological`, and `CircuitLevel`, which are the three most commonly used error models in the literature.

If we want to create our own error model, we just need to specify the faulty gates in `ErrorModel`. Suppose we want to create an error model with faulty initialization, which presents a Pauli X error with probability p, and Pauli Z error with probability p/2 for some p. 
QTop ErrorModel objects take as input lambda functions, not probabilities. This way, you don't have to create a new object each time you want to run a simulation with a new probability.

```python
def probs():
  return lambda p: [float(p), 0, float(p)/2]

errors = probs()
```




Then we instantiate our model:

```python
model = error_models.ErrorModel(initialize = errors)
```

We can simulate our code under this error model with the `CodeCycle` method:

For simplicity, we will set `p = 0.1`.


```python
code = code.CodeCycle(model, p)
```

after this, errors in our code, i.e. non-trivial eigenvalue measurements, are indicated by stars. The size of the star scales with the magnitude of the error.

![alt text](/visualizations/before_decoding.png)

If we want to find all 'red' measurement qudits with Pauli Z errors, we write

```python
s = code.Syndrome('red', 'Z')
print s
```

## Decoders

At this point, we have a code with errors given by the action of our error model. We want to return the code to its codespace by correcting for those errors in an efficient and effective way. To do this, we need to apply homologically trivial correction chains.

Every decoding object has a matching algorithm - the algorithm that matches up elements of the overall error syndrome. In the case of the qubit surface code, an effective algorithm is Edmunds' Blossom algorithm for minimum weight perfect matching, which splits the syndrome into localized pairs of measurement qubits.

Alternatively, renormalization group clustering identifies all maximally disjoint neutral clusters at a given distance scale, and then pairs off elements of the same charge within these clusters. One element in this pair has its excitation charge 'transported' to the other one. This is continued until these neutral clusters are annihilated. Then, the scale is increased and the process starts again. This ensues until the syndrome is empty.

For the purpose of demonstration, we will use our novel Generalized Color Clustering (GCC) decoder. All decoders lie in the `decoders` module:

```python
decoder = gcc.GCC_decoder()
```

Then we perform error correction by applying this decoder:

```python
code = decoder(code)
```

This leads to the corrected lattice represented below:

![alt text](/visualizations/after_decoding.png)

Finally, we can assess our decoding by checking to see if any logical errors have occurred. 

```python
print code.hasLogicalError()
```

Our above decoding was successful in returning the code to the codespace, as can be seen by the presence of only combinations of stabilizers in the corrected lattice. This means that there was no logical error. If, instead we had started with the pre-decoding code

![alt text](/visualizations/before_error.png)

our decoding procedure would result in a collection of excitations that anticommutes with the logical operators supported on the code. This is visualized by a string of data excitations connecting boundaries of all three types.

![alt text](/visualizations/after_error.png)


## Putting it all Together

In this brief tutorial, we've introduced the fundamentals of codes, error models, and decoders. These are the three basic ingredients of a `simulation`. Then, `sim = simulation(2, "666 Color Code", [model, "Test Model"], [decoder, "GCC"])` instantiates a simulation object. The first two input arguments are used to construct the code itself. We can run a simulation by calling our simulation object, with a code depth `L`, and a physical error rate `p`. We succeed if we decode our code without any logical errors, in which case `sim(L, p)` returns `True `. Otherwise, we fail.

We can `run` this simulation by specifying a range of probabilities and code depths, and the number of trials to conduct for each depth-probability combination.

```python
import numpy as np
depths = np.linspace(5,55,6)
probabilities = np.logspace(-3,-0.5,10)
trials = 1000
run(sim, depths, probabilities, trials)
```
This will store a time-stamped pickled file of the logical error rate for each depth-probability combination. If this is stored in the file `data.pickle`, then we can plot out results and obtain a threshold estimate with uncertainties by typing `plot_results.py data.pickle output.png mode` in the command line, where `mode` specifies the type of the plot and must be either `linear`, `loglog` or `semilog`.


---









---

