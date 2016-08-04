# QTop
## Version 0.1 - 22 July 2016

QTop is an open-source python module for simulation and visualization of 
topological quantum codes. QTop is object-oriented and easy to read,
facilitating the addition and assessment of new topological computing
architectures, as well as the development of novel decoders.

QTop allows for the simulation of topologies with arbitrary code depth,
qudit dimension, and error models. Both planar and toric geometries are 
supported. Currently, QTop features Kitaev quantum double models, as well as
color codes in all 3-regular planar tilings. In the future, I hope to add
3-dimensional color codes, gauge color codes, and more exotic topoligical
systems.


Copyright (c) 2016 Jacob Marks, jacob DOT marks AT yale DOT edu.

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
QTop's base classes, including `qudit`, `Stabilizer`, and `Code`, are all included in the 'common' module. In this first installment of the tutorial, you will be introduced to these classes and their methods. By the end, you will be able to create instances of topological quantum error correcting codes.

The simplest object one can work with in QTop is a qudit. A qudit is defined by its position, charge and type. By default, a qudit has trivial X and Z charges, and is of the type `'data'`. To instantiate a measure-X qudit with X charge 1 and Z charge 0 at position (5,2), type the command

```python
position = (5,2)
qudit = qudit(position, charge = {'X':1, 'Z':0}, type = 'X')
```

A stabilizer is stored in the code by the position of its central measurement qudit, and contains information about the positions of the surrounding data qudits, as well as the order in which they must be measured during the error correction code cycle. In addition, the stabilizer has measurement type of the central qudit. 

As a simple example, let's say we want to create a stabilizer with 6 equally spaced data qudits, and central qudit as given above. To do this, we first generate the surrounding data qudits:

```python
data = generateStabilizerData(position, scale = 1, num_sides = 6, angle = 0)
```

This generates the following array, giving the positions of the data qudits:

```python
[(6.0, 2.0), (5.5, 2.866), (4.5, 2.866), (4.0, 2.0), (4.5, 1.134), (5.5, 1.134)]
```

Creating a stabilizer with this data is as simple as:

```python
stabilizer = Stabilizer('X', data, order = False)
```

By default, the measurement order is counterclockwise, but one can program more complex code cycle orderings.

In QTop, `Code` is the base class for topological error correcting codes. From this, surface and color codes are derived as subclasses. In addition, this base class can easily be extended to 3D color codes, gauge color codes, and more exotic topological codes. `Code` contains all information about the data and measurement qudits, the dimension of each qudit, and the code depth, i.e. the smallest number of physical qubits that must present errors for there to be a logical error. In addition, `Code` stores the network connections on the primal and dual lattices. This becomes relevant later, because weightings in the decoding procedure are determined based on shortest path length on the dual lattice.

Let's say we want to instantiate a kitaev surface code of code depth 11, and qudit dimension 3. The relevant subclass of `Code` is `KSC`. 

```python
code = KSC(depth = 11, dimension = 3)
```

We can view this code with the plotting methods

```python
code.plot_primal(1, 'Kitaev Surface Code')
plt.show()
```

which produces the following plot:

![alt text](/visualizations/Kitaev_Surface_Code.png)

Circles are associated with physical qubits inside the lattice, while diamonds are external elements, which are quasi-excitations that are utilized in the decoding process, but don't correspond to anything physical. Black represents data, while blue and red respectively represent measure-X and measure-Z qudits. Edges are connections between data qudits.

If instead, we wanted to view only the dual lattice, we could type:

```python
code.dual(1, 'Kitaev Surface Code')
plt.show()
```

## Error Models

Once we have our code, we need to choose an error model under which to simulate the code's time evolution. In QTop, error models are represented by objects containing the Pauli X, Y and Z error probabilities associated with each gate. All error model objects are instances of the base class `ErrorModel`, or a subclass thereof. Inherent subclasses include `CodeCapacity`, `Phenomenological`, and `CircuitLevel`.

If we want to create our own error model, we just need to specify the gates with non-zero error probabilities as inputs to `ErrorModel`.









---

