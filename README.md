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
QTop's base classes, including `Qubit`, `Stabilizer`, and `Code`, are all included in the 'common' module. In this first installment of the tutorial, you will be introduced to these classes and their methods. By the end, you will be able to create instances of topological quantum error correcting codes.

The simplest object one can work with in QTop is a qubit. A qubit is defined by its position, charge and type. By default, a qubit has trivial X and Z charges, and is of the type `'data'`. To instantiate a measure-X qubit with X charge 1 and Z charge 0 at position (5,2), type the command

```python
position = (5,2)
qubit = Qubit(position, charge = {'X':1, 'Z':0}, type = 'X')
```

A stabilizer is stored in the code by the position of its central measurement qubit, and contains information about the positions of the surrounding data qubits, as well as the order in which they must be measured during the error correction code cycle. In addition, the stabilizer has measurement type of the central qubit. 

As a simple example, let's say we want to create a stabilizer with 6 equally spaced data qubits, and central qubit as given above. To do this, we first generate the surrounding data qubits:

```python
data = generateStabilizerData(position, scale = 1, num_sides = 6, angle = 0)
```

This generates the following array, giving the positions of the data qubits:

```python
[(6.0, 2.0), (5.5, 2.866), (4.5, 2.866), (4.0, 2.0), (4.5, 1.134), (5.5, 1.134)]
```

Creating a stabilizer with this data is as simple as:

```python
stabilizer = Stabilizer('X', data, order = False)
```

By default, the measurement order is counterclockwise. However, one can program more exotic code cycle orderings.











---

