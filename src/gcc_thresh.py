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

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

L = np.array([0,2,3,5,13,25])
thresh = np.array([0,.13,.15,.19,.21,.23])
plt.plot(L, thresh, '.')

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

xs = np.linspace(0,25,26)
popt, pcov = curve_fit(func, L, thresh)

plt.plot(xs, func(xs, *popt), label="Fitted Curve")


title = "Threshold vs Qudit dimension"
plt.title(str(title))
plt.xlabel("Qudit dimension d")
plt.ylabel("Threshold")
plt.savefig('../plots/gcc_thresh.png')
plt.show()
