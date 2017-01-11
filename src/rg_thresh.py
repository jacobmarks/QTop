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

L = np.array([2,3,8,13,25,50])
thresh = np.array([.12,.157,.204,.213,.223,.227])
plt.plot(L, thresh, '.', label="Empirical Data")
# plt.plot(L, thresh)

def func(x, a, b, c):
    return a - float(b)/(c + x)
# def func(x, a, b, c):
#     return a * np.exp(-b * x) + c

# def func(x, a, b, c):
#     return float(a) /(1 + np.exp(-b* (x - c)))
# def func(x, a, b, c):
#     return a * ( 1 - np.exp(-b * x)) + c

xs = np.linspace(2,60,100)
popt, pcov = curve_fit(func, L, thresh)

plt.plot(xs, func(xs, *popt), label="Fitted Curve")

ys = [popt[0]] * 100
thr = round(popt[0],3)
plt.plot(xs, ys, 'r--', label="Plateau at " + str(thr))


title = "Threshold vs Qudit dimension"
plt.title(str(title))
plt.xlabel("Qudit dimension d")
plt.ylabel("Threshold")
plt.legend(loc=4)
plt.savefig('../plots/rg_thresh.png')
plt.show()
