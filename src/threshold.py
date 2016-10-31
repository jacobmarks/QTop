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

import numpy as np
from scipy.optimize import curve_fit

def form(X, p_th, A,B,C,D,u,v):
    return A + B*(X[0] - p_th)* np.power(X[1], float(1)/v)+ C*np.power((X[0] - p_th) * np.power(X[1], float(1)/v),2)+ D*np.power(X[1], float(-1)/u)