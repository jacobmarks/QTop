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
import matplotlib.pyplot as plt

L = [2,3,5,13,25]
thresh = [0.12,.15,.16,.20,.23]
plt.plot(L, thresh)
title = "Threshold vs Qudit dimension"
plt.title(str(title))
plt.xlabel("Qudit dimension d")
plt.ylabel("Threshold")
plt.savefig('../plots/rg_thresh.png')
plt.show()
