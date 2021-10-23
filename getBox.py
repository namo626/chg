#!/usr/bin/env python3

import numpy as np

# Number of nodes
N = 2249093
filename = fort.14

print("Loading fort.14 ...")
arr = np.loadtxt(filename, skiprows=2, max_rows=N)
print("Minimum of fort.14:")
print(np.min(arr,0))
print("Maximum of fort.14:")
print(np.max(arr,0))
