"""
Convert LD matrix to h5 format
Author: Gennady Khvorykh, info [at] inzilico [dot] com
"""

import h5py
import numpy as np
import sys, os
from helpers import show_time_elapsed
import time

# Check input provided
if len(sys.argv) != 3:
    print("""
    Usage: python3 ld2h5-02.py <input> <output>
        input: path/to/ld-matrix
        output: path/to/ld-matrix.h5
    """)
    sys.exit(0)

# Get command line arguments
x = sys.argv[1]
y = sys.argv[2]
ts = time.time()

# Check input 
if not os.path.isfile(x):
    print(x, "doesn't exist")
    sys.exit(1)

# Load ld-matrix into Numpy array
print("Loading", x)
arr = np.loadtxt(x, dtype = np.float32)
print("Shape: ", arr.shape)

# Save ld-matrix into h5 file
print("Writing", y)
with h5py.File(y, "w") as f:
    f.create_dataset("r2", data = arr)

show_time_elapsed(ts)