#!/usr/bin/env python3
## Partition genome with LD matrix by determinant
## Author: Gennady Khvorykh, info@inzilico.com
## Date created: 13.09.2022

import sys, os
import numpy as np
import h5py
import pandas as pd
import argparse
import time

# Initiate
parser = argparse.ArgumentParser(description="Partition genome by determinant")
parser.add_argument("--h5", help="/path/to/ld-matrix.h5", required=True)
parser.add_argument(
    "-p", "--pref", help="/path/to/prefix to save output files", required=True
)
parser.add_argument("-b", "--bim", help="/path/to/filename.bim of Plink file")
parser.add_argument("-t", "--trm2", action="store_true", help="switch to count Tr(M^2)")
args = parser.parse_args()

h5 = args.h5
pref = args.pref
bim = args.bim
trm2 = args.trm2

min_det = 0.001

t1 = time.time()

# Check input
for x in [h5, bim]:
    if x:
        if not os.path.isfile(x):
            print(x, "wasn't found")
            sys.exit(1)

# Open h5 file
h = h5py.File(h5, "r")
ds = h["r2"]
d = ds["block0_values"]

# Load bim data
if bim:
    bm = pd.read_csv(bim, sep="\s+", header=None)

# Open output files
f1 = open(f"{pref}.txt", "w")
if bim:
    f2 = open(f"{pref}.csv", "w")
    print("chr,pos,rs,cl", file=f2)

# Initiate variables for loop
N = d.shape[0]
b, m = 1, 0

# Show input
print("LD matrix:", h5)
print("Size:", N)
if bim:
    print("bim file:", bim)
    if N != bm.shape[0]:
        print(f"Different number of lines in {h5} and {bim}: {N} {bm.shape[0]}")
        sys.exit(1)

# Loop over an axis of the matrix
for n in range(N):
    if m == n:
        if bim:
            print(bm.iloc[n, 0], bm.iloc[n, 3], bm.iloc[n, 1], b, sep=",", file=f2)
        continue
    # Subset a matrix
    d1 = d[m:n, m:n]
    # Count determinant
    sign, logdet = np.linalg.slogdet(d1)

    # Check determinant
    if logdet < np.log(min_det):
        if sign == 0:
            dc = d[m : n - 1, m : n - 1]
            sign, logdet = np.linalg.slogdet(dc)
            if bim:
                print(
                    bm.iloc[n - 1, 0],
                    bm.iloc[n - 1, 3],
                    bm.iloc[n - 1, 1],
                    b,
                    sep=",",
                    file=f2,
                )
            print(b, m, n - 1, sign, logdet, file=f1)
            b += 1
            m = n
        else:
            # Save to file
            if bim:
                print(bm.iloc[n, 0], bm.iloc[n, 3], bm.iloc[n, 1], b, sep=",", file=f2)
            print(b, m, n, sign, logdet, file=f1)
            b += 1
            m = n + 1
            dc = d1
    else:
        # Save to file
        if bim:
            print(bm.iloc[n, 0], bm.iloc[n, 3], bm.iloc[n, 1], b, sep=",", file=f2)
        print(b, m, n, sign, logdet, file=f1)

# Close files
f1.close()
if bim:
    f2.close()
h.close()

# Show time elapsed
dur = time.time() - t1
dur = time.strftime("%H:%M:%S", time.gmtime(dur))
print("Time spent: ", dur)
