#!/usr/bin/env python3
## Partition matrix into submatrices by the determinant.
## The input file has square matrix with 1 on diagonal under h5py format.
## The output file has two columns: index and labels of blocks.
## Author: Gennady Khvorykh, info@inzilico.com
## Date created: 11.10.2024

import sys, os
import numpy as np
import h5py
import pandas as pd
import argparse
import time
from helpers import check_input_files, check_output_dir, attach_matrix, show_time_elapsed


def parser_args():
    # Set command line arguments
    parser = argparse.ArgumentParser(description="Partition square matrix by determinant")
    parser.add_argument("-i", "--input", help="/path/to/input.h5 with matrix under h5py format", required=True)
    parser.add_argument("-o", "--output", help="/path/to/output.txt to save block labels", required=True)
    parser.add_argument("-d", "--min_det", help="Minimal determinant value (default = 0.001)", 
                        default=0.001, type=float)
    # Get arguments
    return parser.parse_args()

def partition_by_determinant(matrix, min_det: float):
    """
    Partition square matrix by determinant into submatrices (blocks).
    Parameters
    ----------
    matrix : h5py object
        Square matrix with 1 on diagonal
    min_det : float
        Minimal determinant value
    Returns
    -------
    list
        Labels of matrix elements
    """

    # Initiate variables
    i1, i2 = 0, 1
    b = 1
    labels = [[0, b]]
    n = matrix.shape[0]
    
    # Run over main diogonal
    while(i2 < n):
        # Get submatrix
        submatrix = matrix[i1:i2, i1:i2]
        # Count determinant
        sign, logdet = np.linalg.slogdet(submatrix)
        if (sign <= 0) | (logdet <= np.log(min_det)):
            labels.append([i2, b])
            # Start new block
            i2 += 1
            if i2 >= n: break
            b += 1 
            labels.append([i2, b])
            i1 = i2
            i2 += 1
        else:
            labels.append([i2, b])
            i2 += 1
    
    return labels, b        
        
def main():

    # Get arguments
    args = parser_args()

    # Initiate variables
    input_file = args.input
    output_file = args.output
    min_det = args.min_det
    ts = time.time()

    # Check input files
    check_input_files([input_file])

    # Create output folder if not exists
    check_output_dir(output_file)

    # Show input
    print("Input:", input_file)
    print("Output:", output_file)
    print("Minimal determinant:", min_det)

    # Open input file, attach and process matrix 
    with h5py.File(input_file, "r") as f:
        matrix = attach_matrix(f)
        labels, nblocks = partition_by_determinant(matrix, min_det=min_det)
    
    # Save labels of blocks
    np.savetxt(output_file, labels, delimiter=" ", fmt="%d")
    
    # Show results
    print("Number of blocks:", nblocks)
    
    show_time_elapsed(ts)


if __name__ == "__main__":
    main()