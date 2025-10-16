"""
    Count features of blocks.
    Author: Gennady Khvorykh, info@inzilico.com
    Created: October 8, 2024
"""

import sys, os
import pandas as pd
import numpy as np
import argparse
import h5py
from helpers import attach_matrix, check_input_files, check_output_dir, show_time_elapsed
import time

def parse_args():
    parser = argparse.ArgumentParser(description = 'Count features of blocks')
    parser.add_argument('-i', '--input', help = 'path/to/input.txt with two columns (index and labels)', 
                        required = True) 
    parser.add_argument('-o', '--output', help = 'path/to/output.txt', required = True)
    parser.add_argument('--h5', help = 'path/to/file.h5 HDF5 file with matrix', required=True)
    parser.add_argument('--r2', help = 'Minimal r2 value to estimate density of blocks', 
                        type = float, default = 0.7)
    parser.add_argument('-b', '--bim_file', help = 'path/to/filename.bim of Plink file')
    return parser.parse_args()
            
def count_features(row, matrix, r2):
    # Subset a block
    block = matrix[row.start:row.end, row.start:row.end]
    # Get above diagonal elements
    v = block[np.triu_indices(block.shape[0], k=1)]
    # Count density
    if len(v) == 0:
        rho = "NA"
        mean = "NA"
    else:
        rho = round(sum(1 for i in v if i >= r2) / len(v), 4)
        mean = round(np.mean(v), 4)
    return rho, mean         

def main():
    
    args = parse_args()
    
    # Initiate variables
    input_file = args.input
    output_file = args.output
    h5_file = args.h5
    bim_file = args.bim_file
    r2 = args.r2
    ts = time.time()
    
    # Check input
    check_input_files([input_file, h5_file, bim_file])
    
    # Check output folder exists
    check_output_dir(output_file)
    
    # Show input
    print("Input:", input_file)
    print("Output:", output_file)
    print("HDF5 file:", h5_file)
    if bim_file: print("bim file:", bim_file)
    print("Minimal r2:", r2)
        
    # Load indexes and labels
    dt1 = pd.read_csv(input_file, sep = " ", header = None)
    
    # Get min and max indexes and size of blocks 
    dt2 = dt1.groupby(1).agg({0: ['min', 'max'], 1: 'count'})
    dt2.columns = ["start", "end", "size"]
    print("Number of blocks:", dt2.shape[0])
   
    # Count other features
    with h5py.File(h5_file, "r") as f:
        # Load matrix
        matrix = attach_matrix(f)
        # Count rho and mean for each block
        dt3 = dt2.apply(lambda row: count_features(row, matrix, r2), axis = 1, result_type="expand")
    
    # Name columns
    dt3.columns = ["rho", "mean"]   
    
    if bim_file:
        # Load bim file
        bim = pd.read_csv(bim_file, sep = "\t", header = None)
        if bim.shape[0] != dt1.shape[0]:
            print("Different number of lines in bim and input files: ", bim.shape[0], dt1.shape[0])
            sys.exit(1)
        bim.columns = ["chr", "rs", "cm", "pos", "a1", "a2"]
        dt4 = dt2.apply(lambda row: bim.loc[row.end, "pos"] - bim.loc[row.start, "pos"], axis = 1, 
                        result_type="expand")
        dt4.name = "length"
        # Create output data frame
        out = pd.concat([dt2, dt3, dt4], axis = 1)
    else: 
        # Create output data frame
        out = pd.concat([dt2, dt3], axis = 1)
        
    # Save output
    out.to_csv(output_file, sep = " ", header = True, index_label="block")
    
    show_time_elapsed(ts)

if __name__ == '__main__':
    main()