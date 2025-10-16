"""
Check and clean from NaN values in square matrix.
Author: Gennady Khvorykh, info [at] inzilico [dot] com
Created: 2025-02-20
"""

import sys
import h5py
import argparse
import numpy as np
import pandas as pd
from helpers import check_input_files

def get_cla():
    # Define the command line arguments
    parser = argparse.ArgumentParser(description='Check and clean from NaN values in square matrix.')
    parser.add_argument('-p', '--prefix', type=str, help='Prefix of input and output files', required=True)
    return parser.parse_args()

def main():
    # Get command line arguments
    args = get_cla()
    
    # Initiate variables
    prefix = args.prefix
    h5_file = prefix + ".ld.h5"
    bim_file = prefix + ".bim"
    cleaned_file_name = prefix + ".cleaned.ld.h5"
    snp_file = prefix + ".snp"

    # Check input files
    check_input_files([h5_file, bim_file])

    # Connect to HDF5 file
    with h5py.File(h5_file, 'r') as file_obj:
        # Check dataset is a square matrix
        dataset_name = "r2"
        dataset = file_obj[dataset_name]
        N, M = dataset.shape

        if N != M:
            print("Dataset is not a square matrix")
            sys.exit(1)
        else:
            print(f"Original dataset shape: {N} x {M}")

        # Remove rows and columns where all elements are NaN
        dataset_array = dataset[:]  # Load the entire dataset into memory
        

    # Get logical vector for non-NaN rows   
    non_nans = ~np.isnan(dataset_array).all(axis=1)

    # Compare the number of non-NaN rows with the total number of rows
    if non_nans.sum() == N:
        print("No rows with all NaN values found")
        # Import the bim data 
        df = pd.read_csv(bim_file, sep="\t", header=None)
        # Save the second column 
        df.iloc[:, 1].to_csv(snp_file, index=False, header=False)
        print(f"The original list of SNPs saved to {snp_file}")
        sys.exit(0)
    else:
        # Remove rows and columns with all NaN values from the dataset array
        cleaned_dataset = dataset_array[non_nans][:,non_nans]

        # Print information about removed rows and columns
        removed_count = N - cleaned_dataset.shape[0]
        if removed_count > 0:
            print(f"Removed {removed_count} rows and {removed_count} columns with all NaN values")
            print(f"New dataset shape: {cleaned_dataset.shape[0]} x {cleaned_dataset.shape[1]}")
        else:
            print("No rows or columns with all NaN values found")

        # Check remaining nan values in dataset
        nan_values = np.isnan(cleaned_dataset)
        nan_count = np.sum(nan_values)
        if nan_count > 0:
            print(f"{nan_count} NaN values found in the cleaned dataset")
        else:
            print("No NaN values found in the cleaned dataset")

        # Save cleaned dataset to a new HDF5 file
        with h5py.File(cleaned_file_name, 'w') as file_obj:
            cleaned_dataset_dset = file_obj.create_dataset(dataset_name, 
                                                        cleaned_dataset.shape, 
                                                        dtype=cleaned_dataset.dtype)
            cleaned_dataset_dset[:] = cleaned_dataset

        print(f"Cleaned dataset saved to {cleaned_file_name}")
        # Import the bim data
        df = pd.read_csv(bim_file, sep="\t", header=None)
        # Subset the rows without NaN values
        df = df.iloc[non_nans]
        # Save the second column
        df.iloc[:, 1].to_csv(snp_file, index=False, header=False)
        print(f"Cleaned list of SNPs saved to {snp_file}")

if __name__ == '__main__':
    main()