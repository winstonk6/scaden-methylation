#!/usr/bin/env python3
import argparse
parser = argparse.ArgumentParser(prog='array_to_sc', description='Converts methylation array to a simulated scMethyl-seq dataset.')
parser.add_argument('filename', help='Methylation array text file')
parser.add_argument('-r', '--repeats', type=int, default=1000, help='Number of times each beta value is sampled from a binomial distribution.')
parser.add_argument('-o', '--outname', default='./data',
                   help='Set the name and location of the outputs, e.g., "./data" will create the files "./data_counts.txt" and "./data_celltypes.txt"')
parser.add_argument('-s', '--seed', type=int, default=100, help='Numpy random generator seed')
args = parser.parse_args()

import numpy as np
import pandas as pd

rng = np.random.default_rng(args.seed)
filename = args.filename
repeats = args.repeats
outname = args.outname

# Read basis matrix
sig = pd.read_table(filename, index_col=0)

# Simulate single cell data by sampling binomial distribution
counts = np.empty(shape=(sig.shape[0], repeats*sig.shape[1]), dtype=np.int8)

for row in range(sig.shape[0]):
    for col in range(sig.shape[1]):
        counts[row, (col*repeats):(col*repeats+repeats)] = rng.binomial(1, sig.iat[row,col], repeats)

counts_df = pd.DataFrame(counts.T, index=[n for n in range(repeats*sig.shape[1])], columns=sig.index)
counts_df.to_csv(f'{outname}_counts.txt', sep='\t')

celltypes = pd.Series(np.repeat(sig.columns, repeats), name='Celltype')
celltypes.to_csv(f'{outname}_celltypes.txt', sep='\t', index=False)

print(f"Created files '{outname}_counts.txt' and '{outname}_celltypes.txt'")
