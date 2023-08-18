#!/usr/bin/env python3
import time
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import tarfile
import gzip
import io

parser = argparse.ArgumentParser(prog='scMap-mem', description='Map snmC-seq to bulk loci. Only does exact matching.')
parser.add_argument('array', 
                    help='Reference array (e.g., Illumina 450K) in BED format (without headers or column names, containing only chrom, chromStart, chromEnd, name, and strand)')
parser.add_argument('singlecell', help='Path to folder/tarball containing the single cell data')
parser.add_argument('-o', '--out', help='Name and path of output file')
parser.add_argument('-s', '--strand', action='store_true', help='Enforce strand specificity')
parser.add_argument('-m', '--mc_classes', default=None, 
                    help='Count the frequency of mc-classes for all positions, and store to the results to the specified output file.')
parser.add_argument('-c', '--coverage', default=None, 
                    help='Count the total number of positions in a cell, then append the result to the specified output file.')
parser.add_argument('-a', '--add_header', action='store_true', 
                    help='Use if the chromosome files do not have a row for the column names (chr, pos, strand, mc_class, mc_count, total, methylated)')
args = parser.parse_args()


# Read each chromosome file and return a df with the required columns
def process_df(file, array_df: pd.DataFrame, strand: bool, mc_classes: bool = False, coverage: bool = False, add_header: bool = False) -> pd.DataFrame:
    if add_header:
        df = pd.read_csv(file, sep='\t', names=['chr', 'pos', 'strand', 'mc_class', 'mc_count', 'total', 'methylated'], na_filter=False, dtype={'chr': str})
    else:
        df = pd.read_csv(file, sep='\t', na_filter=False, dtype={'chr': str})
    
    if coverage:
        global line_count
        line_count += df.shape[0]
    
    # Add 'chr' to chromosome names if names don't already have it
    if df.at[0, 'chr'][:3] != 'chr':
        df['chr'] = 'chr' + df['chr']
    
    if strand:
        df = pd.merge(df, array_df, how='inner', left_on=['chr', 'pos', 'strand'], right_on=['chrom', 'chromStart', 'strand2'])
    else:
        df = pd.merge(df, array_df, how='inner', left_on=['chr', 'pos'], right_on=['chrom', 'chromStart'])
    
    df['pos2'] = df['pos'] + 1
    df['frac'] = df['mc_count'] / df['total']
    
    # Get mc_classes frequency after filtering 
    if mc_classes:
        global mc_counts
        mc_counts = pd.concat([mc_counts, df['mc_class']])
    return df[['chr', 'pos', 'pos2', 'name', 'strand', 'frac']]


# Open chromosome tsv file and append to sc_unsorted_df.
def chromosome_append(chromosome_file, array_df: pd.DataFrame, sc_unsorted_df: pd.DataFrame, 
                      strand: bool = True, mc_classes: bool = False, coverage: bool = False, add_header: bool = False) -> pd.DataFrame:
    try:
        with gzip.open(chromosome_file, 'rb') as chromosome:
            df = process_df(file=chromosome, strand=strand, array_df=array_df, mc_classes=mc_classes, coverage=coverage, add_header=add_header)
    except: 
        # If just a tsv file
        df = process_df(file=chromosome_file, strand=strand, array_df=array_df, mc_classes=mc_classes, coverage=coverage, add_header=add_header)
    return pd.concat([sc_unsorted_df, df])
    
    

def handle_tarfile(tfile, array_df: pd.DataFrame, sc_unsorted_df: pd.DataFrame, strand: bool, mc_classes: bool = False, 
                   coverage: bool = False, add_header: bool = False, fileobj: bool = False) -> pd.DataFrame:
    # Read tarfile from file or memory
    if fileobj:
        tar = tarfile.open(fileobj=io.BytesIO(tfile))
    else:
        tar = tarfile.open(tfile)
    
    df = sc_unsorted_df
    for member in tar:
        if member.isdir():
            # Skip directories since python gets all contained files.
            continue
        mem = tar.extractfile(member)
        if tarfile.is_tarfile(mem):
            # Recursively open tarfiles
            file_content_byte = tar.extractfile(mem).read()
            df = handle_tarfile(tfile=file_content_byte, array_df=array_df, sc_unsorted_df=df, 
                                strand=strand, mc_classes=mc_classes, coverage=coverage, add_header=add_header, 
                                fileobj=True)
        else:
            df = chromosome_append(chromosome_file=mem, array_df=array_df, sc_unsorted_df=df, 
                                   strand=strand, mc_classes=mc_classes, coverage=coverage, add_header=add_header)
    tar.close()
    return df
    

def main():
    path = Path(args.singlecell).expanduser()
    array = args.array
    array_df = pd.read_csv(array, sep='\t', header=None, names=['chrom', 'chromStart', 'chromEnd', 'name', 'strand2'])
    if args.out:
        out_real = Path(args.out).expanduser()
    else:
        out_real = Path(args.singlecell + '.bed')

    mc_classes = args.mc_classes != None # True or False for whether or not to do counting. If args.mc_classes is provided, count is True
    if mc_classes:
        # Create variable to store all mc_classes
        global mc_counts
        mc_counts = pd.Series([], dtype=str)
        count_out = Path(args.mc_classes).expanduser()
    
    coverage = args.coverage != None
    if coverage:
        # Create variable to store total line count
        global line_count
        line_count = 0
        coverage_out = Path(args.coverage).expanduser()
        if not coverage_out.exists():
            with open(coverage_out, 'w') as f:
                f.write('cell\tcoverage\n')
    
    scUnsorted = pd.DataFrame([])

    print(f'Parsing data in {path} ...')
    # If we have a pool_x folder containing one cell's 23 chromosome.tsv.gz
    if path.is_dir():
        for chromosome in path.iterdir():
            if not chromosome.is_dir():
                scUnsorted = chromosome_append(chromosome_file=chromosome, array_df=array_df, 
                                               sc_unsorted_df=scUnsorted, strand=args.strand, mc_classes=mc_classes, 
                                               coverage=coverage, add_header=args.add_header)

    # If we have a tar or tar.gz file containing the pool_x folder or tsv files
    elif tarfile.is_tarfile(path):
        scUnsorted = handle_tarfile(tfile=path, array_df=array_df, sc_unsorted_df=scUnsorted,
                                    strand=args.strand, mc_classes=mc_classes, coverage=coverage,
                                    add_header=args.add_header)

    # If we have a single tsv file containing the data for all chromosomes
    elif path.is_file():
        scUnsorted = chromosome_append(chromosome_file=path, array_df=array_df, sc_unsorted_df=scUnsorted, 
                                       strand=args.strand, mc_classes=mc_classes, coverage=coverage, 
                                       add_header=args.add_header)
    
    print(f'Extracted all single-cell data from {path}')
    
    
    # Get frequency of mc_classes and write to file
    if mc_classes:
        mc_counts.value_counts().to_csv(count_out, sep='\t', index_label = 'mc_class', header=['count'])
        print(f'Saved mc_classes count to {args.mc_classes}')
    
    # Count total number of rows in the cell's files
    if coverage:
        with open(coverage_out, 'a') as outfile:
            outfile.write(f'{path.stem}\t{line_count}\n')
        print(f'Saved coverage count to {args.coverage}')
    
    # Process single_cell data to tsv file containing mapped locations.
    print('Mapped single-cell data to reference array')
    scUnsorted.sort_values(by=['chr', 'pos']).to_csv(out_real, sep='\t', header=False, index=False)
    print(f'Saved output to {out_real}')
    

if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print('Total time: ', round(end - start, 3))
