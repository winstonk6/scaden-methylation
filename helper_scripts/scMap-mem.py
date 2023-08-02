#!/usr/bin/env python3
import time
import subprocess
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import tarfile
import gzip
import io

parser = argparse.ArgumentParser(prog='scMap-mem', description='Map snmC-seq to bulk loci. Only does exact matching.')
parser.add_argument('array', help='Reference array (e.g., Illumina 450K) in BED format (without headers or column names, containing only chrom, chromStart, chromEnd, name, and strand)')
parser.add_argument('singlecell', help='Path to folder/tarball containing the single cell data')
parser.add_argument('-o', '--out', default='./combined.bed', help='Name and path of output file')
parser.add_argument('-s', '--strand', action='store_true', help='Enforce strand specificity')
parser.add_argument('-t', '--tmp', default='/tmp/', help='Directory to store temp files in. Default = /tmp/')
parser.add_argument('-m', '--mc_classes', default=None, help='Count the frequency of mc-classes for all positions, and store to the results to the specified output file.')
parser.add_argument('-c', '--coverage', default=None, help='Count the total number of positions in a cell, then append the result to the specified output file.')
parser.add_argument('-a', '--add_header', action='store_true', help='Use if the chromosome files do not have a row for the column names (chr, pos, strand, mc_class, mc_count, total, methylated)')
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
    
    

def handle_tarfile(tfile, array_df: pd.DataFrame, sc_unsorted_df: pd.DataFrame, 
                   strand: bool, mc_classes: bool = False, coverage: bool = False, add_header: bool = False, fileobj: bool = False) -> pd.DataFrame:
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
            file_content_byte = tar.extractfile(member.name).read()
            df = handle_tarfile(tfile=file_content_byte, array_df=array_df, sc_unsorted_df=df, 
                                strand=strand, mc_classes=mc_classes, coverage=coverage, add_header=add_header, fileobj=True)
        else:
            df = chromosome_append(chromosome_file=mem, array_df=array_df, sc_unsorted_df=df, 
                                   strand=strand, mc_classes=mc_classes, coverage=coverage, add_header=add_header)
    tar.close()
    return df
    

def score_strand_calc(strand_frac_df: pd.DataFrame) -> pd.Series:
    """
    Calculate scores when enforcing strand specificity.
    
    Args:
        strand_frac_df: pd.DataFrame
            Data frame with 3 columns: strand, scStrand, frac
    
    Returns: pd.Series
        Contains scores converted to float. Rows where the strand doesn't match are NA values.
    """
    df = strand_frac_df.reset_index()
    scores = np.empty_like(df['frac'])
    for n in df.index:
        # If only one value and strands match
        if df.at[n,'strand'] == df.at[n,'scStrand']:
            scores[n] = float(df.at[n,'frac'])
            
        # If multiple values
        elif df.at[n,'strand'] in df.at[n,'scStrand']:
            strands = df.at[n,'scStrand'].split(';')
            score_list = df.at[n,'frac'].split(';')
            
            # Calculate score only for the matching strands
            match_scores = []
            for s in range(len(strands)):
                if df.at[n,'strand'] == strands[s]:
                    match_scores.append(float(score_list[s]))
            scores[n] = sum(match_scores)/len(match_scores)

        # If strand doesn't match: NA
        else:
            scores[n] = np.nan
    
    return pd.Series(scores, index=df['index'])


def frac_avg(val: str) -> float:
    """
    Calculate a score. Used in process_bedmap when strand specificity is not enforced
    """
    if ';' in val:
        vals = val.split(';')
        avg = sum([float(n) for n in vals])/len(vals)
    else:
        avg = float(val)
    return avg


def process_bedmap(bed_frac_output: str, bed_strand_output: str = None, strand: bool = False) -> pd.DataFrame:
    """
    Get the final data frame containing chrom, chromStart, chromEnd, cpg, strand, score.
    
    Args:
        bed_frac_output: str
            Location of the bed file with columns [chrom, chromStart, chromEnd, name, strand, score]
        
        bed_strand_output: str
            Location of the bed file with columns 
    """
    # Load files
    df = pd.read_csv(bed_frac_output, sep='\t', header=None, names=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'frac'], dtype={'frac': str}).dropna()
    if strand:
        strandMatch = pd.read_csv(bed_strand_output, header=None, sep='\t', names=['name', 'scStrand']).dropna()
        df['scStrand'] = strandMatch['scStrand']
        del strandMatch
        # Process frac column into score column containing single values + NA values (representing strand mismatch)
        df['score'] = score_strand_calc(df[['strand', 'scStrand', 'frac']])
    
    else:
        df['score'] = df['frac'].apply(frac_avg)
    
    return df[['chrom', 'chromStart', 'chromEnd', 'name', 'strand', 'score']].dropna()
    
    
def main():
    tmp = Path(args.tmp).expanduser()
    path = Path(args.singlecell).expanduser()
    array = args.array
    array_df = pd.read_csv(array, sep='\t', header=None, names=['chrom', 'chromStart', 'chromEnd', 'name', 'strand2'])
    out_real = Path(args.out).expanduser()

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

    # If we have a pool_x folder containing one cell's 23 chromosome.tsv.gz
    if path.is_dir():
        for chromosome in path.iterdir():
            if not chromosome.is_dir():
                scUnsorted = chromosome_append(chromosome_file=chromosome, array_df=array_df, sc_unsorted_df=scUnsorted,
                                               strand=args.strand, mc_classes=mc_classes, coverage=coverage, add_header=args.add_header)

    # If we have a tar or tar.gz file containing the pool_x folder or other tar files
    elif tarfile.is_tarfile(path):
        scUnsorted = handle_tarfile(tfile=path, array_df=array_df, sc_unsorted_df=scUnsorted,
                                    strand=args.strand, mc_classes=mc_classes, coverage=coverage, add_header=args.add_header)

    # If we have a tsv file containing the data for all chromosomes
    elif path.is_file():
        scUnsorted = chromosome_append(chromosome_file=path, array_df=array_df, sc_unsorted_df=scUnsorted, 
                                       strand=args.strand, mc_classes=mc_classes, coverage=coverage, add_header=args.add_header)
    
    else:
        raise TypeError(f"'{args.singlecell}' is not a directory or a tarball.")
    
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

