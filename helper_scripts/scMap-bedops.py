import time
import subprocess
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import tarfile
import gzip
import io

parser = argparse.ArgumentParser(prog='scMap', description='Map snmC-seq to bulk loci. Requires BEDOPS (https://bedops.readthedocs.io) to be installed separately.')
parser.add_argument('array', help='Reference array (e.g., Illumina 450K) in BED format (without headers or column names, containing only chrom, chromStart, chromEnd, name, and strand)')
parser.add_argument('singlecell', help='Path to folder/tarball containing the single cell data')
parser.add_argument('-o', '--out', default='./combined.bed', help='Name and path of output file')
parser.add_argument('-s', '--strand', action='store_true', help='Enforce strand specificity')
parser.add_argument('-t', '--tmp', default='/tmp/', help='Directory to store temp files in. Default = /tmp/')
parser.add_argument('-m', '--mc_classes', default=None, help='Count the frequency of mc-classes for all positions, and store to the results to the specified output file.')
parser.add_argument('-i', '--interval', default=0, type=int, help='Interval to use for the single cell data for mapping.')
parser.add_argument('-c', '--coverage', default=None, help='Count the total number of positions in a cell, then append the result to the specified output file.')
args = parser.parse_args()


# Read each chromosome file and return a df with the required columns
def process_df(file, array_df: pd.DataFrame, interval: int, strand: bool, mc_classes: bool = False, coverage: bool = False) -> pd.DataFrame:
    df = pd.read_csv(file, sep='\t', na_filter=False, dtype={'chr': str})
    
    if coverage:
        global line_count
        line_count += df.shape[0]
    
    df['chr'] = 'chr' + df['chr']
    
    if interval == 0:
        if strand:
            df = pd.merge(df, array_df, how='inner', left_on=['chr', 'pos', 'strand'], right_on=['chrom', 'chromStart', 'strand2'])
        else:
            df = pd.merge(df, array_df, how='inner', left_on=['chr', 'pos'], right_on=['chrom', 'chromStart'])

    df['pos2'] = df['pos'] + interval
    df['frac'] = df['mc_count'] / df['total']
    
    # Get mc_classes frequency after filtering 
    if mc_classes:
        global mc_counts
        mc_counts = pd.concat([mc_counts, df['mc_class']])
    
    if interval == 0:
        return df[['chr', 'pos', 'pos2', 'name', 'strand', 'frac']]
    else:
        if strand:
            return df[['chr', 'pos', 'pos2', 'frac', 'strand']]
        else:
            return df[['chr', 'pos', 'pos2', 'frac']]


# Open chromosome tsv file and append to outpath
def chromosome_append(chromosome_file, outpath, array_df: pd.DataFrame, interval: int = 1, strand: bool = True, mc_classes: bool = False, coverage: bool = False) -> None:
    try:
        with gzip.open(chromosome_file, 'rb') as chromosome:
            df = process_df(file=chromosome, interval=interval, strand=strand, array_df=array_df, mc_classes=mc_classes, coverage=coverage)
    except: 
        # If just a tsv file
        df = process_df(file=chromosome_file, interval=interval, strand=strand, array_df=array_df, mc_classes=mc_classes, coverage=coverage)
    df.to_csv(outpath, sep='\t', mode='a', header=False, index=False)
    

def handle_tarfile(tfile, outpath, array_df, interval: int, strand: bool, fileobj: bool = False, mc_classes: bool = False, coverage: bool = False) -> None:
    # Read tarfile from file or memory
    if fileobj:
        tar = tarfile.open(fileobj=io.BytesIO(tfile))
    else:
        tar = tarfile.open(tfile)

    for member in tar:
        if member.isdir():
            # Skip directories since python gets all contained files.
            continue
        mem = tar.extractfile(member)
        if tarfile.is_tarfile(mem):
            # Recursively open tarfiles
            file_content_byte = tar.extractfile(member.name).read()
            handle_tarfile(tfile=file_content_byte, outpath=outpath, array_df=array_df, interval=interval, strand=strand, fileobj=True, mc_classes=mc_classes, coverage=coverage)
        else:
            chromosome_append(chromosome_file=mem, outpath=outpath, array_df=array_df, interval=interval, strand=strand, mc_classes=mc_classes, coverage=coverage)
    tar.close()
    

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
        count_out = Path(args.mc_classes).expanduser()
        # Create variable to store all mc_classes
        global mc_counts
        mc_counts = pd.Series([], dtype=str)
    
    coverage = args.coverage != None
    if coverage:
        # Create variable to store total line count
        global line_count
        line_count = 0
        coverage_out = Path(args.coverage).expanduser()
        if not coverage_out.exists():
            with open(coverage_out, 'w') as f:
                f.write('cell\tcoverage\n')
    
    # Temp files
    scUnsorted = tmp/f'{path.name}.scUnsorted.tmp'
    scSorted = tmp/f'{path.name}.scSorted.tmp'
    arraySorted = tmp/f'{path.name}.arraySorted.tmp'
    bmapOut = tmp/f'{path.name}.bmapOut.tmp'
    scStrand = tmp/f'{path.name}.scStrand.tmp'
    strandCol = tmp/f'{path.name}.strandCol.tmp'
    
    # Create empty file
    with open(scUnsorted, 'w') as f:
        f.write('')

    # If we have a pool_x folder containing one cell's 23 chromosome.tsv.gz
    if path.is_dir():
        for chromosome in path.iterdir():
            if not chromosome.is_dir():
                chromosome_append(chromosome_file=chromosome, outpath=scUnsorted, array_df=array_df, interval=args.interval, strand=args.strand, mc_classes=mc_classes, coverage=coverage)

    # If we have a tar or tar.gz file containing the pool_x folder or other tar files
    elif tarfile.is_tarfile(path):
        handle_tarfile(tfile=path, outpath=scUnsorted, array_df=array_df, interval=args.interval, strand=args.strand, mc_classes=mc_classes, coverage=coverage)

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
            #rows = pd.read_csv(scUnsorted, sep='\t', header=None).shape[0] # Count using pandas
            #wc = subprocess.run(['wc', '-l', str(scUnsorted)], capture_output=True) # Count lines using wc -l
            #num_lines = wc.stdout.decode().split(' ')[0]
            outfile.write(f'{path.stem}\t{line_count}\n')
    
    
    # Process single_cell data to tsv file containing mapped locations.
    if args.interval == 0:
        print('Mapped single-cell data to reference array')
        pd.read_csv(scUnsorted, na_filter=False, sep='\t', header=None).sort_values(by=[0, 1]).to_csv(out_real, sep='\t', header=False, index=False)
        scUnsorted.unlink() # Remove temp files
        print(f'Saved output to {out_real}')
    
    else:
        # Sort 450K file
        with open(arraySorted, 'w') as outfile:
            subprocess.run(['sort-bed', str(array)], stdout=outfile)

        # Sort single cell file
        with open(scSorted, 'w') as outfile:
            subprocess.run(['sort-bed', str(scUnsorted)], stdout=outfile)
        scUnsorted.unlink()

        # Run bedmap to match chromosome location to fractions
        # Contains chr, chrStart, chrEnd, name, strand, frac
        subprocess.run(f'bedmap --echo --echo-map-id --delim "\t" {arraySorted} {scSorted} > {bmapOut}', shell=True)

        # Match location to strand
        if args.strand:
            with open(scStrand, 'w') as outfile:
                subprocess.run(['cut', '-f', '1-3,5', str(scSorted)], stdout=outfile) 
                # Bed file with only chrom + strand info

            subprocess.run(f'bedmap --echo --echo-map-id --delim "\t" {arraySorted} {scStrand} | cut -f 4,6 > {strandCol}', shell=True) 
            # 2 col file with name and strandMatches

            scStrand.unlink()

        print('Mapped single-cell data to reference array')
        scSorted.unlink()
        arraySorted.unlink()
        
        # Remove na rows and calculate averages
        if args.strand:
            mapped = process_bedmap(bed_frac_output = bmapOut, bed_strand_output = strandCol, strand=True)
            strandCol.unlink()
        else:
            mapped = process_bedmap(bed_frac_output = bmapOut)

        mapped.to_csv(out_real, sep='\t', header=False, index=False)
        print(f'Saved output to {out_real}')
        bmapOut.unlink()
    

if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print('Total time: ', round(end - start, 3))
