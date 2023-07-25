import glob
import argparse
from pathlib import Path
import pandas as pd
import anndata as ad
import gc


parser = argparse.ArgumentParser(prog='cells_to_matrix', description='Combines scores of BED files for single cell data into a matrix. BED files should have 6 columns: chrom, chromStart, chromEnd, name, strand, and score.')
parser.add_argument('data', help='Directory or Unix style pattern matching containing BED files, in quotes. e.g. "~/data/" or "~/data/*.bed"')
parser.add_argument('-o', '--outfile', default='./cell_matrix.txt', help='Name and location of output file')
parser.add_argument('-f', '--format', default='tsv', help='Output file format for the matrix. Options: csv, tsv, h5ad')
args = parser.parse_args()

out = Path(args.outfile).expanduser()
data = Path(args.data).expanduser()
if data.is_dir():
    files = [file for file in data.iterdir() if not file.is_dir()]
    if len(files) == 0:
        raise ValueError(f"'{data}' does not have any files")
elif '*' in str(data) or '?' in str(data) or '[' in str(data) or ']' in str(data):
    files = [Path(file) for file in glob.glob(str(data)) if not Path(file).is_dir()]
    if len(files) == 0:
        raise ValueError(f"No files match the pattern '{data}'")
else:
    raise ValueError(f"'{data}' is not a directory or a Unix file pattern")

    
cellname0 = files[0].stem
df = pd.read_csv(files[0], sep='\t', header=None, names=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', cellname0])[['name', cellname0]]

if len(files) > 1:
    for file in files[1:]:
        cellname = file.stem
        cell = pd.read_csv(file, sep='\t', header=None, names=['chrom', 'chromStart', 'chromEnd', 'name', 'strand', cellname])[['name', cellname]]
        df = df.merge(cell, how='outer', on='name')
        del cell
        gc.collect()

df.set_index('name', inplace=True)

if args.format == 'tsv':
    df.to_csv(args.outfile, sep='\t')
elif args.format == 'csv':
    df.to_csv(args.outfile)
elif args.format == 'h5ad':
    ad.AnnData(df).write(args.outfile)
