import numpy as np
import pandas as pd
import anndata
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--prediction_file', help='File containing the bulk counts data that we want to make predictions on.')
parser.add_argument('--training_data', help='h5ad file generated from scaden simulate')
parser.add_argument('--processed_path', help='Name of the file that the processed data will be saved to. Must end with .h5ad')
parser.add_argument('--cells', help='Parameter taken from "scaden simulate --cells ..."')
parser.add_argument('--var_cutoff', default=0.1)

args = parser.parse_args()

data_path = args.prediction_file
training_data = args.training_data
processed_path = args.processed_path
cells_per_sample = int(args.cells)
var_cutoff = float(args.var_cutoff)


def get_signature_genes(input_path, sig_genes_complete, var_cutoff=0.1):
    """
    Get overlap between signature genes and available genes
    :param input_path:
    :param sig_genes_complete:
    :return: new sig_genes
    """
    data = pd.read_table(input_path, index_col=0)
    keep = data.var(axis=1) > var_cutoff
    data = data.loc[keep]
    available_genes = list(data.index)
    new_sig_genes = list(set(available_genes).intersection(sig_genes_complete))
    #n_sig_genes = len(new_sig_genes)
    new_sig_genes.sort() # Change: sort gene names for easier testing
    return new_sig_genes

def preprocess_h5ad_data(raw_input_path, processed_path, num_cells, sig_genes=None):
    """
    Preprocess raw input data for the model
    :param raw_input_path:
    :param group_small:
    :param signature_genes:
    :return:
    """
    raw_input = anndata.read_h5ad(raw_input_path)

    # Select features go use
    raw_input = raw_input[:, sig_genes]

    # CHANGE: no scaling, counts converted to fractions
    #raw_input.X = sample_scaling(raw_input.X, scaling_option)
    raw_input.X = raw_input.X / num_cells
    
    # CHANGE: sort training data ground truth table columns
    raw_input.obs.sort_index(axis=1, inplace=True)

    raw_input.write(processed_path)

    
def processingv2(data_path, training_data, processed_path, cells_per_sample, var_cutoff=0.1):
    """
    Process a training dataset to contain only the genes also available in the prediction data
    :param data_path: path to prediction data
    :param training_data: path to training data (h5ad file)
    :param processed_path: name of processed file
    :return:
    """   
    # Get the common genes (signature genes)
    raw_input = anndata.read_h5ad(training_data)
    sig_genes_complete = list(raw_input.var_names)
    sig_genes = get_signature_genes(input_path=data_path, sig_genes_complete=sig_genes_complete, var_cutoff=var_cutoff)
    
    # Pre-process data with new signature genes
    preprocess_h5ad_data(raw_input_path=training_data,
                         processed_path=processed_path,
                         num_cells=cells_per_sample,
                         sig_genes=sig_genes)
    
processingv2(data_path, training_data, processed_path, cells_per_sample, var_cutoff)