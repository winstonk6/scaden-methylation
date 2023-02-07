import numpy as np
import pandas as pd
from anndata import read_h5ad
#import argparse
from sklearn.preprocessing import MinMaxScaler
import logging

"""
parser = argparse.ArgumentParser()
parser.add_argument('--prediction_file', help='File containing the bulk counts data (i.e. testing set) that we want to make predictions on.')
parser.add_argument('--training_data', help='h5ad file generated from scaden simulate (containing the training data)')
parser.add_argument('--processed_path', help='Name of the file that the processed data will be saved to. Must end with .h5ad')
parser.add_argument('--cells', help='Parameter taken from "scaden simulate --cells ..."')
parser.add_argument('-sc', '--scaling', default='fraction', help='Option to use for scaling. Default = Fraction')
parser.add_argument('--var_cutoff', default=0.1)

args = parser.parse_args()

testing_data = args.prediction_file # testing data
training_data = args.training_data
processed_path = args.processed_path
cells_per_sample = int(args.cells)
var_cutoff = float(args.var_cutoff)
scaling = args.scaling

"""
logger = logging.getLogger(__name__)

def sample_scaling(x, scaling_option):
    """
    Apply scaling of data
    :param x:
    :param scaling_option:
    :return:
    """
    if scaling_option:
        # Bring in log space
        x = np.log2(x + 1)

        # Normalize data
        mms = MinMaxScaler(feature_range=(0, 1), copy=True)

        # it scales features so transpose is needed
        x = mms.fit_transform(x.T).T

    return x

def get_signature_genes(testing_data, sig_genes_complete, var_cutoff=0.1):
    """
    Get overlap between signature genes (training) and available genes (testing)
    :param testing_data:
    :param sig_genes_complete:
    :return: new sig_genes
    """
    data = pd.read_table(testing_data, index_col=0)
    keep = data.var(axis=1) > var_cutoff
    data = data.loc[keep]
    available_genes = list(data.index)
          
    new_sig_genes = list(set(available_genes).intersection(sig_genes_complete))
    new_sig_genes.sort() # Change: sort gene names for easier testing
    
    if len(new_sig_genes) == 0:
        raise ValueError(f'No genes found in both the training data and the prediction file that have a variance above {var_cutoff}.\nMake sure that the data sets have genes in common, then try lowering the value of --var_cutoff.')
    
    return new_sig_genes

def preprocess_h5ad_data(training_data, processed_path, num_cells, scaling, sig_genes=None):
    """
    Preprocess raw input data for the model
    :param training_data:
    :param group_small:
    :param signature_genes:
    :return:
    """
    logger.info("Pre-processing raw data ...")
    training = read_h5ad(training_data)

    # Select features go use
    logger.info("Subsetting genes ...")
    training = training[:, sig_genes]

    # CHANGE: add scaling option fraction which converts counts to fractions
    logger.info("Scaling using " + str(scaling))
    if scaling == "log" or scaling == "log_min_max":
        training.X = sample_scaling(training.X, scaling)
    elif scaling == "frac" or scaling == "fraction":
        training.X = training.X / num_cells
    
    # CHANGE: sort training data ground truth table columns
    training.obs.sort_index(axis=1, inplace=True)

    logger.info("Writing to disk ...")
    training.write(processed_path)
    logger.info("Data pre-processing done.")
    logger.info(f"Created processed file: [cyan]{processed_path}[/]")

    
def processing(testing_data, training_data, processed_path, cells_per_sample, scaling, var_cutoff=0.1):
    """
    Process a training dataset to contain only the genes also available in the prediction data
    :param testing_data: path to prediction data (TESTING)
    :param training_data: path to training data (h5ad file)
    :param processed_path: name of processed file
    :return:
    """   
    # Get the common genes (signature genes)
    training = read_h5ad(training_data)
    sig_genes_complete = list(training.var_names)
    sig_genes = get_signature_genes(testing_data=testing_data, sig_genes_complete=sig_genes_complete, var_cutoff=var_cutoff)
    
    # Pre-process data with new signature genes
    preprocess_h5ad_data(training_data=training_data,
                         processed_path=processed_path,
                         num_cells=cells_per_sample,
                         scaling=scaling,
                         sig_genes=sig_genes)
    
#processingv2(testing_data, training_data, processed_path, cells_per_sample, scaling, var_cutoff)