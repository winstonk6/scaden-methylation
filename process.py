import numpy as np
import pandas as pd
from anndata import read_h5ad
from sklearn.preprocessing import MinMaxScaler
import logging

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
    try:
        data = pd.read_table(testing_data, index_col=0)
    except UnicodeDecodeError:
        data = read_h5ad(testing_data).to_df()
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

    logger.info("Scaling using " + str(scaling))
    if scaling == "log" or scaling == "log_min_max":
        training.X = sample_scaling(training.X, scaling)
    elif scaling == "frac" or scaling == "fraction":  # frac_notna scaling done in bulk_simulator
        training.X = training.X / num_cells

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

