import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
import json
import logging

logger = logging.getLogger(__name__)


# Determine whether a pd.Series has no NA values. # If an NA value is found, return false
def df_no_na(series):
    isna = pd.isna(series)
    for n in isna:
        if n == True:
            return False
    return True


# Remove rows containing at least one NA value from a data frame.
def get_non_na_rows(df):
    return df[df.apply(df_no_na, axis=1)]


# Calculate correlation coefficient between each column of both data frames
def correlation_matrix(pred, actual):
    pred_mat = pred.to_numpy()
    actual_mat = actual.to_numpy()
    notna = pd.notna(actual_mat)
    
    corrcoefmatrix = np.empty((len(pred.columns), len(actual.columns)))
    for i in range(len(pred.columns)):
        for j in range(len(actual.columns)):
            ind = notna[:,j]
            if sum(ind) == 0:
                corrcoefmatrix[i, j] = np.nan
            else:
                corrcoefmatrix[i, j] = np.corrcoef(pred_mat[:,i][ind], actual_mat[:,j][ind])[0,1]

    corr_mat = pd.DataFrame(corrcoefmatrix, index=pred.columns, columns=actual.columns)
    return corr_mat


# Get max value index for each column, along with the value, and get MSE per pairing
def corr_max_table(corr_mat, pred, actual):
    notna = pd.notna(actual)
    idx = corr_mat.idxmax()
    mse = []
    for n in range(len(idx)):
        if pd.isna(idx[n]) == True:
            mse.append(np.nan)
        else:
            p_col = f'{idx[n]}'
            a_col = f'{idx.index[n]}'
            ind = notna[a_col]            
            mse.append(mean_squared_error(pred[p_col][ind], actual[a_col][ind]))

    table = pd.DataFrame(data = {
                             'Predictions best correlated celltype': idx, 
                             'Correlation coefficient': corr_mat.max(axis=0),
                             'Mean squared error': mse
                         }, 
                         index=idx.index)
    return table


def corr_max_abs_table(corr_mat, pred, actual):
    # Get max absolute value index for each column, the max abs value, and the MSE
    idx = corr_mat.apply(np.abs).idxmax()
    try:
        mse = [mean_squared_error(pred[f'{idx[n]}'], actual[f'{idx.index[n]}']) for n in range(len(idx))]
    except:
        raise ValueError("NA values found in ")
    table = pd.DataFrame(data = {
                             'Highest correlated celltype': idx, 
                             'Correlation coefficient': [corr_mat.at[idx[n], corr_mat.columns[n]] for n in range(len(idx))], 
                             'Mean squared error': mse,
                         }, 
                         index=idx.index)
    return table


# When calculating the correlations, first remove the NA values corresponding samples and figure out the index first (c.csv)
def evaluate(predictions_filename, ground_truth_filename, report_filepath, training_time):
    actual = pd.read_table(ground_truth_filename, index_col=0)
    actual.sort_index(axis=1, inplace=True)
    
    pred = pd.read_csv(predictions_filename, sep='\t', index_col=0)
    pred.sort_index(axis=1, inplace=True)
    
    # Normalize values by the sum of the row
    actualsum = actual.apply(np.sum, axis=1)
    actual = actual / np.array([[n]*len(actual.columns) for n in actualsum])
    
    pred = get_non_na_rows(pred)
    actual = actual.loc[pred.index]
    
    logger.info("Evaluating predictions ...")
    corr_mat = correlation_matrix(pred=pred, actual=actual)
    corr_mat_max = corr_max_table(corr_mat=corr_mat, pred=pred, actual=actual)
    #corr_mat_max_abs = corr_max_abs_table(corr_mat=corr_mat, pred=pred, actual=actual)
    
    results = {"training_time": training_time, "results": corr_mat_max.to_dict()}
    
    with open(report_filepath, 'w', encoding='utf-8') as report_file:
        json.dump(results, report_file, ensure_ascii=False, indent=4)
    
    logger.info(f"[bold]Created evaluation results: [green]{report_filepath}[/]")
    return