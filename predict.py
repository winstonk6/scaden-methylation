"""
scaden Main functionality

Contains code to
- process a training datasets
- train a model
- perform predictions

"""

# Imports
import logging
#import tensorflow as tf
from anndata import read_h5ad
#from scaden.model.architectures import architectures # CHANGE: add code from architectures to this file.
#from scaden.model.scaden import Scaden
from class_scaden import Scaden # CHANGE: uses modified Scaden class to give option to scale by fraction instead of log

logger = logging.getLogger(__name__)

"""
PARAMETERS
"""
"""
scaden.model.architectures contains the three architectures used for the model ensemble in scaden
The architectures are stored in the dict 'architectures', each element of which is a tuple that consists
of the hidden_unit sizes (element 0) and the respective dropout rates (element 1)
"""

architectures = {'m256':    ([256, 128, 64, 32],    [0, 0, 0, 0]),
                 'm512':    ([512, 256, 128, 64],   [0, 0.3, 0.2, 0.1]),
                 'm1024':   ([1024, 512, 256, 128], [0, 0.6, 0.3, 0.1])}

# ==========================================#

# Extract architectures
M256_HIDDEN_UNITS = architectures["m256"][0]
M512_HIDDEN_UNITS = architectures["m512"][0]
M1024_HIDDEN_UNITS = architectures["m1024"][0]
M256_DO_RATES = architectures["m256"][1]
M512_DO_RATES = architectures["m512"][1]
M1024_DO_RATES = architectures["m1024"][1]

# ==========================================#


def prediction(model_dir, data_path, out_name, cells, scaling, seed=0):
    """
    Perform prediction using a trained scaden ensemble
    :param model_dir: the directory containing the models
    :param data_path: the path to the gene expression file
    :param out_name: name of the output prediction file
    :return:
    """

    # Small model predictions
    cdn256 = Scaden(
        model_dir=model_dir + "/m256",
        model_name="m256",
        seed=seed,
        hidden_units=M256_HIDDEN_UNITS,
        do_rates=M256_DO_RATES,
        cells=cells,
        scaling=scaling,
    )
    # Predict ratios
    preds_256 = cdn256.predict(input_path=data_path)

    # Mid model predictions
    cdn512 = Scaden(
        model_dir=model_dir + "/m512",
        model_name="m512",
        seed=seed,
        hidden_units=M512_HIDDEN_UNITS,
        do_rates=M512_DO_RATES,
        cells=cells,
        scaling=scaling,
    )
    # Predict ratios
    preds_512 = cdn512.predict(input_path=data_path)

    # Large model predictions
    cdn1024 = Scaden(
        model_dir=model_dir + "/m1024",
        model_name="m1024",
        seed=seed,
        hidden_units=M1024_HIDDEN_UNITS,
        do_rates=M1024_DO_RATES,
        cells=cells,
        scaling=scaling,
    )
    # Predict ratios
    preds_1024 = cdn1024.predict(input_path=data_path)

    # Average predictions
    preds = (preds_256 + preds_512 + preds_1024) / 3
    preds.to_csv(out_name, sep="\t")
    logger.info(f"[bold]Created cell composition predictions: [green]{out_name}[/]")
