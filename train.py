"""
scaden Main functionality

Contains code to
- process a training datasets
- train a model
- perform predictions

"""

# Imports
import logging
import pandas as pd
import numpy as np
from class_scaden import Scaden


logger = logging.getLogger(__name__)

"""
PARAMETERS
"""
# ==========================================#
architectures = {'m256':    ([256, 128, 64, 32],    [0, 0, 0, 0]),
                 'm512':    ([512, 256, 128, 64],   [0, 0.3, 0.2, 0.1]),
                 'm1024':   ([1024, 512, 256, 128], [0, 0.6, 0.3, 0.1])}

# Extract architectures
M256_HIDDEN_UNITS = architectures["m256"][0]
M512_HIDDEN_UNITS = architectures["m512"][0]
M1024_HIDDEN_UNITS = architectures["m1024"][0]
M256_DO_RATES = architectures["m256"][1]
M512_DO_RATES = architectures["m512"][1]
M1024_DO_RATES = architectures["m1024"][1]

# ==========================================#


def training(
    data_path, train_datasets, model_dir, batch_size, learning_rate, num_steps, seed=0, loss_values=None, loss_curve=None
):
    """
    Perform training of three a scaden model ensemble consisting of three different models
    :param model_dir:
    :param batch_size:
    :param learning_rate:
    :param num_steps:
    :return:
    """
    # Convert training datasets
    if train_datasets == "":
        train_datasets = []
    else:
        train_datasets = train_datasets.split(",")
        logger.info(f"Training on: [cyan]{train_datasets}")
    
    # Training of M256 model
    logger.info("[cyan]Training M256 Model ... [/]")
    cdn256 = Scaden(
        model_dir=model_dir + "/m256",
        model_name="m256",
        batch_size=batch_size,
        learning_rate=learning_rate,
        num_steps=num_steps,
        seed=seed,
        hidden_units=M256_HIDDEN_UNITS,
        do_rates=M256_DO_RATES,
    )
    cdn256.train(input_path=data_path, train_datasets=train_datasets)
    loss256 = cdn256.loss_curve
    del cdn256

    # Training of M512 model
    logger.info("[cyan]Training M512 Model ... [/]")
    cdn512 = Scaden(
        model_dir=model_dir + "/m512",
        model_name="m512",
        batch_size=batch_size,
        learning_rate=learning_rate,
        num_steps=num_steps,
        seed=seed,
        hidden_units=M512_HIDDEN_UNITS,
        do_rates=M512_DO_RATES,
    )
    cdn512.train(input_path=data_path, train_datasets=train_datasets)
    loss512 = cdn512.loss_curve
    del cdn512

    # Training of M1024 model
    logger.info("[cyan]Training M1024 Model ... [/]")
    cdn1024 = Scaden(
        model_dir=model_dir + "/m1024",
        model_name="m1024",
        batch_size=batch_size,
        learning_rate=learning_rate,
        num_steps=num_steps,
        seed=seed,
        hidden_units=M1024_HIDDEN_UNITS,
        do_rates=M1024_DO_RATES,
    )
    cdn1024.train(input_path=data_path, train_datasets=train_datasets)
    loss1024 = cdn1024.loss_curve
    del cdn1024
    
    logger.info("[green]Training finished.")
    
    loss = pd.DataFrame({'step': list(range(num_steps))*3, 
                         'model': np.repeat(['m256', 'm512', 'm1024'], num_steps),
                         'loss': np.concatenate((loss256, loss512, loss1024))
                        })
    if loss_values is not None:
        loss.to_csv(loss_values, index=False, sep='\t')
        logger.info(f'Loss values written to: {loss_values}')
    
    if loss_curve is not None:
        import seaborn as sns
        g = sns.lineplot(data=loss, x='step', y='loss', hue='model')
        g.figure.savefig(loss_curve)
        logger.info(f'Loss curve saved to: {loss_curve}')
        
        