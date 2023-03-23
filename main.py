#!/usr/bin/env python3
import click

@click.command()

# Pipeline controls and logging
@click.option('--no_sim', is_flag=True, help='Skip the creation of simulated training samples if you have already created a training set.')

@click.option('--no_proc', is_flag=True, 
                    help='Skip the preprocessing step if you have already processed the training data. If using this flag, then --proc is required.')

@click.option('--no_pred', is_flag=True, help='Skip the creation of simulated training samples if you have already created a training set.')

@click.option('--no_eval', is_flag=True, help='Skip the evaluation of the model predictions.')

@click.option('--config', default='test', help='Name of configuration.')

@click.option("--reference", "-r", default='1', help="Name of the scMethyl dataset.")

# scaden simulate
@click.option("--out", "-o", default="./", help="Directory to store output files in")

@click.option("--data", "-d", default=".", help="Path to scRNA-seq dataset(s)")

@click.option("--cells", "-c", default=100, type=int, help="Number of cells per sample [default: 100]")

@click.option("--n_samples", "-n", default=1000, type=int, help="Number of samples to simulate [default: 1000]")

@click.option("--pattern", default="*_counts.txt", help="File pattern to recognize your processed scRNA-seq count files",)

@click.option("--unknown", "-u", multiple=True, default=["unknown"], 
                    help="Specifiy cell types to merge into the unknown category. Specify this flag for every cell type you want "
                    "to merge in unknown. [default: unknown]")

@click.option("--prefix", "-p", default="data", help="Prefix to append to training .h5ad file [default: data]")

@click.option("--data_format", "-f", default="txt", help="Data format of scRNA-seq data, can be 'txt' or 'h5ad' [default: 'txt']")


# scaden process
@click.option('--pred', help='Bulk data file (i.e. testing set) that we want to perform deconvolution on.')

@click.option('-proc', '--processed_path', help='Name of the file that the processed data will be saved to. Must end with .h5ad')

@click.option('--var_cutoff', default=0.1, type=float,
                    help='Filter out genes with a variance less than the specified cutoff. A low value is recommended, '
                    'this should only remove genes that are obviously uninformative. [default: 0.1]')

@click.option('-sc', '--scaling', default='fraction', 
                    help='Change scaling option for preprocessing the training data. If something other than the provided options is '
                    'used, then no scaling will be done. [default: fraction] '
                    'Options: None (No scaling), log / log_min_max (log2, then scale to the range 0,1), '
                    'frac / fraction (Divide values by the number of cells)')


# scaden train
@click.option("--train_datasets", default="", help="Comma-separated list of datasets used for training. Uses all by default.")

@click.option("-m", "--model_dir", default="./", help="Path to store the model in")

@click.option("-b", "--batch_size", default=128, type=int, help="Batch size to use for training. [default: 128]")

@click.option("-l", "--learning_rate", default=0.0001, type=float, help="Learning rate used for training. [default: 0.0001]")

@click.option("--steps", default=5000, type=int, help="Number of training steps. [default: 5000]")

@click.option("--seed", default=0, type=int, help="Set random seed")


# scaden predict
@click.option('-po', '--prediction_outname', default="scaden_predictions.txt", help="Name of predictions file.")

@click.option('--prediction_scaling', default="fraction", help="Change scaling option for the preprocessing done when making predictions. Uses the same options as --scaling.")

# Evaluate scaden predictions
@click.option('-g', '--ground_truth', help="Name of file containing the ground truth cell proportions.")
        

def cli(config, no_sim, no_proc, no_pred, no_eval, 
         out, data, cells, n_samples, pattern, unknown, prefix, data_format, 
         pred, processed_path, var_cutoff, scaling, 
         train_datasets, model_dir, batch_size, learning_rate, steps, seed, 
         prediction_outname, prediction_scaling, 
         ground_truth, reference):
    """
    Run all scaden commands in a single program
    """
    import time
    import json

    import numpy as np
    import pandas as pd
    from anndata import read_h5ad
    from sklearn.preprocessing import MinMaxScaler
    
    import sys
    import logging
    import rich.logging
    import os

    os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.addHandler(
        rich.logging.RichHandler(
            level=logging.INFO,
            console=rich.console.Console(file=sys.stderr),
            show_time=False,
            markup=True,
        )
    )
    
    # Create paths if not specified by user
    training_data = out + prefix + ".h5ad" # h5ad file generated from scaden simulate

    if processed_path is None:
        processed_path = out + "processed.h5ad"

    if prediction_outname is None:
        prediction_outname = out + "scaden_predictions.txt"
    
    # Simulate samples
    if no_sim != True:
        from simulate import simulation
        # Create artificial bulk RNA-seq data from scRNA-seq dataset(s)
        simulation(
            simulate_dir = out,
            data_dir = data,
            sample_size = cells,
            num_samples = n_samples,
            pattern = pattern,
            unknown_celltypes = unknown,
            out_prefix = prefix,
            fmt = data_format
        )

    # Preprocess training data
    if no_proc != True:
        from process import processing
        processing(
            testing_data = pred, 
            training_data = training_data, 
            processed_path = processed_path,
            cells_per_sample = cells,
            scaling = scaling,
            var_cutoff = var_cutoff
        )
    
    # Train model
    from train import training
    start = time.time()
    training(
        data_path = processed_path,
        train_datasets = train_datasets,
        model_dir = model_dir,
        batch_size = batch_size,
        learning_rate = learning_rate,
        num_steps = steps,
        seed = seed,
    )
    end = time.time()
    training_time = round(end - start, 3)
    logger.info(f"[bold]Training time: [green]{training_time}[/]")
    
    # Create config.json file
    config_filepath = out + f'config_{config}.json'
    with open(config_filepath, 'w', encoding='utf-8') as config_file:
        config_json = {
            "config": config, 
            "Data": {
                "cells": cells, 
                "n_samples": n_samples, 
                "var_cutoff": var_cutoff,
                "scaling": scaling,
                "reference": reference
            }, 
            "Model": {
                "seed": seed, 
                "steps": steps, 
                "batch_size": batch_size, 
                "learning_rate": learning_rate
            }
        }
        json.dump(config_json, config_file, ensure_ascii=False, indent=4)
    
    logger.info(f"[bold]Created config file: [green]{config_filepath}[/]")

    # Make predictions
    if no_pred != True:
        from predict import prediction
        prediction(
            model_dir = model_dir, 
            data_path = pred, 
            out_name = prediction_outname, 
            cells = cells, 
            scaling = prediction_scaling, 
            seed = seed
        )
        
    # Evaluate predictions
    if no_eval != True:
        from evaluate import evaluate
        report_filepath = out + f'report_{config}.json'
        evaluate(
            predictions_filename = prediction_outname, 
            ground_truth_filename = ground_truth,
            report_filepath = report_filepath,
            training_time = training_time
        )

def main():
    text = """
     ____                _            
    / ___|  ___ __ _  __| | ___ _ __  
    \___ \ / __/ _` |/ _` |/ _ \ '_ \ 
     ___) | (_| (_| | (_| |  __/ | | |
    |____/ \___\__,_|\__,_|\___|_| |_|
    """
    click.echo(click.style(text, fg="blue"))
    cli()
        
if __name__ == "__main__":
    main()


