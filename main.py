#!/usr/bin/env python3
import click

@click.command()

# Pipeline controls and logging
@click.option('--load', help='Load parameters from YAML file.')

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


def cli(load, no_sim, no_proc, no_pred, no_eval, config, reference,
        out, data, cells, n_samples, pattern, unknown, prefix, data_format, 
        pred, processed_path, var_cutoff, scaling, 
        train_datasets, model_dir, batch_size, learning_rate, steps, seed, 
        prediction_outname, prediction_scaling, 
        ground_truth):
    
    # Move args to namespace to allow overwriting variables from strings
    from types import SimpleNamespace
    a = SimpleNamespace(**locals())
    
    # Overwrite args with yaml file if provided
    if load:
        import yaml
        with open(load, 'r') as fname:
            params = yaml.safe_load(fname)
        for p in params:
            if p in locals().keys():           
                setattr(a, p, params[p])
            else:
                raise ValueError(f"Unknown key '{p}' in YAML file. Keys name must be the same as the long form parameters.")
    
    # Run scaden commands
    import sys
    import os
    import time
    import json
    import logging
    import rich.logging
    
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
    training_data = a.out + a.prefix + ".h5ad" # h5ad file generated from scaden simulate

    if a.processed_path is None:
        a.processed_path = a.out + "processed.h5ad"

    if a.prediction_outname is None:
        a.prediction_outname = a.out + "scaden_predictions.txt"
    
    # Simulate samples
    if a.no_sim != True:
        from simulate import simulation
        # Create artificial bulk RNA-seq data from scRNA-seq dataset(s)
        simulation(
            simulate_dir = a.out,
            data_dir = a.data,
            sample_size = a.cells,
            num_samples = a.n_samples,
            pattern = a.pattern,
            unknown_celltypes = a.unknown,
            out_prefix = a.prefix,
            fmt = a.data_format
        )

    # Preprocess training data
    if a.no_proc != True:
        from process import processing
        processing(
            testing_data = a.pred, 
            training_data = training_data, 
            processed_path = a.processed_path,
            cells_per_sample = a.cells,
            scaling = a.scaling,
            var_cutoff = a.var_cutoff
        )
    
    # Train model
    from train import training
    start = time.time()
    training(
        data_path = a.processed_path,
        train_datasets = a.train_datasets,
        model_dir = a.model_dir,
        batch_size = a.batch_size,
        learning_rate = a.learning_rate,
        num_steps = a.steps,
        seed = a.seed,
    )
    end = time.time()
    training_time = round(end - start, 3)
    logger.info(f"[bold]Training time: [green]{training_time}[/]")
    
    # Create config.json file
    config_filepath = a.out + f'config_{a.config}.json'
    with open(config_filepath, 'w', encoding='utf-8') as config_file:
        config_json = {
            "config": a.config, 
            "Data": {
                "cells": a.cells, 
                "n_samples": a.n_samples, 
                "var_cutoff": a.var_cutoff,
                "scaling": a.scaling,
                "reference": a.reference
            }, 
            "Model": {
                "seed": a.seed, 
                "steps": a.steps, 
                "batch_size": a.batch_size, 
                "learning_rate": a.learning_rate
            }
        }
        json.dump(config_json, config_file, ensure_ascii=False, indent=4)
    
    logger.info(f"[bold]Created config file: [green]{config_filepath}[/]")

    # Make predictions
    if a.no_pred != True:
        from predict import prediction
        prediction(
            model_dir = a.model_dir, 
            data_path = a.pred, 
            out_name = a.prediction_outname, 
            cells = a.cells, 
            scaling = a.prediction_scaling, 
            seed = a.seed
        )
        
    # Evaluate predictions
    if a.no_eval != True:
        from evaluate import evaluate
        report_filepath = a.out + f'report_{a.config}.json'
        evaluate(
            predictions_filename = a.prediction_outname, 
            ground_truth_filename = a.ground_truth,
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

