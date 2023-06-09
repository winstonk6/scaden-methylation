#!/usr/bin/env python3
import click
import sys
import os
import json
import logging
import rich.logging
from pathlib import Path

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

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

training_time = 0


def sc_simulate(out, data, cells, n_samples, pattern, unknown, prefix, data_format, seed):
    from simulate import simulation
    simulation(
        simulate_dir = out,
        data_dir = data,
        sample_size = cells,
        num_samples = n_samples,
        pattern = pattern,
        unknown_celltypes = unknown,
        out_prefix = prefix,
        fmt = data_format,
        seed=seed
    )


def sc_process(pred, training_data, processed_path, cells, scaling, var_cutoff):
    from process import processing
    processing(
        testing_data = pred, 
        training_data = training_data, 
        processed_path = processed_path,
        cells_per_sample = cells,
        scaling = scaling,
        var_cutoff = var_cutoff
    )


def sc_train(processed_path, train_datasets, model_dir, batch_size, learning_rate, steps, seed, loss_values, loss_curve):
    from train import training
    import time
    start = time.time()
    training(
        data_path = processed_path,
        train_datasets = train_datasets,
        model_dir = model_dir,
        batch_size = batch_size,
        learning_rate = learning_rate,
        num_steps = steps,
        seed = seed,
        loss_values = loss_values,
        loss_curve = loss_curve
    )
    end = time.time()
    global training_time
    training_time = round(end - start, 3)


def sc_predict(model_dir, pred, prediction_outname, cells, prediction_scaling, seed):
    from predict import prediction
    prediction(
        model_dir = model_dir, 
        data_path = pred, 
        out_name = prediction_outname, 
        cells = cells, 
        scaling = prediction_scaling, 
        seed = seed
    )

    
def sc_evaluate(out, config, prediction_outname, ground_truth):
    from evaluate import evaluate
    report_filepath = out + f'report_{config}.json'
    evaluate(
        predictions_filename = prediction_outname, 
        ground_truth_filename = ground_truth,
        report_filepath = report_filepath,
        training_time = training_time
    )    
    

@click.command()

# Pipeline controls and logging
@click.option('-load', help='Load parameters from YAML file')
@click.option('-simulate', is_flag=True, help='Run scaden simulate')
@click.option('-process', is_flag=True, help='Run scaden process')
@click.option('-train', is_flag=True, help='Run scaden train')
@click.option('-predict', is_flag=True, help='Run scaden predict')
@click.option('-evaluate', is_flag=True, help='Run evaluation')

@click.option('--no_sim', is_flag=True, help='Skip the creation of simulated training samples if you have already created a training set')
@click.option('--no_proc', is_flag=True, help='Skip the preprocessing step if you have already processed the training data. If using this flag, then --proc is required.')
@click.option('--no_pred', is_flag=True, help='Skip the creation of simulated training samples if you have already created a training set')
@click.option('--no_eval', is_flag=True, help='Skip the evaluation of the model predictions')

@click.option('--config', default='test', help='Name of configuration')
@click.option('--reference', default='1', help='Name of the scMethyl dataset')

# scaden simulate
@click.option('--out', '-o', default='./', help='Directory to store output files in')
@click.option('--data', '-d', default='.', help='Path to scRNA-seq dataset(s)')
@click.option('--cells', '-c', default=100, type=int, help='Number of cells per sample [default: 100]')
@click.option('--n_samples', '-n', default=1000, type=int, help='Number of samples to simulate [default: 1000]')
@click.option('--pattern', default='*_counts.txt', help='File pattern to recognize your processed scRNA-seq count files')
@click.option('--unknown', '-u', multiple=True, default=['unknown'], 
                    help='Specifiy cell types to merge into the unknown category. Specify this flag for every cell type you want '
                    'to merge in unknown. [default: unknown]')
@click.option('--prefix', '-p', default='data', help='Prefix to append to training .h5ad file [default: data]')
@click.option('--data_format', '-f', default='txt', help="Data format of scRNA-seq data, can be 'txt' or 'h5ad' [default: 'txt']")

# scaden process
@click.option('--pred', help='Bulk data file (i.e. testing set) that we want to perform deconvolution on')
@click.option('--training_data', help="Training dataset to be processed. Only use this if you are running 'scaden process' by itself.")
@click.option('--processed_path', help='Name of the file that the processed data will be saved to. Must end with .h5ad')
@click.option('--var_cutoff', default=0.1, type=float,
                    help='Filter out genes with a variance less than the specified cutoff. A low value is recommended, '
                    'this should only remove genes that are obviously uninformative. [default: 0.1]')
@click.option('--scaling', default='fraction', 
                    help='Change scaling option for preprocessing the training data. If something other than the provided options is '
                    'used, then no scaling will be done. [default: fraction] '
                    'Options: None (No scaling), log / log_min_max (log2, then scale to the range 0,1), '
                    'frac / fraction (Divide values by the number of cells)')

# scaden train
@click.option('--train_datasets', default='', help='Comma-separated list of datasets used for training. Uses all by default.')
@click.option('--model_dir', default='./', help='Path to store the model in')
@click.option('--batch_size', default=128, type=int, help='Batch size to use for training. [default: 128]')
@click.option('--learning_rate', default=0.0001, type=float, help='Learning rate used for training. [default: 0.0001]')
@click.option('--steps', default=5000, type=int, help='Number of training steps. [default: 5000]')
@click.option('--seed', default=0, type=int, help='Set random seed')
@click.option('--loss_values', help='Name of file to save text file of loss values')
@click.option('--loss_curve', help='Name of file to save line plot figure of loss values')

# scaden predict
@click.option('--prediction_outname', default='scaden_predictions.txt', help='Name of predictions file')
@click.option('--prediction_scaling', default='fraction', help='Change scaling option for the preprocessing done when making predictions. Uses the same options as --scaling.')

# Evaluate scaden predictions
@click.option('--ground_truth', help='Name of file containing the ground truth cell proportions')

def cli(load, simulate, process, train, predict, evaluate,
        no_sim, no_proc, no_pred, no_eval, config, reference,
        out, data, cells, n_samples, pattern, unknown, prefix, data_format, 
        pred, training_data, processed_path, var_cutoff, scaling, 
        train_datasets, model_dir, batch_size, learning_rate, steps, seed, loss_values, loss_curve,
        prediction_outname, prediction_scaling, 
        ground_truth):
    """
    Runs the scaden commands or the entire pipeline.
    """
    
    """
    Handle args
    """
    param_names = list(locals().keys())
    # Move args to namespace to allow overwriting variables from strings
    from types import SimpleNamespace
    a = SimpleNamespace(**locals())
    
    # Overwrite args with yaml file if provided
    if load:
        import yaml
        with open(load, 'r') as fname:
            params = yaml.safe_load(fname)
        logger.info(f'Using params from: {load}')
        for p in params:
            if p in param_names:           
                setattr(a, p, params[p])
            else:
                click.echo(f"Error: Unknown key '{p}' in YAML file. Keys name must be the same as the long form parameters.", err=True)
                sys.exit()
    
    # Ensure that if any mode is being used, only one mode is provided.
    modes = [a.simulate, a.process, a.train, a.predict, a.evaluate]
    if sum(modes) > 1:
        import numpy as np
        commands = np.array(['simulate', 'process', 'train', 'predict', 'evaluate'])
        click.echo(f"Error: Multiple commands {commands[np.where(modes)]} used. You should only provide one command at a time, or don't use any commands to go through the entire pipeline.", err=True)
        sys.exit()
    
    # Create paths if not specified by user
    if a.training_data is None:
        a.training_data = a.out + a.prefix + '.h5ad' # h5ad file generated from scaden simulate
    
    if a.processed_path is None:
        a.processed_path = a.out + 'processed.h5ad'

    if a.prediction_outname is None:
        a.prediction_outname = a.out + 'scaden_predictions.txt'
    
    # Parameter validation
    param_errors = {}
    for n in ['out', 'data', 'pred', 'model_dir', 'ground_truth']: # input/output directories and files
        if getattr(a, n) is not None:
            if not Path(getattr(a, n)).expanduser().exists():
                param_errors[n] = getattr(a, n)
    
    for n in ['processed_path', 'training_data', 'loss_values', 'loss_curve', 'prediction_outname']: # Generated files
        if getattr(a, n) is not None:
            if not Path(getattr(a, n)).expanduser().parent.exists():
                param_errors[n] = getattr(a, n)

    if len(param_errors) > 0:
        click.echo('Error: parameter directories/files or the parent directories do not exist:', err=True)
        for n in param_errors.keys():
            click.echo(f"'{n}':\t{param_errors[n]}", err=True)
        sys.exit()
    
    """
    SCADEN MODES
    """
    # scaden simulate
    if a.simulate:
        sc_simulate(a.out, a.data, a.cells, a.n_samples, a.pattern, a.unknown, a.prefix, a.data_format, a.seed+1)
        sys.exit()
    
    # scaden process
    if a.process:
        sc_process(a.pred, a.training_data, a.processed_path, a.cells, a.scaling, a.var_cutoff)
        sys.exit()
    
    # scaden train
    if a.train:
        sc_train(a.processed_path, a.train_datasets, a.model_dir, a.batch_size, a.learning_rate, a.steps, a.seed, a.loss_values, a.loss_curve)
        sys.exit()
    
    # scaden predict
    if a.predict:
        sc_predict(a.model_dir, a.pred, a.prediction_outname, a.cells, a.prediction_scaling, a.seed)
        sys.exit()
    
    # scaden evaluate
    if a.evaluate:
        sc_evaluate(a.out, a.config, a.prediction_outname, a.ground_truth)
        sys.exit()
    
    
    """
    PIPELINE
    """
    # Go through entire pipeline if not using one mode
    if a.no_sim != True:
        sc_simulate(a.out, a.data, a.cells, a.n_samples, a.pattern, a.unknown, a.prefix, a.data_format, a.seed)

    if a.no_proc != True:
        sc_process(a.pred, a.training_data, a.processed_path, a.cells, a.scaling, a.var_cutoff)

    sc_train(a.processed_path, a.train_datasets, a.model_dir, a.batch_size, a.learning_rate, a.steps, a.seed, a.loss_values, a.loss_curve)
    global training_time
    logger.info(f'[bold]Training time: [green]{training_time}[/]')

    # Create config.json file
    config_filepath = a.out + f'config_{a.config}.json'
    with open(config_filepath, 'w', encoding='utf-8') as config_file:
        config_json = {
            'config': a.config, 
            'Data': {
                'cells': a.cells, 
                'n_samples': a.n_samples, 
                'var_cutoff': a.var_cutoff,
                'scaling': a.scaling,
                'reference': a.reference
            }, 
            'Model': {
                'seed': a.seed, 
                'steps': a.steps, 
                'batch_size': a.batch_size, 
                'learning_rate': a.learning_rate
            }
        }
        json.dump(config_json, config_file, ensure_ascii=False, indent=4)
    logger.info(f'[bold]Created config file: [green]{config_filepath}[/]')
    
    
    if a.no_pred != True:
        sc_predict(a.model_dir, a.pred, a.prediction_outname, a.cells, a.prediction_scaling, a.seed)

    if a.no_eval != True:
        sc_evaluate(a.out, a.config, a.prediction_outname, a.ground_truth)


def main():
    text = """
     ____                _            
    / ___|  ___ __ _  __| | ___ _ __  
    \___ \ / __/ _` |/ _` |/ _ \ '_ \ 
     ___) | (_| (_| | (_| |  __/ | | |
    |____/ \___\__,_|\__,_|\___|_| |_|
    """
    click.echo(click.style(text, fg='blue'))
    cli()
    

if __name__ == '__main__':
    main()

