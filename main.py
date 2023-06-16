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


@click.command()

# Pipeline controls and logging
@click.option('-load', help='YAML file from which parameters are loaded')
@click.option('-all', is_flag=True, help='Run all steps of the pipeline (simulate, process, train, predict, evaluate)')
@click.option('-simulate', is_flag=True, help='Run scaden simulate')
@click.option('-process', is_flag=True, help='Run scaden process')
@click.option('-train', is_flag=True, help='Run scaden train')
@click.option('-predict', is_flag=True, help='Run scaden predict')
@click.option('-evaluate', is_flag=True, help='Run evaluation')

@click.option('--config', default='test', show_default=True, help='Name of configuration')
@click.option('--reference', default='1', show_default=True, help='Name of the scMethyl dataset')
@click.option('--log_params', is_flag=True, default =True, show_default=True, help='Create a json file recording the data and model hyperparameters')
@click.option('--seed', default=0, type=int, show_default=True, help='Set random seed')

# scaden simulate
@click.option('--out', '-o', default='./', show_default=True, help='Directory to store output files in')
@click.option('--data', '-d', default='.', show_default=True, help='Path to scRNA-seq dataset(s)')
@click.option('--cells', '-c', default=100, type=int, show_default=True, help='Number of cells per sample')
@click.option('--n_samples', '-n', default=1000, type=int, show_default=True, help='Number of samples to simulate')
@click.option('--pattern', default='*_counts.txt', show_default=True, help='File pattern to recognize your processed scRNA-seq count files')
@click.option('--unknown', '-u', multiple=True, default=['unknown'], show_default=True,
                    help='Specifiy cell types to merge into the unknown category. Specify this flag for every cell type you want '
                    'to merge in unknown.')
@click.option('--prefix', '-p', default='data', show_default=True, help='Prefix to append to training .h5ad file')
@click.option('--data_format', '-f', default='txt', show_default=True, help="Data format of scRNA-seq data, can be 'txt' or 'h5ad'")

# scaden process
@click.option('--pred', help='Bulk data file (i.e. testing set) that we want to perform deconvolution on')
@click.option('--training_data', help="Training dataset to be processed. Only use this if you are running 'scaden process' by itself.")
@click.option('--processed_path', help='Name of the file that the processed data will be saved to. Must end with .h5ad')
@click.option('--var_cutoff', default=0.1, show_default=True, type=float,
                    help='Filter out genes with a variance less than the specified cutoff. A low value is recommended, '
                    'this should only remove genes that are obviously uninformative.')
@click.option('--scaling', default='fraction', show_default=True,
                    help='Change scaling option for preprocessing the training data. If something other than the provided options is '
                    'used, then no scaling will be done. [default: fraction] '
                    'Options: None (No scaling), log / log_min_max (log2, then scale to the range 0,1), '
                    'frac / fraction (Divide values by the number of cells)')

# scaden train
@click.option('--train_datasets', default='', help='Comma-separated list of datasets used for training. Uses all by default.')
@click.option('--model_dir', default='./', show_default=True, help='Path to store the model in')
@click.option('--batch_size', default=128, type=int, show_default=True, help='Batch size to use for training.')
@click.option('--learning_rate', default=0.0001, type=float, show_default=True, help='Learning rate used for training.')
@click.option('--steps', default=5000, type=int, show_default=True, help='Number of training steps.')
@click.option('--loss_values', help='Name of file to save text file of loss values')
@click.option('--loss_curve', help='Name of file to save line plot figure of loss values')

# scaden predict
@click.option('--prediction_outname', default='scaden_predictions.txt', show_default=True, help='Name of predictions file')
@click.option('--prediction_scaling', default='fraction', show_default=True, help='Change scaling option for the preprocessing done when making predictions. Uses the same options as --scaling.')

# Evaluate scaden predictions
@click.option('--ground_truth', show_default=True, help='Name of file containing the ground truth cell proportions')

def cli(load, all, simulate, process, train, predict, evaluate, config, reference,
        out, data, cells, n_samples, pattern, unknown, prefix, data_format, log_params,
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
    
    # Create paths if not specified by user
    if a.training_data is None:
        a.training_data = Path(a.out) / (a.prefix + '.h5ad') # h5ad file generated from scaden simulate
    
    if a.processed_path is None:
        a.processed_path = Path(a.out) / 'processed.h5ad'

    if a.prediction_outname is None:
        a.prediction_outname = Path(a.out) / 'scaden_predictions.txt'
    
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
    
    if sum([a.all, a.simulate, a.process, a.train, a.predict, a.evaluate]) == 0:
        click.echo('Error: at least one of the following flags must be provided: all, simulate, process, train, predict, evaluate', err=True)
        sys.exit()
    
    if a.all:
        a.simulate = True
        a.process = True
        a.train = True
        a.predict = True
        a.evaluate = True
    
    """
    SCADEN MODES
    """
    # scaden simulate
    if a.simulate:
        from simulate import simulation
        simulation(
            simulate_dir = a.out,
            data_dir = a.data,
            sample_size = a.cells,
            num_samples = a.n_samples,
            pattern = a.pattern,
            unknown_celltypes = a.unknown,
            out_prefix = a.prefix,
            fmt = a.data_format,
            seed = a.seed
        )
    
    # scaden process
    if a.process:
        from process import processing
        processing(
            testing_data = a.pred, 
            training_data = a.training_data, 
            processed_path = a.processed_path,
            cells_per_sample = a.cells,
            scaling = a.scaling,
            var_cutoff = a.var_cutoff
        )
    
    # scaden train
    training_time = 0
    if a.train:
        from train import training
        import time
        start = time.time()
        training(
            data_path = a.processed_path,
            train_datasets = a.train_datasets,
            model_dir = a.model_dir,
            batch_size = a.batch_size,
            learning_rate = a.learning_rate,
            num_steps = a.steps,
            seed = a.seed,
            loss_values = a.loss_values,
            loss_curve = a.loss_curve
        )
        end = time.time()
        training_time = round(end - start, 3)
        logger.info(f'[bold]Training time: [green]{training_time} seconds[/]')

        # Create config.json file
        config_filepath = Path(a.out) / f'config_{a.config}.json'
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
    
    # scaden predict
    if a.predict:
        from predict import prediction
        prediction(
            model_dir = a.model_dir, 
            data_path = a.pred, 
            out_name = a.prediction_outname, 
            cells = a.cells, 
            scaling = a.prediction_scaling, 
            seed = a.seed
        )
    
    # scaden evaluate
    if a.evaluate:
        from evaluate import evaluate
        report_filepath = Path(a.out) / f'report_{a.config}.json'
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
    click.echo(click.style(text, fg='blue'))
    cli()
    

if __name__ == '__main__':
    main()

