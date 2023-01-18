#!/bin/bash
module load python
source ~/venv/scaden/bin/activate

show_help() {
cat << EOF

Run scaden simulate, scaden process, and scaden train. 
Note:
    The sc-RNA dataset files must be formatted as *_counts.txt and *_celltypes.txt
    The --prediction_file argument is required. This file must be generated beforehand and separately.
    Directory names must end with /
    The model can only be trained using data stored in a single h5ad file.

Usage
    --config                   Name of configuration. Will create a text file containing the
                               parameters.
    --no_sim                   Skip the 'scaden simulate' step that creates a training set.
    --no_proc                  Skip the 'scaden process' step that preprocesses the training set.
    
    scaden simulate
      -o, --out TEXT           Directory to store output files in
      -d, --data TEXT          Path to scRNA-seq dataset(s)
      -c, --cells INTEGER      Number of cells per sample [default: 100]
      -n, --n_samples INTEGER  Number of samples to simulate [default: 1000]
      --pattern TEXT           (NOTE: --pattern is not used. Format scRNA-seq files as
                                the default *_counts.txt and *_celltypes.txt)
                               File pattern to recognize your processed scRNA-seq
                               count files
      -u, --unknown TEXT       Specifiy cell types to merge into the unknown
                               category. Specify this flag for every cell type you
                               want to merge in unknown. [default: unknown]
      -p, --prefix TEXT        Prefix to append to training .h5ad file [default:
                               data]
      -f, --data-format TEXT   Data format of scRNA-seq data, can be 'txt' or
                               'h5ad' [default: 'txt']

    scaden process
      --pred, --prediction_file  (REQUIRED) Bulk data file you want to perform deconvolution on.
      --proc, --processed_path   Path of processed file. Must end with .h5ad
      --var_cutoff FLOAT         Filter out genes with a variance less than the
                                 specified cutoff. A low cutoff is recommended,this
                                 should only remove genes that are obviously
                                 uninformative.

    scaden train
      --train_datasets TEXT      (NOTE: --train_datasets is unused.
                                  Script will use preprocessed data file automatically.)
                                 Comma-separated list of datasets used for training.
                                 Uses all by default.
      -m, --model_dir TEXT       Path to store the model in
      -b, --batch_size INTEGER   Batch size to use for training. Default: 128
      -l, --learning_rate FLOAT  Learning rate used for training. Default: 0.0001
      --steps INTEGER            Number of training steps
      --seed INTEGER

EOF
}

err() {
    printf '%s\n' "$1" >&2
    exit 1
}

config="test"
no_sim=false
no_proc=false

out="./"
data="."
cells=100
n_samples=1000
pattern='"*_counts.txt"'
unknown="unknown"
prefix="data"
dataformat="txt"

prediction_file=
processed_path=
var_cutoff="0.1"

train_datasets=""
model_dir="./"
batch_size=128
learning_rate="0.0001"
steps=5000
seed=$RANDOM

# If flag is provided, override default.
while :; do
    case $1 in
        -h|-\?|--help)
            show_help    # Display a usage synopsis.
            exit
            ;;
        -o|--out)
            if [ "$2" ]; then
                out=$2
                shift
            else
                err 'ERROR: "--file" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        -d|--data)
            if [ "$2" ]; then
                data=$2
                shift
            else
                err 'ERROR: "--data" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        -c|--cells)
            if [ "$2" ]; then
                cells=$2
                shift
            else
                err 'ERROR: "--cells" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        -n|--n_samples)
            if [ "$2" ]; then
                n_samples=$2
                shift
            else
                err 'ERROR: "--n_samples" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        --pattern)
            if [ "$2" ]; then
                pattern=$2
                shift
            else
                err 'ERROR: "--pattern" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        -u|--unknown)
            if [ "$2" ]; then
                unknown=$2
                shift
            else
                err 'ERROR: "--unknown" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        -p|--prefix)
            if [ "$2" ]; then
                prefix=$2
                shift
            else
                err 'ERROR: "--prefix" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        -f|--data-format)
            if [ "$2" ]; then
                dataformat=$2
                shift
            else
                err 'ERROR: "--data-format" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        --proc|--processed_path)
            if [ "$2" ]; then
                processed_path=$2
                shift
            else
                err 'ERROR: "--file" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        --var_cutoff)
            if [ "$2" ]; then
                var_cutoff=$2
                shift
            else
                err 'ERROR: "--var_cutoff" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        --train_datasets)
            if [ "$2" ]; then
                train_datasets=$2
                shift
            else
                err 'ERROR: "--train_datasets" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        -m|--model_dir)
            if [ "$2" ]; then
                model_dir=$2
                shift
            else
                err 'ERROR: "--model_dir" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        -b|--batch_size)
            if [ "$2" ]; then
                batch_size=$2
                shift
            else
                err 'ERROR: "--batch_size" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        -l|--learning_rate)
            if [ "$2" ]; then
                learning_rate=$2
                shift
            else
                err 'ERROR: "--learning_rate" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        --steps)
            if [ "$2" ]; then
                steps=$2
                shift
            else
                err 'ERROR: "--steps" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        --seed)
            if [ "$2" ]; then
                seed=$2
                shift
            else
                err 'ERROR: "--seed" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        --pred|--prediction_file)
            if [ "$2" ]; then
                prediction_file=$2
                shift
            else
                err 'ERROR: "--prediction_data" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        --config)
            if [ "$2" ]; then
                config=$2
                shift
            else
                err 'ERROR: "--config" requires a non-empty option argument.'
                exit 1
            fi
            ;;
        --no_sim)
            no_sim=true
            ;;
        --no_proc)
            no_proc=true
            ;;
        --)              # End of all options.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)               # Default case: No more options, so break out of the loop.
            break
    esac

    shift
done

# Training data h5ad file path
trainingdata=${out}${prefix}".h5ad"

# Preprocessing step: if no path is specified, then use the default name 'processed.h5ad' and put in the same directory as the unprocessed training data.
if ! [ $processed_path ]; then
    processed_path=${out}"processed.h5ad"
fi

# Simulate bulk data from sc-RNA counts dataset
if [ $no_sim = false ]; then
    scaden simulate -o $out -d $data -c $cells -n $n_samples -u $unknown -p $prefix -f $dataformat
    mv ${prefix}".h5ad" $out
fi

# Preprocess bulk data with modified preprocessing function
if [ $no_proc = false ]; then
    python /mnt/ceph/users/wkoh1/scaden/pipeline/process_scaden.py --prediction_file $prediction_file --training_data $trainingdata --processed_path $processed_path --cells $cells --var_cutoff $var_cutoff
fi

# Train model with provided parameters
scaden train $processed_path --model_dir $model_dir --batch_size $batch_size --learning_rate $learning_rate --steps $steps --seed $seed

# Record data + pipeline parameters
cat > config_${config}.json << EOF
{
  "config": "${config}",
  "Data": {
    "cells": ${cells},
    "n_samples": ${n_samples}
  },
  "Model": {
    "seed": ${seed},
    "steps": ${steps},
    "batch_size": ${batch_size},
    "learning_rate": ${learning_rate}
  }
}
EOF
