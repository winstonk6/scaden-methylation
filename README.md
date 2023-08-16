# scaden-methylation
Deep learning based cell composition analysis with Scaden using methylation data.

## Overview
[Scaden](https://github.com/KevinMenden/scaden) (made by Menden et al.) is a deep learning based tool that predicts the cell type proportions of tissues from RNA-seq by  training on a separate scRNA-seq dataset.

Scaden-methylation modifies the original scaden in a few ways:
- Can work with methylation data in addition to RNA-seq data.
- Combines all the scaden commands (simulate, process, train, predict) into one command to simplify the pipeline.
- Adds the option of evaluating the performance of the trained model by calculating correlations and MSE with the ground truth proportions.
- Adds and simplifies logging of parameters, performance, etc.

The pipeline has 5 steps:
1. Simulate training data
2. Preprocess training data
3. Train model
4. Predict on testing data
5. Evaluate predictions

## Requirements
- Python 3
- Install required python packages using `python -m pip install -r requirements.txt` or `conda install --file requirements.txt`
- Follow the instructions [here](https://www.tensorflow.org/install/pip) to install TensorFlow.

## Usage
See the Jupyter Notebook [scaden_pipeline_example.ipynb](scaden_pipeline_example.ipynb) for a more detailed explanation and example.

### TL;DR:
1. Find or make a single-cell dataset to use as a reference for training a model.

2. Find a dataset of bulk samples to perform deconvolution on. This should be a text or h5ad file with rows as genes and samples as columns.

3. Make a YAML file with the necessary parameters for your task.

4. Run `python main.py -load params.yaml` 

You can also add arguments through the command line, e.g., `python main.py -simulate --cells 100 ...`

### Single-cell dataset
The single cell dataset should be split into two files in the same directory:
  - One file should end in "`_celltypes.txt`". This file will have one column called "Celltype" with the celltype labels for each cell.
    ```
    Celltype
    type1
    type2
    type9
    type4
    ```
  - The other file should end in "`_counts.txt`" and have the table of methylation values.
    ```
        cpg1    cpg2    cpg100  cpg230
    0    0.0     1.0     0.5     0.0
    1    1.0     0.25    0.0     0.0
    2    0.0     0.0     1.0     0.0
    3    0.75    0.125   0.2     1.0
    ```
Both files should have the same prefix (e.g., `neuron_celltypes.txt` and `neuron_counts.txt`)

### Example YAML files:
  - Use a trained model to make predictions:
    ```
    predict: True
    seed: 100

    # Directory containing the trained model
    model_dir: /mnt/home/wkoh1/ceph/scaden/basis2/configs/new/model/

    # Dataset to do deconvolution, with cpgs as rows and samples as columns
    pred: /mnt/ceph/users/wkoh1/scaden/Scaden_testset/GSE182379/m.csv

    # Scaling to apply to the dataset values. 
    # Options: None, log, frac (log: log2 then scale to the range [0,1], frac: divide by some number N, requires the line below to be uncommented) 
    prediction_scaling: None
    #cells: N

    # Name of the output file
    prediction_outname: /mnt/home/wkoh1/ceph/scaden_predictions.tsv
    ```
  
  - Train a model:
    ```
    simulate: True
    process: True
    train: True
    seed: 100

    # Simulation parameters
    # Directory where the single cell dataset files (*_counts.txt and *_celltypes.txt) are located
    data: /mnt/home/wkoh1/ceph/GSE130711/scaden/original/

    # Number of samples to simulate in the dataset
    n_samples: 10

    # Number of cells per simulated sample
    cells: 10

    # Directory to store the output files
    out: /mnt/home/wkoh1/ceph/scaden/scaden_example/

    # Optional: name of the file containing the training data. This file be stored in the out directory.
    prefix: training

    # Bulk data file (or testing set) that we want to perform deconvolution on. 
    pred: /mnt/home/wkoh1/ceph/GSE130711/scaden/original/test_data.tsv

    # Scaling to apply to the training dataset values. 
    # Options: None, log, frac, frac_notna
    scaling: frac

    # Remove cpgs from the training data that have a variance below this cutoff. Use 0 to keep all cpgs.
    var_cutoff: 0

    # Directory to store the model files
    model_dir: /mnt/home/wkoh1/ceph/scaden/scaden_example/model/

    # Model hyperparameters
    batch_size: 256
    learning_rate: 0.0001
    steps: 10
    ```


# All Options
**Pipeline controls and logging**
```
  -load TEXT                 YAML file from which parameters are loaded
  -v, --verify               Verify that parameters are valid, then exit
                             without doing any steps.
  -all                       Run all steps of the pipeline (simulate, process,
                             train, predict, evaluate)
  -simulate                  Run scaden simulate
  -process                   Run scaden process
  -train                     Run scaden train
  -predict                   Run scaden predict
  -evaluate                  Run evaluation
  --config TEXT              Name of configuration  [default: test]
  --reference TEXT           Name of the dataset  [default: None]
  --log_params               Create a json file recording the data and model
                             hyperparameters
  --seed INTEGER             Set random seed  [default: 0]
```

**Simulate**
```
  -o, --out TEXT             Directory to store output files in  [default: ./]
  -d, --data TEXT            Path to scRNA-seq dataset(s)  [default: .]
  -c, --cells INTEGER        Number of cells per sample  [default: 100]
  -n, --n_samples INTEGER    Number of samples to simulate  [default: 1000]
  --pattern TEXT             File pattern to recognize your processed scRNA-
                             seq count files  [default: *_counts.txt]
  -u, --unknown TEXT         Specifiy cell types to merge into the unknown
                             category. Specify this flag for every cell type
                             you want to merge in unknown.  [default: unknown]
  -p, --prefix TEXT          Prefix to append to training .h5ad file
                             [default: data]
  -f, --data_format TEXT     Data format of scRNA-seq data, can be 'txt' or
                             'h5ad'  [default: txt]
```

**Preprocess**
```
  --pred TEXT                Bulk data file (or testing set) that we want to
                             perform deconvolution on. Should be a text or
                             h5ad file with rows as genes and samples as
                             columns.
  --training_data TEXT       Training dataset to be processed. Only use this
                             if you are running 'scaden process' by itself.
  --processed_path TEXT      Name of the file that the processed data will be
                             saved to. Must end with .h5ad
  --var_cutoff FLOAT         Filter out genes with a variance less than the
                             specified cutoff. A low value is recommended,
                             this should only remove genes that are obviously
                             uninformative.  [default: 0.1]
  --scaling TEXT             Change scaling option for preprocessing the
                             training data. If something other than the
                             provided options is used, then no scaling will be
                             done. Options: None (No scaling), log /
                             log_min_max (log2, then scale to the range 0,1),
                             frac / fraction (Divide values by the number of
                             cells), frac_notna (Divide values by the number
                             of non-NA values)  [default: fraction]
```

**Train**
```
  --train_datasets TEXT      Comma-separated list of datasets used for
                             training. Uses all by default.
  --model_dir TEXT           Path to store the model in  [default: ./]
  --batch_size INTEGER       Batch size to use for training.  [default: 128]
  --learning_rate FLOAT      Learning rate used for training.  [default:
                             0.0001]
  --steps INTEGER            Number of training steps.  [default: 5000]
  --loss_values TEXT         Name of file to save text file of loss values
  --loss_curve TEXT          Name of file to save line plot figure of loss
                             values
```

**Predict**
```
  --prediction_outname TEXT  Name of predictions file  [default:
                             scaden_predictions.txt]
  --prediction_scaling TEXT  Change scaling option for the preprocessing done
                             when making predictions. Uses the same options as
                             --scaling.  [default: fraction]
```

**Evaluate**
```
  --ground_truth TEXT        Name of file containing the ground truth cell
                             proportions
```

