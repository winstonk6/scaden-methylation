# scaden-methylation
Deep learning based cell composition analysis with Scaden using methylation data.

## Overview
[Scaden](https://github.com/KevinMenden/scaden) (made by Menden et al.) is a deep learning based tool that predicts the cell type proportions of tissues from RNA-seq by  training on a separate scRNA-seq dataset.

Scaden-methylation modifies the original scaden in a few ways:
- Uses methylation data rather than RNA-seq data.
- Can combine all the scaden commands (simulate, process, train, predict) into one program to simplify the process of going through the entire pipeline.
- Saves the training loss values.
- Adds the option of evaluating the performance of the trained model by calculating the Pearson correlation coefficients with the ground truth proportions.
- Outputs a json file that logs the model hyperparameters and performance.

The pipeline has 5 steps:
1. Simulate training data
2. Preprocess training data
3. Train model
4. Predict on testing data
5. Evaluate predictions

## Requirements
- Python 3
- Install required python packages using `python -m pip install -r requirements.txt` or `conda install --file requirements.txt`

## Usage
An example can be found in the Jupyter Notebook `scaden_pipeline_example.ipynb`.

To start the pipeline, you need a scMethyl-seq dataset and a methylation array that you want to perform deconvolution on. If you have a scMethyl-seq dataset, it must be formatted as explained [here](https://scaden.readthedocs.io/en/latest/usage.html#scrna-seq-data-processing).
You can also use a methylation array to simulate a scMethyl-seq dataset with the proper formatting using the script `array_to_sc.py`.
```
usage: array_to_sc [-h] [-r REPEATS] [-o OUTNAME] [-s SEED] filename

Converts methylation array to a simulated scMethyl-seq dataset.

positional arguments:
  filename                        Methylation array text file

optional arguments:
  -r, --repeats INTEGER           Number of times each beta value is 
                                  sampled from a binomial distribution.
  -o, --outname TEXT              Set the name and location of the outputs, 
                                  e.g., "./data" will create the files 
                                  "./data_counts.txt" and "./data_celltypes.txt"
  -s, --seed INTEGER              Numpy random generator seed
```

Once you have the necessary files, run `main.py` to go through the entire scaden pipeline and provide the necessary arguments as a YAML file or as command line arguments. Using a YAML file is recommended, since each step has its own set of parameters and the invocation can be quite long.

You can also run each command individually by adding the appropriate flag on the command line (e.g., `-simulate`) or in the YAML file (`simulate: True`)

**Pipeline controls and logging**
```
  -load TEXT                 Load parameters from YAML file
  -simulate                  Run scaden simulate
  -process                   Run scaden process
  -train                     Run scaden train
  -predict                   Run scaden predict
  -evaluate                  Run evaluation
  --no_sim                   Skip the creation of simulated training samples
                             if you have already created a training set
  --no_proc                  Skip the preprocessing step if you have already
                             processed the training data. If using this flag,
                             then --proc is required.
  --no_pred                  Skip the creation of simulated training samples
                             if you have already created a training set
  --no_eval                  Skip the evaluation of the model predictions
  --config TEXT              Name of configuration
  --reference TEXT           Name of the scMethyl dataset
  --seed INTEGER             Set random seed for simulation and training.
```

**Simulate**
```
  -o, --out TEXT             Directory to store output files in
  -d, --data TEXT            Path to scRNA-seq dataset(s)
  -c, --cells INTEGER        Number of cells per sample [default: 100]
  -n, --n_samples INTEGER    Number of samples to simulate [default: 1000]
  --pattern TEXT             File pattern to recognize your processed scRNA-
                             seq count files
  -u, --unknown TEXT         Specifiy cell types to merge into the unknown
                             category. Specify this flag for every cell type
                             you want to merge in unknown. [default: unknown]
  -p, --prefix TEXT          Prefix to append to training .h5ad file [default:
                             data]
  -f, --data_format TEXT     Data format of scRNA-seq data, can be 'txt' or
                             'h5ad' [default: 'txt']
```

**Preprocess**
```
  --pred TEXT                Bulk data file (i.e. testing set) that we want to
                             perform deconvolution on
  --training_data TEXT       Training dataset to be processed. Only use this
                             if you are running 'scaden process' by itself.
  --processed_path TEXT      Name of the file that the processed data will be
                             saved to. Must end with .h5ad
  --var_cutoff FLOAT         Filter out genes with a variance less than the
                             specified cutoff. A low value is recommended,
                             this should only remove genes that are obviously
                             uninformative. [default: 0.1]
  --scaling TEXT             Change scaling option for preprocessing the
                             training data. If something other than the
                             provided options is used, then no scaling will be
                             done. [default: fraction] Options: None (No
                             scaling), log / log_min_max (log2, then scale to
                             the range 0,1), frac / fraction (Divide values by
                             the number of cells)
```

**Train**
```
  --train_datasets TEXT      Comma-separated list of datasets used for
                             training. Uses all by default.
  --model_dir TEXT           Path to store the model in
  --batch_size INTEGER       Batch size to use for training. [default: 128]
  --learning_rate FLOAT      Learning rate used for training. [default:
                             0.0001]
  --steps INTEGER            Number of training steps. [default: 5000]
  --loss_values TEXT         Name of file to save text file of loss values
  --loss_curve TEXT          Name of file to save line plot figure of loss
                             values
```

**Predict**
```
  --prediction_outname TEXT  Name of predictions file
  --prediction_scaling TEXT  Change scaling option for the preprocessing done
                             when making predictions. Uses the same options as
                             --scaling.
```

**Evaluate**
```
  --ground_truth TEXT        Name of file containing the ground truth cell
                             proportions
```

## Examples
Using a YAML file:
`python main.py -load params.yaml`

`params.yaml` file:
```
data: ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/
out: ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/
model_dir: ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/model/
pred: ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/test_bulk_data.txt
ground_truth: ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/test_actual.txt
prediction_outname: ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/test_predictions.txt
scaling: frac
prediction_scaling: None
var_cutoff: 0
seed: 42
cells: 200
n_samples: 2000
steps: 6000
batch_size: 256
learning_rate: 0.0001
config: 42_c200_n2000_s6000_b256_l0.0001_r200_g0.1_c0.2
reference: 200
```


Example invocation using command line (with line breaks):
```
python main.py \
--data ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/ \
--out ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/ \
--model_dir ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/model/ \
--pred ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/test_bulk_data.txt \
--ground_truth ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/test_actual.txt \
--prediction_outname ~/ceph/scaden/simulated_basis_test/g0.1_c0.2/test_predictions.txt \
--scaling frac \
--prediction_scaling None \
--var_cutoff 0 \
--seed 42 \
--cells 200 \
--n_samples 2000 \
--steps 6000 \
--batch_size 256 \
--learning_rate 0.0001 \
--config 42_c200_n2000_s6000_b256_l0.0001_r200_g0.1_c0.2 \
--reference 200
```

Running one command at a time:
```
python main.py -predict -load params.yaml
```
