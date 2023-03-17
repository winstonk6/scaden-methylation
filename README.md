# scaden-methylation
Deep learning based cell composition analysis with Scaden using methylation data.

## Overview
Forked from the [original scaden](https://github.com/KevinMenden/scaden) made by Menden et al. and modified to:
- Use methylation data rather than scRNA-seq data.
- Combines all the scaden commands (simulate, process, train, predict) into one program to simplify the process of going through the entire pipeline.
- Adds the option of evaluating the performance of the trained model by calculating the Spearman correlation coefficients with the ground truth proportions.
- Outputs a json file that logs the model hyperparameters and performance.

## Requirements
- Python 3
- The same required packages as the original scaden (which can be downloaded [here](https://github.com/KevinMenden/scaden#installation-guide)).

## Usage
Run main.py to go through the entire scaden pipeline.
1. Simulate training data
2. Preprocess training data
3. Train model
4. Predict on testing data
5. Evaluate predictions

Each step of the scaden pipeline works the same as in the original. Each step has its own set of parameters, so the invocation can be quite long.

**Pipeline controls**
```
  --config TEXT                   Name of configuration.
  --no_sim                        Skip the creation of simulated training
                                  samples if you have already created a
                                  training set.
  --no_proc                       Skip the preprocessing step if you have
                                  already processed the training data. If
                                  using this flag, then --proc is required.
  --no_pred                       Skip the creation of simulated training
                                  samples if you have already created a
                                  training set.
  --no_eval                       Skip the evaluation of the model
                                  predictions.
  -r, --reference TEXT            Number of times each cell type in the basis
                                  is sampled to generate a simulated sc-
                                  methylation dataset.
```

**Simulate**
```
  -o, --out TEXT                  Directory to store output files in
  -d, --data TEXT                 Path to scMethyl-seq dataset(s)
  -c, --cells INTEGER             Number of cells per sample [default: 100]
  -n, --n_samples INTEGER         Number of samples to simulate [default:
                                  1000]
  --pattern TEXT                  File pattern to recognize your processed
                                  scMethyl-seq count files
  -u, --unknown TEXT              Specifiy cell types to merge into the
                                  unknown category. Specify this flag for
                                  every cell type you want to merge in
                                  unknown. [default: unknown]
  -p, --prefix TEXT               Prefix to append to training .h5ad file
                                  [default: data]
  -f, --data_format TEXT          Data format of scMethyl-seq data, can be 'txt'
                                  or 'h5ad' [default: 'txt']
```

**Preprocess**
```
  --pred TEXT                     Bulk data file (i.e. testing set) that we
                                  want to perform deconvolution on.
  -proc, --processed_path TEXT    Name of the file that the processed data
                                  will be saved to. Must end with .h5ad
  --var_cutoff FLOAT              Filter out genes with a variance less than
                                  the specified cutoff. A low value is
                                  recommended, this should only remove genes
                                  that are obviously uninformative. [default:
                                  0.1]
  -sc, --scaling TEXT             Change scaling option for preprocessing the
                                  training data. If something other than the
                                  provided options is used, then no scaling
                                  will be done. [default: fraction] Options:
                                  None (No scaling), log / log_min_max (log2,
                                  then scale to the range 0,1), frac /
                                  fraction (Divide values by the number of
                                  cells)
```

**Train**
```
  --train_datasets TEXT           Comma-separated list of datasets used for
                                  training. Uses all by default.
  -m, --model_dir TEXT            Path to store the model in
  -b, --batch_size INTEGER        Batch size to use for training. [default:
                                  128]
  -l, --learning_rate FLOAT       Learning rate used for training. [default:
                                  0.0001]
  --steps INTEGER                 Number of training steps. [default: 5000]
  --seed INTEGER                  Set random seed
```

**Predict**
```
  -po, --prediction_outname TEXT  Name of predictions file.
  --prediction_scaling TEXT       Change scaling option for the preprocessing
                                  done when making predictions. Uses the same
                                  options as --scaling.
```

**Evaluate**
```
  -g, --ground_truth TEXT         Name of file containing the ground truth
                                  cell proportions.
```

Example invocation (with line breaks):
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
