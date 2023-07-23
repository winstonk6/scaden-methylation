"""
Cell Deconvolutional Network (scaden) class
"""
import os
import logging
import sys
import gc
import tensorflow as tf
import numpy as np
import pandas as pd
from anndata import read_h5ad
import collections
from sklearn.preprocessing import MinMaxScaler
from rich.progress import Progress, BarColumn

logger = logging.getLogger(__name__)


def sample_scaling(x, scaling_option):
    """
    Apply scaling of data
    """

    if scaling_option:
        # Bring in log space
        x = np.log2(x + 1)

        # Normalize data
        mms = MinMaxScaler(feature_range=(0, 1), copy=True)

        # it scales features so transpose is needed
        x = mms.fit_transform(x.T).T

    return x


@tf.function
def compute_loss(logits, targets):
    """
    Compute L1 loss
    """
    loss = tf.reduce_mean(input_tensor=tf.math.square(logits - targets))
    return loss


@tf.function
def compute_accuracy(logits, targets, pct_cut=0.05):
    """
    Compute prediction accuracy
    """
    equality = tf.less_equal(
        tf.math.abs(tf.math.subtract(logits, targets)), pct_cut
    )
    accuracy = tf.reduce_mean(input_tensor=tf.cast(equality, tf.float32))
    return accuracy


@tf.function
def correlation_coefficient(logits, targets):
    """
    Calculate the pearson correlation coefficient
    """
    mx = tf.reduce_mean(input_tensor=logits)
    my = tf.reduce_mean(input_tensor=targets)
    xm, ym = logits - mx, targets - my
    r_num = tf.reduce_sum(input_tensor=tf.multiply(xm, ym))
    r_den = tf.sqrt(
        tf.multiply(
            tf.reduce_sum(input_tensor=tf.square(xm)),
            tf.reduce_sum(input_tensor=tf.square(ym)),
        )
    )
    r = tf.divide(r_num, r_den)
    r = tf.maximum(tf.minimum(r, 1.0), -1.0)
    return r


class Scaden(object):
    """
    scaden class
    """

    def __init__(
            self,
            model_dir,
            model_name,
            batch_size=128,
            learning_rate=0.0001,
            num_steps=1000,
            seed=0,
            hidden_units=(256, 128, 64, 32),
            do_rates=(0, 0, 0, 0),
            cells=1000,
            scaling='fraction'
    ):

        self.model_dir = model_dir
        self.batch_size = batch_size
        self.model_name = model_name
        self.beta1 = 0.9
        self.beta2 = 0.999
        self.learning_rate = learning_rate
        self.data = None
        self.data_iter = None
        self.n_classes = None
        self.labels = None
        self.x = None
        self.y = None
        self.x_data = None
        self.y_data = None
        self.num_steps = num_steps
        self.scaling = scaling
        self.sig_genes = None
        self.sample_names = None
        self.hidden_units = hidden_units
        self.do_rates = do_rates
        self.cells = cells
        self.model = None
        self.global_step = None
        self.loss = None
        self.loss_curve = np.empty(num_steps)

        # Set seeds for reproducibility
        tf.random.set_seed(seed)
        os.environ["TF_DETERMINISTIC_OPS"] = "1"
        np.random.seed(seed)

    def scaden_model(self, n_classes):
        """Create the Scaden model"""

        model = tf.keras.Sequential()
        model.add(tf.keras.layers.Dense(self.hidden_units[0], activation=tf.nn.relu))
        model.add(tf.keras.layers.Dropout(self.do_rates[0]))
        model.add(tf.keras.layers.Dense(self.hidden_units[1], activation=tf.nn.relu))
        model.add(tf.keras.layers.Dropout(self.do_rates[1]))
        model.add(tf.keras.layers.Dense(self.hidden_units[2], activation=tf.nn.relu))
        model.add(tf.keras.layers.Dropout(self.do_rates[2]))
        model.add(tf.keras.layers.Dense(self.hidden_units[3], activation=tf.nn.relu))
        model.add(tf.keras.layers.Dropout(self.do_rates[3]))
        model.add(tf.keras.layers.Dense(n_classes, activation=tf.nn.softmax))

        return model

    def visualization(self, logits, targets, classes):
        """
        Create evaluation metrics
        """
        # add evaluation metrics
        rmse = tf.compat.v1.metrics.root_mean_squared_error(logits, targets)[1]
        pcor = correlation_coefficient(logits, targets)
        eval_metrics = {"rmse": rmse, "pcor": pcor}

        for i in range(logits.shape[1]):
            eval_metrics[
                "mre_" + str(classes[i])
                ] = tf.compat.v1.metrics.mean_relative_error(
                targets[:, i], logits[:, i], targets[:, i]
            )[
                0
            ]
            eval_metrics[
                "mae_" + str(classes[i])
                ] = tf.compat.v1.metrics.mean_absolute_error(
                targets[:, i], logits[:, i], targets[:, i]
            )[
                0
            ]
            eval_metrics["pcor_" + str(classes[i])] = correlation_coefficient(
                targets[:, i], logits[:, i]
            )

        eval_metrics["mre_total"] = tf.compat.v1.metrics.mean_relative_error(
            targets, logits, targets
        )[1]

        eval_metrics["mae_total"] = tf.compat.v1.metrics.mean_relative_error(
            targets, logits, targets
        )[1]

        eval_metrics["accuracy01"] = compute_accuracy(
            logits, targets, pct_cut=0.01
        )
        eval_metrics["accuracy05"] = compute_accuracy(
            logits, targets, pct_cut=0.05
        )
        eval_metrics["accuracy1"] = compute_accuracy(logits, targets, pct_cut=0.1)

        # Create summary scalars
        for key, value in eval_metrics.items():
            tf.compat.v1.summary.scalar(key, value)

        tf.compat.v1.summary.scalar("loss", self.loss)

        merged_summary_op = tf.compat.v1.summary.merge_all()

        return merged_summary_op

    def load_h5ad_file(self, input_path, batch_size, datasets=()):
        """
        Load input data from a h5ad file and divide into training and test set
        :param input_path: path to h5ad file
        :param batch_size: batch size to use for training
        :param datasets: a list of datasets to extract from the file
        :return: Dataset object
        """
        try:
            raw_input = read_h5ad(input_path)
        except OSError:
            logger.error(
                "Could not load training data file! Is it a .h5ad file generated with `scaden process`?"
            )
            sys.exit(1)

        # Subset dataset if --train_datasets is given
        if len(datasets) > 0:
            all_ds = collections.Counter(raw_input.obs["ds"])

            # Check that given datasets are all actually available
            for ds in datasets:
                if ds not in all_ds:
                    logger.warning(
                        f"The dataset '[cyan]{ds}[/cyan]' could not be found in the training data! Is the name correct?"
                    )

            for ds in all_ds:
                if ds not in datasets:
                    raw_input = raw_input[raw_input.obs["ds"] != ds].copy()

        # Create training dataset
        ratios = [raw_input.obs[ctype] for ctype in raw_input.uns["cell_types"]]
        self.x_data = raw_input.X.astype(np.float32)
        self.y_data = np.array(ratios, dtype=np.float32).transpose()
        self.data = tf.data.Dataset.from_tensor_slices((self.x_data, self.y_data))
        self.data = self.data.shuffle(1000).repeat().batch(batch_size=batch_size)
        self.data_iter = iter(self.data)

        # Extract celltype and feature info
        self.labels = raw_input.uns["cell_types"]
        self.sig_genes = list(raw_input.var_names)

    def load_prediction_file(self, input_path, sig_genes, scaling=None):
        """
        Load a file to perform prediction on it
        :param input_path: path to input file
        :param sig_genes: the signature genes to use
        :param scaling: which scaling to perform
        :return: Dataset object
        """
        # Load data
        try:
            data = pd.read_table(input_path, sep="\t", index_col=0)
        except UnicodeDecodeError:
            try:
                data = read_h5ad(input_path).to_df()
            except OSError:
                logger.error('Unsupported file format. The prediction file must be a text file or a h5ad file.')
                sys.exit(1)

        sample_names = list(data.columns)

        # check for duplicates
        data_index = list(data.index)
        if not (len(data_index) == len(set(data_index))):
            logger.warning(
                "Scaden Warning: Your mixture file contains duplicate genes! The first occurring gene will be used for "
                "every duplicate."
            )
            data = data.loc[~data.index.duplicated(keep="first")]

        try:
            data = data.loc[sig_genes]  # sig_genes: from training data set
        except KeyError:
            logger.warning(
                'Genes in the prediction file do not match the genes that the model was trained on. '
                'Filling in the missing genes with zeroes.'
            )

            sig_genes_ind = pd.Index(sig_genes)
            available = sig_genes_ind[sig_genes_ind.isin(data.index)]
            missing = sig_genes_ind[~sig_genes_ind.isin(data.index)]

            zeros_df = pd.DataFrame(np.zeros(shape=(len(missing), len(data.columns))), index=missing,
                                    columns=data.columns)
            data = pd.concat([data.loc[available], zeros_df])

        data = data.T

        # Scaling
        if self.scaling == "log" or self.scaling == "log_min_max":
            data = sample_scaling(data, scaling_option=scaling)

        elif self.scaling == "fraction" or self.scaling == "frac":
            data = data / self.cells

        self.data = data

        return sample_names

    def build_model(self, input_path, train_datasets, mode="train"):
        """
        Build the model graph
        """
        self.global_step = tf.Variable(0, name="global_step", trainable=False)

        # Load training data
        if mode == "train":
            self.load_h5ad_file(
                input_path=input_path,
                batch_size=self.batch_size,
                datasets=train_datasets,
            )

        # Load prediction data
        if mode == "predict":
            self.sample_names = self.load_prediction_file(
                input_path=input_path,
                sig_genes=self.sig_genes,
                scaling=self.scaling,
            )

        # Build the model or load if available
        self.n_classes = len(self.labels)

        try:
            self.model = tf.keras.models.load_model(self.model_dir, compile=False)
            logger.info(f"Loaded pre-trained model: [cyan]{self.model_name}")
        except OSError:
            self.model = self.scaden_model(n_classes=self.n_classes)

    def train(self, input_path, train_datasets):
        """
        Train the model
        """

        # Define the optimizer
        optimizer = tf.keras.optimizers.Adam(learning_rate=self.learning_rate)

        # Build model graph
        self.build_model(
            input_path=input_path, train_datasets=train_datasets, mode="train"
        )

        # Training loop
        progress_bar = Progress(
            "[bold blue]{task.description}",
            "[bold cyan]Step: {task.fields[step]}, Loss: {task.fields[loss]}",
            BarColumn(bar_width=None),
        )

        training_progress = progress_bar.add_task(
            self.model_name, total=self.num_steps, step=0, loss=1
        )

        @tf.function
        def train_step(x_batch_train, y_batch_train):
            # See: https://www.tensorflow.org/guide/keras/writing_a_training_loop_from_scratch
            with tf.GradientTape() as tape:
                self.logits = self.model(x_batch_train, training=True)  # Forward pass
                loss_value = compute_loss(self.logits, y_batch_train)  # Loss values for this batch

            grads = tape.gradient(loss_value, self.model.trainable_weights)  # Get gradients of loss wrt the weights.

            optimizer.apply_gradients(zip(grads, self.model.trainable_weights))  # Update the weights of the model.
            return loss_value

        with progress_bar:
            # Training loop
            for step in range(self.num_steps):
                x, y = self.data_iter.get_next()

                # Do training step and record loss
                loss = train_step(x, y)
                self.loss_curve[step] = loss

                progress_bar.update(
                    training_progress, advance=1, step=step, loss=f"{loss:.4f}"
                )

                # Collect garbage after 100 steps - otherwise runs out of memory
                if step % 100 == 0:
                    gc.collect()

        # Save the trained model
        self.model.save(self.model_dir)
        pd.DataFrame(self.labels).to_csv(
            os.path.join(self.model_dir, "celltypes.txt"), sep="\t"
        )
        pd.DataFrame(self.sig_genes).to_csv(
            os.path.join(self.model_dir, "genes.txt"), sep="\t"
        )

    def predict(self, input_path):
        """
        Perform prediction with a pre-trained model
        :param input_path: prediction data path
        """
        # Load signature genes and celltype labels
        sig_genes = pd.read_table(self.model_dir + "/genes.txt", index_col=0)
        self.sig_genes = list(sig_genes["0"])
        labels = pd.read_table(self.model_dir + "/celltypes.txt", index_col=0)
        self.labels = list(labels["0"])

        # Build model graph
        self.build_model(input_path=input_path, train_datasets=(), mode="predict")

        predictions = self.model.predict(self.data)

        pred_df = pd.DataFrame(
            predictions, columns=self.labels, index=self.sample_names
        )
        return pred_df
