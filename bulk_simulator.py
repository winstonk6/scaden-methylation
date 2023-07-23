import logging
import glob
import os
import sys
import gc

import pandas as pd
import anndata as ad
import numpy as np

from rich.progress import BarColumn, Progress

logger = logging.getLogger(__name__)


class BulkSimulator:
    """
    BulkSimulator class for the simulation of artificial bulk samples
    from scRNA-seq datasets

    Args
        sample_size: int, default = 100
            Number of cells per sample.

        num_samples: int, default = 1000
            Number of samples to simulate.

        data_path: str, default = './'
            Path to the data directory that contains the single cell datasets (as counts and celltypes files).

        out_dir: str, default = './'
            Output directory where the simulated samples will be saved.

        pattern: str, default = '*_counts.txt'
            File pattern to recognize the counts files

        unknown_celltypes: list, default = None
            List of celltype names to merge into the unknown class

        fmt: str, default = 'txt'
            Format of the single cell dataset files, can be txt or h5ad

        seed: {None, int}, default = None
            Seed to initialize the random number generator.

    """

    def __init__(
            self,
            sample_size: int = 100,
            num_samples: int = 1000,
            data_path: str = './',
            out_dir: str = './',
            pattern: str = '*_counts.txt',
            unknown_celltypes: list | None = None,
            fmt: str = 'txt',
            seed: str | None = None
    ):

        if unknown_celltypes is None:
            unknown_celltypes = ['unknown']

        self.sample_size = sample_size
        self.num_samples = num_samples // 2
        self.data_path = data_path
        self.out_dir = out_dir
        self.pattern = pattern
        self.unknown_celltypes = unknown_celltypes
        self.format = fmt
        self.datasets = []
        self.dataset_files = []
        self.rng = np.random.default_rng(seed)

    def simulate(self) -> None:
        """ Simulate artificial bulk datasets """

        # List available datasets
        if not self.data_path.endswith('/'):
            self.data_path += '/'
        files = glob.glob(os.path.join(self.data_path, self.pattern))
        files = [os.path.basename(x) for x in files]
        self.datasets = [x.replace(self.pattern.replace('*', ''), '') for x in files]
        self.dataset_files = [
            os.path.join(self.out_dir, x + '.h5ad') for x in self.datasets
        ]
        if len(self.datasets) == 0:
            logging.error(
                'No datasets found! Have you specified the pattern correctly?'
            )
            sys.exit(1)

        logger.info('Datasets: [cyan]' + str(self.datasets) + '[/]')

        # Loop over datasets and simulate bulk data
        for dataset in self.datasets:
            gc.collect()
            logger.info(f'[bold u]Simulating data from {dataset}')
            self.simulate_dataset(dataset)

        logger.info('[bold green]Finished data simulation!')

    def merge_unknown_celltypes(self, y: pd.DataFrame) -> pd.DataFrame:
        """
        Modify the celltypes by merging all specified unknown celltypes together into the category 'Unknown'.

        Args
            y: pd.DataFrame
                Contains all the cell types in the celltypes file

        Returns
            pd.DataFrame
        """
        celltypes = list(y['Celltype'])
        new_celltypes = [
            'Unknown' if x in self.unknown_celltypes else x for x in celltypes
        ]
        y['Celltype'] = new_celltypes
        return y

    def load_dataset(self, dataset: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Load the counts and celltypes of a single cell dataset

        Args
            dataset: str
                Name of the dataset

        Returns
            x: pandas.DataFrame
                Contains the counts

            y: pandas.DataFrame
                Contains the cell types
        """
        pattern = self.pattern.replace('*', '')
        logger.info(f'Loading [cyan]{dataset}[/] dataset ...')
        dataset_counts = dataset + pattern
        dataset_celltypes = dataset + '_celltypes.txt'

        # Load data in .txt format
        if self.format == 'txt':
            # Try to load celltypes
            try:
                y = pd.read_table(os.path.join(self.data_path, dataset_celltypes))
                # Check if file has Celltype column
                if 'Celltype' not in y.columns:
                    logger.error(
                        f"No 'Celltype' column found in {dataset}_celltypes.txt! Please make sure to include this "
                        f"column. "
                    )
                    sys.exit(1)
            except FileNotFoundError:
                logger.error(
                    f'No cell types file found for [cyan]{dataset}[/]. It should be called [cyan]{dataset}_celltypes'
                    f'.txt.'
                )
                sys.exit(1)

            # Try to load data file
            try:
                x = pd.read_table(
                    os.path.join(self.data_path, dataset_counts),
                    index_col=0,
                    dtype=np.float32,
                )
            except FileNotFoundError:
                logger.error(
                    f'No counts file found for [cyan]{dataset}[/]. Was looking for file [cyan]{dataset_counts}[/]'
                )
                sys.exit(1)

            # Check that celltypes and count file have the same number of cells
            if y.shape[0] != x.shape[0]:
                logger.error(
                    f'Different number of cells in {dataset}_celltypes and {dataset_counts}! Make sure the data has '
                    f'been processed correctly. '
                )
                sys.exit(1)

        # Load data in .h5ad format
        elif self.format == 'h5ad':
            try:
                data_h5ad = ad.read_h5ad(os.path.join(self.data_path, dataset_counts))
            except FileNotFoundError:
                logger.error(
                    f'No h5ad file found for [cyan]{dataset}[/]. Was looking for file [cyan]{dataset_counts}'
                )
                sys.exit(1)
            # cell types
            try:
                y = pd.DataFrame(data_h5ad.obs.Celltype)
                y.reset_index(inplace=True, drop=True)
            except AttributeError:
                logger.error(f'Celltype attribute not found for [cyan]{dataset}')
                sys.exit(1)
            # counts
            x = data_h5ad.to_df()
            del data_h5ad
        else:
            logger.error(f'Unsupported file format {self.format}!')
            sys.exit(1)

        return x, y

    def simulate_dataset(self, dataset: str) -> None:
        """
        Simulate bulk data from a single scRNA-seq dataset and creates a h5ad file

        Args
            dataset: str
                Name of dataset
        """

        # load the dataset
        data_x, data_y = self.load_dataset(dataset)

        # Merge unknown celltypes
        logger.info(f'Merging unknown cell types: {self.unknown_celltypes}')
        data_y = self.merge_unknown_celltypes(data_y)

        logger.info(f'Subsampling [bold cyan]{dataset}[/] ...')

        # Extract celltypes
        celltypes = pd.Series(data_y['Celltype'].unique())

        # Create simulated samples
        sim_x, sim_y = self.create_subsample_dataset(
            data_x, data_y, celltypes=celltypes
        )
        sim_x.sort_index(axis=1, inplace=True)
        sim_y['ds'] = pd.Series(np.repeat(dataset, sim_y.shape[0]), index=sim_y.index)

        ann_data = ad.AnnData(
            X=sim_x.to_numpy(),
            obs=sim_y,
            var=pd.DataFrame(columns=[], index=list(sim_x)),
        )
        ann_data.uns['unknown'] = self.unknown_celltypes
        ann_data.uns['cell_types'] = celltypes.tolist()

        ann_data.write(os.path.join(self.out_dir, dataset + '.h5ad'))

    def create_subsample_dataset(self, data_x: pd.DataFrame, data_y: pd.DataFrame, celltypes: pd.Series) -> \
            tuple[pd.DataFrame, pd.DataFrame]:
        """
        For one dataset, generate many artificial bulk samples with known fractions.
        This function will create normal and sparse (missing some cell types) samples

        Args
            data_x: pandas.DataFrame
                Counts data, with samples as rows and genes as columns

            data_y: pandas.DataFrame
                Corresponding celltypes for the counts data as one column

            celltypes: list
                All celltypes in the celltypes file

        Returns
            sim_x: pandas.DataFrame
                Simulated counts

            sim_y: pandas.DataFrame
                Simulated fractions, with samples as rows and celltypes as columns
        """
        sim_y = []
        sim_x = np.empty(shape=(self.num_samples * 2, data_x.shape[1]), dtype=np.float32)

        progress_bar = Progress(
            "[bold blue]{task.description}",
            "[bold cyan]{task.fields[samples]}",
            BarColumn(bar_width=None),
        )
        with progress_bar:
            normal_samples_progress = progress_bar.add_task(
                "Normal samples", total=self.num_samples, samples=0
            )
            sparse_samples_progress = progress_bar.add_task(
                "Sparse samples", total=self.num_samples, samples=0
            )

            # Create normal samples
            for i in range(self.num_samples):
                progress_bar.update(normal_samples_progress, advance=1, samples=i + 1)
                sample, label = self.create_subsample(data_x, data_y, celltypes)
                sim_y.append(label)
                sim_x[i, :] = sample

            # Create sparse samples
            for i in range(self.num_samples):
                progress_bar.update(sparse_samples_progress, advance=1, samples=i + 1)
                sample, label = self.create_subsample(data_x, data_y, celltypes, sparse=True)
                # sim_x.append(sample)
                sim_y.append(label)
                sim_x[self.num_samples + i, :] = sample

        sim_x = pd.DataFrame(sim_x, columns=data_x.columns)
        sim_y = pd.concat(sim_y, axis=1).T
        return sim_x, sim_y

    def create_subsample(self, data_x: pd.DataFrame, data_y: pd.DataFrame, celltypes: pd.Series, sparse: bool = False) \
            -> tuple[np.ndarray, pd.Series]:
        """
        Generate one bulk sample with random fractions of celltypes.
        If sparse is set to true, the sample will be composed of a random subset of all celltypes

        Args
            data_x: pandas.DataFrame
                Counts data, with samples as rows and genes as columns

            data_y: pandas.DataFrame
                Corresponding celltypes for the counts data as one column

            celltypes: pd.Series
                All celltypes in the celltypes file

            sparse: bool, default = False
                Create sparse samples

        Returns
            df_samp:
                Counts for the simulated samples

            fracs_complete:
                Cell type fractions for the simulated samples
        """
        if not sparse:
            available_celltypes = np.array(celltypes)
        else:
            # Randomly choose which cell types to sample in order to create a sparse sample
            num_keep = self.rng.integers(low=1, high=len(celltypes))
            keep_index = self.rng.choice(a=len(celltypes), size=num_keep, replace=False)
            available_celltypes = np.array(celltypes)[keep_index]

        # Create fractions for available celltypes
        fractions = self.create_fractions(num_celltypes=available_celltypes.shape[0])

        # Calculate the number of cells to use for each cell type based on the total number of cells
        celltype_count = np.round(np.multiply(fractions, self.sample_size)).astype(int)
        total_cells = np.sum(celltype_count)

        # Record fractions of all cell types
        not_kept = celltypes[~celltypes.isin(available_celltypes)]
        not_kept_fractions = np.zeros(shape=(len(not_kept)))
        fracs_complete = pd.Series(np.concatenate((fractions, not_kept_fractions)),
                                   index=pd.concat([pd.Series(available_celltypes), not_kept]))

        simulated_sample = np.empty(shape=(total_cells, len(data_x.columns)), dtype=np.float32)
        start_row = 0  # counter for row to fill simulated_sample array
        for n, celltype in enumerate(available_celltypes):
            cells_subset = data_x[data_y['Celltype'] == celltype]

            # Randomly pick rows from the subset with equal probability
            choose_cells = self.rng.integers(low=0, high=cells_subset.shape[0], size=celltype_count[n])
            simulated_sample[start_row: start_row + celltype_count[n], :] = cells_subset.iloc[choose_cells].to_numpy()
            start_row += celltype_count[n]  # add to counter

        # Aggregate counts
        df_samp = np.nansum(simulated_sample, axis=0)

        return df_samp, fracs_complete

    def create_fractions(self, num_celltypes: int) -> np.ndarray:
        """
        Create random fractions that add up to 1

        Args
            num_celltypes: int
                Number of fractions to create

        Returns
            fractions: numpy.ndarray
                Array of random fractions of length num_celltypes
        """
        fractions = self.rng.random(size=num_celltypes)
        fractions_sum = np.sum(fractions)
        fractions = np.divide(fractions, fractions_sum)
        return fractions

    @staticmethod
    def merge_datasets(data_dir: str = './', out_name: str = 'data.h5ad', files: list | None = None) -> None:
        """
        Merge multiple datasets of simulated samples into one h5ad file.

        Args:
            out_name: str, default = 'data.h5ad'
                Name of the merged h5ad file

            data_dir: str, default = './'
                Directory to look for datasets

            files: list | None, default = None
        """
        non_celltype_obs = ['ds', 'batch']
        if not files:
            files = glob.glob(os.path.join(data_dir, '*.h5ad'))

        logger.info(f'Merging datasets: {files} into [bold cyan]{out_name}')

        # load first file
        adata = ad.read_h5ad(files[0])

        for i in range(1, len(files)):
            adata = adata.concatenate(ad.read_h5ad(files[i]), uns_merge='same')

        combined_celltypes = list(adata.obs.columns)
        combined_celltypes = [
            x for x in combined_celltypes if x not in non_celltype_obs
        ]
        for ct in combined_celltypes:
            adata.obs[ct].fillna(0, inplace=True)

        adata.uns['cell_types'] = combined_celltypes
        adata.write(os.path.join(data_dir, out_name))
