# Spontaneous movement

## How to run this analysis

The data should be available under the `raw` directory in this repository root (create a symbolic link).

1. Install conda and mamba. Create the _snakemake_ (`env/snakemake.yml`) conda environment
   (required to run the analysis) and the _spontaneous-movement-mne_ (`env/mne.yml`) conda environment
   (required to generate the reports).

```shell
mamba env create -f env/snakemake.yml
mamba env create -f env/mne.yml
```

2. Activate the snakemake environment.

```shell
conda activate snakemake
```

3. Run ``snakemake``:

```shell
snakemake --conda-frontend mamba --use-conda --cores 1 -p all
```

## TODO

- [x] Go through the movement detection and update / make sure everything works as expected.
- [ ] Create reports for all analysis
- [ ] Update the pipeline to use .csv whenever possible and update integration with R
- [x] Generate figures in the pipeline
- [ ] Create reports for all performed analyses
- [ ] Does the *frequency* of action potentials increase during movement/rest? during various phases of the movement onset?
- [ ] wavelet transform

- [ ] average per cell (figure 2, figure 3b, figure 4b,c)
- [ ] redo fig 4 -> fig 2
