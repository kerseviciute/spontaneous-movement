# Spontaneous movement

## Run the analysis

```shell
snakemake --cores 1 --conda-frontend mamba --use-conda --rerun-triggers mtime
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
