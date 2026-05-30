# Resilience of the human T-cell compartment after early-life thymectomy

This repository contains scripts used to generate T-cell receptor sequencing data and figures from the paper "Resilience of the human T-cell compartment after early-life thymectomy" (Figs. 3-4 and Supplementary Figs. 6-10).

Sample metadata linking FASTQ file identifiers to sample information can be found in the `data/metadata` directory.

## To reproduce the figures and analyses in the paper

### Prepare data

Download all `fastq.gz` files from [Dataverse](https://doi.org/10.34894/HQCNGM) and place them in `./data/fastq/`.

### Steps

0. Some samples were sequenced two or more times. To retain as much information as possible, first merge the `fastq.gz` files from the same samples.

   Run:

   ```bash
   python ./src/0.mergeFASTQ/mergeFASTQ_GenomescanIDs.py
   ```

1. TCR alpha and beta alignment uses [MiXCR](https://mixcr.com/mixcr/about/). The pipeline can be run with [Snakemake](https://snakemake.readthedocs.io/en/stable/) from the `src/1.align/` directory.

   The `config.yaml` file contains machine-specific parameters, such as the available memory for running jobs.

   ```bash
   cd src/1.align/
   snakemake
   ```

2. UMAP figures (Fig. 1D and Fig. 3A) can be generated using `./src/2.UMAP/TxAgeing_UMAPfigures_code.R`.

3. Fig. 4 can be generated using `src/3.clone_per_cell_or_UMI/clone_per_cell_plots.py`.

4. Clonotype distribution bar plots per sample (Supplementary Fig. 8) can be generated using `src/4.clone_distribution_bar/clone_dist.R`.

5. Sample metric bar plots for read, UMI, and cell counts (Supplementary Figs. 6 and 7) can be generated using `src/5.explore/plot_sample_metrics.py`.

6. Sample metric scatter plots (Supplementary Figs. 9 and 10) can be generated using `src/5.explore/plot_sample_scatterplots.py`.
