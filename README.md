# Integration of transcriptomics data for _Dunaliella salina_, and _Nannochloropsis gaditana_.

## Introduction

This repository contains the code, data, and results used to integrate transcriptomics data into the genome-scale metabolic models of _Dunaliella salina_,
and _Nannochloropsis gaditana_.

## Data

The data folder contains the gene normalized counts (getmm.tsv), the GEM for each species, the growth media, and a map of reactions to pathways in json format.

## Code

The code used in this study is available in the `code` directory. The code is organized as follows:

- 'topological_validation': Jupyter notebook used to validate the topological properties of the model.
- 'basic_simulations': Jupyter notebook used to simulate the model under different conditions.
- 'light_evaluation': Jupyter notebook used to evaluate the light reactions of the model.
- 'light_absorption': Python script used to evaluate the photosynthetic properties of the models at different light intensities.
- 'utils': Python script containing utility functions used in the study.


To reproduce most results, the packages listed in `requirements.txt` are required.

## Results

The results of this study are available in the `results` directory. For each case study, the results are organized as follows:

- 'preprocessing': Contains the density plots of the transcriptomics data.
- 'deg': Contains the differentially expressed genes (DEGs) and gene set enrichment analysis (GSEA) for each condition.
- 'integration': Results of the integration of transcriptomics data into the GEM, with each method and threshold.
- 'dfa': Contains the results of the differential flux analysis (DFA) for each condition-control pair.

The samples obtained from flux sampling are available upon request.