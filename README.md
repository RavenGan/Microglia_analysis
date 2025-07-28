# Microglia_analysis

This repository contains R scripts for processing and visualizing single-cell RNA-seq (scRNA-seq) data of microglia at embryonic day 14 (E14). The workflow includes raw data processing, quality control, and visualization.

## Data Source

The data analyzed in this project are publicly available from the Gene Expression Omnibus (GEO):

- **Accession:** [GSE121654](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121654)  
- **Reference:**  
  Hammond, T. R., Dufort, C., Dissing-Olesen, L., Giera, S., Young, A., Wysoker, A., ... & Stevens, B. (2019).  
  *Single-cell RNA sequencing of microglia throughout the mouse lifespan and in the injured brain reveals complex cell-state changes*.  
  **Immunity**, 50(1), 253–271. https://doi.org/10.1016/j.immuni.2018.11.004

## Directory Structure

- `./code/`: Contains all R scripts used in the analysis.

## Scripts Overview

1. **`Dat_microglia_E14_raw.R`**  
   Merges all raw input data and generates the count matrices and metadata required for downstream analysis.

2. **`Dat_microglia_E14.R`**  
   Performs preprocessing steps on the scRNA-seq data, including normalization, scaling, dimensionality reduction, clustering, and annotation.

3. **`Dat_microglia_E14_HeatmapVioUmap.R`**  
   Generates key visualizations such as heatmaps, violin plots, and feature plots to highlight gene expression patterns and cluster characteristics.

4. **`Dat_microglia_E14_Violin.R`**  
   Produces additional violin plots for in-depth comparison of gene expression across clusters or conditions.
