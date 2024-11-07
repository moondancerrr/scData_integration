# Single-cell RNA-seq Data Integration Workflow

## Description

This workflow provides a complete data integration and analysis pipeline for single-cell RNA-seq data using a heart cell atlas dataset. The process includes normalization, batch correction, and integration using various methods (scVI, scANVI, BBKNN, and Seurat). Metrics are then calculated to benchmark and evaluate the quality of data integration across different methods.

## Requirements

### Software

- **Python** (version 3.8+)
- **R** (version 4.0+)
  
### Python Packages:

- `scanpy`
- `scvi-tools`
- `bbknn`
- `scib`
- `numpy`
- `pandas`
- `matplotlib`
- `rpy2`
- `anndata2ri`

You can install these packages via pip:

```bash
pip install scanpy scvi-tools bbknn scib numpy pandas matplotlib rpy2 anndata2ri
```

### R Packages:

- `Seurat`

You can install the Seurat package in R:

```R
install.packages("Seurat")
```

## Running the Workflow
1. Load and Preprocess Data

The script begins by loading the Heart Cell Atlas dataset and performing preprocessing steps, including gene filtering, normalization, and selection of highly variable genes (HVGs).

2. Integration Methods

- scVI Integration: Applies scVI for integration without cell labels.
- scANVI Integration: Extends scVI with cell labels for semi-supervised learning.
- BBKNN Integration: Utilizes batch-balanced k-nearest neighbors for batch effect correction.
- Seurat Integration: Leverages Seurat's Mutual Nearest Neighbors (MNN) in R for integration.

3. Evaluation and Visualization

The workflow calculates integration metrics (e.g., ASW, PCR) and visualizes results through UMAPs, bar plots, and scatter plots, allowing you to evaluate each integration method for your dataset.

## Running the Script

To execute the script, simply run it in a Python environment configured with the necessary Python and R dependencies. Each section of the script performs specific tasks as outlined in the comments.

```bash
python integration_workflow.py ```
