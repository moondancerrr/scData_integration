import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

# Python packages
import scanpy as sc
import scvi
import bbknn
import scib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

# R interface
from rpy2.robjects import pandas2ri
from rpy2.robjects import r
import rpy2.robjects as ro
import rpy2.rinterface_lib.callbacks
import anndata2ri

# Activate converters for AnnData and pandas
pandas2ri.activate()
anndata2ri.activate()

# Load the Seurat library in R
ro.r('''library(Seurat)''')

# Load the heart cell atlas dataset
adata = scvi.data.heart_cell_atlas_subsampled()
￼
# Filter out genes with low counts and outliers
sc.pp.filter_genes(adata, min_counts=3)

# Preserve raw counts in a layer and normalize data
adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify highly variable genes (HVGs)
# Perform feature selection, to reduce the number of features (genes in this case)
# For scVI, it's recommended anywhere from 1,000 to 10,000 HVGs, but it will be context-dependent.
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key="cell_source",
)

# Save normalized counts in a layer
adata.layers["logcounts"] = adata.X

# Define the keys for cell type and batch
label_key = "cell_type"
batch_key = "cell_source"

# Display value counts for each batch
adata.obs[batch_key].value_counts()

# compute PCA and UMAP
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Define custom colors for batches and visualize UMAP
adata.uns[batch_key + "_colors"] = [
      "#001524",
      "#15616d",
      "#ffecd1",
      "#ff7d00"
      ] # Set custom colours for batches
sc.pl.umap(adata, color=[label_key, batch_key], wspace=1)

# Plot bar plot for number of batches with highly variable genes
n_batches = adata.var["highly_variable_nbatches"].value_counts()
ax = n_batches.plot(kind="bar")
plt.show()

### Unintegrated data ###
# Create a new AnnData object with only HVGs as an untegrated version of the data
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

### Variational autoencoder (VAE) based integration  without cell labels ###
# Create a copy of HVG AnnData for SCVI
adata_scvi = adata_hvg.copy()

# Set up SCVI with batch and donor as categorical covariates￼
# Alert scvi-tools to the locations of various matrices inside the anndata
# There is a “cell_source” categorical covariate, and within each “cell_source”, multiple “donors”, “gender” and “age_group”. There are also two continuous covariates we’d like to correct for: “percent_mito” and “percent_ribo”
scvi.model.SCVI.setup_anndata(adata_scvi,
                              layer="counts", 
                              categorical_covariate_keys=["cell_source", "donor"],
                              continuous_covariate_keys=["percent_mito", "percent_ribo"],
                              )

# Train SCVI model
model_scvi = scvi.model.SCVI(adata_scvi)

max_epochs_scvi = np.min([round((20000 / adata.n_obs) * 400), 400])

model_scvi.train()

# Get latent representation which is a multi-dimensional embedding where the batch effects have been removed that can be used in a similar way to how we use PCA dimensions when analysing a single dataset
adata_scvi.obsm["X_scVI"] = model_scvi.get_latent_representation()
￼
# Calculate a new UMAP embedding but instead of finding nearest neighbors in PCA space, we start with the corrected representation from scVI
sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
sc.tl.umap(adata_scvi)
sc.pl.umap(adata_scvi, color=[label_key, batch_key], wspace=1)
￼
### VAE integration using cell labels ###

# Normally we would need to run scVI first but we have already done that here
# model_scvi = scvi.model.SCVI(adata_scvi) etc.
# Create a SCANVI model from SCVI model, providing label information
model_scanvi = scvi.model.SCANVI.from_scvi_model(
    model_scvi, labels_key=label_key, unlabeled_category="unlabelled"
)

model_scanvi.view_anndata_setup()

# Train SCANVI model
# This scANVI model object is very similar to what we saw before for scVI. So the training now is much fewer than before as we are just refining the scVI model, rather than training a whole network from scratch.
max_epochs_scanvi = int(np.min([10, np.max([2, round(max_epochs_scvi / 3.0)])]))
model_scanvi.train(max_epochs=max_epochs_scanvi)

# Obtain SCANVI latent representation and compute UMAP
adata_scanvi = adata_scvi.copy()
adata_scanvi.obsm["X_scANVI"] = model_scanvi.get_latent_representation()
sc.pp.neighbors(adata_scanvi, use_rep="X_scANVI")
sc.tl.umap(adata_scanvi)
sc.pl.umap(adata_scanvi, color=[label_key, batch_key], wspace=1)

### Graph based Integration ###
# Set neighbors within batch
neighbors_within_batch = 25 if adata_hvg.n_obs > 100000 else 3

# Copy HVG data and apply BBKNN
adata_bbknn = adata_hvg.copy()

adata_bbknn.X = adata_bbknn.layers["logcounts"].copy()
sc.pp.pca(adata_bbknn)

bbknn.bbknn(
    adata_bbknn, batch_key=batch_key, neighbors_within_batch=neighbors_within_batch
)

# Compute UMAP for BBKNN
sc.tl.umap(adata_bbknn)
sc.pl.umap(adata_bbknn, color=[label_key, batch_key], wspace=1)

### Linear embedding integration using Mutual Nearest Neighbors (MNN) ###
# Prepare data for Seurat integration
adata_seurat = adata_hvg.copy()
# Convert categorical columns to strings
adata_seurat.obs[batch_key] = adata_seurat.obs[batch_key].astype(str)
adata_seurat.obs[label_key] = adata_seurat.obs[label_key].astype(str)

# Delete uns as this can contain arbitrary objects which are difficult to convert
del adata_seurat.uns

# Activate the automatic conversion
anndata2ri.activate()

# Transfer the AnnData object `adata_seurat` to R
ro.globalenv["adata_seurat"] = adata_seurat

# Seurat integration functions require a list of objects added by SplitObject() function
ro.r('''
library(Seurat)
seurat <- as.Seurat(adata_seurat, counts = "counts", data = "logcounts")
batch_list <- SplitObject(seurat, split.by = "cell_source")
# Use this list to find anchors for each pair of datasets.
anchors <- FindIntegrationAnchors(batch_list, anchor.features = rownames(seurat))
# Use the anchors to compute a transformation that maps one dataset onto another
integrated <- IntegrateData(anchors)
# Extract the integrated expression matrix
integrated_expr <- GetAssayData(integrated)
# Make sure the rows and columns are in the same order as the original object
integrated_expr <- integrated_expr[rownames(seurat), colnames(seurat)]
# Transpose the matrix to AnnData format
integrated_expr <- t(integrated_expr)
print(integrated_expr[1:10, 1:10])
''')

# Transfer integrated_expr back to Python
integrated_expr = ro.r("integrated_expr")
# Store the corrected expression matrix as a layer in our AnnData object
adata_seurat.X = integrated_expr
adata_seurat.layers["seurat"] = integrated_expr

# Set custom colors for Seurat batch and visualize UMAP
adata_seurat.uns[batch_key + "_colors"] = [
       "#001524",
       "#15616d",
       "#ffecd1",
       "#ff7d00",
       ]

sc.tl.pca(adata_seurat)
sc.pp.neighbors(adata_seurat)
sc.tl.umap(adata_seurat)
sc.pl.umap(adata_seurat, color=[label_key, batch_key], wspace=1)

# Compare all methods to evaluate the quality of integration 
metrics_scvi = scib.metrics.metrics_fast(
    adata, adata_scvi, batch_key, label_key, embed="X_scVI"
)
metrics_scanvi = scib.metrics.metrics_fast(
    adata, adata_scanvi, batch_key, label_key, embed="X_scANVI"
)
metrics_bbknn = scib.metrics.metrics_fast(adata, adata_bbknn, batch_key, label_key)
metrics_seurat = scib.metrics.metrics_fast(adata, adata_seurat, batch_key, label_key)
metrics_hvg = scib.metrics.metrics_fast(adata, adata_hvg, batch_key, label_key)

# Concatenate metrics results
metrics = pd.concat(
    [metrics_scvi, metrics_scanvi, metrics_bbknn, metrics_seurat, metrics_hvg],
    axis="columns",
)
# Set methods as column names
metrics = metrics.set_axis(
    ["scVI", "scANVI", "BBKNN", "Seurat", "Unintegrated"], axis="columns"
)
# Select only the fast metrics
metrics = metrics.loc[
    [
        "ASW_label",
        "ASW_label/batch",
        "PCR_batch",
        "isolated_label_silhouette",
        "graph_conn",
        "hvg_overlap",
    ],
    :,
]
# Transpose so that metrics are columns and methods are rows
metrics = metrics.T
# Remove the HVG overlap metric because it's not relevant to embedding outputs
metrics = metrics.drop(columns=["hvg_overlap"])

metrics.style.background_gradient(cmap="Blues")
￼
metrics_scaled = (metrics - metrics.min()) / (metrics.max() - metrics.min())
metrics_scaled.style.background_gradient(cmap="Blues")

metrics_scaled["Batch"] = metrics_scaled[
    ["ASW_label/batch", "PCR_batch", "graph_conn"]
].mean(axis=1)
metrics_scaled["Bio"] = metrics_scaled[["ASW_label", "isolated_label_silhouette"]].mean(
    axis=1
)
metrics_scaled.style.background_gradient(cmap="Blues")

fig, ax = plt.subplots()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
metrics_scaled.plot.scatter(
    x="Batch",
    y="Bio",
    c=range(len(metrics_scaled)),
    ax=ax,
)

for k, v in metrics_scaled[["Batch", "Bio"]].iterrows():
    ax.annotate(
        k,
        v,
        xytext=(6, -3),
        textcoords="offset points",
        family="sans-serif",
        fontsize=12,
    )

plt.show()

# Get an overall score for each method we can combine the two summary scores
# The scIB paper suggests a weighting of 40% batch correction and 60% biological conservation
metrics_scaled["Overall"] = 0.4 * metrics_scaled["Batch"] + 0.6 * metrics_scaled["Bio"]
metrics_scaled.style.background_gradient(cmap="Blues")
metrics_scaled.plot.bar(y="Overall")
plt.show()￼

