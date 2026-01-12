# I did the following in the console, using the R kernel: 

# rm(list = ls())
# .rs.restartR(clean = T)

## Setup -----------------------------------------------------------------------
# Actual Python code begins here:

import anndata as ad
import scvi # scVI for deep learning methods
import scanpy as sc # scanpy for processing functions
import pandas as pd
import torch
import matplotlib.pyplot as plt

adata = sc.read_h5ad("sgroi-tnbc_filtered_immune_subset.h5ad")
sc.pl.umap(adata, color = "celltype_level2")
plt.subplots_adjust(right=0.65)
plt.show()

## Subsetting to T and NK cells ------------------------------------------------
adata = adata[adata.obs.celltype_level2.isin(["T and NK"]),].copy()

## Subsetting to genes that matter ---------------------------------------------
# Basically, we need to ignore genes that are likely to be in other cell types 
# and which will only be expressed in this setting as contamination.
# Claude Sonnet 4.5 helped me identify these genes, but Blake can likely confirm.
genes_oi = ~adata.var_names.isin([
  
  # Structural/Stromal
  "COL1A1", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "MEG3", "CXCL14"
  "LUM", "DCN", "SPARC", "SPARCL1", "FN1", "VCAN", "FAP",
  "SERPINE1", "TGFBI", "SPON2", "RGS5", "SEMA3A", "ACTA2", "CAV1",
  
  # Endothelial
  "PLVAP", "FLT1", "IGFBP7", "JAM2",
  
  # Epithelial/Tumor
  "EPCAM", "KRT77", "KRT80", "ESR1", "PGR", "FOXA1",
  "MUC5AC", "TFF3", "REG4", "CDX1", "CDX2",
  "CEACAM6", "CEACAM8", "SOX9", "SOX17",
  "ERBB2", "ERBB3", "ERBB4", "EGFR", "MET", "AXL",
  
  # Immunoglobulins/B cell
  "IGHE", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM",
  "IGKC", "IGLC3", "JCHAIN",
  "MS4A1", "CD19", "PAX5", "SDC1", "MZB1", "CD79A", "CD79B",
  "TCL1A", "BANK1", "FCMR", "FCER2", "CD37", "FCRL3",
  
  # Other
  "XBP1", "CFC1", "PDGFRA", 
  "MMP1", "MMP2", "MMP9", "MMP12", "MMP14"
  
])
adata = adata[:, genes_oi].copy()

## Subsetting to high quality cells --------------------------------------------
import numpy as np
adata.obs["counts_post_subset"] = np.array(adata.layers["counts"].sum(axis=1)).flatten()
cells_keep = adata.obs["counts_post_subset"] >= 20
np.mean(cells_keep) # ~99%
adata = adata[cells_keep,].copy()

## Modeling --------------------------------------------------------------------
scvi.model.SCVI.setup_anndata(
  adata = adata,
  layer = "counts",
  batch_key = "core_global"
)

model4 = scvi.model.SCVI(
    adata,
    gene_likelihood="nb"
)

model4.train(
    check_val_every_n_epoch = 5,  # Validate less frequently (was 1)
    plan_kwargs = {
        "lr": 1e-3,  # Slightly higher learning rate can speed convergence
    },
    batch_size = 4096,
    max_epochs = 400,
    early_stopping = True,
    early_stopping_patience = 20,
    early_stopping_monitor="elbo_validation",
    accelerator = "mps" # This is the apple silicon GPU. Use "gpu" for CUDA.
)
# RUNTIME: 3:38

!mkdir 11_tcell_subclust_scVI_outs
model4.save("11_tcell_subclust_scVI_outs/scVI_model4", overwrite=True)
model4 = scvi.model.SCVI.load(adata = adata, dir_path = "11_tcell_subclust_scVI_outs/scVI_model4")

train_test_results = model4.history["elbo_train"]
train_test_results["elbo_validation"] = model4.history["elbo_validation"]
train_test_results.plot(logy=True)
plt.show()

SCVI_LATENT_KEY = "X_scVI"
latent = model4.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent

adata.write_h5ad("sgroi-tnbc_filtered_tcell_subset.h5ad")
globals().clear()

