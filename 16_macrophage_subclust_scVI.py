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

sc.pl.umap(adata, color = ["TGFBI", "CD68", "CD14", "CD163"], layer = "lognorm", ncols = 2, frameon = False)

## Subsetting to Macrophages ---------------------------------------------------
adata = adata[adata.obs.celltype_level2.isin(["Macrophage1", "Macrophage2"]),].copy()

## Subsetting to genes that matter ---------------------------------------------
genes_oi = ~adata.var_names.isin([
  
  # Structural/Stromal
  "COL1A1", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "MEG3", "CXCL14",
  "LUM", "DCN", "SPARC", "SPARCL1", "FN1", "VCAN", "FAP",
  "SERPINE1", "SPON2", "RGS5", "SEMA3A", "ACTA2", "CAV1", # I removed TGFBI from this list, since macrophages can express it.
  
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

model5 = scvi.model.SCVI(
    adata,
    gene_likelihood="nb"
)

model5.train(
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
# RUNTIME: 2:57

!mkdir 16_macrophage_subclust_scVI_outs
model5.save("16_macrophage_subclust_scVI_outs/scVI_model5", overwrite=True)
model5 = scvi.model.SCVI.load(adata = adata, dir_path = "16_macrophage_subclust_scVI_outs/scVI_model5")

train_test_results = model5.history["elbo_train"]
train_test_results["elbo_validation"] = model5.history["elbo_validation"]
train_test_results.plot(logy=True)
plt.show()

SCVI_LATENT_KEY = "X_scVI"
latent = model5.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent

adata.write_h5ad("sgroi-tnbc_filtered_macrophage_subset.h5ad")
globals().clear()

