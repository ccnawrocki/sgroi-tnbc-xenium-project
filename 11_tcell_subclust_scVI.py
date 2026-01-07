# I did the following in the console, using the R kernel: 

# rm(list = ls())
# .rs.restartR(clean = T)
# reticulate::use_condaenv("ting-crc-to-liver-mets", conda = "/opt/homebrew/Caskroom/miniforge/base/bin/conda")

## Setup -----------------------------------------------------------------------
# Actual Python code begins here:

import anndata as ad
import scvi # scVI for deep learning methods
import scanpy as sc # scanpy for processing functions
import pandas as pd
import torch
import matplotlib.pyplot as plt

adata = sc.read_h5ad("sgroi-tnbc_filtered.h5ad")
sc.pl.umap(adata, color = "celltype_level1")
plt.subplots_adjust(right=0.75)
plt.show()

## Subsetting to immune cells --------------------------------------------------
adata = adata[adata.obs.celltype_level1.isin(["Lymphoid"]),].copy()

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
# RUNTIME: 4:50

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

