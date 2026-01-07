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


## Subsetting to immune cells --------------------------------------------------
adata = adata[adata.obs.celltype_level1.isin(["Lymphoid", "Myeloid", "pDC", "mregDC", "B", "Plasma"]),].copy()

## Modeling --------------------------------------------------------------------
# scvi.model.SCVI.setup_anndata(
#   adata = adata,
#   layer = "counts",
#   batch_key = "core_global"
# )
# 
# model3 = scvi.model.SCVI(
#     adata,
#     gene_likelihood="nb"
# )
# 
# model3.train(
#     check_val_every_n_epoch = 5,  # Validate less frequently (was 1)
#     plan_kwargs = {
#         "lr": 1e-3,  # Slightly higher learning rate can speed convergence
#     },
#     batch_size = 4096,
#     max_epochs = 400,
#     early_stopping = True,
#     early_stopping_patience = 20,
#     early_stopping_monitor="elbo_validation",
#     accelerator = "mps" # This is the apple silicon GPU. Use "gpu" for CUDA.
# )
# RUNTIME: 10:32

# !mkdir 09_immune_subclust_scVI_outs
# model3.save("09_immune_subclust_scVI_outs/scVI_model3", overwrite=True)
# model3 = scvi.model.SCVI.load(adata = adata, dir_path = "09_immune_subclust_scVI_outs/scVI_model3")

train_test_results = model3.history["elbo_train"]
train_test_results["elbo_validation"] = model3.history["elbo_validation"]
train_test_results.plot(logy=True)
plt.show()

SCVI_LATENT_KEY = "X_scVI"
latent = model3.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent

adata.write_h5ad("sgroi-tnbc_filtered_immune_subset.h5ad")
globals().clear()

