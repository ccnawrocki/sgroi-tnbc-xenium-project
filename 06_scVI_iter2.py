# I did the following in the console, using the R kernel: 

# rm(list = ls())
# .rs.restartR(clean = T)
# reticulate::use_condaenv("ting-crc-to-liver-mets", conda = "/opt/homebrew/Caskroom/miniforge/base/bin/conda")

## Setup -----------------------------------------------------------------------
# Actual Python code begins here:

import anndata as ad
import scvi # scVI for deep learning methods
import scanpy as sc # scanpy for processing functions
# import rpy2.robjects as ro # rpy2 lets us make R function calls
import pandas as pd
import torch
import matplotlib.pyplot as plt

adata = sc.read_h5ad("sgroi-tnbc_filtered.h5ad")


##### ITERATION 2 --------------------------------------------------------------
# -- We will keep all samples together.
# -- We will integrate the data, using core as the batch variable.
# -- We will use scVI & scanpy for processing steps.

## scVI modeling
# adata.layers["counts"] = adata.X.copy().astype(int)
# adata.layers["counts"]
# 
# scvi.model.SCVI.setup_anndata(
#   adata = adata,
#   layer = "counts",
#   batch_key = "core_global"
# )
# 
# model2 = scvi.model.SCVI(
#     adata,
#     gene_likelihood="nb"
# )
# 
# model2.train(
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
# RUNTIME: 29:43

# !mkdir 06_scVI_iter2_outs
# model2.save("06_scVI_iter2_outs/scVI_model2", overwrite=True)
model2 = scvi.model.SCVI.load(adata = adata, dir_path = "06_scVI_iter2_outs/scVI_model2")

train_test_results = model2.history["elbo_train"]
train_test_results["elbo_validation"] = model2.history["elbo_validation"]
train_test_results.plot(logy=True)
plt.show()

SCVI_LATENT_KEY = "X_scVI"
latent = model2.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent

adata.write_h5ad("sgroi-tnbc_filtered.h5ad")
globals().clear()

