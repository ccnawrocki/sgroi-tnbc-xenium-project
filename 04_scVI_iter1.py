## Using a pre-made conda environment in our project ---------------------------
# I did the following in the console, using the R kernel: 

# rm(list = ls())
# .rs.restartR(clean = T)
# reticulate::conda_list() # Shows all the environments and their paths
# reticulate::use_condaenv("ting-crc-to-liver-mets", conda = "/opt/homebrew/Caskroom/miniforge/base/bin/conda")

## Setup -----------------------------------------------------------------------
# Actual Python code begins here:

import anndata as ad
import scanpy as sc # scanpy for processing functions
import rpy2.robjects as ro # rpy2 lets us make R function calls
import pandas as pd
import torch
import matplotlib.pyplot as plt
import scvi # scVI for deep learning methods

adata = sc.read_h5ad("sgroi-tnbc.h5ad")


##### ITERATION 1 --------------------------------------------------------------
# -- We will keep all samples together.
# -- We will integrate the data, using core as the batch variable.
# -- We will use scVI & scanpy for processing steps.
# -- Next, we will identify noise clusters and filter them out.

## scVI modeling
adata.layers["counts"] = adata.X.copy().astype(int)
adata.layers["counts"]

adata.obs.loc[:,"core_global"] = adata.obs.slide.str.cat(adata.obs.core, sep="_").astype("category")

scvi.model.SCVI.setup_anndata(
  adata = adata,
  layer = "counts",
  batch_key = "core_global"
)

model1 = scvi.model.SCVI(
    adata,
    gene_likelihood="nb"
)

# Claude helped me with the nitty-gritty for speeding up the training:
# model1.train(
#     check_val_every_n_epoch = 5,  # Validate less frequently (was 1)
#     plan_kwargs={
#         "lr": 1e-3,  # Slightly higher learning rate can speed convergence
#     },
#     batch_size=4096,
#     max_epochs = 400,
#     early_stopping = True,
#     early_stopping_patience = 20,
#     early_stopping_monitor="elbo_validation",
#     accelerator = "mps" # This is the apple silicon GPU. Use "gpu" for CUDA.
# )
# Epoch 392/400:  98%|█████████▊| 392/400 [31:29<00:38,  4.82s/it, v_num=1, train_loss=273]
# Monitored metric elbo_validation did not improve in the last 20 records. Best score: 271.935. Signaling Trainer to stop.

# !mkdir 04_scVI_iter1_outs
# model1.save("04_scVI_iter1_outs/scVI_model1", overwrite=True)
model1 = scvi.model.SCVI.load(adata = adata, dir_path = "04_scVI_iter1_outs/scVI_model1")

train_test_results = model1.history["elbo_train"]
train_test_results["elbo_validation"] = model1.history["elbo_validation"]
train_test_results.plot(logy=True)
plt.show()

SCVI_LATENT_KEY = "X_scVI"
latent = model1.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent

adata.write_h5ad("sgroi-tnbc.h5ad")
globals().clear()

