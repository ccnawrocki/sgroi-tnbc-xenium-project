# In the console, using the R kernel: 
# rm(list = ls())
# .rs.restartR(clean = T)
# reticulate::py_require("anndata")

import anndata as ad
import glob

!ls *.h5ad
adata_paths = glob.glob("*.h5ad")

adata = {}
for p in adata_paths:
  nm = p.replace(".h5ad", "")
  adata[nm] = ad.read_h5ad(p)

adata = ad.concat(adata, index_unique = "_")
adata.shape
# (1353425, 480)

# I manually screened cores that were unsuitable for this study. 
# If a core was a control core, I considered it unsuitable.
# If a core had very little tissue and was also the only core for a patient, 
# then I considered it unsuitable.

# I remove the unsuitable cores here: 
tokeep = ~(((adata.obs["slide"] == "tma1") & (adata.obs["core"].isin(["A1"]))) | ((adata.obs["slide"] == "tma2") & (adata.obs["core"].isin(["A1", "B1"]))) | ((adata.obs["slide"] == "tma3") & (adata.obs["core"].isin(["A1", "B1"]))) | ((adata.obs["slide"] == "tma4") & (adata.obs["core"].isin(["A1", "B1", "B2"]))))
adata = adata[tokeep,]
adata.shape
# (1086571, 480)

# We end up with >1M cells, before QC.
adata.write_h5ad("sgroi-tnbc.h5ad")

# Now that I have all the data organized, I will move the individual anndata objects to storage. 

globals().clear()

