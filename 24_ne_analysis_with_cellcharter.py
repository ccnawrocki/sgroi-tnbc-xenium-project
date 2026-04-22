!which python
#/Users/ccn22/Documents/projects/xenium/sgroi-tnbc-project/.venv/bin/python

import anndata as ad
import scanpy as sc
import numpy as np

adata = ad.read_h5ad("sgroi-tnbc_filtered.h5ad")

adata.obsm 
#AxisArrays with keys: X_scVI, X_umap

#!uv pip install cellcharter
import cellcharter as cc

# Need to add the spatial coordinates to obsm/spatial
adata.obsm["spatial"] = np.array(adata.obs[["x_global", "y_global"]])

# Make some of this data categorical for cellcharter
adata.obs["celltype_final_amended_1"] = adata.obs[["celltype_final_amended_1"]].astype("category")
adata.obs["banksy_niche"] = adata.obs[["banksy_niche"]].astype("category")
adata.obs["HOXB13_group"] = adata.obs[["HOXB13_group"]].astype("category")
adata.obs["patient"] = adata.obs[["patient"]].astype("category")
adata.obs["core_global"] = adata.obs[["core_global"]].astype("category")

# Constructing the neighbors
#!uv pip install squidpy
import squidpy as sq
sq.gr.spatial_neighbors(adata, library_key="core_global", delaunay=True)
cc.gr.remove_long_links(adata, distance_percentile=98)

# Quick plot
import matplotlib.pyplot as plt
sq.pl.spatial_scatter(
    adata, 
    shape = None,
    color=['banksy_niche'], 
    library_key='core_global', 
    img=None, 
    title=["tma4_E3"],
    size=0.25,
    connectivity_key='spatial_connectivities',
    edges_width=0.2,
    #legend_loc=None,
    library_id=["tma4_E3"]
    )
plt.show()


###### NICHE COMPARISON ########################################################

# HOXB13+
hoxpos = adata[(adata.obs['HOXB13_group'] == 'HOXB13+') & (adata.obs['banksy_niche'] != 'singleton')]
cc.gr.nhood_enrichment(
    hoxpos,
    cluster_key='banksy_niche', 
    only_inter=False
)
cc.pl.nhood_enrichment(
    hoxpos,
    cluster_key='banksy_niche',
    annotate=True,
    figsize=(6,6),
    fontsize=10,
)
plt.show()

# HOXB13-
hoxneg = adata[(adata.obs['HOXB13_group'] == 'HOXB13-') & (adata.obs['banksy_niche'] != 'singleton')]
cc.gr.nhood_enrichment(
    hoxneg,
    cluster_key='banksy_niche', 
    only_inter=False
)
cc.pl.nhood_enrichment(
    hoxneg,
    cluster_key='banksy_niche',
    annotate=True,
    figsize=(6,6),
    fontsize=10,
)
plt.show()

# Difference
tmp = adata[(adata.obs['banksy_niche'] != 'singleton')]
cc.gr.diff_nhood_enrichment(
    tmp, 
    cluster_key='banksy_niche',
    only_inter=False,
    condition_key='HOXB13_group',
    library_key='patient',
    pvalues=True,
    n_jobs=4,
    n_perms=1000, 
    batch_size=250
)
cc.pl.diff_nhood_enrichment(
    tmp,
    cluster_key='banksy_niche',
    condition_key='HOXB13_group',
    condition_groups=['HOXB13+', 'HOXB13-'],
    annotate=True,
    figsize=(6,6),
    significance=0.05,
    fontsize=10
)
plt.show()
#plt.savefig("24_ne_analysis_with_cellcharter_outs/banksy_cc_pos_minus_neg.pdf", dpi=96)

# Saving
#!mkdir 24_ne_analysis_with_cellcharter_outs
hoxpos.uns['banksy_niche_nhood_enrichment']['enrichment'].to_csv("24_ne_analysis_with_cellcharter_outs/banksy_cc_hoxpos.csv")
hoxneg.uns['banksy_niche_nhood_enrichment']['enrichment'].to_csv("24_ne_analysis_with_cellcharter_outs/banksy_cc_hoxneg.csv")
tmp.uns['banksy_niche_HOXB13_group_diff_nhood_enrichment']["HOXB13+_HOXB13-"]['enrichment'].to_csv("24_ne_analysis_with_cellcharter_outs/banksy_cc_pos_minus_neg_effects.csv")
tmp.uns['banksy_niche_HOXB13_group_diff_nhood_enrichment']["HOXB13+_HOXB13-"]['pvalue'].to_csv("24_ne_analysis_with_cellcharter_outs/banksy_cc_pos_minus_neg_pvalues.csv")

# Shape metrics
cc.gr.connected_components(tmp, cluster_key='banksy_niche', min_cells=20)
cc.tl.boundaries(tmp, min_hole_area_ratio=0.1)

cc.tl.linearity(tmp)
cc.tl.curl(tmp)
cc.tl.elongation(tmp)
cc.tl.purity(tmp, library_key="core_global")

cc.pl.shape_metrics(
    tmp,
    cluster_key='banksy_niche',
    condition_key='HOXB13_group',
    condition_groups=['HOXB13+', 'HOXB13-'],
    cluster_groups=["6"], 
    metrics=["linearity", "elongation", "curl", "purity"],
    figsize=(10,6),
    fontsize=6, 
    #save="24_ne_analysis_with_cellcharter_outs/banksy_niche6_shapes_cc_pos_vs_neg.pdf"
)
plt.show()

# Saving
import pandas as pd
shape_results = pd.DataFrame(tmp.uns["shape_component"])
tmp.obs.groupby("component")[["banksy_niche", "core_global"]].first().merge(
    shape_results, left_index=True, right_index=True
).to_csv("24_ne_analysis_with_cellcharter_outs/banksy_cc_niche_components.csv", index=True)


###### CELLTYPE COMPARISON #####################################################

tmp = adata[(adata.obs['banksy_niche'] != 'singleton')]
tmp.obs["celltype"] = tmp.obs["celltype_level1"].astype("object")
tmp.obs.loc[(tmp.obs["celltype_level2"] != "NA"), "celltype"] = tmp.obs.loc[(tmp.obs["celltype_level2"] != "NA"), "celltype_level2"].astype("object")
tmp.obs["celltype"] = tmp.obs["celltype"].astype("category")

# HOXB13+
hoxpos = tmp[(tmp.obs['HOXB13_group'] == 'HOXB13+') & (tmp.obs['banksy_niche'] != 'singleton')]
cc.gr.nhood_enrichment(
    hoxpos,
    cluster_key='celltype', 
    only_inter=False
)
cc.pl.nhood_enrichment(
    hoxpos,
    cluster_key='celltype',
    annotate=True,
    figsize=(6,6),
    fontsize=10,
)
plt.show()

# HOXB13-
hoxneg = tmp[(tmp.obs['HOXB13_group'] == 'HOXB13-') & (tmp.obs['banksy_niche'] != 'singleton')]
cc.gr.nhood_enrichment(
    hoxneg,
    cluster_key='celltype', 
    only_inter=False
)
cc.pl.nhood_enrichment(
    hoxneg,
    cluster_key='celltype',
    annotate=True,
    figsize=(6,6),
    fontsize=10,
)
plt.show()

# Difference
cc.gr.diff_nhood_enrichment(
    tmp, 
    cluster_key='celltype',
    only_inter=False,
    condition_key='HOXB13_group',
    library_key='patient',
    pvalues=True,
    n_jobs=8,
    n_perms=1000, 
    batch_size=150
)
cc.pl.diff_nhood_enrichment(
    tmp,
    cluster_key='celltype',
    condition_key='HOXB13_group',
    condition_groups=['HOXB13+', 'HOXB13-'],
    annotate=False,
    figsize=(6,6),
    significance=0.05,
    fontsize=10
)
plt.show()
#plt.savefig("24_ne_analysis_with_cellcharter_outs/celltype_cc_pos_minus_neg.pdf", dpi=96)

# Saving
hoxpos.uns['celltype_nhood_enrichment']['enrichment'].to_csv("24_ne_analysis_with_cellcharter_outs/celltype_cc_hoxpos.csv")
hoxneg.uns['celltype_nhood_enrichment']['enrichment'].to_csv("24_ne_analysis_with_cellcharter_outs/celltype_cc_hoxneg.csv")
tmp.uns['celltype_HOXB13_group_diff_nhood_enrichment']["HOXB13+_HOXB13-"]['enrichment'].to_csv("24_ne_analysis_with_cellcharter_outs/celltype_cc_pos_minus_neg_effects.csv")
tmp.uns['celltype_HOXB13_group_diff_nhood_enrichment']["HOXB13+_HOXB13-"]['pvalue'].to_csv("24_ne_analysis_with_cellcharter_outs/celltype_cc_pos_minus_neg_pvalues.csv")


################################################################################
################################################################################
################################################################################

###### EXPERIMENTAL ############################################################

# I will try using cellcharter to identify niches, instead of Banksy.
globals().clear()

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import squidpy as sq
import cellcharter as cc

adata = ad.read_h5ad("sgroi-tnbc_filtered.h5ad")
adata.obsm["spatial"] = np.array(adata.obs[["x_global", "y_global"]])
adata.obs["core_global"] = adata.obs[["core_global"]].astype("category")
sq.gr.spatial_neighbors(adata, library_key="core_global", delaunay=True)
cc.gr.remove_long_links(adata, distance_percentile=98)

cc.gr.aggregate_neighbors(adata, n_layers=3, sample_key="core_global", use_rep="X_scVI")

# import warnings
# warnings.filterwarnings('ignore')
# 
# model_params = {
#   'random_state': 2001,
#   'trainer_params': {'accelerator': 'cpu', 'enable_progress_bar': False},
#   #'batch_size': 262144, # <-- probably an acceptable size, if data does not fit into memory
#   }
# autok = cc.tl.ClusterAutoK(
#     n_clusters=(2,20),
#     convergence_tol=0.001, 
#     max_runs=10,
#     model_params=model_params,
# )
# autok.fit(adata, use_rep='X_cellcharter')
# autok.save("24_ne_analysis_with_cellcharter_outs/autok_model")

import os
os.remove("24_ne_analysis_with_cellcharter_outs/autok_model/best_models/.DS_Store") # This is jokes!
autok = cc.tl.ClusterAutoK.load('24_ne_analysis_with_cellcharter_outs/autok_model')

import matplotlib.pyplot as plt
plt.close('all')
cc.pl.autok_stability(autok)
plt.grid(True)
plt.show()
#plt.savefig("24_ne_analysis_with_cellcharter_outs/autok_model_stability.pdf", dpi=96)

niches = autok.predict(adata, use_rep='X_cellcharter', k=6).astype("category")
adata = ad.read_h5ad("sgroi-tnbc_filtered.h5ad")
adata.obs["cellcharter_niche"] = niches
adata.write_h5ad("sgroi-tnbc_filtered.h5ad")
adata

adata.obsm["spatial"] = np.array(adata.obs[["x_global", "y_global"]])
sq.pl.spatial_scatter(
    adata, 
    shape = None,
    color=['cellcharter_niche'], 
    library_key='core_global', 
    img=None, 
    title=["tma2_E1", "tma4_F3"],
    size=0.25,
    #connectivity_key='spatial_connectivities',
    #edges_width=0.2,
    legend_loc=None,
    library_id=["tma2_E1", "tma4_F3"]
    )
plt.show()

