rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("sgroi-tnbc_filtered.h5ad", mode = "r+")
adata$obs |> dplyr::glimpse()

# Add globally-defined spatial coordinates
adata$obs$x_global <- dplyr::case_when(adata$obs$slide == "tma2" ~ adata$obs$x_centroid + 15000, 
                                       adata$obs$slide == "tma4" ~ adata$obs$x_centroid + 15000, 
                                       T ~ adata$obs$x_centroid)
adata$obs$y_global <- dplyr::case_when(adata$obs$slide == "tma3" ~ adata$obs$y_centroid + 20000, 
                                       adata$obs$slide == "tma4" ~ adata$obs$y_centroid + 20000, 
                                       T ~ adata$obs$y_centroid)
plot(adata$obs$x_global, adata$obs$y_global, pch = ".", asp = 1)

# Note that for use with squidpy and many python-based tools, you will need to 
# transfer these coordinates into obsm/spatial

# Clean up metadata a bit
adata$obs <- dplyr::relocate(adata$obs, core_global, .after = core)
adata$obs <- dplyr::relocate(adata$obs, x_global, .after = y_centroid)
adata$obs <- dplyr::relocate(adata$obs, y_global, .after = x_global)
adata$obs |> dplyr::glimpse()

# Add HieraType results and ensure that cell-typing versions make sense
hres <- readRDS(file = "18_hieratype_outs/hieratype_custom_lymphoid_pipeline_results.RDS")
adata$obs$celltype_hieratype <- NA
adata$obs[hres$post_probs$lymphoidmajor$cell_ID,]$celltype_hieratype <- hres$post_probs$lymphoidmajor$celltype_granular
adata$obs$celltype_final_amended_1 <- dplyr::case_when(!is.na(adata$obs$celltype_hieratype) ~ adata$obs$celltype_hieratype, 
                                                       T ~ adata$obs$celltype_final)
adata$obs <- dplyr::relocate(adata$obs, celltype_hieratype, celltype_final_amended_1, .after = celltype_final)
adata$obs |> dplyr::glimpse()
adata$obs$celltype_final_amended_1 |> unique()

# Saving again
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered.h5ad", mode = "w")

# Notes: 
# - X and layers/counts are the same, so one could probably delete the latter to save space
# - I used the naming convention _amended_1 because we will probably make more amendments

# Picking up where I left off, after doing Myeloid sub-typing.
rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("sgroi-tnbc_filtered.h5ad", mode = "r+")
adata$obs |> dplyr::glimpse()

# Add HieraType results and ensure that cell-typing versions make sense
hres <- readRDS(file = "18_hieratype_outs/hieratype_custom_myeloid_pipeline_results.RDS")
adata$obs$celltype_hieratype <- as.character(adata$obs$celltype_hieratype)
adata$obs[hres$post_probs$myeloidmajor$cell_ID,]$celltype_hieratype <- hres$post_probs$myeloidmajor$celltype_granular
adata$obs$celltype_final_amended_2 <- dplyr::case_when(adata$obs$celltype_hieratype == "NA" ~ as.character(adata$obs$celltype_final), 
                                                       T ~ adata$obs$celltype_hieratype)
adata$obs <- dplyr::relocate(adata$obs, celltype_final_amended_2, .after = celltype_final_amended_1)
adata$obs |> dplyr::glimpse()
adata$obs$celltype_final_amended_2 <- as.factor(adata$obs$celltype_final_amended_2)
adata$obs$celltype_final_amended_2 |> levels()
BPCells::plot_embedding(source = adata$obs$celltype_final_amended_2, embedding = adata$obsm$X_umap, rasterize = T, labels_discrete = F)

# Saving again
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered.h5ad", mode = "w")

