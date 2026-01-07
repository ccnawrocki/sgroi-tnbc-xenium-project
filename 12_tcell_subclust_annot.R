rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

library(BPCells)
library(Matrix)

# Reading the data
adata = anndataR::read_h5ad(path = "sgroi-tnbc_filtered_tcell_subset.h5ad", mode = "r+")

# UMAP using the scVI dimensions
UM <- uwot::umap(X = adata$obsm$X_scVI, n_neighbors = 50, min_dist = 0.1, metric = "euclidean", nn_method = "hnsw", spread = 1, 
                 fast_sgd = F, n_threads = 8, verbose = T, n_epochs = 500)

# Adding to the anndata
adata$obsm$X_umap <- UM

# Viz
plot_embedding(source = adata$obs$patient, embedding = adata$obsm$X_umap, rasterize = T, labels_discrete = F)

# Normalization
cts <- adata$layers$counts |>  as("CsparseMatrix")
scaling_factor <- 1000
norm_factors <- Matrix::Diagonal(x = scaling_factor/adata$obs$transcript_counts, names=rownames(adata$layers$counts))
norm <- ((norm_factors %*% adata$layers$counts) |> log1p())/log(2)

# Adding to the anndata
adata$layers$lognorm <- norm

# Saving space
remove(norm)
remove(cts)

# Clustering
snn <- knn_hnsw(adata$obsm$X_scVI, k = 50, metric = "euclidean", ef = 1500) |> # Find approximate nearest neighbors
  knn_to_snn_graph() # Convert to a SNN graph
iterative_clustering <- list()
for (i in seq(0.05, 1, 0.05)) {
  iterative_clustering[[as.character(i)]] <- cluster_graph_leiden(snn, resolution = i)
}
iterative_clustering <- dplyr::bind_cols(iterative_clustering)
colnames(iterative_clustering) <- paste("sub_leiden_res", colnames(iterative_clustering), sep = "_")

# Plotting
# dir.create("12_tcell_subclust_annot_outs")
# pdf(file = "12_tcell_subclust_annot_outs/umap_by_clusters.pdf", width = 6, height = 6)
for (nm in colnames(iterative_clustering)) {
  p <- plot_embedding(iterative_clustering[[nm]], adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
    ggplot2::labs(title = nm)
  print(p)
}
# dev.off()

# Adding to the AnnData
adata$obs <- cbind(adata$obs, iterative_clustering)

# Annotating
# 0.7 looks reasonable, so we will continue with that.
library(ggplot2)
# pdf("12_tcell_subclust_annot_outs/umap_by_leiden_res_0.7_tcells.pdf", width = 6, height = 6)
plot_embedding(source = adata$obs$sub_leiden_res_0.7, embedding = adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
  labs(x = "UMAP1", y = "UMAP2")
# dev.off()

# Canonical markers
# pdf("12_tcell_subclust_annot_outs/umap_by_rna_markers_tcells.pdf", width = 12, height = 12)
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("CD4", 
               "CD8A",
               "CD8B",
               "CD3E", 
               "CD2", 
               "CD19", 
               "FOXP3", 
               "CTLA4", 
               "GNLY", 
               "NKG7", 
               "GZMB", 
               "GZMA"
  ), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
) & 
  labs(x = "UMAP1", y = "UMAP2")
# dev.off()

# Protein
protein_mat <- adata$obs[,19:26] |> t() |> as.matrix() |> as("CsparseMatrix")
# pdf("12_tcell_subclust_annot_outs/umap_by_protein_markers_tcells.pdf", width = 12, height = 12)
plot_embedding(
  source = protein_mat,
  embedding = adata$obsm$X_umap,
  features = c("CD45", 
               "CD4", 
               "CD3E", 
               "CD8A", 
               "CD20"
  ),
  rasterize = T, 
  colors_continuous = viridis::viridis(n = 71, option = "C")
) & 
  labs(x = "UMAP1", y = "UMAP2")
# dev.off()

# Need to find cluster markers: 
library(data.table)
psb <- presto::collapse_counts(counts_mat = Matrix::t(adata$layers$counts) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                               meta_data = adata$obs,
                               get_norm = F,
                               varnames = c("sub_leiden_res_0.7", "patient"))
psb$meta_data$logUMI <- log(psb$counts_mat |> colSums())

library(presto)
library(lme4)
library(purrr)
library(dplyr)
# presto_res <- presto.presto(
#   y ~ 1 + (1|sub_leiden_res_0.7) + (1|patient) + (1|sub_leiden_res_0.7:patient) + offset(logUMI),
#   psb$meta_data,
#   psb$counts_mat,
#   size_varname = "logUMI",
#   effects_cov = c("sub_leiden_res_0.7"),
#   ncore = 8,
#   min_sigma = 0.05,
#   family = "poisson",
#   nsim = 1000
# )
# saveRDS(object = presto_res, file = "12_tcell_subclust_annot_outs/presto_model_psb.RDS")

presto_res <- readRDS(file = "12_tcell_subclust_annot_outs/presto_model_psb.RDS")

contrasts_mat <- make_contrast.presto(presto_res, "sub_leiden_res_0.7")
effects_marginal <- contrasts.presto(presto_res, contrasts_mat, one_tailed = T) |> 
  dplyr::mutate(cluster = contrast) |> 
  dplyr::mutate(
    logFC = sign(beta) * log2(exp(abs(beta))), # convert stats to log2
    SD = log2(exp(sigma)),
    zscore = logFC / SD
  ) |> 
  dplyr::select(cluster, feature, logFC, SD, zscore, pvalue) |> 
  dplyr::arrange(pvalue)
effects_marginal$fdr <- p.adjust(effects_marginal$pvalue, method = "BH")

library(ggplot2)
p <- effects_marginal |> 
  ggplot() + 
  scattermore::geom_scattermore(mapping = aes(x = logFC, y = -log10(fdr)), pointsize = 3) + 
  geom_hline(yintercept = -log10(0.01)) + 
  geom_vline(xintercept = 1.5) + 
  ggrepel::geom_text_repel(data = effects_marginal |> filter(fdr < 0.01 & logFC > 1.5), 
                           mapping = aes(x = logFC, y = -log10(fdr), label = feature), 
                           color = "red", size = 4, max.overlaps = 30, box.padding = 0.25) + 
  theme_classic() + 
  labs(title = "Subcluster Markers") +
  facet_wrap(.~cluster, ncol = 4, scales = "free")
p
# ggsave(filename = "12_tcell_subclust_annot_outs/leiden_res_0.7_tcells_cluster_markers_volcano_plots.pdf", height = 12, width = 12)

# Determining the cluster identities: 
# Immediately obvious is that clusters 2, 3, 6, 8, and 13 look strange. We will 
# look at these in the explorer, after we annotate the rest.

# 1: Treg 1 -- Maybe a proliferative subtype. Clear expression of FOXP3, IL2RA, and CTLA4. Also expression of MKI67, HIST1H2BC, TUBA1B, CDK1.
# 4: NK -- Clear and unique expression of GNLY. Also positive for NKG7 and granzymes. Not positive for CD8. 
# 5: CD4+ T -- Memory subtype? LTB is a strong marker.  
# 7: CD4+ T -- Exhausted subtype? CXCL13, PDCD1, CD40LG, and IFNG are all markers. 
# 9: Treg 2 -- Clear expression of FOXP3, IL2RA, and CTLA4.
# 10: CD8+ T cells -- Tissue-resident memory? ITGAE is a strong marker. 
# 11: Treg 3 -- Clear expression of FOXP3, IL2RA, and CTLA4.
# 12: CD8+ T -- Not exactly sure the type. 

table(adata$obs$sub_leiden_res_0.7, adata$obs$slide)
adata$obs |> filter(slide == "tma4") |>
  select(cell_id, sub_leiden_res_0.7) |>
  rename(group = sub_leiden_res_0.7) |>
  write.csv(file = "tmp.csv", row.names = F)

# 2: Seems to show endothelial markers.On the explorer, hard to say... but does not look like endothelial cells. Probably contaminated lymphocytes.
# 3: Seems to show myeloid markers. On the explorer, hard to say...
# 6: Seems to show tumor markers. On the explorer, this appears to be TILs. 
# 8: Seems to show stromal markers. I looked at this one on the explorer. It appears to be CD8+ T cells that are contaminated.
# 13: CD8+ T -- Clearly not B cells, but often next to them. This explains the strange markers... it is contamination.

unlink("tmp.csv")

# Labeling level 3 of cell-typing
celltype_l3_map <- c(
  "1" = "c1_Treg_1",
  "2" = "c2_ambiguous",
  "3" = "c3_ambiguous",
  "4" = "c4_NK",
  "5" = "c5_CD4+T_1",
  "6" = "c6_ambiguous",
  "7" = "c7_CD4+T_2",
  "8" = "c8_ambiguous",
  "9" = "c9_Treg_2",
  "10" = "c10_CD8+T_1",
  "11" = "c11_Treg_3",
  "12" =  "c12_CD8+T_2",
  "13" = "c13_ambiguous"
)
adata$obs$celltype_level3 <- plyr::mapvalues(x = adata$obs$sub_leiden_res_0.7, from = names(celltype_l3_map), to = celltype_l3_map)
# pdf("12_tcell_subclust_annot_outs/umap_by_celltypes_level3_t_cells.pdf", width = 6, height = 6)
plot_embedding(adata$obs$celltype_level3, adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
  labs(title = "Cell Types (level 3)", x = "UMAP1", y = "UMAP2")
# dev.off()

# Saving 
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered_tcell_subset.h5ad", mode = "w")



# We will set up a uv venv here
# $ uv venv
# $ source .venv/bin/activate
# $ uv pip install scvi-tools

# Telling renv to use that uv venv
# renv::use_python(python = ".venv/bin/python")

rm(list = ls())
.rs.restartR(clean = T)

library(Matrix)
library(BPCells)
library(ggplot2)


## scVI normalized expression  -------------------------------------------------
SCVI <- reticulate::import("scvi")
ad <- reticulate::import("anndata")
adata <- ad$read_h5ad("sgroi-tnbc_filtered_tcell_subset.h5ad")

model <- SCVI$model$SCVI$load("11_tcell_subclust_scVI_outs/scVI_model4")
denoised <- model$get_normalized_expression(adata = adata, library_size = 1000)


## Plotting some canonical markers ---------------------------------------------
plot_embedding(
  source = t(denoised),
  embedding = adata$obsm["X_umap"],
  features = c(
    "CD4", 
    "CD8A",
    "CD8B",
    # "CD3E", 
    # "CD2", 
    "CD19", 
    "FOXP3", 
    "CTLA4", 
    "GNLY", 
    "NKG7", 
    "GZMB", 
    "GZMA", 
    "IFNG", 
    "TNF"
  ),
  rasterize = T, 
  colors_continuous = viridis::viridis(n = 71), 
  quantile_range = c(0.01, 0.99)
) & 
  labs(x = "UM1", y = "UM2")





