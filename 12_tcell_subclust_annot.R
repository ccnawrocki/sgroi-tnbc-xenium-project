rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

library(BPCells)
library(Matrix)

# Reading the data
adata = anndataR::read_h5ad(path = "sgroi-tnbc_filtered_tcell_subset.h5ad", mode = "r+")

# UMAP using the scVI dimensions
UM <- uwot::umap(X = adata$obsm$X_scVI, n_neighbors = 30, min_dist = 0.25, metric = "euclidean", nn_method = "hnsw", spread = 2, 
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
snn <- knn_hnsw(adata$obsm$X_scVI, k = 30, metric = "euclidean", ef = 1500) |> # Find approximate nearest neighbors
  knn_to_snn_graph() # Convert to a SNN graph
iterative_clustering <- list()
for (i in seq(0.05, 1, 0.05)) {
  iterative_clustering[[as.character(i)]] <- cluster_graph_leiden(snn, resolution = i)
}
iterative_clustering <- dplyr::bind_cols(iterative_clustering)
colnames(iterative_clustering) <- paste("sub_sub_leiden_res", colnames(iterative_clustering), sep = "_")

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

# Canonical markers
library(ggplot2)
# pdf("12_tcell_subclust_annot_outs/umap_by_rna_markers_tcells.pdf", width = 12, height = 12)
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("CD4", 
               "CD8A",
               "CD8B",
               "TNF", 
               "PRF1", 
               "FOXP3", 
               "CTLA4", 
               "GNLY", 
               "NKG7", 
               "GZMB", 
               "GZMA", 
               "IFNG"
  ), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
) & 
  labs(x = "UMAP1", y = "UMAP2")
# dev.off()

# Annotating
# 0.35 looks reasonable, so we will continue with that.
# pdf("12_tcell_subclust_annot_outs/umap_by_leiden_res_0.35_tcells.pdf", width = 6, height = 6)
plot_embedding(source = adata$obs$sub_sub_leiden_res_0.35, embedding = adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
  labs(x = "UMAP1", y = "UMAP2", title = "Resolution = 0.35")
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

# Saving
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered_tcell_subset.h5ad", mode = "w")

# Need to find cluster markers: 
rm(list = ls())
.rs.restartR(clean = T)

library(Matrix)
library(BPCells)
library(ggplot2)

# Getting scVI normalized expression  
SCVI <- reticulate::import("scvi")
ad <- reticulate::import("anndata")

adata <- ad$read_h5ad("sgroi-tnbc_filtered_tcell_subset.h5ad")

model <- SCVI$model$SCVI$load("11_tcell_subclust_scVI_outs/scVI_model4")
denoised <- model$get_normalized_expression(adata = adata, library_size = 1000)
denoised_scaled <- apply(X = denoised, MARGIN = 2, FUN = scale)
rownames(denoised_scaled) <- rownames(denoised)

# Making profiles for each cluster
# renv::install("NanoString-BioStats/InSituType")
prof <- InSituType::Estep(counts = denoised_scaled, 
                          assay_type = "protein",
                          clust = paste0("c", adata$obs$sub_sub_leiden_res_0.35), 
                          neg = array(data = 0, dim = adata$shape[[1]]) # We do not have negative probe info
)$profiles

# renv::install("pheatmap", prompt = F)
# pdf(file = "12_tcell_subclust_annot_outs/heatmap_z-score_scaled_profiles.pdf", width = 6, height = 8)
pheatmap::pheatmap(prof, 
                   color = viridis::plasma(n = 101),
                   fontsize_col = 10, show_rownames = F,
                   treeheight_row = 8, treeheight_col = 8, 
                   main = "Mean Z-score Profiles", 
                   scale = "none")
# dev.off()

# limma-trend + quantile normalization (no voom)
mm <- model.matrix(~0+sub_sub_leiden_res_0.35, 
                   data = adata$obs |> 
                     dplyr::mutate(sub_sub_leiden_res_0.35 = as.character(sub_sub_leiden_res_0.35) |> factor(levels = 1:6)))
colnames(mm) <- paste0("c", 1:6)
lfit <- limma::lmFit(limma::normalizeQuantiles(t(log2(denoised))), mm)

contrast_list <- list()
for (clst in levels(adata$obs$sub_sub_leiden_res_0.35)) {
  contrast_list[[clst]] <- mm[adata$obs$sub_sub_leiden_res_0.35 == clst,] |> colMeans()
}

GetClusterMarks <- function(clst) {
  clst_oi <- contrast_list[[clst]]
  other_clsts <- Reduce(f = "+", x = contrast_list[names(contrast_list) != clst])/(length(contrast_list[names(contrast_list) != clst]))
  clst_vs_other <- clst_oi-other_clsts
  tmp <- limma::contrasts.fit(lfit, clst_vs_other)
  res <- limma::eBayes(tmp, trend = T, robust = T) |> limma::topTable(number = Inf)
  colnames(res) <- c("log2FC", "AveExpr", "t", "p", "p_adj", "B")
  res$gene <- rownames(res)
  rownames(res) <- 1:nrow(res)
  res$cluster <- clst
  return(res)
}

library(magrittr)
cluster_marks_limma <- lapply(as.list(unique(adata$obs$sub_sub_leiden_res_0.35)), GetClusterMarks) |> dplyr::bind_rows()
cluster_marks_limma %<>% dplyr::group_by(cluster) %<>% dplyr::arrange(cluster, desc(log2FC), p_adj)

top_marks <- dplyr::group_by(cluster_marks_limma, cluster) |> 
  dplyr::mutate(rank = order(log2FC, decreasing = T)) |> 
  dplyr::group_by(cluster) |> 
  dplyr::top_n(n = -10, wt = rank)
top_marks <- tidyr::pivot_wider(data = top_marks, id_cols = rank, names_from = cluster, values_from = gene) |> 
  dplyr::arrange(rank)

# write.csv(x = top_marks, file = "12_tcell_subclust_annot_outs/cluster_markers_by_limma-trend.csv", row.names = F)

# pdf(file = "12_tcell_subclust_annot_outs/cluster_markers_by_limma-trend_overall_MA_plot.pdf", width = 6, height = 4)
plot(x = cluster_marks_limma$AveExpr, y = cluster_marks_limma$log2FC, pch = 16, cex = 0.5, xlab = "AveExpr", ylab = "log2FC", main = "MA Plot")
abline(h = 0, lwd = 4, col = "red")
# dev.off()

# pdf(file = "12_tcell_subclust_annot_outs/cluster_markers_by_limma-trend_heatmap.pdf", width = 6, height = 10)
pheatmap::pheatmap(prof[unique(top_marks[,-1] |> as.matrix() |> as.vector()), ], 
                   color = viridis::plasma(n = 101), 
                   breaks = seq(-1, 1, length.out = 102),
                   fontsize_col = 10, fontsize_row = 10, 
                   cellheight = 8, cellwidth = 10, 
                   treeheight_row = 8, treeheight_col = 8, 
                   main = "Markers by limma-trend"
)
# dev.off()

# pdf(file = "12_tcell_subclust_annot_outs/canonical_markers_heatmap.pdf", width = 6, height = 10)
marks <- c(
  "CD2", "CD3D", "CD3E", "CD247", 
  "CD4", "CD8A", "CD8B",
  "CD69", "CD38", 
  "PDCD1", "HAVCR2", "LAG3", "TIGIT", "ENTPD1", 
  "CD160", 
  "GZMA", "GZMB", "GZMH", "GZMK",
  "PRF1", "NKG7", "GNLY", "IFNG", "TNF", 
  "IL2RA", "FOXP3", "CTLA4", 
  "CCR7", "SELL", "IL7R", "CD27", 
  "TBX21", "GATA3", "RORC", "EOMES", "TCF7"
)
pheatmap::pheatmap(prof[marks, ], 
                   color = viridis::plasma(n = 101), 
                   breaks = seq(-1, 1, length.out = 102),
                   fontsize_col = 10, fontsize_row = 10, 
                   cellheight = 15, cellwidth = 15, 
                   treeheight_row = 8, treeheight_col = 8, 
                   main = "Canonical Markers"
)
# dev.off()

adata = anndataR::read_h5ad(path = "sgroi-tnbc_filtered_tcell_subset.h5ad", mode = "r+")
library(data.table)
psb <- presto::collapse_counts(counts_mat = Matrix::t(adata$layers$counts) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                               meta_data = adata$obs,
                               get_norm = F,
                               varnames = c("sub_sub_leiden_res_0.35", "patient"), 
                               min_cells_per_group = 20)
psb$meta_data$logUMI <- log(psb$counts_mat |> colSums())

library(presto)
library(lme4)
library(purrr)
library(dplyr)
# presto_res <- presto.presto(
#   y ~ 1 + (1|sub_sub_leiden_res_0.35) + (1|sub_sub_leiden_res_0.35:patient) + (1|patient) + offset(logUMI),
#   psb$meta_data,
#   psb$counts_mat,
#   size_varname = "logUMI",
#   effects_cov = c("sub_sub_leiden_res_0.35"),
#   ncore = 8,
#   min_sigma = 0.05,
#   family = "poisson",
#   nsim = 1000
# )
# saveRDS(object = presto_res, file = "12_tcell_subclust_annot_outs/presto_model_psb.RDS")

presto_res <- readRDS(file = "12_tcell_subclust_annot_outs/presto_model_psb.RDS")

contrasts_mat <- make_contrast.presto(presto_res, "sub_sub_leiden_res_0.35")
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
                           color = "red", size = 3, max.overlaps = 30, box.padding = 0.25) + 
  theme_classic() + 
  labs(title = "Subcluster Markers") +
  facet_wrap(.~cluster, ncol = 3, scales = "free")
p
# ggsave(filename = "12_tcell_subclust_annot_outs/sub_sub_leiden_res_0.35_cluster_markers_volcano_plots.pdf", height = 12, width = 10)

MARKS <- effects_marginal |> filter(fdr < 0.01 & logFC > 1.5) |> arrange(cluster, desc(zscore)) |> pull(feature) |> unique()
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names), 
                  features = MARKS, 
                  groups = adata$obs$sub_sub_leiden_res_0.35) + 
  scale_color_viridis_c(oob = scales::squish, limits = c(-1, 1), option = "B", direction = -1) + 
  labs(y = "Cluster", title = "DE Markers") + 
  theme(axis.title.x = element_blank())
# ggsave(filename = "12_tcell_subclust_annot_outs/bubble_plot_for_DE_markers.pdf", height = 4, width = 12)
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names), 
                  features = marks, 
                  groups = adata$obs$sub_sub_leiden_res_0.35) + 
  scale_color_viridis_c(oob = scales::squish, limits = c(-1, 1), option = "B", direction = -1) + 
  labs(y = "Cluster", title = "Canonical Markers") + 
  theme(axis.title.x = element_blank())
# ggsave(filename = "12_tcell_subclust_annot_outs/bubble_plot_for_canonical_markers.pdf", height = 4, width = 12)

# Determining the cluster identities: 
celltype_l3_map <- c(
  "1" = "Treg",
  "2" = "T.CD4.naive",
  "3" = "NK",
  "4" = "T.CD8.activated.1",
  "5" = "T.exhausted",
  "6" = "T.CD8.activated.2"
)
adata$obs$celltype_level3 <- plyr::mapvalues(x = adata$obs$sub_sub_leiden_res_0.35, from = names(celltype_l3_map), to = celltype_l3_map)
# pdf("12_tcell_subclust_annot_outs/umap_by_celltypes_level3_t_cells.pdf", width = 6, height = 6)
plot_embedding(adata$obs$celltype_level3, adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
  labs(title = "Cell Types (level 3)", x = "UMAP1", y = "UMAP2")
# dev.off()

# Saving 
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered_tcell_subset.h5ad", mode = "w")

# Finally, assign the lower-quality T cells to one of these clusters.
### TBD ###





