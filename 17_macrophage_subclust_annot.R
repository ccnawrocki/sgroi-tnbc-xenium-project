rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

library(BPCells)
library(Matrix)

# Reading the data
adata = anndataR::read_h5ad(path = "sgroi-tnbc_filtered_macrophage_subset.h5ad", mode = "r+")

# UMAP using the scVI dimensions
UM <- uwot::umap(X = adata$obsm$X_scVI, n_neighbors = 30, min_dist = 0.25, metric = "euclidean", nn_method = "hnsw", spread = 2, 
                 fast_sgd = F, n_threads = 8, verbose = T, n_epochs = 500)

# Adding to the anndata
adata$obsm$X_umap <- UM

# Viz
plot_embedding(source = adata$obs$patient, embedding = adata$obsm$X_umap, rasterize = T, labels_discrete = F)
plot_embedding(source = adata$obs$celltype_level2, embedding = adata$obsm$X_umap, rasterize = T, labels_discrete = F)

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
# dir.create("17_macrophage_subclust_annot_outs")
# pdf(file = "17_macrophage_subclust_annot_outs/umap_by_clusters.pdf", width = 6, height = 6)
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
# pdf("17_macrophage_subclust_annot_outs/umap_by_rna_markers_macrophages.pdf", width = 12, height = 12)
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("IL4", 
               "IL13",
               "CXCL13", 
               "CXCL1", 
               "CXCL2", 
               "CCL18", 
               "CCL20", 
               "CCL26", 
               "IL10", 
               "TGFBI", 
               "TGFB1"
  ), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
) & 
  labs(x = "UMAP1", y = "UMAP2")
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("IFNG", 
               "TNF", 
               "IL1B", 
               "IL6", 
               "CCL8", 
               "CCL19", 
               "CCL20", 
               "CXCL5", 
               "CXCL9", 
               "CXCL10"
  ), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
) & 
  labs(x = "UMAP1", y = "UMAP2")
# dev.off()

# Annotating
# 0.2 looks reasonable, so we will continue with that.
# pdf("17_macrophage_subclust_annot_outs/umap_by_leiden_res_0.2_macrophages.pdf", width = 6, height = 6)
plot_embedding(source = adata$obs$sub_sub_leiden_res_0.2, embedding = adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
  labs(x = "UMAP1", y = "UMAP2", title = "Resolution = 0.2")
# dev.off()

# Saving
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered_macrophage_subset.h5ad", mode = "w")

# Need to find cluster markers: 
rm(list = ls())
.rs.restartR(clean = T)

library(Matrix)
library(BPCells)
library(ggplot2)

# Getting scVI normalized expression  
SCVI <- reticulate::import("scvi")
ad <- reticulate::import("anndata")

adata <- ad$read_h5ad("sgroi-tnbc_filtered_macrophage_subset.h5ad")

model <- SCVI$model$SCVI$load("16_macrophage_subclust_scVI_outs/scVI_model5")
denoised <- model$get_normalized_expression(adata = adata, library_size = 1000)
denoised_scaled <- apply(X = denoised, MARGIN = 2, FUN = scale)
rownames(denoised_scaled) <- rownames(denoised)

# Making profiles for each cluster
prof <- InSituType::Estep(counts = denoised_scaled, 
                          assay_type = "protein",
                          clust = paste0("c", adata$obs$sub_sub_leiden_res_0.2), 
                          neg = array(data = 0, dim = adata$shape[[1]]) # We do not have negative probe info
)$profiles

# pdf(file = "17_macrophage_subclust_annot_outs/heatmap_z-score_scaled_profiles.pdf", width = 6, height = 8)
pheatmap::pheatmap(prof, 
                   color = viridis::plasma(n = 101),
                   fontsize_col = 10, show_rownames = F,
                   treeheight_row = 8, treeheight_col = 8, 
                   main = "Mean Z-score Profiles", 
                   scale = "none")
# dev.off()

# limma-trend + quantile normalization (no voom)
mm <- model.matrix(~0+sub_sub_leiden_res_0.2, 
                   data = adata$obs |> 
                     dplyr::mutate(sub_sub_leiden_res_0.2 = as.character(sub_sub_leiden_res_0.2) |> factor(levels = 1:3)))
colnames(mm) <- paste0("c", 1:3)
lfit <- limma::lmFit(limma::normalizeQuantiles(t(log2(denoised))), mm)

contrast_list <- list()
for (clst in levels(adata$obs$sub_sub_leiden_res_0.2)) {
  contrast_list[[clst]] <- mm[adata$obs$sub_sub_leiden_res_0.2 == clst,] |> colMeans()
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
cluster_marks_limma <- lapply(as.list(unique(adata$obs$sub_sub_leiden_res_0.2)), GetClusterMarks) |> dplyr::bind_rows()
cluster_marks_limma %<>% dplyr::group_by(cluster) %<>% dplyr::arrange(cluster, desc(log2FC), p_adj)

top_marks <- dplyr::group_by(cluster_marks_limma, cluster) |> 
  dplyr::mutate(rank = order(log2FC, decreasing = T)) |> 
  dplyr::group_by(cluster) |> 
  dplyr::top_n(n = -10, wt = rank)
top_marks <- tidyr::pivot_wider(data = top_marks, id_cols = rank, names_from = cluster, values_from = gene) |> 
  dplyr::arrange(rank)

# pdf(file = "17_macrophage_subclust_annot_outs/cluster_markers_by_limma-trend_overall_MA_plot.pdf", width = 6, height = 4)
plot(x = cluster_marks_limma$AveExpr, y = cluster_marks_limma$log2FC, pch = 16, cex = 0.5, xlab = "AveExpr", ylab = "log2FC", main = "MA Plot")
abline(h = 0, lwd = 4, col = "red")
# dev.off()

# pdf(file = "17_macrophage_subclust_annot_outs/cluster_markers_by_limma-trend_heatmap.pdf", width = 6, height = 10)
pheatmap::pheatmap(prof[unique(top_marks[,-1] |> as.matrix() |> as.vector()), ], 
                   color = viridis::plasma(n = 101), 
                   breaks = seq(-1, 1, length.out = 102),
                   fontsize_col = 10, fontsize_row = 10, 
                   cellheight = 20, cellwidth = 20, 
                   treeheight_row = 8, treeheight_col = 8, 
                   main = "Markers by limma-trend"
)
# dev.off()

# pdf(file = "17_macrophage_subclust_annot_outs/canonical_markers_heatmap.pdf", width = 6, height = 10)
marks <- c(
  "IL4", 
  "IL13",
  "CXCL13", 
  "CXCL1", 
  "CXCL2", 
  "CCL18", 
  "CCL20", 
  "CCL26", 
  "IL10", 
  "TGFBI", 
  "TGFB1",
  
  "IFNG", 
  "TNF", 
  "IL1B", 
  "IL6", 
  "CCL8", 
  "CCL19", 
  "CCL20", 
  "CXCL5", 
  "CXCL9", 
  "CXCL10",
  
  "S100A9"
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

adata <- anndataR::read_h5ad(path = "sgroi-tnbc_filtered_macrophage_subset.h5ad", mode = "r+")
library(data.table)
psb <- presto::collapse_counts(counts_mat = Matrix::t(adata$layers$counts) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                               meta_data = adata$obs,
                               get_norm = F,
                               varnames = c("sub_sub_leiden_res_0.2", "patient"), 
                               min_cells_per_group = 20)
psb$meta_data$logUMI <- log(psb$counts_mat |> colSums())

library(presto)
library(lme4)
library(purrr)
library(dplyr)
# presto_res <- presto.presto(
#   y ~ 1 + (1|sub_sub_leiden_res_0.2) + (1|sub_sub_leiden_res_0.2:patient) + (1|patient) + offset(logUMI),
#   psb$meta_data,
#   psb$counts_mat,
#   size_varname = "logUMI",
#   effects_cov = c("sub_sub_leiden_res_0.2"),
#   ncore = 8,
#   min_sigma = 0.05,
#   family = "poisson",
#   nsim = 1000
# )
# saveRDS(object = presto_res, file = "17_macrophage_subclust_annot_outs/presto_model_psb.RDS")

presto_res <- readRDS(file = "17_macrophage_subclust_annot_outs/presto_model_psb.RDS")

contrasts_mat <- make_contrast.presto(presto_res, "sub_sub_leiden_res_0.2")
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
  geom_vline(xintercept = 1) + 
  ggrepel::geom_text_repel(data = effects_marginal |> filter(fdr < 0.01 & logFC > 1), 
                           mapping = aes(x = logFC, y = -log10(fdr), label = feature), 
                           color = "red", size = 2, max.overlaps = 30, box.padding = 0.25) + 
  theme_classic() + 
  labs(title = "Subcluster Markers") +
  facet_wrap(.~cluster, ncol = 3, scales = "free")
p
# ggsave(filename = "17_macrophage_subclust_annot_outs/sub_sub_leiden_res_0.2_cluster_markers_volcano_plots.pdf", height = 6, width = 10)

MARKS <- effects_marginal |> filter(fdr < 0.01 & logFC > 1) |> arrange(cluster, desc(zscore)) |> pull(feature) |> unique()
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names), 
                  features = MARKS, 
                  groups = adata$obs$sub_sub_leiden_res_0.2) + 
  scale_color_viridis_c(oob = scales::squish, limits = c(-1, 1), option = "B", direction = -1) + 
  labs(y = "Cluster", title = "DE Markers") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 4))
# ggsave(filename = "17_macrophage_subclust_annot_outs/bubble_plot_for_DE_markers.pdf", height = 4, width = 16)
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names), 
                  features = marks, 
                  groups = adata$obs$sub_sub_leiden_res_0.2) + 
  scale_color_viridis_c(oob = scales::squish, limits = c(-1, 1), option = "B", direction = -1) + 
  labs(y = "Cluster", title = "Canonical Markers") + 
  theme(axis.title.x = element_blank())
# ggsave(filename = "17_macrophage_subclust_annot_outs/bubble_plot_for_canonical_markers.pdf", height = 4, width = 10)

top_marks <- dplyr::filter(effects_marginal, fdr < 0.01 & logFC > 1) |> 
  dplyr::group_by(cluster) |> 
  dplyr::mutate(rank = order(zscore, decreasing = T)) |> 
  dplyr::group_by(cluster) |> 
  dplyr::top_n(n = -10, wt = rank)
top_marks <- tidyr::pivot_wider(data = top_marks, id_cols = rank, names_from = cluster, values_from = feature) |> 
  dplyr::arrange(rank)
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names), 
                  features = top_marks[,-1] |> unlist() |> unique(), 
                  groups = adata$obs$sub_sub_leiden_res_0.2, 
                  group_order = colnames(top_marks[,-1])) + 
  scale_color_viridis_c(oob = scales::squish, limits = c(-1, 2), option = "B", direction = -1) + 
  labs(y = "Cluster", title = "DE Markers") + 
  theme(axis.title.x = element_blank())
# ggsave(filename = "17_macrophage_subclust_annot_outs/bubble_plot_for_DE_markers_alt_layout.pdf", height = 4, width = 10)

# Determining the cluster identities: 
### TBD



