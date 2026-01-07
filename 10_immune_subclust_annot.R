rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

library(BPCells)
library(Matrix)

# Reading the data
adata = anndataR::read_h5ad(path = "sgroi-tnbc_filtered_immune_subset.h5ad", mode = "r+")

# UMAP using the scVI dimensions
UM <- uwot::umap(X = adata$obsm$X_scVI, n_neighbors = 50, min_dist = 0.1, metric = "euclidean", nn_method = "hnsw", spread = 1, 
                 fast_sgd = F, n_threads = 8, verbose = T, n_epochs = 500)

# Adding to the anndata
adata$obsm$X_umap <- UM

# Viz
# pdf("10_immune_subclust_annot_outs/umap_by_celltypes_level1_immune.pdf", width = 6, height = 6)
plot_embedding(source = adata$obs$celltype_level1, embedding = adata$obsm$X_umap, rasterize = T, labels_discrete = F)
# dev.off()

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
# dir.create("10_immune_subclust_annot_outs")
# pdf(file = "10_immune_subclust_annot_outs/umap_by_clusters.pdf", width = 6, height = 6)
for (nm in colnames(iterative_clustering)) {
  p <- plot_embedding(iterative_clustering[[nm]], adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
    ggplot2::labs(title = nm)
  print(p)
}
# dev.off()

# Adding to the AnnData
adata$obs <- cbind(adata$obs, iterative_clustering)

# Annotating
# 0.45 looks reasonable, so we will continue with that.
# pdf("10_immune_subclust_annot_outs/umap_by_leiden_res_0.45_immune.pdf", width = 6, height = 6)
plot_embedding(source = adata$obs$sub_leiden_res_0.45, embedding = adata$obsm$X_umap, rasterize = T, labels_discrete = F)
# dev.off()

# Canonical markers
# pdf("10_immune_subclust_annot_outs/umap_by_markers_immune.pdf", width = 12, height = 12)
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("CD4", 
               "CD8A", 
               "CD3E", 
               "CD2", 
               "CD19", 
               "LAMP3", 
               "FOXP3", 
               "JCHAIN", 
               "CD163"
  ), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
) & 
  labs(x = "UMAP1", y = "UMAP2")
# dev.off()

# Need to find cluster markers: 
library(data.table)
psb <- presto::collapse_counts(counts_mat = Matrix::t(adata$layers$counts) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                               meta_data = adata$obs,
                               get_norm = F,
                               varnames = c("sub_leiden_res_0.45", "patient"))
psb$meta_data$logUMI <- log(psb$counts_mat |> colSums())

library(presto)
library(lme4)
library(purrr)
library(dplyr)
# presto_res <- presto.presto(
#   y ~ 1 + (1|sub_leiden_res_0.45) + (1|sub_leiden_res_0.45:patient) + (1|patient) + offset(logUMI),
#   psb$meta_data,
#   psb$counts_mat,
#   size_varname = "logUMI",
#   effects_cov = c("sub_leiden_res_0.45"),
#   ncore = 8,
#   min_sigma = 0.05,
#   family = "poisson",
#   nsim = 1000
# )
# saveRDS(object = presto_res, file = "10_immune_subclust_annot_outs/presto_model_psb.RDS")

presto_res <- readRDS(file = "10_immune_subclust_annot_outs/presto_model_psb.RDS")

contrasts_mat <- make_contrast.presto(presto_res, "sub_leiden_res_0.45")
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
                           color = "red", size = 2, max.overlaps = 30, box.padding = 0.25) + 
  theme_classic() + 
  labs(title = "Subcluster Markers") +
  facet_wrap(.~cluster, ncol = 4, scales = "free")
p
# ggsave(filename = "10_immune_subclust_annot_outs/sub_leiden_res_0.45_cluster_markers_volcano_plots.pdf", height = 12, width = 10)

# Determining the cluster identities: 
effects_marginal[effects_marginal$cluster == 1 & effects_marginal$logFC > 1.5, ] |> arrange(fdr) |> knitr::kable() # TAM
effects_marginal[effects_marginal$cluster == 2 & effects_marginal$logFC > 2, ] |> arrange(fdr) # mregDC
effects_marginal[effects_marginal$cluster == 3 & effects_marginal$logFC > 1.5, ] |> arrange(fdr) |> knitr::kable() # TAM
effects_marginal[effects_marginal$cluster == 4 & effects_marginal$logFC > 2, ] |> arrange(fdr) # Plasmablast
effects_marginal[effects_marginal$cluster == 5 & effects_marginal$logFC > 2, ] |> arrange(fdr) |> knitr::kable() # Tc and NK, probably activated
effects_marginal[effects_marginal$cluster == 6 & effects_marginal$logFC > 2, ] |> arrange(fdr) |> knitr::kable() # Treg
effects_marginal[effects_marginal$cluster == 7 & effects_marginal$logFC > 2, ] |> arrange(fdr) |> knitr::kable() # also Tc and NK, probably exhausted
effects_marginal[effects_marginal$cluster == 8 & effects_marginal$logFC > 2, ] |> arrange(fdr) |> knitr::kable() # Memory T
effects_marginal[effects_marginal$cluster == 9 & effects_marginal$logFC > 2, ] |> arrange(fdr) # Plasma
effects_marginal[effects_marginal$cluster == 10 & effects_marginal$logFC > 2, ] |> arrange(fdr) # pDC
effects_marginal[effects_marginal$cluster == 11 & effects_marginal$logFC > 2, ] |> arrange(fdr) # B

# Labeling level 2 of cell-typing
celltype_l2_map <- c(
  "1" = "TAM",
  "2" = "mregDC",
  "3" = "TAM",
  "4" = "Plasmablast",
  "5" = "Activated Tc and NK",
  "6" = "Treg",
  "7" = "Exhausted Tc and NK",
  "8" = "Memory T",
  "9" = "Plasma",
  "10" = "pDC", 
  "11" = "B"
)
adata$obs$celltype_level2 <- plyr::mapvalues(x = adata$obs$sub_leiden_res_0.45, from = names(celltype_l2_map), to = celltype_l2_map)
# pdf("10_immune_subclust_annot_outs/umap_by_celltypes_level2_immune.pdf", width = 6, height = 6)
plot_embedding(adata$obs$celltype_level2, adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
  labs(title = "Cell Types (level 2)", x = "UMAP1", y = "UMAP2")
# dev.off()

# Saving 
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered_immune_subset.h5ad", mode = "w")

