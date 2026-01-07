rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

library(BPCells)
library(Matrix)

## QC & PP ITER 2 --------------------------------------------------------------
# We are now re-running the pipeline without the filtered cells.

# Reading the data
adata = anndataR::read_h5ad(path = "sgroi-tnbc_filtered.h5ad", mode = "r+")

# UMAP using the scVI dimensions
UM <- uwot::umap(X = adata$obsm$X_scVI, n_neighbors = 50, min_dist = 0.1, metric = "euclidean", nn_method = "hnsw", spread = 1, 
                 fast_sgd = F, n_threads = 8, verbose = T, n_epochs = 500)

# Quick viz
umemb <- UM |> as.data.frame()
colnames(umemb) <- c("umap_1", "umap_2")
umemb$core_global <- adata$obs$core_global
umemb$sample_type <- adata$obs$sample_type

tinyplot::plt(umap_2 ~ umap_1 | core_global, data = umemb, pal = "Polychrome 36", pch = ".", legend = legend(pt.cex = 9))

# Adding to the anndata
adata$obsm$X_umap <- UM

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

# Checking out the results with the normalized data
# dir.create("07_qc_and_pp_iter2_outs")
# pdf(file = "07_qc_and_pp_iter2_outs/umap_by_markers.pdf", width = 10, height = 6)
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("EPCAM",
               "FN1",
               "CD68",
               "CIITA",
               "SPARCL1", 
               "VCAN"
               ),
  rasterize = T, 
  colors_continuous = viridis::viridis(n = 71)
)
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("CD4", 
               "CD8A", 
               "CD3E", 
               "CD19", 
               "CD2", 
               "NKG7"
               ), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
)
# dev.off()

# Leiden clustering (for many resolutions this time)
snn <- knn_hnsw(adata$obsm$X_scVI, k = 50, metric = "euclidean", ef = 1500) |> # Find approximate nearest neighbors
  knn_to_snn_graph() # Convert to a SNN graph
iterative_clustering <- list()
for (i in seq(0.2, 0.8, 0.05)) {
  iterative_clustering[[as.character(i)]] <- cluster_graph_leiden(snn, resolution = i)
}
iterative_clustering <- dplyr::bind_cols(iterative_clustering)
colnames(iterative_clustering) <- paste("leiden_res", colnames(iterative_clustering), sep = "_")

# pdf(file = "07_qc_and_pp_iter2_outs/umap_by_clusters.pdf", width = 6, height = 6)
for (nm in colnames(iterative_clustering)) {
  p <- plot_embedding(iterative_clustering[[nm]], adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
    ggplot2::labs(title = nm)
  print(p)
}
# dev.off()

# Adding to the AnnData and re-saving
adata$obs <- cbind(adata$obs, iterative_clustering)
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered.h5ad", mode = "w")

