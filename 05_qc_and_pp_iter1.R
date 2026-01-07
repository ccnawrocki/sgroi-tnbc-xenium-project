rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

# renv::install(packages = c("tinyplot", "viridis"), prompt = F)

library(BPCells)
library(Matrix)

## QC & PP ITER 1 --------------------------------------------------------------
# We want to identify the noise clusters and filter them out.
# Next, we will re-run the pipeline without these cells.

# Reading the data
adata = anndataR::read_h5ad(path = "sgroi-tnbc.h5ad", mode = "r+")

# UMAP using the scVI dimensions
UM <- uwot::umap(X = adata$obsm$X_scVI, n_neighbors = 50, min_dist = 0.1, metric = "euclidean", nn_method = "hnsw", spread = 1, 
                 fast_sgd = F, n_threads = 8, verbose = T, n_epochs = 500)

# Quick viz
umemb <- UM |> as.data.frame()
colnames(umemb) <- c("umap_1", "umap_2")
umemb$core_global <- adata$obs$core_global

tinyplot::plt(umap_2 ~ umap_1 | core_global, data = umemb, pal = "Polychrome 36", pch = ".", legend = legend(pt.cex = 9), asp = 1)

# Adding to the anndata
adata$obsm$X_umap <- UM

# Normalization
scaling_factor <- 1000
norm_factors <- Matrix::Diagonal(x = scaling_factor/adata$obs$transcript_counts, names=rownames(adata$layers$counts))
norm <- ((norm_factors %*% adata$layers$counts) |> log1p())/log(2)

# Adding to the anndata
adata$layers$lognorm <- norm

# Saving space
remove(norm)

# Checking out the results with the normalized data
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("EPCAM",
               "FN1",
               "CD68", 
               "SPARCL1"),
  rasterize = T, 
  colors_continuous = viridis::viridis(n = 71)
)
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("CD4", 
               "CD8A", 
               "TCF7", 
               "CD3E"), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
)

plot(adata$obsm$X_umap, pch = ".",
     col = viridis::viridis_pal(option = "C")(101)[
       pmin(101, 1 + round(100 * adata$obs$CD3E / quantile(adata$obs$CD3E, 0.99)))
     ]
     )
plot(adata$obsm$X_umap, pch = ".",
     col = viridis::viridis_pal(option = "C")(101)[
       pmin(101, 1 + round(100 * adata$obs$CD20 / quantile(adata$obs$CD20, 0.99)))
     ]
)
plot(adata$obsm$X_umap, pch = ".",
     col = viridis::viridis_pal(option = "C")(101)[
       pmin(101, 1 + round(100 * adata$obs$CD45 / quantile(adata$obs$CD45, 0.99)))
     ]
)
plot(adata$obsm$X_umap, pch = ".",
     col = viridis::viridis_pal(option = "C")(101)[
       pmin(101, 1 + round(100 * adata$obs$CD8A / quantile(adata$obs$CD8A, 0.99)))
     ]
)
plot(adata$obsm$X_umap, pch = ".",
     col = viridis::viridis_pal(option = "C")(101)[
       pmin(101, 1 + round(100 * adata$obs$`E-Cadherin` / quantile(adata$obs$`E-Cadherin`, 0.99)))
     ]
)

# Leiden clustering
clusts <- knn_hnsw(adata$obsm$X_scVI, k = 50, metric = "euclidean", ef = 1500) |> # Find approximate nearest neighbors
  knn_to_snn_graph() |> # Convert to a SNN graph
  cluster_graph_leiden(resolution = 0.5) # Perform graph-based clustering

# dir.create("05_qc_and_pp_iter1_outs")
# pdf(file = "05_qc_and_pp_iter1_outs/umap_by_clusters.pdf", width = 6, height = 6)
plot_embedding(clusts, adata$obsm$X_umap, rasterize = T)
# dev.off()

adata$obs$initial_cluster <- clusts

umemb <- adata$obsm$X_umap |> as.data.frame()
colnames(umemb) <- c("umap_1", "umap_2")
umemb$transcript_counts <- adata$obs$transcript_counts

tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(slide == "tma4"), 
              col = discrete_palette("stallion"),
              pch = ".", legend = legend(pt.cex = 9), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(core_global == "tma4_F3"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.25,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(core_global == "tma1_D1"), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.25,
              legend = legend(pt.cex = 2), asp = 1)

# QC
# renv::install("bioc::scuttle", prompt = F)
plot_embedding(adata$obs$initial_cluster, adata$obsm$X_umap, rasterize = T)

# In most scRNA-seq pipelines nowadays, QC thresholds are set for each sample
# before they are integrated together. Let's apply this principle here for each
# core.

# Using scuttle to get the metrics
library(scuttle)
cts <- adata$layers$counts |> Matrix::t()
dimnames(cts) <- list(adata$var_names, adata$obs_names)
qc_metrics <- scuttle::perCellQCMetrics(x = cts)
summary(qc_metrics$sum)
summary(qc_metrics$detected)

# Finding core-wise (batch-wise) outliers for counts, which would be affected by
# the batch.
all_outliers <- isOutlier(qc_metrics$sum, type = "both", log = T, batch = adata$obs$core_global, 3)
summary(all_outliers)
(attributes(all_outliers)$thresholds["lower",]) |> sort() |> knitr::kable()

# Seems like tma1_A2 and tma4_A3 may need to go... 
tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(core_global == "tma1_A2"), 
              pal = "Polychrome 36",
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | initial_cluster, 
              data = adata$obs |> dplyr::filter(core_global == "tma4_A3"), 
              pal = "Polychrome 36",
              pch = 16, cex = 0.5,
              legend = legend(pt.cex = 2), asp = 1)

# Adding to the AnnData
adata$obs$counts_outlier_scuttle <- as.vector(all_outliers)
tinyplot::plt(y_centroid ~ x_centroid | counts_outlier_scuttle, 
              data = adata$obs |> dplyr::filter(core_global == "tma1_A2"), 
              pch = 16, cex = 0.5, legend = legend(pt.cex = 2), asp = 1)

# Looks reasonable
par(las = 2)
outliers_scuttle <- table(adata$obs$core_global, adata$obs$counts_outlier_scuttle) 
sweep(x = outliers_scuttle, STATS = rowSums(outliers_scuttle), MARGIN = 1, FUN = "/") |> plot(main = "counts outliers")

# Area filtering
areas <- adata$obs$cell_area
par(mar=c(4, 4, 4, 4))
hist(x = log2(areas), breaks = 50)
ol <- isOutlier(log2(areas), log = F, type = "both", nmads = 3)
th <- attr(ol, "threshold")
abline(v = th, col = "blue")

# Adding to the AnnData
adata$obs$area_outlier_scuttle <- as.vector(ol)

# Examining the embedding again: 
plot_embedding(adata$obs$counts_outlier_scuttle, adata$obsm$X_umap, rasterize = T, labels_discrete = F)
plot_embedding(adata$obs$area_outlier_scuttle, adata$obsm$X_umap, rasterize = T, labels_discrete = F)

# Checking cluster-wise:
props <- adata$obs |> dplyr::group_by(initial_cluster) |> 
  dplyr::summarise(prop_ol = mean(counts_outlier_scuttle)) |>
  dplyr::arrange(desc(prop_ol)) |>  data.frame()
Ns <- adata$obs |> dplyr::group_by(initial_cluster) |> 
  dplyr::tally() |>
  dplyr::arrange(desc(n)) |>  data.frame()

dplyr::inner_join(x = props, y = Ns, by = "initial_cluster")
#    initial_cluster     prop_ol      n
# 1               10 0.997124371   2782 --> filter whole thing
# 2                7 0.788327790   3307 --> filter whole thing
# 3                9 0.473278916  12986 --> subject to normal filtering
# 4                2 0.030418956  86262
# 5                5 0.022576119  54881
# 6                3 0.020638769  47435
# 7                8 0.019378629 127460
# 8                6 0.018980916 164481
# 9               12 0.013921114  31463
# 10               4 0.010541207 291238
# 11              11 0.008788689   2617
# 12               1 0.006017491 256253
# 13              13 0.003144654   5406

# renv::install("modelr", prompt = F)
library(presto)
markers <- presto::wilcoxauc(X = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names), 
                             y = adata$obs$initial_cluster)

# 10
markers |> dplyr::filter(group == "10", logFC > 0)
plot_embedding((adata$obs$transcript_counts == 0), adata$obsm$X_umap, rasterize = T, labels_discrete = F)

# 7
markers |> dplyr::filter(group == "7", logFC > 0)

# 9
markers |> dplyr::filter(group == "9", logFC > 0) 

# I do not see any systematic trends on the area outlier axis
props <- adata$obs |> dplyr::group_by(initial_cluster) |> 
  dplyr::summarise(prop_ol = mean(area_outlier_scuttle)) |>
  dplyr::arrange(desc(prop_ol)) |>  data.frame()
dplyr::inner_join(x = props, y = Ns, by = "initial_cluster") 

# Same for nuclei
adata$obs$nuclei_outlier <- ifelse(adata$obs$nucleus_count > 2, yes = T, no = F)
props <- adata$obs |> dplyr::group_by(initial_cluster) |> 
  dplyr::summarise(prop_ol = mean(nuclei_outlier)) |>
  dplyr::arrange(desc(prop_ol)) |>  data.frame()
dplyr::inner_join(x = props, y = Ns, by = "initial_cluster") 

# Adding the var data
ens <- rhdf5::h5read(file = "tma1.h5ad", name = "/var") |> dplyr::bind_rows()
adata$var[["gene"]] <- NA
adata$var[ens$name, "gene"] <- ens$name
adata$var[["ensembl_id"]] <- NA
adata$var[ens$name, "ensembl_id"] <- ens$id

# Saving one last time
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc.h5ad", mode = "w")

# Final filtering
adata <- anndataR::read_h5ad("sgroi-tnbc.h5ad", mode = "r+")
cellstokeep <- !(adata$obs$counts_outlier_scuttle | # Counts outliers
                   adata$obs$area_outlier_scuttle | # Area outliers
                   adata$obs$nuclei_outlier | # Nuclei outliers
                   (adata$obs$initial_cluster %in% c(10, 7)) # Obvious noise clusters
                 ) | 
  (adata$obs$initial_cluster == "") # Protected clusters (white-listed)
mean(cellstokeep) # 0.9677766

# Final filtering
reticulate::py_require("anndata")
ad <- reticulate::import("anndata")
adata <- ad$read_h5ad("sgroi-tnbc.h5ad")
adata_filtered <- adata[cellstokeep, ]
adata_filtered$write_h5ad("sgroi-tnbc_filtered.h5ad")

# A couple summary outputs
adata_filtered$obs |> 
  group_by(slide, patient) |> 
  reframe(core = unique(core)) |> 
  group_by(slide, patient) |> 
  reframe(cores = stringr::str_flatten(unique(core), collapse = ", ")) |> 
  arrange(slide) # |> 
  # write.csv(file = "05_qc_and_pp_iter1_outs/core_summary.csv")

# pdf(file = "05_qc_and_pp_iter1_outs/QC_summary.pdf", width = 10, height = 6)
par(las = 2, mar = c(1, 5, 2, 5), yaxt = "n")
table(adata2$obs$core_global, ifelse(test = rownames(adata$obs) %in% rownames(adata_filtered$obs), yes = "retained", no = "filtered")) |> 
  plot(main = "QC Summary")
# dev.off()

