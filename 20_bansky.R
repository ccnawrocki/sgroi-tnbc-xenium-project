rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("sgroi-tnbc_filtered.h5ad", mode = "r+")

adata$obs$x_global <- dplyr::case_when(adata$obs$slide == "tma2" ~ adata$obs$x_centroid + 15000, 
                                       adata$obs$slide == "tma4" ~ adata$obs$x_centroid + 15000, 
                                       T ~ adata$obs$x_centroid)
adata$obs$y_global <- dplyr::case_when(adata$obs$slide == "tma3" ~ adata$obs$y_centroid + 20000, 
                                       adata$obs$slide == "tma4" ~ adata$obs$y_centroid + 20000, 
                                       T ~ adata$obs$y_centroid)
plot(adata$obs$x_global, adata$obs$y_global, pch = ".")

# Make "spe"
library(SpatialExperiment)
spe <- SpatialExperiment(assays = list(lognorm = adata$layers$lognorm |> # To start, we will use the log-normalized data
                                         magrittr::set_colnames(adata$var_names) |> 
                                         magrittr::set_rownames(adata$obs_names) |> 
                                         Matrix::t() |> as("CsparseMatrix")),
                         colData = adata$obs, 
                         spatialCoords = adata$obs |> dplyr::select(x_global, y_global),
                         spatialCoordsNames = c("x_global", "y_global"), 
                         sample_id = adata$obs$core_global |> as.character()
)
rm(adata)

# We want to do: banksy --> pca --> harmony --> clustering + umap

#renv::install("bioc::Banksy")
library(Banksy)
lambda <- 0.8 # Recommended for tissue niches
k_geom <- c(15, 30) # Recommended by the package authors

# Computes the Banksy data
spe <- Banksy::computeBanksy(spe, assay_name = "lognorm", compute_agf = T, k_geom = k_geom)

# Now we have the original data (lognorm), 
# the (15-NN-derived) median neighborhood expression (H0), 
# and the (30-NN-derived) magnitude of the expression gradient (H1).

# The next step is to scale each matrix (best to do this for each sample separately), 
# then transform them with the lambda value, then run PCA on them. 

# However, this is memory-intense for this dataset. So I will switch to the HPC 
# for this. 

################################################################################

# adata <- anndataR::read_h5ad("sgroi.h5ad", mode = "r+")
# adata$obs$x_global <- dplyr::case_when(adata$obs$slide == "tma2" ~ adata$obs$x_centroid + 15000, 
#                                        adata$obs$slide == "tma4" ~ adata$obs$x_centroid + 15000, 
#                                        T ~ adata$obs$x_centroid)
# adata$obs$y_global <- dplyr::case_when(adata$obs$slide == "tma3" ~ adata$obs$y_centroid + 20000, 
#                                        adata$obs$slide == "tma4" ~ adata$obs$y_centroid + 20000, 
#                                        T ~ adata$obs$y_centroid)
# plot(adata$obs$x_global, adata$obs$y_global, pch = ".")
# 
# # Make "spe"
# library(SpatialExperiment)
# spe <- SpatialExperiment(assays = list(lognorm = adata$layers$lognorm |> # To start, we will use the log-normalized data
#                                          magrittr::set_colnames(adata$var_names) |> 
#                                          magrittr::set_rownames(adata$obs_names) |> 
#                                          Matrix::t() |> as("CsparseMatrix")),
#                          colData = adata$obs, 
#                          spatialCoords = adata$obs |> dplyr::select(x_global, y_global),
#                          spatialCoordsNames = c("x_global", "y_global"), 
#                          sample_id = adata$obs$core_global |> as.character()
# )
# rm(adata)
# 
# # We want to do: banksy --> pca --> harmony --> clustering + umap
# 
# #renv::install("bioc::Banksy")
# library(Banksy)
# lambda <- 0.8 # Recommended for tissue niches
# k_geom <- c(15, 30) # Recommended by the package authors
# 
# # Computes the Banksy data
# spe <- Banksy::computeBanksy(spe, assay_name = "lognorm", compute_agf = T, k_geom = k_geom)
# 
# # The next step is to scale each matrix (best to do this for each sample separately), 
# # then transform them with the lambda value, then run PCA on them. 
# ## -- This is well-illustrated with Figure 1a
# ## -- This can all be achieved with: runBanksyPCA()
# 
# # However, this is memory-intense for this dataset. So we could simply get the 
# # matrix, then use BPCells to do the PCA and clustering. 
# 
# # # Get the Banksy matrix on which to run PCA
# # bmat <- getBanksyMatrix(se = spe, M = 1, lambda = lambda, assay_name = "lognorm", scale = T, group = "sample_id")
# 
# # # We will just run it as usual though... since we are on the HPC with plenty of memory
# # #renv::install("bioc::harmony")
# # # Performing PCA on the Banksy matrix
# # spe <- runBanksyPCA(spe, npcs = 50, M = 1, lambda = lambda, group = "sample_id", assay_name = "lognorm")
# # # spe@int_colData@listData[["reducedDims"]]@listData[["PCA_M1_lam0.8"]] |> head()
# # # saveRDS(object = spe@int_colData@listData[["reducedDims"]]@listData[["PCA_M1_lam0.8"]], file = "BanksyPCA.RDS")
# 
# PCA <- readRDS("BanksyPCA.RDS")
# attr(PCA, which = "percentVar") |> plot()
# 
# # H <- harmony::RunHarmony(data_mat = PCA[,1:15], meta_data = colData(spe), vars_use = "sample_id", 
# #                          plot_convergence = T, verbose = T)
# # saveRDS(H, "BanksyH.RDS")

################################################################################

# Now, I will do clustering, after switching back. This should not be memory intensive
# but it will probably take awhile. 

H <- readRDS(file = "BanksyH.RDS")

# Leiden clustering
library(BPCells)
snn_ngr <- BPCells::knn_hnsw(data = H, k = 50, metric = "euclidean", ef = 1500) |> 
  BPCells::knn_to_snn_graph() 
clusts <- cluster_graph_leiden(mat = snn_ngr, resolution = 0.35)
table(clusts) # 4 singeltons, but that is okay

# UMAP
UM <- uwot::umap(X = H, n_neighbors = 50, min_dist = 0.1, metric = "euclidean", nn_method = "hnsw", spread = 1, 
                 fast_sgd = F, n_threads = 8, verbose = T, n_epochs = 500)

# Viz
d <- colData(spe)
d$UMAP_1 <- UM[,1]
d$UMAP_2 <- UM[,2]
d$niche <- dplyr::case_when(clusts %in% 1:11 ~ clusts, 
                            T ~ "singleton") |> as.factor() # Deals with those 4 singletons
d <- as.data.frame(d)

tinyplot::plt(UMAP_2 ~ UMAP_1 | niche, data = d, pch = ".", col = BPCells::discrete_palette(name = "stallion")[1:12], 
              legend = list(pch = 16, pt.cex = 2))
tinyplot::plt(y_centroid ~ x_centroid | niche, data = d[d$core_global == "tma4_F3",], pch = 16, cex = 0.25, col = BPCells::discrete_palette(name = "stallion")[1:12], 
              legend = list(pch = 16, pt.cex = 2), asp = 1, main = "tma4, coreF3")
tinyplot::plt(y_centroid ~ x_centroid | celltype_level1, data = d[d$core_global == "tma4_F3",], pch = 16, cex = 0.25, col = BPCells::discrete_palette(name = "stallion")[1:12], 
              legend = list(pch = 16, pt.cex = 2), asp = 1)

# Adding to the anndata
adata <- anndataR::read_h5ad("sgroi-tnbc_filtered.h5ad", mode = "r+")
adata$obs$banksy_niche <- plyr::mapvalues(x = adata$obs_names, from = rownames(d), to = d$niche |> as.character()) |> as.factor()

tinyplot::plt(y_centroid ~ x_centroid | banksy_niche, data = adata$obs[adata$obs$core_global == "tma4_F3",], pch = 16, cex = 0.25, col = BPCells::discrete_palette(name = "stallion")[1:12], 
              legend = list(pch = 16, pt.cex = 2), asp = 1, main = "tma4, coreF3")

# Saving
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered.h5ad", mode = "w")


#### Downstream Analysis ------------------------------------------------------
rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("sgroi-tnbc_filtered.h5ad", mode = "r+")

# First, are any of the niches enriched in the HOXB13+ group?
niche_props <- table(adata$obs$banksy_niche, adata$obs$patient)[1:11,] |> # exclude the 4 singletons
  prop.table(margin = 2) |> as.data.frame.matrix() # turn into proportions
groups <- dplyr::group_by(adata$obs, patient) |> dplyr::summarise(g = unique(HOXB13_group)) # get the grouping variable
rownames(groups) <- groups$patient
compstats <- presto::wilcoxauc(X = niche_props, y = groups[colnames(niche_props),]$g) |> # Use presto to do the tests
  dplyr::filter(group == "HOXB13+") |> 
  dplyr::mutate(padj = p.adjust(p = pval, method = "BH"))
compstats <- dplyr::select(compstats, feature, group, logFC, statistic, pval, padj)
colnames(compstats) <- c("niche", "group", "diff", "statistic", "pval", "fdr")
compstats$group <- ifelse(test = compstats$diff > 0, yes = "HOXB13+", no = "HOXB13-")
compstats <- dplyr::arrange(compstats, pval)
compstats

# Is the above code doing what we think? Let's recreate with base R:
outs <- apply(X = niche_props, MARGIN = 1, FUN = function(z) {
  dtmp <- data.frame(x = z, y = groups[colnames(niche_props),]$g)
  testout <- wilcox.test(dtmp$x ~ dtmp$y, correct = F) 
  data.frame(stat = testout$statistic, pval = testout$p.value)
  }
)
outs <- dplyr::bind_rows(outs, .id = "niche")
outs$diff <- apply(X = niche_props, MARGIN = 1, FUN = function(z) {
  dtmp <- data.frame(x = z, y = groups[colnames(niche_props),]$g)
  means <- tapply(dtmp$x, dtmp$y, mean)
  means[1]-means[2]
}, simplify = T
)

# Yep--it leads to the same results. P-values are slightly different because 
# presto uses an approximation. Nonetheless, the interpretations are the same.

#dir.create(path = "20_banksy_outs")
#system(command = "mv BanksyH.RDS 20_banksy_outs/BanksyH.RDS")
#openxlsx::write.xlsx(compstats, "20_banksy_outs/banksy_niches_abundance_differences.xlsx")

# The organization of the niches looks like it varies... maybe we can test this 
# with cellcharter or something like that? 
par(mar = c(1, 1, 1, 1))
d <- adata$obs
d <- d[d$banksy_niche != "singleton",]
d$banksy_niche <- factor(x = d$banksy_niche, levels = 1:11)
tinyplot::plt(y_centroid ~ x_centroid | banksy_niche, data = d[d$core_global == "tma4_D2",], 
              palette = "Alphabet", pch = 16, cex = 0.5, asp = 1, legend = list(pt.cex = 2), 
              frame = F, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, main = "tma4_D2")
tinyplot::plt(y_centroid ~ x_centroid | banksy_niche, data = d[d$core_global == "tma2_D3",], 
              palette = "Alphabet", pch = 16, cex = 0.5, asp = 1, legend = list(pt.cex = 2), 
              frame = F, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, main = "tma2_D3")

