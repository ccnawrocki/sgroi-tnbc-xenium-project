rm(list = ls())
.rs.restartR(clean = T)
.libPaths()

library(BPCells)
library(Matrix)

adata <- anndataR::read_h5ad("sgroi-tnbc_filtered.h5ad", mode = "r+")

# After looking at our clustering results, 0.35 seems like a good resolution to 
# begin with.
plot_embedding(adata$obs$leiden_res_0.35, adata$obsm$X_umap, rasterize = T, labels_discrete = F)

## Finding Cluster Markers
library(data.table)
psb <- presto::collapse_counts(counts_mat = Matrix::t(adata$layers$counts) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                               meta_data = adata$obs,
                               get_norm = F,
                               varnames = c("leiden_res_0.35", "patient"))
psb$meta_data$logUMI <- log(psb$counts_mat |> colSums())
 
library(presto)
library(lme4)
library(purrr)
library(dplyr)
# presto_res <- presto.presto(
#   y ~ 1 + (1|leiden_res_0.35) + (1|leiden_res_0.35:patient) + (1|patient) + offset(logUMI),
#   psb$meta_data,
#   psb$counts_mat,
#   size_varname = "logUMI",
#   effects_cov = c("leiden_res_0.35"),
#   ncore = 8,
#   min_sigma = 0.05,
#   family = "poisson",
#   nsim = 1000
# )
# dir.create("08_cell-typing_outs")
# saveRDS(object = presto_res, file = "08_cell-typing_outs/presto_model_psb.RDS")

presto_res <- readRDS(file = "08_cell-typing_outs/presto_model_psb.RDS")

contrasts_mat <- make_contrast.presto(presto_res, "leiden_res_0.35")
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
  labs(title = "Cluster Markers") +
  facet_wrap(.~cluster, ncol = 4, scales = "free")
p
# dir.create("08_cell-typing_outs")
# ggsave(filename = "08_cell-typing_outs/leiden_res_0.35_cluster_markers_volcano_plots.pdf", height = 12, width = 10)

## Plotting Markers
## 1-2 ## 
# These can easily be identified as epithelial.

## 3 ##
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  features = c("ITGAX", "CD163", "CD68", "CSF1R", "TREM2"), groups = adata$obs$leiden_res_0.35, 
                  colors = viridis::plasma(n = 71)
)
# This seems to be a myeloid cluster. It is positive for CD4 and CIITA as well, 
# which makes sense.

## 4 ##
# This can easily be identified as epithelial.

## 5 ##
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  features = c("CD3E", "CD2", "TIGIT", "CTLA4", "GZMA"), groups = adata$obs$leiden_res_0.35, 
                  colors = viridis::plasma(n = 71)
)
# Seems to be a lymphoid cluster. Probably T and NK cells.

## 6 ##
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  features = c("FLT1", "PLVAP", "SOX17", "RGS5", "ITGA1"), groups = adata$obs$leiden_res_0.35, 
                  colors = viridis::plasma(n = 71)
)
# Endothelial cells.

## 7 ##
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  features = c("COL6A3", "MEG3", "DCN", "LUM", "COL1A1"), groups = adata$obs$leiden_res_0.35, 
                  colors = viridis::plasma(n = 71)
)
# Stromal cells. Note that MEG3 expression is high in many of these cells, which
# may be interesting. 

## 8 ## 
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  features = c("CD19", "BANK1", "MS4A1", "CD79A", "PAX5"), groups = adata$obs$leiden_res_0.35, 
                  colors = viridis::plasma(n = 71)
                  )
# B cells.

## 9 ##
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  features = c("FKBP11", "IGHG2", "IGHGP", "JCHAIN", "IGLC3"), groups = adata$obs$leiden_res_0.35, 
                  colors = viridis::plasma(n = 71)
)
# Plasma cells and plasmablasts.

## 10 ##
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  features = c("CPA3", "MS4A2", "KIT", "HPGDS", "GATA2"), groups = adata$obs$leiden_res_0.35, 
                  colors = viridis::plasma(n = 71)
)
# This cluster definitely includes MAST cells, since CPA3 and MS4A2 are strong
# markers. However, KIT and GATA2 may suggest that some of these cells are 
# maybe other basophils. 

plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("CPA3", "MS4A2", "KIT", "HPGDS", "GATA2"), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
)
# For now, we will label as basophils.

## 11 ##
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  features = c("LILRA4", "PLD4", "GZMB", "IRF8", "MPEG1"), groups = adata$obs$leiden_res_0.35, 
                  colors = viridis::plasma(n = 71)
)
# Plasmacytoid dendritic cells.

## 12 ##
BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  features = c("IL12B", "CCR7", "CD80", "TRAF1", "LAMP3"), groups = adata$obs$leiden_res_0.35, 
                  colors = viridis::plasma(n = 71)
)
# mregDCs, evidenced by LAMP3 and CCR7 being markers.

# Labeling level 1 of cell-typing
celltype_l1_map <- c(
  "1" = "Epithelial",
  "2" = "Epithelial",
  "3" = "Myeloid",
  "4" = "Epithelial",
  "5" = "Lymphoid",
  "6" = "Endothelial",
  "7" = "Stromal",
  "8" = "B",
  "9" = "Plasma",
  "10" = "Basophil", 
  "11" = "pDC", 
  "12" = "mregDC"
)
adata$obs$celltype_level1 <- plyr::mapvalues(x = adata$obs$leiden_res_0.35, from = names(celltype_l1_map), to = celltype_l1_map)
plot_embedding(adata$obs$celltype_level1, adata$obsm$X_umap, rasterize = T, labels_discrete = F)

# Validating with protein data
protein_mat <- adata$obs[,19:26] |> t() |> as.matrix() |> as("CsparseMatrix")
plot_embedding(
  source = protein_mat,
  embedding = adata$obsm$X_umap,
  features = c("CD45", "E-Cadherin", "alphaSMA", "CD4", "CD3E", "Vimentin", "CD8A", "CD20"),
  rasterize = T, 
  colors_continuous = viridis::viridis(n = 71)
)
# One part of that basophil cluster looks like it has really strong alphaSMA 
# signal. This does not make much sense, biologically.

# In addition, I looked at the explorer and found that these cells are clearly 
# luminal secretory cells. 
# table(adata$obs$leiden_res_0.8, adata$obs$slide)
# adata$obs |> filter(slide == "tma1") |> 
#   select(cell_id, leiden_res_0.8) |> 
#   rename(group = leiden_res_0.8) |> 
#   write.csv(file = "tmp.csv", row.names = F)
# adata$obs |> filter(slide == "tma3") |> 
#   select(cell_id, leiden_res_0.8) |> 
#   rename(group = leiden_res_0.8) |> 
#   write.csv(file = "tmp.csv", row.names = F)
# unlink("tmp.csv")

# Furthermore, according to the HPA, KIT can be a marker for these cells. These
# cells should be ID4 positive as well. 
plot_embedding(
  source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
  embedding = adata$obsm$X_umap,
  features = c("KIT", "ID4"), 
  rasterize = T,
  colors_continuous = viridis::viridis(n = 71)
)

# Amending the basophil cluster
adata$obs$celltype_level1 <- ifelse(test = (adata$obs$leiden_res_0.8 == 19 & adata$obs$leiden_res_0.35 == 10), yes = "Ductal", no = as.character(adata$obs$celltype_level1))

# Summary plots
# pdf("08_cell-typing_outs/umap_by_celltype_level1.pdf", width = 6, height = 6)
plot_embedding(adata$obs$celltype_level1, adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
  labs(title = "Cell Types (level 1)", x = "UMAP1", y = "UMAP2")
# dev.off()

BPCells::plot_dot(source = Matrix::t(adata$layers$lognorm) |> as("CsparseMatrix") |> magrittr::set_rownames(value = adata$var_names),
                  groups = adata$obs$celltype_level1, 
                  group_order = c("Epithelial", "Ductal", "Stromal", "Endothelial",
                                  "Myeloid", "Basophil", "mregDC", 
                                  "pDC", "Plasma", "B", "Lymphoid"), 
                  features = c("EPCAM", "ERBB2", "ERBB3", "SOX9", "CENPF", "CEACAM1",
                               "KIT", "ID4", 
                               "COL6A3", "MEG3", "DCN", "LUM", "COL1A1",
                               "FLT1", "PLVAP", "SOX17", "RGS5", "ITGA1",
                               "ITGAX", "CD163", "CD68", "CSF1R", "TREM2",
                               "CPA3", "MS4A2", 
                               "HPGDS", "GATA2",
                               "IL12B", "CCR7", "CD80", "TRAF1", "LAMP3",
                               "LILRA4", "PLD4", "GZMB", "IRF8", "MPEG1",
                               "FKBP11", "IGHG2", "IGHGP", "JCHAIN", "IGLC3",
                               "CD19", "BANK1", "MS4A1", "CD79A", "PAX5",
                               "CD3E", "CD2", "TIGIT", "CTLA4", "GZMA"
                               )
                  # , colors = viridis::plasma(n = 71)
) + 
  scale_color_viridis_c(limits = c(-1, 2), oob = scales::squish, option = "B", direction = -1) + 
  labs(title = "RNA Markers for Cell Types (level 1)", y = "Cell Type", x = "Gene")
# ggsave(filename = "08_cell-typing_outs/bubble_plot_for_celltype_level1_markers.pdf", height = 6, width = 10)

# Looking at some cores: 
tinyplot::plt(y_centroid ~ x_centroid | celltype_level1, 
              data = adata$obs |> dplyr::filter(core_global == "tma4_F3") |> dplyr::mutate(y_centroid = -1*y_centroid), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.25,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | celltype_level1, 
              data = adata$obs |> dplyr::filter(core_global == "tma1_E1") |> dplyr::mutate(y_centroid = -1*y_centroid), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.25,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | celltype_level1, 
              data = adata$obs |> dplyr::filter(core_global == "tma2_E1") |> dplyr::mutate(y_centroid = -1*y_centroid), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.25,
              legend = legend(pt.cex = 2), asp = 1)
tinyplot::plt(y_centroid ~ x_centroid | celltype_level1, 
              data = adata$obs |> dplyr::filter(core_global == "tma1_F1") |> dplyr::mutate(y_centroid = -1*y_centroid), 
              col = discrete_palette("stallion"),
              pch = 16, cex = 0.25,
              legend = legend(pt.cex = 2), asp = 1)

# Saving
anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered.h5ad", mode = "w")

