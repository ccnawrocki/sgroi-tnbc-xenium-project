rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("sgroi-tnbc_filtered.h5ad", mode = "r+")

# Let's just look a bit more broadly at the niches
par(mar = c(1, 1, 1, 1))
tinyplot::plt(y_global ~ x_global | banksy_niche, data = adata$obs, pch = ".", cex = 0.02, palette = "Alphabet", asp = 1, 
              legend = list(pt.cex = 10), 
              frame = F, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
d <- adata$obs |> dplyr::summarise(X = median(x_global), Y = median(y_global), core = unique(core), .by = core_global)
text(x = d$X, y = d$Y, label = d$core_global)

D <- adata$obs
D$y_centroid <- -1*D$y_centroid
tinyplot::plt(y_centroid ~ x_centroid | banksy_niche, data = D[D$core_global == "tma2_E1",], 
              palette = "Alphabet", pch = 16, cex = 0.5, asp = 1, legend = list(pt.cex = 2), 
              frame = F, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, main = "tma2_E1")
tinyplot::plt(y_centroid ~ x_centroid | celltype_final_amended_1, data = D[D$core_global == "tma2_E1",], 
              palette = "Alphabet", pch = 16, cex = 0.4, asp = 1, legend = list(pt.cex = 2), 
              frame = F, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, main = "tma2_E1")

# niche 6 = TLS niche?
# niche 4 = innate immune niche? 

tinyplot::plt(y_centroid ~ x_centroid | banksy_niche, data = D[D$core_global == "tma3_E1",], 
              palette = "Alphabet", pch = 16, cex = 0.5, asp = 1, legend = list(pt.cex = 2), 
              frame = F, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, main = "tma3_E1")
tinyplot::plt(y_centroid ~ x_centroid | celltype_final_amended_1, data = D[D$core_global == "tma3_E1",], 
              palette = "Alphabet", pch = 16, cex = 0.4, asp = 1, legend = list(pt.cex = 2), 
              frame = F, xaxt = "n", yaxt = "n", xlab = NA, ylab = NA, main = "tma3_E1")

table(adata$obs$core_global, adata$obs$celltype_final_amended_1)

# At a glance, it appears that some samples have lots of T.CD4.fh cells. This 
# could be an interesting thing to focus on. 

# Maybe the niches themselves, and their combinations, are important. Does 
# niche 6 confer a survival benefit or response benefit? What about when 
# there is a combo of 6 and 4? 

# Characterizing the niches by their cell types

# First, overall proportions
props <- table(adata$obs$celltype_final_amended_1, adata$obs$banksy_niche) |> 
  prop.table(margin = 2) |> as.data.frame.matrix() |> dplyr::select(-singleton) 
pheatmap::pheatmap(props, scale = "row")

# Better: average proportion
D$celltype_final_amended_1 <- as.factor(D$celltype_final_amended_1)
D <- D[D$banksy_niche != "singleton",]
core_props <- dplyr::group_by(D, core_global, banksy_niche, celltype_final_amended_1, .drop = F) |> 
  dplyr::tally() |> 
  dplyr::group_by(core_global, banksy_niche, .drop = T) |> 
  dplyr::mutate(total = sum(n)) |> 
  dplyr::mutate(prop = n/total)
core_props <- na.omit(core_props)

mean_props <- lapply(X = as.character(1:11), FUN = function(nn) {
  niche <- core_props[core_props$banksy_niche == nn,] |> 
    tidyr::pivot_wider(values_from = "prop", names_from = "core_global", 
                       id_cols = c("celltype_final_amended_1"), values_fill = 0) |> 
    tibble::column_to_rownames("celltype_final_amended_1") |> 
    as.matrix()
  data.frame(mean.prop = rowMeans(niche)) |> t() |> as.data.frame() |> magrittr::set_rownames(NULL)
})
mean_props <- dplyr::bind_rows(mean_props, .id = "banksy_niche")

scaled <- sweep(mean_props[,-1], 2, apply(mean_props[,-1], 2, min), FUN = "-")
scaled <- sweep(scaled, 2, apply(scaled, 2, max), FUN = "/")

dir.create("22_banksy_downstream_analysis_outs")
pheatmap::pheatmap(scaled, scale = "none", cellwidth = 15, cellheight = 15, 
                   color = viridis::inferno(n=101), 
                   #filename = "22_banksy_downstream_analysis_outs/celltype_enrichment_heatmap.png",
                   #width = 8, height = 6
                   )

# Characterizing the niches by their marker genes
library(presto)
library(data.table)
library(lme4)
library(lmerTest)
library(singlecellmethods)
psb <- collapse_counts(counts_mat = adata$layers$counts |> magrittr::set_colnames(adata$var_names) |> magrittr::set_rownames(adata$obs_names) |> Matrix::t() |> as('CsparseMatrix'), 
                       meta_data = adata$obs, varnames = c("core_global", "banksy_niche"), min_cells_per_group = 20, get_norm = T)
wilout <- wilcoxauc(X = psb$exprs_norm, y = psb$meta_data$banksy_niche)
top_marks <- wilout |> dplyr::group_by(group) |> dplyr::top_n(n = 10, wt = logFC) |> 
  dplyr::arrange(group, desc(logFC))
lfcs <- tidyr::pivot_wider(data = wilout[wilout$feature %in% top_marks$feature,], 
                           id_cols = "feature", names_from = "group", values_from = "logFC") |> 
  tibble::column_to_rownames("feature") |> as.matrix()
pheatmap::pheatmap(lfcs, scale = "none", cellwidth = 15, cellheight = 8, 
                   color = viridis::inferno(n=101), breaks = seq(-1, 2, length.out = 101), 
                   treeheight_col = 25, fontsize_row = 8, 
                   name = "log2FC",
                   #filename = "22_banksy_downstream_analysis_outs/gene_enrichment_heatmap.png",
                   #width = 6, height = 12
                   )

# Next step: we could dig deeper into the individual cell types. For example, 
# macrophages in niche 4 vs all other macrophages. Or, T cells in niche 11 vs 
# T cells in niche 6. 

# Characterizing which niches interface with one another. 








