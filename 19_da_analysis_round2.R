rm(list = ls())
.rs.restartR(clean = T)

# dir.create("19_da_analysis_round2")

### Using HieraType Results ----------------------------------------------------

pdf(file = "19_da_analysis_round2/boxplots_with_points_new_T_cell_typing.pdf", width = 6, height = 8)

hres <- readRDS(file = "18_hieratype_outs/hieratype_custom_lymphoid_pipeline_results.RDS")
adata <- anndataR::read_h5ad(path = "sgroi-tnbc_filtered.h5ad", mode = "r+")

adata$obs$celltype_h <- adata$obs$celltype_final
adata$obs[hres$post_probs$lymphoidmajor$cell_ID,]$celltype_h <- hres$post_probs$lymphoidmajor$celltype_granular

d <- adata$obs
rm(hres); rm(adata)

dtmp <- table(d$patient, d$celltype_h)[,c("T.CD8.activated", "Treg")]
dtmp <- as.data.frame.matrix(dtmp)
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = d$patient, to = as.character(d$HOXB13_group))
dtmp$ratio <- (dtmp$T.CD8/dtmp$Treg)
dtmp$patient <- rownames(dtmp)

# barplot(ratio ~ patient, data = dtmp, 
#         col = ifelse(test = dtmp$group == "HOXB13+", yes = "lightgreen", no = "darkblue"))
boxplot(ratio ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "N Activated CD8+ T Cells / N Tregs", 
        frame = F, 
        boxwex = 0.5)
stripchart(ratio ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(ratio ~ group, data = dtmp) # p = 0.5251


dtmp <- table(d$patient, d$celltype_h)
dtmp <- dtmp[,-which(colnames(dtmp) == "Epithelial")] |> prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = d$patient, to = as.character(d$HOXB13_group))
dtmp$proportion <- dtmp$T.CD8.activated
dtmp$patient <- rownames(dtmp)

# barplot(proportion ~ patient, data = dtmp, 
#         col = ifelse(test = dtmp$group == "HOXB13+", yes = "lightgreen", no = "darkblue"))
boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "N Activated CD8+ T Cells / Total Non-tumor Cells", 
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # p = 0.3123


dtmp <- table(d$patient, d$celltype_h) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = d$patient, to = as.character(d$HOXB13_group))
dtmp$ratio <- dtmp$T.CD8.activated/dtmp$T.CD8.exhausted
dtmp$patient <- rownames(dtmp)

# barplot(ratio ~ patient, data = dtmp, 
#         col = ifelse(test = dtmp$group == "HOXB13+", yes = "lightgreen", no = "darkblue"))
boxplot(ratio ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "N Activated CD8+ T Cells / N Exhausted CD8+ T Cells", 
        frame = F, 
        boxwex = 0.5)
stripchart(ratio ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(ratio ~ group, data = dtmp) # p = 0.8741

dev.off()


#### Proportions expressing key genes ####

rm(list = ls())
.rs.restartR(clean = T)

# Combine the CD8 T cell clusters, then analyze the following subsets as a % of the overall CD8 T cell population in each sample and compare  HOXB13+ to HOXB13- samples

adata <- anndataR::read_h5ad(path = "sgroi-tnbc_filtered.h5ad", mode = "r+")
hres <- readRDS(file = "18_hieratype_outs/hieratype_custom_lymphoid_pipeline_results.RDS")
adata$obs$celltype_h <- adata$obs$celltype_final
adata$obs[hres$post_probs$lymphoidmajor$cell_ID,]$celltype_h <- hres$post_probs$lymphoidmajor$celltype_granular

idx <- adata$obs_names[grepl(adata$obs$celltype_h, pattern = "T.CD8")]
cts <- adata$layers$counts |> magrittr::set_rownames(adata$obs_names) |> magrittr::set_colnames(adata$var_names)
cts <- cts[idx,]
cts <- (cts > 0)

pdf(file = "19_da_analysis_round2/boxplots_key_marker_genes_proportion_positive_for_CD8_T_cells_new.pdf", width = 6, height = 8)

meta <- adata$obs[idx,]
meta[,c("CXCL13+", "LAG3+", "ENTPD1+", "PDCD1+", "CTLA4+", "TIGIT+", "HAVCR2+", "TNF+", "IFNG+", "PRF1+", "GZMB+", "GZMA+", "ITGAE+", "IL2RA+")] <- 
  cts[,c("CXCL13", "LAG3", "ENTPD1", "PDCD1", "CTLA4", "TIGIT", "HAVCR2", "TNF", "IFNG", "PRF1", "GZMB", "GZMA", "ITGAE", "IL2RA")]

# CXCl13+ and LAG3+ CD8 T cells  (Late Exhausted-subset 1)
dtmp <- table(meta$patient, meta[,c("CXCL13+", "LAG3+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "CXCL13+ LAG3+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.3042

# ENTPD1+, PDCD1 (PD-1)  CD8 T cells (Early exhausted – subset 2)
dtmp <- table(meta$patient, meta[,c("ENTPD1+", "PDCD1+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "ENTPD1+ PDCD1+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.2988

# CTLA4+, TIGIT+,  ENTPD1+ (CD39), LAG3+, TIM3+  CD8 T cells ( Exhausted subset 3)
dtmp <- table(meta$patient, meta[,c("CTLA4+", "TIGIT+", "ENTPD1+", "LAG3+", "HAVCR2+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "CTLA4+ TIGIT+ ENTPD1+ LAG3+ TIM3+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.3491

# TNF+, IFN-G+ CD8 T cells (activated CD8 T cell-subset 1)
dtmp <- table(meta$patient, meta[,c("TNF+", "IFNG+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "TNF+ IFNG+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.04363

# TNF+, IFN-G+, PRFN1+, GZMB+, GZMA+ (activated CD8 T cell-subset 2)
dtmp <- table(meta$patient, meta[,c("TNF+", "IFNG+", "PRF1+", "GZMB+", "GZMA+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "TNF+ IFNG+ PRF1+ GZMB+ GZMA+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.1191

# TNF+, IFN-G+, PRFN1+, GZMB+, GZMA+ ITGAE+ (activated CD8 T cell-subset 3)
dtmp <- table(meta$patient, meta[,c("TNF+", "IFNG+", "PRF1+", "GZMB+", "GZMA+", "ITGAE+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "TNF+ IFNG+ PRF1+ GZMB+ GZMA+ ITGAE+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.09941

# TNF+, IFN-G+, PRFN1+, GZMB+, GZMA+ IL2RA+ (activated CD8 T cell-subset 14)
dtmp <- table(meta$patient, meta[,c("TNF+", "IFNG+", "PRF1+", "GZMB+", "GZMA+", "IL2RA+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "TNF+ IFNG+ PRF1+ GZMB+ GZMA+ IL2RA+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.08486

# See the y-axes! Very low signal for some of these genes... 

dev.off()


### Using previous results -----------------------------------------------------

pdf(file = "19_da_analysis_round2/boxplots_with_points_old_T_cell_typing.pdf", width = 6, height = 8)

d$celltype_simple <- ifelse(test = grepl(pattern = "CD8.activated", x = d$celltype_final), yes = "T.CD8.activated", no = d$celltype_final)
dtmp <- table(d$patient, d$celltype_simple)[,c("T.CD8.activated", "Treg")]
dtmp <- as.data.frame.matrix(dtmp)
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = d$patient, to = as.character(d$HOXB13_group))
dtmp$ratio <- (dtmp$T.CD8.activated/dtmp$Treg)
dtmp$patient <- rownames(dtmp)

# barplot(ratio ~ patient, data = dtmp, 
#         col = ifelse(test = dtmp$group == "HOXB13+", yes = "lightgreen", no = "darkblue"))
boxplot(ratio ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "N Activated CD8+ T Cells / N Tregs", 
        frame = F, 
        boxwex = 0.5)
stripchart(ratio ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(ratio ~ group, data = dtmp) # p = 1


dtmp <- table(d$patient, d$celltype_simple)
dtmp <- dtmp[,-which(colnames(dtmp) == "Epithelial")] |> prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = d$patient, to = as.character(d$HOXB13_group))
dtmp$proportion <- dtmp$T.CD8.activated
dtmp$patient <- rownames(dtmp)

# barplot(proportion ~ patient, data = dtmp, 
#         col = ifelse(test = dtmp$group == "HOXB13+", yes = "lightgreen", no = "darkblue"))
boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "N Activated CD8+ T Cells / Total Non-tumor Cells", 
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # p = 0.458


dtmp <- table(d$patient, d$celltype_simple) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = d$patient, to = as.character(d$HOXB13_group))
dtmp$ratio <- dtmp$T.CD8.activated/dtmp$T.CD8.exhausted
dtmp$patient <- rownames(dtmp)

# barplot(ratio ~ patient, data = dtmp, 
#         col = ifelse(test = dtmp$group == "HOXB13+", yes = "lightgreen", no = "darkblue"))
boxplot(ratio ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "N Activated CD8+ T Cells / N Exhausted CD8+ T Cells", 
        sub = "Note: p61 has zero exhausted CD8 T cells and is not shown.",
        frame = F, 
        boxwex = 0.5)
stripchart(ratio ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(ratio ~ group, data = dtmp) # p = 0.5604

# wilcox.test(ratio ~ group, data = dtmp)$statistic
# wilcox.test(ratio ~ group, data = dtmp[is.finite(dtmp$ratio), ])$statistic

dev.off()


#### Proportions expressing key genes ####

rm(list = ls())
.rs.restartR(clean = T)

# Combine the CD8 T cell clusters 4, 5 and 6. Then analyze the following subsets as a % of the overall CD8 T cell population in each sample and compare  HOXB13+ to HOXB13- samples

adata <- anndataR::read_h5ad(path = "sgroi-tnbc_filtered.h5ad", mode = "r+")
adata$obs$celltype_simple <- dplyr::case_when(grepl(pattern = "CD8", x = adata$obs$celltype_final) ~ "T.CD8", 
                                         grepl(pattern = "CD4", x = adata$obs$celltype_final) ~ "T.CD4", 
                                         grepl(pattern = "Macr", x = adata$obs$celltype_final) ~ "Myeloid",
                                         T ~ adata$obs$celltype_final
)
idx <- adata$obs_names[adata$obs$celltype_simple == "T.CD8"]
cts <- adata$layers$counts |> magrittr::set_rownames(adata$obs_names) |> magrittr::set_colnames(adata$var_names)
cts <- cts[idx,]
cts <- (cts > 0)

pdf(file = "19_da_analysis_round2/boxplots_key_marker_genes_proportion_positive_for_CD8_T_cells.pdf", width = 6, height = 8)

meta <- adata$obs[idx,]
meta[,c("CXCL13+", "LAG3+", "ENTPD1+", "PDCD1+", "CTLA4+", "TIGIT+", "HAVCR2+", "TNF+", "IFNG+", "PRF1+", "GZMB+", "GZMA+", "ITGAE+", "IL2RA+")] <- 
  cts[,c("CXCL13", "LAG3", "ENTPD1", "PDCD1", "CTLA4", "TIGIT", "HAVCR2", "TNF", "IFNG", "PRF1", "GZMB", "GZMA", "ITGAE", "IL2RA")]

# CXCl13+ and LAG3+ CD8 T cells  (Late Exhausted-subset 1)
dtmp <- table(meta$patient, meta[,c("CXCL13+", "LAG3+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "CXCL13+ LAG3+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.2252

# ENTPD1+, PDCD1 (PD-1)  CD8 T cells (Early exhausted – subset 2)
dtmp <- table(meta$patient, meta[,c("ENTPD1+", "PDCD1+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "ENTPD1+ PDCD1+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.2574

# CTLA4+, TIGIT+,  ENTPD1+ (CD39), LAG3+, TIM3+  CD8 T cells ( Exhausted subset 3)
dtmp <- table(meta$patient, meta[,c("CTLA4+", "TIGIT+", "ENTPD1+", "LAG3+", "HAVCR2+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "CTLA4+ TIGIT+ ENTPD1+ LAG3+ TIM3+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.509

# TNF+, IFN-G+ CD8 T cells (activated CD8 T cell-subset 1)
dtmp <- table(meta$patient, meta[,c("TNF+", "IFNG+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "TNF+ IFNG+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.3303

# TNF+, IFN-G+, PRFN1+, GZMB+, GZMA+ (activated CD8 T cell-subset 2)
dtmp <- table(meta$patient, meta[,c("TNF+", "IFNG+", "PRF1+", "GZMB+", "GZMA+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "TNF+ IFNG+ PRF1+ GZMB+ GZMA+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.4373

# TNF+, IFN-G+, PRFN1+, GZMB+, GZMA+ ITGAE+ (activated CD8 T cell-subset 3)
dtmp <- table(meta$patient, meta[,c("TNF+", "IFNG+", "PRF1+", "GZMB+", "GZMA+", "ITGAE+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "TNF+ IFNG+ PRF1+ GZMB+ GZMA+ ITGAE+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.5263

# TNF+, IFN-G+, PRFN1+, GZMB+, GZMA+ IL2RA+ (activated CD8 T cell-subset 14)
dtmp <- table(meta$patient, meta[,c("TNF+", "IFNG+", "PRF1+", "GZMB+", "GZMA+", "IL2RA+")] |> apply(MARGIN = 1, FUN = all)) |> 
  prop.table(margin = 1) |> as.data.frame.matrix()
dtmp$group <- plyr::mapvalues(x = rownames(dtmp), from = meta$patient, to = meta$HOXB13_group |> as.character())
dtmp$proportion <- dtmp$`TRUE`

boxplot(proportion ~ group, data = dtmp, 
        col = c("darkblue", "lightgreen"),
        outline = F, 
        main = "TNF+ IFNG+ PRF1+ GZMB+ GZMA+ IL2RA+",
        frame = F, 
        boxwex = 0.5)
stripchart(proportion ~ group, data = dtmp,
           method = "jitter", jitter = 0.15,
           add = T, vertical = T, 
           pch = 21, cex = 1, bg = "grey")
wilcox.test(proportion ~ group, data = dtmp) # 0.1874

# See the y-axes! Very low signal for some of these genes... 

dev.off()

