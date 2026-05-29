rm(list = ls())
.rs.restartR(clean = T)

#renv::install("Nanostring-Biostats/CosMx-Analysis-Scratch-Space/_code/HieraType")
#renv::install("mvtnorm", prompt = F)
#renv::install("bioc::ComplexHeatmap", prompt = F)

#### Testing out HieraType -----------------------------------------------------

# We will start out with the T and NK cells
adata <- anndataR::read_h5ad("sgroi-tnbc_filtered.h5ad", mode = "r+")

# First, we must make a custom pipeline for this gene panel.
library(HieraType)

#!#! NK vs. T cells #!#!
markerslist_lymphoidmajor <- 
  make_markerslist(index_marker = list(
    "NK" = c("NCAM1", "KLRD1", "KLRF1"), 
    "T.cell.Lineage" = c("CD3E", "CD3D", "IL7R", "CD247")
  ), 
  predictors = list(
    "NK" = c(
      "NCAM1", "KLRD1", "KLRF1", "FCGR3A", "PRF1", "GZMB", 
      "GNLY", "KLRC1", "KLRB1", "CD160", "IL2RB", "CD244", 
      "IFNG", "TNF", "CCL5", "CX3CR1", "CD69", "TIGIT", 
      "EOMES", "TBX21", "RUNX3", "GATA3", "CXCR6", "SELL", 
      "NKG7", "SPON2", "FGFBP2", "ZNF683", "CMKLR1", "CTSW"
      ),
    "T.cell.Lineage" = c(
      "CD3D", "CD3E", "IL7R", "CD247", "CD8A", "CD8B", "CD27", "SELL", "CCR7", 
      "CD28", "STAT3", "STAT4", "STAT5B", "RUNX3", "BATF3", "EOMES", "TBX21", 
      "PDCD1", "LAG3", "TIGIT", "HAVCR2", "GZMK", "GZMA", "GZMB", "PRF1", "GNLY", 
      "TCF7", "IFNG", "CCL5", "CX3CR1", "NKG7", "FAS", "FASLG", "IL2", "CD244", 
      "CXCR3", "KLRB1", "TOX", "BCL6", "CXCL13", "CD4", "CD40LG", "GATA3", "RORC", 
      "FOXP3", "IL2RA", "CTLA4", "IL4", "IL5", "IL13", "IL10", "CCR6", "CXCR5", 
      "TNFRSF9", "CD69", "IRF4", "HIF1A", "CXCR6", "IL2RB", "STAT1", "CD2", "CD44", 
      "ID2", "PDCD1LG2", "ENTPD1", "GZMH", "CCR5", "IL2RG", "BCL2", "LCK", 
      "TRAT1", "CD38", "CTSW", "CST7", "ITGAE", "ITGA1", "TAGAP", "RGS1"
    )
  )
)

#!#! Within T cells, CD8 vs. CD4 #!#!
markerslist_cd8_cd4 <- 
  make_markerslist(index_marker = list(
    "T.CD8" = c("CD8A", "CD8B"), 
    "T.CD4" = c("CD4")
  ), 
  predictors = list(
    "T.CD8" = c(
      "CD8A", "CD8B", "CD3E", "CD3D", "CD247", "IL7R", "CD27", "SELL", "CCR7", "TCF7", 
      "CD28", "STAT3", "RUNX3", "EOMES", "TBX21", "PDCD1", "LAG3", "TIGIT", 
      "HAVCR2", "GZMK", "GZMA", "GZMB", "PRF1", "GNLY", "KLRD1", "IFNG", "CCL5", 
      "CX3CR1", "NKG7", "FAS", "FASLG", "IL2", "CD244", "CXCR3", "KLRB1", "TOX", 
      "BCL6", "CXCL13", "STAT5B", "GZMH", "CST7", "CTSW", "FGFBP2", "ITGAE", 
      "ITGA1", "ZNF683", "RGS1", "STAT4", "STAT1", "CD38", "CD69", "TAGAP", "TRAT1",
      "LCK"
    ),
    "T.CD4" = c(
      "CD4", "CD3E", "CD3D", "CD247", "IL7R", "CD27", "SELL", "CCR7", "TCF7", "CD28", "CD40LG", "STAT3", 
      "RUNX3", "BATF3", "GATA3", "TBX21", "RORC", "FOXP3", "IL2RA", "CTLA4", 
      "PDCD1", "LAG3", "TIGIT", "HAVCR2", "IL2", "IFNG", "IL4", "IL5", "IL13", 
      "IL10", "CCL5", "CXCR3", "CCR6", "CXCR5", "BCL6", "EOMES", "CXCL13", "TOX", 
      "IRF4", "HIF1A", "IL2RB", "IL2RG", "CCR4", "TNFRSF9", "CD38", "CD69", 
      "STAT4", "STAT5B", "STAT1", "RGS1"
    )
  )
)

#!#! Within CD4 T cells #!#!
markerslist_cd4custom <- 
  make_markerslist(index_marker = list(
    "Treg" = c("FOXP3", "IL2RA", "CTLA4"), 
    "T.CD4.naive" = c("TCF7", "SELL", "CCR7"), 
    "T.CD4.fh" = c("CXCL13", "CXCR5", "BCL6"), 
    "T.CD4.effector" = c("IFNG", "TBX21", "CXCR3")
  ), 
  predictors = list(
    "Treg" = c(
      "CD3D", "CD3E", "CD4", "FOXP3", "IL2RA", "CTLA4", "LAG3", "CCR4", 
      "TGFBI", "IL10", "TIGIT", "FAS", "IDO1", "CD69", "TNFRSF9", "HIF1A",
      "STAT5B", "ENTPD1", "STAT3", "IL2RB", "IL2RG", "SMAD3", "TGFB1"
      ),
    "T.CD4.naive" = c(
      "CD3D", "CD3E", "CD4", "TCF7", "SELL", "CCR7", "IL7R", "CD28", 
      "CD27", "STAT6", "LTB", "PTPRC", "IL2RB", "TOX", "BCL2", "IL2RG", 
      "STAT3", "JUN", "FOS"
    ), 
    "T.CD4.fh" = c(
      "CD3D", "CD3E", "CD4", "CXCL13", "CXCR5", "BCL6", "PDCD1", "TOX", "IL2RB", 
      "STAT3", "IRF4", "IL10", "IL4", "CCR7", "CD69", "CD40LG", "TIGIT", "LAG3", 
      "CXCR3", "IFNG", "FASLG"
      ), 
    "T.CD4.effector" = c(
      "CD3D", "CD3E", "CD4", "IFNG", "TBX21", "CXCR3", 
      "TNF", "CCL5", "STAT1", "STAT4", "IRF1",
      "GATA3", "IL4", "IL5", "IL13", "CCR4", "STAT6",
      "RORC", "CCR6",
      "GZMB", "GZMK", "PRF1", "GNLY",
      "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CD38", "CD69",
      "ITGAE", "ITGA1", "RGS1", "CX3CR1",
      "CD44", "FAS", "FASLG", "TNFRSF9", "HIF1A"
    )
  )
)

#!#! Within CD8 T cells !#!# 
markerslist_cd8custom <- 
  make_markerslist(index_marker = list(
    "T.CD8.activated" = c("GZMB", "PRF1", "IFNG", "CD38", "GNLY", "KLRD1"),
    "T.CD8.exhausted" = c("TOX", "PDCD1", "HAVCR2", "LAG3", "TIGIT")
  ), 
  predictors = list(
    "T.CD8.activated" = c(
      "CD3D", "CD3E", "CD8A", "CD8B", "GZMB", "PRF1", "IFNG", "CD38", "GNLY", 
      "KLRD1", "GZMA", "GZMH", "GZMK", "NKG7", "CST7", "CTSW", "FGFBP2", "TNF", 
      "CCL5", "FASLG", "CX3CR1", "KLRB1", "CD69", "STAT1", "STAT4", "ITGAE", "ITGA1"
    ),
    "T.CD8.exhausted" = c(
      "CD3D", "CD3E", "CD8A", "CD8B", "TOX", "PDCD1", "HAVCR2", 
      "LAG3", "TIGIT", "BATF3", "IRF4", "RUNX3", "BCL6", "EOMES", 
      "CXCL13", "CD244", "ID2", "ENTPD1", "STAT3")
  )
)

# Now we can make the pipeline
pipeline_lymphoid <- 
  make_pipeline(markerslists = list("lymphoidmajor" = markerslist_lymphoidmajor,
                                    "tmajor" = markerslist_cd8_cd4,
                                    "tcd4minor" = markerslist_cd4custom,
                                    "tcd8minor" = markerslist_cd8custom
  ), 
  priors = list("tmajor" = "lymphoidmajor", 
                "tcd4minor" = "tmajor", 
                "tcd8minor" = "tmajor"
  ), 
  priors_category = list("tmajor" = "T.cell.Lineage", 
                         "tcd4minor" = "T.CD4", 
                         "tcd8minor" = "T.CD8"
  )
)

# Finally, we can run HieraType, using this pipeline. We need the neighborhood
# graph from expression space.
library(BPCells)
snn <- knn_hnsw(adata$obsm$X_scVI, k = 30, metric = "euclidean", ef = 1500) |> # Find approximate nearest neighbors
  knn_to_snn_graph() # Convert to a SNN graph
dimnames(snn) <- list(adata$obs_names, adata$obs_names)

library(Matrix)
TandNKcells <- adata$obs[adata$obs$celltype_level2 == "T and NK" & !is.na(adata$obs$celltype_level2),] |> rownames()
cts <- adata$layers$counts |> as("CsparseMatrix") |> magrittr::set_colnames(value = adata$var_names) |> magrittr::set_rownames(value = adata$obs_names)

lymphoid_typing <- 
  run_pipeline(pipeline = pipeline_lymphoid, 
               counts_matrix = cts[TandNKcells, ], 
               adjacency_matrix = snn[TandNKcells, TandNKcells], 
               celltype_call_threshold = 0.5
)

# dir.create("18_hieratype_outs")
# saveRDS(pipeline_lymphoid, file = "18_hieratype_outs/hieratype_custom_lymphoid_pipeline.RDS")
# saveRDS(lymphoid_typing, file = "18_hieratype_outs/hieratype_custom_lymphoid_pipeline_results.RDS")

# Visualizing
lymphoid_typing <- readRDS("18_hieratype_outs/hieratype_custom_lymphoid_pipeline_results.RDS")
norm <- (adata$layers$lognorm) |> magrittr::set_colnames(rownames(adata$var)) |> magrittr::set_rownames(rownames(adata$obs)) |> as("CsparseMatrix")
norm <- norm[TandNKcells,]
fct <- clusterwise_foldchange_metrics(normed = Matrix::t(norm),
                                      metadata = lymphoid_typing$post_probs$lymphoidmajor,
                                      cluster_column = "celltype_granular")
idxs <- unique(unname(c(unlist(lapply(markerslist_cd4custom, "[[", "index_marker")), 
                        unlist(lapply(markerslist_cd8custom, "[[", "index_marker")), 
                        unlist(lapply(markerslist_lymphoidmajor, "[[", "index_marker")))))
nms <- c("NK", names(markerslist_cd4custom), names(markerslist_cd8custom))
hmsubtype <- marker_heatmap(fct, featsuse = c("CD8A", "CD8B", "CD4", "ENTPD1", "TNF", "GZMA", "GZMK", idxs[is.element(idxs, fct$gene)]), 
                            clusterorder = nms, orient_diagonal = T) + 
  ggplot2::scale_fill_viridis_c(option = "C") + 
  ggplot2::labs(fill = "Mean expression in group", title = "T cells")

print(hmsubtype)
#ggplot2::ggsave("18_hieratype_outs/t_cells_bubble_plot_canonical_markers.pdf", width = 10, height = 4)

table(adata$obs[TandNKcells,]$patient, lymphoid_typing$post_probs$lymphoidmajor$celltype_granular)


#### Using HieraType for Macrophages -------------------------------------------

adata <- anndataR::read_h5ad(path = "sgroi-tnbc_filtered.h5ad", mode = "r+")
adata$var_names

library(HieraType)

## Claude helped me with the marker genes ##

#!#! Macrophage vs. DC (top layer) #!#!
markerslist_myeloidmajor <-
  make_markerslist(index_marker = list(
    "Macrophage" = c("CD68", "CD163", "C1QB", "MARCO", "CSF1R", "VSIG4"),
    "DC"         = c("CD1C", "LAMP3", "ITGAX", "IRF8", "BATF3", "IRF4")
  ),
  predictors = list(
    "Macrophage" = c(
      "CD68", "CD163", "C1QB", "MARCO", "CSF1R", "VSIG4", "MPEG1",
      "AIF1", "TREM2", "SPP1", "APOE", "AXL", "ITGAM",
      "TGFB1", "IL10", "VEGFA", "MMP9", "MMP12",
      "LGALS3", "CTSB", "CTSD", "CTSL", "CTSS",
      "CCL7", "CCL8", "CCL18", "CCL13",
      "ARG1", "NOS2", "CXCL10", "TNF", "IL1B", "IL6",
      "CIITA", "STAT1", "IRF1", "MKI67", "CHIT1", "PLA2G7", "ALOX5AP"
    ),
    "DC" = c(
      "CD1C", "LAMP3", "ITGAX", "CLEC10A", "FCER1A", "IRF8", "BATF3", "IRF4",
      "CLEC12A", "FCGR2B", "LILRB2", "LILRB4", "NOTCH2", "RUNX3",
      "CCR7", "CD83", "CD40", "CD80", "CD86", "CIITA",
      "MARCKSL1", "IDO1", "CXCL9", "CXCL10", "PDCD1LG2",
      "SPIB", "RELB", "AREG", "TGFB1", "CX3CR1", "ITGAM", "CD1A"
    )
  )
  )

#!#! Within DCs: cDC vs. mregDC #!#!
markerslist_dcs <-
  make_markerslist(index_marker = list(
    "cDC"    = c("CD1C", "IRF8", "BATF3", "IRF4", "CLEC10A", "FCER1A"),
    "mregDC" = c("LAMP3", "CCR7", "CD83", "MARCKSL1", "IDO1")
  ),
  predictors = list(
    "cDC" = c(
      "CD1C", "CLEC10A", "FCER1A", "IRF4", "BATF3", "IRF8",
      "NOTCH2", "RUNX3", "CLEC12A", "FCGR2B",
      "LILRB2", "LILRB4", "CD1A", "AREG", "CX3CR1",
      "ITGAX", "ITGAM", "CIITA", "SPIB", "RELB", "CD74"
    ),
    "mregDC" = c(
      "LAMP3", "CCR7", "CD83", "MARCKSL1", "IDO1",
      "CCL19", "CCL21", "CD274", "PDCD1LG2",
      "CXCL9", "CXCL10", "STAT3", "CD40", "CD80", "CD86",
      "CIITA", "AIRE", "TGFB1", "AREG", "IL4R"
    )
  )
  )

#!#! Within Macrophages: M1 vs. M2 #!#!
markerslist_macros <-
  make_markerslist(index_marker = list(
    "Mac.M1" = c("CXCL10", "CXCL9", "NOS2", "STAT1", "IL12B", "IRF1"),
    "Mac.M2" = c("CD163", "ARG1", "VSIG4", "TREM2", "CCL18", "IL10")
  ),
  predictors = list(
    "Mac.M1" = c(
      "CXCL10", "CXCL9", "CXCL11", "CXCL5",
      "TNF", "IL1B", "IL1A", "IL6", "IL12A", "IL12B", "IL18",
      "STAT1", "STAT2", "IRF1", "IRF3", "IRF5", "NOS2", "IDO1",
      "CD80", "CD86", "ICAM1", "ISG15", "IFIT2", "IFIT3", "MX1", "MX2", "OAS2",
      "IFNG", "IFNGR1", "CGAS", "TMEM173", "NFKB1", "NFKB2", "RELB",
      "TNFAIP3", "SOCS1", "SOCS3", "HIF1A", "SERPINE1", "CCL20"
    ),
    "Mac.M2" = c(
      "CD163", "ARG1", "IL10", "TGFB1", "TGFBI", "SMAD3",
      "CCL13", "CCL18", "CCL2", "CCL7", "CCL8", "CCL26",
      "MARCO", "VSIG4", "APOE", "TREM2", "SPP1",
      "FN1", "VEGFA", "MMP9", "MMP12", "MMP14",
      "LGALS3", "CHIT1", "CTSL", "CTSS", "CTSB",
      "IL4R", "STAT3", "STAT6", "SOCS2", "AXL",
      "ALOX5AP", "PLA2G7", "PTGS1", "MYD88", "MPEG1", "C1QB"
    )
  )
  )

#!#! Build pipeline #!#!
pipeline_myeloid <-
  make_pipeline(
    markerslists = list(
      "myeloidmajor" = markerslist_myeloidmajor,
      "dcminor"      = markerslist_dcs,
      "macminor"     = markerslist_macros
    ),
    priors = list(
      "dcminor"  = "myeloidmajor",
      "macminor" = "myeloidmajor"
    ),
    priors_category = list(
      "dcminor"  = "DC",
      "macminor" = "Macrophage"
    )
  )

# library(BPCells)
# snn <- knn_hnsw(adata$obsm$X_scVI, k = 30, metric = "euclidean", ef = 1500) |> # Find approximate nearest neighbors
#   knn_to_snn_graph() # Convert to a SNN graph
# dimnames(snn) <- list(adata$obs_names, adata$obs_names)
# 
# library(Matrix)
# myeloid_cells <- adata$obs[adata$obs$celltype_level2 %in% c("Macrophage1", "Macrophage2", "mDC") & !is.na(adata$obs$celltype_level2),] |> rownames()
# cts <- adata$layers$counts |> as("CsparseMatrix") |> magrittr::set_colnames(value = adata$var_names) |> magrittr::set_rownames(value = adata$obs_names)
# 
# myeloid_typing <- 
#   run_pipeline(pipeline = pipeline_myeloid, 
#                counts_matrix = cts[myeloid_cells, ], 
#                adjacency_matrix = snn[myeloid_cells, myeloid_cells], 
#                celltype_call_threshold = 0.5
#   )

#saveRDS(pipeline_myeloid, file = "18_hieratype_outs/hieratype_custom_myeloid_pipeline.RDS")
#saveRDS(myeloid_typing, file = "18_hieratype_outs/hieratype_custom_myeloid_pipeline_results.RDS")

# Visualizing
myeloid_typing <- readRDS("18_hieratype_outs/hieratype_custom_myeloid_pipeline_results.RDS")
norm <- (adata$layers$lognorm) |> magrittr::set_colnames(rownames(adata$var)) |> magrittr::set_rownames(rownames(adata$obs)) |> as("CsparseMatrix")
norm <- norm[myeloid_cells,]
fct <- clusterwise_foldchange_metrics(normed = Matrix::t(norm),
                                      metadata = myeloid_typing$post_probs$myeloidmajor,
                                      cluster_column = "celltype_granular")
idxs <- unique(unname(c(unlist(lapply(markerslist_macros, "[[", "index_marker")), 
                        unlist(lapply(markerslist_dcs, "[[", "index_marker")), 
                        unlist(lapply(markerslist_myeloidmajor, "[[", "index_marker")))))
nms <- c(names(markerslist_macros), names(markerslist_dcs))
hmsubtype <- marker_heatmap(fct, featsuse = c("SPP1", "TGFB1", "CIITA", "CXCL1", 
                                              "CXCL2", "CXCL5", "CXCL13", "CCL20", 
                                              "CCL26", "CCL19", "CD40", "PDCD1", 
                                              idxs[is.element(idxs, fct$gene)]), 
                            clusterorder = nms, orient_diagonal = T) + 
  ggplot2::scale_fill_viridis_c(option = "C") + 
  ggplot2::labs(fill = "Mean expression in group", title = "Myeloid cells")

print(hmsubtype)
#ggplot2::ggsave("18_hieratype_outs/myeloid_cells_bubble_plot_canonical_markers.pdf", width = 10, height = 4)

table(adata$obs[myeloid_cells,]$patient, myeloid_typing$post_probs$myeloidmajor$celltype_granular)

