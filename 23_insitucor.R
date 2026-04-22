rm(list = ls())
.rs.restartR(clean = T)

adata <- anndataR::read_h5ad("sgroi-tnbc_filtered.h5ad", mode = "r+")
raw <- adata$layers$counts |> magrittr::set_colnames(adata$var_names) |> magrittr::set_rownames(adata$obs_names) |> Matrix::t() |> as('CsparseMatrix')
norm <- Matrix::t(sweep(x = raw, MARGIN = 2, Matrix::colSums(raw), FUN = "/") * 1000)
meta <- adata$obs

neighbs <- InSituType::nearestNeighborGraph(x = meta$x_global, y = meta$y_global, N = 30, subset = meta$core_global)
neighbs@x[neighbs@x > 100] <- 0
neighbs <- Matrix::drop0(neighbs)

##### First, InSituCor for the entire cohort together --------------------------

#renv::install("NanoString-BioStats/InSituCor")
library(InSituCor)
# set.seed(2001)
# res <- insitucor(counts = norm,
#                  conditionon = dplyr::select(meta, celltype_final_amended_1, control_probe_counts, transcript_counts, patient),
#                  celltype = meta$celltype_final_amended_1,
#                  neighbors = neighbs,
#                  tissue = meta$core_global,
#                  xy = dplyr::select(meta, x_global, y_global),
#                  k = 30,
#                  roundcortozero = 0.01,
#                  max_cells = 20000,
#                  # max_module_size = 20, # Not implemented yet!
#                  min_module_size = 3,
#                  resolution = 0.025)
# dir.create("23_insitucor_outs")
# saveRDS(res, "23_insitucor_outs/insitucor_global_results.RDS")
res <- readRDS("23_insitucor_outs/insitucor_global_results.RDS")

library(pheatmap)
hc1 <- hclust(dist(res$condcor))
pheatmap(as.matrix(res$condcor[hc1$order, hc1$order]), cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("darkblue",'blue', "white","red","darkred"))(100),
         breaks = seq(-0.6,0.6,length.out = 101),
         show_rownames = F, show_colnames = F, legend = F, 
         #filename = "23_insitucor_outs/correlation_heatmap.png",
         #width = 6, height = 6, res = 96
         )

#svg(filename = "23_insitucor_outs/insitucor_network.svg", width = 6, height = 7)
#pdf(file = "23_insitucor_outs/insitucor_network.pdf", width = 6, height = 7)
par(mar = c(3, 2, 2, 2) + 0.1)
set.seed(913)
InSituCor::plotCorrelationNetwork(x = (res$condcor * (abs(res$condcor) > 0.1)), 
                                  modules = res$modules, 
                                  show_gene_names = F, 
                                  vertex_size = 4, 
                                  corthresh = 0.2)
rect(ytop = 0.12, ybottom = -0.08, xleft = 0.65, xright = 0.88, col = NA, lwd = 4, border = "red2")
rect(ytop = -0.85, ybottom = -1.08, xleft = 0.05, xright = 0.35, col = NA, lwd = 4, border = "orange2")
text(x = c(0.65 + 0.1, 0.05 + 0.1), y = c(-0.08, -1.08), labels = c("M1 Module", "TLS Module"), pos = 1, col = c("red2", "orange2"))
#dev.off()

res$modules[res$modules$module == "CD79A_PAX5_11",]
pheatmap(res$attributionmats$CD79A_PAX5_11, main = "TLS Module", fontsize_row = 8,
         col = colorRampPalette(c("white", "firebrick"))(100),
         breaks = seq(0,1,length.out=101), name = " ", cellheight = 10, cellwidth = 10, 
         #filename = "23_insitucor_outs/TLS_module_attribution_heatmap.png",
         #width = 6, height = 6
         )

res$modules[res$modules$module == "G0S2_CCL3_9",]
pheatmap(res$attributionmats$G0S2_CCL3_9, main = "M1 Macrophage Module", fontsize_row = 8,
         col = colorRampPalette(c("white", "firebrick"))(100),
         breaks = seq(0,1,length.out=101), name = " ", cellheight = 10, cellwidth = 10, 
         #filename = "23_insitucor_outs/M1_module_attribution_heatmap.png",
         #width = 6, height = 6
         )

# Next steps might be to associate these scores with our banksy niches and/or see
# if they are higher among immune cells in the HOXB13+/- group...





