rm(list = ls())
.rs.restartR(clean = T)

# This is the finalized DE. Dennis needed to confirm the patient groupings. We 
# also want to ensure that some of the results are not impacted by contamination.

# dir.create("14_de_analysis_final_outs")
adata <- anndataR::read_h5ad(path = "sgroi-tnbc_filtered.h5ad", mode = "r+")

# HOXB13 groupings
sample_data <- openxlsx::read.xlsx(xlsxFile = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/xenium/sgroi-tnbc-data/Normalized HOXB13.xlsx")
sample_data$patient <- gsub(pattern = "([a-z])|([A-Z])|( )|(-)", replacement = "", x = sample_data$Sample_Name) |> sprintf(fmt = "%02s") |> paste0("p", ... = _)
sample_data$group <- ifelse(test = sample_data$H_norm > -13.5, yes = "HOXB13+", no = "HOXB13-")

# qPCR on the samples
tmp <- dplyr::arrange(sample_data, H_norm)
barplot(tmp$H_norm, col = ifelse(test = tmp$group == "HOXB13+", yes = "lightgreen", no = "darkblue"))
abline(h = -10, col = "black", lwd = 1, lty = "dashed")
legend(x = 25, y = -15, legend = c("+", "-"), fill = c("lightgreen", "darkblue"), pch = NA, title = "HOXB13 group", adj = c(0.5, 0.5))


# qPCR on the patients
patient_data <- dplyr::group_by(sample_data, patient) |> 
  dplyr::summarise(H_norm = mean(H_norm)) |> 
  dplyr::arrange(H_norm)
patient_data$group <- ifelse(test = patient_data$H_norm > -13.5, yes = "HOXB13+", no = "HOXB13-")
par(las = 1)
mps <- barplot(patient_data$H_norm, 
        col = ifelse(test = patient_data$group == "HOXB13+", yes = "lightgreen", no = "darkblue"), 
        border = "black", 
        main = "HOXB13 Normalized qPCR by Patient"
        )
abline(h = -13.5, col = "black", lwd = 1, lty = "dashed")
mtext(text = patient_data$patient, at = mps, side = 3, cex = 0.7, col = ifelse(test = patient_data$patient %in% levels(adata$obs$patient), yes = "red", no = "black"))
legend(x = 25, y = -15, legend = c("+", "-"), fill = c("lightgreen", "darkblue"), pch = NA, title = "HOXB13 group", adj = c(0.5, 0.5))
legend(x = 20, y = -15, legend = c("yes", "no"), col = c("red", "black"), pch = "a", title = "Xenium", adj = c(0, 0.5))

setdiff(x = levels(adata$obs$patient), y = patient_data$patient)

patient_data_filtered <- patient_data[patient_data$patient %in% levels(adata$obs$patient),]
# pdf("14_de_analysis_final_outs/HOXB13_qPCR_by_patient.pdf", width = 10, height = 6)
mps <- barplot(patient_data_filtered$H_norm, 
               col = ifelse(test = patient_data_filtered$group == "HOXB13+", yes = "lightgreen", no = "darkblue"), 
               border = "black", 
               main = "HOXB13 Normalized qPCR by Patient"
)
abline(h = -13.5, col = "black", lwd = 1, lty = "dashed")
mtext(text = patient_data_filtered$patient, at = mps, side = 3, cex = 0.7, col = "red")
legend(x = 20, y = -15, legend = c("+", "-"), fill = c("lightgreen", "darkblue"), pch = NA, title = "HOXB13 group", adj = c(0.5, 0.5))
# dev.off()

wilcox.test(patient_data$H_norm ~ patient_data$group)
wilcox.test(patient_data_filtered$H_norm ~ patient_data_filtered$group)

# Adding to the adata
adata$obs$HOXB13_group <- plyr::mapvalues(x = adata$obs$patient, from = patient_data_filtered$patient, to = patient_data_filtered$group)


## DE 
# Epithelial cells -------------------------------------------------------------

# Need to identify genes that are likely very contaminated. For every epithelial 
# cell, we will find it's spatial neighbors using Delaunay. Then, among these
# neighbors that are not epithelial, we will calculate the average expression
# of every gene.

source(file = url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/refs/heads/Main/_code/smiDE/R/pre_de_contamination_screen.R"))
coords <- split(x = adata$obs[,c("x_centroid", "y_centroid")], f = adata$obs$core_global)
ids <- split(x = rownames(adata$obs), f = adata$obs$core_global)

library(Matrix)
# renv::install("geometry")
neighbs <- lapply(X = coords, FUN = spatula::getSpatialNeighbors, return_weights = F, dist_thresh_quantile = 0.98)
neighbs <- Matrix::bdiag(neighbs)
dimnames(neighbs) <- list(unlist(ids), unlist(ids))

cts <- adata$layers$counts |> t() |> as("CsparseMatrix") |> magrittr::set_colnames(adata$obs_names) |> magrittr::set_rownames(adata$var_names)
meta <- adata$obs
meta$cell_ID <- rownames(meta)

library(data.table)
contam <- overlap_ratio_metric(assay_matrix = cts, 
                               metadata = meta, 
                               adjacency_matrix = neighbs, 
                               cluster_col = "celltype_level1", 
                               cellid_col = "cell_ID", 
                               verbose = T)
good_targets <- contam[contam$celltype_level1 == "Epithelial" & contam$ratio <= 2,] |> dplyr::pull(target)

idx <- meta |> dplyr::filter(celltype_level1 == "Epithelial") |> dplyr::arrange(patient) |> rownames()
cts <- cts[good_targets, idx]
meta <- meta[idx,]

mm <- model.matrix(~HOXB13_group, data = meta)
eff <- meta$total_counts
ids <- meta$patient

nfit <- nebula::nebula(count = cts, id = ids, pred = mm, offset = eff, method = "LN", covariance = T, ncore = 6)
saveRDS(nfit, "14_de_analysis_final_outs/epithelial_de_nebula_model.RDS")

nfit <- readRDS("14_de_analysis_final_outs/epithelial_de_nebula_model.RDS")

local_wald_test <- function(nfit, .contr) {
  
  # Get lfc values
  lfcs <- rowSums(sweep(x = nfit$summary[, 1:length(.contr)], MARGIN = 2, STATS = .contr, FUN = "*"))
  
  # Get ses, wald stat, and p-values
  ses <- array(data = NA, dim = c(length(lfcs)))
  wss <- array(data = NA, dim = c(length(lfcs)))
  pvals <- array(data = NA, dim = c(length(lfcs)))
  for (i in 1:length(pvals)) {
    cov <- matrix(NA, length(.contr), length(.contr))
    cov[lower.tri(cov, diag=T)] <- as.numeric(nfit$covariance[i,])
    cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
    ses[i] <- sqrt((t(.contr)%*%cov%*%(.contr)))
    wss[i] <- lfcs[i]^2/(t(.contr)%*%cov%*%(.contr))
    pvals[i] <- pchisq(lfcs[i]^2/(t(.contr)%*%cov%*%(.contr)),1,lower.tail=FALSE)
  }
  
  # Construct output
  out <- data.frame("logFC" = lfcs,
                    "log2FC" = lfcs / log(2),
                    "se" = ses,
                    "wald.stat" = wss,
                    "p.value" = pvals, 
                    "target" = nfit$summary$gene, 
                    "convergence" = nfit$convergence
  )
  
  return(out)
  
}

contr <- c(0, -1)

out <- local_wald_test(nfit = nfit, .contr = contr)
out <- out[out$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
out$p.adj <- p.adjust(out$p.value, method = "BH") # FDR

plot(x = out$log2FC, 
     y = -log10(out$p.adj), 
     pch = 16, cex = 0.5, frame = F, 
     col = ifelse(test = (out$p.adj < 0.01 & abs(out$log2FC) > 1), yes = "red", no = "black"), 
     xlab = "log2FC", ylab = "-log10(adj. P-value)")

library(ggplot2)
out$enrichment <- dplyr::case_when(out$log2FC > 1 & out$p.adj < 0.05 ~ "HOXB13+", 
                                   out$log2FC < -1 & out$p.adj < 0.05 ~ "HOXB13-", 
                                   T ~ NA)
ggplot() + 
  geom_point(data = out[out$log2FC > 1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), fill = enrichment), shape = 21, size = 2) +
  geom_point(data = out[out$log2FC < -1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), fill = enrichment), shape = 21, size = 2) +
  scale_fill_manual(values = c("HOXB13+" = "lightgreen", "HOXB13-" = "darkblue")) + 
  geom_point(data = out[abs(out$log2FC) < 1 | out$p.adj >= 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj)), color = "grey", shape = 16, size = 1) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) + 
  ggrepel::geom_text_repel(data = out[abs(out$log2FC) > 1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), label = target), size = 3) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5), 
        legend.position = "top"
  ) + 
  labs(title = "Epithelial Cells", fill = "Enrichment") + 
  guides(color = guide_legend(override.aes = list(size = 4)))
# ggsave(filename = "14_de_analysis_final_outs/epithelial_de_nebula_volcano_plot.pdf", width = 8, height = 8)
# openxlsx::write.xlsx(x = out |> dplyr::arrange(enrichment, desc(wald.stat)), file = "14_de_analysis_final_outs/epithelial_de_nebula_results_file.xlsx")

mean(out$log2FC > 0) # 0.5418848

test <- data.frame(target = nfit$summary$gene, A = nfit$summary$`logFC_(Intercept)`)
rownames(test) <- test$target

# pdf(file = "14_de_analysis_final_outs/epithelial_de_nebula_MA_plot.pdf", width = 7, height = 6)
plot(x = test[out$target,]$A, y = out$log2FC, pch = 21, xlab = "Avg. Expression", ylab = "log2FC", cex = 1, 
     bg = dplyr::case_when(is.na(out$enrichment) ~ "grey", out$enrichment == "HOXB13+" ~ "lightgreen", out$enrichment == "HOXB13-" ~ "darkblue"))
abline(h = 0, col = "red", lwd = 2)
text(x = test[out[!is.na(out$enrichment),]$target,]$A, y = out[!is.na(out$enrichment),]$log2FC, labels = out[!is.na(out$enrichment),]$target, 
     pos = sample(x = 1:4, size = sum(!is.na(out$enrichment)), replace = T), cex = 0.5)
# dev.off()

norm <- adata$layers$lognorm |> t() |> as("CsparseMatrix") |> magrittr::set_colnames(adata$obs_names) |> magrittr::set_rownames(adata$var_names)
# pdf("14_de_analysis_final_outs/epithelial_IDO1_norm_expr_boxplot.pdf", width = 4, height = 8)
boxplot(norm["IDO1", idx] ~ meta$HOXB13_group, outline = F, frame = F, col = c("lightgreen", "darkblue"), 
        xlab = NA, ylab = "Normalized Expression", main = "IDO1 in Epithelial Cells")
# dev.off()

test <- data.frame(ido1 = norm["IDO1", idx], group = meta$HOXB13_group)
ggplot(data = test) + 
  scattermore::geom_scattermore(mapping = aes(x = group, y = ido1), position = position_jitter(width = 0.1, height = 0)) + 
  stat_summary(mapping = aes(x = group, y = ido1, color = group), geom = "crossbar", fun = "mean", width = 0.25) + 
  scale_color_manual(values = c("HOXB13+" = "lightgreen", "HOXB13-" = "darkblue")) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5), 
        legend.position = "top", 
        axis.title.x = element_blank()
  ) + 
  labs(title = "IDO1 Expression in Epithelial Cells", color = "Group", y = "Normalized Expression")
# ggsave("14_de_analysis_final_outs/epithelial_IDO1_norm_expr_means.pdf", width = 4, height = 8)


# T.CD4.naive cells ------------------------------------------------------------
source(file = url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/refs/heads/Main/_code/smiDE/R/pre_de_contamination_screen.R"))
coords <- split(x = adata$obs[,c("x_centroid", "y_centroid")], f = adata$obs$core_global)
ids <- split(x = rownames(adata$obs), f = adata$obs$core_global)

library(Matrix)
neighbs <- lapply(X = coords, FUN = spatula::getSpatialNeighbors, return_weights = F, dist_thresh_quantile = 0.98)
neighbs <- Matrix::bdiag(neighbs)
dimnames(neighbs) <- list(unlist(ids), unlist(ids))

cts <- adata$layers$counts |> t() |> as("CsparseMatrix") |> magrittr::set_colnames(adata$obs_names) |> magrittr::set_rownames(adata$var_names)
meta <- adata$obs
meta$cell_ID <- rownames(meta)

library(data.table)
contam <- overlap_ratio_metric(assay_matrix = cts, 
                               metadata = meta, 
                               adjacency_matrix = neighbs, 
                               cluster_col = "celltype_final", 
                               cellid_col = "cell_ID", 
                               verbose = T)
good_targets <- contam[contam$celltype_final == "T.CD4.naive" & contam$ratio <= 2,] |> dplyr::pull(target)

idx <- meta |> dplyr::filter(celltype_final == "T.CD4.naive") |> dplyr::arrange(patient) |> rownames()
cts <- cts[good_targets, idx]
meta <- meta[idx,]

mm <- model.matrix(~HOXB13_group, data = meta)
eff <- meta$total_counts
ids <- meta$patient

nfit <- nebula::nebula(count = cts, id = ids, pred = mm, offset = eff, method = "LN", covariance = T, ncore = 4)
saveRDS(nfit, "14_de_analysis_final_outs/t.cd4.naive_de_nebula_model.RDS")

nfit <- readRDS("14_de_analysis_final_outs/t.cd4.naive_de_nebula_model.RDS")

out <- local_wald_test(nfit = nfit, .contr = contr)
out <- out[out$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
out$p.adj <- p.adjust(out$p.value, method = "BH") # FDR

plot(x = out$log2FC, 
     y = -log10(out$p.adj), 
     pch = 16, cex = 0.5, frame = F, 
     col = ifelse(test = (out$p.adj < 0.05), yes = "red", no = "black"), 
     xlab = "log2FC", ylab = "-log10(adj. P-value)")

library(ggplot2)
out$enrichment <- dplyr::case_when(out$log2FC > 0 & out$p.adj < 0.05 ~ "HOXB13+", 
                                   out$log2FC < 0 & out$p.adj < 0.05 ~ "HOXB13-", 
                                   T ~ NA)
ggplot() + 
  geom_point(data = out[out$log2FC > 0 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), fill = enrichment), shape = 21, size = 2) +
  geom_point(data = out[out$log2FC < 0 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), fill = enrichment), shape = 21, size = 2) +
  scale_fill_manual(values = c("HOXB13+" = "lightgreen", "HOXB13-" = "darkblue")) + 
  geom_point(data = out[out$p.adj >= 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj)), color = "grey", shape = 16, size = 1) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) + 
  ggrepel::geom_text_repel(data = out[abs(out$log2FC) > 0 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), label = target), size = 3) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5), 
        legend.position = "top"
  ) + 
  labs(title = "Naive CD4+ T Cells", fill = "Enrichment") + 
  guides(color = guide_legend(override.aes = list(size = 4)))
# ggsave(filename = "14_de_analysis_final_outs/t.cd4.naive_de_nebula_volcano_plot.pdf", width = 8, height = 8)
# openxlsx::write.xlsx(x = out |> dplyr::arrange(enrichment, desc(wald.stat)), file = "14_de_analysis_final_outs/t.cd4.naive_de_nebula_results_file.xlsx")

mean(out$log2FC > 0) # 0.4272727


# Macrophages ------------------------------------------------------------------
source(file = url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/refs/heads/Main/_code/smiDE/R/pre_de_contamination_screen.R"))
coords <- split(x = adata$obs[,c("x_centroid", "y_centroid")], f = adata$obs$core_global)
ids <- split(x = rownames(adata$obs), f = adata$obs$core_global)

library(Matrix)
neighbs <- lapply(X = coords, FUN = spatula::getSpatialNeighbors, return_weights = F, dist_thresh_quantile = 0.98)
neighbs <- Matrix::bdiag(neighbs)
dimnames(neighbs) <- list(unlist(ids), unlist(ids))

cts <- adata$layers$counts |> t() |> as("CsparseMatrix") |> magrittr::set_colnames(adata$obs_names) |> magrittr::set_rownames(adata$var_names)
meta <- adata$obs
meta$cell_ID <- rownames(meta)

library(data.table)
contam <- overlap_ratio_metric(assay_matrix = cts, 
                               metadata = meta, 
                               adjacency_matrix = neighbs, 
                               cluster_col = "celltype_final", 
                               cellid_col = "cell_ID", 
                               verbose = T)
good_targets <- contam[contam$celltype_final %in% c("Macrophage1", "Macrophage2") & contam$ratio <= 2,] |> dplyr::pull(target) |> unique()

idx <- meta |> dplyr::filter(celltype_final %in% c("Macrophage1", "Macrophage2")) |> dplyr::arrange(patient) |> rownames()
cts <- cts[good_targets, idx]
meta <- meta[idx,]

mm <- model.matrix(~HOXB13_group, data = meta)
eff <- meta$total_counts
ids <- meta$patient

nfit <- nebula::nebula(count = cts, id = ids, pred = mm, offset = eff, method = "LN", covariance = T, ncore = 4)
saveRDS(nfit, "14_de_analysis_final_outs/macrophage_de_nebula_model.RDS")

nfit <- readRDS("14_de_analysis_final_outs/macrophage_de_nebula_model.RDS")

out <- local_wald_test(nfit = nfit, .contr = contr)
out <- out[out$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
out$p.adj <- p.adjust(out$p.value, method = "BH") # FDR

plot(x = out$log2FC, 
     y = -log10(out$p.adj), 
     pch = 16, cex = 0.5, frame = F, 
     col = ifelse(test = (out$p.adj < 0.05), yes = "red", no = "black"), 
     xlab = "log2FC", ylab = "-log10(adj. P-value)")

library(ggplot2)
out$enrichment <- dplyr::case_when(out$log2FC > 0 & out$p.adj < 0.05 ~ "HOXB13+", 
                                   out$log2FC < 0 & out$p.adj < 0.05 ~ "HOXB13-", 
                                   T ~ NA)
ggplot() + 
  geom_point(data = out[out$log2FC > 0 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), fill = enrichment), shape = 21, size = 2) +
  geom_point(data = out[out$log2FC < 0 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), fill = enrichment), shape = 21, size = 2) +
  scale_fill_manual(values = c("HOXB13+" = "lightgreen", "HOXB13-" = "darkblue")) + 
  geom_point(data = out[out$p.adj >= 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj)), color = "grey", shape = 16, size = 1) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) + 
  ggrepel::geom_text_repel(data = out[abs(out$log2FC) > 0 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), label = target), size = 3) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5), 
        legend.position = "top"
  ) + 
  labs(title = "Macrophages", fill = "Enrichment") + 
  guides(color = guide_legend(override.aes = list(size = 4)))
# ggsave(filename = "14_de_analysis_final_outs/macrophage_de_nebula_volcano_plot.pdf", width = 8, height = 8)
# openxlsx::write.xlsx(x = out |> dplyr::arrange(enrichment, desc(wald.stat)), file = "14_de_analysis_final_outs/macrophage_de_nebula_results_file.xlsx")

mean(out$log2FC > 0) # 0.5183099


# CD8 T cells ------------------------------------------------------------------
source(file = url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/refs/heads/Main/_code/smiDE/R/pre_de_contamination_screen.R"))
coords <- split(x = adata$obs[,c("x_centroid", "y_centroid")], f = adata$obs$core_global)
ids <- split(x = rownames(adata$obs), f = adata$obs$core_global)

library(Matrix)
neighbs <- lapply(X = coords, FUN = spatula::getSpatialNeighbors, return_weights = F, dist_thresh_quantile = 0.98)
neighbs <- Matrix::bdiag(neighbs)
dimnames(neighbs) <- list(unlist(ids), unlist(ids))

cts <- adata$layers$counts |> t() |> as("CsparseMatrix") |> magrittr::set_colnames(adata$obs_names) |> magrittr::set_rownames(adata$var_names)
meta <- adata$obs
meta$cell_ID <- rownames(meta)
meta$celltype_simple <- dplyr::case_when(grepl(pattern = "CD8", x = meta$celltype_final) ~ "T.CD8", 
                                         grepl(pattern = "CD4", x = meta$celltype_final) ~ "T.CD4", 
                                         grepl(pattern = "Macr", x = meta$celltype_final) ~ "Myeloid",
                                         T ~ meta$celltype_final
                                         )

library(data.table)
contam <- overlap_ratio_metric(assay_matrix = cts, 
                               metadata = meta, 
                               adjacency_matrix = neighbs, 
                               cluster_col = "celltype_simple", 
                               cellid_col = "cell_ID", 
                               verbose = T)
good_targets <- contam[contam$celltype_simple == "T.CD8" & contam$ratio <= 1.5,] |> dplyr::pull(target) |> unique()

idx <- meta |> dplyr::filter(celltype_simple == "T.CD8") |> dplyr::arrange(patient) |> rownames()
cts <- cts[good_targets, idx]
meta <- meta[idx,]

mm <- model.matrix(~HOXB13_group, data = meta)
eff <- meta$total_counts
ids <- meta$patient

nfit <- nebula::nebula(count = cts, id = ids, pred = mm, offset = eff, method = "LN", covariance = T, ncore = 4)
saveRDS(nfit, "14_de_analysis_final_outs/cd8t_de_nebula_model.RDS")

nfit <- readRDS("14_de_analysis_final_outs/cd8t_de_nebula_model.RDS")

out <- local_wald_test(nfit = nfit, .contr = contr)
out <- out[out$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
out$p.adj <- p.adjust(out$p.value, method = "BH") # FDR

plot(x = out$log2FC, 
     y = -log10(out$p.adj), 
     pch = 16, cex = 0.5, frame = F, 
     col = ifelse(test = (out$p.adj < 0.05), yes = "red", no = "black"), 
     xlab = "log2FC", ylab = "-log10(adj. P-value)")

library(ggplot2)
out$enrichment <- dplyr::case_when(out$log2FC > 0 & out$p.adj < 0.05 ~ "HOXB13+", 
                                   out$log2FC < 0 & out$p.adj < 0.05 ~ "HOXB13-", 
                                   T ~ NA)
ggplot() + 
  geom_point(data = out[out$log2FC > 0 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), fill = enrichment), shape = 21, size = 2) +
  geom_point(data = out[out$log2FC < 0 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), fill = enrichment), shape = 21, size = 2) +
  scale_fill_manual(values = c("HOXB13+" = "lightgreen", "HOXB13-" = "darkblue")) + 
  geom_point(data = out[out$p.adj >= 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj)), color = "grey", shape = 16, size = 1) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) + 
  ggrepel::geom_text_repel(data = out[abs(out$log2FC) > 0 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), label = target), size = 3) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5), 
        legend.position = "top"
  ) + 
  labs(title = "CD8+ T Cells", fill = "Enrichment") + 
  guides(color = guide_legend(override.aes = list(size = 4)))
# ggsave(filename = "14_de_analysis_final_outs/cd8t_de_nebula_volcano_plot.pdf", width = 8, height = 8)
# openxlsx::write.xlsx(x = out |> dplyr::arrange(enrichment, desc(wald.stat)), file = "14_de_analysis_final_outs/cd8t_de_nebula_results_file.xlsx")

mean(out$log2FC > 0) # 0.4690265


## TBD 
# - Differentially expressed genes in CD8 T cells adjacent to IDO1 hi vs lo tumor cells?
# - Differentially expressed genes in Macrophages adjacent to IDO1 hi vs lo tumor cells?

