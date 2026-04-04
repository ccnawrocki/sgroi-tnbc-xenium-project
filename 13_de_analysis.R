rm(list = ls())
.rs.restartR(clean = T)

# This is a preliminary round of DE, not finalized! 

# renv::install("nebula")
# renv::install("ggthemes")
# dir.create("13_de_analysis_outs")
adata <- anndataR::read_h5ad(path = "sgroi-tnbc_filtered.h5ad", mode = "r+")

### Adding the HOXB13 groupings
sample_data <- openxlsx::read.xlsx("../../sgroi-tnbc-data/For Cole TNBC 2021-2025 cohort clinpath and HOX data.xlsx")
sample_data$`TNBC.#` <- paste0("p", sprintf(sample_data$`TNBC.#`, fmt = "%02d"))
sample_data$HOXB13_group <- ifelse(test = grepl(pattern = "Pos", x = sample_data$HOXB13), yes = "positive", no = "negative")
sample_data <- dplyr::select(sample_data, c(1, 3, 40))
adata$obs$HOXB13_group <- plyr::mapvalues(x = adata$obs$patient, from = sample_data$`TNBC.#`, to = sample_data$HOXB13_group)


### DE 
# Epithelial cells -------------------------------------------------------------
library(Matrix)
idx <- adata$obs |> dplyr::filter(celltype_level1 == "Epithelial") |> dplyr::arrange(patient) |> rownames()
cts <- adata$layers$counts |> t() |> as("CsparseMatrix") |> magrittr::set_colnames(adata$obs_names) |> magrittr::set_rownames(adata$var_names)
cts <- cts[,idx]
meta <- adata$obs[idx,]

# mm <- model.matrix(~HOXB13_group, data = meta)
# eff <- meta$total_counts
# ids <- meta$patient
# 
# nfit <- nebula::nebula(count = cts, id = ids, pred = mm, offset = eff, method = "LN", covariance = T, ncore = 4)
# saveRDS(nfit, "13_de_analysis_outs/epithelial_de_nebula_model.RDS")

nfit <- readRDS("13_de_analysis_outs/epithelial_de_nebula_model.RDS")

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
out$enrichment <- dplyr::case_when(out$log2FC > 1 & out$p.adj < 0.05 ~ "positive", 
                                   out$log2FC < -1 & out$p.adj < 0.05 ~ "negative", 
                                   T ~ NA)
ggplot() + 
  geom_point(data = out[out$log2FC > 1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), color = enrichment), shape = 16, size = 2) +
  geom_point(data = out[out$log2FC < -1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), color = enrichment), shape = 16, size = 2) +
  scale_color_manual(values = c("positive" = "steelblue1", "negative" = "darkblue")) + 
  geom_point(data = out, mapping = aes(x = log2FC, y = -log10(p.adj)), color = "grey", shape = 16, size = 0.5) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) + 
  ggrepel::geom_text_repel(data = out[abs(out$log2FC) > 1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), label = target)) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5), 
        legend.position = "top"
  ) + 
  labs(title = "Epithelial Cells", color = "Enrichment") + 
  guides(color = guide_legend(override.aes = list(size = 4)))

mean(out$log2FC > 0) # 0.6092437

library(data.table)
psb <- presto::collapse_counts(counts_mat = cts, meta_data = meta, 
                               varnames = c("core_global", "patient", "HOXB13_group"), 
                               get_norm = F, how = "sum", min_cells_per_group = 20) # Only 20+ cell samples!

psb$meta$patient |> dplyr::n_distinct() # 26
adata$obs$patient |> dplyr::n_distinct() # 26

mm <- model.matrix(~HOXB13_group, data = psb$meta)
cont <- c(0, -1)

psb$meta$N |> barplot(ylab = "N cells", main = "Pseudo-bulk Samples")
psb$counts_mat |> colSums() |> barplot(ylab = "Library Size", main = "Pseudo-bulk Samples", xaxt = "n")

library(edgeR)
y <- psb$counts_mat |> DGEList() |> calcNormFactors(method = "TMM")
v <- voomLmFit(counts = y,
               design = mm, 
               sample.weights = T, 
               plot = T, 
               save.plot = T,
               normalize.method = "none",
               keep.EList = T, 
               block = psb$meta_data$patient
)
y <- y[v$voom.xy$x > 5.5, , keep.lib.sizes = T]
v <- voomLmFit(counts = y,
               design = mm, 
               sample.weights = T, 
               plot = T, 
               save.plot = T,
               normalize.method = "none",
               keep.EList = T
)
out <- v |> contrasts.fit(contrasts = cont) |> eBayes() |> topTable(number = Inf)
plot(x = out$logFC, 
     y = -log10(out$adj.P.Val), 
     pch = 16, cex = 0.5, frame = F, 
     col = ifelse(test = (out$adj.P.Val < 0.01 & abs(out$logFC) > 1), yes = "red", no = "black"), 
     xlab = "log2FC", ylab = "-log10(adj. P-value)")

out$enrichment <- dplyr::case_when(out$logFC > 1 & out$adj.P.Val < 0.05 ~ "positive", 
                                   out$logFC < -1 & out$adj.P.Val < 0.05 ~ "negative", 
                                   T ~ NA)
out$target <- rownames(out)
ggplot() + 
  geom_point(data = out[out$logFC > 1 & out$adj.P.Val < 0.05,], mapping = aes(x = logFC, y = -log10(adj.P.Val), color = enrichment), shape = 16, size = 2) +
  geom_point(data = out[out$logFC < -1 & out$adj.P.Val < 0.05,], mapping = aes(x = logFC, y = -log10(adj.P.Val), color = enrichment), shape = 16, size = 2) +
  scale_color_manual(values = c("positive" = "steelblue1", "negative" = "darkblue")) + 
  geom_point(data = out, mapping = aes(x = logFC, y = -log10(adj.P.Val)), color = "grey", shape = 16, size = 0.5) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) + 
  ggrepel::geom_text_repel(data = out[abs(out$logFC) > 1 & out$adj.P.Val < 0.05,], mapping = aes(x = logFC, y = -log10(adj.P.Val), label = target)) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5), 
        legend.position = "top"
  ) + 
  labs(title = "Epithelial Cells", color = "Enrichment") + 
  guides(color = guide_legend(override.aes = list(size = 4)))

mean(out$logFC > 0) # 0.5021008

# Fibroblasts ------------------------------------------------------------------
idx <- adata$obs |> dplyr::filter(celltype_level1 == "Stromal") |> dplyr::arrange(patient) |> rownames()
cts <- adata$layers$counts |> t() |> as("CsparseMatrix") |> magrittr::set_colnames(adata$obs_names) |> magrittr::set_rownames(adata$var_names)
cts <- cts[,idx]
meta <- adata$obs[idx,]

# mm <- model.matrix(~HOXB13_group, data = meta)
# eff <- meta$total_counts
# ids <- meta$patient
# 
# nfit <- nebula::nebula(count = cts, id = ids, pred = mm, offset = eff, method = "LN", covariance = T, ncore = 4)
# saveRDS(nfit, "13_de_analysis_outs/stromal_de_nebula_model.RDS")

nfit <- readRDS("13_de_analysis_outs/stromal_de_nebula_model.RDS")

out <- local_wald_test(nfit = nfit, .contr = contr)
out <- out[out$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
out$p.adj <- p.adjust(out$p.value, method = "BH") # FDR

plot(x = out$log2FC, 
     y = -log10(out$p.adj), 
     pch = 16, cex = 0.5, frame = F, 
     col = ifelse(test = (out$p.adj < 0.01 & abs(out$log2FC) > 1), yes = "red", no = "black"), 
     xlab = "log2FC", ylab = "-log10(adj. P-value)")

out$enrichment <- dplyr::case_when(out$log2FC > 1 & out$p.adj < 0.05 ~ "positive", 
                                   out$log2FC < -1 & out$p.adj < 0.05 ~ "negative", 
                                   T ~ NA)
ggplot() + 
  geom_point(data = out[out$log2FC > 1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), color = enrichment), shape = 16, size = 2) +
  geom_point(data = out[out$log2FC < -1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), color = enrichment), shape = 16, size = 2) +
  scale_color_manual(values = c("positive" = "steelblue1", "negative" = "darkblue")) + 
  geom_point(data = out, mapping = aes(x = log2FC, y = -log10(p.adj)), color = "grey", shape = 16, size = 0.5) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) + 
  ggrepel::geom_text_repel(data = out[abs(out$log2FC) > 1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), label = target)) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5), 
        legend.position = "top"
  ) + 
  labs(title = "Stromal", color = "Enrichment") + 
  guides(color = guide_legend(override.aes = list(size = 4)))

mean(out$logFC > 0) # 0.6892779

# Endothelial ------------------------------------------------------------------
idx <- adata$obs |> dplyr::filter(celltype_level1 == "Endothelial") |> dplyr::arrange(patient) |> rownames()
cts <- adata$layers$counts |> t() |> as("CsparseMatrix") |> magrittr::set_colnames(adata$obs_names) |> magrittr::set_rownames(adata$var_names)
cts <- cts[,idx]
meta <- adata$obs[idx,]

# mm <- model.matrix(~HOXB13_group, data = meta)
# eff <- meta$total_counts
# ids <- meta$patient
# 
# nfit <- nebula::nebula(count = cts, id = ids, pred = mm, offset = eff, method = "LN", covariance = T, ncore = 4)
# saveRDS(nfit, "13_de_analysis_outs/endothelial_de_nebula_model.RDS")

nfit <- readRDS("13_de_analysis_outs/endothelial_de_nebula_model.RDS")

out <- local_wald_test(nfit = nfit, .contr = contr)
out <- out[out$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
out$p.adj <- p.adjust(out$p.value, method = "BH") # FDR

plot(x = out$log2FC, 
     y = -log10(out$p.adj), 
     pch = 16, cex = 0.5, frame = F, 
     col = ifelse(test = (out$p.adj < 0.01 & abs(out$log2FC) > 1), yes = "red", no = "black"), 
     xlab = "log2FC", ylab = "-log10(adj. P-value)")

out$enrichment <- dplyr::case_when(out$log2FC > 1 & out$p.adj < 0.05 ~ "positive", 
                                   out$log2FC < -1 & out$p.adj < 0.05 ~ "negative", 
                                   T ~ NA)
ggplot() + 
  geom_point(data = out[out$log2FC > 1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), color = enrichment), shape = 16, size = 2) +
  geom_point(data = out[out$log2FC < -1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), color = enrichment), shape = 16, size = 2) +
  scale_color_manual(values = c("positive" = "steelblue1", "negative" = "darkblue")) + 
  geom_point(data = out, mapping = aes(x = log2FC, y = -log10(p.adj)), color = "grey", shape = 16, size = 0.5) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) + 
  ggrepel::geom_text_repel(data = out[abs(out$log2FC) > 1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), label = target)) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5), 
        legend.position = "top"
  ) + 
  labs(title = "Endothelial Cells", color = "Enrichment") + 
  guides(color = guide_legend(override.aes = list(size = 4)))

# Lymphoid ---------------------------------------------------------------------
idx <- adata$obs |> dplyr::filter(celltype_level1 == "Lymphoid") |> dplyr::arrange(patient) |> rownames()
cts <- adata$layers$counts |> t() |> as("CsparseMatrix") |> magrittr::set_colnames(adata$obs_names) |> magrittr::set_rownames(adata$var_names)
cts <- cts[,idx]
meta <- adata$obs[idx,]

# mm <- model.matrix(~HOXB13_group, data = meta)
# eff <- meta$total_counts
# ids <- meta$patient
# 
# nfit <- nebula::nebula(count = cts, id = ids, pred = mm, offset = eff, method = "LN", covariance = T, ncore = 4)
# saveRDS(nfit, "13_de_analysis_outs/lymphoid_de_nebula_model.RDS")

nfit <- readRDS("13_de_analysis_outs/lymphoid_de_nebula_model.RDS")

out <- local_wald_test(nfit = nfit, .contr = contr)
out <- out[out$convergence %in% c(1, -10),] # To be safe, do not trust genes whose model did not converge
out$p.adj <- p.adjust(out$p.value, method = "BH") # FDR

plot(x = out$log2FC, 
     y = -log10(out$p.adj), 
     pch = 16, cex = 0.5, frame = F, 
     col = ifelse(test = (out$p.adj < 0.01 & abs(out$log2FC) > 1), yes = "red", no = "black"), 
     xlab = "log2FC", ylab = "-log10(adj. P-value)")

out$enrichment <- dplyr::case_when(out$log2FC > 1 & out$p.adj < 0.05 ~ "positive", 
                                   out$log2FC < -1 & out$p.adj < 0.05 ~ "negative", 
                                   T ~ NA)
ggplot() + 
  geom_point(data = out[out$log2FC > 1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), color = enrichment), shape = 16, size = 2) +
  geom_point(data = out[out$log2FC < -1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), color = enrichment), shape = 16, size = 2) +
  scale_color_manual(values = c("positive" = "steelblue1", "negative" = "darkblue")) + 
  geom_point(data = out, mapping = aes(x = log2FC, y = -log10(p.adj)), color = "grey", shape = 16, size = 0.5) + 
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) + 
  ggrepel::geom_text_repel(data = out[abs(out$log2FC) > 1 & out$p.adj < 0.05,], mapping = aes(x = log2FC, y = -log10(p.adj), label = target)) + 
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 12, color = "black"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"), 
        plot.title = element_text(hjust = 0.5, size = 16), plot.subtitle = element_text(hjust = 0.5), 
        legend.position = "top"
  ) + 
  labs(title = "T and NK Cells", color = "Enrichment") + 
  guides(color = guide_legend(override.aes = list(size = 4)))

