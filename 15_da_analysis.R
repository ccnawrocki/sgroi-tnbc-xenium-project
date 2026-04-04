rm(list = ls())
.rs.restartR(clean = T)

# dir.create("15_da_analysis_final_outs")
adata <- anndataR::read_h5ad(path = "sgroi-tnbc_filtered.h5ad", mode = "r+")

# Reconciling all cell-typing results ------------------------------------------
rhdf5::h5ls("sgroi-tnbc_filtered_immune_subset.h5ad")
imm <- rhdf5::h5read(file = "sgroi-tnbc_filtered_immune_subset.h5ad", name = "obs")[c("celltype_level2", "_index")]
imm <- data.frame(cell = imm$`_index`, celltype = 
                    plyr::mapvalues(x = imm$celltype_level2$codes, 
                                    from = seq_along(imm$celltype_level2$categories)-1, 
                                    to = imm$celltype_level2$categories)
                  )
tcs <- rhdf5::h5read(file = "sgroi-tnbc_filtered_tcell_subset.h5ad", name = "obs")[c("celltype_level3", "_index")]
tcs <- data.frame(cell = tcs$`_index`, celltype = 
                    plyr::mapvalues(x = tcs$celltype_level3$codes, 
                                    from = seq_along(tcs$celltype_level3$categories)-1, 
                                    to = tcs$celltype_level3$categories)
)

adata$obs$celltype_level2 <- plyr::mapvalues(x = adata$obs_names, from = imm$cell, to = imm$celltype)
adata$obs$celltype_level2 <- ifelse(test = grepl(pattern = "-", x = adata$obs$celltype_level2), yes = NA, no = adata$obs$celltype_level2)
adata$obs$celltype_level3 <- plyr::mapvalues(x = adata$obs_names, from = tcs$cell, to = tcs$celltype)
adata$obs$celltype_level3 <- ifelse(test = grepl(pattern = "-", x = adata$obs$celltype_level3), yes = NA, no = adata$obs$celltype_level3)

adata$obs$celltype_final <- dplyr::case_when(
  !is.na(adata$obs$celltype_level3) ~ adata$obs$celltype_level3, 
  !is.na(adata$obs$celltype_level2) ~ adata$obs$celltype_level2,
  T ~ adata$obs$celltype_level1
)

adata$obs$celltype_final <- ifelse(test = (adata$obs$celltype_final == "T and NK"), yes = "T.low.quality", no = adata$obs$celltype_final)

library(BPCells)
library(ggplot2)
plot_embedding(adata$obs$celltype_final, adata$obsm$X_umap, rasterize = T, labels_discrete = F) + 
  labs(x = "UMAP1", y = "UMAP2", title = "Final Celltype")
# anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered.h5ad", mode = "w")

# renv::install("Nanostring-Biostats/CosMx-Analysis-Scratch-Space/_code/HieraType")

# Using the Hieratype's plotting function
fctbl_unsup <- 
  HieraType::clusterwise_foldchange_metrics(normed  = Matrix::t(adata$layers$lognorm) |> 
                                              magrittr::set_rownames(rownames(adata$var)) |> 
                                              magrittr::set_colnames(rownames(adata$obs)) |> 
                                              as("CsparseMatrix"), 
                                            metadata = adata$obs, 
                                            cluster_column = "celltype_final"
  )
hm_unsup <- HieraType::marker_heatmap(fctbl_unsup, 
                                      extras = c("CD3D", "CD4", "FOXP3", "CD8B", "CD8A")
)
print(hm_unsup) # Looks very reasonable


# Inflammation scores ----------------------------------------------------------
nontumor <- (adata$obs$celltype_final != "Epithelial")
immune <- !(adata$obs$celltype_final %in% c("Epithelial", "Ductal", "Endothelial", "Stromal"))

inflammation <- data.frame(nontumor = nontumor, immune = immune, patient = adata$obs$patient)
scores <- dplyr::group_by(inflammation, patient) |> dplyr::summarise(infl_score = sum(immune)/sum(nontumor))

sample_data <- openxlsx::read.xlsx(xlsxFile = "~/Library/CloudStorage/OneDrive-MassGeneralBrigham/xenium/sgroi-tnbc-data/Normalized HOXB13.xlsx")
sample_data$patient <- gsub(pattern = "([a-z])|([A-Z])|( )|(-)", replacement = "", x = sample_data$Sample_Name) |> sprintf(fmt = "%02s") |> paste0("p", ... = _)
sample_data$group <- ifelse(test = sample_data$H_norm > -13.5, yes = "HOXB13+", no = "HOXB13-")
scores$group <- plyr::mapvalues(x = scores$patient, from = sample_data$patient, to = sample_data$group)
scores <- dplyr::arrange(scores, infl_score)

# pdf("15_da_analysis_final_outs/inflammation_scores_by_patient.pdf", width = 10, height = 6)
mps <- barplot(scores$infl_score, col = ifelse(test = scores$group == "HOXB13+", yes = "lightgreen", no = "darkblue"), 
               ylim = c(0, 1), axes = F,
               main = "Inflammation Score by Patient", ylab = "(N immune cells) / (N non-tumor cells)")
abline(h = 0.5, col = "black", lwd = 1, lty = "dashed")
legend(x = 0.5, y = 0.95, legend = c("+", "-"), fill = c("lightgreen", "darkblue"), pch = NA, title = "HOXB13 group", adj = c(0.5, 0.5))
mtext(text = scores$patient, at = mps, side = 1, cex = 0.7, col = ifelse(test = scores$infl_score > 0.5, yes = "red3", no = "steelblue"))
axis(side = 2, at = seq(0, 1, 0.1), labels = paste0(seq(0, 1, 0.1)*100, "%"), las = 1)
text(x = c(1, 1), y = c(0.45, 0.55), labels = c("uninflamed", "inflamed"), col = c("steelblue", "red3"))
# dev.off()

scores$H_norm <- plyr::mapvalues(x = scores$patient, from = sample_data$patient, to = sample_data$H_norm) |> as.numeric()
# pdf("15_da_analysis_final_outs/inflammation_scores_vs_HOXB13_qPCR.pdf", width = 8, height = 6)
plot(scores$H_norm, scores$infl_score, pch = 16, 
     xlab = "HOXB13 Normalized qPCR", ylab = "Inflammation Score", main = "Relationship Between HOXB13 qPCR and Inflammation Score")
abline(reg = lm(infl_score ~ H_norm, scores), lwd = 1, lty = "dashed")
text(22, 0.95, paste("rho", "=", cor(scores$H_norm, scores$infl_score, method = "spearman") |> round(3)))
# dev.off()

# Add the HOXB13 group and the inflammation group to the anndata
adata$obs$HOXB13_group <- plyr::mapvalues(x = adata$obs$patient, from = scores$patient, to = scores$group |> as.character())
scores$infl_group <- ifelse(test = scores$infl_score > 0.5, yes = "inflamed", no = "uninflamed")
adata$obs$inflammation_group <- plyr::mapvalues(x = adata$obs$patient, from = scores$patient, to = scores$infl_group |> as.character())
# anndataR::write_h5ad(object = adata, path = "sgroi-tnbc_filtered.h5ad", mode = "w")


# CD8 to Treg ratios -----------------------------------------------------------
cd8s <- grepl(pattern = "CD8", x = adata$obs$celltype_final)
tregs <- grepl(pattern = "Treg", x = adata$obs$celltype_final)

ts <- data.frame(cd8s = cd8s, tregs = tregs, patient = adata$obs$patient)
ratios <- dplyr::group_by(ts, patient) |> dplyr::summarise(cd8_to_treg = sum(cd8s)/sum(tregs))
ratios <- dplyr::arrange(ratios, cd8_to_treg)
ratios <- ratios[ratios$patient %in% scores[scores$infl_group == "inflamed",]$patient,]
ratios$group <- plyr::mapvalues(x = ratios$patient, from = sample_data$patient, to = sample_data$group)

# pdf("15_da_analysis_final_outs/Tcd8_to_Treg_ratio_for_inflamed_patients.pdf", width = 8, height = 6)
mps <- barplot(ratios$cd8_to_treg, col = ifelse(test = ratios$group == "HOXB13+", yes = "lightgreen", no = "darkblue"),
               main = "CD8 T cell to Treg Ratio by Patient", ylab = "(N CD8 T cells) / (N Tregs)")
legend(x = 0.5, y = 5, legend = c("+", "-"), fill = c("lightgreen", "darkblue"), pch = NA, title = "HOXB13 group", adj = c(0.5, 0.5))
mtext(text = ratios$patient, at = mps, side = 1, cex = 0.7, col = "red3")
# dev.off()


# CD8 activated and exhausted --------------------------------------------------
# - (T.CD8.activated.classical + T.CD8.activated.proliferative) / total non-tumor cells (for each patient colored by Hoxb13 high/low like your recent graphs)
act <- grepl(pattern = "CD8.activated", x = adata$obs$celltype_final)
nontumor <- (adata$obs$celltype_final != "Epithelial")

df <- data.frame(cd8_act = act, nontumor = nontumor, patient = adata$obs$patient)
ratios <- dplyr::group_by(df, patient) |> dplyr::summarise(cd8_act_to_nontumor = sum(cd8_act)/sum(nontumor))
ratios <- dplyr::arrange(ratios, cd8_act_to_nontumor)
ratios$group <- plyr::mapvalues(x = ratios$patient, from = sample_data$patient, to = sample_data$group)

# pdf("15_da_analysis_final_outs/Tcd8_act_to_nontumor_ratio_by_patient.pdf", width = 8, height = 6)
mps <- barplot(ratios$cd8_act_to_nontumor, col = ifelse(test = ratios$group == "HOXB13+", yes = "lightgreen", no = "darkblue"),
               main = "Activated CD8 T Cell to Non-Tumor Cell Ratio by Patient", ylab = "(N activated CD8 T cells) / (N non-tumor cells)")
legend(x = 0.5, y = 0.2, legend = c("+", "-"), fill = c("lightgreen", "darkblue"), pch = NA, title = "HOXB13 group", adj = c(0.5, 0.5))
mtext(text = ratios$patient, at = mps, side = 1, cex = 0.7, col = "black")
# dev.off()

# - (T.CD8.activated.classical + T.CD8.activated.proliferative) / Tregs specifically (for each patient colored by Hoxb13 high/low like your recent graphs)
act <- grepl(pattern = "CD8.activated", x = adata$obs$celltype_final)
tregs <- (adata$obs$celltype_final == "Treg")

df <- data.frame(cd8_act = act, treg = tregs, patient = adata$obs$patient)
ratios <- dplyr::group_by(df, patient) |> dplyr::summarise(cd8_act_to_treg = sum(cd8_act)/sum(treg))
ratios <- dplyr::arrange(ratios, cd8_act_to_treg)
ratios$group <- plyr::mapvalues(x = ratios$patient, from = sample_data$patient, to = sample_data$group)

# pdf("15_da_analysis_final_outs/Tcd8_act_to_Treg_ratio_by_patient.pdf", width = 8, height = 6)
mps <- barplot(ratios$cd8_act_to_treg, col = ifelse(test = ratios$group == "HOXB13+", yes = "lightgreen", no = "darkblue"),
               main = "Activated CD8 T Cell to Treg Ratio by Patient", ylab = "(N activated CD8 T cells) / (N Tregs)")
legend(x = 0.5, y = 3, legend = c("+", "-"), fill = c("lightgreen", "darkblue"), pch = NA, title = "HOXB13 group", adj = c(0.5, 0.5))
mtext(text = ratios$patient, at = mps, side = 1, cex = 0.7, col = "black")
# dev.off()

# - (T.CD8.activated.classical + T.CD8.activated.proliferative) / T.CD8.exhausted (for each patient colored by Hoxb13 high/low like your recent graphs)
act <- grepl(pattern = "CD8.activated", x = adata$obs$celltype_final)
exh <- (adata$obs$celltype_final == "T.CD8.exhausted")

df <- data.frame(cd8_act = act, cd8_exh = exh, patient = adata$obs$patient)
ratios <- dplyr::group_by(df, patient) |> dplyr::summarise(cd8_act_to_exh = sum(cd8_act)/sum(cd8_exh))
ratios <- dplyr::arrange(ratios, cd8_act_to_exh)
ratios$group <- plyr::mapvalues(x = ratios$patient, from = sample_data$patient, to = sample_data$group)
ratios <- ratios[ratios$cd8_act_to_exh < Inf,]

# pdf("15_da_analysis_final_outs/Tcd8_act_to_Tcd8_exh_ratio_by_patient.pdf", width = 8, height = 6)
mps <- barplot(ratios$cd8_act_to_exh, col = ifelse(test = ratios$group == "HOXB13+", yes = "lightgreen", no = "darkblue"),
               main = "Activated CD8 T Cell to Exhausted CD8 T Cell Ratio by Patient", ylab = "(N activated CD8 T cells) / (N exhausted CD8 T cells)", 
               sub = "Note: p61 has zero exhausted CD8 T cells and is not shown.")
legend(x = 0.5, y = 25, legend = c("+", "-"), fill = c("lightgreen", "darkblue"), pch = NA, title = "HOXB13 group", adj = c(0.5, 0.5))
mtext(text = ratios$patient, at = mps, side = 1, cex = 0.7, col = "black")
# dev.off()

library(ggplot2)
df <- dplyr::group_by(df, patient) |> dplyr::summarise(activated = sum(cd8_act), exhausted = sum(cd8_exh))
df <- tidyr::pivot_longer(data = df, cols = 2:3, names_to = "CD8 T", values_to = "value")
df <- dplyr::group_by(df, patient) |> dplyr::mutate(total = sum(value))
df$prop <- (df$value/df$total)
df$patient <- factor(x = df$patient, levels = df |> dplyr::filter(`CD8 T` == "activated") |> dplyr::arrange(prop) |> dplyr::pull(patient))
df$group <- plyr::mapvalues(x = df$patient, from = sample_data$patient, to = sample_data$group)
ggplot() + 
  geom_bar(data = df, mapping = aes(x = patient, y = value, fill = `CD8 T`), stat = "identity", position = "dodge") + 
  scale_fill_manual(values = c("red3", "pink3")) +
  ggnewscale::new_scale_fill() +
  geom_rect(data = df, mapping = aes(x = patient, width = 1, y = -100, height = 100, fill = group)) +
  scale_fill_manual(values = c("HOXB13-" = "darkblue", "HOXB13+" = "lightgreen")) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(vjust = 1, size = 12, color = ifelse(test = levels(df$patient) %in% sample_data[sample_data$group == "HOXB13-",]$patient, yes = "darkblue", no = "lightgreen"), face = "bold"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"),
        axis.title.x = element_blank(),
        plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5), 
        legend.position = "top"
  ) + 
  labs(y = "count", title = "CD8 T Cells by Patient")
# ggsave(filename = "15_da_analysis_final_outs/Tcd8_act_and_Tcd8_exh_counts_by_patient.pdf", width = 10, height = 8)
ggplot() + 
  geom_bar(data = df, mapping = aes(x = patient, y = prop, fill = `CD8 T`), stat = "identity") + 
  scale_fill_manual(values = c("red3", "pink3")) +
  ggnewscale::new_scale_fill() +
  geom_rect(data = df, mapping = aes(x = patient, width = 1, y = -0.01, height = 0.01, fill = group)) +
  scale_fill_manual(values = c("HOXB13-" = "darkblue", "HOXB13+" = "lightgreen")) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  ggthemes::theme_hc() + 
  theme(axis.text.x = element_text(vjust = 1, size = 12, color = ifelse(test = levels(df$patient) %in% sample_data[sample_data$group == "HOXB13-",]$patient, yes = "darkblue", no = "lightgreen"), face = "bold"), 
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks = element_line(linewidth = 0.25, color = "black"), 
        axis.ticks.length = unit(0.25, units = "cm"),
        axis.title.x = element_blank(),
        plot.subtitle = element_text(hjust = 0.5), plot.title = element_text(hjust = 0.5), 
        legend.position = "top"
        ) + 
  labs(y = "proportion", title = "CD8 T Cells by Patient")
# ggsave(filename = "15_da_analysis_final_outs/Tcd8_act_and_Tcd8_exh_proportions_by_patient.pdf", width = 10, height = 8)

# Plasmablasts -----------------------------------------------------------------
# - Plasmablasts / total non-tumor cells (for each patient colored by Hoxb13 high/low like your recent graphs)
plsm <- (adata$obs$celltype_final == "Plasmablast")
nontumor <- (adata$obs$celltype_final != "Epithelial")

df <- data.frame(plasmablast = plsm, nontumor = nontumor, patient = adata$obs$patient)
ratios <- dplyr::group_by(df, patient) |> dplyr::summarise(plsm_to_nontumor = sum(plasmablast)/sum(nontumor))
ratios <- dplyr::arrange(ratios, plsm_to_nontumor)
ratios$group <- plyr::mapvalues(x = ratios$patient, from = sample_data$patient, to = sample_data$group)

# pdf("15_da_analysis_final_outs/plasmablast_to_nontumor_ratio_by_patient.pdf", width = 8, height = 6)
mps <- barplot(ratios$plsm_to_nontumor, col = ifelse(test = ratios$group == "HOXB13+", yes = "lightgreen", no = "darkblue"),
               main = "Plasmablast to Non-Tumor Cell Ratio by Patient", ylab = "(N plasmablasts) / (N non-tumor cells)")
legend(x = 0.5, y = 0.06, legend = c("+", "-"), fill = c("lightgreen", "darkblue"), pch = NA, title = "HOXB13 group", adj = c(0.5, 0.5))
mtext(text = ratios$patient, at = mps, side = 1, cex = 0.7, col = "black")
# dev.off()

# - (B+Plasma+Plasmablasts) / total non-tumor cells (for each patient colored by Hoxb13 high/low like your recent graphs)
plsm_all <- (adata$obs$celltype_final %in% c("Plasmablast", "Plasma", "B"))
nontumor <- (adata$obs$celltype_final != "Epithelial")

df <- data.frame(lineageB = plsm_all, nontumor = nontumor, patient = adata$obs$patient)
ratios <- dplyr::group_by(df, patient) |> dplyr::summarise(lineageB_to_nontumor = sum(lineageB)/sum(nontumor))
ratios <- dplyr::arrange(ratios, lineageB_to_nontumor)
ratios$group <- plyr::mapvalues(x = ratios$patient, from = sample_data$patient, to = sample_data$group)

# pdf("15_da_analysis_final_outs/B_cell_lineage_to_nontumor_ratio_by_patient.pdf", width = 8, height = 6)
mps <- barplot(ratios$lineageB_to_nontumor, col = ifelse(test = ratios$group == "HOXB13+", yes = "lightgreen", no = "darkblue"),
               main = "B Lineage Cell to Non-Tumor Cell Ratio by Patient", ylab = "(N B lineage cells) / (N non-tumor cells)")
legend(x = 0.5, y = 0.2, legend = c("+", "-"), fill = c("lightgreen", "darkblue"), pch = NA, title = "HOXB13 group", adj = c(0.5, 0.5))
mtext(text = ratios$patient, at = mps, side = 1, cex = 0.7, col = "black")
# dev.off()

