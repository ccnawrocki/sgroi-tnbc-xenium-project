rm(list = ls())
.rs.restartR(clean = T)

#dir.create("25_da_analysis_round3_outs")
adata <- anndataR::read_h5ad(path = "sgroi-tnbc_filtered.h5ad")
counts <- table(adata$obs$patient, adata$obs$celltype_final_amended_2) 
props <- counts |> 
  prop.table(margin = 1) |> 
  as.data.frame.matrix()
ptgroups <- (table(adata$obs$patient, adata$obs$HOXB13_group) > 0)[,"HOXB13+"] |> ifelse(test = _, yes = "HOXB13+", no = "HOXB13-")

all(rownames(props) == names(ptgroups))

# Simplest approach
wtests <- presto::wilcoxauc(X = t(props), y = ptgroups)

# Probably the more principled approach
modtest <- function(ct) {
  dtmp <- data.frame("prop" = props[,ct] |> as.numeric(), "count" = counts[,ct] |> as.numeric(), "total" = rowSums(counts),
                     "group" = ptgroups |> as.factor())
  #mod <- glm(formula = cbind(count, total - count) ~ group, data = dtmp, family = quasibinomial(link = "logit"))
  mod <- MASS::glm.nb(formula = count ~ group + offset(log(total)), data = dtmp)
  return(summary(mod)[["coefficients"]][2,,drop=F])
}

ctoi <- levels(adata$obs$celltype_final_amended_2) |> setNames(levels(adata$obs$celltype_final_amended_2))
res <- lapply(X = ctoi, FUN = modtest)
res <- lapply(X = res, FUN = as.data.frame)
res <- dplyr::bind_rows(res, .id = "celltype")
res$contrast <- "HOXB13+ minus HOXB13-"
colnames(res) <- c("celltype", "logOR", "logOR_se", "zstat", "pval", "contrast")
rownames(res) <- NULL
res$log2OR <- res$logOR / log(2)
res$log2OR_se <- res$logOR_se / log(2)
res$padj <- p.adjust(p = res$pval, method = "BH")

# edgeR might even be better... 
library(edgeR)
counts_mat <- table(adata$obs$patient, adata$obs$celltype_final_amended_2) |> 
  as.matrix() |> t()  # cell types as rows, patients as columns

group <- ptgroups[colnames(counts_mat)] |> as.factor()
mm <- model.matrix(~group)

y <- DGEList(counts = counts_mat, 
             lib.size = colSums(counts_mat),
             group = group)

# Estimate dispersion with empirical Bayes shrinkage across cell types
y <- estimateDisp(y, mm)

# Fit and test
fit <- glmQLFit(y, mm)
fitres <- glmQLFTest(fit, coef = "groupHOXB13+")
fitres <- topTags(fitres, n = Inf) |> as.data.frame()

#pdf(file = "25_da_analysis_round3_outs/boxplots_with_points_all_celltypes_proportion_of_patient.pdf", width = 6, height = 8)
for (ct in levels(adata$obs$celltype_final_amended_2)) {
  boxplot(props[,ct] ~ ptgroups, 
          col = c("darkblue", "lightgreen"),
          outline = F, 
          main = ct, 
          frame = F, 
          boxwex = 0.5, 
          xlab = NA, ylab = "Proportion of Patient's Cells")
  stripchart(props[,ct] ~ ptgroups,
             method = "jitter", jitter = 0.15,
             add = T, vertical = T, 
             pch = 21, cex = 1, bg = "grey")
}
#dev.off()

# I will stick with my negative binomial modeling from above, since I do not 
# completely understand the quasi framework in edgeR right now.
library(ggplot2)
res <- dplyr::arrange(res, desc(log2OR)) |> dplyr::mutate(celltype = factor(celltype, levels = celltype))
res$enrichment <- dplyr::case_when(res$log2OR > 0 & res$pval < 0.05 ~ "HOXB13+", 
                                   res$log2OR < 0 & res$pval < 0.05 ~ "HOXB13-", 
                                   T ~ "neither")
ggplot() + 
  geom_errorbar(data = res, aes(xmin = log2OR-1.96*log2OR_se, xmax = log2OR+1.96*log2OR_se, y = celltype), color = "black", width = 0.5) +
  geom_point(data = res, aes(x = log2OR, y = celltype, fill = enrichment), shape = 21, color = "black", size = 4) + 
  geom_text(data = res, aes(x = (log2OR+1.96*log2OR_se)+0.25, y = celltype, label = paste0("p = ", round(pval, digits = 3), "\n", "FDR = ", round(padj, digits = 3))), color = "black", size = 2.5, hjust = 0, fontface = "bold") + 
  scale_fill_manual(values = c("HOXB13-" = "darkblue","HOXB13+" = "lightgreen", "neither" = "black")) + 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  ggthemes::theme_hc() + 
  theme(panel.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5), 
        axis.title.y = element_blank(), axis.text = element_text(color = "black", size = 12), 
        axis.ticks.length = unit(2.5, "mm"), axis.line = element_line(color = "black")) + 
  labs(title = "Differential Abundance\nTesting Summary") + 
  scale_x_continuous(limits = c(-4, 5), breaks = seq(-4, 5, 1))
#ggsave(filename = "25_da_analysis_round3_outs/forest_plot_all_celltypes_proportion_of_patient.pdf", height = 8.5, width = 8.5)

