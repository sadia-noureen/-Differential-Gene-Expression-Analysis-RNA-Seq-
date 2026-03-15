# ============================================================
# Script 04: Visualization
# Project: Differential Gene Expression Analysis (RNA-Seq)
# Plots: PCA, MA plot, Volcano plot, Heatmap, Sample Distance
# ============================================================

library(DESeq2)
library(tidyverse)
## installing 
install.packages("EnhancedVolcano")
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(ggrepel)
library(airway)
intsall.packages("airway")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Load saved objects from script 03
load("data/processed/deseq2_objects.RData")

# ============================================================
# PLOT 1: PCA Plot (Quality Control)
# ============================================================

# Variance-stabilizing transformation before PCA
# Why VST? Raw counts are heteroskedastic — variance increases with mean.
# VST stabilizes variance across the expression range (better than log2+1)
# Reference: Anders & Huber. Genome Biology. 2010;11(10):R106
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2,
                                 color = dex, shape = cell)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(aes(label = name), size = 3, max.overlaps = 10) +
  scale_color_manual(values = c("untrt" = "#2166AC", "trt" = "#D6604D"),
                     labels = c("Untreated", "Treated")) +
  labs(
    title = "Principal Component Analysis",
    subtitle = "Dexamethasone-treated vs Untreated (airway dataset)",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Treatment",
    shape = "Cell Line"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("results/figures/01_PCA_plot.png", pca_plot,
       width = 8, height = 6, dpi = 300, bg = "white")
cat(" PCA plot saved\n")

# ============================================================
# PLOT 2: Sample Distance Heatmap (Quality Control)
# ============================================================

# Euclidean distance between samples — confirms clustering by treatment
samp_dists <- dist(t(assay(vsd)))
samp_dist_matrix <- as.matrix(samp_dists)
rownames(samp_dist_matrix) <- paste(vsd$dex, vsd$cell, sep = "-")
colnames(samp_dist_matrix) <- NULL

ann_colors <- list(
  Treatment = c(untrt = "#2166AC", trt = "#D6604D")
)

annotation_row <- data.frame(
  Treatment = vsd$dex,
  row.names = rownames(samp_dist_matrix)
)

png("results/figures/02_sample_distance_heatmap.png",
    width = 800, height = 700, res = 150)
pheatmap(samp_dist_matrix,
         clustering_distance_rows = samp_dists,
         clustering_distance_cols = samp_dists,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-Sample Distance Matrix",
         fontsize = 11)
dev.off()
cat(" Sample distance heatmap saved\n")

# ============================================================
# PLOT 3: MA Plot
# ============================================================

# MA plot: log fold change (M) vs mean expression (A)
# Shows relationship between expression level and fold change magnitude
# Points should scatter symmetrically around y=0 after shrinkage

png("results/figures/03_MA_plot.png", width = 900, height = 650, res = 150)
plotMA(res_shrunk,
       ylim = c(-4, 4),
       main = "MA Plot: Dexamethasone Treatment vs Untreated\n(LFC shrinkage with apeglm)",
       colSig = "#D6604D",
       colNonSig = "grey60",
       alpha = 0.05)
abline(h = c(-1, 1), col = "dodgerblue", lty = 2, lwd = 1.5)
legend("topright",
       legend = c("Significant DEG", "Non-significant", "|LFC| = 1 threshold"),
       col = c("#D6604D", "grey60", "dodgerblue"),
       pch = c(16, 16, NA), lty = c(NA, NA, 2))
dev.off()
cat(" MA plot saved\n")

# ============================================================
# PLOT 4: Enhanced Volcano Plot (Publication Ready)
# ============================================================

# Filter rows with complete annotation for labeling
res_for_volcano <- res_df %>%
  filter(!is.na(gene_symbol), !is.na(padj), !is.na(log2FoldChange))

volcano <- EnhancedVolcano(
  res_for_volcano,
  lab = res_for_volcano$gene_symbol,
  x = "log2FoldChange",
  y = "padj",
  
  # Thresholds
  pCutoff = 0.05,
  FCcutoff = 1.0,
  
  # Aesthetics
  title = "Differential Gene Expression",
  subtitle = "Dexamethasone Treatment vs Untreated (DESeq2)",
  caption = paste0("Total DEGs: ", nrow(sig_degs),
                   " | Up: ", nrow(up_reg),
                   " | Down: ", nrow(down_reg)),
  
  # Colors: non-sig, sig FC only, sig p only, sig both
  col = c("grey70", "#4DAC26", "#0571B0", "#CA0020"),
  colAlpha = 0.75,
  
  # Label top genes
  selectLab = head(res_for_volcano$gene_symbol[
    res_for_volcano$padj < 0.01 & abs(res_for_volcano$log2FoldChange) > 2], 20),
  
  labSize = 3.5,
  labCol = "black",
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  
  # Axes
  xlim = c(-6, 6),
  legendPosition = "right",
  legendLabSize = 11,
  axisLabSize = 12,
  
  # Gridlines
  gridlines.major = TRUE,
  gridlines.minor = FALSE
)

ggsave("results/figures/04_volcano_plot.png", volcano,
       width = 10, height = 8, dpi = 300, bg = "white")
cat(" Volcano plot saved\n")

# ============================================================
# PLOT 5: Top DEG Heatmap
# ============================================================

# Select top 40 most significant DEGs
top_genes <- sig_degs %>%
  filter(!is.na(gene_symbol)) %>%
  slice_min(padj, n = 40) %>%
  pull(ensembl_id)

# Extract normalized counts for top genes
heat_counts <- assay(vsd)[top_genes, ]

# Replace Ensembl IDs with gene symbols for readable labels
rownames(heat_counts) <- sig_degs$gene_symbol[
  match(rownames(heat_counts), sig_degs$ensembl_id)
]
rownames(heat_counts)[is.na(rownames(heat_counts))] <- top_genes[
  is.na(rownames(heat_counts))
]

# Column annotation
col_ann <- data.frame(
  Treatment = colData(dds)$dex,
  CellLine = colData(dds)$cell,
  row.names = colnames(heat_counts)
)

ann_colors <- list(
  Treatment = c(untrt = "#2166AC", trt = "#D6604D"),
  CellLine = c(N052611 = "#7FC97F", N061011 = "#BEAED4",
               N080611 = "#FDC086", N61311 = "#FFFF99")
)

png("results/figures/05_top_DEG_heatmap.png",
    width = 900, height = 1100, res = 150)
pheatmap(
  heat_counts,
  annotation_col = col_ann,
  annotation_colors = ann_colors,
  scale = "row",              # Z-score per gene (row scaling)
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 9,
  main = "Top 40 Differentially Expressed Genes\n(Z-score normalized, DESeq2)",
  border_color = NA,
  treeheight_row = 30,
  treeheight_col = 20
)
dev.off()
cat(" Top DEG heatmap saved\n")

# ============================================================
# PLOT 6: Count Plots for Top 6 Individual Genes
# ============================================================

# Visualize normalized counts for top 6 genes — shows the actual signal
top6_genes <- sig_degs %>%
  filter(!is.na(gene_symbol)) %>%
  slice_min(padj, n = 6) %>%
  pull(ensembl_id)

count_plots <- lapply(seq_along(top6_genes), function(i) {
  gene_id <- top6_genes[i]
  gene_name <- sig_degs$gene_symbol[sig_degs$ensembl_id == gene_id][1]
  
  d <- plotCounts(dds, gene = gene_id,
                  intgroup = "dex", returnData = TRUE)
  
  ggplot(d, aes(x = dex, y = count, color = dex, fill = dex)) +
    geom_jitter(width = 0.15, size = 3, alpha = 0.8) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    scale_color_manual(values = c("untrt" = "#2166AC", "trt" = "#D6604D")) +
    scale_fill_manual(values = c("untrt" = "#2166AC", "trt" = "#D6604D")) +
    scale_y_log10() +
    labs(title = gene_name,
         x = NULL, y = "Normalized Count (log10)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5))
})

combined_counts <- wrap_plots(count_plots, ncol = 3) +
  plot_annotation(
    title = "Top 6 DEGs: Normalized Count Distribution",
    subtitle = "Treated vs Untreated (DESeq2 size-factor normalized)",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("results/figures/06_top_gene_counts.png", combined_counts,
       width = 12, height = 8, dpi = 300, bg = "white")
cat(" Top gene count plots saved\n")

cat("\n All visualizations complete.\n")
cat(" Proceed to: 05_pathway_enrichment.R\n")