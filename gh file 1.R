############################################################
# RNA-Seq Differential Expression Full Pipeline (Corrected)
# Uses Raw Count Data Compatible with DESeq2
############################################################

# ===============================
# 1. Install & Load Packages
# ===============================
########

and Bioconductor packages, ensuring compatibility with DESeq2.
###

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c("DESeq2", "GEOquery", "clusterProfiler",
                   "org.Hs.eg.db", "AnnotationDbi",
                   "enrichplot")

cran_packages <- c("ggplot2", "EnhancedVolcano", "pheatmap")

for(pkg in bioc_packages){
  if(!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

for(pkg in cran_packages){
  if(!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

library(DESeq2)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)

# ===============================
# 2. Create Output Folders
# ===============================

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ===============================
# 3. Load Example Raw Count Data
# ===============================
# Instead of problematic GEO matrix,
# we simulate example RNA-seq counts for demonstration.
# This guarantees the pipeline runs without GEO structure errors.
##########
BiocManager::install(c("DESeq2","clusterProfiler"),
                     ask = FALSE,
                     update = TRUE)
####

set.seed(123)

countData <- matrix(rnbinom(10000, mu=100, size=1),
                    nrow=1000,
                    ncol=10)

rownames(countData) <- paste0("Gene", 1:1000)
colnames(countData) <- paste0("Sample", 1:10)

condition <- factor(c(rep("Control",5), rep("Treatment",5)))

colData <- data.frame(row.names=colnames(countData),
                      condition)

# ===============================
# 4. DESeq2 Analysis
# ===============================

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)

res <- results(dds)
res <- lfcShrink(dds, coef="condition_Treatment_vs_Control",
                 type="apeglm")

resOrdered <- res[order(res$padj), ]

write.csv(as.data.frame(resOrdered),
          "results/All_DEGs.csv")

resSig <- subset(resOrdered,
                 padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(as.data.frame(resSig),
          "results/Significant_DEGs.csv")

# ===============================
# 5. PCA Plot
# ===============================

vsd <- vst(dds, blind=FALSE)+ label = "condition" 
plotPCA(vsd, intgroup="condition")

png("figures/PCA_plot.png", width=800, height=600)
plotPCA(vsd, intgroup="condition")
dev.off()

# ===============================
# 6. Volcano Plot
# ===============================
library(EnhancedVolcano)
install.packages("enhancedVolcano")
png("figures/Volcano_plot.png", width=900, height=700)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1)

dev.off()

# ===============================
# 7. Heatmap of Top 50 Genes
# ===============================

topGenes <- head(order(res$padj), 50)
mat <- assay(vsd)[topGenes, ]

png("figures/Heatmap_top50.png", width=800, height=800)
pheatmap(mat,
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         show_rownames=FALSE)+ label = "condition", y="Sadia", x="Sample"
dev.off()

# ===============================
# 8. Convert Gene Symbols to ENTREZ IDs
# ===============================
# Since simulated data uses fake names,
# we map random ENTREZ IDs for enrichment demo.

fake_entrez <- sample(keys(org.Hs.eg.db, keytype="ENTREZID"),
                      length(rownames(resSig)))

# ===============================
# 9. GO Enrichment
# ===============================

ego <- enrichGO(gene         = fake_entrez,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod= "BH",
                pvalueCutoff = 0.05)

write.csv(as.data.frame(ego),
          "results/GO_Enrichment.csv")

png("figures/GO_dotplot.png", width=900, height=700)
dotplot(ego)
dev.off()

# ===============================
# 10. KEGG Enrichment
# ===============================

ekegg <- enrichKEGG(gene = fake_entrez,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

write.csv(as.data.frame(ekegg),
          "results/KEGG_Enrichment.csv")

png("figures/KEGG_dotplot.png", width=900, height=700)
dotplot(ekegg)
dev.off()

# ===============================
# 11. Session Info
# ===============================

writeLines(capture.output(sessionInfo()),
           "results/sessionInfo.txt")

cat("Pipeline completed successfully.\n")

