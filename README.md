# -Differential-Gene-Expression-Analysis-RNA-Seq-
 Differential Gene Expression Analysis (RNA-Seq)
Computational Biology & Text Mining Projects in R
A collection of two independent analytical projects implemented in R, covering RNA-seq differential gene expression analysis and literary text mining using sentiment and visual analysis pipelines.

Project 1: Differential Gene Expression Analysis (RNA-Seq)
Overview
This project performs a full RNA-seq differential gene expression (DGE) analysis pipeline using publicly available data retrieved via GEOquery, processed through DESeq2, and visualized using ggplot2, pheatmap, and EnhancedVolcano. Pathway enrichment is conducted using clusterProfiler with KEGG annotations.
Workflow
Raw Count Data (GEO) --> DESeq2 Normalization --> DGE Analysis --> Visualization --> Pathway Enrichment
Key Analyses
Principal Component Analysis (PCA): Variance stabilized expression data (VST) was used to assess sample-level clustering. PC1 and PC2 together explained 26% of total variance. Notable: one Control sample appears as an outlier on PC1 (~45), which warrants inspection before final reporting.
Heatmap (Top 50 DEGs): Hierarchical clustering of the top 50 genes ranked by adjusted p-value reveals two broad expression clusters separating Control and Treatment groups.
Volcano Plot: Differentially expressed genes were filtered at |log2FC| > 1 and padj < 0.05. Significantly upregulated and downregulated genes are highlighted.
KEGG Pathway Enrichment: Over-representation analysis was performed using enrichKEGG() from clusterProfiler, with results exported as a dot plot and CSV.
Dependencies
# Bioconductor
BiocManager::install(c("DESeq2", "GEOquery", "clusterProfiler", "EnhancedVolcano"))

# CRAN
install.packages(c("ggplot2", "pheatmap", "dplyr", "tidyr"))
File Structure
Differential Gene Expression Analysis (RNA-Seq)/
├── adavance.R                  # Main analysis script
├── figures/
│   ├── PCA_plot.png
│   ├── Heatmap_top50.png
│   ├── Volcano_plot.png
│   └── KEGG_dotplot.png
├── results/
│   ├── KEGG_Enrichment.csv
│   └── sessionInfo.txt
└── README.md
Usage
Open adavance.R in RStudio and run sections sequentially. Ensure you have a stable internet connection for GEOquery data retrieval and KEGG database access.
# Example: Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = coldata,
                               design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
Notes and Known Issues
The PCA plot shows a potential outlier in the Control group (high PC1 value). This sample should be verified in the metadata before drawing biological conclusions.
fake_entrez IDs are used in the KEGG enrichment section, which means enrichment results from this run are not biologically valid. Replace with real Entrez gene IDs mapped from your DEG results using bitr() from clusterProfiler.
Sample size is small (n=12 total from metadata), which limits statistical power and should be explicitly acknowledged in any manuscript.
Citation
Love, M.I., Huber, W. & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550. https://doi.org/10.1186/s13059-014-0550-8
Yu, G., Wang, L.G., Han, Y. & He, Q.Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology, 16(5), 284–287. https://doi.org/10.1089/omi.2011.0118

Project 2: Text Mining & Sentiment Analysis — The Tinder-box
Overview
This project applies text mining and sentiment analysis to Hans Christian Andersen's fairy tale The Tinder-box using the tidytext and ggpage packages in R. The analysis visualizes word-level sentiment progression across the narrative using a rolling average approach, and maps longer words (9+ characters) across the full text layout using ggpage.
Key Analyses
Word-length visualization (ggpage): Highlights words of 9 or more characters in red across the full text layout, revealing the distribution of lexically complex words throughout the story.
Sentiment trajectory (AFINN + rolling average): Joins tokenized words to the AFINN sentiment lexicon and computes a rolling average (zoo::rollmean) across 500-word windows to track emotional arc across the narrative.
Animated sentiment plot (gganimate): Produces a transition_time() animated visualization of the sentiment rolling average across narrative positions.
Dependencies
install.packages(c("tidytext", "ggpage", "gganimate", "zoo", "dplyr", "ggplot2", "paletteer"))
Known Error and Fix
Error encountered:
Error in `filter()`:
! In argument: `lexicon == "nrc"`.
Caused by error:
! object 'lexicon' not found
Cause: This error occurs because get_sentiments("afinn") returns a tibble where the column is named word and value, not lexicon. When you try to filter by lexicon == "nrc" on an AFINN object, the column does not exist.
Fix:
# Wrong
sentiment_types <- sentiments %>%
  filter(lexicon == "nrc") %>%
  pull(sentiment) %>%
  unique()

# Correct — use get_sentiments() directly
nrc_sentiments <- get_sentiments("nrc")
sentiment_types <- nrc_sentiments %>%
  pull(sentiment) %>%
  unique()
File Structure
Learning with Ammar/
├── ggpage.R                    # Main text mining script
├── All_Figures_Complete.R      # Consolidated figure generation
└── README.md
Citation
Silge, J. & Robinson, D. (2017). Text Mining with R: A Tidy Approach. O'Reilly Media. https://www.tidytextmining.com
Hvitfeldt, E. (2019). ggpage: Create Page Layout Visualizations in R. CRAN. https://cran.r-project.org/package=ggpage

Session Info
Both projects were developed in R version 4.5.2 on Windows. Full session info is logged in results/sessionInfo.txt for the RNA-seq project.

Author
Sadia Noureen - Muzzamil Hussain - Abdullah Afzal Alvi
Researcher | Computational Biology & Bioinformatics

License
This repository is for educational and research purposes. Data sourced from NCBI GEO is subject to its respective data use policies.
