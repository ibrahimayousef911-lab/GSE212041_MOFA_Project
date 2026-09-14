# ==============================================================================
# Project: GSE212041 Longitudinal Multi-Omics Analysis
# Title: Integrating Neutrophil Transcriptomics and Plasma Proteomics Using MOFA2
#
# Dataset:
#   NCBI GEO: GSE212041
#
# Country:
#   Egypt
#
# Author:
#   Ibrahim Ashraf
#
# Analysis:
#   RNA-seq / Differential Expression / Longitudinal Transcriptomics /
#   Plasma Proteomics / Pathway Analysis / GSEA / PPI / MOFA2
#
# R version:
#   4.6.1
#
# Purpose:
#   Reproducible analysis of longitudinal neutrophil transcriptomics and
#   plasma proteomics with multi-omics integration using MOFA2.
#
# Important:
#   This repository is a cleaned organization of the original project script.
#   Final analytical methods/results are retained; repeated troubleshooting,
#   duplicate exploratory blocks, and transient debugging are not intended
#   as independent analyses.
# ==============================================================================

# STEP 16 - Differential Expression Analysis
# ============================================================

# ------------------------------------------------------------
# STEP 16A - Check DESeq2 Object
# ------------------------------------------------------------

# Check the number of samples
ncol(dds)

# Check the number of genes
nrow(dds)

# Check the experimental design
design(dds)

# Check sample metadata
head(as.data.frame(colData(dds)))

# Check patient categories
table(colData(dds)$patient_category)

# Check time points
table(colData(dds)$time_point)

# Check acuity
table(colData(dds)$acuity)

# ------------------------------------------------------------
# STEP 16B - Check DESeq2 Factor Levels
# ------------------------------------------------------------

levels(colData(dds)$patient_category)

table(colData(dds)$patient_category)

# ------------------------------------------------------------
# STEP 16C - Define Primary DE Comparison
# ------------------------------------------------------------

# Set COVID- symptomatic as the reference level
colData(dds)$patient_category <- relevel(
  colData(dds)$patient_category,
  ref = "COVID- symptomatic"
)

# Confirm factor levels
levels(colData(dds)$patient_category)

# Confirm sample counts
table(colData(dds)$patient_category)

# Confirm design
design(dds)

# ------------------------------------------------------------
# STEP 16D - Prepare Primary DESeq2 Comparison
# ------------------------------------------------------------

# Keep only COVID+ and COVID- symptomatic samples
keep_primary <- colData(dds)$patient_category %in% c(
  "COVID- symptomatic",
  "COVID+"
)

dds_primary <- dds[, keep_primary]

# Drop unused factor level
colData(dds_primary)$patient_category <- droplevels(
  colData(dds_primary)$patient_category
)

# Confirm sample numbers
ncol(dds_primary)

table(colData(dds_primary)$patient_category)

# Confirm design
design(dds_primary)
# ------------------------------------------------------------
# STEP 16E - Additional Gene Filtering for DE Analysis
# ------------------------------------------------------------

# Count matrix for the primary comparison
primary_counts <- counts(dds_primary)

# Keep genes expressed in at least 10% of the samples
min_samples_de <- ceiling(0.10 * ncol(dds_primary))

keep_de_genes <- rowSums(
  primary_counts >= 10
) >= min_samples_de

# Filter the DESeq2 object
dds_primary <- dds_primary[keep_de_genes, ]

# Check filtering results
cat("Genes before DE filtering:", nrow(primary_counts), "\n")
cat("Genes after DE filtering:", nrow(dds_primary), "\n")
cat("Genes removed:", sum(!keep_de_genes), "\n")

# Check final dimensions
dim(dds_primary)
# ------------------------------------------------------------
# STEP 16F - Save Primary DESeq2 Object Before Analysis
# ------------------------------------------------------------

if (!dir.exists("Results/Differential_Expression")) {
  dir.create(
    "Results/Differential_Expression",
    recursive = TRUE
  )
}

saveRDS(
  dds_primary,
  "Results/Differential_Expression/dds_primary_filtered.rds"
)

file.exists(
  "Results/Differential_Expression/dds_primary_filtered.rds"
)

dds_primary <- DESeq(dds_primary)


# ------------------------------------------------------------
# STEP 16G - Extract Differential Expression Results
# ------------------------------------------------------------

res_covid_vs_symptomatic <- results(
  dds_primary,
  contrast = c(
    "patient_category",
    "COVID+",
    "COVID- symptomatic"
  )
)

# Convert to data frame
res_covid_vs_symptomatic <- as.data.frame(
  res_covid_vs_symptomatic
)

# Add gene IDs
res_covid_vs_symptomatic$gene_id <- rownames(
  res_covid_vs_symptomatic
)

# Reorder columns
res_covid_vs_symptomatic <- res_covid_vs_symptomatic[
  ,
  c(
    "gene_id",
    setdiff(
      colnames(res_covid_vs_symptomatic),
      "gene_id"
    )
  )
]

# Sort by adjusted p-value
res_covid_vs_symptomatic <- res_covid_vs_symptomatic[
  order(res_covid_vs_symptomatic$padj),
]

# Check dimensions
dim(res_covid_vs_symptomatic)

# Show first results
head(res_covid_vs_symptomatic)
# ------------------------------------------------------------
# STEP 16H - Save Differential Expression Results
# ------------------------------------------------------------

write.csv(
  res_covid_vs_symptomatic,
  "Results/Differential_Expression/COVID_vs_Symptomatic_all_genes.csv",
  row.names = FALSE
)

saveRDS(
  res_covid_vs_symptomatic,
  "Results/Differential_Expression/COVID_vs_Symptomatic_all_genes.rds"
)

# Verify files
file.exists(
  "Results/Differential_Expression/COVID_vs_Symptomatic_all_genes.csv"
)

file.exists(
  "Results/Differential_Expression/COVID_vs_Symptomatic_all_genes.rds"
)

# ============================================================
# STEP 17 - Significant Differentially Expressed Genes
# ============================================================

# ------------------------------------------------------------
# STEP 17A - Define DEG Thresholds
# ------------------------------------------------------------

padj_cutoff <- 0.05
log2fc_cutoff <- 1

# ------------------------------------------------------------
# STEP 17B - Remove Genes Without Valid Statistics
# ------------------------------------------------------------

res_valid <- res_covid_vs_symptomatic[
  !is.na(res_covid_vs_symptomatic$padj) &
    !is.na(res_covid_vs_symptomatic$log2FoldChange),
]

# ------------------------------------------------------------
# STEP 17C - Identify Significant DEGs
# ------------------------------------------------------------

sig_degs <- res_valid[
  res_valid$padj < padj_cutoff &
    abs(res_valid$log2FoldChange) >= log2fc_cutoff,
]

# ------------------------------------------------------------
# STEP 17D - Separate Upregulated and Downregulated Genes
# ------------------------------------------------------------

upregulated_genes <- sig_degs[
  sig_degs$log2FoldChange >= log2fc_cutoff,
]

downregulated_genes <- sig_degs[
  sig_degs$log2FoldChange <= -log2fc_cutoff,
]

# ------------------------------------------------------------
# STEP 17E - Check DEG Counts
# ------------------------------------------------------------

cat("Total genes tested:", nrow(res_valid), "\n")
cat("Significant DEGs:", nrow(sig_degs), "\n")
cat("Upregulated genes:", nrow(upregulated_genes), "\n")
cat("Downregulated genes:", nrow(downregulated_genes), "\n")
# ------------------------------------------------------------
# STEP 17F - Save Significant DEG Tables
# ------------------------------------------------------------

write.csv(
  sig_degs,
  "Results/Differential_Expression/COVID_vs_Symptomatic_significant_DEGs.csv",
  row.names = FALSE
)

write.csv(
  upregulated_genes,
  "Results/Differential_Expression/COVID_vs_Symptomatic_upregulated.csv",
  row.names = FALSE
)

write.csv(
  downregulated_genes,
  "Results/Differential_Expression/COVID_vs_Symptomatic_downregulated.csv",
  row.names = FALSE
)
# ============================================================
# STEP 18 - Differential Expression Visualization
# ============================================================

# ------------------------------------------------------------
# STEP 18A - Prepare Volcano Plot
# ------------------------------------------------------------

# Create Figures directory
if (!dir.exists("Results/Figures")) {
  dir.create("Results/Figures", recursive = TRUE)
}

# Create plotting data
volcano_df <- res_valid

# Add significance classification
volcano_df$Significance <- "Not Significant"

volcano_df$Significance[
  volcano_df$padj < padj_cutoff &
    volcano_df$log2FoldChange >= log2fc_cutoff
] <- "Upregulated"

volcano_df$Significance[
  volcano_df$padj < padj_cutoff &
    volcano_df$log2FoldChange <= -log2fc_cutoff
] <- "Downregulated"

# Calculate -log10 adjusted p-value
volcano_df$neg_log10_padj <- -log10(
  volcano_df$padj
)

# Check categories
table(volcano_df$Significance)

saveRDS(
  sig_degs,
  "Results/Differential_Expression/COVID_vs_Symptomatic_significant_DEGs.rds"
)

# Verify
file.exists(
  "Results/Differential_Expression/COVID_vs_Symptomatic_significant_DEGs.csv"
)

file.exists(
  "Results/Differential_Expression/COVID_vs_Symptomatic_upregulated.csv"
)

file.exists(
  "Results/Differential_Expression/COVID_vs_Symptomatic_downregulated.csv"
)
# ------------------------------------------------------------
# STEP 18B - Volcano Plot
# ------------------------------------------------------------

library(ggplot2)

volcano_plot <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = Significance
  )
) +
  geom_point(
    size = 1.8,
    alpha = 0.7
  ) +
  geom_vline(
    xintercept = c(-log2fc_cutoff, log2fc_cutoff),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed"
  ) +
  labs(
    title = "Differential Expression: COVID+ vs COVID- symptomatic",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significance"
  ) +
  theme_classic()

print(volcano_plot)
# ------------------------------------------------------------
# STEP 18C - Save Volcano Plot
# ------------------------------------------------------------

ggsave(
  filename = "Results/Figures/Volcano_COVID_vs_Symptomatic.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

file.exists(
  "Results/Figures/Volcano_COVID_vs_Symptomatic.png"
)
# ------------------------------------------------------------
# STEP 18D - Prepare Top DEG Heatmap
# ------------------------------------------------------------

# Select top 50 significant DEGs
top_n <- 50

top_degs <- sig_degs[
  order(sig_degs$padj),
]

top_degs <- head(top_degs, top_n)

# Extract gene IDs
top_gene_ids <- top_degs$gene_id

# Extract normalized counts
heatmap_counts <- normalized_counts[
  top_gene_ids,
  ,
  drop = FALSE
]

# Log2 transform
heatmap_log2 <- log2(
  heatmap_counts + 1
)

# Check dimensions
dim(heatmap_log2)

# Check missing values
sum(is.na(heatmap_log2))

# Check selected genes
head(top_gene_ids)
# ------------------------------------------------------------
# STEP 18E - Top 50 DEG Heatmap
# ------------------------------------------------------------

# Load heatmap package
library(pheatmap)

# Create sample annotations
annotation_col <- data.frame(
  Patient_Category = metadata_deseq$patient_category,
  Time_Point = metadata_deseq$time_point,
  Acuity = metadata_deseq$acuity
)

# Make sure annotation row names match sample names
rownames(annotation_col) <- colnames(heatmap_log2)

# Scale each gene across samples
heatmap_scaled <- t(
  scale(t(heatmap_log2))
)

# Check dimensions
dim(heatmap_scaled)

# Check annotation dimensions
dim(annotation_col)

# Draw heatmap
heatmap_top50 <- pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 7,
  main = "Top 50 DEGs: COVID+ vs COVID- symptomatic"
)
# ------------------------------------------------------------
# STEP 18F - Save Top 50 DEG Heatmap
# ------------------------------------------------------------

pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 7,
  main = "Top 50 DEGs: COVID+ vs COVID- symptomatic",
  filename = "Results/Figures/Heatmap_Top50_DEGs_COVID_vs_Symptomatic.png",
  width = 12,
  height = 8
)

file.exists(
  "Results/Figures/Heatmap_Top50_DEGs_COVID_vs_Symptomatic.png"
)

# ============================================================
