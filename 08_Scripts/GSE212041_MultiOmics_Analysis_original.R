# ============================================================
# STEP 1 - DOWNLOAD GEO METADATA
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install GEOquery if needed
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

# Load GEOquery
library(GEOquery)

# Define GEO study ID
geo_id <- "GSE212041"

# Download GEO metadata from NCBI
gse <- getGEO(
  geo_id,
  GSEMatrix = TRUE
)

# Check the number of GEO objects
length(gse)

# Extract the GEO object
gse_object <- gse[[1]]

# Extract sample metadata
metadata <- pData(gse_object)

# Check metadata dimensions
dim(metadata)

# Show metadata column names
colnames(metadata)

# Show the first samples
head(metadata)


# ============================================================
# STEP 2 - CREATE CLEAN SAMPLE METADATA
# Project: GSE212041 Multi-Omics Analysis
# ============================================================

# Use GEO metadata from GPL24676
metadata <- pData(gse[[2]])

# Create sample ID from the GEO title
metadata$sample_id <- metadata$title

# Rename important clinical variables
metadata$patient_category <- metadata$`patient category:ch1`
metadata$covid_status     <- metadata$`covid-19 status:ch1`
metadata$time_point       <- metadata$`time point:ch1`
metadata$acuity           <- metadata$`acuity.max:ch1`
metadata$cell_type        <- metadata$`cell type:ch1`

# Keep the main variables needed for the project
metadata_clean <- metadata[, c(
  "sample_id",
  "geo_accession",
  "patient_category",
  "covid_status",
  "time_point",
  "acuity",
  "cell_type"
)]

# Save clean metadata
write.csv(
  metadata_clean,
  "03_Metadata/GSE212041_metadata_clean.csv",
  row.names = FALSE
)

# Check the result
head(metadata_clean)

# Check group counts
table(metadata_clean$patient_category)

# Check time points
table(metadata_clean$time_point)

# ============================================================
# STEP 3 - DEFINE AND CHECK RNA-SEQ FILES
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

# Define RNA-seq count matrix
count_file <- "GSE212041_Neutrophil_RNAseq_Count_Matrix.txt.gz"

# Define RNA-seq TPM matrix
tpm_file <- "GSE212041_Neutrophil_RNAseq_TPM_Matrix.txt.gz"

# Check that the count matrix exists
file.exists(count_file)

# Check that the TPM matrix exists
file.exists(tpm_file)

# ============================================================
# STEP 4 - LOAD RNA-SEQ COUNT MATRIX
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

counts <- read.table(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Check matrix dimensions
dim(counts)

# ============================================================
# STEP 5 - CHECK RNA-SEQ / METADATA SAMPLE MATCHING
# ============================================================

rna_samples <- colnames(counts)[-(1:2)]

metadata_samples <- metadata_clean$sample_id

missing_from_metadata <- setdiff(
  rna_samples,
  metadata_samples
)

length(missing_from_metadata)

missing_from_metadata

# ============================================================
# STEP 6 - CHECK THE 16 ADDITIONAL RNA-SEQ SAMPLES
# ============================================================

metadata_gpl18573 <- pData(gse[[1]])

# Show the available metadata columns
colnames(metadata_gpl18573)

# Show the sample titles and clinical information
metadata_gpl18573[, c(
  "geo_accession",
  "title",
  "patient category:ch1",
  "covid-19 status:ch1",
  "time point:ch1"
)]

# ============================================================
# STEP 7 - PREPARE ADDITIONAL SAMPLE METADATA
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

# Create sample ID from GEO title
metadata_gpl18573$sample_id <- metadata_gpl18573$title

# Rename important clinical variables
metadata_gpl18573$patient_category <-
  metadata_gpl18573$`patient category:ch1`

metadata_gpl18573$covid_status <-
  metadata_gpl18573$`covid-19 status:ch1`

metadata_gpl18573$time_point <-
  metadata_gpl18573$`time point:ch1`

metadata_gpl18573$acuity <-
  metadata_gpl18573$`acuity.max:ch1`

metadata_gpl18573$cell_type <-
  metadata_gpl18573$`cell type:ch1`

# Keep the same columns as the main metadata
metadata_extra <- metadata_gpl18573[, c(
  "sample_id",
  "geo_accession",
  "patient_category",
  "covid_status",
  "time_point",
  "acuity",
  "cell_type"
)]

# Check the 16 additional samples
metadata_extra

# ============================================================
# STEP 8 - MERGE COMPLETE SAMPLE METADATA
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

# Combine the main metadata with the additional metadata
metadata_all <- rbind(
  metadata_clean,
  metadata_extra
)

# Check total number of samples
nrow(metadata_all)

# Check for duplicated sample IDs
sum(duplicated(metadata_all$sample_id))

# Check metadata sample IDs against RNA sample IDs
length(setdiff(
  rna_samples,
  metadata_all$sample_id
))

# Check RNA samples that are not in metadata
setdiff(
  rna_samples,
  metadata_all$sample_id
)
# ============================================================
# STEP 9 - SAVE COMPLETE SAMPLE METADATA
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

# Save the complete metadata for all 781 RNA-seq samples
write.csv(
  metadata_all,
  "03_Metadata/GSE212041_metadata_all_781_samples.csv",
  row.names = FALSE
)

# Check that the file was created
file.exists(
  "03_Metadata/GSE212041_metadata_all_781_samples.csv"
)
# ============================================================
# STEP 10 - PREPARE RNA-SEQ COUNT MATRIX
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

# Extract gene IDs
gene_ids <- counts$Gene.ID

# Extract gene symbols
gene_symbols <- counts$Symbol

# Extract RNA-seq expression counts
count_matrix <- counts[, -(1:2)]

# Convert to matrix
count_matrix <- as.matrix(count_matrix)

# Assign gene IDs as row names
rownames(count_matrix) <- gene_ids

# Check dimensions
dim(count_matrix)

# Check sample names
colnames(count_matrix)[1:10]

# Check gene IDs
rownames(count_matrix)[1:10]

# ============================================================
# STEP 11 - RNA-SEQ QUALITY CONTROL
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

# Calculate total sequencing counts per sample
library_sizes <- colSums(count_matrix)

# Summary of library sizes
summary(library_sizes)

# Identify samples with very low library size
low_library_samples <- names(library_sizes[
  library_sizes < 1e6
])

# Number of low-library samples
length(low_library_samples)

# Display low-library samples
low_library_samples

# ============================================================
# STEP 12 - VISUALIZE LIBRARY SIZE DISTRIBUTION
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

library(ggplot2)

library_size_df <- data.frame(
  sample_id = names(library_sizes),
  library_size = as.numeric(library_sizes)
)

ggplot(library_size_df, aes(x = library_size)) +
  geom_histogram(bins = 40) +
  scale_x_log10() +
  labs(
    title = "RNA-seq Library Size Distribution",
    x = "Total Counts per Sample (log10 scale)",
    y = "Number of Samples"
  ) +
  theme_minimal()

# ============================================================
# STEP 13 - FILTER LOW-EXPRESSION GENES
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

# Define minimum expression threshold
min_count <- 10

# Define minimum number of samples
min_samples <- 10

# Keep genes expressed at sufficient levels
keep_genes <- rowSums(
  count_matrix >= min_count
) >= min_samples

# Create filtered count matrix
count_matrix_filtered <- count_matrix[keep_genes, ]

# Check filtering results
cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering:", nrow(count_matrix_filtered), "\n")
cat("Genes removed:", sum(!keep_genes), "\n")

# Check final dimensions
dim(count_matrix_filtered)

# ============================================================
# STEP 14 - DESeq2 Preparation and Normalization
# ============================================================

# ------------------------------------------------------------
# STEP 14A - Verify Final Metadata
# ------------------------------------------------------------

metadata_all <- read.csv(
  "GSE212041_metadata_final_781.csv",
  stringsAsFactors = FALSE
)

all(
  colnames(count_matrix_filtered) ==
    metadata_all$rna_sample_id
)

nrow(metadata_all)

sum(is.na(metadata_all$rna_sample_id))


# ------------------------------------------------------------
# STEP 14B - Create DESeq2 Object
# ------------------------------------------------------------

library(DESeq2)

metadata_deseq <- metadata_all

rownames(metadata_deseq) <- metadata_deseq$rna_sample_id

metadata_deseq$patient_category <- factor(
  metadata_deseq$patient_category
)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = metadata_deseq,
  design = ~ patient_category
)


# ------------------------------------------------------------
# STEP 14C - Estimate Size Factors
# ------------------------------------------------------------

dds <- estimateSizeFactors(dds)

head(sizeFactors(dds))


# ------------------------------------------------------------
# STEP 14D - Extract Normalized Counts
# ------------------------------------------------------------

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

dim(normalized_counts)


# ------------------------------------------------------------
# STEP 14E - Save Normalized Counts
# ------------------------------------------------------------

saveRDS(
  normalized_counts,
  "normalized_counts_28117x781.rds"
)

# ============================================================
# STEP 15 - Exploratory Data Analysis
# ============================================================

# ------------------------------------------------------------
# STEP 15A - Check Sample Groups
# ------------------------------------------------------------

table(metadata_deseq$patient_category)

table(metadata_deseq$time_point)

table(metadata_deseq$acuity)

# ------------------------------------------------------------
# STEP 15B - PCA Preparation
# ------------------------------------------------------------

# Log2-transform normalized counts
log2_counts <- log2(normalized_counts + 1)

# Check dimensions
dim(log2_counts)

# Check for missing values
sum(is.na(log2_counts))

# Check for infinite values
sum(is.infinite(log2_counts))

# ------------------------------------------------------------
# STEP 15C - Principal Component Analysis (PCA)
# ------------------------------------------------------------

# Transpose the expression matrix
pca_input <- t(log2_counts)

# Run PCA
pca_result <- prcomp(
  pca_input,
  center = TRUE,
  scale. = FALSE
)

# Check PCA dimensions
dim(pca_result$x)

# ------------------------------------------------------------
# STEP 15D - PCA Plot by Patient Category
# ------------------------------------------------------------

# Create PCA data frame
pca_df <- data.frame(
  Sample = rownames(pca_result$x),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Patient_Category = metadata_deseq$patient_category
)

# Calculate variance explained
PC1_var <- round(
  100 * summary(pca_result)$importance[2, 1],
  2
)

PC2_var <- round(
  100 * summary(pca_result)$importance[2, 2],
  2
)

# Plot PCA
library(ggplot2)

pca_plot <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Patient_Category
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.7
  ) +
  labs(
    title = "PCA of GSE212041 Neutrophil RNA-seq",
    x = paste0("PC1 (", PC1_var, "%)"),
    y = paste0("PC2 (", PC2_var, "%)"),
    color = "Patient Category"
  ) +
  theme_classic()

print(pca_plot)

# ------------------------------------------------------------
# STEP 15E - PCA Plot by Time Point
# ------------------------------------------------------------

pca_df$Time_Point <- metadata_deseq$time_point

pca_plot_time <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Time_Point
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.7
  ) +
  labs(
    title = "PCA of GSE212041 Neutrophil RNA-seq by Time Point",
    x = paste0("PC1 (", PC1_var, "%)"),
    y = paste0("PC2 (", PC2_var, "%)"),
    color = "Time Point"
  ) +
  theme_classic()

print(pca_plot_time)

# ------------------------------------------------------------
# STEP 15F - Save PCA Plots
# ------------------------------------------------------------

ggsave(
  filename = "PCA_by_Patient_Category.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "PCA_by_Time_Point.png",
  plot = pca_plot_time,
  width = 8,
  height = 6,
  dpi = 300
)
file.exists("PCA_by_Patient_Category.png")
file.exists("PCA_by_Time_Point.png")

# ------------------------------------------------------------
# STEP 15G - PCA Plot by Acuity
# ------------------------------------------------------------

pca_df$Acuity <- metadata_deseq$acuity

pca_plot_acuity <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Acuity
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.7
  ) +
  labs(
    title = "PCA of GSE212041 Neutrophil RNA-seq by Acuity",
    x = paste0("PC1 (", PC1_var, "%)"),
    y = paste0("PC2 (", PC2_var, "%)"),
    color = "Acuity"
  ) +
  theme_classic()

print(pca_plot_acuity)
# ------------------------------------------------------------
# STEP 15H - Save PCA by Acuity
# ------------------------------------------------------------

ggsave(
  filename = "PCA_by_Acuity.png",
  plot = pca_plot_acuity,
  width = 8,
  height = 6,
  dpi = 300
)

file.exists("PCA_by_Acuity.png")
# Create results directory
if (!dir.exists("Results")) {
  dir.create("Results")
}

if (!dir.exists("Results/PCA")) {
  dir.create("Results/PCA", recursive = TRUE)
}

ggsave(
  filename = "Results/PCA/PCA_by_Patient_Category.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "Results/PCA/PCA_by_Time_Point.png",
  plot = pca_plot_time,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "Results/PCA/PCA_by_Acuity.png",
  plot = pca_plot_acuity,
  width = 8,
  height = 6,
  dpi = 300
)
# ============================================================
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
# STEP 19 - Functional Enrichment Analysis
# ============================================================

# ------------------------------------------------------------
# STEP 19A - Gene ID Conversion
# ------------------------------------------------------------

# Install/load required Bioconductor packages if needed

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE)
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
}

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Clean Ensembl IDs
# ------------------------------------------------------------

all_gene_ids <- sig_degs$gene_id
up_gene_ids <- upregulated_genes$gene_id
down_gene_ids <- downregulated_genes$gene_id

# Remove Ensembl version numbers if present
all_gene_ids <- sub("\\..*$", "", all_gene_ids)
up_gene_ids <- sub("\\..*$", "", up_gene_ids)
down_gene_ids <- sub("\\..*$", "", down_gene_ids)

# Remove duplicates
all_gene_ids <- unique(all_gene_ids)
up_gene_ids <- unique(up_gene_ids)
down_gene_ids <- unique(down_gene_ids)

# ------------------------------------------------------------
# Convert Ensembl IDs → Entrez IDs
# ------------------------------------------------------------

all_id_map <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_id_map <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_id_map <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# ------------------------------------------------------------
# Check conversion
# ------------------------------------------------------------

cat("Total significant DEGs:", length(all_gene_ids), "\n")
cat("DEGs mapped to Entrez:", nrow(all_id_map), "\n")

cat("Upregulated genes:", length(up_gene_ids), "\n")
cat("Upregulated mapped:", nrow(up_id_map), "\n")

cat("Downregulated genes:", length(down_gene_ids), "\n")
cat("Downregulated mapped:", nrow(down_id_map), "\n")

# Show first mappings
head(all_id_map)

# ============================================================
# STEP 19A - Check Enrichment Packages
# ============================================================

library(BiocManager)

BiocManager::valid()

requireNamespace("GO.db", quietly = TRUE)
requireNamespace("org.Hs.eg.db", quietly = TRUE)
requireNamespace("clusterProfiler", quietly = TRUE)
requireNamespace("DOSE", quietly = TRUE)
requireNamespace("enrichplot", quietly = TRUE)
# ============================================================
# STEP 19B - Install Enrichment Analysis Packages
# ============================================================

BiocManager::install(
  "clusterProfiler",
  ask = FALSE,
  update = FALSE
)

BiocManager::install(
  "DOSE",
  ask = FALSE,
  update = FALSE
)

BiocManager::install(
  "enrichplot",
  ask = FALSE,
  update = FALSE
)
# ============================================================
# STEP 19C - Ensembl to Entrez ID Mapping
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# Get Ensembl IDs
all_gene_ids <- unique(sig_degs$gene_id)
up_gene_ids <- unique(upregulated_genes$gene_id)
down_gene_ids <- unique(downregulated_genes$gene_id)

# Remove Ensembl version numbers if present
all_gene_ids <- sub("\\..*$", "", all_gene_ids)
up_gene_ids <- sub("\\..*$", "", up_gene_ids)
down_gene_ids <- sub("\\..*$", "", down_gene_ids)

# Convert Ensembl → Entrez
all_id_map <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_id_map <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_id_map <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Create unique Entrez gene lists
all_entrez <- unique(all_id_map$ENTREZID)
up_entrez <- unique(up_id_map$ENTREZID)
down_entrez <- unique(down_id_map$ENTREZID)

# Check mapping
cat("Significant DEGs:", length(all_gene_ids), "\n")
cat("Mapped Entrez IDs:", length(all_entrez), "\n\n")

cat("Upregulated genes:", length(up_gene_ids), "\n")
cat("Mapped Up Entrez IDs:", length(up_entrez), "\n\n")

cat("Downregulated genes:", length(down_gene_ids), "\n")
cat("Mapped Down Entrez IDs:", length(down_entrez), "\n")

# Preview mappings
head(all_id_map)
# ============================================================
# STEP 19D - GO Biological Process Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Create background gene list
# ------------------------------------------------------------

background_ensembl <- unique(res_valid$gene_id)

background_ensembl <- sub(
  "\\..*$",
  "",
  background_ensembl
)

background_map <- bitr(
  background_ensembl,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_entrez <- unique(
  background_map$ENTREZID
)

cat(
  "Background genes tested:",
  length(background_ensembl),
  "\n"
)

cat(
  "Background Entrez IDs:",
  length(background_entrez),
  "\n"
)

# ------------------------------------------------------------
# GO-BP enrichment - All significant DEGs
# ------------------------------------------------------------

ego_all <- enrichGO(
  gene = all_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP enrichment - Upregulated genes
# ------------------------------------------------------------

ego_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP enrichment - Downregulated genes
# ------------------------------------------------------------

ego_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# Check results
# ------------------------------------------------------------

cat(
  "GO-BP terms - All DEGs:",
  nrow(as.data.frame(ego_all)),
  "\n"
)

cat(
  "GO-BP terms - Upregulated:",
  nrow(as.data.frame(ego_up)),
  "\n"
)

cat(
  "GO-BP terms - Downregulated:",
  nrow(as.data.frame(ego_down)),
  "\n"
)

# ============================================================
# STEP 19D - Prepare GO Background
# ============================================================

background_ensembl <- unique(
  res_covid_vs_symptomatic$gene_id
)

background_ensembl <- sub(
  "\\..*$",
  "",
  background_ensembl
)

background_map <- bitr(
  background_ensembl,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_entrez <- unique(
  background_map$ENTREZID
)

cat(
  "Background Ensembl genes:",
  length(background_ensembl),
  "\n"
)

cat(
  "Background Entrez genes:",
  length(background_entrez),
  "\n"
)
# ============================================================
# STEP 19E - GO Biological Process Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# GO-BP: All Significant DEGs
# ------------------------------------------------------------

ego_all <- enrichGO(
  gene = all_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP: Upregulated Genes
# ------------------------------------------------------------

ego_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP: Downregulated Genes
# ------------------------------------------------------------

ego_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# Check Results
# ------------------------------------------------------------

cat(
  "GO-BP terms - All DEGs:",
  nrow(as.data.frame(ego_all)),
  "\n"
)

cat(
  "GO-BP terms - Upregulated:",
  nrow(as.data.frame(ego_up)),
  "\n"
)

cat(
  "GO-BP terms - Downregulated:",
  nrow(as.data.frame(ego_down)),
  "\n"
)
# ============================================================
# STEP 19E - Recreate Entrez Gene Lists
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# Significant DEG lists
all_gene_ids <- unique(sig_degs$gene_id)
up_gene_ids <- unique(upregulated_genes$gene_id)
down_gene_ids <- unique(downregulated_genes$gene_id)

# Remove Ensembl version numbers
all_gene_ids <- sub("\\..*$", "", all_gene_ids)
up_gene_ids <- sub("\\..*$", "", up_gene_ids)
down_gene_ids <- sub("\\..*$", "", down_gene_ids)

# Ensembl -> Entrez
all_id_map <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_id_map <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_id_map <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Unique Entrez IDs
all_entrez <- unique(all_id_map$ENTREZID)
up_entrez <- unique(up_id_map$ENTREZID)
down_entrez <- unique(down_id_map$ENTREZID)

# Check
cat("All DEGs:", length(all_gene_ids), "\n")
cat("Mapped all:", length(all_entrez), "\n")

cat("Upregulated:", length(up_gene_ids), "\n")
cat("Mapped up:", length(up_entrez), "\n")

cat("Downregulated:", length(down_gene_ids), "\n")
cat("Mapped down:", length(down_entrez), "\n")

# ============================================================
# STEP 19D - Recreate Entrez IDs for Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# Ensembl IDs
all_gene_ids <- unique(sig_degs$gene_id)
up_gene_ids <- unique(upregulated_genes$gene_id)
down_gene_ids <- unique(downregulated_genes$gene_id)

# Remove Ensembl version numbers
all_gene_ids <- sub("\\..*$", "", all_gene_ids)
up_gene_ids <- sub("\\..*$", "", up_gene_ids)
down_gene_ids <- sub("\\..*$", "", down_gene_ids)

# Ensembl -> Entrez
all_id_map <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_id_map <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_id_map <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Unique Entrez IDs
all_entrez <- unique(all_id_map$ENTREZID)
up_entrez <- unique(up_id_map$ENTREZID)
down_entrez <- unique(down_id_map$ENTREZID)

# Check
cat("All DEGs:", length(all_gene_ids), "\n")
cat("Mapped all:", length(all_entrez), "\n")
cat("Upregulated:", length(up_gene_ids), "\n")
cat("Mapped up:", length(up_entrez), "\n")
cat("Downregulated:", length(down_gene_ids), "\n")
cat("Mapped down:", length(down_entrez), "\n")
# ============================================================
# STEP 19E - GO Biological Process Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# GO-BP: All Significant DEGs
# ------------------------------------------------------------

ego_all <- enrichGO(
  gene = all_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP: Upregulated Genes
# ------------------------------------------------------------

ego_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP: Downregulated Genes
# ------------------------------------------------------------

ego_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# Check Results
# ------------------------------------------------------------

cat(
  "GO-BP terms - All DEGs:",
  nrow(as.data.frame(ego_all)),
  "\n"
)

cat(
  "GO-BP terms - Upregulated:",
  nrow(as.data.frame(ego_up)),
  "\n"
)

cat(
  "GO-BP terms - Downregulated:",
  nrow(as.data.frame(ego_down)),
  "\n"
)
# ============================================================
# STEP 19F - Inspect and Save GO-BP Results
# ============================================================

# Create GO results directory
if (!dir.exists("Results/Enrichment/GO")) {
  dir.create("Results/Enrichment/GO", recursive = TRUE)
}

# Convert enrichment results to data frames
go_all_df <- as.data.frame(ego_all)
go_up_df <- as.data.frame(ego_up)
go_down_df <- as.data.frame(ego_down)

# Check number of enriched terms
cat("GO-BP - All DEGs:", nrow(go_all_df), "\n")
cat("GO-BP - Upregulated:", nrow(go_up_df), "\n")
cat("GO-BP - Downregulated:", nrow(go_down_df), "\n")

# Show top 10 terms
cat("\n===== TOP 10 GO-BP: UPREGULATED =====\n")
print(
  go_up_df[
    order(go_up_df$p.adjust),
    c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "Count")
  ][1:min(10, nrow(go_up_df)), ]
)

cat("\n===== TOP 10 GO-BP: DOWNREGULATED =====\n")
print(
  go_down_df[
    order(go_down_df$p.adjust),
    c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "Count")
  ][1:min(10, nrow(go_down_df)), ]
)

# Save tables
write.csv(
  go_all_df,
  "Results/Enrichment/GO/GO_BP_all_DEGs.csv",
  row.names = FALSE
)

write.csv(
  go_up_df,
  "Results/Enrichment/GO/GO_BP_upregulated.csv",
  row.names = FALSE
)

write.csv(
  go_down_df,
  "Results/Enrichment/GO/GO_BP_downregulated.csv",
  row.names = FALSE
)

# Save RDS objects
saveRDS(
  ego_all,
  "Results/Enrichment/GO/GO_BP_all_DEGs.rds"
)

saveRDS(
  ego_up,
  "Results/Enrichment/GO/GO_BP_upregulated.rds"
)

saveRDS(
  ego_down,
  "Results/Enrichment/GO/GO_BP_downregulated.rds"
)

# Verify files
cat("\n===== FILE CHECK =====\n")

cat(
  "All GO:",
  file.exists("Results/Enrichment/GO/GO_BP_all_DEGs.csv"),
  "\n"
)

cat(
  "Up GO:",
  file.exists("Results/Enrichment/GO/GO_BP_upregulated.csv"),
  "\n"
)

cat(
  "Down GO:",
  file.exists("Results/Enrichment/GO/GO_BP_downregulated.csv"),
  "\n"
)
# ============================================================
# STEP 19G - GO-BP Visualization
# ============================================================

library(enrichplot)
library(ggplot2)

# Create GO figures directory
if (!dir.exists("Results/Enrichment/GO/Figures")) {
  dir.create(
    "Results/Enrichment/GO/Figures",
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# UPREGULATED GO-BP
# ------------------------------------------------------------

p_go_up <- dotplot(
  ego_up,
  showCategory = 20
) +
  ggtitle("GO Biological Process - Upregulated Genes")

print(p_go_up)

ggsave(
  "Results/Enrichment/GO/Figures/GO_BP_Upregulated_Dotplot.png",
  p_go_up,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# DOWNREGULATED GO-BP
# ------------------------------------------------------------

p_go_down <- dotplot(
  ego_down,
  showCategory = 14
) +
  ggtitle("GO Biological Process - Downregulated Genes")

print(p_go_down)

ggsave(
  "Results/Enrichment/GO/Figures/GO_BP_Downregulated_Dotplot.png",
  p_go_down,
  width = 10,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# ALL DEGs GO-BP
# ------------------------------------------------------------

p_go_all <- dotplot(
  ego_all,
  showCategory = 20
) +
  ggtitle("GO Biological Process - All Significant DEGs")

print(p_go_all)

ggsave(
  "Results/Enrichment/GO/Figures/GO_BP_All_DEGs_Dotplot.png",
  p_go_all,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Verify figures
# ------------------------------------------------------------

cat(
  "Upregulated figure:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_BP_Upregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "Downregulated figure:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_BP_Downregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "All DEGs figure:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_BP_All_DEGs_Dotplot.png"
  ),
  "\n"
)

# ============================================================
# STEP 19 - KEGG Pathway Enrichment Analysis
# ============================================================

library(clusterProfiler)
library(ggplot2)

# ------------------------------------------------------------
# STEP 19A - Create KEGG results directory
# ------------------------------------------------------------

if (!dir.exists("Results/Enrichment/KEGG")) {
  dir.create(
    "Results/Enrichment/KEGG",
    recursive = TRUE
  )
}

if (!dir.exists("Results/Enrichment/KEGG/Figures")) {
  dir.create(
    "Results/Enrichment/KEGG/Figures",
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# STEP 19B - Check input gene lists
# ------------------------------------------------------------

cat("All significant Entrez genes:",
    length(all_entrez), "\n")

cat("Upregulated Entrez genes:",
    length(up_entrez), "\n")

cat("Downregulated Entrez genes:",
    length(down_entrez), "\n")

cat("Background Entrez genes:",
    length(background_entrez), "\n")


# ------------------------------------------------------------
# STEP 19C - KEGG enrichment: ALL DEGs
# ------------------------------------------------------------

ekegg_all <- enrichKEGG(
  gene = all_entrez,
  universe = background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500
)

# ------------------------------------------------------------
# STEP 19D - KEGG enrichment: UPREGULATED
# ------------------------------------------------------------

ekegg_up <- enrichKEGG(
  gene = up_entrez,
  universe = background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500
)

# ------------------------------------------------------------
# STEP 19E - KEGG enrichment: DOWNREGULATED
# ------------------------------------------------------------

ekegg_down <- enrichKEGG(
  gene = down_entrez,
  universe = background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 500
)

# ------------------------------------------------------------
# STEP 19F - Convert to data frames
# ------------------------------------------------------------

kegg_all_df <- as.data.frame(ekegg_all)
kegg_up_df <- as.data.frame(ekegg_up)
kegg_down_df <- as.data.frame(ekegg_down)

# ------------------------------------------------------------
# STEP 19G - Check number of pathways
# ------------------------------------------------------------

cat("\n===== KEGG RESULTS =====\n")

cat(
  "KEGG pathways - All DEGs:",
  nrow(kegg_all_df),
  "\n"
)

cat(
  "KEGG pathways - Upregulated:",
  nrow(kegg_up_df),
  "\n"
)

cat(
  "KEGG pathways - Downregulated:",
  nrow(kegg_down_df),
  "\n"
)

# ------------------------------------------------------------
# STEP 19H - Show top pathways
# ------------------------------------------------------------

if (nrow(kegg_all_df) > 0) {
  
  cat("\n===== TOP KEGG - ALL DEGs =====\n")
  
  print(
    kegg_all_df[
      order(kegg_all_df$p.adjust),
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ][1:min(15, nrow(kegg_all_df)), ]
  )
}


if (nrow(kegg_up_df) > 0) {
  
  cat("\n===== TOP KEGG - UPREGULATED =====\n")
  
  print(
    kegg_up_df[
      order(kegg_up_df$p.adjust),
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ][1:min(15, nrow(kegg_up_df)), ]
  )
}


if (nrow(kegg_down_df) > 0) {
  
  cat("\n===== TOP KEGG - DOWNREGULATED =====\n")
  
  print(
    kegg_down_df[
      order(kegg_down_df$p.adjust),
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ][1:min(15, nrow(kegg_down_df)), ]
  )
}

# ------------------------------------------------------------
# STEP 19I - Save KEGG tables
# ------------------------------------------------------------

write.csv(
  kegg_all_df,
  "Results/Enrichment/KEGG/KEGG_all_DEGs.csv",
  row.names = FALSE
)

write.csv(
  kegg_up_df,
  "Results/Enrichment/KEGG/KEGG_upregulated.csv",
  row.names = FALSE
)

write.csv(
  kegg_down_df,
  "Results/Enrichment/KEGG/KEGG_downregulated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# STEP 19J - Save RDS objects
# ------------------------------------------------------------

saveRDS(
  ekegg_all,
  "Results/Enrichment/KEGG/KEGG_all_DEGs.rds"
)

saveRDS(
  ekegg_up,
  "Results/Enrichment/KEGG/KEGG_upregulated.rds"
)

saveRDS(
  ekegg_down,
  "Results/Enrichment/KEGG/KEGG_downregulated.rds"
)

# ------------------------------------------------------------
# STEP 19K - Verify files
# ------------------------------------------------------------

cat("\n===== KEGG FILE CHECK =====\n")

cat(
  "All KEGG:",
  file.exists(
    "Results/Enrichment/KEGG/KEGG_all_DEGs.csv"
  ),
  "\n"
)

cat(
  "Up KEGG:",
  file.exists(
    "Results/Enrichment/KEGG/KEGG_upregulated.csv"
  ),
  "\n"
)

cat(
  "Down KEGG:",
  file.exists(
    "Results/Enrichment/KEGG/KEGG_downregulated.csv"
  ),
  "\n"
)
# ============================================================
# STEP 20 - KEGG Enrichment Figures
# ============================================================

library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# STEP 20A - Create Figures directory
# ------------------------------------------------------------

if (!dir.exists("Results/Enrichment/KEGG/Figures")) {
  dir.create(
    "Results/Enrichment/KEGG/Figures",
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# STEP 20B - KEGG All DEGs
# ------------------------------------------------------------

p_kegg_all <- dotplot(
  ekegg_all,
  showCategory = min(20, nrow(as.data.frame(ekegg_all)))
) +
  ggtitle("KEGG Pathway Enrichment - All Significant DEGs")

print(p_kegg_all)

ggsave(
  "Results/Enrichment/KEGG/Figures/KEGG_All_DEGs_Dotplot.png",
  p_kegg_all,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# STEP 20C - KEGG Upregulated
# ------------------------------------------------------------

if (nrow(as.data.frame(ekegg_up)) > 0) {
  
  p_kegg_up <- dotplot(
    ekegg_up,
    showCategory = min(20, nrow(as.data.frame(ekegg_up)))
  ) +
    ggtitle("KEGG Pathway Enrichment - Upregulated Genes")
  
  print(p_kegg_up)
  
  ggsave(
    "Results/Enrichment/KEGG/Figures/KEGG_Upregulated_Dotplot.png",
    p_kegg_up,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ------------------------------------------------------------
# STEP 20D - KEGG Downregulated
# ------------------------------------------------------------

if (nrow(as.data.frame(ekegg_down)) > 0) {
  
  p_kegg_down <- dotplot(
    ekegg_down,
    showCategory = min(20, nrow(as.data.frame(ekegg_down)))
  ) +
    ggtitle("KEGG Pathway Enrichment - Downregulated Genes")
  
  print(p_kegg_down)
  
  ggsave(
    "Results/Enrichment/KEGG/Figures/KEGG_Downregulated_Dotplot.png",
    p_kegg_down,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ------------------------------------------------------------
# STEP 20E - Verify Figures
# ------------------------------------------------------------

cat("\n===== KEGG FIGURE CHECK =====\n")

cat(
  "All KEGG figure:",
  file.exists(
    "Results/Enrichment/KEGG/Figures/KEGG_All_DEGs_Dotplot.png"
  ),
  "\n"
)

cat(
  "Up KEGG figure:",
  file.exists(
    "Results/Enrichment/KEGG/Figures/KEGG_Upregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "Down KEGG figure:",
  file.exists(
    "Results/Enrichment/KEGG/Figures/KEGG_Downregulated_Dotplot.png"
  ),
  "\n"
)
# ============================================================
# STEP 21A - Check ReactomePA
# ============================================================

cat(
  "ReactomePA:",
  requireNamespace("ReactomePA", quietly = TRUE),
  "\n"
)
# ============================================================
# STEP 21A - Install ReactomePA
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("ReactomePA", ask = FALSE, update = FALSE)

# ============================================================
# STEP 22 - Restore Gene IDs and Background
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Define gene IDs
# ------------------------------------------------------------

all_gene_ids <- unique(sig_degs$gene_id)
up_gene_ids <- unique(upregulated_genes$gene_id)
down_gene_ids <- unique(downregulated_genes$gene_id)

# ------------------------------------------------------------
# Convert Ensembl -> Entrez
# ------------------------------------------------------------

all_id_map <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_id_map <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_id_map <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# ------------------------------------------------------------
# Unique Entrez IDs
# ------------------------------------------------------------

all_entrez <- unique(all_id_map$ENTREZID)
up_entrez <- unique(up_id_map$ENTREZID)
down_entrez <- unique(down_id_map$ENTREZID)

# ------------------------------------------------------------
# Background
# ------------------------------------------------------------

background_ensembl <- unique(res_covid_vs_symptomatic$gene_id)

background_id_map <- bitr(
  background_ensembl,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_entrez <- unique(
  background_id_map$ENTREZID
)

# ------------------------------------------------------------
# Check
# ------------------------------------------------------------

cat("\n===== GENE ID CHECK =====\n")

cat(
  "Significant DEGs:",
  length(all_gene_ids),
  "\n"
)

cat(
  "Mapped all Entrez:",
  length(all_entrez),
  "\n"
)

cat(
  "Upregulated:",
  length(up_gene_ids),
  "\n"
)

cat(
  "Mapped up Entrez:",
  length(up_entrez),
  "\n"
)

cat(
  "Downregulated:",
  length(down_gene_ids),
  "\n"
)

cat(
  "Mapped down Entrez:",
  length(down_entrez),
  "\n"
)

cat(
  "Background Ensembl:",
  length(background_ensembl),
  "\n"
)

cat(
  "Background Entrez:",
  length(background_entrez),
  "\n"
)
# ============================================================
# STEP 23 - Install ReactomePA
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  "ReactomePA",
  ask = FALSE,
  update = FALSE
)
# ============================================================
# STEP 23A - Install missing Reactome dependency
# ============================================================

BiocManager::install(
  "reactome.db",
  ask = FALSE,
  update = FALSE,
  force = TRUE
)
# ============================================================
# STEP 24 - Reactome Pathway Enrichment
# ============================================================

library(ReactomePA)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# 1. Reactome enrichment - All significant DEGs
# ------------------------------------------------------------

reactome_all <- enrichPathway(
  gene = all_entrez,
  universe = background_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 2. Reactome enrichment - Upregulated genes
# ------------------------------------------------------------

reactome_up <- enrichPathway(
  gene = up_entrez,
  universe = background_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 3. Reactome enrichment - Downregulated genes
# ------------------------------------------------------------

reactome_down <- enrichPathway(
  gene = down_entrez,
  universe = background_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 4. Check results
# ------------------------------------------------------------

cat("\n===== REACTOME RESULTS =====\n")
cat("Reactome All pathways:", nrow(as.data.frame(reactome_all)), "\n")
cat("Reactome Up pathways:", nrow(as.data.frame(reactome_up)), "\n")
cat("Reactome Down pathways:", nrow(as.data.frame(reactome_down)), "\n")

# ============================================================
# STEP 25 - Save Reactome Results
# ============================================================

# Create output directories
dir.create(
  "Results/Enrichment/Reactome",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/Enrichment/Reactome/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Save tables
# ------------------------------------------------------------

write.csv(
  as.data.frame(reactome_all),
  "Results/Enrichment/Reactome/Reactome_all_DEGs.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(reactome_up),
  "Results/Enrichment/Reactome/Reactome_upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(reactome_down),
  "Results/Enrichment/Reactome/Reactome_downregulated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Save RDS objects
# ------------------------------------------------------------

saveRDS(
  reactome_all,
  "Results/Enrichment/Reactome/Reactome_all_DEGs.rds"
)

saveRDS(
  reactome_up,
  "Results/Enrichment/Reactome/Reactome_upregulated.rds"
)

saveRDS(
  reactome_down,
  "Results/Enrichment/Reactome/Reactome_downregulated.rds"
)

# ------------------------------------------------------------
# Verify files
# ------------------------------------------------------------

cat("\n===== REACTOME FILE CHECK =====\n")

cat(
  "All CSV:",
  file.exists("Results/Enrichment/Reactome/Reactome_all_DEGs.csv"),
  "\n"
)

cat(
  "Up CSV:",
  file.exists("Results/Enrichment/Reactome/Reactome_upregulated.csv"),
  "\n"
)

cat(
  "Down CSV:",
  file.exists("Results/Enrichment/Reactome/Reactome_downregulated.csv"),
  "\n"
)

cat(
  "All RDS:",
  file.exists("Results/Enrichment/Reactome/Reactome_all_DEGs.rds"),
  "\n"
)

cat(
  "Up RDS:",
  file.exists("Results/Enrichment/Reactome/Reactome_upregulated.rds"),
  "\n"
)

cat(
  "Down RDS:",
  file.exists("Results/Enrichment/Reactome/Reactome_downregulated.rds"),
  "\n"
)

# ============================================================
# STEP 26 - Reactome Dotplots
# ============================================================

library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# 1. All significant DEGs
# ------------------------------------------------------------

reactome_all_plot <- dotplot(
  reactome_all,
  showCategory = 20
) +
  ggtitle("Reactome Pathway Enrichment - All Significant DEGs")

ggsave(
  "Results/Enrichment/Reactome/Figures/Reactome_All_DEGs_Dotplot.png",
  reactome_all_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 2. Upregulated genes
# ------------------------------------------------------------

reactome_up_plot <- dotplot(
  reactome_up,
  showCategory = 20
) +
  ggtitle("Reactome Pathway Enrichment - Upregulated Genes")

ggsave(
  "Results/Enrichment/Reactome/Figures/Reactome_Upregulated_Dotplot.png",
  reactome_up_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 3. Downregulated genes
# ------------------------------------------------------------

if (nrow(as.data.frame(reactome_down)) > 0) {
  
  reactome_down_plot <- dotplot(
    reactome_down,
    showCategory = 20
  ) +
    ggtitle("Reactome Pathway Enrichment - Downregulated Genes")
  
  ggsave(
    "Results/Enrichment/Reactome/Figures/Reactome_Downregulated_Dotplot.png",
    reactome_down_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  
  cat(
    "\nNo significant Reactome pathways for downregulated genes.\n"
  )
}

# ------------------------------------------------------------
# 4. Verify figures
# ------------------------------------------------------------

cat("\n===== REACTOME FIGURE CHECK =====\n")

cat(
  "All:",
  file.exists(
    "Results/Enrichment/Reactome/Figures/Reactome_All_DEGs_Dotplot.png"
  ),
  "\n"
)

cat(
  "Up:",
  file.exists(
    "Results/Enrichment/Reactome/Figures/Reactome_Upregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "Down:",
  file.exists(
    "Results/Enrichment/Reactome/Figures/Reactome_Downregulated_Dotplot.png"
  ),
  "\n"
)
# ============================================================
# STEP 27 - Top Reactome Pathways
# ============================================================

# Top 20 Reactome pathways - All DEGs
top_reactome_all <- as.data.frame(reactome_all) |>
  dplyr::arrange(p.adjust) |>
  dplyr::select(
    ID,
    Description,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    geneID
  ) |>
  head(20)

# Top 20 Reactome pathways - Upregulated
top_reactome_up <- as.data.frame(reactome_up) |>
  dplyr::arrange(p.adjust) |>
  dplyr::select(
    ID,
    Description,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    geneID
  ) |>
  head(20)

# Print results
cat("\n===== TOP REACTOME - ALL DEGs =====\n")
print(top_reactome_all)

cat("\n===== TOP REACTOME - UPREGULATED =====\n")
print(top_reactome_up)
# ============================================================
# STEP 28 - Save Top Reactome Pathways
# ============================================================

write.csv(
  top_reactome_all,
  "Results/Enrichment/Reactome/Top20_Reactome_All_DEGs.csv",
  row.names = FALSE
)

write.csv(
  top_reactome_up,
  "Results/Enrichment/Reactome/Top20_Reactome_Upregulated.csv",
  row.names = FALSE
)

cat("\n===== TOP REACTOME TABLE CHECK =====\n")

cat(
  "All:",
  file.exists(
    "Results/Enrichment/Reactome/Top20_Reactome_All_DEGs.csv"
  ),
  "\n"
)

cat(
  "Up:",
  file.exists(
    "Results/Enrichment/Reactome/Top20_Reactome_Upregulated.csv"
  ),
  "\n"
)
# ============================================================
# STEP 29 - GO Molecular Function Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# 1. All significant DEGs
# ------------------------------------------------------------

ego_mf_all <- enrichGO(
  gene = all_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 2. Upregulated genes
# ------------------------------------------------------------

ego_mf_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 3. Downregulated genes
# ------------------------------------------------------------

ego_mf_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 4. Check results
# ------------------------------------------------------------

cat("\n===== GO MOLECULAR FUNCTION RESULTS =====\n")

cat(
  "GO-MF All pathways:",
  nrow(as.data.frame(ego_mf_all)),
  "\n"
)

cat(
  "GO-MF Up pathways:",
  nrow(as.data.frame(ego_mf_up)),
  "\n"
)

cat(
  "GO-MF Down pathways:",
  nrow(as.data.frame(ego_mf_down)),
  "\n"
)
# ============================================================
# STEP 30 - Save GO Molecular Function Results
# ============================================================

dir.create(
  "Results/Enrichment/GO/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Save CSV tables
# ------------------------------------------------------------

write.csv(
  as.data.frame(ego_mf_all),
  "Results/Enrichment/GO/GO_MF_all_DEGs.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_mf_up),
  "Results/Enrichment/GO/GO_MF_upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_mf_down),
  "Results/Enrichment/GO/GO_MF_downregulated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Save RDS objects
# ------------------------------------------------------------

saveRDS(
  ego_mf_all,
  "Results/Enrichment/GO/GO_MF_all_DEGs.rds"
)

saveRDS(
  ego_mf_up,
  "Results/Enrichment/GO/GO_MF_upregulated.rds"
)

saveRDS(
  ego_mf_down,
  "Results/Enrichment/GO/GO_MF_downregulated.rds"
)

# ------------------------------------------------------------
# Verify
# ------------------------------------------------------------

cat("\n===== GO-MF FILE CHECK =====\n")

cat(
  "All CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_all_DEGs.csv"
  ),
  "\n"
)

cat(
  "Up CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_upregulated.csv"
  ),
  "\n"
)

cat(
  "Down CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_downregulated.csv"
  ),
  "\n"
)

cat(
  "All RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_all_DEGs.rds"
  ),
  "\n"
)

cat(
  "Up RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_upregulated.rds"
  ),
  "\n"
)

cat(
  "Down RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_downregulated.rds"
  ),
  "\n"
)
# ============================================================
# STEP 31 - GO Molecular Function Dotplots
# ============================================================

library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# All significant DEGs
# ------------------------------------------------------------

go_mf_all_plot <- dotplot(
  ego_mf_all,
  showCategory = 20
) +
  ggtitle("GO Molecular Function - All Significant DEGs")

ggsave(
  "Results/Enrichment/GO/Figures/GO_MF_All_DEGs_Dotplot.png",
  go_mf_all_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Upregulated
# ------------------------------------------------------------

go_mf_up_plot <- dotplot(
  ego_mf_up,
  showCategory = 20
) +
  ggtitle("GO Molecular Function - Upregulated Genes")

ggsave(
  "Results/Enrichment/GO/Figures/GO_MF_Upregulated_Dotplot.png",
  go_mf_up_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Downregulated
# ------------------------------------------------------------

if (nrow(as.data.frame(ego_mf_down)) > 0) {
  
  go_mf_down_plot <- dotplot(
    ego_mf_down,
    showCategory = 20
  ) +
    ggtitle("GO Molecular Function - Downregulated Genes")
  
  ggsave(
    "Results/Enrichment/GO/Figures/GO_MF_Downregulated_Dotplot.png",
    go_mf_down_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  
  cat(
    "\nNo significant GO-MF terms for downregulated genes.\n"
  )
}

# ------------------------------------------------------------
# Verify figures
# ------------------------------------------------------------

cat("\n===== GO-MF FIGURE CHECK =====\n")

cat(
  "All:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_MF_All_DEGs_Dotplot.png"
  ),
  "\n"
)

cat(
  "Up:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_MF_Upregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "Down:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_MF_Downregulated_Dotplot.png"
  ),
  "\n"
)

# ============================================================
# STEP 32 - GO Cellular Component Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# 1. All significant DEGs
# ------------------------------------------------------------

ego_cc_all <- enrichGO(
  gene = all_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 2. Upregulated genes
# ------------------------------------------------------------

ego_cc_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 3. Downregulated genes
# ------------------------------------------------------------

ego_cc_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 4. Check results
# ------------------------------------------------------------

cat("\n===== GO CELLULAR COMPONENT RESULTS =====\n")

cat(
  "GO-CC All pathways:",
  nrow(as.data.frame(ego_cc_all)),
  "\n"
)

cat(
  "GO-CC Up pathways:",
  nrow(as.data.frame(ego_cc_up)),
  "\n"
)

cat(
  "GO-CC Down pathways:",
  nrow(as.data.frame(ego_cc_down)),
  "\n"
)
# ============================================================
# STEP 33 - Save GO Cellular Component Results
# ============================================================

# Save CSV tables

write.csv(
  as.data.frame(ego_cc_all),
  "Results/Enrichment/GO/GO_CC_all_DEGs.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_cc_up),
  "Results/Enrichment/GO/GO_CC_upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_cc_down),
  "Results/Enrichment/GO/GO_CC_downregulated.csv",
  row.names = FALSE
)

# Save RDS objects

saveRDS(
  ego_cc_all,
  "Results/Enrichment/GO/GO_CC_all_DEGs.rds"
)

saveRDS(
  ego_cc_up,
  "Results/Enrichment/GO/GO_CC_upregulated.rds"
)

saveRDS(
  ego_cc_down,
  "Results/Enrichment/GO/GO_CC_downregulated.rds"
)

# Verify files

cat("\n===== GO-CC FILE CHECK =====\n")

cat(
  "All CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_all_DEGs.csv"
  ),
  "\n"
)

cat(
  "Up CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_upregulated.csv"
  ),
  "\n"
)

cat(
  "Down CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_downregulated.csv"
  ),
  "\n"
)

cat(
  "All RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_all_DEGs.rds"
  ),
  "\n"
)

cat(
  "Up RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_upregulated.rds"
  ),
  "\n"
)

cat(
  "Down RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_downregulated.rds"
  ),
  "\n"
)
# ============================================================
# STEP 34 - GO Cellular Component Dotplots
# ============================================================

library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# All significant DEGs
# ------------------------------------------------------------

go_cc_all_plot <- dotplot(
  ego_cc_all,
  showCategory = 20
) +
  ggtitle("GO Cellular Component - All Significant DEGs")

ggsave(
  "Results/Enrichment/GO/Figures/GO_CC_All_DEGs_Dotplot.png",
  go_cc_all_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Upregulated
# ------------------------------------------------------------

go_cc_up_plot <- dotplot(
  ego_cc_up,
  showCategory = 20
) +
  ggtitle("GO Cellular Component - Upregulated Genes")

ggsave(
  "Results/Enrichment/GO/Figures/GO_CC_Upregulated_Dotplot.png",
  go_cc_up_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Downregulated
# ------------------------------------------------------------

if (nrow(as.data.frame(ego_cc_down)) > 0) {
  
  go_cc_down_plot <- dotplot(
    ego_cc_down,
    showCategory = 20
  ) +
    ggtitle("GO Cellular Component - Downregulated Genes")
  
  ggsave(
    "Results/Enrichment/GO/Figures/GO_CC_Downregulated_Dotplot.png",
    go_cc_down_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  
  cat(
    "\nNo significant GO-CC terms for downregulated genes.\n"
  )
}

# ------------------------------------------------------------
# Verify figures
# ------------------------------------------------------------

cat("\n===== GO-CC FIGURE CHECK =====\n")

cat(
  "All:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_CC_All_DEGs_Dotplot.png"
  ),
  "\n"
)

cat(
  "Up:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_CC_Upregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "Down:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_CC_Downregulated_Dotplot.png"
  ),
  "\n"
)

# ============================================================
# STEP 35 - Prepare Ranked Gene List for GSEA
# ============================================================

# Make sure DESeq2 results are available
head(res_covid_vs_symptomatic)

# ------------------------------------------------------------
# Create ranked statistics
# ------------------------------------------------------------

gsea_df <- as.data.frame(res_covid_vs_symptomatic)

# Keep genes with valid statistics
gsea_df <- gsea_df[
  !is.na(gsea_df$stat) &
    !is.na(gsea_df$gene_id),
]

# Remove duplicated Ensembl IDs
gsea_df <- gsea_df[
  !duplicated(gsea_df$gene_id),
]

# Create ranked vector
gene_ranking <- gsea_df$stat

names(gene_ranking) <- gsea_df$gene_id

# Sort decreasing
gene_ranking <- sort(
  gene_ranking,
  decreasing = TRUE
)

# ------------------------------------------------------------
# Check ranking
# ------------------------------------------------------------

cat("\n===== GSEA RANKING CHECK =====\n")

cat(
  "Total ranked genes:",
  length(gene_ranking),
  "\n"
)

cat(
  "Highest statistic:",
  max(gene_ranking),
  "\n"
)

cat(
  "Lowest statistic:",
  min(gene_ranking),
  "\n"
)

cat(
  "NA values:",
  sum(is.na(gene_ranking)),
  "\n"
)

cat("\nTop 10 genes:\n")
print(head(gene_ranking, 10))

cat("\nBottom 10 genes:\n")
print(tail(gene_ranking, 10))
# ============================================================
# STEP 36 - Convert GSEA Ranking: Ensembl -> Entrez
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Convert Ensembl IDs to Entrez IDs
# ------------------------------------------------------------

gene_map_gsea <- bitr(
  names(gene_ranking),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# ------------------------------------------------------------
# Remove duplicated mappings
# Keep one Ensembl gene per Entrez ID
# ------------------------------------------------------------

gene_map_gsea <- gene_map_gsea[
  !duplicated(gene_map_gsea$ENTREZID),
]

# ------------------------------------------------------------
# Match ranking statistics to mapped genes
# ------------------------------------------------------------

gene_map_gsea$stat <- gene_ranking[
  match(
    gene_map_gsea$ENSEMBL,
    names(gene_ranking)
  )
]

# ------------------------------------------------------------
# Create Entrez-ranked vector
# ------------------------------------------------------------

gene_ranking_entrez <- gene_map_gsea$stat

names(gene_ranking_entrez) <- gene_map_gsea$ENTREZID

gene_ranking_entrez <- sort(
  gene_ranking_entrez,
  decreasing = TRUE
)

# ------------------------------------------------------------
# Check
# ------------------------------------------------------------

cat("\n===== ENTREZ GSEA RANKING CHECK =====\n")

cat(
  "Original Ensembl genes:",
  length(gene_ranking),
  "\n"
)

cat(
  "Mapped Entrez genes:",
  length(gene_ranking_entrez),
  "\n"
)

cat(
  "Mapping percentage:",
  round(
    100 * length(gene_ranking_entrez) /
      length(gene_ranking),
    2
  ),
  "%\n"
)

cat(
  "NA values:",
  sum(is.na(gene_ranking_entrez)),
  "\n"
)

cat(
  "Duplicated Entrez IDs:",
  sum(duplicated(names(gene_ranking_entrez))),
  "\n"
)

cat("\nTop 10 Entrez genes:\n")
print(head(gene_ranking_entrez, 10))

cat("\nBottom 10 Entrez genes:\n")
print(tail(gene_ranking_entrez, 10))

# ============================================================
# STEP 37 - GSEA: GO Biological Process
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# Run GO-BP GSEA
# ------------------------------------------------------------

gsea_go_bp <- gseGO(
  geneList = gene_ranking_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE
)

# ------------------------------------------------------------
# Check result
# ------------------------------------------------------------

cat("\n===== GO-BP GSEA RESULT =====\n")

cat(
  "Number of enriched pathways:",
  nrow(as.data.frame(gsea_go_bp)),
  "\n"
)

if (nrow(as.data.frame(gsea_go_bp)) > 0) {
  
  cat("\nTop 10 enriched GO-BP pathways:\n")
  
  print(
    as.data.frame(gsea_go_bp)[
      1:min(10, nrow(as.data.frame(gsea_go_bp))),
      c(
        "ID",
        "Description",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
} else {
  
  cat("\nNo significant GO-BP pathways detected.\n")
}

# ============================================================
# STEP 38 - Save GO-BP GSEA Results + Figures
# ============================================================

# Create output directories
dir.create(
  "Results/Enrichment/GSEA",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/Enrichment/GSEA/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Convert result to data frame
# ------------------------------------------------------------

gsea_go_bp_df <- as.data.frame(gsea_go_bp)

# ------------------------------------------------------------
# Save complete results
# ------------------------------------------------------------

write.csv(
  gsea_go_bp_df,
  "Results/Enrichment/GSEA/GO_BP_GSEA.csv",
  row.names = FALSE
)

saveRDS(
  gsea_go_bp,
  "Results/Enrichment/GSEA/GO_BP_GSEA.rds"
)

# ------------------------------------------------------------
# Save top 20 pathways
# ------------------------------------------------------------

top20_gsea_go_bp <- gsea_go_bp_df[
  order(gsea_go_bp_df$p.adjust),
]

top20_gsea_go_bp <- head(
  top20_gsea_go_bp,
  20
)

write.csv(
  top20_gsea_go_bp,
  "Results/Enrichment/GSEA/Top20_GO_BP_GSEA.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Dotplot
# ------------------------------------------------------------

p_gsea_go_bp <- dotplot(
  gsea_go_bp,
  showCategory = 20,
  title = "GO Biological Process GSEA"
)

ggsave(
  "Results/Enrichment/GSEA/Figures/GO_BP_GSEA_Dotplot.png",
  p_gsea_go_bp,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Enrichment plots for selected top pathways
# ------------------------------------------------------------

top_pathways <- gsea_go_bp_df$ID[
  order(gsea_go_bp_df$p.adjust)
]

top_pathways <- head(
  top_pathways,
  3
)

for (i in seq_along(top_pathways)) {
  
  pathway_id <- top_pathways[i]
  
  p <- gseaplot2(
    gsea_go_bp,
    geneSetID = pathway_id,
    title = gsea_go_bp_df$Description[
      match(pathway_id, gsea_go_bp_df$ID)
    ]
  )
  
  ggsave(
    paste0(
      "Results/Enrichment/GSEA/Figures/",
      "GO_BP_GSEA_",
      i,
      "_",
      pathway_id,
      ".png"
    ),
    p,
    width = 10,
    height = 7,
    dpi = 300
  )
}

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

cat("\n===== GSEA SAVE CHECK =====\n")

cat(
  "GO-BP GSEA CSV exists:",
  file.exists(
    "Results/Enrichment/GSEA/GO_BP_GSEA.csv"
  ),
  "\n"
)

cat(
  "GO-BP GSEA RDS exists:",
  file.exists(
    "Results/Enrichment/GSEA/GO_BP_GSEA.rds"
  ),
  "\n"
)

cat(
  "GO-BP GSEA Dotplot exists:",
  file.exists(
    "Results/Enrichment/GSEA/Figures/GO_BP_GSEA_Dotplot.png"
  ),
  "\n"
)

cat(
  "Number of GSEA pathways saved:",
  nrow(gsea_go_bp_df),
  "\n"
)

# ============================================================
# STEP 39 - Reactome GSEA
# ============================================================

library(ReactomePA)
library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# Run Reactome GSEA
# ------------------------------------------------------------

gsea_reactome <- gsePathway(
  geneList = gene_ranking_entrez,
  organism = "human",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE
)

# ------------------------------------------------------------
# Check result
# ------------------------------------------------------------

gsea_reactome_df <- as.data.frame(gsea_reactome)

cat("\n===== REACTOME GSEA RESULT =====\n")

cat(
  "Number of enriched Reactome pathways:",
  nrow(gsea_reactome_df),
  "\n"
)

if (nrow(gsea_reactome_df) > 0) {
  
  cat("\nTop 10 Reactome pathways:\n")
  
  print(
    gsea_reactome_df[
      1:min(10, nrow(gsea_reactome_df)),
      c(
        "ID",
        "Description",
        "setSize",
        "enrichmentScore",
        "NES",
        "pvalue",
        "p.adjust"
      )
    ]
  )
  
} else {
  
  cat("\nNo significant Reactome pathways detected.\n")
}

# ============================================================
# STEP 40 - Save Reactome GSEA Results + Figures
# ============================================================

# ------------------------------------------------------------
# Save complete Reactome GSEA results
# ------------------------------------------------------------

write.csv(
  gsea_reactome_df,
  "Results/Enrichment/GSEA/Reactome_GSEA.csv",
  row.names = FALSE
)

saveRDS(
  gsea_reactome,
  "Results/Enrichment/GSEA/Reactome_GSEA.rds"
)

# ------------------------------------------------------------
# Top 20 Reactome pathways
# ------------------------------------------------------------

top20_gsea_reactome <- gsea_reactome_df[
  order(gsea_reactome_df$p.adjust),
]

top20_gsea_reactome <- head(
  top20_gsea_reactome,
  20
)

write.csv(
  top20_gsea_reactome,
  "Results/Enrichment/GSEA/Top20_Reactome_GSEA.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Reactome GSEA Dotplot
# ------------------------------------------------------------

p_gsea_reactome <- dotplot(
  gsea_reactome,
  showCategory = 20,
  title = "Reactome GSEA"
)

ggsave(
  "Results/Enrichment/GSEA/Figures/Reactome_GSEA_Dotplot.png",
  p_gsea_reactome,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

cat("\n===== REACTOME GSEA SAVE CHECK =====\n")

cat(
  "Reactome GSEA CSV exists:",
  file.exists(
    "Results/Enrichment/GSEA/Reactome_GSEA.csv"
  ),
  "\n"
)

cat(
  "Reactome GSEA RDS exists:",
  file.exists(
    "Results/Enrichment/GSEA/Reactome_GSEA.rds"
  ),
  "\n"
)

cat(
  "Reactome GSEA Dotplot exists:",
  file.exists(
    "Results/Enrichment/GSEA/Figures/Reactome_GSEA_Dotplot.png"
  ),
  "\n"
)

cat(
  "Number of Reactome pathways saved:",
  nrow(gsea_reactome_df),
  "\n"
)

# ============================================================
# STEP 41-A - Restore Previous DE Results
# ============================================================

res_covid_vs_symptomatic <- readRDS(
  "Results/Differential_Expression/COVID_vs_Symptomatic_all_genes.rds"
)

sig_degs <- readRDS(
  "Results/Differential_Expression/COVID_vs_Symptomatic_significant_DEGs.rds"
)

upregulated_genes <- read.csv(
  "Results/Differential_Expression/COVID_vs_Symptomatic_upregulated.csv",
  stringsAsFactors = FALSE
)

downregulated_genes <- read.csv(
  "Results/Differential_Expression/COVID_vs_Symptomatic_downregulated.csv",
  stringsAsFactors = FALSE
)

cat("\n===== DE RESULTS RESTORED =====\n")

cat(
  "Complete DESeq2 results:",
  nrow(res_covid_vs_symptomatic),
  "\n"
)

cat(
  "Significant DEGs:",
  nrow(sig_degs),
  "\n"
)

cat(
  "Upregulated:",
  nrow(upregulated_genes),
  "\n"
)

cat(
  "Downregulated:",
  nrow(downregulated_genes),
  "\n"
)
# ============================================================
# STEP 41-B - Restore GSEA Ranked Gene Lists
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Create Ensembl-ranked list
# ------------------------------------------------------------

gsea_df <- as.data.frame(res_covid_vs_symptomatic)

gsea_df <- gsea_df[
  !is.na(gsea_df$stat) &
    !is.na(gsea_df$gene_id),
]

gsea_df <- gsea_df[
  !duplicated(gsea_df$gene_id),
]

gene_ranking <- gsea_df$stat

names(gene_ranking) <- gsea_df$gene_id

gene_ranking <- sort(
  gene_ranking,
  decreasing = TRUE
)

# ------------------------------------------------------------
# Convert Ensembl -> Entrez
# ------------------------------------------------------------

gene_map_gsea <- bitr(
  names(gene_ranking),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Remove duplicate Entrez IDs
gene_map_gsea <- gene_map_gsea[
  !duplicated(gene_map_gsea$ENTREZID),
]

# Match statistics
gene_map_gsea$stat <- gene_ranking[
  match(
    gene_map_gsea$ENSEMBL,
    names(gene_ranking)
  )
]

# Create Entrez-ranked vector
gene_ranking_entrez <- gene_map_gsea$stat

names(gene_ranking_entrez) <- gene_map_gsea$ENTREZID

gene_ranking_entrez <- sort(
  gene_ranking_entrez,
  decreasing = TRUE
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

cat("\n===== GSEA RANKING RESTORED =====\n")

cat(
  "Ensembl ranked genes:",
  length(gene_ranking),
  "\n"
)

cat(
  "Entrez ranked genes:",
  length(gene_ranking_entrez),
  "\n"
)

cat(
  "Highest statistic:",
  max(gene_ranking_entrez),
  "\n"
)

cat(
  "Lowest statistic:",
  min(gene_ranking_entrez),
  "\n"
)

cat(
  "NA values:",
  sum(is.na(gene_ranking_entrez)),
  "\n"
)

cat(
  "Duplicated Entrez IDs:",
  sum(duplicated(names(gene_ranking_entrez))),
  "\n"
)
# ============================================================
# STEP 41 - Leading-edge Genes
# ============================================================

# ------------------------------------------------------------
# Select top Reactome pathways
# ------------------------------------------------------------

selected_reactome <- c(
  "R-HSA-69481",
  "R-HSA-69620",
  "R-HSA-69306"
)

reactome_selected <- gsea_reactome_df[
  gsea_reactome_df$ID %in% selected_reactome,
]

reactome_selected <- reactome_selected[
  order(reactome_selected$p.adjust),
]

# ------------------------------------------------------------
# Display pathway statistics
# ------------------------------------------------------------

cat("\n===== SELECTED REACTOME PATHWAYS =====\n\n")

print(
  reactome_selected[
    ,
    c(
      "ID",
      "Description",
      "setSize",
      "enrichmentScore",
      "NES",
      "pvalue",
      "p.adjust"
    )
  ]
)

# ------------------------------------------------------------
# Extract core enrichment genes
# ------------------------------------------------------------

leading_edge_list <- list()

for (i in seq_len(nrow(reactome_selected))) {
  
  pathway_id <- reactome_selected$ID[i]
  
  pathway_name <- reactome_selected$Description[i]
  
  genes <- unlist(
    strsplit(
      reactome_selected$core_enrichment[i],
      "/"
    )
  )
  
  leading_edge_list[[pathway_id]] <- genes
  
  cat("\n============================================\n")
  cat("Pathway:", pathway_name, "\n")
  cat("Pathway ID:", pathway_id, "\n")
  cat("Leading-edge genes:", length(genes), "\n")
  cat("============================================\n")
  
  print(genes)
}

# ============================================================
# STEP 42 — Convert Leading-edge Entrez IDs to Gene Symbols
# ============================================================

library(org.Hs.eg.db)
library(AnnotationDbi)

# ------------------------------------------------------------
# Combine all leading-edge genes
# ------------------------------------------------------------

leading_edge_df <- data.frame(
  EntrezID = unique(unlist(leading_edge_list)),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Map Entrez IDs -> Gene Symbols
# ------------------------------------------------------------

gene_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = leading_edge_df$EntrezID,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "ENTREZID"
)

gene_mapping <- gene_mapping[
  !duplicated(gene_mapping$ENTREZID),
]

# ------------------------------------------------------------
# Merge mapping
# ------------------------------------------------------------

leading_edge_df <- merge(
  leading_edge_df,
  gene_mapping,
  by.x = "EntrezID",
  by.y = "ENTREZID",
  all.x = TRUE
)

# ------------------------------------------------------------
# Add pathway membership
# ------------------------------------------------------------

leading_edge_df$G2M_Checkpoints <- leading_edge_df$EntrezID %in%
  leading_edge_list[["R-HSA-69481"]]

leading_edge_df$Cell_Cycle_Checkpoints <- leading_edge_df$EntrezID %in%
  leading_edge_list[["R-HSA-69620"]]

leading_edge_df$DNA_Replication <- leading_edge_df$EntrezID %in%
  leading_edge_list[["R-HSA-69306"]]

# Number of pathways containing each gene
leading_edge_df$N_Pathways <- rowSums(
  leading_edge_df[
    ,
    c(
      "G2M_Checkpoints",
      "Cell_Cycle_Checkpoints",
      "DNA_Replication"
    )
  ]
)

# ------------------------------------------------------------
# Sort: genes shared by all pathways first
# ------------------------------------------------------------

leading_edge_df <- leading_edge_df[
  order(
    -leading_edge_df$N_Pathways,
    leading_edge_df$SYMBOL
  ),
]

# ------------------------------------------------------------
# Display summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 42 — LEADING-EDGE GENE SUMMARY\n")
cat("============================================\n\n")

cat("Unique leading-edge Entrez IDs:",
    nrow(leading_edge_df), "\n")

cat(
  "Mapped Gene Symbols:",
  sum(!is.na(leading_edge_df$SYMBOL)),
  "\n"
)

cat("\nGenes shared by ALL 3 pathways:\n")

shared_all <- leading_edge_df$SYMBOL[
  leading_edge_df$N_Pathways == 3 &
    !is.na(leading_edge_df$SYMBOL)
]

print(shared_all)

cat("\nNumber shared by all 3:",
    length(shared_all), "\n")

cat("\nTop leading-edge genes by pathway overlap:\n")

print(
  leading_edge_df[
    ,
    c(
      "EntrezID",
      "SYMBOL",
      "G2M_Checkpoints",
      "Cell_Cycle_Checkpoints",
      "DNA_Replication",
      "N_Pathways"
    )
  ],
  row.names = FALSE
)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

dir.create(
  "Results/Enrichment/GSEA/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  leading_edge_df,
  "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.csv",
  row.names = FALSE
)

write.csv(
  data.frame(GeneSymbol = shared_all),
  "Results/Enrichment/GSEA/Shared_Leading_Edge_Genes_Top3.csv",
  row.names = FALSE
)

saveRDS(
  leading_edge_df,
  "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.rds"
)

cat("\n============================================\n")
cat("FILES SAVED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "Leading-edge table: ",
  file.exists(
    "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.csv"
  ),
  "\n"
)

cat(
  "Shared genes table: ",
  file.exists(
    "Results/Enrichment/GSEA/Shared_Leading_Edge_Genes_Top3.csv"
  ),
  "\n"
)

# ============================================================
# STEP 43 — Core Leading-edge Genes + DESeq2 Statistics
# ============================================================

core_genes <- shared_all

# ------------------------------------------------------------
# Extract DESeq2 results
# ------------------------------------------------------------

core_de <- res_covid_vs_symptomatic[
  res_covid_vs_symptomatic$SYMBOL %in% core_genes |
    rownames(res_covid_vs_symptomatic) %in% core_genes,
]

# ------------------------------------------------------------
# If SYMBOL column is unavailable, map rownames
# ------------------------------------------------------------

if (nrow(core_de) == 0) {
  
  core_entrez <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = core_genes,
    columns = c("SYMBOL"),
    keytype = "SYMBOL"
  )
  
  core_de <- res_covid_vs_symptomatic[
    rownames(res_covid_vs_symptomatic) %in%
      core_entrez$SYMBOL,
  ]
}

# ------------------------------------------------------------
# Add Gene Symbol from rownames when necessary
# ------------------------------------------------------------

if (!"SYMBOL" %in% colnames(core_de)) {
  core_de$SYMBOL <- rownames(core_de)
}

# ------------------------------------------------------------
# Select useful columns
# ------------------------------------------------------------

core_de <- core_de[
  ,
  intersect(
    c(
      "SYMBOL",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    ),
    colnames(core_de)
  ),
  drop = FALSE
]

# ------------------------------------------------------------
# Remove duplicates
# ------------------------------------------------------------

core_de <- core_de[
  !duplicated(core_de$SYMBOL),
]

# ------------------------------------------------------------
# Add regulation status
# ------------------------------------------------------------

core_de$Regulation <- ifelse(
  !is.na(core_de$padj) &
    core_de$padj < 0.05 &
    core_de$log2FoldChange > 1,
  "Upregulated",
  ifelse(
    !is.na(core_de$padj) &
      core_de$padj < 0.05 &
      core_de$log2FoldChange < -1,
    "Downregulated",
    "Not significant"
  )
)

# ------------------------------------------------------------
# Sort by log2FC
# ------------------------------------------------------------

core_de <- core_de[
  order(
    -core_de$log2FoldChange
  ),
]

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 43 — CORE LEADING-EDGE + DESEQ2\n")
cat("============================================\n\n")

cat(
  "Core genes:",
  length(core_genes),
  "\n"
)

cat(
  "Genes found in DESeq2:",
  nrow(core_de),
  "\n"
)

cat("\nRegulation summary:\n")
print(table(core_de$Regulation))

cat("\nCore genes with DESeq2 statistics:\n\n")

print(
  core_de,
  row.names = FALSE
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  core_de,
  "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.csv",
  row.names = FALSE
)

saveRDS(
  core_de,
  "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.rds"
)

cat("\n============================================\n")
cat("FILES SAVED\n")
cat("============================================\n")

cat(
  "CSV:",
  file.exists(
    "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.csv"
  ),
  "\n"
)

cat(
  "RDS:",
  file.exists(
    "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.rds"
  ),
  "\n"
)
# ============================================================
# STEP 43 FIX — Map DESeq2 Ensembl IDs to Gene Symbols
# ============================================================

library(org.Hs.eg.db)
library(AnnotationDbi)

# ------------------------------------------------------------
# Check DESeq2 structure
# ------------------------------------------------------------

cat("\n============================================\n")
cat("DESeq2 RESULT STRUCTURE\n")
cat("============================================\n\n")

cat("Number of DESeq2 genes:",
    nrow(res_covid_vs_symptomatic), "\n\n")

cat("DESeq2 columns:\n")
print(colnames(res_covid_vs_symptomatic))

cat("\nFirst 5 rownames:\n")
print(head(rownames(res_covid_vs_symptomatic)))

# ------------------------------------------------------------
# Create Ensembl ID column
# ------------------------------------------------------------

de_ids <- rownames(res_covid_vs_symptomatic)

# Remove version numbers if present
# Example: ENSG00000123456.12 -> ENSG00000123456

de_ensembl <- sub(
  "\\..*$",
  "",
  de_ids
)

# ------------------------------------------------------------
# Map Ensembl -> Gene Symbol
# ------------------------------------------------------------

de_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(de_ensembl),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
)

# Keep one mapping per Ensembl ID
de_mapping <- de_mapping[
  !duplicated(de_mapping$ENSEMBL),
]

# ------------------------------------------------------------
# Add Ensembl + Symbol
# ------------------------------------------------------------

de_table <- as.data.frame(res_covid_vs_symptomatic)

de_table$ENSEMBL <- de_ensembl

de_table$SYMBOL <- de_mapping$SYMBOL[
  match(
    de_table$ENSEMBL,
    de_mapping$ENSEMBL
  )
]

# ------------------------------------------------------------
# Extract the 54 core genes
# ------------------------------------------------------------

core_de <- de_table[
  !is.na(de_table$SYMBOL) &
    de_table$SYMBOL %in% core_genes,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Remove duplicate symbols
# ------------------------------------------------------------

core_de <- core_de[
  !duplicated(core_de$SYMBOL),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Add regulation status
# ------------------------------------------------------------

core_de$Regulation <- ifelse(
  !is.na(core_de$padj) &
    core_de$padj < 0.05 &
    core_de$log2FoldChange > 1,
  "Upregulated",
  ifelse(
    !is.na(core_de$padj) &
      core_de$padj < 0.05 &
      core_de$log2FoldChange < -1,
    "Downregulated",
    "Not significant"
  )
)

# ------------------------------------------------------------
# Sort by log2FC
# ------------------------------------------------------------

core_de <- core_de[
  order(
    -core_de$log2FoldChange
  ),
]

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 43 FIX — CORE GENE RESULTS\n")
cat("============================================\n\n")

cat(
  "Core genes from GSEA:",
  length(core_genes),
  "\n"
)

cat(
  "Core genes found in DESeq2:",
  nrow(core_de),
  "\n"
)

cat("\nRegulation summary:\n")
print(table(core_de$Regulation))

cat("\nCore genes + DESeq2 statistics:\n\n")

print(
  core_de[
    ,
    c(
      "ENSEMBL",
      "SYMBOL",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj",
      "Regulation"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

# ------------------------------------------------------------
# Save corrected results
# ------------------------------------------------------------

write.csv(
  core_de,
  "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.csv",
  row.names = FALSE
)

saveRDS(
  core_de,
  "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.rds"
)

cat("\n============================================\n")
cat("FILES SAVED\n")
cat("============================================\n")

cat(
  "CSV:",
  file.exists(
    "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.csv"
  ),
  "\n"
)

cat(
  "RDS:",
  file.exists(
    "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.rds"
  ),
  "\n"
)
# ============================================================
# STEP 44 — Core Cell-Cycle Signature
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# Select significant core genes
# ------------------------------------------------------------

core_sig <- core_de[
  core_de$Regulation == "Upregulated",
  ,
  drop = FALSE
]

# Sort by log2FC
core_sig <- core_sig[
  order(core_sig$log2FoldChange),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Print
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 44 — CORE CELL-CYCLE SIGNATURE\n")
cat("============================================\n\n")

cat(
  "Significant core genes:",
  nrow(core_sig),
  "\n\n"
)

print(
  core_sig[
    ,
    c(
      "SYMBOL",
      "baseMean",
      "log2FoldChange",
      "pvalue",
      "padj"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

# ------------------------------------------------------------
# Save table
# ------------------------------------------------------------

write.csv(
  core_sig,
  "Results/Enrichment/GSEA/Core_Significant_Genes.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Create plot
# ------------------------------------------------------------

p_core_signature <- ggplot(
  core_sig,
  aes(
    x = log2FoldChange,
    y = SYMBOL
  )
) +
  geom_point(
    size = 4
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  labs(
    title = "Core Cell-Cycle / DNA-Replication Signature",
    subtitle = "Significant leading-edge genes from top Reactome GSEA pathways",
    x = "log2 Fold Change",
    y = "Gene"
  ) +
  theme_minimal(base_size = 13)

print(p_core_signature)

# ------------------------------------------------------------
# Save plot
# ------------------------------------------------------------

ggsave(
  "Results/Enrichment/GSEA/Figures/Core_Cell_Cycle_Signature.png",
  p_core_signature,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FILES SAVED\n")
cat("============================================\n")

cat(
  "Table:",
  file.exists(
    "Results/Enrichment/GSEA/Core_Significant_Genes.csv"
  ),
  "\n"
)

cat(
  "Plot:",
  file.exists(
    "Results/Enrichment/GSEA/Figures/Core_Cell_Cycle_Signature.png"
  ),
  "\n"
)
# ============================================================
# STEP 45 — Core Gene–Pathway Network
# ============================================================

cat("\n============================================================\n")
cat("STEP 45 — CORE GENE–PATHWAY NETWORK\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load required objects
# ------------------------------------------------------------

core_sig <- read.csv(
  "Results/Enrichment/GSEA/Core_Significant_Genes.csv",
  stringsAsFactors = FALSE
)

leading_edge <- read.csv(
  "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Check column names
# ------------------------------------------------------------

cat("Core signature columns:\n")
print(colnames(core_sig))

cat("\nLeading-edge columns:\n")
print(colnames(leading_edge))

# ------------------------------------------------------------
# 3. Keep only the 12 significant core genes
# ------------------------------------------------------------

core_genes <- unique(core_sig$SYMBOL)

cat("\nNumber of core genes:", length(core_genes), "\n")
print(core_genes)

# ------------------------------------------------------------
# 4. Identify pathway column
# ------------------------------------------------------------

pathway_col <- intersect(
  c("Description", "Pathway", "pathway", "ID"),
  colnames(leading_edge)
)

gene_col <- intersect(
  c("SYMBOL", "gene_symbol", "Gene", "gene"),
  colnames(leading_edge)
)

cat("\nDetected pathway column:", pathway_col[1], "\n")
cat("Detected gene column:", gene_col[1], "\n")

# ------------------------------------------------------------
# 5. Build gene–pathway network
# ------------------------------------------------------------

network_df <- leading_edge[
  leading_edge[[gene_col[1]]] %in% core_genes,
  c(pathway_col[1], gene_col[1])
]

colnames(network_df) <- c(
  "Pathway",
  "Gene"
)

network_df <- unique(network_df)

# ------------------------------------------------------------
# 6. Add gene regulation information
# ------------------------------------------------------------

network_df <- merge(
  network_df,
  core_sig[, c("SYMBOL", "log2FoldChange", "padj", "Regulation")],
  by.x = "Gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

network_df <- network_df[
  order(network_df$Pathway, -network_df$log2FoldChange),
]

# ------------------------------------------------------------
# 7. Print network
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE–PATHWAY NETWORK\n")
cat("============================================================\n\n")

print(network_df)

# ------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------

cat("\nNumber of network edges:", nrow(network_df), "\n")
cat("Number of genes:", length(unique(network_df$Gene)), "\n")
cat("Number of pathways:", length(unique(network_df$Pathway)), "\n")

# ------------------------------------------------------------
# 9. Save network table
# ------------------------------------------------------------

dir.create(
  "Results/Network",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  network_df,
  "Results/Network/Core_Gene_Pathway_Network.csv",
  row.names = FALSE
)

saveRDS(
  network_df,
  "Results/Network/Core_Gene_Pathway_Network.rds"
)

# ------------------------------------------------------------
# 10. Validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n")

cat(
  "CSV:",
  file.exists(
    "Results/Network/Core_Gene_Pathway_Network.csv"
  ),
  "\n"
)

cat(
  "RDS:",
  file.exists(
    "Results/Network/Core_Gene_Pathway_Network.rds"
  ),
  "\n"
)

cat("\nSTEP 45 completed successfully.\n")
# ============================================================
# STEP 45 — FIX
# Core Gene–Pathway Network
# ============================================================

cat("\n============================================================\n")
cat("STEP 45 — CORE GENE–PATHWAY NETWORK\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

core_sig <- read.csv(
  "Results/Enrichment/GSEA/Core_Significant_Genes.csv",
  stringsAsFactors = FALSE
)

leading_edge <- read.csv(
  "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Core genes
# ------------------------------------------------------------

core_genes <- unique(core_sig$SYMBOL)

cat("Number of core genes:", length(core_genes), "\n")
print(core_genes)

# ------------------------------------------------------------
# 3. Reactome pathway columns
# ------------------------------------------------------------

pathway_cols <- c(
  "G2M_Checkpoints",
  "Cell_Cycle_Checkpoints",
  "DNA_Replication"
)

cat("\nPathway columns:\n")
print(pathway_cols)

# Check that all columns exist
cat("\nAll pathway columns exist:",
    all(pathway_cols %in% colnames(leading_edge)),
    "\n")

# ------------------------------------------------------------
# 4. Create gene–pathway edges
# ------------------------------------------------------------

network_list <- lapply(pathway_cols, function(pw) {
  
  x <- leading_edge[
    leading_edge$SYMBOL %in% core_genes &
      !is.na(leading_edge[[pw]]) &
      leading_edge[[pw]] != FALSE &
      leading_edge[[pw]] != 0,
    c("SYMBOL", pw)
  ]
  
  if (nrow(x) == 0) {
    return(NULL)
  }
  
  data.frame(
    Gene = x$SYMBOL,
    Pathway = pw,
    stringsAsFactors = FALSE
  )
})

network_df <- do.call(rbind, network_list)

# Remove duplicates
network_df <- unique(network_df)

# ------------------------------------------------------------
# 5. Add DESeq2 information
# ------------------------------------------------------------

network_df <- merge(
  network_df,
  core_sig[, c(
    "SYMBOL",
    "log2FoldChange",
    "padj",
    "Regulation"
  )],
  by.x = "Gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

# Sort
network_df <- network_df[
  order(network_df$Pathway, -network_df$log2FoldChange),
]

# ------------------------------------------------------------
# 6. Clean pathway names
# ------------------------------------------------------------

network_df$Pathway <- factor(
  network_df$Pathway,
  levels = pathway_cols
)

network_df$Pathway <- as.character(network_df$Pathway)

network_df$Pathway[
  network_df$Pathway == "G2M_Checkpoints"
] <- "G2/M Checkpoints"

network_df$Pathway[
  network_df$Pathway == "Cell_Cycle_Checkpoints"
] <- "Cell Cycle Checkpoints"

network_df$Pathway[
  network_df$Pathway == "DNA_Replication"
] <- "DNA Replication"

# ------------------------------------------------------------
# 7. Print network
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE–PATHWAY NETWORK\n")
cat("============================================================\n\n")

print(network_df)

# ------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("NETWORK SUMMARY\n")
cat("============================================================\n\n")

cat("Number of edges:",
    nrow(network_df),
    "\n")

cat("Number of unique genes:",
    length(unique(network_df$Gene)),
    "\n")

cat("Number of pathways:",
    length(unique(network_df$Pathway)),
    "\n")

cat("\nEdges per pathway:\n")
print(table(network_df$Pathway))

# ------------------------------------------------------------
# 9. Save results
# ------------------------------------------------------------

dir.create(
  "Results/Network",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  network_df,
  "Results/Network/Core_Gene_Pathway_Network.csv",
  row.names = FALSE
)

saveRDS(
  network_df,
  "Results/Network/Core_Gene_Pathway_Network.rds"
)

# ------------------------------------------------------------
# 10. Validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n")

cat(
  "CSV:",
  file.exists(
    "Results/Network/Core_Gene_Pathway_Network.csv"
  ),
  "\n"
)

cat(
  "RDS:",
  file.exists(
    "Results/Network/Core_Gene_Pathway_Network.rds"
  ),
  "\n"
)

cat("\nSTEP 45 FIX completed successfully.\n")


# ============================================================
# STEP 46 — PREPARE CYTOSCAPE NETWORK FILES
# ============================================================

cat("\n============================================================\n")
cat("STEP 46 — CYTOSCAPE NETWORK PREPARATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load network
# ------------------------------------------------------------

network_df <- readRDS(
  "Results/Network/Core_Gene_Pathway_Network.rds"
)

# ------------------------------------------------------------
# 2. Create EDGE table
# ------------------------------------------------------------

edges <- network_df[, c(
  "Gene",
  "Pathway"
)]

colnames(edges) <- c(
  "source",
  "target"
)

edges$interaction <- "core_gene_in_pathway"

# ------------------------------------------------------------
# 3. Create GENE nodes
# ------------------------------------------------------------

gene_nodes <- unique(
  network_df[, c(
    "Gene",
    "log2FoldChange",
    "padj",
    "Regulation"
  )]
)

gene_nodes$NodeType <- "Gene"

colnames(gene_nodes)[1] <- "Node"

# ------------------------------------------------------------
# 4. Create PATHWAY nodes
# ------------------------------------------------------------

pathway_nodes <- data.frame(
  Node = unique(network_df$Pathway),
  log2FoldChange = NA,
  padj = NA,
  Regulation = NA,
  NodeType = "Pathway",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 5. Combine nodes
# ------------------------------------------------------------

nodes <- rbind(
  gene_nodes,
  pathway_nodes
)

# Remove duplicate nodes if any
nodes <- unique(nodes)

# ------------------------------------------------------------
# 6. Save Cytoscape files
# ------------------------------------------------------------

dir.create(
  "Results/Network/Cytoscape",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  edges,
  "Results/Network/Cytoscape/Core_Network_Edges.csv",
  row.names = FALSE
)

write.csv(
  nodes,
  "Results/Network/Cytoscape/Core_Network_Nodes.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n\n")

cat("Number of edges:", nrow(edges), "\n")
cat("Number of nodes:", nrow(nodes), "\n")

cat("\nNode types:\n")
print(table(nodes$NodeType))

cat("\nEdges file:",
    file.exists(
      "Results/Network/Cytoscape/Core_Network_Edges.csv"
    ),
    "\n")

cat("Nodes file:",
    file.exists(
      "Results/Network/Cytoscape/Core_Network_Nodes.csv"
    ),
    "\n")

cat("\nSTEP 46 completed successfully.\n")

# ============================================================
# STEP 47 — FIX
# CORE GENE–PATHWAY NETWORK VISUALIZATION
# ============================================================

cat("\n============================================================\n")
cat("STEP 47 — NETWORK VISUALIZATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(igraph)
library(ggraph)
library(ggplot2)

# ------------------------------------------------------------
# 2. Load network files
# ------------------------------------------------------------

edges <- read.csv(
  "Results/Network/Cytoscape/Core_Network_Edges.csv",
  stringsAsFactors = FALSE
)

nodes <- read.csv(
  "Results/Network/Cytoscape/Core_Network_Nodes.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 3. Create TRUE bipartite node type
# ------------------------------------------------------------

nodes$type <- nodes$NodeType == "Pathway"

# ------------------------------------------------------------
# 4. Create igraph object
# ------------------------------------------------------------

graph <- graph_from_data_frame(
  d = edges[, c("source", "target")],
  vertices = nodes,
  directed = FALSE
)

# ------------------------------------------------------------
# 5. Add readable node labels
# ------------------------------------------------------------

V(graph)$type_label <- ifelse(
  V(graph)$type,
  "Pathway",
  "Gene"
)

# ------------------------------------------------------------
# 6. Add node sizes
# ------------------------------------------------------------

logfc_values <- nodes$log2FoldChange[
  match(V(graph)$name, nodes$Node)
]

V(graph)$node_size <- ifelse(
  V(graph)$type_label == "Gene",
  4 + abs(logfc_values) * 4,
  9
)

# ------------------------------------------------------------
# 7. Create network plot
# ------------------------------------------------------------

network_plot <- ggraph(
  graph,
  layout = "bipartite"
) +
  geom_edge_link(
    linewidth = 0.7,
    alpha = 0.6
  ) +
  geom_node_point(
    aes(
      size = node_size,
      shape = type_label
    )
  ) +
  geom_node_text(
    aes(
      label = name
    ),
    repel = TRUE,
    size = 4
  ) +
  scale_shape_manual(
    values = c(
      Gene = 16,
      Pathway = 17
    )
  ) +
  scale_size_identity() +
  labs(
    title = "Core Gene–Pathway Network",
    subtitle = "12 significant core genes connected to top Reactome pathways",
    shape = "Node type"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 12
    ),
    legend.position = "bottom"
  )

# ------------------------------------------------------------
# 8. Display
# ------------------------------------------------------------

print(network_plot)

# ------------------------------------------------------------
# 9. Save figure
# ------------------------------------------------------------

dir.create(
  "Results/Network/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "Results/Network/Figures/Core_Gene_Pathway_Network.png",
  network_plot,
  width = 12,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n\n")

cat("Nodes:", vcount(graph), "\n")
cat("Edges:", ecount(graph), "\n")

cat(
  "Gene nodes:",
  sum(V(graph)$type_label == "Gene"),
  "\n"
)

cat(
  "Pathway nodes:",
  sum(V(graph)$type_label == "Pathway"),
  "\n"
)

cat(
  "Plot:",
  file.exists(
    "Results/Network/Figures/Core_Gene_Pathway_Network.png"
  ),
  "\n"
)

cat("\nSTEP 47 FIX completed successfully.\n")
# ============================================================
# STEP 48 — IMPROVED CORE NETWORK FIGURE
# ============================================================

cat("\n============================================================\n")
cat("STEP 48 — IMPROVED NETWORK FIGURE\n")
cat("============================================================\n\n")

library(igraph)
library(ggraph)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load network
# ------------------------------------------------------------

network_df <- readRDS(
  "Results/Network/Core_Gene_Pathway_Network.rds"
)

# ------------------------------------------------------------
# 2. Create edges
# ------------------------------------------------------------

edges <- network_df[, c("Gene", "Pathway")]
colnames(edges) <- c("from", "to")

# ------------------------------------------------------------
# 3. Create nodes
# ------------------------------------------------------------

gene_nodes <- unique(
  network_df[, c(
    "Gene",
    "log2FoldChange",
    "padj",
    "Regulation"
  )]
)

colnames(gene_nodes)[1] <- "name"
gene_nodes$NodeType <- "Gene"

pathway_nodes <- data.frame(
  name = unique(network_df$Pathway),
  log2FoldChange = NA,
  padj = NA,
  Regulation = NA,
  NodeType = "Pathway",
  stringsAsFactors = FALSE
)

nodes <- rbind(
  gene_nodes,
  pathway_nodes
)

# ------------------------------------------------------------
# 4. Create graph
# ------------------------------------------------------------

graph <- graph_from_data_frame(
  edges,
  vertices = nodes,
  directed = FALSE
)

# TRUE = pathway
V(graph)$type <- V(graph)$NodeType == "Pathway"

# ------------------------------------------------------------
# 5. Node size
# ------------------------------------------------------------

V(graph)$size_value <- ifelse(
  V(graph)$NodeType == "Gene",
  5 + abs(V(graph)$log2FoldChange) * 3,
  10
)

# ------------------------------------------------------------
# 6. Improved plot
# ------------------------------------------------------------

network_plot_improved <- ggraph(
  graph,
  layout = "bipartite"
) +
  
  geom_edge_link(
    linewidth = 0.5,
    alpha = 0.35
  ) +
  
  geom_node_point(
    aes(
      size = size_value,
      shape = NodeType
    )
  ) +
  
  geom_node_text(
    aes(
      label = name
    ),
    repel = TRUE,
    size = 3.8,
    fontface = ifelse(
      V(graph)$NodeType == "Pathway",
      "bold",
      "plain"
    )
  ) +
  
  scale_size_identity() +
  
  scale_shape_manual(
    values = c(
      Gene = 16,
      Pathway = 17
    )
  ) +
  
  labs(
    title = "Core Gene–Pathway Network",
    subtitle = paste(
      "12 significant core genes across 3 Reactome pathways"
    ),
    shape = "Node type"
  ) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5
    ),
    legend.position = "bottom",
    plot.margin = margin(
      20, 40, 20, 40
    )
  )

# ------------------------------------------------------------
# 7. Display
# ------------------------------------------------------------

print(network_plot_improved)

# ------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------

ggsave(
  "Results/Network/Figures/Core_Gene_Pathway_Network_Improved.png",
  network_plot_improved,
  width = 14,
  height = 10,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n\n")

cat("Nodes:", vcount(graph), "\n")
cat("Edges:", ecount(graph), "\n")

cat(
  "Improved plot:",
  file.exists(
    "Results/Network/Figures/Core_Gene_Pathway_Network_Improved.png"
  ),
  "\n"
)

cat("\nSTEP 48 completed successfully.\n")

# ============================================================
# STEP 49 — FINAL CYTOSCAPE FILE VALIDATION
# ============================================================

cat("\n============================================================\n")
cat("STEP 49 — CYTOSCAPE FILE VALIDATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Load files
# ------------------------------------------------------------

edges_cyto <- read.csv(
  "Results/Network/Cytoscape/Core_Network_Edges.csv",
  stringsAsFactors = FALSE
)

nodes_cyto <- read.csv(
  "Results/Network/Cytoscape/Core_Network_Nodes.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Basic dimensions
# ------------------------------------------------------------

cat("Edges:", nrow(edges_cyto), "\n")
cat("Nodes:", nrow(nodes_cyto), "\n")

# ------------------------------------------------------------
# Node types
# ------------------------------------------------------------

cat("\nNode types:\n")
print(table(nodes_cyto$NodeType))

# ------------------------------------------------------------
# Check genes
# ------------------------------------------------------------

gene_nodes <- nodes_cyto[
  nodes_cyto$NodeType == "Gene",
]

cat("\nGene nodes:", nrow(gene_nodes), "\n")

cat("\nCore genes:\n")
print(gene_nodes$Node)

# ------------------------------------------------------------
# Check pathways
# ------------------------------------------------------------

pathway_nodes <- nodes_cyto[
  nodes_cyto$NodeType == "Pathway",
]

cat("\nPathway nodes:", nrow(pathway_nodes), "\n")

cat("\nPathways:\n")
print(pathway_nodes$Node)

# ------------------------------------------------------------
# Check edge integrity
# ------------------------------------------------------------

cat("\nAll edge sources are valid:",
    all(edges_cyto$source %in% nodes_cyto$Node),
    "\n")

cat("All edge targets are valid:",
    all(edges_cyto$target %in% nodes_cyto$Node),
    "\n")

# ------------------------------------------------------------
# Check duplicate edges
# ------------------------------------------------------------

cat(
  "Duplicate edges:",
  sum(duplicated(edges_cyto[, c("source", "target")])),
  "\n"
)

# ------------------------------------------------------------
# Final validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL VALIDATION\n")
cat("============================================================\n\n")

cat(
  "Edges file:",
  file.exists(
    "Results/Network/Cytoscape/Core_Network_Edges.csv"
  ),
  "\n"
)

cat(
  "Nodes file:",
  file.exists(
    "Results/Network/Cytoscape/Core_Network_Nodes.csv"
  ),
  "\n"
)

cat("\nSTEP 49 completed successfully.\n")
# STEP 52 — Prepare Core Genes for STRING

string_genes <- core_genes$SYMBOL

cat("Genes prepared for STRING:", length(string_genes), "\n")
print(string_genes)

# STEP 53 — STRING/PPI analysis

if (!requireNamespace("STRINGdb", quietly = TRUE)) {
  BiocManager::install("STRINGdb", ask = FALSE, update = FALSE)
}

library(STRINGdb)

string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400,
  input_directory = ""
)

string_mapped <- string_db$map(
  data.frame(gene = string_genes),
  "gene",
  removeUnmappedRows = TRUE
)

cat("Input genes:", length(string_genes), "\n")
cat("Mapped genes:", nrow(string_mapped), "\n")

print(string_mapped[, c("gene", "STRING_id")])

# STEP 54 — Extract STRING interactions

string_ids <- string_mapped$STRING_id

ppi <- string_db$get_interactions(string_ids)

cat("Number of PPI interactions:", nrow(ppi), "\n")
print(head(ppi))


# STEP 55 — Save STRING PPI network

dir.create("Results/Network/STRING", recursive = TRUE, showWarnings = FALSE)

# Add DESeq2 information to STRING nodes
string_nodes <- string_mapped

string_nodes <- merge(
  string_nodes,
  core_genes[, c("SYMBOL", "log2FoldChange", "padj", "Regulation")],
  by.x = "gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

# Save interactions
write.csv(
  ppi,
  "Results/Network/STRING/Core_Genes_STRING_PPI.csv",
  row.names = FALSE
)

# Save nodes
write.csv(
  string_nodes,
  "Results/Network/STRING/Core_Genes_STRING_Nodes.csv",
  row.names = FALSE
)

cat("PPI edges saved:", file.exists(
  "Results/Network/STRING/Core_Genes_STRING_PPI.csv"
), "\n")

cat("PPI nodes saved:", file.exists(
  "Results/Network/STRING/Core_Genes_STRING_Nodes.csv"
), "\n")

cat("Edges:", nrow(ppi), "\n")
cat("Nodes:", nrow(string_nodes), "\n")

# STEP 56 — Calculate PPI degree

library(igraph)

ppi_graph <- graph_from_data_frame(
  ppi[, c("from_name", "to_name", "combined_score")],
  directed = FALSE,
  vertices = string_nodes[, c("gene", "log2FoldChange", "padj", "Regulation")]
)

V(ppi_graph)$degree <- degree(ppi_graph)

hub_genes <- data.frame(
  Gene = V(ppi_graph)$name,
  Degree = V(ppi_graph)$degree
)

hub_genes <- hub_genes[order(-hub_genes$Degree), ]

print(hub_genes)
# STEP 57 — Save hub gene analysis

write.csv(
  hub_genes,
  "Results/Network/STRING/PPI_Hub_Genes.csv",
  row.names = FALSE
)

cat(
  "Hub genes file saved:",
  file.exists("Results/Network/STRING/PPI_Hub_Genes.csv"),
  "\n"
)

cat("Top hub gene:", hub_genes$Gene[1], "\n")
cat("Top degree:", hub_genes$Degree[1], "\n")

# STEP 58 — Prepare STRING edges for Cytoscape

cyto_string_edges <- ppi[, c(
  "from_name",
  "to_name",
  "combined_score"
)]

colnames(cyto_string_edges) <- c(
  "source",
  "target",
  "combined_score"
)

write.csv(
  cyto_string_edges,
  "Results/Network/STRING/Core_Genes_STRING_Cytoscape_Edges.csv",
  row.names = FALSE
)

cat(
  "Cytoscape STRING edges saved:",
  file.exists(
    "Results/Network/STRING/Core_Genes_STRING_Cytoscape_Edges.csv"
  ),
  "\n"
)

cat("Edges:", nrow(cyto_string_edges), "\n")
print(head(cyto_string_edges))

# STEP 59 — Prepare STRING nodes for Cytoscape

cyto_string_nodes <- merge(
  hub_genes,
  core_genes[, c("SYMBOL", "log2FoldChange", "padj", "Regulation")],
  by.x = "Gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

cyto_string_nodes$NodeType <- "Gene"

colnames(cyto_string_nodes)[1] <- "shared.name"

cyto_string_nodes <- cyto_string_nodes[, c(
  "shared.name",
  "Degree",
  "log2FoldChange",
  "padj",
  "Regulation",
  "NodeType"
)]

write.csv(
  cyto_string_nodes,
  "Results/Network/STRING/Core_Genes_STRING_Cytoscape_Nodes.csv",
  row.names = FALSE
)

cat(
  "Cytoscape STRING nodes saved:",
  file.exists(
    "Results/Network/STRING/Core_Genes_STRING_Cytoscape_Nodes.csv"
  ),
  "\n"
)

cat("Nodes:", nrow(cyto_string_nodes), "\n")
print(cyto_string_nodes)

# STEP 60 — Final STRING network validation

edge_genes <- unique(c(
  cyto_string_edges$source,
  cyto_string_edges$target
))

node_genes <- cyto_string_nodes$shared.name

cat("Total nodes:", length(node_genes), "\n")
cat("Total edges:", nrow(cyto_string_edges), "\n")
cat("Unique genes in edges:", length(edge_genes), "\n")

cat(
  "All edge genes exist in nodes:",
  all(edge_genes %in% node_genes),
  "\n"
)

cat(
  "Duplicate edges:",
  sum(duplicated(
    apply(
      cyto_string_edges[, c("source", "target")],
      1,
      function(x) paste(sort(x), collapse = "_")
    )
  )),
  "\n"
)
# STEP 61 — Extended PPI network centrality

V(ppi_graph)$degree <- degree(ppi_graph)
V(ppi_graph)$betweenness <- betweenness(
  ppi_graph,
  directed = FALSE,
  normalized = TRUE
)
V(ppi_graph)$closeness <- closeness(
  ppi_graph,
  normalized = TRUE
)

ppi_centrality <- data.frame(
  Gene = V(ppi_graph)$name,
  Degree = V(ppi_graph)$degree,
  Betweenness = V(ppi_graph)$betweenness,
  Closeness = V(ppi_graph)$closeness
)

ppi_centrality <- ppi_centrality[
  order(-ppi_centrality$Degree,
        -ppi_centrality$Betweenness),
]

print(ppi_centrality)

# STEP 62 — Save PPI centrality analysis

ppi_centrality <- merge(
  ppi_centrality,
  core_genes[, c("SYMBOL", "log2FoldChange", "padj", "Regulation")],
  by.x = "Gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

write.csv(
  ppi_centrality,
  "Results/Network/STRING/PPI_Centrality_Analysis.csv",
  row.names = FALSE
)

cat(
  "Centrality file saved:",
  file.exists(
    "Results/Network/STRING/PPI_Centrality_Analysis.csv"
  ),
  "\n"
)

cat("Genes analyzed:", nrow(ppi_centrality), "\n")

# STEP 63 — STRING PPI Network Figure

set.seed(123)

# Add edge weights
E(ppi_graph)$weight <- E(ppi_graph)$combined_score / 1000

# Node size based on degree
node_size <- 8 + (V(ppi_graph)$degree * 2)

# Node labels
node_labels <- V(ppi_graph)$name

png(
  "Results/Network/STRING/Core_Genes_STRING_PPI_Network.png",
  width = 2400,
  height = 2000,
  res = 300
)

plot(
  ppi_graph,
  layout = layout_with_fr(ppi_graph),
  vertex.size = node_size,
  vertex.label = node_labels,
  vertex.label.cex = 0.9,
  vertex.label.color = "black",
  vertex.frame.color = "black",
  edge.width = 1 + (E(ppi_graph)$weight * 3),
  edge.color = "grey70",
  main = "STRING PPI Network of Core Genes"
)

dev.off()

cat(
  "PPI network figure saved:",
  file.exists(
    "Results/Network/STRING/Core_Genes_STRING_PPI_Network.png"
  ),
  "\n"
)
# STEP 64 — Link PPI hubs to GSEA core pathways

leading_edge <- read.csv(
  "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.csv",
  stringsAsFactors = FALSE
)

# Keep only the 12 significant core genes
hub_core <- ppi_centrality[
  ppi_centrality$Gene %in% core_genes$SYMBOL,
]

# Pathway membership
pathway_cols <- c(
  "G2M_Checkpoints",
  "Cell_Cycle_Checkpoints",
  "DNA_Replication"
)

hub_core_pathways <- leading_edge[
  leading_edge$SYMBOL %in% hub_core$Gene,
  c("SYMBOL", pathway_cols)
]

# Merge centrality with pathway membership
hub_core_pathways <- merge(
  hub_core,
  hub_core_pathways,
  by.x = "Gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

# Save integrated table
write.csv(
  hub_core_pathways,
  "Results/Network/STRING/Hub_Genes_GSEA_Integration.csv",
  row.names = FALSE
)

cat(
  "Integrated file saved:",
  file.exists(
    "Results/Network/STRING/Hub_Genes_GSEA_Integration.csv"
  ),
  "\n"
)

cat("Genes integrated:", nrow(hub_core_pathways), "\n")

print(hub_core_pathways)

# STEP 65 — Final integrated core signature

final_core_signature <- hub_core_pathways

final_core_signature$N_Core_Pathways <- rowSums(
  final_core_signature[, c(
    "G2M_Checkpoints",
    "Cell_Cycle_Checkpoints",
    "DNA_Replication"
  )]
)

final_core_signature <- final_core_signature[
  order(
    -final_core_signature$Degree,
    -final_core_signature$log2FoldChange
  ),
]

write.csv(
  final_core_signature,
  "Results/Network/STRING/Final_Integrated_Core_Signature.csv",
  row.names = FALSE
)

cat(
  "Final signature saved:",
  file.exists(
    "Results/Network/STRING/Final_Integrated_Core_Signature.csv"
  ),
  "\n"
)

cat(
  "Genes:",
  nrow(final_core_signature),
  "\n"
)

cat(
  "Genes shared across all 3 pathways:",
  sum(final_core_signature$N_Core_Pathways == 3),
  "\n"
)

print(
  final_core_signature[, c(
    "Gene",
    "Degree",
    "Betweenness",
    "log2FoldChange",
    "padj",
    "N_Core_Pathways"
  )]
)
# STEP 66 — Final Integrated Core Signature Figure

library(ggplot2)

p_final_signature <- ggplot(
  final_core_signature,
  aes(
    x = log2FoldChange,
    y = Degree,
    size = Betweenness,
    label = Gene
  )
) +
  geom_point() +
  geom_text(
    vjust = -0.8,
    size = 3.5
  ) +
  scale_size_continuous(
    range = c(4, 10)
  ) +
  labs(
    title = "Integrated Core Gene Signature",
    subtitle = "Expression change and STRING network centrality",
    x = "log2 Fold Change",
    y = "STRING Degree",
    size = "Betweenness"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

ggsave(
  "Results/Network/STRING/Final_Integrated_Core_Signature.png",
  plot = p_final_signature,
  width = 9,
  height = 7,
  dpi = 300
)

cat(
  "Final integrated figure saved:",
  file.exists(
    "Results/Network/STRING/Final_Integrated_Core_Signature.png"
  ),
  "\n"
)
# STEP 67 — Final Master Core Gene Table

master_core_table <- final_core_signature[, c(
  "Gene",
  "log2FoldChange",
  "padj",
  "Degree",
  "Betweenness",
  "Closeness",
  "N_Core_Pathways"
)]

master_core_table <- master_core_table[
  order(-master_core_table$Degree,
        -master_core_table$log2FoldChange),
]

write.csv(
  master_core_table,
  "Results/Network/STRING/MASTER_Core_Gene_Table.csv",
  row.names = FALSE
)

cat(
  "Master table saved:",
  file.exists(
    "Results/Network/STRING/MASTER_Core_Gene_Table.csv"
  ),
  "\n"
)

cat("Genes:", nrow(master_core_table), "\n")

print(master_core_table)

# STEP 68 — Inspect current project metadata

cat("Metadata dimensions:", dim(metadata_all), "\n\n")

cat("Metadata columns:\n")
print(names(metadata_all))

cat("\nFirst 6 rows:\n")
print(head(metadata_all))
# STEP 69 — Extract metadata from current DESeq2 object

metadata_primary <- as.data.frame(colData(dds_primary))

cat("Metadata dimensions:", dim(metadata_primary), "\n\n")

cat("Metadata columns:\n")
print(names(metadata_primary))

cat("\nFirst 6 rows:\n")
print(head(metadata_primary))
# STEP 70 — Inspect clinical metadata values

cat("Patient category:\n")
print(table(metadata_primary$patient_category, useNA = "ifany"))

cat("\nCOVID status:\n")
print(table(metadata_primary$covid_status, useNA = "ifany"))

cat("\nTime point:\n")
print(table(metadata_primary$time_point, useNA = "ifany"))

cat("\nAcuity:\n")
print(table(metadata_primary$acuity, useNA = "ifany"))

cat("\nCell type:\n")
print(table(metadata_primary$cell_type, useNA = "ifany"))
# STEP 71 — Extract and validate patient IDs

metadata_primary$patient_id <- sub("_.*$", "", rownames(metadata_primary))

cat("Number of unique patients:", 
    length(unique(metadata_primary$patient_id)), "\n\n")

cat("Samples per patient — first 20 patients:\n")
print(head(table(metadata_primary$patient_id), 20))

cat("\nExample longitudinal patients:\n")
print(
  head(
    metadata_primary[
      metadata_primary$patient_id %in%
        names(sort(table(metadata_primary$patient_id), decreasing = TRUE))[1:10],
      c("patient_id", "patient_category", "time_point", "acuity")
    ],
    30
  )
)
# STEP 72 — Longitudinal coverage per patient

patient_time_table <- table(
  metadata_primary$patient_id,
  metadata_primary$time_point
)

cat("Patient × Time-point matrix dimensions:\n")
print(dim(patient_time_table))

cat("\nNumber of patients with each time point:\n")
print(colSums(patient_time_table > 0))

cat("\nNumber of patients with each longitudinal pattern:\n")

patient_patterns <- apply(
  patient_time_table,
  1,
  function(x) paste(names(x)[x > 0], collapse = " + ")
)

print(sort(table(patient_patterns), decreasing = TRUE))

# STEP 73 — Check whether acuity is constant within each patient

acuity_by_patient <- tapply(
  metadata_primary$acuity,
  metadata_primary$patient_id,
  function(x) length(unique(x))
)

cat("Patients with constant acuity:",
    sum(acuity_by_patient == 1), "\n")

cat("Patients with changing acuity:",
    sum(acuity_by_patient > 1), "\n")

cat("\nPatients with changing acuity:\n")

changing_acuity_patients <- names(acuity_by_patient[acuity_by_patient > 1])

print(
  metadata_primary[
    metadata_primary$patient_id %in% changing_acuity_patients,
    c("patient_id", "patient_category", "time_point", "acuity")
  ]
)
# STEP 74 — Define complete longitudinal cohort

complete_longitudinal_patients <- names(
  which(
    rowSums(patient_time_table[, c("D0", "D3", "D7")] > 0) == 3
  )
)

metadata_longitudinal <- metadata_primary[
  metadata_primary$patient_id %in% complete_longitudinal_patients,
]

cat("Complete longitudinal patients:",
    length(complete_longitudinal_patients), "\n")

cat("Samples in longitudinal cohort:",
    nrow(metadata_longitudinal), "\n")

cat("\nTime-point distribution:\n")
print(table(metadata_longitudinal$time_point))

cat("\nPatient category distribution:\n")
print(table(
  metadata_longitudinal$patient_id,
  metadata_longitudinal$patient_category
))
# STEP 75 — Check duplicate samples within patient × time point

patient_time_counts <- as.data.frame(
  table(
    metadata_longitudinal$patient_id,
    metadata_longitudinal$time_point
  )
)

names(patient_time_counts) <- c(
  "patient_id",
  "time_point",
  "n_samples"
)

duplicates <- patient_time_counts[
  patient_time_counts$n_samples > 1,
]

cat("Patient × time-point combinations with >1 sample:",
    nrow(duplicates), "\n\n")

cat("Duplicate combinations:\n")
print(duplicates)

cat("\nCorresponding metadata:\n")

if (nrow(duplicates) > 0) {
  
  duplicate_ids <- paste(
    duplicates$patient_id,
    duplicates$time_point,
    sep = "_"
  )
  
  print(
    metadata_longitudinal[
      paste(
        metadata_longitudinal$patient_id,
        metadata_longitudinal$time_point,
        sep = "_"
      ) %in% duplicate_ids,
      c(
        "patient_id",
        "sample_id",
        "geo_accession",
        "time_point",
        "acuity",
        "rna_sample_id"
      )
    ]
  )
}

# STEP 76 — Compare duplicate D7 samples

duplicate_d7 <- c("154_D7A", "154_D7B", "43_D7A", "43_D7B")

duplicate_d7_counts <- counts(dds_primary, normalized = TRUE)[
  ,
  intersect(duplicate_d7, colnames(dds_primary))
]

cat("Duplicate D7 samples found:\n")
print(colnames(duplicate_d7_counts))

cat("\nSample correlations:\n")
print(round(cor(duplicate_d7_counts, method = "spearman"), 3))

# STEP 77 — Inspect metadata for duplicate D7 samples

duplicate_samples <- c(
  "154_D7A",
  "154_D7B",
  "43_D7A",
  "43_D7B"
)

print(
  metadata_primary[
    rownames(metadata_primary) %in% duplicate_samples,
    ,
    drop = FALSE
  ]
)
# ============================================================
# STEP 73 — Validate completed analysis before moving forward
# ============================================================

required_files <- c(
  "Results/Differential_Expression/COVID_vs_Symptomatic_all_genes.rds",
  "Results/Differential_Expression/COVID_vs_Symptomatic_significant_DEGs.csv",
  "Results/Enrichment/GSEA/GO_BP_GSEA.rds",
  "Results/Enrichment/GSEA/Reactome_GSEA.rds",
  "Results/Enrichment/GSEA/Core_Significant_Genes.csv",
  "Results/Network/STRING/Core_Genes_STRING_PPI.csv",
  "Results/Network/STRING/PPI_Centrality_Analysis.csv",
  "Results/Network/STRING/Final_Integrated_Core_Signature.csv",
  "Results/Network/STRING/MASTER_Core_Gene_Table.csv"
)

validation <- data.frame(
  File = required_files,
  Exists = file.exists(required_files)
)

print(validation)

cat("\nAll required files available:",
    all(validation$Exists), "\n")

# ============================================================
# STEP 74 — Inspect final integrated core signature
# ============================================================

core_signature <- read.csv(
  "Results/Network/STRING/Final_Integrated_Core_Signature.csv",
  stringsAsFactors = FALSE
)

cat("Number of core genes:", nrow(core_signature), "\n\n")

print(core_signature)

cat("\nRegulation:\n")
print(table(core_signature$Regulation))

cat("\nShared Reactome pathways:\n")
print(table(core_signature$N_Core_Pathways))

cat("\nDegree summary:\n")
print(summary(core_signature$Degree))

# ============================================================
# STEP 75 — Rank Core Genes
# ============================================================

core_ranked <- core_signature

# Rank by network degree first, then statistical significance
core_ranked <- core_ranked[
  order(
    -core_ranked$Degree,
    core_ranked$padj,
    -core_ranked$log2FoldChange
  ),
]

# Add rank
core_ranked$Rank <- seq_len(nrow(core_ranked))

# Reorder columns
core_ranked <- core_ranked[, c(
  "Rank",
  "Gene",
  "Degree",
  "Betweenness",
  "Closeness",
  "log2FoldChange",
  "padj",
  "Regulation",
  "N_Core_Pathways"
)]

print(core_ranked)

# Save ranked table
write.csv(
  core_ranked,
  "Results/Network/STRING/Ranked_Core_Genes.csv",
  row.names = FALSE
)

saveRDS(
  core_ranked,
  "Results/Network/STRING/Ranked_Core_Genes.rds"
)

cat("\nRanked core genes saved successfully.\n")

# ============================================================
# STEP 76 — Core Gene log2FC Visualization
# ============================================================

library(ggplot2)

core_plot <- ggplot(
  core_ranked,
  aes(
    x = reorder(Gene, log2FoldChange),
    y = log2FoldChange
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Core Gene Signature",
    subtitle = "COVID+ vs Symptomatic Controls",
    x = "Core Gene",
    y = "log2 Fold Change"
  ) +
  theme_minimal(base_size = 12)

print(core_plot)

ggsave(
  "Results/Network/STRING/Core_Gene_log2FC_Barplot.png",
  plot = core_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat(
  "\nPlot saved:",
  file.exists(
    "Results/Network/STRING/Core_Gene_log2FC_Barplot.png"
  ),
  "\n"
)

# ============================================================
# STEP 77 — Final Core Signature Annotation
# ============================================================

core_annotated <- core_ranked

core_annotated$Functional_Group <- NA_character_

core_annotated$Functional_Group[
  core_annotated$Gene %in% c("CCNA1", "CCNA2")
] <- "Cell cycle regulation"

core_annotated$Functional_Group[
  core_annotated$Gene %in% c("CDC6", "CDC45", "MCM2", "MCM4", "MCM10", "ORC1")
] <- "DNA replication"

core_annotated$Functional_Group[
  core_annotated$Gene %in% c("H2BC5", "H2BC7", "H2BC9", "H2BC17")
] <- "Chromatin / Histone"

print(core_annotated[, c(
  "Rank",
  "Gene",
  "Functional_Group",
  "Degree",
  "log2FoldChange",
  "padj",
  "N_Core_Pathways"
)])

# Save final annotated table
write.csv(
  core_annotated,
  "Results/Network/STRING/Final_Core_Signature_Annotated.csv",
  row.names = FALSE
)

saveRDS(
  core_annotated,
  "Results/Network/STRING/Final_Core_Signature_Annotated.rds"
)

cat(
  "\nAnnotated core signature saved:",
  file.exists(
    "Results/Network/STRING/Final_Core_Signature_Annotated.csv"
  ),
  "\n"
)
# Olink file
olink_file <- "/mnt/data/Olink Proteomics.xlsx"

# Read Olink data
olink_data <- readxl::read_excel(olink_file)

# Create subject ID and time point
olink_data$subject_id <- sub("_.*$", "", olink_data$`Public ID`)
olink_data$time_point <- sub("^.*_", "", olink_data$`Public ID`)

# Clinical metadata
clinical_file <- "/mnt/data/Clinical Metadata.xlsx"

clinical_metadata <- readxl::read_excel(
  clinical_file,
  sheet = "Subject-level metadata"
)

# Check matching
cat("Olink samples:", nrow(olink_data), "\n")
cat("Olink subjects:", length(unique(olink_data$subject_id)), "\n")

cat("Clinical subjects:", nrow(clinical_metadata), "\n")

cat(
  "Subject overlap:",
  length(
    intersect(
      unique(olink_data$subject_id),
      as.character(clinical_metadata$`Public ID`)
    )
  ),
  "\n"
)

cat("\nOlink time points:\n")
print(table(olink_data$time_point))
library(dplyr)
library(stringr)

olink_check <- olink_data %>%
  mutate(
    Patient_ID = as.integer(str_extract(`Public ID`, "^[0-9]+")),
    Timepoint = str_extract(`Public ID`, "(?<=_)D[037E]")
  )

# Check the extracted IDs and timepoints
olink_check %>%
  select(`Public ID`, Patient_ID, Timepoint) %>%
  head(15)

cat("Unique Olink patients:", n_distinct(olink_check$Patient_ID), "\n")
cat("Unique Olink timepoints:", paste(sort(unique(olink_check$Timepoint)), collapse = ", "), "\n")

cat("Patient IDs matching Clinical Metadata:",
    sum(unique(olink_check$Patient_ID) %in% clinical_ids), "\n")
clinical_long_check <- clinical_subject %>%
  mutate(
    Patient_ID = as.integer(`Public ID`)
  ) %>%
  select(Patient_ID, D0_draw, D3_draw, D7_draw, DE_draw)

olink_check %>%
  count(Timepoint)

missing_patients <- setdiff(
  clinical_long_check$Patient_ID,
  unique(olink_check$Patient_ID)
)

cat("Patients in Clinical Metadata but absent from Olink:",
    length(missing_patients), "\n")

missing_patients

olink_check %>%
  select(`Public ID`, Patient_ID, Timepoint) %>%
  left_join(
    clinical_subject %>%
      mutate(Patient_ID = as.integer(`Public ID`)) %>%
      select(Patient_ID, D0_draw, D3_draw, D7_draw, DE_draw),
    by = "Patient_ID"
  ) %>%
  head(20)
olink_clinical <- olink_check %>%
  left_join(
    clinical_subject %>%
      mutate(Patient_ID = as.integer(`Public ID`)) %>%
      select(-`Public ID`),
    by = "Patient_ID"
  )

dim(olink_clinical)

sum(is.na(olink_clinical$COVID))

head(
  olink_clinical %>%
    select(`Public ID`, Patient_ID, Timepoint, COVID, `Acuity max`),
  10
)

protein_cols <- names(olink_data)[3:ncol(olink_data)]

protein_matrix <- as.matrix(
  olink_clinical[, protein_cols]
)

cat("Number of proteins:", ncol(protein_matrix), "\n")
cat("Number of samples:", nrow(protein_matrix), "\n")
cat("Total missing values:", sum(is.na(protein_matrix)), "\n")
cat("Missing percentage:",
    round(mean(is.na(protein_matrix)) * 100, 2), "%\n")

protein_summary <- data.frame(
  Protein = protein_cols,
  Mean = colMeans(protein_matrix),
  SD = apply(protein_matrix, 2, sd),
  Min = apply(protein_matrix, 2, min),
  Max = apply(protein_matrix, 2, max)
)

summary(protein_summary$Mean)
summary(protein_summary$SD)

cat("Proteins with SD = 0:",
    sum(protein_summary$SD == 0), "\n")

head(
  protein_summary[order(protein_summary$SD, decreasing = TRUE), ],
  10
)
olink_clinical[non_numeric] <- lapply(
  olink_clinical[non_numeric],
  function(x) as.numeric(x)
)

protein_matrix <- as.matrix(
  olink_clinical[, protein_cols]
)

cat("Number of proteins:", ncol(protein_matrix), "\n")
cat("Number of samples:", nrow(protein_matrix), "\n")
cat("Total missing values:", sum(is.na(protein_matrix)), "\n")

non_numeric_after <- protein_cols[
  !sapply(olink_clinical[, protein_cols], is.numeric)
]

cat("Non-numeric protein columns after conversion:",
    length(non_numeric_after), "\n")

protein_summary <- data.frame(
  Protein = protein_cols,
  Mean = colMeans(protein_matrix, na.rm = TRUE),
  SD = apply(protein_matrix, 2, sd, na.rm = TRUE),
  Min = apply(protein_matrix, 2, min, na.rm = TRUE),
  Max = apply(protein_matrix, 2, max, na.rm = TRUE)
)

summary(protein_summary$Mean)
summary(protein_summary$SD)

cat("Proteins with SD = 0:",
    sum(protein_summary$SD == 0, na.rm = TRUE), "\n")

cat("Proteins with any missing values:",
    sum(colSums(is.na(protein_matrix)) > 0), "\n")

head(
  protein_summary[order(protein_summary$SD, decreasing = TRUE), ],
  10
)
missing_summary <- data.frame(
  Protein = protein_cols,
  Missing_N = colSums(is.na(protein_matrix)),
  Missing_Percent = colMeans(is.na(protein_matrix)) * 100
)

missing_summary %>%
  filter(Missing_N > 0)

sample_missing <- data.frame(
  Public_ID = olink_clinical$`Public ID`,
  Missing_N = rowSums(is.na(protein_matrix)),
  Missing_Percent = rowMeans(is.na(protein_matrix)) * 100
)

summary(sample_missing$Missing_N)

sample_missing %>%
  filter(Missing_N > 0)

sample_qc <- data.frame(
  Public_ID = olink_clinical$`Public ID`,
  Timepoint = olink_clinical$Timepoint,
  Median_NPX = apply(protein_matrix, 1, median, na.rm = TRUE),
  Mean_NPX = apply(protein_matrix, 1, mean, na.rm = TRUE)
)

summary(sample_qc$Median_NPX)
summary(sample_qc$Mean_NPX)

sample_qc %>%
  arrange(Median_NPX) %>%
  head(10)

sample_qc %>%
  arrange(desc(Median_NPX)) %>%
  head(10)

sample_qc %>%
  group_by(Timepoint) %>%
  summarise(
    N = n(),
    Median = median(Median_NPX),
    Mean = mean(Median_NPX),
    SD = sd(Median_NPX),
    Min = min(Median_NPX),
    Max = max(Median_NPX)
  )
clinical_subject %>%
  filter(`Public ID` %in% c(59, 321, 344)) %>%
  select(
    `Public ID`,
    COVID,
    `D0_draw`, `D3_draw`, `D7_draw`, `DE_draw`,
    `Acuity 0`, `Acuity 3`, `Acuity 7`, `Acuity 28`, `Acuity max`
  )
duplicate_pairs <- list(
  "59_D3" = c("59_D3", "59_D3.1"),
  "321_D7" = c("321_D7.1", "321_D7.2"),
  "344_DE" = c("344_DE.1", "344_DE.2")
)

duplicate_correlations <- lapply(names(duplicate_pairs), function(x) {
  ids <- duplicate_pairs[[x]]
  
  v1 <- as.numeric(protein_matrix[match(ids[1], olink_clinical$`Public ID`), ])
  v2 <- as.numeric(protein_matrix[match(ids[2], olink_clinical$`Public ID`), ])
  
  data.frame(
    Pair = x,
    Sample_1 = ids[1],
    Sample_2 = ids[2],
    Spearman_Correlation = cor(v1, v2, method = "spearman", use = "complete.obs")
  )
}) %>%
  bind_rows()

duplicate_correlations

duplicate_summary <- olink_check %>%
  count(Patient_ID, Timepoint, name = "n") %>%
  filter(n > 1)

duplicate_summary

replicate_ids <- c(
  "59_D3", "59_D3.1",
  "321_D7.1", "321_D7.2",
  "344_DE.1", "344_DE.2"
)

replicate_ids <- c(
  "59_D3", "59_D3.1",
  "321_D7.1", "321_D7.2",
  "344_DE.1", "344_DE.2"
)

replicate_group <- c(
  "59_D3", "59_D3",
  "321_D7", "321_D7",
  "344_DE", "344_DE"
)

olink_clinical_clean <- olink_clinical %>%
  mutate(
    Replicate_Group = if_else(
      `Public ID` %in% replicate_ids,
      replicate_group[match(`Public ID`, replicate_ids)],
      `Public ID`
    )
  )

olink_clinical_clean <- olink_clinical_clean %>%
  group_by(Replicate_Group) %>%
  summarise(
    Patient_ID = first(Patient_ID),
    Timepoint = first(Timepoint),
    across(
      everything(),
      ~ if (is.numeric(.x)) mean(.x, na.rm = TRUE) else  first(.x)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    `Public ID` = Replicate_Group
  ) %>%
  select(
    `Public ID`,
    Patient_ID,
    Timepoint,
    everything()
  ) %>%
  select(-Replicate_Group)

dim(olink_clinical_clean)
library(dplyr)

olink_clinical_clean %>%
  count(Patient_ID, Timepoint, name = "n") %>%
  filter(n > 1)

protein_matrix_clean <- as.matrix(
  olink_clinical_clean[, protein_cols]
)

sample_qc_clean <- data.frame(
  Public_ID = olink_clinical_clean$`Public ID`,
  Patient_ID = olink_clinical_clean$Patient_ID,
  Timepoint = olink_clinical_clean$Timepoint,
  Median_NPX = apply(protein_matrix_clean, 1, median, na.rm = TRUE),
  Mean_NPX = apply(protein_matrix_clean, 1, mean, na.rm = TRUE),
  SD_NPX = apply(protein_matrix_clean, 1, sd, na.rm = TRUE)
)

summary(sample_qc_clean[, c("Median_NPX", "Mean_NPX", "SD_NPX")])

library(ggplot2)

ggplot(sample_qc_clean, aes(x = Timepoint, y = Median_NPX)) +
  geom_boxplot(outlier.size = 1.5) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  theme_minimal() +
  labs(
    title = "Olink Proteomics QC: Sample Median NPX",
    x = "Timepoint",
    y = "Median NPX"Z
  )

qc_plot <- ggplot(sample_qc_clean, aes(x = Timepoint, y = Median_NPX)) +
  geom_boxplot(outlier.size = 1.5) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  theme_minimal() +
  labs(
    title = "Olink Proteomics QC: Sample Median NPX",
    x = "Timepoint",
    y = "Median NPX"
  )

qc_plot

ggsave(
  filename = file.path(results_dir, "Olink_QC_Median_NPX_by_Timepoint.png"),
  plot = qc_plot,
  width = 9,
  height = 6,
  dpi = 300
)

protein_matrix_clean <- as.matrix(
  olink_clinical_clean[, protein_cols]
)

protein_matrix_clean <- apply(
  protein_matrix_clean,
  2,
  as.numeric
)

protein_matrix_clean <- t(scale(t(protein_matrix_clean)))

pca_proteomics <- prcomp(
  t(protein_matrix_clean),
  center = FALSE,
  scale. = FALSE
)

pca_df <- data.frame(
  Sample_ID = olink_clinical_clean$`Public ID`,
  Patient_ID = olink_clinical_clean$Patient_ID,
  Timepoint = olink_clinical_clean$Timepoint,
  PC1 = pca_proteomics$x[, 1],
  PC2 = pca_proteomics$x[, 2]
)

pca_variance <- (pca_proteomics$sdev^2) /
  sum(pca_proteomics$sdev^2)

pca_df$PC1_percent <- pca_variance[1] * 100
pca_df$PC2_percent <- pca_variance[2] * 100

head(pca_df)

cat(
  "PC1 variance explained:",
  round(pca_variance[1] * 100, 2), "%\n"
)

cat(
  "PC2 variance explained:",
  round(pca_variance[2] * 100, 2), "%\n"
)

protein_matrix_imputed <- protein_matrix_clean

for (j in seq_len(ncol(protein_matrix_imputed))) {
  missing_idx <- is.na(protein_matrix_imputed[, j])
  
  if (any(missing_idx)) {
    protein_matrix_imputed[missing_idx, j] <- median(
      protein_matrix_imputed[, j],
      na.rm = TRUE
    )
  }
}

cat("Missing values before imputation: 18\n")
cat(
  "Missing values after imputation:",
  sum(is.na(protein_matrix_imputed)),
  "\n"
)

cat(
  "Infinite values:",
  sum(is.infinite(protein_matrix_imputed)),
  "\n"
)

saveRDS(
  protein_matrix_imputed,
  file.path(results_dir, "Olink_Protein_Matrix_Imputed.rds")
)

pca_proteomics <- prcomp(
  protein_matrix_imputed,
  center = FALSE,
  scale. = FALSE
)

pca_variance <- (pca_proteomics$sdev^2) /
  sum(pca_proteomics$sdev^2)

pca_proteomics_df <- data.frame(
  Sample_ID = olink_clinical_clean$`Public ID`,
  Patient_ID = olink_clinical_clean$Patient_ID,
  Timepoint = olink_clinical_clean$Timepoint,
  PC1 = pca_proteomics$x[, 1],
  PC2 = pca_proteomics$x[, 2]
)

write.csv(
  pca_proteomics_df,
  file.path(results_dir, "Olink_PCA_Scores.csv"),
  row.names = FALSE
)

saveRDS(
  pca_proteomics_df,
  file.path(results_dir, "Olink_PCA_Scores.rds")
)

pca_variance_df <- data.frame(
  PC = paste0("PC", seq_along(pca_variance)),
  Variance_Explained = pca_variance * 100
)

write.csv(
  pca_variance_df,
  file.path(results_dir, "Olink_PCA_Variance_Explained.csv"),
  row.names = FALSE
)

saveRDS(
  pca_proteomics,
  file.path(results_dir, "Olink_PCA_Object.rds")
)

cat(
  "PC1 variance explained:",
  round(pca_variance[1] * 100, 2),
  "%\n"
)

cat(
  "PC2 variance explained:",
  round(pca_variance[2] * 100, 2),
  "%\n"
)

cat(
  "Total PC1 + PC2:",
  round(sum(pca_variance[1:2]) * 100, 2),
  "%\n"
)

cat(
  "PCA samples:",
  nrow(pca_proteomics_df),
  "\n"
)

cat(
  "PCA proteins:",
  ncol(pca_proteomics),
  "\n"
)
protein_matrix_scaled <- scale(
  protein_matrix_imputed,
  center = TRUE,
  scale = TRUE
)

pca_proteomics <- prcomp(
  protein_matrix_scaled,
  center = FALSE,
  scale. = FALSE
)

pca_variance <- (pca_proteomics$sdev^2) /
  sum(pca_proteomics$sdev^2)

pca_proteomics_df <- data.frame(
  Sample_ID = olink_clinical_clean$`Public ID`,
  Patient_ID = olink_clinical_clean$Patient_ID,
  Timepoint = olink_clinical_clean$Timepoint,
  PC1 = pca_proteomics$x[, 1],
  PC2 = pca_proteomics$x[, 2]
)

pca_proteomics_df$PC1_percent <- pca_variance[1] * 100
pca_proteomics_df$PC2_percent <- pca_variance[2] * 100

pca_variance_df <- data.frame(
  PC = paste0("PC", seq_along(pca_variance)),
  Variance_Explained = pca_variance * 100
)

saveRDS(
  pca_proteomics,
  file.path(results_dir, "Olink_PCA_Object.rds")
)

saveRDS(
  pca_proteomics_df,
  file.path(results_dir, "Olink_PCA_Scores.rds")
)

write.csv(
  pca_proteomics_df,
  file.path(results_dir, "Olink_PCA_Scores.csv"),
  row.names = FALSE
)

saveRDS(
  pca_variance_df,
  file.path(results_dir, "Olink_PCA_Variance_Explained.rds")
)

write.csv(
  pca_variance_df,
  file.path(results_dir, "Olink_PCA_Variance_Explained.csv"),
  row.names = FALSE
)

cat("Samples:", nrow(protein_matrix_imputed), "\n")
cat("Proteins:", ncol(protein_matrix_imputed), "\n")
cat("PC1 variance explained:", round(pca_variance[1] * 100, 2), "%\n")
cat("PC2 variance explained:", round(pca_variance[2] * 100, 2), "%\n")
cat("PC1 + PC2:", round(sum(pca_variance[1:2]) * 100, 2), "%\n")
sum(pca_variance[1:2]) * 100

pca_plot_timepoint <- ggplot(
  pca_proteomics_df,
  aes(x = PC1, y = PC2, color = Timepoint)
) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    title = "Olink Proteomics PCA by Timepoint",
    x = paste0("PC1 (", round(pca_variance[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(pca_variance[2] * 100, 2), "%)"),
    color = "Timepoint"
  )

pca_plot_timepoint

ggsave(
  filename = file.path(
    results_dir,
    "Olink_PCA_by_Timepoint.png"
  ),
  plot = pca_plot_timepoint,
  width = 9,
  height = 7,
  dpi = 300
)

saveRDS(
  pca_plot_timepoint,
  file.path(
    results_dir,
    "Olink_PCA_by_Timepoint.rds"
  )
)

cat(
  "PCA plot saved:",
  file.exists(
    file.path(results_dir, "Olink_PCA_by_Timepoint.png")
  ),
  "\n"
  
)
pca_proteomics_df$COVID <- olink_clinical_clean$COVID

table(pca_proteomics_df$COVID, useNA = "ifany")

pca_plot_covid <- ggplot(
  pca_proteomics_df,
  aes(x = PC1, y = PC2, color = COVID)
) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    title = "Olink Proteomics PCA by COVID Status",
    x = paste0("PC1 (", round(pca_variance[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(pca_variance[2] * 100, 2), "%)"),
    color = "COVID Status"
  )

pca_plot_covid

ggsave(
  filename = file.path(
    results_dir,
    "Olink_PCA_by_COVID_Status.png"
  ),
  plot = pca_plot_covid,
  width = 9,
  height = 7,
  dpi = 300
)

saveRDS(
  pca_plot_covid,
  file.path(
    results_dir,
    "Olink_PCA_by_COVID_Status.rds"
  )
)

cat(
  "PCA COVID plot saved:",
  file.exists(
    file.path(results_dir, "Olink_PCA_by_COVID_Status.png")
  ),
  "\n"
)
pca_proteomics_df$COVID <- factor(
  pca_proteomics_df$COVID,
  levels = c(0, 1),
  labels = c("COVID-", "COVID+")
)

table(pca_proteomics_df$COVID, useNA = "ifany")

pca_plot_covid <- ggplot(
  pca_proteomics_df,
  aes(x = PC1, y = PC2, color = COVID)
) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    title = "Olink Proteomics PCA by COVID Status",
    x = paste0("PC1 (", round(pca_variance[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(pca_variance[2] * 100, 2), "%)"),
    color = "COVID Status"
  )

pca_plot_covid

ggsave(
  filename = file.path(
    results_dir,
    "Olink_PCA_by_COVID_Status.png"
  ),
  plot = pca_plot_covid,
  width = 9,
  height = 7,
  dpi = 300
)

saveRDS(
  pca_plot_covid,
  file.path(
    results_dir,
    "Olink_PCA_by_COVID_Status.rds"
  )
)

cat(
  "PCA COVID plot saved:",
  file.exists(
    file.path(results_dir, "Olink_PCA_by_COVID_Status.png")
  ),
  "\n"
)
olink_d0 <- olink_clinical_clean %>%
  filter(Timepoint == "D0") %>%
  mutate(
    COVID = factor(
      COVID,
      levels = c(0, 1),
      labels = c("COVID-", "COVID+")
    )
  )

cat("D0 samples:", nrow(olink_d0), "\n")

cat("\nCOVID status at D0:\n")
print(table(olink_d0$COVID, useNA = "ifany"))

cat("\nCOVID status × Acuity max:\n")
print(table(
  olink_d0$COVID,
  olink_d0$`Acuity max`,
  useNA = "ifany"
))

cat("\nUnique patients at D0:", length(unique(olink_d0$Patient_ID)), "\n")

d0_summary <- olink_d0 %>%
  count(COVID, name = "N")

print(d0_summary)

write.csv(
  d0_summary,
  file.path(
    results_dir,
    "Olink_D0_COVID_Group_Summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  d0_summary,
  file.path(
    results_dir,
    "Olink_D0_COVID_Group_Summary.rds"
  )
)

d0_protein_matrix <- as.matrix(
  olink_d0[, protein_cols]
)

storage.mode(d0_protein_matrix) <- "numeric"

d0_protein_matrix <- t(d0_protein_matrix)

design_proteomics <- model.matrix(
  ~ COVID,
  data = olink_d0
)

cat(
  "Protein matrix:",
  nrow(d0_protein_matrix),
  "proteins ×",
  ncol(d0_protein_matrix),
  "samples\n"
)

cat(
  "Design matrix:",
  nrow(design_proteomics),
  "samples ×",
  ncol(design_proteomics),
  "coefficients\n"
)

cat(
  "Sample dimensions match:",
  ncol(d0_protein_matrix) == nrow(design_proteomics),
  "\n"
)

cat("\nDesign columns:\n")
print(colnames(design_proteomics))

fit_proteomics <- lmFit(
  d0_protein_matrix,
  design_proteomics
)

fit_proteomics <- eBayes(
  fit_proteomics,
  trend = TRUE,
  robust = TRUE
)

proteomics_de <- topTable(
  fit_proteomics,
  coef = "COVIDCOVID+",
  number = Inf,
  sort.by = "P"
)

proteomics_de$Protein <- rownames(proteomics_de)

proteomics_de <- proteomics_de[, c(
  "Protein",
  setdiff(names(proteomics_de), "Protein")
)]

proteomics_de$Regulation <- ifelse(
  proteomics_de$logFC > 0,
  "Up",
  "Down"
)

proteomics_de$Significant <- (
  proteomics_de$adj.P.Val < 0.05 &
    abs(proteomics_de$logFC) >= 0.5
)

proteomics_sig <- proteomics_de %>%
  filter(Significant)

proteomics_up <- proteomics_sig %>%
  filter(Regulation == "Up")

proteomics_down <- proteomics_sig %>%
  filter(Regulation == "Down")

cat("\nProteins analyzed:", nrow(proteomics_de), "\n")
cat("Significant proteins:", nrow(proteomics_sig), "\n")
cat("Upregulated:", nrow(proteomics_up), "\n")
cat("Downregulated:", nrow(proteomics_down), "\n")

cat("\nTop 20 proteins:\n")

print(
  proteomics_de %>%
    select(
      Protein,
      logFC,
      AveExpr,
      t,
      P.Value,
      adj.P.Val,
      Regulation
    ) %>%
    head(20)
)

write.csv(
  proteomics_de,
  file.path(
    results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  proteomics_sig,
  file.path(
    results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_significant.csv"
  ),
  row.names = FALSE
)

write.csv(
  proteomics_up,
  file.path(
    results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_upregulated.csv"
  ),
  row.names = FALSE
)

write.csv(
  proteomics_down,
  file.path(
    results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_downregulated.csv"
  ),
  row.names = FALSE
)

saveRDS(
  proteomics_de,
  file.path(
    results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_all.rds"
  )
)

saveRDS(
  proteomics_sig,
  file.path(
    results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_significant.rds"
  )
)

saveRDS(
  fit_proteomics,
  file.path(
    results_dir,
    "Olink_limma_D0_fit.rds"
  )
)

cat(
  "\nResults saved:",
  file.exists(
    file.path(
      results_dir,
      "Differential_Proteomics_D0_COVID_vs_Control_all.csv"
    )
  ),
  "\n"
)

library(ggplot2)

proteomics_de$Plot_Group <- "Not significant"

proteomics_de$Plot_Group[
  proteomics_de$adj.P.Val < 0.05 &
    proteomics_de$logFC >= 0.5
] <- "Upregulated"

proteomics_de$Plot_Group[
  proteomics_de$adj.P.Val < 0.05 &
    proteomics_de$logFC <= -0.5
] <- "Downregulated"

proteomics_de$Plot_Group <- factor(
  proteomics_de$Plot_Group,
  levels = c(
    "Downregulated",
    "Not significant",
    "Upregulated"
  )
)

proteomics_volcano <- ggplot(
  proteomics_de,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = Plot_Group
  )
) +
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  geom_vline(
    xintercept = c(-0.5, 0.5),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  theme_minimal() +
  labs(
    title = "Differential Proteomics at D0",
    subtitle = "COVID+ vs symptomatic COVID-",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Protein group"
  )

proteomics_volcano

ggsave(
  filename = file.path(
    results_dir,
    "Differential_Proteomics_D0_Volcano.png"
  ),
  plot = proteomics_volcano,
  width = 10,
  height = 8,
  dpi = 300
)

saveRDS(
  proteomics_volcano,
  file.path(
    results_dir,
    "Differential_Proteomics_D0_Volcano.rds"
  )
)

saveRDS(
  proteomics_up,
  file.path(
    results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_upregulated.rds"
  )
)

saveRDS(
  proteomics_down,
  file.path(
    results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_downregulated.rds"
  )
)

cat(
  "Volcano saved:",
  file.exists(
    file.path(
      results_dir,
      "Differential_Proteomics_D0_Volcano.png"
    )
  ),
  "\n"
)

cat(
  "Upregulated RDS saved:",
  file.exists(
    file.path(
      results_dir,
      "Differential_Proteomics_D0_COVID_vs_Control_upregulated.rds"
    )
  ),
  "\n"
)

cat(
  "Downregulated RDS saved:",
  file.exists(
    file.path(
      results_dir,
      "Differential_Proteomics_D0_COVID_vs_Control_downregulated.rds"
    )
  ),
  "\n"
)


olink_annotation <- data.frame(
  OlinkID = protein_cols,
  stringsAsFactors = FALSE
)

cat("Number of Olink IDs:", nrow(olink_annotation), "\n")

cat("\nFirst 20 Olink IDs:\n")
print(head(olink_annotation, 20))

cat("\nNumber of unique Olink IDs:\n")
print(length(unique(olink_annotation$OlinkID)))

cat("\nExample significant proteins:\n")
print(
  proteomics_de %>%
    filter(Significant) %>%
    select(
      Protein,
      logFC,
      P.Value,
      adj.P.Val,
      Regulation
    ) %>%
    head(10)
)

library(readxl)
library(dplyr)

supp_file <- "C:/Users/ibrah/OneDrive/Desktop/GSE212041_MultiOmics_Project/Plasma proteomics reveals tissue-specific cell death and mediators of cell-cell interactions in severe COVID-19 patients. Filbin et al/Supplemental-Table-2-Olink-Assays-NPX-values_v2.xlsx"

olink_mapping <- read_excel(
  supp_file,
  sheet = "2A-Olink-Assay",
  skip = 1
)

olink_mapping <- olink_mapping %>%
  select(
    OlinkID,
    Assay,
    UniProt,
    MissingFreq,
    LOD,
    Panel,
    Panel_Version
  ) %>%
  mutate(OlinkID = as.character(OlinkID))

cat("Number of mapping rows:", nrow(olink_mapping), "\n")
cat("Number of unique Olink IDs:", length(unique(olink_mapping$OlinkID)), "\n")

proteomics_de_annotated <- proteomics_de %>%
  mutate(OlinkID = as.character(Protein)) %>%
  left_join(
    olink_mapping,
    by = "OlinkID"
  )

cat("\nTotal DE proteins:", nrow(proteomics_de_annotated), "\n")

cat(
  "DE proteins with successful mapping:",
  sum(!is.na(proteomics_de_annotated$Assay)),
  "\n"
)

cat(
  "DE proteins without mapping:",
  sum(is.na(proteomics_de_annotated$Assay)),
  "\n"
)

cat("\nTop significant annotated proteins:\n")

print(
  proteomics_de_annotated %>%
    filter(Significant) %>%
    select(
      OlinkID,
      Assay,
      UniProt,
      Panel,
      logFC,
      P.Value,
      adj.P.Val,
      Regulation
    ) %>%
    head(20)
)

proteomics_results_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics"

dir.create(
  proteomics_results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  olink_mapping,
  file.path(
    proteomics_results_dir,
    "Olink_Protein_Annotation_1429.csv"
  ),
  row.names = FALSE
)

saveRDS(
  olink_mapping,
  file.path(
    proteomics_results_dir,
    "Olink_Protein_Annotation_1429.rds"
  )
)

write.csv(
  proteomics_de_annotated,
  file.path(
    proteomics_results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_Annotated_all.csv"
  ),
  row.names = FALSE
)

saveRDS(
  proteomics_de_annotated,
  file.path(
    proteomics_results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_Annotated_all.rds"
  )
)

proteomics_sig_annotated <- proteomics_de_annotated %>%
  filter(Significant)

write.csv(
  proteomics_sig_annotated,
  file.path(
    proteomics_results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_Annotated_significant.csv"
  ),
  row.names = FALSE
)

saveRDS(
  proteomics_sig_annotated,
  file.path(
    proteomics_results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_Annotated_significant.rds"
  )
)

cat("Annotated DE proteins:", nrow(proteomics_de_annotated), "\n")
cat("Significant annotated proteins:", nrow(proteomics_sig_annotated), "\n")
cat("Upregulated:", sum(proteomics_sig_annotated$Regulation == "Up"), "\n")
cat("Downregulated:", sum(proteomics_sig_annotated$Regulation == "Down"), "\n")

cat("\nFiles saved:\n")

print(
  file.exists(
    c(
      file.path(proteomics_results_dir, "Olink_Protein_Annotation_1429.csv"),
      file.path(proteomics_results_dir, "Olink_Protein_Annotation_1429.rds"),
      file.path(proteomics_results_dir, "Differential_Proteomics_D0_COVID_vs_Control_Annotated_all.csv"),
      file.path(proteomics_results_dir, "Differential_Proteomics_D0_COVID_vs_Control_Annotated_all.rds"),
      file.path(proteomics_results_dir, "Differential_Proteomics_D0_COVID_vs_Control_Annotated_significant.csv"),
      file.path(proteomics_results_dir, "Differential_Proteomics_D0_COVID_vs_Control_Annotated_significant.rds")
    )
  )
)

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

proteomics_up <- proteomics_sig_annotated %>%
  filter(Regulation == "Up")

proteomics_down <- proteomics_sig_annotated %>%
  filter(Regulation == "Down")

up_mapping <- bitr(
  unique(proteomics_up$UniProt),
  fromType = "UNIPROT",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

down_mapping <- bitr(
  unique(proteomics_down$UniProt),
  fromType = "UNIPROT",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

background_mapping <- bitr(
  unique(proteomics_de_annotated$UniProt),
  fromType = "UNIPROT",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

cat("Significant Up proteins:", nrow(proteomics_up), "\n")
cat("Significant Down proteins:", nrow(proteomics_down), "\n")
cat("Total tested proteins:", nrow(proteomics_de_annotated), "\n")

cat("\nUp mapped to Entrez:", nrow(up_mapping), "\n")
cat("Down mapped to Entrez:", nrow(down_mapping), "\n")
cat("Background mapped to Entrez:", nrow(background_mapping), "\n")

cat("\nUnique Up Entrez:", length(unique(up_mapping$ENTREZID)), "\n")
cat("Unique Down Entrez:", length(unique(down_mapping$ENTREZID)), "\n")
cat("Unique Background Entrez:", length(unique(background_mapping$ENTREZID)), "\n")

cat("\nExample Up mapping:\n")
print(head(up_mapping, 10))

cat("\nExample Down mapping:\n")
print(head(down_mapping, 10))

proteomics_up_entrez <- unique(up_mapping$ENTREZID)

proteomics_background_entrez <- unique(background_mapping$ENTREZID)

ego_proteomics_up <- enrichGO(
  gene = proteomics_up_entrez,
  universe = proteomics_background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

cat("GO-BP terms for upregulated proteins:", nrow(as.data.frame(ego_proteomics_up)), "\n")

if (nrow(as.data.frame(ego_proteomics_up)) > 0) {
  print(
    as.data.frame(ego_proteomics_up) %>%
      select(
        ID,
        Description,
        GeneRatio,
        BgRatio,
        pvalue,
        p.adjust,
        qvalue,
        Count
      ) %>%
      head(20)
  )
} else {
  cat("No significant GO-BP terms found.\n")
}
ego_proteomics_up_df <- as.data.frame(ego_proteomics_up)

cat("GO-BP significant terms:", nrow(ego_proteomics_up_df), "\n\n")

print(
  ego_proteomics_up_df %>%
    dplyr::select(
      ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      Count
    ) %>%
    head(20)
)

proteomics_enrichment_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/Enrichment"

dir.create(
  proteomics_enrichment_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ego_proteomics_up_df <- as.data.frame(ego_proteomics_up)

write.csv(
  ego_proteomics_up_df,
  file.path(
    proteomics_enrichment_dir,
    "GO_BP_Proteomics_Upregulated.csv"
  ),
  row.names = FALSE
)

saveRDS(
  ego_proteomics_up,
  file.path(
    proteomics_enrichment_dir,
    "GO_BP_Proteomics_Upregulated.rds"
  )
)

up_dotplot <- dotplot(
  ego_proteomics_up,
  showCategory = min(20, nrow(ego_proteomics_up_df))
) +
  ggplot2::ggtitle(
    "GO Biological Process – Upregulated Plasma Proteins"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5
    )
  )

ggsave(
  filename = file.path(
    proteomics_enrichment_dir,
    "GO_BP_Proteomics_Upregulated_Dotplot.png"
  ),
  plot = up_dotplot,
  width = 11,
  height = 8,
  dpi = 300
)

saveRDS(
  up_dotplot,
  file.path(
    proteomics_enrichment_dir,
    "GO_BP_Proteomics_Upregulated_Dotplot.rds"
  )
)

cat("GO-BP terms saved:", nrow(ego_proteomics_up_df), "\n")

cat(
  "CSV exists:",
  file.exists(
    file.path(
      proteomics_enrichment_dir,
      "GO_BP_Proteomics_Upregulated.csv"
    )
  ),
  "\n"
)

cat(
  "RDS exists:",
  file.exists(
    file.path(
      proteomics_enrichment_dir,
      "GO_BP_Proteomics_Upregulated.rds"
    )
  ),
  "\n"
)

cat(
  "PNG exists:",
  file.exists(
    file.path(
      proteomics_enrichment_dir,
      "GO_BP_Proteomics_Upregulated_Dotplot.png"
    )
  ),
  "\n"
)

ego_proteomics_down <- enrichGO(
  gene = unique(down_mapping$ENTREZID),
  universe = unique(background_mapping$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_proteomics_down_df <- as.data.frame(ego_proteomics_down)

cat(
  "GO-BP terms for downregulated proteins:",
  nrow(ego_proteomics_down_df),
  "\n\n"
)

if (nrow(ego_proteomics_down_df) > 0) {
  print(
    ego_proteomics_down_df %>%
      dplyr::select(
        ID,
        Description,
        GeneRatio,
        BgRatio,
        pvalue,
        p.adjust,
        qvalue,
        Count
      ) %>%
      head(20)
  )
} else {
  cat("No significant GO-BP terms found.\n")
}

ego_proteomics_down_check <- enrichGO(
  gene = unique(down_mapping$ENTREZID),
  universe = unique(background_mapping$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

ego_proteomics_down_check_df <- as.data.frame(
  ego_proteomics_down_check
)

cat(
  "Total GO-BP terms tested:",
  nrow(ego_proteomics_down_check_df),
  "\n"
)

cat(
  "Terms with nominal p < 0.05:",
  sum(ego_proteomics_down_check_df$pvalue < 0.05),
  "\n"
)

cat(
  "Terms with adjusted p < 0.10:",
  sum(ego_proteomics_down_check_df$p.adjust < 0.10),
  "\n"
)

cat("\nTop 20 GO-BP terms by nominal p-value:\n")

print(
  ego_proteomics_down_check_df %>%
    dplyr::arrange(pvalue) %>%
    dplyr::select(
      ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      Count
    ) %>%
    head(20)
)

write.csv(
  ego_proteomics_down_df,
  file.path(
    proteomics_enrichment_dir,
    "GO_BP_Proteomics_Downregulated.csv"
  ),
  row.names = FALSE
)

saveRDS(
  ego_proteomics_down,
  file.path(
    proteomics_enrichment_dir,
    "GO_BP_Proteomics_Downregulated.rds"
  )
)

write.csv(
  ego_proteomics_down_check_df,
  file.path(
    proteomics_enrichment_dir,
    "GO_BP_Proteomics_Downregulated_Diagnostic_All_Terms.csv"
  ),
  row.names = FALSE
)

saveRDS(
  ego_proteomics_down_check,
  file.path(
    proteomics_enrichment_dir,
    "GO_BP_Proteomics_Downregulated_Diagnostic_All_Terms.rds"
  )
)

cat(
  "Official significant GO-BP terms:",
  nrow(ego_proteomics_down_df),
  "\n"
)

cat(
  "Diagnostic nominal p < 0.05:",
  sum(ego_proteomics_down_check_df$pvalue < 0.05),
  "\n"
)

cat(
  "Diagnostic adjusted p < 0.10:",
  sum(ego_proteomics_down_check_df$p.adjust < 0.10),
  "\n"
)

cat("\nFiles saved:\n")

print(
  file.exists(
    c(
      file.path(
        proteomics_enrichment_dir,
        "GO_BP_Proteomics_Downregulated.csv"
      ),
      file.path(
        proteomics_enrichment_dir,
        "GO_BP_Proteomics_Downregulated.rds"
      ),
      file.path(
        proteomics_enrichment_dir,
        "GO_BP_Proteomics_Downregulated_Diagnostic_All_Terms.csv"
      ),
      file.path(
        proteomics_enrichment_dir,
        "GO_BP_Proteomics_Downregulated_Diagnostic_All_Terms.rds"
      )
    )
  )
)


ego_proteomics_up_mf <- enrichGO(
  gene = unique(up_mapping$ENTREZID),
  universe = unique(background_mapping$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_proteomics_up_mf_df <- as.data.frame(ego_proteomics_up_mf)

cat(
  "GO-MF terms for upregulated proteins:",
  nrow(ego_proteomics_up_mf_df),
  "\n\n"
)

if (nrow(ego_proteomics_up_mf_df) > 0) {
  print(
    ego_proteomics_up_mf_df %>%
      dplyr::select(
        ID,
        Description,
        GeneRatio,
        BgRatio,
        pvalue,
        p.adjust,
        qvalue,
        Count
      ) %>%
      head(20)
  )
} else {
  cat("No significant GO-MF terms found.\n")
}


mf_up_dir <- file.path(results_dir, "Enrichment")
dir.create(mf_up_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  ego_proteomics_up_mf_df,
  file.path(mf_up_dir, "GO_MF_Proteomics_Upregulated.csv"),
  row.names = FALSE
)

saveRDS(
  ego_proteomics_up_mf,
  file.path(mf_up_dir, "GO_MF_Proteomics_Upregulated.rds")
)

ego_proteomics_up_mf_diag <- enrichGO(
  gene = unique(up_mapping$ENTREZID),
  universe = unique(background_mapping$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

ego_proteomics_up_mf_diag_df <- as.data.frame(
  ego_proteomics_up_mf_diag
)

write.csv(
  ego_proteomics_up_mf_diag_df,
  file.path(
    mf_up_dir,
    "GO_MF_Proteomics_Upregulated_Diagnostic_All_Terms.csv"
  ),
  row.names = FALSE
)

saveRDS(
  ego_proteomics_up_mf_diag,
  file.path(
    mf_up_dir,
    "GO_MF_Proteomics_Upregulated_Diagnostic_All_Terms.rds"
  )
)

cat(
  "Significant GO-MF terms:",
  nrow(ego_proteomics_up_mf_df),
  "\n"
)

cat(
  "Total GO-MF terms tested:",
  nrow(ego_proteomics_up_mf_diag_df),
  "\n"
)

cat(
  "Nominal p < 0.05:",
  sum(ego_proteomics_up_mf_diag_df$pvalue < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Adjusted p < 0.05:",
  sum(ego_proteomics_up_mf_diag_df$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Files saved:",
  file.exists(file.path(mf_up_dir, "GO_MF_Proteomics_Upregulated.csv")),
  file.exists(file.path(mf_up_dir, "GO_MF_Proteomics_Upregulated.rds")),
  file.exists(file.path(
    mf_up_dir,
    "GO_MF_Proteomics_Upregulated_Diagnostic_All_Terms.csv"
  )),
  file.exists(file.path(
    mf_up_dir,
    "GO_MF_Proteomics_Upregulated_Diagnostic_All_Terms.rds"
  )),
  "\n"
)

ego_proteomics_up_cc <- enrichGO(
  gene = unique(up_mapping$ENTREZID),
  universe = unique(background_mapping$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_proteomics_up_cc_df <- as.data.frame(ego_proteomics_up_cc)

cat(
  "GO-CC terms for upregulated proteins:",
  nrow(ego_proteomics_up_cc_df),
  "\n\n"
)

if (nrow(ego_proteomics_up_cc_df) > 0) {
  print(
    ego_proteomics_up_cc_df %>%
      dplyr::select(
        ID,
        Description,
        GeneRatio,
        BgRatio,
        pvalue,
        p.adjust,
        qvalue,
        Count
      ) %>%
      head(20)
  )
} else {
  cat("No significant GO-CC terms found.\n")
}

cc_up_dir <- file.path(results_dir, "Enrichment")
dir.create(cc_up_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  ego_proteomics_up_cc_df,
  file.path(cc_up_dir, "GO_CC_Proteomics_Upregulated.csv"),
  row.names = FALSE
)

saveRDS(
  ego_proteomics_up_cc,
  file.path(cc_up_dir, "GO_CC_Proteomics_Upregulated.rds")
)

ego_proteomics_up_cc_diag <- enrichGO(
  gene = unique(up_mapping$ENTREZID),
  universe = unique(background_mapping$ENTREZID),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

ego_proteomics_up_cc_diag_df <- as.data.frame(
  ego_proteomics_up_cc_diag
)

write.csv(
  ego_proteomics_up_cc_diag_df,
  file.path(
    cc_up_dir,
    "GO_CC_Proteomics_Upregulated_Diagnostic_All_Terms.csv"
  ),
  row.names = FALSE
)

saveRDS(
  ego_proteomics_up_cc_diag,
  file.path(
    cc_up_dir,
    "GO_CC_Proteomics_Upregulated_Diagnostic_All_Terms.rds"
  )
)

cat(
  "Significant GO-CC terms:",
  nrow(ego_proteomics_up_cc_df),
  "\n"
)

cat(
  "Total GO-CC terms tested:",
  nrow(ego_proteomics_up_cc_diag_df),
  "\n"
)

cat(
  "Nominal p < 0.05:",
  sum(ego_proteomics_up_cc_diag_df$pvalue < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Adjusted p < 0.05:",
  sum(ego_proteomics_up_cc_diag_df$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Files saved:",
  file.exists(file.path(cc_up_dir, "GO_CC_Proteomics_Upregulated.csv")),
  file.exists(file.path(cc_up_dir, "GO_CC_Proteomics_Upregulated.rds")),
  file.exists(file.path(
    cc_up_dir,
    "GO_CC_Proteomics_Upregulated_Diagnostic_All_Terms.csv"
  )),
  file.exists(file.path(
    cc_up_dir,
    "GO_CC_Proteomics_Upregulated_Diagnostic_All_Terms.rds"
  )),
  "\n"
)
ekegg_proteomics_up <- enrichKEGG(
  gene = unique(up_mapping$ENTREZID),
  universe = unique(background_mapping$ENTREZID),
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ekegg_proteomics_up_df <- as.data.frame(ekegg_proteomics_up)

cat(
  "KEGG terms for upregulated proteins:",
  nrow(ekegg_proteomics_up_df),
  "\n\n"
)

if (nrow(ekegg_proteomics_up_df) > 0) {
  print(
    ekegg_proteomics_up_df %>%
      dplyr::select(
        ID,
        Description,
        GeneRatio,
        BgRatio,
        pvalue,
        p.adjust,
        qvalue,
        Count
      ) %>%
      head(20)
  )
} else {
  cat("No significant KEGG pathways found.\n")
}

kegg_prot_dir <- file.path(results_dir, "Enrichment", "KEGG")
dir.create(kegg_prot_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  ekegg_proteomics_up_df,
  file.path(
    kegg_prot_dir,
    "KEGG_Proteomics_Upregulated.csv"
  ),
  row.names = FALSE
)

saveRDS(
  ekegg_proteomics_up,
  file.path(
    kegg_prot_dir,
    "KEGG_Proteomics_Upregulated.rds"
  )
)

if (nrow(ekegg_proteomics_up_df) > 0) {
  
  kegg_up_plot <- dotplot(
    ekegg_proteomics_up,
    showCategory = min(20, nrow(ekegg_proteomics_up_df))
  ) +
    ggtitle(
      "KEGG Enrichment of Upregulated Proteins"
    ) +
    theme_minimal()
  
  ggsave(
    file.path(
      kegg_prot_dir,
      "KEGG_Proteomics_Upregulated_Dotplot.png"
    ),
    kegg_up_plot,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  saveRDS(
    kegg_up_plot,
    file.path(
      kegg_prot_dir,
      "KEGG_Proteomics_Upregulated_Dotplot.rds"
    )
  )
}

cat(
  "KEGG significant pathways:",
  nrow(ekegg_proteomics_up_df),
  "\n"
)

cat(
  "Files saved:",
  file.exists(file.path(
    kegg_prot_dir,
    "KEGG_Proteomics_Upregulated.csv"
  )),
  file.exists(file.path(
    kegg_prot_dir,
    "KEGG_Proteomics_Upregulated.rds"
  )),
  file.exists(file.path(
    kegg_prot_dir,
    "KEGG_Proteomics_Upregulated_Dotplot.png"
  )),
  file.exists(file.path(
    kegg_prot_dir,
    "KEGG_Proteomics_Upregulated_Dotplot.rds"
  )),
  "\n"
)

ekegg_proteomics_down <- enrichKEGG(
  gene = unique(down_mapping$ENTREZID),
  universe = unique(background_mapping$ENTREZID),
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

ekegg_proteomics_down_df <- as.data.frame(
  ekegg_proteomics_down
)

cat(
  "KEGG terms for downregulated proteins:",
  nrow(ekegg_proteomics_down_df),
  "\n\n"
)

if (nrow(ekegg_proteomics_down_df) > 0) {
  print(
    ekegg_proteomics_down_df %>%
      dplyr::select(
        ID,
        Description,
        GeneRatio,
        BgRatio,
        pvalue,
        p.adjust,
        qvalue,
        Count
      ) %>%
      head(20)
  )
} else {
  cat("No significant KEGG pathways found.\n")
}

kegg_prot_down_dir <- file.path(results_dir, "Enrichment", "KEGG")
dir.create(kegg_prot_down_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  ekegg_proteomics_down_df,
  file.path(
    kegg_prot_down_dir,
    "KEGG_Proteomics_Downregulated.csv"
  ),
  row.names = FALSE
)

saveRDS(
  ekegg_proteomics_down,
  file.path(
    kegg_prot_down_dir,
    "KEGG_Proteomics_Downregulated.rds"
  )
)

cat(
  "KEGG significant downregulated pathways:",
  nrow(ekegg_proteomics_down_df),
  "\n"
)

cat(
  "Files saved:",
  file.exists(file.path(
    kegg_prot_down_dir,
    "KEGG_Proteomics_Downregulated.csv"
  )),
  file.exists(file.path(
    kegg_prot_down_dir,
    "KEGG_Proteomics_Downregulated.rds"
  )),
  "\n"
)
proteomics_ranked <- proteomics_de %>%
  dplyr::select(Protein, t) %>%
  dplyr::filter(
    !is.na(t),
    is.finite(t)
  ) %>%
  dplyr::left_join(
    olink_mapping %>%
      dplyr::select(OlinkID, UniProt, Assay, Panel) %>%
      dplyr::distinct(OlinkID, .keep_all = TRUE),
    by = c("Protein" = "OlinkID")
  ) %>%
  dplyr::filter(
    !is.na(UniProt),
    UniProt != ""
  ) %>%
  dplyr::distinct(UniProt, .keep_all = TRUE) %>%
  dplyr::arrange(desc(t))

proteomics_geneList <- proteomics_ranked$t
names(proteomics_geneList) <- proteomics_ranked$UniProt

cat(
  "Proteins in GSEA ranking:",
  length(proteomics_geneList),
  "\n"
)

cat(
  "Maximum t-statistic:",
  max(proteomics_geneList),
  "\n"
)

cat(
  "Minimum t-statistic:",
  min(proteomics_geneList),
  "\n"
)

cat(
  "NA values:",
  sum(is.na(proteomics_geneList)),
  "\n"
)

cat(
  "Duplicate UniProt IDs:",
  sum(duplicated(names(proteomics_geneList))),
  "\n"
)

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

results_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics"

proteomics_de <- readRDS(
  file.path(
    results_dir,
    "Differential_Proteomics_D0_COVID_vs_Control_all.rds"
  )
)

olink_mapping <- readRDS(
  file.path(
    results_dir,
    "Olink_Protein_Annotation_1429.rds"
  )
)

proteomics_ranked <- proteomics_de %>%
  dplyr::select(Protein, t) %>%
  dplyr::filter(
    !is.na(t),
    is.finite(t)
  ) %>%
  dplyr::left_join(
    olink_mapping %>%
      dplyr::select(OlinkID, UniProt) %>%
      dplyr::distinct(OlinkID, .keep_all = TRUE),
    by = c("Protein" = "OlinkID")
  ) %>%
  dplyr::filter(
    !is.na(UniProt),
    UniProt != ""
  ) %>%
  dplyr::distinct(UniProt, .keep_all = TRUE) %>%
  dplyr::arrange(desc(t))

proteomics_geneList <- proteomics_ranked$t
names(proteomics_geneList) <- proteomics_ranked$UniProt

proteomics_uniprot_entrez <- clusterProfiler::bitr(
  names(proteomics_geneList),
  fromType = "UNIPROT",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

proteomics_uniprot_entrez <- proteomics_uniprot_entrez %>%
  dplyr::distinct(UNIPROT, .keep_all = TRUE)

proteomics_ranked_entrez <- proteomics_ranked %>%
  dplyr::left_join(
    proteomics_uniprot_entrez,
    by = c("UniProt" = "UNIPROT")
  ) %>%
  dplyr::filter(
    !is.na(ENTREZID)
  ) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
  dplyr::arrange(desc(t))

proteomics_geneList_entrez <- proteomics_ranked_entrez$t
names(proteomics_geneList_entrez) <- proteomics_ranked_entrez$ENTREZID

saveRDS(
  proteomics_geneList_entrez,
  file.path(
    results_dir,
    "Proteomics_GSEA_Ranking_Entrez.rds"
  )
)

write.csv(
  proteomics_ranked_entrez,
  file.path(
    results_dir,
    "Proteomics_GSEA_Ranking_Entrez.csv"
  ),
  row.names = FALSE
)

cat(
  "Original UniProt ranking:",
  length(proteomics_geneList),
  "\n"
)

cat(
  "Mapped Entrez ranking:",
  length(proteomics_geneList_entrez),
  "\n"
)

cat(
  "NA values:",
  sum(is.na(proteomics_geneList_entrez)),
  "\n"
)

cat(
  "Duplicate Entrez IDs:",
  sum(duplicated(names(proteomics_geneList_entrez))),
  "\n"
)

cat(
  "Maximum t:",
  max(proteomics_geneList_entrez),
  "\n"
)

cat(
  "Minimum t:",
  min(proteomics_geneList_entrez),
  "\n"
)

cat(
  "Ranking saved:",
  file.exists(
    file.path(results_dir, "Proteomics_GSEA_Ranking_Entrez.rds")
  ),
  "\n"
)
proteomics_gsea_go_bp <- clusterProfiler::gseGO(
  geneList = proteomics_geneList_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

proteomics_gsea_go_bp_df <- as.data.frame(
  proteomics_gsea_go_bp
)

cat(
  "GO-BP GSEA pathways:",
  nrow(proteomics_gsea_go_bp_df),

gsea_go_bp_dir <- file.path(
    results_dir,
    "GSEA"
  )
  
  dir.create(
    gsea_go_bp_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    proteomics_gsea_go_bp,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_GO_BP.rds"
    )
  )
  
  write.csv(
    proteomics_gsea_go_bp_df,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_GO_BP.csv"
    ),
    row.names = FALSE
  )
  
  proteomics_gsea_go_bp_sig <- proteomics_gsea_go_bp_df %>%
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    ) %>%
    dplyr::arrange(p.adjust)
  
  write.csv(
    proteomics_gsea_go_bp_sig,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_GO_BP_significant.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    proteomics_gsea_go_bp_sig[
      1:min(20, nrow(proteomics_gsea_go_bp_sig)),
    ],
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_GO_BP_top20.csv"
    ),
    row.names = FALSE
  )
  
  cat(
    "GSEA object saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_GO_BP.rds"
      )
    ),
    "\n"
  )
  
  cat(
    "Full CSV saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_GO_BP.csv"
      )
    ),
    "\n"
  )
  
  cat(
    "Significant CSV saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_GO_BP_significant.csv"
      )
    ),
    "\n"
  )
  
  cat(
    "Significant pathways saved:",
    nrow(proteomics_gsea_go_bp_sig),
    "\n"
  )    
  library(enrichplot)
  library(ggplot2)
  
  gsea_go_bp_plot <- enrichplot::dotplot(
    proteomics_gsea_go_bp,
    showCategory = 20,
    color = "p.adjust",
    size = "setSize"
  ) +
    ggplot2::ggtitle(
      "Proteomics GSEA — GO Biological Process"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
  
  print(gsea_go_bp_plot)
  
  ggsave(
    filename = file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_GO_BP_Dotplot_Top20.png"
    ),
    plot = gsea_go_bp_plot,
    width = 11,
    height = 9,
    dpi = 300
  )
  
  saveRDS(
    gsea_go_bp_plot,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_GO_BP_Dotplot_Top20.rds"
    )
  )
  
  cat(
    "PNG saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_GO_BP_Dotplot_Top20.png"
      )
    ),
    "\n"
  )
  
  cat(
    "RDS saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_GO_BP_Dotplot_Top20.rds"
      )
    ),
    "\n"
  ) 
  proteomics_gsea_go_bp_top20 <- proteomics_gsea_go_bp_df %>%
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    ) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 20)
  
  proteomics_gsea_go_bp_top20_ids <- proteomics_gsea_go_bp_top20$ID
  
  proteomics_gsea_go_bp_top20_ids  
  
  proteomics_gsea_leading_edge <- proteomics_gsea_go_bp_df %>%
    dplyr::filter(
      ID %in% proteomics_gsea_go_bp_top20_ids
    ) %>%
    dplyr::select(
      ID,
      Description,
      NES,
      p.adjust,
      core_enrichment
    ) %>%
    tidyr::separate_rows(
      core_enrichment,
      sep = "/"
    ) %>%
    dplyr::rename(
      EntrezID = core_enrichment
    ) %>%
    dplyr::left_join(
      proteomics_ranked_entrez %>%
        dplyr::select(
          EntrezID = ENTREZID,
          UniProt,
          Protein,
          t
        ),
      by = "EntrezID"
    ) %>%
    dplyr::arrange(
      p.adjust,
      desc(abs(t))
    )
  
  write.csv(
    proteomics_gsea_leading_edge,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_GO_BP_Top20_LeadingEdge.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    proteomics_gsea_leading_edge,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_GO_BP_Top20_LeadingEdge.rds"
    )
  )
  
  cat(
    "Leading-edge rows:",
    nrow(proteomics_gsea_leading_edge),
    "\n"
  )
  
  cat(
    "Unique proteins:",
    dplyr::n_distinct(
      proteomics_gsea_leading_edge$UniProt,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "Pathways represented:",
    dplyr::n_distinct(
      proteomics_gsea_leading_edge$ID
    ),
    "\n"
  )
  
  cat(
    "CSV saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_GO_BP_Top20_LeadingEdge.csv"
      )
    ),
    "\n"
  )
  
  cat(
    "RDS saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_GO_BP_Top20_LeadingEdge.rds"
      )
    ),
    "\n"
  )
  proteomics_leading_edge_core <- proteomics_gsea_leading_edge %>%
    dplyr::filter(
      !is.na(UniProt),
      UniProt != ""
    ) %>%
    dplyr::group_by(
      UniProt
    ) %>%
    dplyr::summarise(
      Protein = dplyr::first(Protein),
      Pathway_Count = dplyr::n_distinct(ID),
      Mean_NES = mean(NES, na.rm = TRUE),
      Max_abs_t = max(abs(t), na.rm = TRUE),
      Mean_t = mean(t, na.rm = TRUE),
      Pathways = paste(
        unique(Description),
        collapse = " | "
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      desc(Pathway_Count),
      desc(Max_abs_t)
    )
  
  write.csv(
    proteomics_leading_edge_core,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_GO_BP_Core_LeadingEdge_Proteins.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    proteomics_leading_edge_core,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_GO_BP_Core_LeadingEdge_Proteins.rds"
    )
  )
  
  cat(
    "Unique leading-edge proteins:",
    nrow(proteomics_leading_edge_core),
    "\n"
  )
  
  cat(
    "Proteins in >=5 pathways:",
    sum(
      proteomics_leading_edge_core$Pathway_Count >= 5
    ),
    "\n"
  )
  
  cat(
    "Proteins in >=10 pathways:",
    sum(
      proteomics_leading_edge_core$Pathway_Count >= 10
    ),
    "\n"
  )
  
  cat(
    "Maximum pathway recurrence:",
    max(
      proteomics_leading_edge_core$Pathway_Count
    ),
    "\n"
  )
  
  print(
    proteomics_leading_edge_core[
      1:min(20, nrow(proteomics_leading_edge_core)),
      c(
        "UniProt",
        "Protein",
        "Pathway_Count",
        "Mean_NES",
        "Max_abs_t",
        "Mean_t"
      )
    ]
  )
  
  cat(
    "CSV saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_GO_BP_Core_LeadingEdge_Proteins.csv"
      )
    ),
    "\n"
  )
  
  cat(
    "RDS saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_GO_BP_Core_LeadingEdge_Proteins.rds"
      )
    ),
    "\n"
  )
  
  
  proteomics_core_de_overlap <- proteomics_leading_edge_core %>%
    dplyr::filter(
      Pathway_Count >= 5
    ) %>%
    dplyr::left_join(
      proteomics_de %>%
        dplyr::select(
          Protein,
          logFC,
          adj.P.Val,
          Regulation,
          Significant
        ),
      by = "Protein"
    ) %>%
    dplyr::left_join(
      olink_mapping %>%
        dplyr::select(
          OlinkID,
          Assay,
          UniProt
        ) %>%
        dplyr::distinct(OlinkID, .keep_all = TRUE),
      by = c("Protein" = "OlinkID"),
      suffix = c("", "_mapping")
    ) %>%
    dplyr::mutate(
      UniProt = dplyr::coalesce(
        UniProt,
        UniProt_mapping
      )
    ) %>%
    dplyr::select(
      Protein,
      Assay,
      UniProt,
      Pathway_Count,
      Mean_NES,
      Max_abs_t,
      Mean_t,
      logFC,
      adj.P.Val,
      Regulation,
      Significant,
      Pathways
    ) %>%
    dplyr::arrange(
      desc(Significant),
      adj.P.Val,
      desc(Pathway_Count)
    )
  
  write.csv(
    proteomics_core_de_overlap,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_Core_LeadingEdge_DE_Overlap.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    proteomics_core_de_overlap,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_Core_LeadingEdge_DE_Overlap.rds"
    )
  )
  
  cat(
    "Core proteins (>=5 pathways):",
    nrow(proteomics_core_de_overlap),
    "\n"
  )
  
  cat(
    "Significant DE proteins among core:",
    sum(
      proteomics_core_de_overlap$Significant %in% TRUE,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "Upregulated DE proteins among core:",
    sum(
      proteomics_core_de_overlap$Regulation == "Up",
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "Downregulated DE proteins among core:",
    sum(
      proteomics_core_de_overlap$Regulation == "Down",
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "CSV saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_Core_LeadingEdge_DE_Overlap.csv"
      )
    ),
    "\n"
  )
  
  cat(
    "RDS saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_Core_LeadingEdge_DE_Overlap.rds"
      )
    ),
    "\n"
  )
  
  print(
    proteomics_core_de_overlap %>%
      dplyr::filter(
        Significant %in% TRUE
      ) %>%
      dplyr::select(
        Protein,
        Assay,
        UniProt,
        Pathway_Count,
        logFC,
        adj.P.Val,
        Regulation
      ) %>%
      dplyr::arrange(
        adj.P.Val
      ) %>%
      head(30)
  )
  
  
  proteomics_core_de_sig <- proteomics_core_de_overlap %>%
    dplyr::filter(
      Significant %in% TRUE
    ) %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        logFC > 0 ~ "Up",
        logFC < 0 ~ "Down",
        TRUE ~ "Zero"
      )
    ) %>%
    dplyr::arrange(
      Direction,
      adj.P.Val
    )
  
  write.csv(
    proteomics_core_de_sig,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_Core_LeadingEdge_Significant_DE_36.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    proteomics_core_de_sig,
    file.path(
      gsea_go_bp_dir,
      "Proteomics_GSEA_Core_LeadingEdge_Significant_DE_36.rds"
    )
  )
  
  cat(
    "Significant core proteins:",
    nrow(proteomics_core_de_sig),
    "\n"
  )
  
  cat(
    "Significant Up:",
    sum(
      proteomics_core_de_sig$Direction == "Up"
    ),
    "\n"
  )
  
  cat(
    "Significant Down:",
    sum(
      proteomics_core_de_sig$Direction == "Down"
    ),
    "\n"
  )
  
  cat(
    "CSV saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_Core_LeadingEdge_Significant_DE_36.csv"
      )
    ),
    "\n"
  )
  
  cat(
    "RDS saved:",
    file.exists(
      file.path(
        gsea_go_bp_dir,
        "Proteomics_GSEA_Core_LeadingEdge_Significant_DE_36.rds"
      )
    ),
    "\n"
  )
  
  print(
    proteomics_core_de_sig %>%
      dplyr::select(
        Protein,
        Assay,
        UniProt,
        Pathway_Count,
        logFC,
        adj.P.Val,
        Direction
      ) %>%
      dplyr::arrange(
        adj.P.Val
      ),
    n = Inf
  )
  
  
  library(enrichplot)
  library(clusterProfiler)
  library(ggplot2)
  
  gsea_map <- proteomics_gsea_go_bp
  
  top30_ids <- gsea_map@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 30) %>%
    dplyr::pull(ID)
  
  gsea_map@result <- gsea_map@result[
    gsea_map@result$ID %in% top30_ids,
  ]
  
  cat("Pathways selected:", nrow(gsea_map@result), "\n")
  
  # Calculate pairwise similarity first
  gsea_map_sim <- pairwise_termsim(gsea_map)
  
  cat(
    "Similarity matrix dimensions:",
    paste(dim(gsea_map_sim@termsim), collapse = " x "),
    "\n"
  )
  
  # Draw enrichment map
  p_emap <- emapplot(
    gsea_map_sim,
    showCategory = 30
  )
  
  print(p_emap)
  
  ggsave(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/GSEA/Proteomics_GSEA_GO_BP_EnrichmentMap_Top30.png",
    p_emap,
    width = 14,
    height = 11,
    dpi = 300
  )
  
  saveRDS(
    p_emap,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/GSEA/Proteomics_GSEA_GO_BP_EnrichmentMap_Top30.rds"
  )
  
  cat(
    "PNG exists:",
    file.exists(
      "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/GSEA/Proteomics_GSEA_GO_BP_EnrichmentMap_Top30.png"
    ),
    "\n"
  )
 
  library(msigdbr)
  library(fgsea)
  library(dplyr)
  
  # Load Reactome pathways for Homo sapiens
  reactome_sets <- msigdbr(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:REACTOME"
  )
  
  cat("Reactome gene-set rows:", nrow(reactome_sets), "\n")
  cat(
    "Reactome pathways:",
    dplyr::n_distinct(reactome_sets$gs_name),
    "\n"
  )
  
  # Convert to a list format required by fgsea
  reactome_list <- split(
    reactome_sets$entrez_gene,
    reactome_sets$gs_name
  )
  
  # Remove duplicated genes within each pathway
  reactome_list <- lapply(
    reactome_list,
    unique
  )
  
  # Run Reactome GSEA
  set.seed(123)
  
  proteomics_gsea_reactome_fgsea <- fgsea(
    pathways = reactome_list,
    stats = proteomics_gsea_rank,
    minSize = 10,
    maxSize = 500,
    nperm = 10000
  )
  
  # Sort results by adjusted p-value
  reactome_result <- proteomics_gsea_reactome_fgsea %>%
    as.data.frame() %>%
    arrange(padj)
  
  cat(
    "Reactome pathways tested:",
    nrow(reactome_result),
    "\n"
  )
  
  cat(
    "Significant pathways:",
    sum(reactome_result$padj < 0.05, na.rm = TRUE),
    "\n"
  )
  
  # Display the top 20 pathways
  print(
    reactome_result %>%
      select(
        pathway,
        size,
        ES,
        NES,
        pval,
        padj
      ) %>%
      head(20)
  )

  # Prepare Reactome gene sets for fgsea
  
  reactome_sets_clean <- reactome_sets %>%
    dplyr::filter(
      !is.na(ncbi_gene),
      ncbi_gene != ""
    ) %>%
    dplyr::mutate(
      ncbi_gene = as.numeric(ncbi_gene)
    ) %>%
    dplyr::filter(
      !is.na(ncbi_gene)
    )
  
  reactome_list <- split(
    reactome_sets_clean$ncbi_gene,
    reactome_sets_clean$gs_name
  )
  
  # Remove duplicated genes within each pathway
  reactome_list <- lapply(
    reactome_list,
    unique
  )
  
  # Keep pathways within the desired size range
  reactome_list <- reactome_list[
    sapply(reactome_list, length) >= 10 &
      sapply(reactome_list, length) <= 500
  ]
  
  cat(
    "Reactome pathways after filtering:",
    length(reactome_list),
    "\n"
  )
  
  cat(
    "Smallest pathway size:",
    min(sapply(reactome_list, length)),
    "\n"
  )
  
  cat(
    "Largest pathway size:",
    max(sapply(reactome_list, length)),
    "\n"
  )
  
  cat(
    "Genes in ranking:",
    length(proteomics_gsea_rank),
    "\n"
  )
  # Run Reactome GSEA using fgsea
  
  set.seed(123)
  
  proteomics_gsea_reactome_fgsea <- fgsea(
    pathways = reactome_list,
    stats = proteomics_gsea_rank,
    minSize = 10,
    maxSize = 500,
    nperm = 10000
  )
  
  # Convert results to data frame and sort by adjusted p-value
  reactome_result <- proteomics_gsea_reactome_fgsea %>%
    as.data.frame() %>%
    dplyr::arrange(padj)
  
  cat(
    "Reactome pathways tested:",
    nrow(reactome_result),
    "\n"
  )
  
  cat(
    "Significant pathways:",
    sum(reactome_result$padj < 0.05, na.rm = TRUE),
    "\n"
  )
  
  cat("\nTop 20 Reactome pathways:\n")
  
  print(
    reactome_result %>%
      dplyr::select(
        pathway,
        size,
        ES,
        NES,
        pval,
        padj
      ) %>%
      head(20)
  )  
  # Save Reactome GSEA results safely
  
  gsea_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/GSEA"
  
  dir.create(
    gsea_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Save the original fgsea object
  saveRDS(
    proteomics_gsea_reactome_fgsea,
    file.path(
      gsea_dir,
      "Proteomics_GSEA_Reactome_fgsea.rds"
    )
  )
  
  # Convert list-columns to character for CSV export
  reactome_result_csv <- reactome_result %>%
    dplyr::mutate(
      dplyr::across(
        where(is.list),
        ~ vapply(
          .x,
          function(x) paste(x, collapse = ";"),
          character(1)
        )
      )
    )
  
  # Save full results
  write.csv(
    reactome_result_csv,
    file.path(
      gsea_dir,
      "Proteomics_GSEA_Reactome.csv"
    ),
    row.names = FALSE
  )
  
  # Significant pathways
  reactome_significant_csv <- reactome_result_csv %>%
    dplyr::filter(padj < 0.05)
  
  write.csv(
    reactome_significant_csv,
    file.path(
      gsea_dir,
      "Proteomics_GSEA_Reactome_significant.csv"
    ),
    row.names = FALSE
  )
  
  # Top 20 pathways
  reactome_top20_csv <- reactome_result_csv %>%
    dplyr::slice_head(n = 20)
  
  write.csv(
    reactome_top20_csv,
    file.path(
      gsea_dir,
      "Proteomics_GSEA_Reactome_top20.csv"
    ),
    row.names = FALSE
  )
  
  # Validation
  cat(
    "RDS exists:",
    file.exists(
      file.path(gsea_dir, "Proteomics_GSEA_Reactome_fgsea.rds")
    ),
    "\n"
  )
  
  cat(
    "Full CSV exists:",
    file.exists(
      file.path(gsea_dir, "Proteomics_GSEA_Reactome.csv")
    ),
    "\n"
  )
  
  cat(
    "Significant CSV exists:",
    file.exists(
      file.path(gsea_dir, "Proteomics_GSEA_Reactome_significant.csv")
    ),
    "\n"
  )
  
  cat(
    "Significant pathways saved:",
    nrow(reactome_significant_csv),
    "\n"
  )
  
  
  library(ggplot2)
  library(dplyr)
  
  gsea_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/GSEA"
  
  # Select the 14 significant Reactome pathways
  reactome_plot_data <- reactome_result %>%
    filter(padj < 0.05) %>%
    arrange(NES) %>%
    mutate(
      pathway_clean = gsub("REACTOME_", "", pathway),
      pathway_clean = gsub("_", " ", pathway_clean),
      pathway_clean = tools::toTitleCase(tolower(pathway_clean))
    )
  
  # Create a horizontal NES bar plot
  p_reactome <- ggplot(
    reactome_plot_data,
    aes(
      x = NES,
      y = reorder(pathway_clean, NES)
    )
  ) +
    geom_col() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    labs(
      title = "Reactome GSEA — Proteomics",
      subtitle = "Significant pathways (BH-adjusted p < 0.05)",
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    theme_minimal(base_size = 12)
  
  print(p_reactome)
  
  # Save plot
  ggsave(
    file.path(
      gsea_dir,
      "Proteomics_GSEA_Reactome_Significant_NES.png"
    ),
    p_reactome,
    width = 14,
    height = 9,
    dpi = 300
  )
  
  # Save plot object
  saveRDS(
    p_reactome,
    file.path(
      gsea_dir,
      "Proteomics_GSEA_Reactome_Significant_NES.rds"
    )
  )
  
  cat(
    "Reactome plot PNG exists:",
    file.exists(
      file.path(
        gsea_dir,
        "Proteomics_GSEA_Reactome_Significant_NES.png"
      )
    ),
    "\n"
  )
  
  cat(
    "Number of significant pathways plotted:",
    nrow(reactome_plot_data),
    "\n"
  )
  # Inspect available clinical variables for covariate adjustment
  
  olink_clinical_clean <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/Olink_Clinical_Clean_781_samples.rds"
  )
  
  cat("Samples:", nrow(olink_clinical_clean), "\n")
  cat("Columns:", ncol(olink_clinical_clean), "\n\n")
  
  cat("Available clinical variables:\n")
  print(names(olink_clinical_clean))
  
  cat("\nCOVID distribution:\n")
  print(table(olink_clinical_clean$COVID, useNA = "ifany"))
  
  # Identify clinical variables available for covariate adjustment
  
  clinical_cols <- names(olink_clinical_clean)[
    !grepl("^OID[0-9]+$", names(olink_clinical_clean))
  ]
  
  cat("Number of non-protein variables:", length(clinical_cols), "\n\n")
  
  cat("Clinical variables:\n")
  print(clinical_cols)
  
  cat("\nVariable summaries:\n")
  
  for (v in clinical_cols) {
    cat("\n---", v, "---\n")
    print(table(olink_clinical_clean[[v]], useNA = "ifany"))
  } 
 
  
  # Validate clinical covariates before building the adjusted model
  
  candidate_covariates <- c(
    "Age cat",
    "BMI cat",
    "HEART",
    "LUNG",
    "KIDNEY",
    "DIABETES",
    "HTN",
    "IMMUNO"
  )
  
  for (v in candidate_covariates) {
    x <- olink_clinical_clean[[v]]
    
    cat("\n---", v, "---\n")
    cat("Class:", class(x)[1], "\n")
    cat("Missing:", sum(is.na(x)), "\n")
    cat("Unique values:", length(unique(na.omit(x))), "\n")
    print(table(x, useNA = "ifany"))
  }
  # Build and validate the adjusted limma design matrix
  
  d0_adjusted <- olink_clinical_clean %>%
    dplyr::filter(Timepoint == "D0") %>%
    dplyr::filter(COVID %in% c("COVID-", "COVID+")) %>%
    dplyr::mutate(
      COVID = factor(COVID, levels = c("COVID-", "COVID+")),
      `Age cat` = factor(`Age cat`),
      `BMI cat` = factor(`BMI cat`),
      HEART = factor(HEART),
      LUNG = factor(LUNG),
      KIDNEY = factor(KIDNEY),
      DIABETES = factor(DIABETES),
      HTN = factor(HTN),
      IMMUNO = factor(IMMUNO)
    )
  
  cat("D0 samples:", nrow(d0_adjusted), "\n")
  cat("COVID-:", sum(d0_adjusted$COVID == "COVID-"), "\n")
  cat("COVID+:", sum(d0_adjusted$COVID == "COVID+"), "\n")
  
  adjusted_design <- model.matrix(
    ~ COVID +
      `Age cat` +
      `BMI cat` +
      HEART +
      LUNG +
      KIDNEY +
      DIABETES +
      HTN +
      IMMUNO,
    data = d0_adjusted
  )
  
  cat("\nDesign dimensions:\n")
  print(dim(adjusted_design))
  
  cat("\nDesign rank:\n")
  cat(qr(adjusted_design)$rank, "\n")
  
  cat("\nNumber of design columns:\n")
  cat(ncol(adjusted_design), "\n")
  
  cat("\nDesign column names:\n")
  print(colnames(adjusted_design)) 
  cat("Timepoint values:\n")
  print(table(olink_clinical_clean$Timepoint, useNA = "ifany"))
  
  cat("\nCOVID values:\n")
  print(table(olink_clinical_clean$COVID, useNA = "ifany"))
  
  cat("\nFirst 20 Timepoint values:\n")
  print(head(olink_clinical_clean$Timepoint, 20))
  
  cat("\nFirst 20 COVID values:\n")
  print(head(olink_clinical_clean$COVID, 20))
  # Prepare D0 cohort and validate the adjusted design
  
  d0_adjusted <- olink_clinical_clean %>%
    dplyr::filter(Timepoint == "D0") %>%
    dplyr::filter(COVID %in% c(0, 1)) %>%
    dplyr::mutate(
      COVID = factor(COVID, levels = c(0, 1),
                     labels = c("Control", "COVID+")),
      `Age cat` = factor(`Age cat`),
      `BMI cat` = factor(`BMI cat`),
      HEART = factor(HEART),
      LUNG = factor(LUNG),
      KIDNEY = factor(KIDNEY),
      DIABETES = factor(DIABETES),
      HTN = factor(HTN),
      IMMUNO = factor(IMMUNO)
    )
  
  cat("D0 samples:", nrow(d0_adjusted), "\n")
  cat("Control:", sum(d0_adjusted$COVID == "Control"), "\n")
  cat("COVID+:", sum(d0_adjusted$COVID == "COVID+"), "\n")
  
  adjusted_design <- model.matrix(
    ~ COVID +
      `Age cat` +
      `BMI cat` +
      HEART +
      LUNG +
      KIDNEY +
      DIABETES +
      HTN +
      IMMUNO,
    data = d0_adjusted
  )
  
  cat("\nDesign dimensions:\n")
  print(dim(adjusted_design))
  
  cat("\nDesign rank:", qr(adjusted_design)$rank, "\n")
  cat("Number of columns:", ncol(adjusted_design), "\n")
  
  cat("\nRank check:\n")
  if (qr(adjusted_design)$rank == ncol(adjusted_design)) {
    cat("FULL RANK: No complete linear dependency detected.\n")
  } else {
    cat("NOT FULL RANK: Potential collinearity detected.\n")
  }
  
  cat("\nDesign columns:\n")
  print(colnames(adjusted_design))
  
  # Run covariate-adjusted differential proteomics
  
  d0_protein_matrix_adjusted <- t(
    as.matrix(d0_adjusted[, protein_cols, drop = FALSE])
  )
  
  rownames(d0_protein_matrix_adjusted) <- protein_cols
  colnames(d0_protein_matrix_adjusted) <- d0_adjusted$`Public ID`
  
  storage.mode(d0_protein_matrix_adjusted) <- "numeric"
  
  fit_adjusted <- limma::lmFit(
    d0_protein_matrix_adjusted,
    adjusted_design
  )
  
  fit_adjusted <- limma::eBayes(
    fit_adjusted,
    trend = TRUE,
    robust = TRUE
  )
  
  proteomics_adjusted <- limma::topTable(
    fit_adjusted,
    coef = "COVIDCOVID+",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  
  proteomics_adjusted$Protein <- rownames(proteomics_adjusted)
  rownames(proteomics_adjusted) <- NULL
  
  proteomics_adjusted$Regulation <- ifelse(
    proteomics_adjusted$adj.P.Val < 0.05 &
      proteomics_adjusted$logFC >= 0.5,
    "Up",
    ifelse(
      proteomics_adjusted$adj.P.Val < 0.05 &
        proteomics_adjusted$logFC <= -0.5,
      "Down",
      "NS"
    )
  )
  
  proteomics_adjusted$Significant <-
    proteomics_adjusted$Regulation != "NS"
  
  proteomics_adjusted <- proteomics_adjusted[, c(
    "Protein",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B",
    "Regulation",
    "Significant"
  )]
  
  cat("Proteins analyzed:", nrow(proteomics_adjusted), "\n")
  cat("Significant:", sum(proteomics_adjusted$Significant), "\n")
  cat("Up:", sum(proteomics_adjusted$Regulation == "Up"), "\n")
  cat("Down:", sum(proteomics_adjusted$Regulation == "Down"), "\n")
  cat("NS:", sum(proteomics_adjusted$Regulation == "NS"), "\n")
  
  cat("\nTop 20 adjusted proteins:\n")
  print(
    head(
      proteomics_adjusted[
        order(proteomics_adjusted$adj.P.Val),
      ],
      20
    )
  )
  
  
  # Save adjusted differential proteomics results
  
  saveRDS(
    fit_adjusted,
    file.path(
      proteomics_results_path,
      "Differential_Proteomics_D0_COVID_vs_Control_Adjusted_limma_fit.rds"
    )
  )
  
  saveRDS(
    proteomics_adjusted,
    file.path(
      proteomics_results_path,
      "Differential_Proteomics_D0_COVID_vs_Control_Adjusted_all.rds"
    )
  )
  
  write.csv(
    proteomics_adjusted,
    file.path(
      proteomics_results_path,
      "Differential_Proteomics_D0_COVID_vs_Control_Adjusted_all.csv"
    ),
    row.names = FALSE
  )
  
  proteomics_adjusted_sig <- proteomics_adjusted[
    proteomics_adjusted$Significant,
  ]
  
  proteomics_adjusted_up <- proteomics_adjusted[
    proteomics_adjusted$Regulation == "Up",
  ]
  
  proteomics_adjusted_down <- proteomics_adjusted[
    proteomics_adjusted$Regulation == "Down",
  ]
  
  saveRDS(
    proteomics_adjusted_sig,
    file.path(
      proteomics_results_path,
      "Differential_Proteomics_D0_COVID_vs_Control_Adjusted_significant.rds"
    )
  )
  
  write.csv(
    proteomics_adjusted_sig,
    file.path(
      proteomics_results_path,
      "Differential_Proteomics_D0_COVID_vs_Control_Adjusted_significant.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    proteomics_adjusted_up,
    file.path(
      proteomics_results_path,
      "Differential_Proteomics_D0_COVID_vs_Control_Adjusted_upregulated.rds"
    )
  )
  
  write.csv(
    proteomics_adjusted_up,
    file.path(
      proteomics_results_path,
      "Differential_Proteomics_D0_COVID_vs_Control_Adjusted_upregulated.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    proteomics_adjusted_down,
    file.path(
      proteomics_results_path,
      "Differential_Proteomics_D0_COVID_vs_Control_Adjusted_downregulated.rds"
    )
  )
  
  write.csv(
    proteomics_adjusted_down,
    file.path(
      proteomics_results_path,
      "Differential_Proteomics_D0_COVID_vs_Control_Adjusted_downregulated.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    adjusted_design,
    file.path(
      proteomics_results_path,
      "Proteomics_D0_Adjusted_Design_Matrix.rds"
    )
  )
  
  saveRDS(
    d0_adjusted,
    file.path(
      proteomics_results_path,
      "Proteomics_D0_Adjusted_Cohort.rds"
    )
  )
  
  cat("Adjusted differential proteomics results saved successfully.\n")
  cat("All proteins:", nrow(proteomics_adjusted), "\n")
  cat("Significant:", nrow(proteomics_adjusted_sig), "\n")
  cat("Up:", nrow(proteomics_adjusted_up), "\n")
  cat("Down:", nrow(proteomics_adjusted_down), "\n")
  
  # Compare unadjusted and adjusted differential proteomics
  
  # Restore unadjusted results if needed
  if (!exists("proteomics_de")) {
    proteomics_de <- readRDS(
      file.path(
        proteomics_results_path,
        "Differential_Proteomics_D0_COVID_vs_Control_all.rds"
      )
    )
  }
  
  # Add significance definition to unadjusted results
  proteomics_unadjusted <- proteomics_de
  
  proteomics_unadjusted$Regulation <- ifelse(
    proteomics_unadjusted$adj.P.Val < 0.05 &
      proteomics_unadjusted$logFC >= 0.5,
    "Up",
    ifelse(
      proteomics_unadjusted$adj.P.Val < 0.05 &
        proteomics_unadjusted$logFC <= -0.5,
      "Down",
      "NS"
    )
  )
  
  proteomics_unadjusted$Significant <-
    proteomics_unadjusted$Regulation != "NS"
  
  # Prepare comparison table
  comparison <- merge(
    proteomics_unadjusted[, c(
      "Protein",
      "logFC",
      "adj.P.Val",
      "Regulation",
      "Significant"
    )],
    proteomics_adjusted[, c(
      "Protein",
      "logFC",
      "adj.P.Val",
      "Regulation",
      "Significant"
    )],
    by = "Protein",
    suffixes = c("_Unadjusted", "_Adjusted")
  )
  
  # Classify robustness
  comparison$Status <- ifelse(
    comparison$Significant_Unadjusted &
      comparison$Significant_Adjusted,
    "Significant_Both",
    ifelse(
      comparison$Significant_Unadjusted &
        !comparison$Significant_Adjusted,
      "Lost_After_Adjustment",
      ifelse(
        !comparison$Significant_Unadjusted &
          comparison$Significant_Adjusted,
        "New_After_Adjustment",
        "NS_Both"
      )
    )
  )
  
  # Direction consistency
  comparison$Direction_Consistent <- comparison$Regulation_Unadjusted ==
    comparison$Regulation_Adjusted
  
  comparison$Direction_Change <- ifelse(
    comparison$Significant_Unadjusted &
      comparison$Significant_Adjusted &
      !comparison$Direction_Consistent,
    TRUE,
    FALSE
  )
  
  # LogFC change
  comparison$logFC_Change <- comparison$logFC_Adjusted -
    comparison$logFC_Unadjusted
  
  # Summary
  cat("Total proteins compared:", nrow(comparison), "\n\n")
  
  cat("Unadjusted significant:",
      sum(comparison$Significant_Unadjusted), "\n")
  
  cat("Adjusted significant:",
      sum(comparison$Significant_Adjusted), "\n\n")
  
  cat("Significant in both:",
      sum(comparison$Status == "Significant_Both"), "\n")
  
  cat("Lost after adjustment:",
      sum(comparison$Status == "Lost_After_Adjustment"), "\n")
  
  cat("New after adjustment:",
      sum(comparison$Status == "New_After_Adjustment"), "\n")
  
  cat("NS in both:",
      sum(comparison$Status == "NS_Both"), "\n\n")
  
  robust <- comparison[
    comparison$Status == "Significant_Both" &
      comparison$Direction_Consistent,
  ]
  
  cat("Robust significant proteins:",
      nrow(robust), "\n")
  
  cat("\nRobust Up:", sum(
    robust$Regulation_Unadjusted == "Up"
  ), "\n")
  
  cat("Robust Down:", sum(
    robust$Regulation_Unadjusted == "Down"
  ), "\n")
  
  cat("\nDirection changes among significant in both:",
      sum(comparison$Direction_Change), "\n")
  
  cat("\nTop robust proteins by adjusted significance:\n")
  
  robust_top <- robust[
    order(robust$adj.P.Val_Adjusted),
  ]
  
  print(
    head(robust_top, 30)
  )
  
  # Save complete comparison
  saveRDS(
    comparison,
    file.path(
      proteomics_results_path,
      "Proteomics_D0_Unadjusted_vs_Adjusted_Comparison.rds"
    )
  )
  
  write.csv(
    comparison,
    file.path(
      proteomics_results_path,
      "Proteomics_D0_Unadjusted_vs_Adjusted_Comparison.csv"
    ),
    row.names = FALSE
  )
  
  # Save robust proteins
  saveRDS(
    robust,
    file.path(
      proteomics_results_path,
      "Proteomics_D0_Robust_Significant_Proteins.rds"
    )
  )
  
  write.csv(
    robust,
    file.path(
      proteomics_results_path,
      "Proteomics_D0_Robust_Significant_Proteins.csv"
    ),
    row.names = FALSE
  )
  
  # Save status-specific tables
  for (status in c(
    "Significant_Both",
    "Lost_After_Adjustment",
    "New_After_Adjustment",
    "NS_Both"
  )) {
    
    tmp <- comparison[
      comparison$Status == status,
    ]
    
    safe_name <- gsub("_", "_", status)
    
    write.csv(
      tmp,
      file.path(
        proteomics_results_path,
        paste0(
          "Proteomics_D0_",
          safe_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )
  }
  
  cat("\nRobustness comparison saved successfully.\n")
  
  # Annotate the 164 robust proteomic proteins
  
  robust_annotated <- robust %>%
    dplyr::left_join(
      olink_mapping,
      by = c("Protein" = "OlinkID")
    )
  
  cat("Robust proteins:", nrow(robust_annotated), "\n")
  
  cat("Mapped to Olink annotation:",
      sum(!is.na(robust_annotated$Assay)), "\n")
  
  cat("Unmapped:",
      sum(is.na(robust_annotated$Assay)), "\n")
  
  cat("\nPanel distribution:\n")
  print(table(robust_annotated$Panel, useNA = "ifany"))
  
  cat("\nTop 20 robust proteins:\n")
  print(
    robust_annotated %>%
      dplyr::arrange(adj.P.Val_Adjusted) %>%
      dplyr::select(
        Protein,
        Assay,
        UniProt,
        Panel,
        logFC_Unadjusted,
        adj.P.Val_Unadjusted,
        logFC_Adjusted,
        adj.P.Val_Adjusted,
        Regulation_Adjusted
      ) %>%
      head(20)
  )
  
  # Save annotated robust signature
  saveRDS(
    robust_annotated,
    file.path(
      proteomics_results_path,
      "Proteomics_D0_Robust_Significant_Proteins_Annotated.rds"
    )
  )
  
  write.csv(
    robust_annotated,
    file.path(
      proteomics_results_path,
      "Proteomics_D0_Robust_Significant_Proteins_Annotated.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nAnnotated robust signature saved successfully.\n")
  
  
  # GO Biological Process enrichment for robust Upregulated proteins
  
  robust_up <- robust_annotated %>%
    dplyr::filter(Regulation_Adjusted == "Up")
  
  cat("Robust Upregulated proteins:", nrow(robust_up), "\n")
  
  # Map UniProt to Entrez
  up_entrez <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(robust_up$UniProt),
    columns = c("ENTREZID"),
    keytype = "UNIPROT"
  )
  
  up_entrez <- up_entrez %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(UNIPROT, .keep_all = TRUE)
  
  up_gene_ids <- unique(up_entrez$ENTREZID)
  
  cat("Mapped Entrez IDs:", length(up_gene_ids), "\n")
  
  # Background: all Olink proteins with valid Entrez mapping
  background_entrez <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(olink_mapping$UniProt),
    columns = c("ENTREZID"),
    keytype = "UNIPROT"
  ) %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(UNIPROT, .keep_all = TRUE)
  
  background_ids <- unique(background_entrez$ENTREZID)
  
  cat("Background Entrez IDs:", length(background_ids), "\n")
  
  # GO Biological Process enrichment
  go_bp_robust_up <- clusterProfiler::enrichGO(
    gene = up_gene_ids,
    universe = background_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  go_bp_robust_up_df <- as.data.frame(go_bp_robust_up)
  
  cat("\nSignificant GO-BP terms:", nrow(go_bp_robust_up_df), "\n")
  
  if (nrow(go_bp_robust_up_df) > 0) {
    print(
      go_bp_robust_up_df %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::select(
          ID,
          Description,
          GeneRatio,
          BgRatio,
          pvalue,
          p.adjust,
          qvalue,
          Count
        ) %>%
        head(20)
    )
  }
  
  # Save results
  dir.create(
    file.path(proteomics_results_path, "Enrichment", "Robust"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    go_bp_robust_up,
    file.path(
      proteomics_results_path,
      "Enrichment",
      "Robust",
      "GO_BP_Robust_Upregulated.rds"
    )
  )
  
  write.csv(
    go_bp_robust_up_df,
    file.path(
      proteomics_results_path,
      "Enrichment",
      "Robust",
      "GO_BP_Robust_Upregulated.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nGO-BP robust Upregulated enrichment saved successfully.\n")
  
  # GO Biological Process enrichment for robust Downregulated proteins
  
  robust_down <- robust_annotated %>%
    dplyr::filter(Regulation_Adjusted == "Down")
  
  cat("Robust Downregulated proteins:", nrow(robust_down), "\n")
  
  # Map UniProt to Entrez
  down_entrez <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(robust_down$UniProt),
    columns = c("ENTREZID"),
    keytype = "UNIPROT"
  )
  
  down_entrez <- down_entrez %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(UNIPROT, .keep_all = TRUE)
  
  down_gene_ids <- unique(down_entrez$ENTREZID)
  
  cat("Mapped Entrez IDs:", length(down_gene_ids), "\n")
  
  # GO Biological Process enrichment
  go_bp_robust_down <- clusterProfiler::enrichGO(
    gene = down_gene_ids,
    universe = background_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  go_bp_robust_down_df <- as.data.frame(go_bp_robust_down)
  
  cat("\nSignificant GO-BP terms:", nrow(go_bp_robust_down_df), "\n")
  
  if (nrow(go_bp_robust_down_df) > 0) {
    print(
      go_bp_robust_down_df %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::select(
          ID,
          Description,
          GeneRatio,
          BgRatio,
          pvalue,
          p.adjust,
          qvalue,
          Count
        ) %>%
        head(20)
    )
  }
  
  # Save results
  saveRDS(
    go_bp_robust_down,
    file.path(
      proteomics_results_path,
      "Enrichment",
      "Robust",
      "GO_BP_Robust_Downregulated.rds"
    )
  )
  
  write.csv(
    go_bp_robust_down_df,
    file.path(
      proteomics_results_path,
      "Enrichment",
      "Robust",
      "GO_BP_Robust_Downregulated.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nGO-BP robust Downregulated enrichment saved successfully.\n")  

  
  # KEGG enrichment for robust Upregulated proteins
  
  kegg_robust_up <- clusterProfiler::enrichKEGG(
    gene = up_gene_ids,
    universe = background_ids,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  kegg_robust_up_df <- as.data.frame(kegg_robust_up)
  
  cat("Robust Upregulated proteins:", nrow(robust_up), "\n")
  cat("Mapped Entrez IDs:", length(up_gene_ids), "\n")
  cat("Significant KEGG pathways:", nrow(kegg_robust_up_df), "\n")
  
  if (nrow(kegg_robust_up_df) > 0) {
    print(
      kegg_robust_up_df %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::select(
          ID,
          Description,
          GeneRatio,
          BgRatio,
          pvalue,
          p.adjust,
          qvalue,
          Count
        ) %>%
        head(20)
    )
  }
  
  # Save KEGG results
  saveRDS(
    kegg_robust_up,
    file.path(
      proteomics_results_path,
      "Enrichment",
      "Robust",
      "KEGG_Robust_Upregulated.rds"
    )
  )
  
  write.csv(
    kegg_robust_up_df,
    file.path(
      proteomics_results_path,
      "Enrichment",
      "Robust",
      "KEGG_Robust_Upregulated.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nKEGG robust Upregulated enrichment saved successfully.\n")

  # Prepare ranked proteomics signal for robust GSEA
  
  # Restore objects if needed
  if (!exists("proteomics_de")) {
    proteomics_de <- readRDS(
      file.path(
        proteomics_results_path,
        "Differential_Proteomics_D0_COVID_vs_Control_all.rds"
      )
    )
  }
  
  if (!exists("olink_mapping")) {
    olink_mapping <- readRDS(
      file.path(
        proteomics_results_path,
        "Olink_Protein_Annotation_1429.rds"
      )
    )
  }
  
  # Join differential statistics with protein annotation
  proteomics_gsea_df <- proteomics_de %>%
    dplyr::left_join(
      olink_mapping %>%
        dplyr::select(OlinkID, Assay, UniProt, Panel),
      by = c("Protein" = "OlinkID")
    ) %>%
    dplyr::filter(
      !is.na(UniProt),
      !is.na(t),
      is.finite(t)
    )
  
  cat("Proteins with valid UniProt and t-statistic:",
      nrow(proteomics_gsea_df), "\n")
  
  # Convert UniProt to Entrez
  gsea_entrez_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(proteomics_gsea_df$UniProt),
    columns = "ENTREZID",
    keytype = "UNIPROT"
  )
  
  gsea_entrez_map <- gsea_entrez_map %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::distinct(UNIPROT, .keep_all = TRUE)
  
  # Add Entrez IDs
  proteomics_gsea_df <- proteomics_gsea_df %>%
    dplyr::left_join(
      gsea_entrez_map,
      by = "UniProt"
    ) %>%
    dplyr::filter(!is.na(ENTREZID))
  
  # Keep one protein per Entrez ID
  proteomics_gsea_df <- proteomics_gsea_df %>%
    dplyr::arrange(dplyr::desc(abs(t))) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
  
  # Create ranked vector
  proteomics_gsea_ranked <- proteomics_gsea_df$t
  names(proteomics_gsea_ranked) <- proteomics_gsea_df$ENTREZID
  
  proteomics_gsea_ranked <- sort(
    proteomics_gsea_ranked,
    decreasing = TRUE
  )
  
  cat("Final ranked Entrez proteins:",
      length(proteomics_gsea_ranked), "\n")
  
  cat("Maximum t-statistic:",
      max(proteomics_gsea_ranked), "\n")
  
  cat("Minimum t-statistic:",
      min(proteomics_gsea_ranked), "\n")
  
  cat("Missing values:",
      sum(is.na(proteomics_gsea_ranked)), "\n")
  
  cat("Duplicate Entrez IDs:",
      sum(duplicated(names(proteomics_gsea_ranked))), "\n")
  
  cat("\nTop 15 ranked proteins:\n")
  print(
    proteomics_gsea_df %>%
      dplyr::arrange(dplyr::desc(t)) %>%
      dplyr::select(
        Protein,
        Assay,
        UniProt,
        ENTREZID,
        Panel,
        t,
        logFC,
        adj.P.Val
      ) %>%
      head(15)
  )
  
  cat("\nBottom 15 ranked proteins:\n")
  print(
    proteomics_gsea_df %>%
      dplyr::arrange(t) %>%
      dplyr::select(
        Protein,
        Assay,
        UniProt,
        ENTREZID,
        Panel,
        t,
        logFC,
        adj.P.Val
      ) %>%
      head(15)
  )
  
  # Save ranked GSEA input
  saveRDS(
    proteomics_gsea_ranked,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_Robust_GSEA_Ranking_Entrez.rds"
    )
  )
  
  write.csv(
    proteomics_gsea_df,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_Robust_GSEA_Ranking_Annotation.csv"
    ),
    row.names = FALSE
  )
  
  dir.create(
    file.path(proteomics_results_path, "GSEA"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  cat("\nRobust GSEA ranking prepared and saved successfully.\n")    
  
  # Fix UniProt-to-Entrez mapping and prepare the ranked GSEA vector
  
  # Rename the annotation key to match the proteomics data
  gsea_entrez_map <- gsea_entrez_map %>%
    dplyr::rename(UniProt = UNIPROT)
  
  cat("UniProt mapping rows:", nrow(gsea_entrez_map), "\n")
  cat("Unique UniProt IDs:", dplyr::n_distinct(gsea_entrez_map$UniProt), "\n")
  
  # Add Entrez IDs
  proteomics_gsea_df <- proteomics_gsea_df %>%
    dplyr::left_join(
      gsea_entrez_map,
      by = "UniProt"
    ) %>%
    dplyr::filter(!is.na(ENTREZID))
  
  cat("Proteins mapped to Entrez:", nrow(proteomics_gsea_df), "\n")
  
  # Keep one protein per Entrez ID
  proteomics_gsea_df <- proteomics_gsea_df %>%
    dplyr::arrange(dplyr::desc(abs(t))) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
  
  # Create ranked vector
  proteomics_gsea_ranked <- proteomics_gsea_df$t
  names(proteomics_gsea_ranked) <- proteomics_gsea_df$ENTREZID
  
  proteomics_gsea_ranked <- sort(
    proteomics_gsea_ranked,
    decreasing = TRUE
  )
  
  cat("Final ranked Entrez proteins:",
      length(proteomics_gsea_ranked), "\n")
  
  cat("Maximum t-statistic:",
      max(proteomics_gsea_ranked), "\n")
  
  cat("Minimum t-statistic:",
      min(proteomics_gsea_ranked), "\n")
  
  cat("Missing values:",
      sum(is.na(proteomics_gsea_ranked)), "\n")
  
  cat("Duplicate Entrez IDs:",
      sum(duplicated(names(proteomics_gsea_ranked))), "\n")
  
  # Create output directory before saving
  dir.create(
    file.path(proteomics_results_path, "GSEA"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # Save ranked GSEA input
  saveRDS(
    proteomics_gsea_ranked,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_Robust_GSEA_Ranking_Entrez.rds"
    )
  )
  
  write.csv(
    proteomics_gsea_df,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_Robust_GSEA_Ranking_Annotation.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nRobust GSEA ranking prepared and saved successfully.\n")  

  # Robust proteomics GSEA using GO Biological Process gene sets
  
  library(msigdbr)
  library(fgsea)
  
  # Retrieve GO Biological Process gene sets
  go_bp_sets <- msigdbr(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:BP"
  )
  
  cat("GO-BP rows:", nrow(go_bp_sets), "\n")
  cat("GO-BP pathways:", dplyr::n_distinct(go_bp_sets$gs_name), "\n")
  
  # Prepare Entrez-based gene sets
  go_bp_list <- split(
    go_bp_sets$ncbi_gene,
    go_bp_sets$gs_name
  )
  
  go_bp_list <- lapply(
    go_bp_list,
    function(x) unique(as.character(x[!is.na(x)]))
  )
  
  # Keep pathways with reasonable size
  go_bp_list <- go_bp_list[
    vapply(go_bp_list, length, integer(1)) >= 10 &
      vapply(go_bp_list, length, integer(1)) <= 500
  ]
  
  cat("GO-BP pathways after size filtering:",
      length(go_bp_list), "\n")
  
  # Run fgsea
  set.seed(123)
  
  gsea_go_bp_robust <- fgsea::fgsea(
    pathways = go_bp_list,
    stats = proteomics_gsea_ranked,
    minSize = 10,
    maxSize = 500,
    nperm = 10000
  )
  
  # Order by adjusted p-value
  gsea_go_bp_robust <- gsea_go_bp_robust[
    order(gsea_go_bp_robust$padj),
  ]
  
  gsea_go_bp_robust_df <- as.data.frame(
    gsea_go_bp_robust
  )
  
  # Significant pathways
  gsea_go_bp_robust_sig <- gsea_go_bp_robust_df %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::arrange(padj)
  
  cat("\nGO-BP pathways tested:",
      nrow(gsea_go_bp_robust_df), "\n")
  
  cat("Significant GO-BP pathways:",
      nrow(gsea_go_bp_robust_sig), "\n")
  
  if (nrow(gsea_go_bp_robust_sig) > 0) {
    print(
      gsea_go_bp_robust_sig %>%
        dplyr::select(
          pathway,
          NES,
          pval,
          padj,
          size
        ) %>%
        head(20)
    )
  }
  
  # Save complete GSEA results
  saveRDS(
    gsea_go_bp_robust,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_Robust_GSEA_GO_BP_fgsea.rds"
    )
  )
  
  write.csv(
    gsea_go_bp_robust_df,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_Robust_GSEA_GO_BP_fgsea.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    gsea_go_bp_robust_sig,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_Robust_GSEA_GO_BP_significant.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nRobust proteomics GO-BP GSEA saved successfully.\n") 
  # Re-run robust proteomics GO-BP GSEA using fgseaMultilevel
  
  set.seed(123)
  
  gsea_go_bp_robust <- fgsea::fgseaMultilevel(
    pathways = go_bp_list,
    stats = proteomics_gsea_ranked,
    minSize = 10,
    maxSize = 500
  )
  
  # Convert to data.frame
  gsea_go_bp_robust_df <- as.data.frame(
    gsea_go_bp_robust
  )
  
  # Remove list-columns so the result can be exported safely
  list_columns <- vapply(
    gsea_go_bp_robust_df,
    is.list,
    logical(1)
  )
  
  cat("List-columns detected:",
      sum(list_columns), "\n")
  
  if (any(list_columns)) {
    cat("Removing list-columns:",
        paste(
          names(gsea_go_bp_robust_df)[list_columns],
          collapse = ", "
        ),
        "\n")
    
    gsea_go_bp_robust_df <- gsea_go_bp_robust_df[
      ,
      !list_columns,
      drop = FALSE
    ]
  }
  
  # Significant pathways
  gsea_go_bp_robust_sig <- gsea_go_bp_robust_df %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::arrange(padj)
  
  cat("\nGO-BP pathways tested:",
      nrow(gsea_go_bp_robust_df), "\n")
  
  cat("Significant GO-BP pathways:",
      nrow(gsea_go_bp_robust_sig), "\n")
  
  if (nrow(gsea_go_bp_robust_sig) > 0) {
    
    cat("\nTop significant GO-BP pathways:\n")
    
    print(
      gsea_go_bp_robust_sig %>%
        dplyr::select(
          pathway,
          NES,
          pval,
          padj,
          size
        ) %>%
        head(20)
    )
    
  } else {
    
    cat("\nNo GO-BP pathways reached padj < 0.05.\n")
    
    cat("\nTop pathways by adjusted p-value:\n")
    
    print(
      gsea_go_bp_robust_df %>%
        dplyr::arrange(padj) %>%
        dplyr::select(
          pathway,
          NES,
          pval,
          padj,
          size
        ) %>%
        head(20)
    )
  }
  
  # Save results
  dir.create(
    file.path(proteomics_results_path, "GSEA"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    gsea_go_bp_robust_df,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_Robust_GSEA_GO_BP_fgseaMultilevel.rds"
    )
  )
  
  write.csv(
    gsea_go_bp_robust_df,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_Robust_GSEA_GO_BP_fgseaMultilevel.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    gsea_go_bp_robust_sig,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_Robust_GSEA_GO_BP_significant.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nRobust GO-BP GSEA completed and saved successfully.\n") 
  
  # Build the truly covariate-adjusted GSEA ranking
  # using the t-statistics from the adjusted limma model
  
  cat("Checking adjusted limma model...\n")
  
  cat("Coefficient names:\n")
  print(colnames(fit_adjusted$coefficients))
  
  # Extract the adjusted-model t-statistics for COVID+
  robust_t <- fit_adjusted$t[, "COVIDCOVID+"]
  
  cat("\nAdjusted t-statistics extracted:", length(robust_t), "\n")
  cat("NA values:", sum(is.na(robust_t)), "\n")
  cat("Maximum t:", max(robust_t, na.rm = TRUE), "\n")
  cat("Minimum t:", min(robust_t, na.rm = TRUE), "\n")
  
  # Match protein IDs
  robust_ranking_df <- data.frame(
    OlinkID = rownames(fit_adjusted$t),
    t = as.numeric(robust_t),
    stringsAsFactors = FALSE
  )
  
  # Add protein annotation
  robust_ranking_df <- robust_ranking_df %>%
    dplyr::left_join(
      olink_mapping,
      by = c("OlinkID" = "OlinkID")
    )
  
  cat("\nProteins before Entrez mapping:",
      nrow(robust_ranking_df), "\n")
  
  cat("Proteins with UniProt:",
      sum(!is.na(robust_ranking_df$UniProt) &
            robust_ranking_df$UniProt != ""),
      "\n")
  
  # UniProt -> Entrez mapping
  robust_entrez_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(
      robust_ranking_df$UniProt[
        !is.na(robust_ranking_df$UniProt) &
          robust_ranking_df$UniProt != ""
      ]
    ),
    columns = c("UNIPROT", "ENTREZID"),
    keytype = "UNIPROT"
  )
  
  robust_entrez_map <- robust_entrez_map %>%
    dplyr::filter(
      !is.na(ENTREZID),
      !is.na(UNIPROT)
    ) %>%
    dplyr::distinct(UNIPROT, ENTREZID)
  
  # Join Entrez IDs
  robust_gsea_df <- robust_ranking_df %>%
    dplyr::left_join(
      robust_entrez_map,
      by = c("UniProt" = "UNIPROT")
    ) %>%
    dplyr::filter(
      !is.na(ENTREZID),
      !is.na(t)
    ) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)
  
  # Create named ranking
  robust_gsea_ranked <- robust_gsea_df$t
  names(robust_gsea_ranked) <- as.character(robust_gsea_df$ENTREZID)
  
  robust_gsea_ranked <- sort(
    robust_gsea_ranked,
    decreasing = TRUE
  )
  
  cat("\nFinal robust GSEA ranking:\n")
  cat("Proteins:", length(robust_gsea_ranked), "\n")
  cat("NA values:", sum(is.na(robust_gsea_ranked)), "\n")
  cat("Duplicate Entrez IDs:",
      sum(duplicated(names(robust_gsea_ranked))), "\n")
  cat("Maximum t:", max(robust_gsea_ranked), "\n")
  cat("Minimum t:", min(robust_gsea_ranked), "\n")
  
  cat("\nTop 10 proteins by adjusted t-statistic:\n")
  print(
    robust_gsea_df %>%
      dplyr::arrange(dplyr::desc(t)) %>%
      dplyr::select(
        OlinkID,
        Assay,
        UniProt,
        ENTREZID,
        t,
        Panel
      ) %>%
      head(10)
  )
  
  cat("\nBottom 10 proteins by adjusted t-statistic:\n")
  print(
    robust_gsea_df %>%
      dplyr::arrange(t) %>%
      dplyr::select(
        OlinkID,
        Assay,
        UniProt,
        ENTREZID,
        t,
        Panel
      ) %>%
      head(10)
  )
  
  # Save the corrected robust GSEA ranking
  dir.create(
    file.path(proteomics_results_path, "GSEA"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    robust_gsea_ranked,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Ranking_Adjusted_t.rds"
    )
  )
  
  write.csv(
    robust_gsea_df,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Ranking_Adjusted_t_Annotation.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nCorrected robust GSEA ranking saved successfully.\n") 
  # Run GO-BP GSEA using the truly covariate-adjusted ranking
  
  library(msigdbr)
  library(fgsea)
  library(dplyr)
  
  cat("Preparing GO-BP gene sets...\n")
  
  go_bp_msigdbr <- msigdbr(
    species = "Homo sapiens",
    collection = "C5",
    subcollection = "GO:BP"
  )
  
  go_bp_list <- split(
    go_bp_msigdbr$ncbi_gene[
      !is.na(go_bp_msigdbr$ncbi_gene)
    ],
    go_bp_msigdbr$gs_name[
      !is.na(go_bp_msigdbr$ncbi_gene)
    ]
  )
  
  go_bp_list <- lapply(
    go_bp_list,
    function(x) unique(as.character(x))
  )
  
  go_bp_list <- go_bp_list[
    sapply(go_bp_list, length) >= 10 &
      sapply(go_bp_list, length) <= 500
  ]
  
  cat("GO-BP pathways prepared:",
      length(go_bp_list), "\n")
  
  cat("Running fgseaMultilevel...\n")
  
  gsea_go_bp_true_robust <- fgsea::fgseaMultilevel(
    pathways = go_bp_list,
    stats = robust_gsea_ranked,
    minSize = 10,
    maxSize = 500
  )
  
  gsea_go_bp_true_robust <- as.data.frame(
    gsea_go_bp_true_robust
  )
  
  gsea_go_bp_true_robust <- gsea_go_bp_true_robust %>%
    dplyr::arrange(padj)
  
  cat("\nGO-BP pathways tested:",
      nrow(gsea_go_bp_true_robust), "\n")
  
  cat("Significant pathways (padj < 0.05):",
      sum(gsea_go_bp_true_robust$padj < 0.05),
      "\n")
  
  cat("\nTop 20 GO-BP pathways:\n")
  
  print(
    gsea_go_bp_true_robust %>%
      dplyr::select(
        pathway,
        NES,
        pval,
        padj,
        size
      ) %>%
      head(20)
  )
  
  # Save the complete result including leadingEdge
  saveRDS(
    gsea_go_bp_true_robust,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_fgseaMultilevel.rds"
    )
  )
  
  # Save CSV without the list-column
  gsea_go_bp_true_robust_csv <- gsea_go_bp_true_robust %>%
    dplyr::select(
      -dplyr::any_of("leadingEdge")
    )
  
  write.csv(
    gsea_go_bp_true_robust_csv,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_fgseaMultilevel.csv"
    ),
    row.names = FALSE
  )
  
  # Save significant pathways
  gsea_go_bp_true_robust_sig <- gsea_go_bp_true_robust %>%
    dplyr::filter(padj < 0.05)
  
  write.csv(
    gsea_go_bp_true_robust_sig %>%
      dplyr::select(
        -dplyr::any_of("leadingEdge")
      ),
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_significant.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nTRUE robust GO-BP GSEA completed and saved successfully.\n")
  
  
  # Extract leading-edge proteins from significant robust GO-BP pathways
  
  gsea_go_bp_true_robust_sig <- gsea_go_bp_true_robust %>%
    dplyr::filter(padj < 0.05)
  
  cat("Significant GO-BP pathways:", 
      nrow(gsea_go_bp_true_robust_sig), "\n")
  
  # Extract leading-edge Entrez IDs
  leading_edge_list <- gsea_go_bp_true_robust_sig$leadingEdge
  
  leading_edge_df <- purrr::map2_dfr(
    leading_edge_list,
    gsea_go_bp_true_robust_sig$pathway,
    ~data.frame(
      pathway = .y,
      ENTREZID = as.character(.x),
      stringsAsFactors = FALSE
    )
  )
  
  cat("Leading-edge rows:", 
      nrow(leading_edge_df), "\n")
  
  cat("Unique leading-edge proteins:", 
      dplyr::n_distinct(leading_edge_df$ENTREZID), 
      "\n")
  
  # Count recurrence across significant pathways
  leading_edge_core <- leading_edge_df %>%
    dplyr::count(
      ENTREZID,
      name = "pathway_count"
    ) %>%
    dplyr::arrange(
      dplyr::desc(pathway_count)
    )
  
  cat("\nTop recurrent leading-edge proteins:\n")
  
  print(
    leading_edge_core %>%
      head(30)
  )
  
  # Add protein annotation
  leading_edge_core_annotated <- leading_edge_core %>%
    dplyr::left_join(
      robust_gsea_df %>%
        dplyr::select(
          OlinkID,
          Assay,
          UniProt,
          ENTREZID,
          t,
          Panel
        ),
      by = "ENTREZID"
    ) %>%
    dplyr::arrange(
      dplyr::desc(pathway_count),
      dplyr::desc(t)
    )
  
  cat("\nTop recurrent leading-edge proteins with annotation:\n")
  
  print(
    leading_edge_core_annotated %>%
      head(30)
  )
  
  # Save detailed leading-edge table
  saveRDS(
    leading_edge_df,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_LeadingEdge.rds"
    )
  )
  
  write.csv(
    leading_edge_df,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_LeadingEdge.csv"
    ),
    row.names = FALSE
  )
  
  # Save recurrent/core proteins
  saveRDS(
    leading_edge_core_annotated,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_Core_LeadingEdge_Proteins.rds"
    )
  )
  
  write.csv(
    leading_edge_core_annotated,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_Core_LeadingEdge_Proteins.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nTRUE robust GO-BP leading-edge analysis completed and saved.\n")  
  # Separate significant GO-BP pathways by enrichment direction
  
  gsea_go_bp_direction <- gsea_go_bp_true_robust_sig %>%
    dplyr::mutate(
      Direction = dplyr::case_when(
        NES > 0 ~ "Positive",
        NES < 0 ~ "Negative",
        TRUE ~ "Neutral"
      )
    )
  
  cat("Positive NES pathways:",
      sum(gsea_go_bp_direction$Direction == "Positive"),
      "\n")
  
  cat("Negative NES pathways:",
      sum(gsea_go_bp_direction$Direction == "Negative"),
      "\n")
  
  cat("\nPositive NES pathways:\n")
  
  print(
    gsea_go_bp_direction %>%
      dplyr::filter(Direction == "Positive") %>%
      dplyr::select(
        pathway,
        NES,
        pval,
        padj,
        size
      ) %>%
      dplyr::arrange(dplyr::desc(NES))
  )
  
  cat("\nNegative NES pathways:\n")
  
  print(
    gsea_go_bp_direction %>%
      dplyr::filter(Direction == "Negative") %>%
      dplyr::select(
        pathway,
        NES,
        pval,
        padj,
        size
      ) %>%
      dplyr::arrange(NES)
  )
  
  # Build direction-specific leading-edge tables
  
  positive_pathways <- gsea_go_bp_direction %>%
    dplyr::filter(Direction == "Positive")
  
  negative_pathways <- gsea_go_bp_direction %>%
    dplyr::filter(Direction == "Negative")
  
  positive_leading_edge <- purrr::map2_dfr(
    positive_pathways$leadingEdge,
    positive_pathways$pathway,
    ~data.frame(
      pathway = .y,
      ENTREZID = as.character(.x),
      stringsAsFactors = FALSE
    )
  )
  
  negative_leading_edge <- purrr::map2_dfr(
    negative_pathways$leadingEdge,
    negative_pathways$pathway,
    ~data.frame(
      pathway = .y,
      ENTREZID = as.character(.x),
      stringsAsFactors = FALSE
    )
  )
  
  positive_core <- positive_leading_edge %>%
    dplyr::count(
      ENTREZID,
      name = "pathway_count"
    ) %>%
    dplyr::left_join(
      robust_gsea_df %>%
        dplyr::select(
          OlinkID,
          Assay,
          UniProt,
          ENTREZID,
          t,
          Panel
        ),
      by = "ENTREZID"
    ) %>%
    dplyr::arrange(
      dplyr::desc(pathway_count),
      dplyr::desc(t)
    )
  
  negative_core <- negative_leading_edge %>%
    dplyr::count(
      ENTREZID,
      name = "pathway_count"
    ) %>%
    dplyr::left_join(
      robust_gsea_df %>%
        dplyr::select(
          OlinkID,
          Assay,
          UniProt,
          ENTREZID,
          t,
          Panel
        ),
      by = "ENTREZID"
    ) %>%
    dplyr::arrange(
      dplyr::desc(pathway_count),
      t
    )
  
  cat("\nTop positive-NES core proteins:\n")
  print(head(positive_core, 20))
  
  cat("\nTop negative-NES core proteins:\n")
  print(head(negative_core, 20))
  
  # Save direction-specific results
  
  saveRDS(
    positive_core,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_Positive_Core.rds"
    )
  )
  
  write.csv(
    positive_core,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_Positive_Core.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    negative_core,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_Negative_Core.rds"
    )
  )
  
  write.csv(
    negative_core,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_GO_BP_Negative_Core.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nDirection-specific GO-BP core analysis completed and saved.\n")  
  
  # Run TRUE robust Reactome GSEA using the covariate-adjusted limma t-statistics
  
  library(msigdbr)
  library(fgsea)
  library(dplyr)
  
  cat("Preparing Reactome gene sets...\n")
  
  reactome_msigdbr <- msigdbr(
    species = "Homo sapiens",
    collection = "C2",
    subcollection = "CP:REACTOME"
  )
  
  reactome_list <- split(
    reactome_msigdbr$ncbi_gene[
      !is.na(reactome_msigdbr$ncbi_gene)
    ],
    reactome_msigdbr$gs_name[
      !is.na(reactome_msigdbr$ncbi_gene)
    ]
  )
  
  reactome_list <- lapply(
    reactome_list,
    function(x) unique(as.character(x))
  )
  
  reactome_list <- reactome_list[
    sapply(reactome_list, length) >= 10 &
      sapply(reactome_list, length) <= 500
  ]
  
  cat("Reactome pathways prepared:",
      length(reactome_list), "\n")
  
  cat("Running fgseaMultilevel...\n")
  
  gsea_reactome_true_robust <- fgsea::fgseaMultilevel(
    pathways = reactome_list,
    stats = robust_gsea_ranked,
    minSize = 10,
    maxSize = 500
  )
  
  gsea_reactome_true_robust <- as.data.frame(
    gsea_reactome_true_robust
  ) %>%
    dplyr::arrange(padj)
  
  cat("\nReactome pathways tested:",
      nrow(gsea_reactome_true_robust), "\n")
  
  cat("Significant Reactome pathways (padj < 0.05):",
      sum(gsea_reactome_true_robust$padj < 0.05),
      "\n")
  
  cat("\nTop 20 Reactome pathways:\n")
  
  print(
    gsea_reactome_true_robust %>%
      dplyr::select(
        pathway,
        NES,
        pval,
        padj,
        size
      ) %>%
      head(20)
  )
  
  # Save complete result including leadingEdge
  saveRDS(
    gsea_reactome_true_robust,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_fgseaMultilevel.rds"
    )
  )
  
  # Save CSV without list-column
  gsea_reactome_true_robust_csv <- gsea_reactome_true_robust %>%
    dplyr::select(
      -dplyr::any_of("leadingEdge")
    )
  
  write.csv(
    gsea_reactome_true_robust_csv,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_fgseaMultilevel.csv"
    ),
    row.names = FALSE
  )
  
  # Save significant pathways
  gsea_reactome_true_robust_sig <- gsea_reactome_true_robust %>%
    dplyr::filter(padj < 0.05)
  
  write.csv(
    gsea_reactome_true_robust_sig %>%
      dplyr::select(
        -dplyr::any_of("leadingEdge")
      ),
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_significant.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nTRUE robust Reactome GSEA completed and saved successfully.\n")
  
  
  library(dplyr)
  
  # Get significant Reactome pathways
  reactome_sig <- gsea_reactome_true_robust %>%
    filter(padj < 0.05)
  
  cat("Significant Reactome pathways:", nrow(reactome_sig), "\n")
  
  # Separate pathways by direction
  reactome_positive <- reactome_sig %>%
    filter(NES > 0)
  
  reactome_negative <- reactome_sig %>%
    filter(NES < 0)
  
  cat("Positive NES pathways:", nrow(reactome_positive), "\n")
  cat("Negative NES pathways:", nrow(reactome_negative), "\n")
  
  # Extract leading-edge genes
  reactome_positive_leading <- reactome_positive %>%
    select(pathway, NES, padj, leadingEdge) %>%
    tidyr::unnest_longer(leadingEdge) %>%
    rename(Entrez = leadingEdge)
  
  reactome_negative_leading <- reactome_negative %>%
    select(pathway, NES, padj, leadingEdge) %>%
    tidyr::unnest_longer(leadingEdge) %>%
    rename(Entrez = leadingEdge)
  
  # Count recurrence across pathways
  positive_core <- reactome_positive_leading %>%
    count(Entrez, name = "pathway_count") %>%
    arrange(desc(pathway_count))
  
  negative_core <- reactome_negative_leading %>%
    count(Entrez, name = "pathway_count") %>%
    arrange(desc(pathway_count))
  
  cat("\nPositive leading-edge proteins:",
      nrow(positive_core), "\n")
  
  cat("Negative leading-edge proteins:",
      nrow(negative_core), "\n")
  
  cat("\nTop positive core proteins:\n")
  print(head(positive_core, 20))
  
  cat("\nTop negative core proteins:\n")
  print(head(negative_core, 20))
  
  # Save complete leading-edge tables
  saveRDS(
    reactome_positive_leading,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_Positive_LeadingEdge.rds"
    )
  )
  
  write.csv(
    reactome_positive_leading,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_Positive_LeadingEdge.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    reactome_negative_leading,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_Negative_LeadingEdge.rds"
    )
  )
  
  write.csv(
    reactome_negative_leading,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_Negative_LeadingEdge.csv"
    ),
    row.names = FALSE
  )
  
  # Save positive and negative core recurrence
  saveRDS(
    positive_core,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_Positive_Core.rds"
    )
  )
  
  write.csv(
    positive_core,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_Positive_Core.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    negative_core,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_Negative_Core.rds"
    )
  )
  
  write.csv(
    negative_core,
    file.path(
      proteomics_results_path,
      "GSEA",
      "Proteomics_TRUE_Robust_GSEA_Reactome_Negative_Core.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nReactome TRUE Robust leading-edge analysis completed and saved successfully.\n")
  
  # Extract leading-edge genes using explicit dplyr namespace
  
  reactome_positive_leading <- reactome_positive %>%
    dplyr::select(pathway, NES, padj, leadingEdge) %>%
    tidyr::unnest_longer(leadingEdge) %>%
    dplyr::rename(Entrez = leadingEdge)
  
  reactome_negative_leading <- reactome_negative %>%
    dplyr::select(pathway, NES, padj, leadingEdge) %>%
    tidyr::unnest_longer(leadingEdge) %>%
    dplyr::rename(Entrez = leadingEdge)
  
  # Count recurrence across pathways
  positive_core <- reactome_positive_leading %>%
    dplyr::count(Entrez, name = "pathway_count") %>%
    dplyr::arrange(dplyr::desc(pathway_count))
  
  negative_core <- reactome_negative_leading %>%
    dplyr::count(Entrez, name = "pathway_count") %>%
    dplyr::arrange(dplyr::desc(pathway_count))
  
  cat("\nPositive leading-edge proteins:", nrow(positive_core), "\n")
  cat("Negative leading-edge proteins:", nrow(negative_core), "\n")
  
  cat("\nTop positive core proteins:\n")
  print(head(positive_core, 20))
  
  cat("\nTop negative core proteins:\n")
  print(head(negative_core, 20))
  
  library(dplyr)
  
  # Convert Reactome core Entrez IDs to character
  positive_core_annotated <- positive_core %>%
    dplyr::left_join(
      robust_gsea_ranked_annotation %>%
        dplyr::select(
          Entrez,
          OlinkID,
          Assay,
          UniProt,
          t
        ),
      by = "Entrez"
    ) %>%
    dplyr::arrange(desc(pathway_count), desc(t))
  
  negative_core_annotated <- negative_core %>%
    dplyr::left_join(
      robust_gsea_ranked_annotation %>%
        dplyr::select(
          Entrez,
          OlinkID,
          Assay,
          UniProt,
          t
        ),
      by = "Entrez"
    ) %>%
    dplyr::arrange(desc(pathway_count), t)
  
  cat("\nPositive core annotation:\n")
  print(head(positive_core_annotated, 20))
  
  cat("\nNegative core annotation:\n")
  print(head(negative_core_annotated, 20))
  
  cat(
    "\nPositive unmapped:",
    sum(is.na(positive_core_annotated$Assay)),
    "\n"
  )
  
  cat(
    "Negative unmapped:",
    sum(is.na(negative_core_annotated$Assay)),
    "\n"
  )  
  library(dplyr)
  
  # Build annotation for the TRUE robust adjusted t-statistic ranking
  
  robust_rank_annotation <- data.frame(
    OlinkID = rownames(fit_adjusted$t),
    t = as.numeric(fit_adjusted$t[, "COVIDCOVID+"]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(
      olink_mapping %>%
        dplyr::select(
          OlinkID,
          Assay,
          UniProt
        ),
      by = "OlinkID"
    )
  
  # Add Entrez IDs from the saved robust GSEA ranking
  robust_rank_annotation <- robust_rank_annotation %>%
    dplyr::left_join(
      readRDS(
        file.path(
          proteomics_results_path,
          "GSEA",
          "Proteomics_TRUE_Robust_GSEA_Ranking_Adjusted_t.rds"
        )
      ) %>%
        dplyr::select(Entrez),
      by = "OlinkID"
    )
  
  # Annotate positive Reactome core
  positive_core_annotated <- positive_core %>%
    dplyr::left_join(
      robust_rank_annotation,
      by = "Entrez"
    ) %>%
    dplyr::arrange(
      dplyr::desc(pathway_count),
      dplyr::desc(t)
    )
  
  # Annotate negative Reactome core
  negative_core_annotated <- negative_core %>%
    dplyr::left_join(
      robust_rank_annotation,
      by = "Entrez"
    ) %>%
    dplyr::arrange(
      dplyr::desc(pathway_count),
      t
    )
  
  cat("\nPositive Reactome core:\n")
  print(head(positive_core_annotated, 20))
  
  cat("\nNegative Reactome core:\n")
  print(head(negative_core_annotated, 20))
  
  cat(
    "\nPositive core unmapped:",
    sum(is.na(positive_core_annotated$Assay)),
    "\n"
  )
  
  cat(
    "Negative core unmapped:",
    sum(is.na(negative_core_annotated$Assay)),
    "\n"
  )
  
  robust_rank_file <- file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_GSEA_Ranking_Adjusted_t.rds"
  )
  
  robust_rank_check <- readRDS(robust_rank_file)
  
  cat("Class:\n")
  print(class(robust_rank_check))
  
  cat("\nStructure:\n")
  str(robust_rank_check)
  
  cat("\nNames:\n")
  print(names(robust_rank_check))
  
  cat("\nDimensions:\n")
  print(dim(robust_rank_check))
  
  
  library(dplyr)
  
  # Read the saved TRUE robust ranking annotation
  robust_annotation_file <- file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_GSEA_Ranking_Adjusted_t_Annotation.csv"
  )
  
  robust_annotation <- read.csv(
    robust_annotation_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat("Annotation dimensions:\n")
  print(dim(robust_annotation))
  
  cat("\nAnnotation columns:\n")
   cat("\nFirst rows:\n")]
print(names(robust_annotation))
library(dplyr)

# Standardize Entrez ID type
robust_annotation <- robust_annotation %>%
  dplyr::mutate(
    Entrez = as.character(ENTREZID)
  )

# Annotate positive Reactome core
positive_core_annotated <- positive_core %>%
  dplyr::mutate(
    Entrez = as.character(Entrez)
  ) %>%
  dplyr::left_join(
    robust_annotation %>%
      dplyr::select(
        Entrez,
        OlinkID,
        Assay,
        UniProt,
        t,
        Panel
      ),
    by = "Entrez"
  ) %>%
  dplyr::arrange(
    dplyr::desc(pathway_count),
    dplyr::desc(t)
  )

# Annotate negative Reactome core
negative_core_annotated <- negative_core %>%
  dplyr::mutate(
    Entrez = as.character(Entrez)
  ) %>%
  dplyr::left_join(
    robust_annotation %>%
      dplyr::select(
        Entrez,
        OlinkID,
        Assay,
        UniProt,
        t,
        Panel
      ),
    by = "Entrez"
  ) %>%
  dplyr::arrange(
    dplyr::desc(pathway_count),
    t
  )

cat("\nPositive Reactome core annotation:\n")
print(head(positive_core_annotated, 20))

cat("\nNegative Reactome core annotation:\n")
print(head(negative_core_annotated, 20))

cat(
  "\nPositive core unmapped:",
  sum(is.na(positive_core_annotated$Assay)),
  "\n"
)

cat(
  "Negative core unmapped:",
  sum(is.na(negative_core_annotated$Assay)),
  "\n"
)
# Save annotated Reactome positive and negative core signatures

positive_core_annotated <- positive_core_annotated %>%
  dplyr::arrange(
    dplyr::desc(pathway_count),
    dplyr::desc(t)
  )

negative_core_annotated <- negative_core_annotated %>%
  dplyr::arrange(
    dplyr::desc(pathway_count),
    t
  )

saveRDS(
  positive_core_annotated,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_GSEA_Reactome_Positive_Core_Annotated.rds"
  )
)

write.csv(
  positive_core_annotated,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_GSEA_Reactome_Positive_Core_Annotated.csv"
  ),
  row.names = FALSE
)

saveRDS(
  negative_core_annotated,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_GSEA_Reactome_Negative_Core_Annotated.rds"
  )
)

write.csv(
  negative_core_annotated,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_GSEA_Reactome_Negative_Core_Annotated.csv"
  ),
  row.names = FALSE
)

cat("Positive annotated core saved:", nrow(positive_core_annotated), "\n")
cat("Negative annotated core saved:", nrow(negative_core_annotated), "\n")

cat(
  "\nUnique positive Olink assays:",
  dplyr::n_distinct(positive_core_annotated$OlinkID),
  "\n"
)

cat(
  "Unique negative Olink assays:",
  dplyr::n_distinct(negative_core_annotated$OlinkID),
  "\n"
)

cat("\nReactome TRUE Robust leading-edge annotation saved successfully.\n")
# Load TRUE robust GO-BP core results

go_positive_core <- readRDS(
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_GSEA_GO_BP_Positive_Core.rds"
  )
)

# Inspect structure
cat("GO-BP positive core:\n")
print(dim(go_positive_core))

cat("\nColumns:\n")
print(names(go_positive_core))

cat("\nTop rows:\n")
print(head(go_positive_core, 10))


# Collapse GO-BP positive core to unique Olink assays

go_positive_assay <- go_positive_core %>%
  arrange(desc(pathway_count), desc(t)) %>%
  group_by(OlinkID) %>%
  summarise(
    GO_BP_pathway_count = max(pathway_count, na.rm = TRUE),
    ENTREZID = paste(sort(unique(ENTREZID)), collapse = ";"),
    Assay = Assay[1],
    UniProt = UniProt[1],
    t = t[1],
    Panel = Panel[1],
    .groups = "drop"
  )

# Prepare Reactome positive core

reactome_positive_assay <- positive_core_annotated %>%
  arrange(desc(pathway_count), desc(t)) %>%
  group_by(OlinkID) %>%
  summarise(
    Reactome_pathway_count = max(pathway_count, na.rm = TRUE),
    Assay = Assay[1],
    UniProt = UniProt[1],
    t = t[1],
    Panel = Panel[1],
    .groups = "drop"
  )

cat("GO-BP unique positive assays:",
    nrow(go_positive_assay), "\n")

cat("Reactome unique positive assays:",
    nrow(reactome_positive_assay), "\n")

cat("\nGO-BP top 10:\n")
print(head(go_positive_assay, 10))

cat("\nReactome top 10:\n")
print(head(reactome_positive_assay, 10))

# Integrate TRUE Robust positive cores from GO-BP and Reactome

positive_integrated_core <- go_positive_assay %>%
  inner_join(
    reactome_positive_assay,
    by = "OlinkID",
    suffix = c("_GO_BP", "_Reactome")
  ) %>%
  mutate(
    mean_t = rowMeans(
      cbind(t_GO_BP, t_Reactome),
      na.rm = TRUE
    )
  ) %>%
  arrange(
    desc(GO_BP_pathway_count),
    desc(Reactome_pathway_count),
    desc(mean_t)
  )

cat(
  "Integrated positive core assays:",
  nrow(positive_integrated_core),
  "\n"
)

cat(
  "GO-BP positive assays:",
  nrow(go_positive_assay),
  "\n"
)

cat(
  "Reactome positive assays:",
  nrow(reactome_positive_assay),
  "\n"
)

cat("\nTop integrated positive proteins:\n")
print(
  positive_integrated_core %>%
    select(
      OlinkID,
      Assay_GO_BP,
      UniProt_GO_BP,
      GO_BP_pathway_count,
      Reactome_pathway_count,
      t_GO_BP,
      Panel_GO_BP
    ) %>%
    head(20)
)

# Inspect the integrated positive core

print(
  positive_integrated_core %>%
    dplyr::select(
      OlinkID,
      Assay_GO_BP,
      UniProt_GO_BP,
      GO_BP_pathway_count,
      Reactome_pathway_count,
      t_GO_BP,
      Panel_GO_BP
    ) %>%
    head(20)
)


# Load TRUE robust GO-BP negative core

go_negative_core <- readRDS(
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_GSEA_GO_BP_Negative_Core.rds"
  )
)

cat("GO-BP negative core rows:",
    nrow(go_negative_core), "\n")

cat("\nColumns:\n")
print(names(go_negative_core))

cat("\nTop rows:\n")
print(head(go_negative_core, 10))


# Collapse GO-BP negative core to unique Olink assays

go_negative_assay <- go_negative_core %>%
  arrange(desc(pathway_count), t) %>%
  group_by(OlinkID) %>%
  summarise(
    GO_BP_pathway_count = max(pathway_count, na.rm = TRUE),
    ENTREZID = paste(sort(unique(ENTREZID)), collapse = ";"),
    Assay = Assay[1],
    UniProt = UniProt[1],
    t = t[1],
    Panel = Panel[1],
    .groups = "drop"
  )

# Prepare Reactome negative core

reactome_negative_assay <- negative_core_annotated %>%
  arrange(desc(pathway_count), t) %>%
  group_by(OlinkID) %>%
  summarise(
    Reactome_pathway_count = max(pathway_count, na.rm = TRUE),
    Assay = Assay[1],
    UniProt = UniProt[1],
    t = t[1],
    Panel = Panel[1],
    .groups = "drop"
  )

cat(
  "GO-BP unique negative assays:",
  nrow(go_negative_assay),
  "\n"
)

cat(
  "Reactome unique negative assays:",
  nrow(reactome_negative_assay),
  "\n"
)

cat("\nGO-BP top 10:\n")
print(head(go_negative_assay, 10))

cat("\nReactome top 10:\n")
print(head(reactome_negative_assay, 10))


# Integrate TRUE Robust negative cores from GO-BP and Reactome

negative_integrated_core <- go_negative_assay %>%
  inner_join(
    reactome_negative_assay,
    by = "OlinkID",
    suffix = c("_GO_BP", "_Reactome")
  ) %>%
  mutate(
    mean_t = rowMeans(
      cbind(t_GO_BP, t_Reactome),
      na.rm = TRUE
    )
  ) %>%
  arrange(
    desc(GO_BP_pathway_count),
    desc(Reactome_pathway_count),
    t_GO_BP
  )

cat(
  "Integrated negative core assays:",
  nrow(negative_integrated_core),
  "\n"
)

cat(
  "GO-BP negative assays:",
  nrow(go_negative_assay),
  "\n"
)

cat(
  "Reactome negative assays:",
  nrow(reactome_negative_assay),
  "\n"
)

cat("\nTop integrated negative proteins:\n")

print(
  negative_integrated_core %>%
    dplyr::select(
      OlinkID,
      Assay_GO_BP,
      UniProt_GO_BP,
      GO_BP_pathway_count,
      Reactome_pathway_count,
      t_GO_BP,
      Panel_GO_BP
    ) %>%
    head(20)
)

# Build and save the final TRUE Robust integrated proteomics core

integrated_proteomics_core <- bind_rows(
  positive_integrated_core %>%
    mutate(Direction = "Positive"),
  negative_integrated_core %>%
    mutate(Direction = "Negative")
) %>%
  mutate(
    OlinkID = OlinkID,
    Assay = Assay_GO_BP,
    UniProt = UniProt_GO_BP,
    Panel = Panel_GO_BP,
    t_adjusted = t_GO_BP
  ) %>%
  dplyr::select(
    Direction,
    OlinkID,
    Assay,
    UniProt,
    Panel,
    t_adjusted,
    GO_BP_pathway_count,
    Reactome_pathway_count
  ) %>%
  arrange(
    factor(Direction, levels = c("Positive", "Negative")),
    desc(GO_BP_pathway_count),
    desc(Reactome_pathway_count),
    desc(t_adjusted)
  )

# Save integrated core
saveRDS(
  integrated_proteomics_core,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_Integrated_Core_28_Assays.rds"
  )
)

write.csv(
  integrated_proteomics_core,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_Integrated_Core_28_Assays.csv"
  ),
  row.names = FALSE
)

# Save direction-specific tables
integrated_positive_core <- integrated_proteomics_core %>%
  filter(Direction == "Positive")

integrated_negative_core <- integrated_proteomics_core %>%
  filter(Direction == "Negative")

saveRDS(
  integrated_positive_core,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_Integrated_Positive_Core_19_Assays.rds"
  )
)

write.csv(
  integrated_positive_core,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_Integrated_Positive_Core_19_Assays.csv"
  ),
  row.names = FALSE
)

saveRDS(
  integrated_negative_core,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_Integrated_Negative_Core_9_Assays.rds"
  )
)

write.csv(
  integrated_negative_core,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_Integrated_Negative_Core_9_Assays.csv"
  ),
  row.names = FALSE
)

# Validation
cat("Total integrated assays:", nrow(integrated_proteomics_core), "\n")
cat("Positive integrated assays:", nrow(integrated_positive_core), "\n")
cat("Negative integrated assays:", nrow(integrated_negative_core), "\n")

cat("\nDuplicate OlinkIDs:",
    sum(duplicated(integrated_proteomics_core$OlinkID)),
    "\n")

cat("\nDirection counts:\n")
print(table(integrated_proteomics_core$Direction))

cat("\nIntegrated core:\n")
print(integrated_proteomics_core)

# Characterize integrated proteomics core by Olink panel

integrated_panel_summary <- integrated_proteomics_core %>%
  count(Direction, Panel, name = "Protein_Count") %>%
  arrange(
    factor(Direction, levels = c("Positive", "Negative")),
    desc(Protein_Count),
    Panel
  )

cat("Integrated core panel distribution:\n")
print(integrated_panel_summary)

cat("\nTotal proteins by panel:\n")
print(
  integrated_proteomics_core %>%
    count(Panel, name = "Protein_Count") %>%
    arrange(desc(Protein_Count))
)

# Save results

saveRDS(
  integrated_panel_summary,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_Integrated_Core_Panel_Summary.rds"
  )
)

write.csv(
  integrated_panel_summary,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_Integrated_Core_Panel_Summary.csv"
  ),
  row.names = FALSE
)

cat("\nPanel summary saved successfully.\n")

# Map the 28-protein integrated core to STRING

library(STRINGdb)

string_db_integrated <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400
)

integrated_string_input <- integrated_proteomics_core %>%
  dplyr::select(
    OlinkID,
    Assay,
    UniProt,
    Direction,
    Panel,
    t_adjusted
  ) %>%
  distinct(UniProt, .keep_all = TRUE)

cat("Integrated proteins before STRING mapping:",
    nrow(integrated_string_input), "\n")

# Map UniProt IDs to STRING

integrated_string_mapped <- string_db_integrated$map(
  integrated_string_input,
  "UniProt",
  removeUnmappedRows = FALSE
)

cat(
  "Rows after STRING mapping:",
  nrow(integrated_string_mapped),
  "\n"
)

cat(
  "Mapped STRING proteins:",
  sum(!is.na(integrated_string_mapped$STRING_id)),
  "\n"
)

cat(
  "Unmapped proteins:",
  sum(is.na(integrated_string_mapped$STRING_id)),
  "\n"
)

cat("\nUnmapped proteins:\n")
print(
  integrated_string_mapped %>%
    filter(is.na(STRING_id)) %>%
    dplyr::select(
      OlinkID,
      Assay,
      UniProt,
      Direction,
      Panel
    )
)

cat("\nMapped proteins:\n")
print(
  integrated_string_mapped %>%
    filter(!is.na(STRING_id)) %>%
    dplyr::select(
      OlinkID,
      Assay,
      UniProt,
      Direction,
      Panel,
      t_adjusted,
      STRING_id
    )
)

# Validate UniProt IDs before direct STRING interaction retrieval

integrated_uniprot <- integrated_proteomics_core %>%
  dplyr::select(
    OlinkID,
    Assay,
    UniProt,
    Direction,
    Panel,
    t_adjusted
  ) %>%
  distinct(UniProt, .keep_all = TRUE)

cat(
  "Unique UniProt IDs:",
  dplyr::n_distinct(integrated_uniprot$UniProt),
  "\n"
)

cat(
  "Total integrated proteins:",
  nrow(integrated_uniprot),
  "\n"
)

cat(
  "Missing UniProt IDs:",
  sum(
    is.na(integrated_uniprot$UniProt) |
      integrated_uniprot$UniProt == ""
  ),
  "\n"
)

cat("\nUniProt IDs:\n")
print(integrated_uniprot$UniProt)

cat("\nDuplicate UniProt IDs:\n")
print(
  integrated_uniprot %>%
    dplyr::count(UniProt) %>%
    dplyr::filter(n > 1)
)


# Test direct STRING mapping using UniProt identifiers

cat("STRINGdb version:\n")
print(packageVersion("STRINGdb"))

cat("\nSTRING database version:\n")
print(string_db_integrated$get_version())

cat("\nSTRING species:\n")
print(string_db_integrated$get_species())

cat("\nTesting one UniProt ID:\n")

test_mapping <- tryCatch(
  {
    string_db_integrated$get_proteins(
      string_db_integrated$get_proteins()
    )
  },
  error = function(e) {
    e
  }
)

print(test_mapping)

# Inspect available STRINGdb methods

cat("Available STRINGdb methods:\n")

print(
  methods(class = "STRINGdb")
)

cat("\nSTRINGdb object fields/methods:\n")

print(
  ls(string_db_integrated)
)
# Inspect STRING protein table and identify UniProt mapping column

cat("STRING protein table dimensions:\n")
print(dim(string_db_integrated$proteins))

cat("\nSTRING protein table columns:\n")
print(names(string_db_integrated$proteins))

cat("\nSTRING database version field:\n")
print(string_db_integrated$version)

cat("\nFirst rows of STRING protein table:\n")
print(head(string_db_integrated$proteins))

# Test STRING API connection with the 28 integrated proteins

library(httr)

integrated_uniprot <- c(
  "P01579","Q10589","P00813","O95786","P21589","Q9H773",
  "Q9UKK9","P29466","P19474","P32456","P0DMV8","P19971",
  "O75356","P38484","Q06830","P09958","O95544","Q7LG56",
  "P30041","Q9UMF0","P12644","P05556","Q13332","P49747",
  "P12830","P08648","P22105","P31997"
)

string_ids <- paste0("9606.", integrated_uniprot, collapse = "%0d")

url <- paste0(
  "https://string-db.org/api/tsv/network?identifiers=",
  string_ids,
  "&species=9606",
  "&required_score=400"
)

response <- GET(url)

cat("HTTP status:", status_code(response), "\n")

if (status_code(response) == 200) {
  string_network <- read.delim(
    text = content(response, "text", encoding = "UTF-8"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  cat("STRING network downloaded successfully.\n")
  cat("Rows:", nrow(string_network), "\n")
  cat("Columns:", ncol(string_network), "\n")
  
  print(head(string_network))
} else {
  cat("STRING API request failed.\n")
  print(content(response, "text"))
}

# Test STRING using gene symbols

integrated_genes <- c(
  "IFNG",
  "BST2",
  "ADA",
  "DDX58",
  "NT5E",
  "DCTPP1",
  "NUDT5",
  "CASP1",
  "TRIM21",
  "GBP2",
  "HSPA1A",
  "TYMP",
  "ENTPD5",
  "IFNGR2",
  "PRDX1",
  "FURIN",
  "NADK",
  "RRM2B",
  "PRDX6",
  "ICAM5",
  "BMP4",
  "ITGB1",
  "PTPRS",
  "COMP",
  "CDH1",
  "ITGA5",
  "TNXB",
  "CEACAM8"
)

response <- httr::POST(
  "https://string-db.org/api/tsv/network",
  body = list(
    identifiers = paste(integrated_genes, collapse = "\r"),
    species = 9606,
    required_score = 400
  ),
  encode = "form"
)

cat("HTTP status:", httr::status_code(response), "\n")

if (httr::status_code(response) == 200) {
  
  string_network <- read.delim(
    text = httr::content(response, "text", encoding = "UTF-8"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  cat("STRING network downloaded successfully.\n")
  cat("Edges:", nrow(string_network), "\n")
  cat("Columns:", ncol(string_network), "\n\n")
  
  print(head(string_network))
  
} else {
  
  cat("STRING API request failed.\n")
  print(httr::content(response, "text", encoding = "UTF-8"))
}

# Keep only interactions among the 28 integrated proteins
# and prepare the PPI edge table

integrated_ppi_edges <- string_network[
  string_network$preferredName_A %in% integrated_genes &
    string_network$preferredName_B %in% integrated_genes,
]

integrated_ppi_edges <- integrated_ppi_edges[
  integrated_ppi_edges$preferredName_A != integrated_ppi_edges$preferredName_B,
]

integrated_ppi_edges <- integrated_ppi_edges[
  order(-integrated_ppi_edges$score),
]

cat("Integrated PPI edges:", nrow(integrated_ppi_edges), "\n")
cat("Unique proteins connected:",
    length(unique(c(
      integrated_ppi_edges$preferredName_A,
      integrated_ppi_edges$preferredName_B
    ))), "\n\n")

cat("Top interactions:\n")
print(
  integrated_ppi_edges[
    c("preferredName_A", "preferredName_B", "score")
  ]
)

# Calculate PPI network centrality measures

library(igraph)

ppi_graph <- graph_from_data_frame(
  integrated_ppi_edges[, c("preferredName_A", "preferredName_B", "score")],
  directed = FALSE
)

V(ppi_graph)$name <- as.character(V(ppi_graph)$name)

centrality_table <- data.frame(
  Protein = V(ppi_graph)$name,
  Degree = degree(ppi_graph),
  Betweenness = betweenness(ppi_graph, normalized = TRUE),
  Closeness = closeness(ppi_graph, normalized = TRUE),
  stringsAsFactors = FALSE
)

centrality_table <- centrality_table[
  order(
    -centrality_table$Degree,
    -centrality_table$Betweenness,
    -centrality_table$Closeness
  ),
]

cat("Number of nodes:", vcount(ppi_graph), "\n")
cat("Number of edges:", ecount(ppi_graph), "\n\n")

cat("PPI centrality ranking:\n")
print(centrality_table)

write.csv(
  centrality_table,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_Integrated_Core_PPI_Centrality.csv"
  ),
  row.names = FALSE
)

saveRDS(
  centrality_table,
  file.path(
    proteomics_results_path,
    "GSEA",
    "Proteomics_TRUE_Robust_Integrated_Core_PPI_Centrality.rds"
  )
)

cat("\nCentrality results saved successfully.\n")

# Build and save the final integrated proteomics PPI network

library(igraph)

# Prepare node annotation
ppi_nodes <- data.frame(
  Protein = V(ppi_graph)$name,
  stringsAsFactors = FALSE
)

ppi_nodes <- ppi_nodes %>%
  left_join(
    integrated_proteomics_core %>%
      dplyr::select(
        Assay,
        OlinkID,
        Direction,
        Panel,
        t_adjusted,
        GO_BP_pathway_count,
        Reactome_pathway_count
      ) %>%
      dplyr::rename(Protein = Assay),
    by = "Protein"
  ) %>%
  left_join(
    centrality_table,
    by = "Protein"
  )

# Add node attributes
V(ppi_graph)$degree <- ppi_nodes$Degree
V(ppi_graph)$betweenness <- ppi_nodes$Betweenness
V(ppi_graph)$closeness <- ppi_nodes$Closeness
V(ppi_graph)$Direction <- ppi_nodes$Direction

# Edge width based on STRING confidence score
E(ppi_graph)$width <- 1 + 5 * E(ppi_graph)$score

# Node size based on degree
V(ppi_graph)$size <- 8 + 3 * V(ppi_graph)$degree

# Output paths
gsea_path <- file.path(proteomics_results_path, "GSEA")

# Save node table
write.csv(
  ppi_nodes,
  file.path(
    gsea_path,
    "Proteomics_TRUE_Robust_Integrated_Core_PPI_Nodes.csv"
  ),
  row.names = FALSE
)

saveRDS(
  ppi_nodes,
  file.path(
    gsea_path,
    "Proteomics_TRUE_Robust_Integrated_Core_PPI_Nodes.rds"
  )
)

# Save edge table
write.csv(
  integrated_ppi_edges,
  file.path(
    gsea_path,
    "Proteomics_TRUE_Robust_Integrated_Core_PPI_Edges.csv"
  ),
  row.names = FALSE
)

saveRDS(
  integrated_ppi_edges,
  file.path(
    gsea_path,
    "Proteomics_TRUE_Robust_Integrated_Core_PPI_Edges.rds"
  )
)

# Save igraph object
saveRDS(
  ppi_graph,
  file.path(
    gsea_path,
    "Proteomics_TRUE_Robust_Integrated_Core_PPI_igraph.rds"
  )
)

# Draw publication-quality network
png(
  filename = file.path(
    gsea_path,
    "Proteomics_TRUE_Robust_Integrated_Core_PPI_Network.png"
  ),
  width = 2400,
  height = 2000,
  res = 300
)

set.seed(123)

plot(
  ppi_graph,
  layout = layout_with_fr(ppi_graph),
  vertex.label = V(ppi_graph)$name,
  vertex.label.cex = 0.75,
  vertex.label.color = "black",
  vertex.size = V(ppi_graph)$size,
  vertex.frame.color = "black",
  edge.width = E(ppi_graph)$width,
  edge.color = "grey60",
  main = "Integrated Proteomics PPI Network"
)

legend(
  "topleft",
  legend = c(
    "Positive integrated core",
    "Negative integrated core"
  ),
  pch = 21,
  pt.cex = 2,
  col = "black",
  pt.bg = c("white", "grey70"),
  bty = "n"
)

dev.off()

cat("Final PPI network completed successfully.\n")
cat("Nodes:", vcount(ppi_graph), "\n")
cat("Edges:", ecount(ppi_graph), "\n")
cat("Node table rows:", nrow(ppi_nodes), "\n")
cat("Edge table rows:", nrow(integrated_ppi_edges), "\n")
cat("Network PNG saved.\n")


# Restore the primary transcriptomics dataset for longitudinal analysis

library(DESeq2)
library(dplyr)

dds_longitudinal <- readRDS(
  "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/dds_primary_filtered.rds"
)

metadata_longitudinal <- as.data.frame(colData(dds_longitudinal))

cat("Transcriptomics object restored successfully.\n")
cat("Samples:", ncol(dds_longitudinal), "\n")
cat("Genes:", nrow(dds_longitudinal), "\n\n")

cat("Metadata columns:\n")
print(colnames(metadata_longitudinal))

cat("\nPatient ID column check:\n")
print("patient_id" %in% colnames(metadata_longitudinal))

cat("\nTimepoint-related columns:\n")
print(
  colnames(metadata_longitudinal)[
    grepl("time|day|visit|patient", colnames(metadata_longitudinal),
          ignore.case = TRUE)
  ]
)

if ("patient_id" %in% colnames(metadata_longitudinal)) {
  cat("\nUnique patients:",
      length(unique(metadata_longitudinal$patient_id)), "\n")
}

cat("\nObject design:\n")
print(design(dds_longitudinal))

# Inspect sample identifiers to recover patient-level pairing information

cat("First 30 sample IDs:\n")
print(head(metadata_longitudinal$sample_id, 30))

cat("\nFirst 30 GEO accessions:\n")
print(head(metadata_longitudinal$geo_accession, 30))

cat("\nTimepoint distribution:\n")
print(table(metadata_longitudinal$time_point, useNA = "ifany"))

cat("\nUnique sample IDs:\n")
cat(length(unique(metadata_longitudinal$sample_id)), "\n")

cat("\nUnique GEO accessions:\n")
cat(length(unique(metadata_longitudinal$geo_accession)), "\n")


# Retrieve original GEO sample metadata for patient-level longitudinal pairing

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}

library(GEOquery)

gse_longitudinal <- getGEO(
  "GSE212041",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

cat("GEO object retrieved successfully.\n")
cat("Number of GEO platforms:", length(gse_longitudinal), "\n\n")

geo_eset <- gse_longitudinal[[1]]

geo_pheno <- pData(geo_eset)

cat("GEO samples:", nrow(geo_pheno), "\n")
cat("GEO metadata columns:\n")
print(colnames(geo_pheno))

cat("\nColumns containing patient/subject/sample information:\n")
print(
  colnames(geo_pheno)[
    grepl(
      "patient|subject|participant|donor|sample",
      colnames(geo_pheno),
      ignore.case = TRUE
    )
  ]
)

# Check both GEO platforms and identify the full sample set

cat("Number of samples in each GEO platform:\n")

for (i in seq_along(gse_longitudinal)) {
  cat(
    "Platform", i,
    ":", nrow(pData(gse_longitudinal[[i]])),
    "samples |",
    "GPL:", annotation(gse_longitudinal[[i]]),
    "\n"
  )
}

cat("\nPlatform 2 metadata columns:\n")
print(colnames(pData(gse_longitudinal[[2]])))

cat("\nPlatform 2 sample count:\n")
print(nrow(pData(gse_longitudinal[[2]])))

cat("\nPlatform 2 relevant characteristics:\n")
print(
  colnames(pData(gse_longitudinal[[2]]))[
    grepl(
      "patient|subject|participant|donor|sample|time|acuity",
      colnames(pData(gse_longitudinal[[2]])),
      ignore.case = TRUE
    )
  ]
)
# Inspect the actual metadata text for patient/sample identifiers

geo_pheno2 <- pData(gse_longitudinal[[2]])

cat("Number of samples:", nrow(geo_pheno2), "\n\n")

cat("First 20 sample titles:\n")
print(
  geo_pheno2[, c(
    "geo_accession",
    "title",
    "source_name_ch1",
    "characteristics_ch1",
    "characteristics_ch1.1",
    "characteristics_ch1.2",
    "characteristics_ch1.3",
    "characteristics_ch1.4"
  )],
  row.names = FALSE
)

cat("\n\nSearching all metadata text for patient/subject/participant identifiers:\n")

metadata_text_cols <- c(
  "title",
  "source_name_ch1",
  "characteristics_ch1",
  "characteristics_ch1.1",
  "characteristics_ch1.2",
  "characteristics_ch1.3",
  "characteristics_ch1.4",
  "description",
  "description.1"
)

for (col in metadata_text_cols) {
  if (col %in% colnames(geo_pheno2)) {
    cat("\n---", col, "---\n")
    vals <- unique(as.character(geo_pheno2[[col]]))
    vals <- vals[
      grepl(
        "patient|subject|participant|donor|individual|ID|id|sample",
        vals,
        ignore.case = TRUE
      )
    ]
    print(head(vals, 30))
  }
}

# Match GEO platform 2 samples to the RNA-seq dataset

geo_pheno2 <- pData(gse_longitudinal[[2]])

dds_sample_ids <- colnames(dds_longitudinal)

geo_ids <- as.character(geo_pheno2$geo_accession)

cat("RNA-seq samples:", length(dds_sample_ids), "\n")
cat("GEO platform 2 samples:", length(geo_ids), "\n")
cat("GEO IDs unique:", length(unique(geo_ids)), "\n\n")

# Samples present in RNA-seq but missing from GEO platform 2
missing_from_geo <- setdiff(dds_sample_ids, geo_ids)

# Samples present in GEO platform 2 but missing from RNA-seq
missing_from_dds <- setdiff(geo_ids, dds_sample_ids)

cat("RNA-seq samples missing from GEO platform 2:",
    length(missing_from_geo), "\n")
print(missing_from_geo)

cat("\nGEO platform 2 samples missing from RNA-seq:",
    length(missing_from_dds), "\n")
print(missing_from_dds)

# Build patient/timepoint information from GEO titles
geo_patient_info <- geo_pheno2 %>%
  dplyr::transmute(
    geo_accession = as.character(geo_accession),
    title = as.character(title),
    Patient_ID = sub("_.*$", "", title),
    Timepoint_from_title = sub("^[^_]+_", "", title)
  )

cat("\nUnique patients in GEO platform 2:",
    length(unique(geo_patient_info$Patient_ID)), "\n")

cat("\nTimepoints from GEO title:\n")
print(table(geo_patient_info$Timepoint_from_title))

cat("\nFirst 20 parsed samples:\n")
print(head(geo_patient_info, 20), row.names = FALSE)

# Match RNA-seq sample IDs to GEO sample titles

geo_patient_info <- geo_pheno2 %>%
  dplyr::transmute(
    geo_accession = as.character(geo_accession),
    title = as.character(title),
    Patient_ID = sub("_.*$", "", title),
    Timepoint_from_title = sub("^[^_]+_", "", title)
  )

# RNA sample IDs
dds_sample_ids <- colnames(dds_longitudinal)

# Exact matching using GEO title
matched_sample_ids <- intersect(dds_sample_ids, geo_patient_info$title)

# Samples in RNA-seq but absent from GEO titles
missing_from_geo_titles <- setdiff(dds_sample_ids, geo_patient_info$title)

# Samples in GEO titles but absent from RNA-seq
missing_from_dds_titles <- setdiff(geo_patient_info$title, dds_sample_ids)

cat("RNA-seq samples:", length(dds_sample_ids), "\n")
cat("GEO titles:", nrow(geo_patient_info), "\n")
cat("Exact matched samples:", length(matched_sample_ids), "\n\n")

cat("RNA-seq samples missing from GEO titles:",
    length(missing_from_geo_titles), "\n")
print(missing_from_geo_titles)

cat("\nGEO titles missing from RNA-seq:",
    length(missing_from_dds_titles), "\n")
print(missing_from_dds_titles)

cat("\nGEO-only titles by timepoint:\n")
print(
  table(
    geo_patient_info$Timepoint_from_title[
      geo_patient_info$title %in% missing_from_dds_titles
    ]
  )
)

cat("\nRNA-only samples by timepoint:\n")
print(
  table(
    sub("^[^_]+_", "", missing_from_geo_titles)
  )
)

# Build the longitudinal metadata table from GEO and RNA-seq metadata

geo_longitudinal_metadata <- geo_patient_info %>%
  dplyr::filter(title %in% matched_sample_ids) %>%
  dplyr::left_join(
    geo_pheno2 %>%
      dplyr::select(
        geo_accession,
        `acuity.max:ch1`,
        `cell type:ch1`,
        `covid-19 status:ch1`,
        `patient category:ch1`,
        `time point:ch1`
      ),
    by = "geo_accession"
  ) %>%
  dplyr::transmute(
    sample_id = title,
    geo_accession = geo_accession,
    Patient_ID = Patient_ID,
    Timepoint = Timepoint_from_title,
    COVID_status_GEO = `covid-19 status:ch1`,
    Patient_category_GEO = `patient category:ch1`,
    Acuity_max_GEO = `acuity.max:ch1`,
    Cell_type_GEO = `cell type:ch1`,
    Timepoint_GEO = `time point:ch1`
  ) %>%
  dplyr::arrange(
    suppressWarnings(as.numeric(Patient_ID)),
    factor(Timepoint, levels = c("D0", "D3", "D7", "DE"))
  )

cat("Longitudinal metadata rows:", nrow(geo_longitudinal_metadata), "\n")
cat("Unique patients:", length(unique(geo_longitudinal_metadata$Patient_ID)), "\n")
cat("Unique GEO accessions:", length(unique(geo_longitudinal_metadata$geo_accession)), "\n\n")

cat("Timepoint distribution:\n")
print(table(geo_longitudinal_metadata$Timepoint))

cat("\nCOVID status distribution:\n")
print(table(geo_longitudinal_metadata$COVID_status_GEO))

cat("\nPatient category distribution:\n")
print(table(geo_longitudinal_metadata$Patient_category_GEO))

cat("\nDuplicated sample IDs:\n")
print(
  geo_longitudinal_metadata$sample_id[
    duplicated(geo_longitudinal_metadata$sample_id)
  ]
)

cat("\nDuplicated patient × timepoint combinations:\n")
duplicate_pt <- geo_longitudinal_metadata %>%
  dplyr::count(Patient_ID, Timepoint, name = "n") %>%
  dplyr::filter(n > 1)

print(duplicate_pt)

cat("\nFirst 20 rows:\n")
print(head(geo_longitudinal_metadata, 20), row.names = FALSE)

# Create the final longitudinal RNA-seq metadata

metadata_longitudinal_final <- metadata_longitudinal %>%
  dplyr::left_join(
    geo_longitudinal_metadata %>%
      dplyr::select(
        sample_id,
        geo_accession,
        Patient_ID,
        Timepoint,
        COVID_status_GEO,
        Patient_category_GEO,
        Acuity_max_GEO,
        Cell_type_GEO,
        Timepoint_GEO
      ),
    by = "sample_id"
  )

cat("Final metadata rows:", nrow(metadata_longitudinal_final), "\n")
cat("Final metadata columns:", ncol(metadata_longitudinal_final), "\n\n")

cat("Samples with Patient_ID:", sum(!is.na(metadata_longitudinal_final$Patient_ID)), "\n")
cat("Samples missing Patient_ID:", sum(is.na(metadata_longitudinal_final$Patient_ID)), "\n\n")

cat("Timepoint distribution:\n")
print(table(metadata_longitudinal_final$Timepoint, useNA = "ifany"))

cat("\nCOVID status distribution:\n")
print(table(metadata_longitudinal_final$COVID_status_GEO, useNA = "ifany"))

cat("\nPatient category distribution:\n")
print(table(metadata_longitudinal_final$Patient_category_GEO, useNA = "ifany"))

cat("\nAcquity max distribution:\n")
print(table(metadata_longitudinal_final$Acuity_max_GEO, useNA = "ifany"))

cat("\nPatient × Timepoint duplicates:\n")
duplicate_final <- metadata_longitudinal_final %>%
  dplyr::filter(!is.na(Patient_ID), !is.na(Timepoint)) %>%
  dplyr::count(Patient_ID, Timepoint, name = "n") %>%
  dplyr::filter(n > 1)

print(duplicate_final)

cat("\nFirst 15 rows:\n")
print(
  metadata_longitudinal_final %>%
    dplyr::select(
      sample_id,
      geo_accession,
      Patient_ID,
      Timepoint,
      COVID_status_GEO,
      Patient_category_GEO,
      Acuity_max_GEO
    ) %>%
    head(15),
  row.names = FALSE
)
# Save the final longitudinal metadata mapping

longitudinal_metadata_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression"

dir.create(
  longitudinal_metadata_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  metadata_longitudinal_final,
  file.path(
    longitudinal_metadata_dir,
    "Longitudinal_RNAseq_Final_Metadata_773_samples.rds"
  )
)

write.csv(
  metadata_longitudinal_final,
  file.path(
    longitudinal_metadata_dir,
    "Longitudinal_RNAseq_Final_Metadata_773_samples.csv"
  ),
  row.names = FALSE
)

# Save the 13 RNA-seq samples without GEO metadata separately

longitudinal_missing_geo <- metadata_longitudinal_final %>%
  dplyr::filter(is.na(Patient_ID))

saveRDS(
  longitudinal_missing_geo,
  file.path(
    longitudinal_metadata_dir,
    "Longitudinal_RNAseq_13_samples_missing_GEO_metadata.rds"
  )
)

write.csv(
  longitudinal_missing_geo,
  file.path(
    longitudinal_metadata_dir,
    "Longitudinal_RNAseq_13_samples_missing_GEO_metadata.csv"
  ),
  row.names = FALSE
)

cat("Final metadata saved successfully.\n")
cat("Full metadata:", nrow(metadata_longitudinal_final), "samples\n")
cat("Mapped samples:", sum(!is.na(metadata_longitudinal_final$Patient_ID)), "\n")
cat("Missing GEO metadata:", sum(is.na(metadata_longitudinal_final$Patient_ID)), "\n")
cat("Missing-metadata file:", nrow(longitudinal_missing_geo), "samples\n")


# Characterize the longitudinal cohort

longitudinal_cohort <- metadata_longitudinal_final %>%
  dplyr::filter(
    !is.na(Patient_ID),
    Timepoint %in% c("D0", "D3", "D7")
  ) %>%
  dplyr::mutate(
    Patient_ID = as.character(Patient_ID)
  )

patient_timepoint_summary <- longitudinal_cohort %>%
  dplyr::distinct(Patient_ID, Timepoint) %>%
  dplyr::count(Patient_ID, name = "n_timepoints")

cat("Longitudinal samples:", nrow(longitudinal_cohort), "\n")
cat("Unique longitudinal patients:",
    dplyr::n_distinct(longitudinal_cohort$Patient_ID), "\n\n")

cat("Samples by timepoint:\n")
print(table(longitudinal_cohort$Timepoint))

cat("\nPatients by number of timepoints:\n")
print(table(patient_timepoint_summary$n_timepoints))

cat("\nPatients with D0 + D3 + D7:\n")
complete_3tp <- patient_timepoint_summary %>%
  dplyr::filter(n_timepoints == 3)

cat(nrow(complete_3tp), "\n")

cat("\nPatients with at least 2 timepoints:\n")
cat(sum(patient_timepoint_summary$n_timepoints >= 2), "\n")

cat("\nPatients with only 1 timepoint:\n")
cat(sum(patient_timepoint_summary$n_timepoints == 1), "\n")

cat("\nCOVID status among longitudinal samples:\n")
print(table(longitudinal_cohort$COVID_status_GEO))

cat("\nPatient category among longitudinal samples:\n")
print(table(longitudinal_cohort$Patient_category_GEO))

cat("\nAcuity among longitudinal samples:\n")
print(table(longitudinal_cohort$Acuity_max_GEO))

cat("\nTimepoint combination counts:\n")

patient_timepoint_wide <- longitudinal_cohort %>%
  dplyr::distinct(Patient_ID, Timepoint) %>%
  dplyr::mutate(present = 1) %>%
  tidyr::pivot_wider(
    names_from = Timepoint,
    values_from = present,
    values_fill = 0
  )

print(
  patient_timepoint_wide %>%
    dplyr::count(D0, D3, D7, name = "patients") %>%
    dplyr::arrange(desc(patients))
)


# Validate patient-level consistency across longitudinal timepoints

patient_consistency <- longitudinal_cohort %>%
  dplyr::group_by(Patient_ID) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_timepoints = dplyr::n_distinct(Timepoint),
    n_covid_status = dplyr::n_distinct(COVID_status_GEO),
    n_patient_category = dplyr::n_distinct(Patient_category_GEO),
    n_acuity = dplyr::n_distinct(Acuity_max_GEO),
    covid_status = paste(sort(unique(COVID_status_GEO)), collapse = ";"),
    patient_category = paste(sort(unique(Patient_category_GEO)), collapse = ";"),
    acuity_values = paste(sort(unique(Acuity_max_GEO)), collapse = ";"),
    .groups = "drop"
  )

cat("Patients:", nrow(patient_consistency), "\n\n")

cat("COVID status inconsistencies:\n")
print(
  patient_consistency %>%
    dplyr::filter(n_covid_status > 1)
)

cat("\nPatient category inconsistencies:\n")
print(
  patient_consistency %>%
    dplyr::filter(n_patient_category > 1)
)

cat("\nAcuity variation across timepoints:\n")
print(
  table(patient_consistency$n_acuity)
)

cat("\nPatients with COVID status = 0:\n")
print(
  patient_consistency %>%
    dplyr::filter(covid_status == "0") %>%
    dplyr::select(
      Patient_ID,
      n_samples,
      n_timepoints,
      covid_status,
      patient_category,
      acuity_values
    )
)

cat("\nPatients with COVID status inconsistencies:",
    sum(patient_consistency$n_covid_status > 1), "\n")

cat("Patients with patient-category inconsistencies:",
    sum(patient_consistency$n_patient_category > 1), "\n")

# Define the final COVID-positive longitudinal cohort

covid_longitudinal <- metadata_longitudinal_final %>%
  dplyr::filter(
    !is.na(Patient_ID),
    Patient_category_GEO == "COVID+",
    Timepoint %in% c("D0", "D3", "D7")
  ) %>%
  dplyr::mutate(
    Patient_ID = as.character(Patient_ID),
    Timepoint = factor(
      Timepoint,
      levels = c("D0", "D3", "D7")
    ),
    COVID_status = factor(
      COVID_status_GEO,
      levels = c(0, 1),
      labels = c("COVID-", "COVID+")
    )
  )

covid_patient_summary <- covid_longitudinal %>%
  dplyr::distinct(Patient_ID, Timepoint) %>%
  dplyr::count(Patient_ID, name = "n_timepoints")

cat("COVID+ longitudinal samples:", nrow(covid_longitudinal), "\n")
cat(
  "COVID+ longitudinal patients:",
  dplyr::n_distinct(covid_longitudinal$Patient_ID),
  "\n\n"
)

cat("Samples by timepoint:\n")
print(table(covid_longitudinal$Timepoint))

cat("\nPatients by number of timepoints:\n")
print(table(covid_patient_summary$n_timepoints))

cat(
  "\nPatients with >=2 timepoints:",
  sum(covid_patient_summary$n_timepoints >= 2),
  "\n"
)

cat(
  "Patients with all 3 timepoints:",
  sum(covid_patient_summary$n_timepoints == 3),
  "\n\n"
)

cat("Acuity max by timepoint:\n")
print(
  table(
    covid_longitudinal$Timepoint,
    covid_longitudinal$Acuity_max_GEO
  )
)

cat("\nUnique patients by timepoint:\n")
print(
  covid_longitudinal %>%
    dplyr::group_by(Timepoint) %>%
    dplyr::summarise(
      patients = dplyr::n_distinct(Patient_ID),
      samples = dplyr::n(),
      .groups = "drop"
    )
)
# Save the final COVID-positive longitudinal cohort

longitudinal_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"

dir.create(
  longitudinal_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  covid_longitudinal,
  file.path(
    longitudinal_dir,
    "COVID_Positive_Longitudinal_Cohort_635_samples.rds"
  )
)

write.csv(
  covid_longitudinal,
  file.path(
    longitudinal_dir,
    "COVID_Positive_Longitudinal_Cohort_635_samples.csv"
  ),
  row.names = FALSE
)

saveRDS(
  covid_patient_summary,
  file.path(
    longitudinal_dir,
    "COVID_Positive_Longitudinal_Patient_Timepoint_Summary.rds"
  )
)

write.csv(
  covid_patient_summary,
  file.path(
    longitudinal_dir,
    "COVID_Positive_Longitudinal_Patient_Timepoint_Summary.csv"
  ),
  row.names = FALSE
)

cat("COVID+ longitudinal cohort saved successfully.\n")
cat("Patients:", dplyr::n_distinct(covid_longitudinal$Patient_ID), "\n")
cat("Samples:", nrow(covid_longitudinal), "\n")
cat("D0:", sum(covid_longitudinal$Timepoint == "D0"), "\n")
cat("D3:", sum(covid_longitudinal$Timepoint == "D3"), "\n")
cat("D7:", sum(covid_longitudinal$Timepoint == "D7"), "\n")
cat(
  "Patients with >=2 timepoints:",
  sum(covid_patient_summary$n_timepoints >= 2),
  "\n"
)
cat(
  "Patients with all 3 timepoints:",
  sum(covid_patient_summary$n_timepoints == 3),
  "\n"
)

# Diagnose the PCA matrix before scaling

pca_matrix_check <- t(pca_matrix_longitudinal)

# Check dimensions
cat("PCA matrix dimensions:", nrow(pca_matrix_check), "samples x", 
    ncol(pca_matrix_check), "genes\n")

# Count problematic genes
gene_sd <- apply(pca_matrix_check, 2, sd, na.rm = TRUE)

zero_sd_genes <- sum(gene_sd == 0, na.rm = TRUE)
nonfinite_sd_genes <- sum(!is.finite(gene_sd))
na_genes <- sum(!is.finite(as.matrix(pca_matrix_check)) |> colSums() > 0)

cat("Zero-SD genes:", zero_sd_genes, "\n")
cat("Non-finite SD genes:", nonfinite_sd_genes, "\n")
cat("Genes containing non-finite values:", na_genes, "\n")

# Identify all genes that cannot be scaled
bad_genes <- which(
  !is.finite(gene_sd) | gene_sd == 0
)

cat("Total problematic genes:", length(bad_genes), "\n")

# Show first problematic gene indices
if (length(bad_genes) > 0) {
  cat("First problematic gene indices:\n")
  print(head(bad_genes, 20))
}

# Check the actual values of the first problematic gene
if (length(bad_genes) > 0) {
  cat("Values of first problematic gene:\n")
  print(pca_matrix_check[, bad_genes[1]])
}

# Check longitudinal sample IDs and their matching with DESeq2

cat("Number of samples in covid_longitudinal:",
    length(covid_longitudinal$sample_id), "\n")

cat("Number of columns in dds_longitudinal:",
    ncol(dds_longitudinal), "\n")

cat("Number of matching sample IDs:",
    sum(covid_longitudinal$sample_id %in% colnames(dds_longitudinal)), "\n")

cat("First COVID longitudinal sample IDs:\n")
print(head(covid_longitudinal$sample_id, 10))

cat("\nFirst DESeq2 column names:\n")
print(head(colnames(dds_longitudinal), 10))

cat("\nExample matching IDs:\n")
print(
  head(
    covid_longitudinal$sample_id[
      covid_longitudinal$sample_id %in% colnames(dds_longitudinal)
    ],
    10
  )
)

# Inspect the longitudinal metadata used to construct DESeq2 sample IDs

cat("Columns in metadata_longitudinal_final:\n")
print(colnames(metadata_longitudinal_final))

cat("\nFirst rows of Patient_ID and Timepoint_from_title:\n")

if (all(c("Patient_ID", "Timepoint_from_title") %in% colnames(metadata_longitudinal_final))) {
  
  print(
    head(
      metadata_longitudinal_final[
        ,
        c("sample_id", "Patient_ID", "Timepoint_from_title")
      ],
      10
    )
  )
  
} else {
  
  cat(
    "\nPatient_ID or Timepoint_from_title is missing from metadata_longitudinal_final.\n"
  )
}

# Reconstruct the correct DESeq2 sample IDs

covid_longitudinal$dds_sample_id <- paste0(
  covid_longitudinal$Patient_ID,
  "_",
  covid_longitudinal$Timepoint
)

cat("First constructed DESeq2 sample IDs:\n")
print(head(covid_longitudinal$dds_sample_id, 10))

# Check matching with DESeq2
matched_ids <- covid_longitudinal$dds_sample_id[
  covid_longitudinal$dds_sample_id %in% colnames(dds_longitudinal)
]

cat("\nCOVID longitudinal samples:", nrow(covid_longitudinal), "\n")
cat("Matching DESeq2 sample IDs:", length(matched_ids), "\n")
cat(
  "Missing DESeq2 sample IDs:",
  sum(!covid_longitudinal$dds_sample_id %in% colnames(dds_longitudinal)),
  "\n"
)

cat("\nTimepoint distribution of matched samples:\n")
print(table(sub("^[^_]+_", "", matched_ids)))

# Build the normalized expression matrix for the COVID-positive longitudinal cohort

longitudinal_sample_ids <- covid_longitudinal$dds_sample_id

normalized_counts_longitudinal <- counts(
  dds_longitudinal,
  normalized = TRUE
)[
  ,
  longitudinal_sample_ids,
  drop = FALSE
]

log2_normalized_longitudinal <- log2(
  normalized_counts_longitudinal + 1
)

cat("Normalized expression matrix:\n")
cat("Genes:", nrow(log2_normalized_longitudinal), "\n")
cat("Samples:", ncol(log2_normalized_longitudinal), "\n")

# Calculate gene variance across the 635 longitudinal samples
gene_variance <- apply(
  log2_normalized_longitudinal,
  1,
  var,
  na.rm = TRUE
)

valid_genes <- which(
  is.finite(gene_variance) &
    gene_variance > 0
)

top_n <- min(2000, length(valid_genes))

top_variable_gene_indices <- valid_genes[
  order(
    gene_variance[valid_genes],
    decreasing = TRUE
  )[1:top_n]
]

pca_matrix_longitudinal <- log2_normalized_longitudinal[
  top_variable_gene_indices,
  ,
  drop = FALSE
]

cat("\nGenes with non-zero variance:",
    length(valid_genes), "\n")

cat("Genes selected for PCA:",
    nrow(pca_matrix_longitudinal), "\n")

cat("Samples selected for PCA:",
    ncol(pca_matrix_longitudinal), "\n")

# Verify that all selected genes can be scaled
pca_sd <- apply(
  pca_matrix_longitudinal,
  1,
  sd,
  na.rm = TRUE
)

cat("Zero-SD genes:",
    sum(pca_sd == 0, na.rm = TRUE), "\n")

cat("Non-finite SD genes:",
    sum(!is.finite(pca_sd)), "\n")

# Create PCA scores with longitudinal metadata

pca_scores_longitudinal <- as.data.frame(
  pca_longitudinal$x[, 1:5]
)

pca_scores_longitudinal$dds_sample_id <- rownames(
  pca_longitudinal$x
)

# Match metadata using the DESeq2 sample ID
metadata_pca <- covid_longitudinal[
  match(
    pca_scores_longitudinal$dds_sample_id,
    covid_longitudinal$dds_sample_id
  ),
  ,
  drop = FALSE
]

# Add key metadata variables
pca_scores_longitudinal$Patient_ID <- metadata_pca$Patient_ID
pca_scores_longitudinal$Timepoint <- metadata_pca$Timepoint
pca_scores_longitudinal$Acuity_max <- metadata_pca$Acuity_max_GEO

# Validation
cat("PCA samples:", nrow(pca_scores_longitudinal), "\n")
cat("Missing Patient IDs:",
    sum(is.na(pca_scores_longitudinal$Patient_ID)), "\n")
cat("Missing Timepoints:",
    sum(is.na(pca_scores_longitudinal$Timepoint)), "\n")

cat("\nTimepoint distribution:\n")
print(table(pca_scores_longitudinal$Timepoint))

cat("\nPatient timepoint distribution:\n")
print(
  table(
    table(pca_scores_longitudinal$Patient_ID)
  )
)

cat("\nFirst PCA score rows:\n")
print(head(
  pca_scores_longitudinal[
    ,
    c(
      "dds_sample_id",
      "Patient_ID",
      "Timepoint",
      "Acuity_max",
      "PC1",
      "PC2"
    )
  ],
  10
))

# Load plotting package
library(ggplot2)

# Define output directory
pca_output_dir <- file.path(
  "C:/Users/ibrah/OneDrive/Documents/Results",
  "Differential_Expression",
  "Longitudinal"
)

dir.create(
  pca_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Convert timepoint to an ordered factor
pca_scores_longitudinal$Timepoint <- factor(
  pca_scores_longitudinal$Timepoint,
  levels = c("D0", "D3", "D7")
)

# Create PCA plot
pca_plot_longitudinal <- ggplot(
  pca_scores_longitudinal,
  aes(
    x = PC1,
    y = PC2,
    color = Timepoint
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.7
  ) +
  labs(
    title = "Longitudinal PCA of COVID-positive samples",
    subtitle = "D0, D3 and D7 | 635 samples | 2,000 most variable genes",
    x = paste0(
      "PC1 (",
      round(pca_variance[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance[2], 2),
      "%)"
    ),
    color = "Timepoint"
  ) +
  theme_minimal(base_size = 13)

print(pca_plot_longitudinal)

# Save PCA figure
ggsave(
  filename = file.path(
    pca_output_dir,
    "COVID_Positive_Longitudinal_PCA_D0_D3_D7.png"
  ),
  plot = pca_plot_longitudinal,
  width = 9,
  height = 7,
  dpi = 300
)

# Save PCA scores
write.csv(
  pca_scores_longitudinal,
  file = file.path(
    pca_output_dir,
    "COVID_Positive_Longitudinal_PCA_Scores.csv"
  ),
  row.names = FALSE
)

saveRDS(
  pca_scores_longitudinal,
  file = file.path(
    pca_output_dir,
    "COVID_Positive_Longitudinal_PCA_Scores.rds"
  )
)

# Save PCA object and variance explained
saveRDS(
  pca_longitudinal,
  file = file.path(
    pca_output_dir,
    "COVID_Positive_Longitudinal_PCA_Object.rds"
  )
)

write.csv(
  data.frame(
    PC = paste0("PC", seq_along(pca_variance)),
    Variance_Explained_Percent = pca_variance
  ),
  file = file.path(
    pca_output_dir,
    "COVID_Positive_Longitudinal_PCA_Variance_Explained.csv"
  ),
  row.names = FALSE
)

cat("PCA plot saved successfully.\n")
cat("PCA scores saved as CSV and RDS.\n")
cat("PCA object saved as RDS.\n")
cat("Variance explained saved as CSV.\n")

# Identify patients with repeated longitudinal samples

patient_timepoint_counts <- table(
  pca_scores_longitudinal$Patient_ID
)

repeated_patients <- names(
  patient_timepoint_counts[
    patient_timepoint_counts >= 2
  ]
)

pca_trajectory <- pca_scores_longitudinal[
  pca_scores_longitudinal$Patient_ID %in% repeated_patients,
  ,
  drop = FALSE
]

# Order samples within each patient
pca_trajectory$Timepoint <- factor(
  pca_trajectory$Timepoint,
  levels = c("D0", "D3", "D7")
)

pca_trajectory <- pca_trajectory[
  order(
    pca_trajectory$Patient_ID,
    pca_trajectory$Timepoint
  ),
  ,
  drop = FALSE
]

# Validate trajectory cohort
cat("Patients with >=2 timepoints:",
    length(unique(pca_trajectory$Patient_ID)), "\n")

cat("Samples in trajectory analysis:",
    nrow(pca_trajectory), "\n")

cat("\nTimepoint distribution:\n")
print(table(pca_trajectory$Timepoint))

cat("\nNumber of timepoints per patient:\n")
print(
  table(
    table(pca_trajectory$Patient_ID)
  )
)

cat("\nComplete D0-D3-D7 patients:\n")

complete_patients <- names(
  table(pca_trajectory$Patient_ID)[
    table(pca_trajectory$Patient_ID) == 3
  ]
)

cat(length(complete_patients), "\n")


# Prepare complete D0-D3-D7 longitudinal trajectories

complete_trajectory <- pca_trajectory[
  pca_trajectory$Patient_ID %in% complete_patients,
  ,
  drop = FALSE
]

complete_trajectory$Timepoint <- factor(
  complete_trajectory$Timepoint,
  levels = c("D0", "D3", "D7")
)

# Create longitudinal trajectory plot for PC1
trajectory_plot_pc1 <- ggplot(
  complete_trajectory,
  aes(
    x = Timepoint,
    y = PC1,
    group = Patient_ID
  )
) +
  geom_line(
    alpha = 0.20,
    linewidth = 0.5
  ) +
  geom_point(
    aes(color = Timepoint),
    size = 2,
    alpha = 0.70
  ) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    linewidth = 1.3
  ) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "point",
    size = 3.5
  ) +
  labs(
    title = "Within-patient longitudinal trajectory of PC1",
    subtitle = "115 patients with complete D0-D3-D7 sampling",
    x = "Timepoint",
    y = paste0("PC1 (", round(pca_variance[1], 2), "% variance)")
  ) +
  theme_minimal(base_size = 13)

print(trajectory_plot_pc1)

# Create longitudinal trajectory plot for PC2
trajectory_plot_pc2 <- ggplot(
  complete_trajectory,
  aes(
    x = Timepoint,
    y = PC2,
    group = Patient_ID
  )
) +
  geom_line(
    alpha = 0.20,
    linewidth = 0.5
  ) +
  geom_point(
    aes(color = Timepoint),
    size = 2,
    alpha = 0.70
  ) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    linewidth = 1.3
  ) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "point",
    size = 3.5
  ) +
  labs(
    title = "Within-patient longitudinal trajectory of PC2",
    subtitle = "115 patients with complete D0-D3-D7 sampling",
    x = "Timepoint",
    y = paste0("PC2 (", round(pca_variance[2], 2), "% variance)")
  ) +
  theme_minimal(base_size = 13)

print(trajectory_plot_pc2)

# Calculate within-patient changes
trajectory_wide <- reshape(
  complete_trajectory[
    ,
    c("Patient_ID", "Timepoint", "PC1", "PC2")
  ],
  timevar = "Timepoint",
  idvar = "Patient_ID",
  direction = "wide"
)

trajectory_wide$PC1_D3_minus_D0 <- (
  trajectory_wide$PC1.D3 -
    trajectory_wide$PC1.D0
)

trajectory_wide$PC1_D7_minus_D0 <- (
  trajectory_wide$PC1.D7 -
    trajectory_wide$PC1.D0
)

trajectory_wide$PC2_D3_minus_D0 <- (
  trajectory_wide$PC2.D3 -
    trajectory_wide$PC2.D0
)

trajectory_wide$PC2_D7_minus_D0 <- (
  trajectory_wide$PC2.D7 -
    trajectory_wide$PC2.D0
)

# Summary of within-patient changes
trajectory_change_summary <- data.frame(
  Comparison = c(
    "D3 - D0",
    "D7 - D0"
  ),
  PC1_Mean_Change = c(
    mean(trajectory_wide$PC1_D3_minus_D0, na.rm = TRUE),
    mean(trajectory_wide$PC1_D7_minus_D0, na.rm = TRUE)
  ),
  PC1_Median_Change = c(
    median(trajectory_wide$PC1_D3_minus_D0, na.rm = TRUE),
    median(trajectory_wide$PC1_D7_minus_D0, na.rm = TRUE)
  ),
  PC2_Mean_Change = c(
    mean(trajectory_wide$PC2_D3_minus_D0, na.rm = TRUE),
    mean(trajectory_wide$PC2_D7_minus_D0, na.rm = TRUE)
  ),
  PC2_Median_Change = c(
    median(trajectory_wide$PC2_D3_minus_D0, na.rm = TRUE),
    median(trajectory_wide$PC2_D7_minus_D0, na.rm = TRUE)
  )
)

# Save trajectory figures
ggsave(
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_Within_Patient_Trajectory_PC1_D0_D3_D7.png"
  ),
  trajectory_plot_pc1,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_Within_Patient_Trajectory_PC2_D0_D3_D7.png"
  ),
  trajectory_plot_pc2,
  width = 9,
  height = 7,
  dpi = 300
)

# Save trajectory data
write.csv(
  complete_trajectory,
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_Complete_115_Patients_Trajectory_Scores.csv"
  ),
  row.names = FALSE
)

write.csv(
  trajectory_wide,
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_Complete_115_Patients_PC_Changes.csv"
  ),
  row.names = FALSE
)

write.csv(
  trajectory_change_summary,
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_PC_Trajectory_Change_Summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  complete_trajectory,
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_Complete_115_Patients_Trajectory_Scores.rds"
  )
)

saveRDS(
  trajectory_wide,
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_Complete_115_Patients_PC_Changes.rds"
  )
)

saveRDS(
  trajectory_change_summary,
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_PC_Trajectory_Change_Summary.rds"
  )
)

cat("Complete longitudinal trajectory analysis completed.\n\n")

cat("Complete patients:", length(complete_patients), "\n")
cat("Samples analyzed:", nrow(complete_trajectory), "\n\n")

cat("Within-patient PC changes:\n")
print(
  round(
    trajectory_change_summary,
    3
  )
)

cat("\nTrajectory figures and data saved successfully.\n")


# Statistical analysis of within-patient longitudinal PCA changes

# Ensure complete trajectories are ordered correctly
complete_trajectory$Timepoint <- factor(
  complete_trajectory$Timepoint,
  levels = c("D0", "D3", "D7")
)

# Convert PCA scores to wide format
trajectory_wide <- reshape(
  complete_trajectory[
    ,
    c("Patient_ID", "Timepoint", "PC1", "PC2")
  ],
  timevar = "Timepoint",
  idvar = "Patient_ID",
  direction = "wide"
)

# Calculate paired changes
trajectory_wide$PC1_D3_D0 <- trajectory_wide$PC1.D3 -
  trajectory_wide$PC1.D0

trajectory_wide$PC1_D7_D0 <- trajectory_wide$PC1.D7 -
  trajectory_wide$PC1.D0

trajectory_wide$PC1_D7_D3 <- trajectory_wide$PC1.D7 -
  trajectory_wide$PC1.D3

trajectory_wide$PC2_D3_D0 <- trajectory_wide$PC2.D3 -
  trajectory_wide$PC2.D0

trajectory_wide$PC2_D7_D0 <- trajectory_wide$PC2.D7 -
  trajectory_wide$PC2.D0

trajectory_wide$PC2_D7_D3 <- trajectory_wide$PC2.D7 -
  trajectory_wide$PC2.D3

# Paired Wilcoxon tests
pc1_wilcox_D3_D0 <- wilcox.test(
  trajectory_wide$PC1.D3,
  trajectory_wide$PC1.D0,
  paired = TRUE,
  exact = FALSE
)

pc1_wilcox_D7_D0 <- wilcox.test(
  trajectory_wide$PC1.D7,
  trajectory_wide$PC1.D0,
  paired = TRUE,
  exact = FALSE
)

pc1_wilcox_D7_D3 <- wilcox.test(
  trajectory_wide$PC1.D7,
  trajectory_wide$PC1.D3,
  paired = TRUE,
  exact = FALSE
)

pc2_wilcox_D3_D0 <- wilcox.test(
  trajectory_wide$PC2.D3,
  trajectory_wide$PC2.D0,
  paired = TRUE,
  exact = FALSE
)

pc2_wilcox_D7_D0 <- wilcox.test(
  trajectory_wide$PC2.D7,
  trajectory_wide$PC2.D0,
  paired = TRUE,
  exact = FALSE
)

pc2_wilcox_D7_D3 <- wilcox.test(
  trajectory_wide$PC2.D7,
  trajectory_wide$PC2.D3,
  paired = TRUE,
  exact = FALSE
)

# Friedman tests across all three timepoints
pc1_friedman <- friedman.test(
  cbind(
    trajectory_wide$PC1.D0,
    trajectory_wide$PC1.D3,
    trajectory_wide$PC1.D7
  )
)

pc2_friedman <- friedman.test(
  cbind(
    trajectory_wide$PC2.D0,
    trajectory_wide$PC2.D3,
    trajectory_wide$PC2.D7
  )
)

# Create statistical summary
trajectory_statistics <- data.frame(
  PC = c(
    "PC1", "PC1", "PC1",
    "PC2", "PC2", "PC2"
  ),
  Comparison = c(
    "D3 - D0",
    "D7 - D0",
    "D7 - D3",
    "D3 - D0",
    "D7 - D0",
    "D7 - D3"
  ),
  Mean_Change = c(
    mean(trajectory_wide$PC1_D3_D0, na.rm = TRUE),
    mean(trajectory_wide$PC1_D7_D0, na.rm = TRUE),
    mean(trajectory_wide$PC1_D7_D3, na.rm = TRUE),
    mean(trajectory_wide$PC2_D3_D0, na.rm = TRUE),
    mean(trajectory_wide$PC2_D7_D0, na.rm = TRUE),
    mean(trajectory_wide$PC2_D7_D3, na.rm = TRUE)
  ),
  Median_Change = c(
    median(trajectory_wide$PC1_D3_D0, na.rm = TRUE),
    median(trajectory_wide$PC1_D7_D0, na.rm = TRUE),
    median(trajectory_wide$PC1_D7_D3, na.rm = TRUE),
    median(trajectory_wide$PC2_D3_D0, na.rm = TRUE),
    median(trajectory_wide$PC2_D7_D0, na.rm = TRUE),
    median(trajectory_wide$PC2_D7_D3, na.rm = TRUE)
  ),
  SD_Change = c(
    sd(trajectory_wide$PC1_D3_D0, na.rm = TRUE),
    sd(trajectory_wide$PC1_D7_D0, na.rm = TRUE),
    sd(trajectory_wide$PC1_D7_D3, na.rm = TRUE),
    sd(trajectory_wide$PC2_D3_D0, na.rm = TRUE),
    sd(trajectory_wide$PC2_D7_D0, na.rm = TRUE),
    sd(trajectory_wide$PC2_D7_D3, na.rm = TRUE)
  ),
  Wilcoxon_P = c(
    pc1_wilcox_D3_D0$p.value,
    pc1_wilcox_D7_D0$p.value,
    pc1_wilcox_D7_D3$p.value,
    pc2_wilcox_D3_D0$p.value,
    pc2_wilcox_D7_D0$p.value,
    pc2_wilcox_D7_D3$p.value
  )
)

# Add Friedman test results
friedman_results <- data.frame(
  PC = c("PC1", "PC2"),
  Friedman_P = c(
    pc1_friedman$p.value,
    pc2_friedman$p.value
  )
)

# Save results
write.csv(
  trajectory_statistics,
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_PCA_Within_Patient_Statistics.csv"
  ),
  row.names = FALSE
)

write.csv(
  friedman_results,
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_PCA_Friedman_Tests.csv"
  ),
  row.names = FALSE
)

saveRDS(
  trajectory_statistics,
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_PCA_Within_Patient_Statistics.rds"
  )
)

saveRDS(
  friedman_results,
  file.path(
    pca_output_dir,
    "COVID_Longitudinal_PCA_Friedman_Tests.rds"
  )
)

# Print results
cat("Within-patient longitudinal PCA statistics completed.\n\n")

cat("Paired Wilcoxon results:\n")
print(
  data.frame(
    PC = trajectory_statistics$PC,
    Comparison = trajectory_statistics$Comparison,
    Mean_Change = round(trajectory_statistics$Mean_Change, 3),
    Median_Change = round(trajectory_statistics$Median_Change, 3),
    P_value = signif(trajectory_statistics$Wilcoxon_P, 4)
  )
)

cat("\nFriedman tests:\n")
print(
  data.frame(
    PC = friedman_results$PC,
    P_value = signif(friedman_results$Friedman_P, 4)
  )
)

cat("\nStatistical results saved successfully.\n")

# Prepare the repeated-patient cohort for gene-level longitudinal analysis

repeated_longitudinal <- pca_trajectory[
  pca_trajectory$Patient_ID %in% repeated_patients,
  ,
  drop = FALSE
]

# Keep only samples with valid longitudinal timepoints
repeated_longitudinal <- repeated_longitudinal[
  repeated_longitudinal$Timepoint %in% c("D0", "D3", "D7"),
  ,
  drop = FALSE
]

# Create ordered factors
repeated_longitudinal$Patient_ID <- factor(
  repeated_longitudinal$Patient_ID
)

repeated_longitudinal$Timepoint <- factor(
  repeated_longitudinal$Timepoint,
  levels = c("D0", "D3", "D7")
)

# Get the corresponding DESeq2 sample IDs
longitudinal_repeated_ids <- repeated_longitudinal$dds_sample_id

# Verify sample matching
matched_repeated_ids <- longitudinal_repeated_ids[
  longitudinal_repeated_ids %in% colnames(dds_longitudinal)
]

cat("Repeated patients:", 
    length(unique(repeated_longitudinal$Patient_ID)), "\n")

cat("Samples in repeated-patient cohort:",
    nrow(repeated_longitudinal), "\n")

cat("Matching DESeq2 samples:",
    length(matched_repeated_ids), "\n")

cat("Missing DESeq2 samples:",
    sum(
      !longitudinal_repeated_ids %in% colnames(dds_longitudinal)
    ), "\n")

cat("\nTimepoint distribution:\n")
print(table(repeated_longitudinal$Timepoint))

cat("\nPatients by number of timepoints:\n")
print(
  table(
    table(repeated_longitudinal$Patient_ID)
  )
)

cat("\nDesign variables:\n")
cat("Unique patients:",
    nlevels(repeated_longitudinal$Patient_ID), "\n")

cat("Timepoint levels:",
    paste(levels(repeated_longitudinal$Timepoint), collapse = ", "),
    "\n")

# Build the repeated-measures DESeq2 dataset

library(DESeq2)

# Use the correctly matched DESeq2 sample IDs
repeated_sample_ids <- repeated_longitudinal$dds_sample_id

# Extract counts for the repeated-patient cohort
counts_repeated <- counts(
  dds_longitudinal,
  normalized = FALSE
)[
  ,
  repeated_sample_ids,
  drop = FALSE
]

# Build sample metadata
metadata_repeated <- as.data.frame(
  colData(dds_longitudinal)[repeated_sample_ids, , drop = FALSE]
)

# Add longitudinal variables from the validated cohort
metadata_repeated$Patient_ID <- factor(
  repeated_longitudinal$Patient_ID[
    match(
      rownames(metadata_repeated),
      repeated_longitudinal$dds_sample_id
    )
  ]
)

metadata_repeated$Timepoint <- factor(
  repeated_longitudinal$Timepoint[
    match(
      rownames(metadata_repeated),
      repeated_longitudinal$dds_sample_id
    )
  ],
  levels = c("D0", "D3", "D7")
)

# Remove unused patient levels
metadata_repeated$Patient_ID <- droplevels(
  metadata_repeated$Patient_ID
)

# Create DESeq2 object
dds_repeated <- DESeqDataSetFromMatrix(
  countData = counts_repeated,
  colData = metadata_repeated,
  design = ~ Patient_ID + Timepoint
)

# Re-estimate size factors for the repeated-patient cohort
dds_repeated <- estimateSizeFactors(dds_repeated)

# Check design information
design_matrix_repeated <- model.matrix(
  ~ Patient_ID + Timepoint,
  data = metadata_repeated
)

cat("Repeated-measures DESeq2 dataset created.\n\n")

cat("Genes:", nrow(dds_repeated), "\n")
cat("Samples:", ncol(dds_repeated), "\n")
cat("Patients:", nlevels(metadata_repeated$Patient_ID), "\n\n")

cat("Timepoint distribution:\n")
print(table(metadata_repeated$Timepoint))

cat("\nDesign matrix dimensions:\n")
cat(
  nrow(design_matrix_repeated),
  "samples x",
  ncol(design_matrix_repeated),
  "coefficients\n"
)

cat("\nDesign matrix rank:\n")
cat(
  qr(design_matrix_repeated)$rank,
  "of",
  ncol(design_matrix_repeated),
  "\n"
)

cat("\nDesign matrix is full rank:",
    qr(design_matrix_repeated)$rank ==
      ncol(design_matrix_repeated),
    "\n")

cat("\nTimepoint coefficients expected:\n")
print(
  colnames(design_matrix_repeated)[
    grepl("Timepoint", colnames(design_matrix_repeated))
  ]
)

# Check limma availability

if (!requireNamespace("limma", quietly = TRUE)) {
  cat("limma is not installed.\n")
} else {
  cat("limma is installed.\n")
  cat("limma version:", as.character(packageVersion("limma")), "\n")
}

if (!requireNamespace("edgeR", quietly = TRUE)) {
  cat("edgeR is not installed.\n")
} else {
  cat("edgeR is installed.\n")
  cat("edgeR version:", as.character(packageVersion("edgeR")), "\n")
}

# Install edgeR if needed

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  "edgeR",
  ask = FALSE,
  update = FALSE
)

library(edgeR)

cat("edgeR installed successfully.\n")
cat("edgeR version:", as.character(packageVersion("edgeR")), "\n")


# Check whether the saved DESeq2 RDS file exists

rds_path <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/dds_primary_filtered.rds"

cat("File exists:", file.exists(rds_path), "\n")
cat("File path:", rds_path, "\n")

# Restore the saved DESeq2 object

library(DESeq2)

dds_longitudinal <- readRDS(rds_path)

cat("Object restored:", exists("dds_longitudinal"), "\n")
cat("Genes:", nrow(dds_longitudinal), "\n")
cat("Samples:", ncol(dds_longitudinal), "\n")

# Restore metadata and rebuild the longitudinal sample identifiers

metadata_longitudinal <- as.data.frame(
  colData(dds_longitudinal)
)

# Recreate the longitudinal cohort metadata from the saved metadata file
metadata_longitudinal_final <- readRDS(
  "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal_RNAseq_Final_Metadata_773_samples.rds"
)

# Keep only COVID-positive D0/D3/D7 samples
covid_longitudinal <- metadata_longitudinal_final[
  metadata_longitudinal_final$patient_category == "COVID+" &
    metadata_longitudinal_final$Timepoint %in% c("D0", "D3", "D7"),
  ,
  drop = FALSE
]

# Recreate the DESeq2 sample IDs
covid_longitudinal$dds_sample_id <- paste0(
  covid_longitudinal$Patient_ID,
  "_",
  covid_longitudinal$Timepoint
)

# Validate matching against the DESeq2 object
matched_ids <- covid_longitudinal$dds_sample_id %in%
  colnames(dds_longitudinal)

cat("COVID-positive longitudinal samples:",
    nrow(covid_longitudinal), "\n")

cat("Matched to dds_longitudinal:",
    sum(matched_ids), "\n")

cat("Missing from dds_longitudinal:",
    sum(!matched_ids), "\n")

cat("\nTimepoint distribution:\n")
print(table(covid_longitudinal$Timepoint))

# Prepare RNA-seq counts for limma-voom

library(edgeR)
library(limma)

counts_repeated <- counts(
  dds_longitudinal,
  normalized = FALSE
)[
  ,
  covid_longitudinal$dds_sample_id,
  drop = FALSE
]

metadata_repeated <- covid_longitudinal

# Convert variables to factors
metadata_repeated$Patient_ID <- factor(
  metadata_repeated$Patient_ID
)

metadata_repeated$Timepoint <- factor(
  metadata_repeated$Timepoint,
  levels = c("D0", "D3", "D7")
)

# Create edgeR object
dge_repeated <- DGEList(
  counts = counts_repeated
)

# Filter lowly expressed genes
keep_genes <- filterByExpr(
  dge_repeated,
  group = metadata_repeated$Timepoint
)

dge_repeated <- dge_repeated[
  keep_genes,
  ,
  keep.lib.sizes = FALSE
]

# TMM normalization
dge_repeated <- calcNormFactors(
  dge_repeated,
  method = "TMM"
)

cat("limma-voom preparation completed.\n\n")

cat("Genes before filtering:",
    nrow(counts_repeated), "\n")

cat("Genes retained after filtering:",
    nrow(dge_repeated), "\n")

cat("Samples:",
    ncol(dge_repeated), "\n")

cat("Patients:",
    nlevels(metadata_repeated$Patient_ID), "\n")

cat("\nTimepoint distribution:\n")
print(table(metadata_repeated$Timepoint))

cat("\nLibrary size summary:\n")
print(summary(dge_repeated$samples$lib.size))

cat("\nTMM normalization factors:\n")
print(summary(dge_repeated$samples$norm.factors))

# Build the repeated-measures design matrix

design_repeated <- model.matrix(
  ~ Patient_ID + Timepoint,
  data = metadata_repeated
)

cat("Design matrix dimensions:",
    nrow(design_repeated), "samples x",
    ncol(design_repeated), "coefficients\n\n")

cat("Design matrix columns:\n")
print(colnames(design_repeated))

cat("\nFull-rank check:\n")
cat(
  "Rank:",
  qr(design_repeated)$rank,
  "\nColumns:",
  ncol(design_repeated),
  "\n"
)

cat("\nTimepoint coefficients:\n")
print(
  grep(
    "^Timepoint",
    colnames(design_repeated),
    value = TRUE
  )
  
  # Apply voom transformation without opening a plotting device
  
  voom_repeated <- voom(
    dge_repeated,
    design = design_repeated,
    plot = FALSE
  )
  
  cat("voom transformation completed successfully.\n\n")
  
  cat("Voom expression matrix dimensions:",
      nrow(voom_repeated$E), "genes x",
      ncol(voom_repeated$E), "samples\n")
  
  cat("Genes:",
      nrow(voom_repeated$E), "\n")
  
  cat("Samples:",
      ncol(voom_repeated$E), "\n")
  
  cat("Missing values in voom matrix:",
      sum(!is.finite(voom_repeated$E)), "\n") 
  # Apply voom transformation without opening a plotting device
  
  voom_repeated <- voom(
    dge_repeated,
    design = design_repeated,
    plot = FALSE
  )
  
  cat("voom transformation completed successfully.\n\n")
  
  cat("Voom expression matrix dimensions:",
      nrow(voom_repeated$E), "genes x",
      ncol(voom_repeated$E), "samples\n")
  
  cat("Genes:",
      nrow(voom_repeated$E), "\n")
  
  cat("Samples:",
      ncol(voom_repeated$E), "\n")
  
  cat("Missing values in voom matrix:",
      sum(!is.finite(voom_repeated$E)), "\n") 

  # Prepare complete-patient paired data for fast longitudinal analysis
  
  library(edgeR)
  library(limma)
  
  # Identify patients with exactly one sample at each timepoint
  patient_timepoint_counts <- table(
    metadata_repeated$Patient_ID,
    metadata_repeated$Timepoint
  )
  
  complete_patient_ids <- rownames(
    patient_timepoint_counts[
      patient_timepoint_counts[, "D0"] == 1 &
        patient_timepoint_counts[, "D3"] == 1 &
        patient_timepoint_counts[, "D7"] == 1,
      ,
      drop = FALSE
    ]
  )
  
  cat("Complete patients:",
      length(complete_patient_ids), "\n")
  
  # Keep only complete patients
  complete_metadata <- metadata_repeated[
    metadata_repeated$Patient_ID %in% complete_patient_ids,
    ,
    drop = FALSE
  ]
  
  # Create sample IDs
  complete_metadata$dds_sample_id <- paste0(
    complete_metadata$Patient_ID,
    "_",
    complete_metadata$Timepoint
  )
  
  # Order patients and timepoints consistently
  complete_metadata$Patient_ID <- factor(
    complete_metadata$Patient_ID,
    levels = sort(complete_patient_ids)
  )
  
  complete_metadata$Timepoint <- factor(
    complete_metadata$Timepoint,
    levels = c("D0", "D3", "D7")
  )
  
  complete_metadata <- complete_metadata[
    order(
      complete_metadata$Patient_ID,
      complete_metadata$Timepoint
    ),
    ,
    drop = FALSE
  ]
  
  # Extract normalized log-CPM values
  logCPM_complete <- cpm(
    dge_repeated,
    log = TRUE,
    prior.count = 1
  )[
    ,
    complete_metadata$dds_sample_id,
    drop = FALSE
  ]
  
  # Create patient-by-timepoint sample index
  sample_index <- matrix(
    NA_character_,
    nrow = length(complete_patient_ids),
    ncol = 3,
    dimnames = list(
      sort(complete_patient_ids),
      c("D0", "D3", "D7")
    )
  )
  
  for (tp in c("D0", "D3", "D7")) {
    idx <- complete_metadata$Timepoint == tp
    sample_index[
      as.character(complete_metadata$Patient_ID[idx]),
      tp
    ] <- complete_metadata$dds_sample_id[idx]
  }
  
  # Create within-patient expression differences
  logCPM_D0 <- logCPM_complete[, sample_index[, "D0"], drop = FALSE]
  logCPM_D3 <- logCPM_complete[, sample_index[, "D3"], drop = FALSE]
  logCPM_D7 <- logCPM_complete[, sample_index[, "D7"], drop = FALSE]
  
  diff_D3_D0 <- logCPM_D3 - logCPM_D0
  diff_D7_D0 <- logCPM_D7 - logCPM_D0
  diff_D7_D3 <- logCPM_D7 - logCPM_D3
  
  cat("\nComplete-patient timepoint distribution:\n")
  print(table(complete_metadata$Timepoint))
  
  cat("\nExpression matrix dimensions:\n")
  cat("Genes:",
      nrow(logCPM_complete), "\n")
  cat("Complete patients:",
      ncol(logCPM_complete) / 3, "\n")
  cat("Samples:",
      ncol(logCPM_complete), "\n")
  
  cat("\nDifference matrix dimensions:\n")
  cat("D3-D0:",
      nrow(diff_D3_D0), "genes x",
      ncol(diff_D3_D0), "patients\n")
  
  cat("D7-D0:",
      nrow(diff_D7_D0), "genes x",
      ncol(diff_D7_D0), "patients\n")
  
  cat("D7-D3:",
      nrow(diff_D7_D3), "genes x",
      ncol(diff_D7_D3), "patients\n")
  
  cat("\nNon-finite values:\n")
  cat("D3-D0:",
      sum(!is.finite(diff_D3_D0)), "\n")
  cat("D7-D0:",
      sum(!is.finite(diff_D7_D0)), "\n")
  cat("D7-D3:",
      sum(!is.finite(diff_D7_D3)), "\n")
  
  
  # Perform paired gene-level differential analysis: D3 vs D0
  
  p_D3_D0 <- rowMeans(diff_D3_D0)
  
  t_D3_D0 <- apply(
    diff_D3_D0,
    1,
    function(x) {
      t.test(x, mu = 0)$statistic
    }
  )
  
  pvalue_D3_D0 <- apply(
    diff_D3_D0,
    1,
    function(x) {
      t.test(x, mu = 0)$p.value
    }
  )
  
  fdr_D3_D0 <- p.adjust(
    pvalue_D3_D0,
    method = "BH"
  )
  
  results_D3_D0 <- data.frame(
    Gene = rownames(diff_D3_D0),
    logFC = p_D3_D0,
    t_statistic = as.numeric(t_D3_D0),
    P.Value = pvalue_D3_D0,
    FDR = fdr_D3_D0,
    stringsAsFactors = FALSE
  )
  
  results_D3_D0 <- results_D3_D0[
    order(results_D3_D0$FDR, -abs(results_D3_D0$logFC)),
  ]
  
  significant_D3_D0 <- subset(
    results_D3_D0,
    FDR < 0.05 & abs(logFC) >= 1
  )
  
  cat("D3 vs D0 paired analysis completed.\n\n")
  
  cat("Genes tested:",
      nrow(results_D3_D0), "\n")
  
  cat("Significant genes (FDR < 0.05 and |logFC| >= 1):",
      nrow(significant_D3_D0), "\n")
  
  cat("Upregulated:",
      sum(significant_D3_D0$logFC > 0), "\n")
  
  cat("Downregulated:",
      sum(significant_D3_D0$logFC < 0), "\n")
  
  cat("\nTop 10 genes:\n")
  print(
    head(
      results_D3_D0,
      10
    )
  )
  
  
  # Save D3 vs D0 longitudinal differential expression results
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write.csv(
    results_D3_D0,
    file.path(
      output_dir,
      "COVID_Longitudinal_D3_vs_D0_all_genes.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    results_D3_D0,
    file.path(
      output_dir,
      "COVID_Longitudinal_D3_vs_D0_all_genes.rds"
    )
  )
  
  write.csv(
    significant_D3_D0,
    file.path(
      output_dir,
      "COVID_Longitudinal_D3_vs_D0_significant_DEGs.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    significant_D3_D0,
    file.path(
      output_dir,
      "COVID_Longitudinal_D3_vs_D0_significant_DEGs.rds"
    )
  )
  
  write.csv(
    subset(significant_D3_D0, logFC > 0),
    file.path(
      output_dir,
      "COVID_Longitudinal_D3_vs_D0_upregulated.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    subset(significant_D3_D0, logFC < 0),
    file.path(
      output_dir,
      "COVID_Longitudinal_D3_vs_D0_downregulated.csv"
    ),
    row.names = FALSE
  )
  
  cat("D3 vs D0 results saved successfully.\n")
  cat("Output directory:\n", output_dir, "\n\n")
  
  cat("Files saved:\n")
  cat("- All genes: CSV + RDS\n")
  cat("- Significant DEGs: CSV + RDS\n")
  cat("- Upregulated DEGs: CSV\n")
  cat("- Downregulated DEGs: CSV\n")
  p_D7_D0 <- rowMeans(diff_D7_D0)
  
  t_D7_D0 <- apply(
    diff_D7_D0,
    1,
    function(x) t.test(x, mu = 0)$statistic
  )
  
  pvalue_D7_D0 <- apply(
    diff_D7_D0,
    1,
    function(x) t.test(x, mu = 0)$p.value
  )
  
  fdr_D7_D0 <- p.adjust(pvalue_D7_D0, method = "BH")
  
  results_D7_D0 <- data.frame(
    Gene = rownames(diff_D7_D0),
    logFC = p_D7_D0,
    t = as.numeric(t_D7_D0),
    PValue = pvalue_D7_D0,
    FDR = fdr_D7_D0
  )
  
  results_D7_D0 <- results_D7_D0[
    order(results_D7_D0$FDR, -abs(results_D7_D0$logFC)),
  ]
  
  significant_D7_D0 <- subset(
    results_D7_D0,
    FDR < 0.05 & abs(logFC) >= 1
  )
  
  cat("D7 vs D0 analysis completed.\n")
  cat("Genes tested:", nrow(results_D7_D0), "\n")
  cat("Significant DEGs:", nrow(significant_D7_D0), "\n")
  cat("Upregulated:", sum(significant_D7_D0$logFC > 0), "\n")
  cat("Downregulated:", sum(significant_D7_D0$logFC < 0), "\n\n")
  
  cat("Top 10 genes:\n")
  print(head(results_D7_D0, 10))  
 
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    results_D7_D0,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_all_genes.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    results_D7_D0,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_all_genes.rds")
  )
  
  write.csv(
    significant_D7_D0,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_significant_DEGs.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    significant_D7_D0,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_significant_DEGs.rds")
  )
  
  write.csv(
    subset(significant_D7_D0, logFC > 0),
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    subset(significant_D7_D0, logFC < 0),
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_downregulated.csv"),
    row.names = FALSE
  )
  
  cat("D7 vs D0 results saved successfully.\n")
  cat("Output directory:\n", output_dir, "\n\n")
  cat("Files saved:\n")
  cat("- All genes: CSV + RDS\n")
  cat("- Significant DEGs: CSV + RDS\n")
  cat("- Upregulated DEGs: CSV\n")
  cat("- Downregulated DEGs: CSV\n")

  p_D7_D3 <- rowMeans(diff_D7_D3)
  
  t_D7_D3 <- apply(
    diff_D7_D3,
    1,
    function(x) t.test(x, mu = 0)$statistic
  )
  
  pvalue_D7_D3 <- apply(
    diff_D7_D3,
    1,
    function(x) t.test(x, mu = 0)$p.value
  )
  
  fdr_D7_D3 <- p.adjust(pvalue_D7_D3, method = "BH")
  
  results_D7_D3 <- data.frame(
    Gene = rownames(diff_D7_D3),
    logFC = p_D7_D3,
    t = as.numeric(t_D7_D3),
    PValue = pvalue_D7_D3,
    FDR = fdr_D7_D3
  )
  
  results_D7_D3 <- results_D7_D3[
    order(results_D7_D3$FDR, -abs(results_D7_D3$logFC)),
  ]
  
  significant_D7_D3 <- subset(
    results_D7_D3,
    FDR < 0.05 & abs(logFC) >= 1
  )
  
  cat("D7 vs D3 analysis completed.\n")
  cat("Genes tested:", nrow(results_D7_D3), "\n")
  cat("Significant DEGs:", nrow(significant_D7_D3), "\n")
  cat("Upregulated:", sum(significant_D7_D3$logFC > 0), "\n")
  cat("Downregulated:", sum(significant_D7_D3$logFC < 0), "\n\n")
  
  cat("Top 10 genes:\n")
  print(head(results_D7_D3, 10))
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    results_D7_D3,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_all_genes.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    results_D7_D3,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_all_genes.rds")
  )
  
  write.csv(
    significant_D7_D3,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_significant_DEGs.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    significant_D7_D3,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_significant_DEGs.rds")
  )
  
  write.csv(
    subset(significant_D7_D3, logFC > 0),
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    subset(significant_D7_D3, logFC < 0),
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_downregulated.csv"),
    row.names = FALSE
  )
  
  cat("D7 vs D3 results saved successfully.\n")
  cat("Output directory:\n", output_dir, "\n\n")
  cat("Files saved:\n")
  cat("- All genes: CSV + RDS\n")
  cat("- Significant DEGs: CSV + RDS\n")
  cat("- Upregulated DEGs: CSV\n")
  cat("- Downregulated DEGs: CSV\n")
  colnames(results_D3_D0_annotated)
  
  cat("\nNumber of rows:", nrow(results_D3_D0_annotated), "\n")
  
  cat(
    "Mapped symbols:",
    sum(!is.na(results_D3_D0_annotated$SYMBOL) &
          results_D3_D0_annotated$SYMBOL != ""),
    "\n"
  )
  
  cat(
    "Unmapped symbols:",
    sum(is.na(results_D3_D0_annotated$SYMBOL) |
          results_D3_D0_annotated$SYMBOL == ""),
    "\n"
  )
  
  cat("\nFirst rows:\n")
  print(head(results_D3_D0_annotated))
  
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  significant_D3_D0_annotated <- subset(
    results_D3_D0_annotated,
    FDR < 0.05 & abs(logFC) >= 1
  )
  
  upregulated_D3_D0_annotated <- subset(
    significant_D3_D0_annotated,
    logFC > 0
  )
  
  downregulated_D3_D0_annotated <- subset(
    significant_D3_D0_annotated,
    logFC < 0
  )
  
  write.csv(
    results_D3_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_all_genes_annotated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    results_D3_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_all_genes_annotated.rds")
  )
  
  write.csv(
    significant_D3_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_significant_DEGs_annotated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    significant_D3_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_significant_DEGs_annotated.rds")
  )
  
  write.csv(
    upregulated_D3_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_upregulated_annotated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    downregulated_D3_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_downregulated_annotated.csv"),
    row.names = FALSE
  )
  
  cat("D3 vs D0 annotated results saved successfully.\n")
  cat("Total genes:", nrow(results_D3_D0_annotated), "\n")
  cat("Mapped symbols:", sum(!is.na(results_D3_D0_annotated$SYMBOL) &
                               results_D3_D0_annotated$SYMBOL != ""), "\n")
  cat("Significant DEGs:", nrow(significant_D3_D0_annotated), "\n")
  cat("Upregulated:", nrow(upregulated_D3_D0_annotated), "\n")
  cat("Downregulated:", nrow(downregulated_D3_D0_annotated), "\n")  

  
  results_D7_D0_annotated <- annotate_longitudinal_results(results_D7_D0)
  
  cat("Annotation completed for D7 vs D0.\n")
  cat("Genes:", nrow(results_D7_D0_annotated), "\n")
  cat(
    "Mapped symbols:",
    sum(!is.na(results_D7_D0_annotated$SYMBOL) &
          results_D7_D0_annotated$SYMBOL != ""),
    "\n"
  )
  cat(
    "Unmapped symbols:",
    sum(is.na(results_D7_D0_annotated$SYMBOL) |
          results_D7_D0_annotated$SYMBOL == ""),
    "\n\n"
  )
  
  cat("Top 10 annotated genes:\n")
  print(
    head(
      results_D7_D0_annotated[
        , c("Gene", "SYMBOL", "GENENAME", "logFC", "P.Value", "FDR")
      ],
      10
    )
  )
  
  
  colnames(results_D7_D0_annotated)
  
  cat("\nNumber of rows:", nrow(results_D7_D0_annotated), "\n")
  
  cat(
    "Mapped symbols:",
    sum(!is.na(results_D7_D0_annotated$SYMBOL) &
          results_D7_D0_annotated$SYMBOL != ""),
    "\n"
  )
  
  cat(
    "Unmapped symbols:",
    sum(is.na(results_D7_D0_annotated$SYMBOL) |
          results_D7_D0_annotated$SYMBOL == ""),
    "\n\n"
  )
  
  cat("First 10 annotated genes:\n")
  print(head(results_D7_D0_annotated, 10))
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  significant_D7_D0_annotated <- subset(
    results_D7_D0_annotated,
    FDR < 0.05 & abs(logFC) >= 1
  )
  
  upregulated_D7_D0_annotated <- subset(
    significant_D7_D0_annotated,
    logFC > 0
  )
  
  downregulated_D7_D0_annotated <- subset(
    significant_D7_D0_annotated,
    logFC < 0
  )
  
  write.csv(
    results_D7_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_all_genes_annotated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    results_D7_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_all_genes_annotated.rds")
  )
  
  write.csv(
    significant_D7_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_significant_DEGs_annotated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    significant_D7_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_significant_DEGs_annotated.rds")
  )
  
  write.csv(
    upregulated_D7_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_upregulated_annotated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    downregulated_D7_D0_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_downregulated_annotated.csv"),
    row.names = FALSE
  )
  
  cat("D7 vs D0 annotated results saved successfully.\n")
  cat("Total genes:", nrow(results_D7_D0_annotated), "\n")
  cat("Mapped symbols:", sum(!is.na(results_D7_D0_annotated$SYMBOL) &
                               results_D7_D0_annotated$SYMBOL != ""), "\n")
  cat("Significant DEGs:", nrow(significant_D7_D0_annotated), "\n")
  cat("Upregulated:", nrow(upregulated_D7_D0_annotated), "\n")
  cat("Downregulated:", nrow(downregulated_D7_D0_annotated), "\n")
  
  
  results_D7_D3_annotated <- annotate_longitudinal_results(results_D7_D3)
  
  cat("Annotation completed for D7 vs D3.\n")
  cat("Genes:", nrow(results_D7_D3_annotated), "\n")
  
  cat(
    "Mapped symbols:",
    sum(!is.na(results_D7_D3_annotated$SYMBOL) &
          results_D7_D3_annotated$SYMBOL != ""),
    "\n"
  )
  
  cat(
    "Unmapped symbols:",
    sum(is.na(results_D7_D3_annotated$SYMBOL) |
          results_D7_D3_annotated$SYMBOL == ""),
    "\n\n"
  )
  
  cat("First 10 annotated genes:\n")
  print(head(results_D7_D3_annotated, 10))
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  significant_D7_D3_annotated <- subset(
    results_D7_D3_annotated,
    FDR < 0.05 & abs(logFC) >= 1
  )
  
  upregulated_D7_D3_annotated <- subset(
    significant_D7_D3_annotated,
    logFC > 0
  )
  
  downregulated_D7_D3_annotated <- subset(
    significant_D7_D3_annotated,
    logFC < 0
  )
  
  write.csv(
    results_D7_D3_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_all_genes_annotated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    results_D7_D3_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_all_genes_annotated.rds")
  )
  
  write.csv(
    significant_D7_D3_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_significant_DEGs_annotated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    significant_D7_D3_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_significant_DEGs_annotated.rds")
  )
  
  write.csv(
    upregulated_D7_D3_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_upregulated_annotated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    downregulated_D7_D3_annotated,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_downregulated_annotated.csv"),
    row.names = FALSE
  )
  
  cat("D7 vs D3 annotated results saved successfully.\n")
  cat("Total genes:", nrow(results_D7_D3_annotated), "\n")
  cat("Mapped symbols:", sum(!is.na(results_D7_D3_annotated$SYMBOL) &
                               results_D7_D3_annotated$SYMBOL != ""), "\n")
  cat("Significant DEGs:", nrow(significant_D7_D3_annotated), "\n")
  cat("Upregulated:", nrow(upregulated_D7_D3_annotated), "\n")
  cat("Downregulated:", nrow(downregulated_D7_D3_annotated), "\n")  

  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  sig_D3_D0_up_symbols <- unique(
    na.omit(upregulated_D3_D0_annotated$SYMBOL)
  )
  
  sig_D3_D0_down_symbols <- unique(
    na.omit(downregulated_D3_D0_annotated$SYMBOL)
  )
  
  sig_D3_D0_all_symbols <- unique(
    na.omit(significant_D3_D0_annotated$SYMBOL)
  )
  
  entrez_D3_D0_up <- bitr(
    sig_D3_D0_up_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  entrez_D3_D0_down <- bitr(
    sig_D3_D0_down_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  entrez_D3_D0_all <- bitr(
    sig_D3_D0_all_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  entrez_D3_D0_up <- unique(entrez_D3_D0_up$ENTREZID)
  entrez_D3_D0_down <- unique(entrez_D3_D0_down$ENTREZID)
  entrez_D3_D0_all <- unique(entrez_D3_D0_all$ENTREZID)
  
  ego_D3_D0_all <- enrichGO(
    gene = entrez_D3_D0_all,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_D3_D0_up <- enrichGO(
    gene = entrez_D3_D0_up,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_D3_D0_down <- enrichGO(
    gene = entrez_D3_D0_down,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_D3_D0_all_df <- as.data.frame(ego_D3_D0_all)
  ego_D3_D0_up_df <- as.data.frame(ego_D3_D0_up)
  ego_D3_D0_down_df <- as.data.frame(ego_D3_D0_down)
  
  cat("D3 vs D0 GO-BP enrichment completed.\n\n")
  
  cat("Input genes:\n")
  cat("All:", length(entrez_D3_D0_all), "\n")
  cat("Up:", length(entrez_D3_D0_up), "\n")
  cat("Down:", length(entrez_D3_D0_down), "\n\n")
  
  cat("Significant GO-BP pathways:\n")
  cat("All DEGs:", nrow(ego_D3_D0_all_df), "\n")
  cat("Upregulated:", nrow(ego_D3_D0_up_df), "\n")
  cat("Downregulated:", nrow(ego_D3_D0_down_df), "\n\n")
  
  cat("Top 10 GO-BP pathways - All DEGs:\n")
  print(
    head(
      ego_D3_D0_all_df[
        order(ego_D3_D0_all_df$p.adjust),
        c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "geneID")
      ],
      10
    )
  )
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    ego_D3_D0_all_df,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_GO_BP_all_DEGs.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ego_D3_D0_up_df,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_GO_BP_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ego_D3_D0_down_df,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_GO_BP_downregulated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    ego_D3_D0_all,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_GO_BP_all_DEGs.rds")
  )
  
  saveRDS(
    ego_D3_D0_up,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_GO_BP_upregulated.rds")
  )
  
  saveRDS(
    ego_D3_D0_down,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_GO_BP_downregulated.rds")
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_D3_D0_all),
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_GO_BP_input_all_Entrez.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_D3_D0_up),
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_GO_BP_input_up_Entrez.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_D3_D0_down),
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_GO_BP_input_down_Entrez.csv"),
    row.names = FALSE
  )
  
  cat("D3 vs D0 GO-BP results saved successfully.\n")
  cat("All DEGs pathways:", nrow(ego_D3_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(ego_D3_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(ego_D3_D0_down_df), "\n")
  cat("All input Entrez IDs:", length(entrez_D3_D0_all), "\n")
  cat("Up input Entrez IDs:", length(entrez_D3_D0_up), "\n")
  cat("Down input Entrez IDs:", length(entrez_D3_D0_down), "\n") 
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  ekegg_D3_D0_all <- enrichKEGG(
    gene = entrez_D3_D0_all,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  ekegg_D3_D0_up <- enrichKEGG(
    gene = entrez_D3_D0_up,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  ekegg_D3_D0_down <- enrichKEGG(
    gene = entrez_D3_D0_down,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  ekegg_D3_D0_all_df <- as.data.frame(ekegg_D3_D0_all)
  ekegg_D3_D0_up_df <- as.data.frame(ekegg_D3_D0_up)
  ekegg_D3_D0_down_df <- as.data.frame(ekegg_D3_D0_down)
  
  cat("D3 vs D0 KEGG enrichment completed.\n")
  cat("All DEGs pathways:", nrow(ekegg_D3_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(ekegg_D3_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(ekegg_D3_D0_down_df), "\n")
  
  cat("\nTop KEGG pathways - All DEGs:\n")
  print(
    ekegg_D3_D0_all_df[
      order(ekegg_D3_D0_all_df$p.adjust),
      c("ID", "Description", "GeneRatio", "p.adjust")
    ][1:min(10, nrow(ekegg_D3_D0_all_df)), ]
  )
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  write.csv(
    ekegg_D3_D0_all_df,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_KEGG_all_DEGs.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ekegg_D3_D0_up_df,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_KEGG_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ekegg_D3_D0_down_df,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_KEGG_downregulated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    ekegg_D3_D0_all,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_KEGG_all_DEGs.rds")
  )
  
  saveRDS(
    ekegg_D3_D0_up,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_KEGG_upregulated.rds")
  )
  
  saveRDS(
    ekegg_D3_D0_down,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_KEGG_downregulated.rds")
  )
  
  cat("D3 vs D0 KEGG results saved successfully.\n")
  cat("All DEGs pathways:", nrow(ekegg_D3_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(ekegg_D3_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(ekegg_D3_D0_down_df), "\n") 
  
  library(ReactomePA)
  library(org.Hs.eg.db)
  
  reactome_D3_D0_all <- enrichPathway(
    gene = entrez_D3_D0_all,
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  reactome_D3_D0_up <- enrichPathway(
    gene = entrez_D3_D0_up,
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  reactome_D3_D0_down <- enrichPathway(
    gene = entrez_D3_D0_down,
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  reactome_D3_D0_all_df <- as.data.frame(reactome_D3_D0_all)
  reactome_D3_D0_up_df <- as.data.frame(reactome_D3_D0_up)
  reactome_D3_D0_down_df <- as.data.frame(reactome_D3_D0_down)
  
  cat("D3 vs D0 Reactome enrichment completed.\n")
  cat("All DEGs pathways:", nrow(reactome_D3_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(reactome_D3_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(reactome_D3_D0_down_df), "\n")
  
  cat("\nTop Reactome pathways - All DEGs:\n")
  
  print(
    reactome_D3_D0_all_df[
      order(reactome_D3_D0_all_df$p.adjust),
      c("ID", "Description", "GeneRatio", "p.adjust")
    ][1:min(10, nrow(reactome_D3_D0_all_df)), ]
  )  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  write.csv(
    reactome_D3_D0_all_df,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_Reactome_all_DEGs.csv"),
    row.names = FALSE
  )
  
  write.csv(
    reactome_D3_D0_up_df,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_Reactome_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    reactome_D3_D0_down_df,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_Reactome_downregulated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    reactome_D3_D0_all,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_Reactome_all_DEGs.rds")
  )
  
  saveRDS(
    reactome_D3_D0_up,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_Reactome_upregulated.rds")
  )
  
  saveRDS(
    reactome_D3_D0_down,
    file.path(output_dir, "COVID_Longitudinal_D3_vs_D0_Reactome_downregulated.rds")
  )
  
  cat("D3 vs D0 Reactome results saved successfully.\n")
  cat("All DEGs pathways:", nrow(reactome_D3_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(reactome_D3_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(reactome_D3_D0_down_df), "\n")
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  sig_D7_D0 <- readRDS(
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_significant_DEGs.rds")
  )
  
  sig_D7_D0_up <- read.csv(
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_upregulated.csv"),
    stringsAsFactors = FALSE
  )
  
  sig_D7_D0_down <- read.csv(
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_downregulated.csv"),
    stringsAsFactors = FALSE
  )
  
  entrez_D7_D0_all <- bitr(
    sig_D7_D0$Gene,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )$ENTREZID
  
  entrez_D7_D0_up <- bitr(
    sig_D7_D0_up$Gene,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )$ENTREZID
  
  entrez_D7_D0_down <- bitr(
    sig_D7_D0_down$Gene,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )$ENTREZID
  
  entrez_D7_D0_all <- unique(entrez_D7_D0_all)
  entrez_D7_D0_up <- unique(entrez_D7_D0_up)
  entrez_D7_D0_down <- unique(entrez_D7_D0_down)
  
  cat("D7 vs D0 Entrez conversion completed.\n")
  cat("Significant DEGs:", nrow(sig_D7_D0), "\n")
  cat("Upregulated DEGs:", nrow(sig_D7_D0_up), "\n")
  cat("Downregulated DEGs:", nrow(sig_D7_D0_down), "\n")
  cat("All Entrez IDs:", length(entrez_D7_D0_all), "\n")
  cat("Up Entrez IDs:", length(entrez_D7_D0_up), "\n")
  cat("Down Entrez IDs:", length(entrez_D7_D0_down), "\n")
  
  
  ego_D7_D0_all <- enrichGO(
    gene = entrez_D7_D0_all,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_D7_D0_up <- enrichGO(
    gene = entrez_D7_D0_up,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_D7_D0_down <- enrichGO(
    gene = entrez_D7_D0_down,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_D7_D0_all_df <- as.data.frame(ego_D7_D0_all)
  ego_D7_D0_up_df <- as.data.frame(ego_D7_D0_up)
  ego_D7_D0_down_df <- as.data.frame(ego_D7_D0_down)
  
  cat("D7 vs D0 GO-BP enrichment completed.\n")
  cat("All DEGs pathways:", nrow(ego_D7_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(ego_D7_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(ego_D7_D0_down_df), "\n")
  
  cat("\nTop GO-BP pathways - All DEGs:\n")
  
  print(
    ego_D7_D0_all_df[
      order(ego_D7_D0_all_df$p.adjust),
      c("ID", "Description", "GeneRatio", "p.adjust")
    ][1:min(10, nrow(ego_D7_D0_all_df)), ]
  )
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  write.csv(
    ego_D7_D0_all_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_GO_BP_all_DEGs.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ego_D7_D0_up_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_GO_BP_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ego_D7_D0_down_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_GO_BP_downregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_D7_D0_all),
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_GO_BP_input_all_Entrez.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_D7_D0_up),
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_GO_BP_input_up_Entrez.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_D7_D0_down),
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_GO_BP_input_down_Entrez.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    ego_D7_D0_all,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_GO_BP_all_DEGs.rds")
  )
  
  saveRDS(
    ego_D7_D0_up,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_GO_BP_upregulated.rds")
  )
  
  saveRDS(
    ego_D7_D0_down,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_GO_BP_downregulated.rds")
  )
  
  cat("D7 vs D0 GO-BP results saved successfully.\n")
  cat("All DEGs pathways:", nrow(ego_D7_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(ego_D7_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(ego_D7_D0_down_df), "\n")
  cat("All input Entrez IDs:", length(entrez_D7_D0_all), "\n")
  cat("Up input Entrez IDs:", length(entrez_D7_D0_up), "\n")
  cat("Down input Entrez IDs:", length(entrez_D7_D0_down), "\n")  
  ekegg_D7_D0_all <- enrichKEGG(
    gene = entrez_D7_D0_all,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  ekegg_D7_D0_up <- enrichKEGG(
    gene = entrez_D7_D0_up,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  ekegg_D7_D0_down <- enrichKEGG(
    gene = entrez_D7_D0_down,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  ekegg_D7_D0_all_df <- as.data.frame(ekegg_D7_D0_all)
  ekegg_D7_D0_up_df <- as.data.frame(ekegg_D7_D0_up)
  ekegg_D7_D0_down_df <- as.data.frame(ekegg_D7_D0_down)
  
  cat("D7 vs D0 KEGG enrichment completed.\n")
  cat("All DEGs pathways:", nrow(ekegg_D7_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(ekegg_D7_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(ekegg_D7_D0_down_df), "\n")
  
  cat("\nTop KEGG pathways - All DEGs:\n")
  
  print(
    ekegg_D7_D0_all_df[
      order(ekegg_D7_D0_all_df$p.adjust),
      c("ID", "Description", "GeneRatio", "p.adjust")
    ][1:min(10, nrow(ekegg_D7_D0_all_df)), ]
  )
  
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  write.csv(
    ekegg_D7_D0_all_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_KEGG_all_DEGs.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ekegg_D7_D0_up_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_KEGG_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ekegg_D7_D0_down_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_KEGG_downregulated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    ekegg_D7_D0_all,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_KEGG_all_DEGs.rds")
  )
  
  saveRDS(
    ekegg_D7_D0_up,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_KEGG_upregulated.rds")
  )
  
  saveRDS(
    ekegg_D7_D0_down,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_KEGG_downregulated.rds")
  )
  
  cat("D7 vs D0 KEGG results saved successfully.\n")
  cat("All DEGs pathways:", nrow(ekegg_D7_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(ekegg_D7_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(ekegg_D7_D0_down_df), "\n") 
  
  reactome_D7_D0_all <- enrichPathway(
    gene = entrez_D7_D0_all,
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  reactome_D7_D0_up <- enrichPathway(
    gene = entrez_D7_D0_up,
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  reactome_D7_D0_down <- enrichPathway(
    gene = entrez_D7_D0_down,
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  reactome_D7_D0_all_df <- as.data.frame(reactome_D7_D0_all)
  reactome_D7_D0_up_df <- as.data.frame(reactome_D7_D0_up)
  reactome_D7_D0_down_df <- as.data.frame(reactome_D7_D0_down)
  
  cat("D7 vs D0 Reactome enrichment completed.\n")
  cat("All DEGs pathways:", nrow(reactome_D7_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(reactome_D7_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(reactome_D7_D0_down_df), "\n")
  
  cat("\nTop Reactome pathways - All DEGs:\n")
  
  print(
    reactome_D7_D0_all_df[
      order(reactome_D7_D0_all_df$p.adjust),
      c("ID", "Description", "GeneRatio", "p.adjust")
    ][1:min(10, nrow(reactome_D7_D0_all_df)), ]
  ) 
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  write.csv(
    reactome_D7_D0_all_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_Reactome_all_DEGs.csv"),
    row.names = FALSE
  )
  
  write.csv(
    reactome_D7_D0_up_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_Reactome_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    reactome_D7_D0_down_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_Reactome_downregulated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    reactome_D7_D0_all,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_Reactome_all_DEGs.rds")
  )
  
  saveRDS(
    reactome_D7_D0_up,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_Reactome_upregulated.rds")
  )
  
  saveRDS(
    reactome_D7_D0_down,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D0_Reactome_downregulated.rds")
  )
  
  cat("D7 vs D0 Reactome results saved successfully.\n")
  cat("All DEGs pathways:", nrow(reactome_D7_D0_all_df), "\n")
  cat("Upregulated pathways:", nrow(reactome_D7_D0_up_df), "\n")
  cat("Downregulated pathways:", nrow(reactome_D7_D0_down_df), "\n")  
  
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  convert_to_entrez <- function(df) {
    genes <- unique(na.omit(df$Gene))
    
    annotation <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = genes,
      keytype = "ENSEMBL",
      columns = c("ENSEMBL", "ENTREZID")
    )
    
    annotation <- annotation[
      !is.na(annotation$ENTREZID) &
        !duplicated(annotation$ENSEMBL),
    ]
    
    unique(annotation$ENTREZID)
  }
  
  entrez_D7_D3_all <- convert_to_entrez(sig_D7_D3)
  entrez_D7_D3_up <- convert_to_entrez(up_D7_D3)
  entrez_D7_D3_down <- convert_to_entrez(down_D7_D3)
  
  cat("D7 vs D3 Entrez mapping:\n")
  cat("Significant DEGs:", nrow(sig_D7_D3),
      "->", length(entrez_D7_D3_all), "Entrez IDs\n")
  cat("Upregulated:", nrow(up_D7_D3),
      "->", length(entrez_D7_D3_up), "Entrez IDs\n")
  cat("Downregulated:", nrow(down_D7_D3),
      "->", length(entrez_D7_D3_down), "Entrez IDs\n")
  
  cat("\nValidation:\n")
  cat("All IDs unique:", length(entrez_D7_D3_all) == length(unique(entrez_D7_D3_all)), "\n")
  cat("Up IDs unique:", length(entrez_D7_D3_up) == length(unique(entrez_D7_D3_up)), "\n")
  cat("Down IDs unique:", length(entrez_D7_D3_down) == length(unique(entrez_D7_D3_down)), "\n")
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  ego_D7_D3_all <- enrichGO(
    gene = entrez_D7_D3_all,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  ego_D7_D3_up <- enrichGO(
    gene = entrez_D7_D3_up,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  ego_D7_D3_down <- enrichGO(
    gene = entrez_D7_D3_down,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  ego_D7_D3_all_df <- as.data.frame(ego_D7_D3_all)
  ego_D7_D3_up_df <- as.data.frame(ego_D7_D3_up)
  ego_D7_D3_down_df <- as.data.frame(ego_D7_D3_down)
  
  write.csv(
    ego_D7_D3_all_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_GO_BP_all_DEGs.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ego_D7_D3_up_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_GO_BP_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ego_D7_D3_down_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_GO_BP_downregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_D7_D3_all),
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_GO_BP_input_all_Entrez.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_D7_D3_up),
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_GO_BP_input_up_Entrez.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(ENTREZID = entrez_D7_D3_down),
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_GO_BP_input_down_Entrez.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    ego_D7_D3_all,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_GO_BP_all_DEGs.rds")
  )
  
  saveRDS(
    ego_D7_D3_up,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_GO_BP_upregulated.rds")
  )
  
  saveRDS(
    ego_D7_D3_down,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_GO_BP_downregulated.rds")
  )
  
  cat("D7 vs D3 GO-BP enrichment completed and saved successfully.\n")
  cat("All DEGs pathways:", nrow(ego_D7_D3_all_df), "\n")
  cat("Upregulated pathways:", nrow(ego_D7_D3_up_df), "\n")
  cat("Downregulated pathways:", nrow(ego_D7_D3_down_df), "\n") 
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  ekegg_D7_D3_all <- enrichKEGG(
    gene = entrez_D7_D3_all,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
  
  ekegg_D7_D3_up <- enrichKEGG(
    gene = entrez_D7_D3_up,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
  
  ekegg_D7_D3_down <- enrichKEGG(
    gene = entrez_D7_D3_down,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
  
  ekegg_D7_D3_all_df <- as.data.frame(ekegg_D7_D3_all)
  ekegg_D7_D3_up_df <- as.data.frame(ekegg_D7_D3_up)
  ekegg_D7_D3_down_df <- as.data.frame(ekegg_D7_D3_down)
  
  write.csv(
    ekegg_D7_D3_all_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_KEGG_all_DEGs.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ekegg_D7_D3_up_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_KEGG_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    ekegg_D7_D3_down_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_KEGG_downregulated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    ekegg_D7_D3_all,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_KEGG_all_DEGs.rds")
  )
  
  saveRDS(
    ekegg_D7_D3_up,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_KEGG_upregulated.rds")
  )
  
  saveRDS(
    ekegg_D7_D3_down,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_KEGG_downregulated.rds")
  )
  
  cat("D7 vs D3 KEGG enrichment completed and saved successfully.\n")
  cat("All DEGs pathways:", nrow(ekegg_D7_D3_all_df), "\n")
  cat("Upregulated pathways:", nrow(ekegg_D7_D3_up_df), "\n")
  cat("Downregulated pathways:", nrow(ekegg_D7_D3_down_df), "\n")  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  reactome_D7_D3_all <- enrichPathway(
    gene = entrez_D7_D3_all,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  reactome_D7_D3_up <- enrichPathway(
    gene = entrez_D7_D3_up,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  reactome_D7_D3_down <- enrichPathway(
    gene = entrez_D7_D3_down,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  reactome_D7_D3_all_df <- as.data.frame(reactome_D7_D3_all)
  reactome_D7_D3_up_df <- as.data.frame(reactome_D7_D3_up)
  reactome_D7_D3_down_df <- as.data.frame(reactome_D7_D3_down)
  
  write.csv(
    reactome_D7_D3_all_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_Reactome_all_DEGs.csv"),
    row.names = FALSE
  )
  
  write.csv(
    reactome_D7_D3_up_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_Reactome_upregulated.csv"),
    row.names = FALSE
  )
  
  write.csv(
    reactome_D7_D3_down_df,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_Reactome_downregulated.csv"),
    row.names = FALSE
  )
  
  saveRDS(
    reactome_D7_D3_all,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_Reactome_all_DEGs.rds")
  )
  
  saveRDS(
    reactome_D7_D3_up,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_Reactome_upregulated.rds")
  )
  
  saveRDS(
    reactome_D7_D3_down,
    file.path(output_dir, "COVID_Longitudinal_D7_vs_D3_Reactome_downregulated.rds")
  )
  
  cat("D7 vs D3 Reactome enrichment completed and saved successfully.\n")
  cat("All DEGs pathways:", nrow(reactome_D7_D3_all_df), "\n")
  cat("Upregulated pathways:", nrow(reactome_D7_D3_up_df), "\n")
  cat("Downregulated pathways:", nrow(reactome_D7_D3_down_df), "\n") 
  
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  prepare_pathway_summary <- function(df, database, comparison) {
    if (nrow(df) == 0) {
      return(data.frame())
    }
    
    result <- data.frame(
      Database = database,
      Comparison = comparison,
      Pathway = df$Description,
      GeneRatio = df$GeneRatio,
      BgRatio = df$BgRatio,
      GeneCount = df$Count,
      PValue = df$pvalue,
      AdjustedPValue = df$p.adjust,
      QValue = df$qvalue,
      stringsAsFactors = FALSE
    )
    
    result
  }
  
  pathway_D3_D0_GO <- prepare_pathway_summary(
    ego_D3_D0_all_df,
    "GO-BP",
    "D3_vs_D0"
  )
  
  pathway_D3_D0_KEGG <- prepare_pathway_summary(
    ekegg_D3_D0_all_df,
    "KEGG",
    "D3_vs_D0"
  )
  
  pathway_D3_D0_Reactome <- prepare_pathway_summary(
    reactome_D3_D0_all_df,
    "Reactome",
    "D3_vs_D0"
  )
  
  pathway_D7_D0_GO <- prepare_pathway_summary(
    ego_D7_D0_all_df,
    "GO-BP",
    "D7_vs_D0"
  )
  
  pathway_D7_D0_KEGG <- prepare_pathway_summary(
    ekegg_D7_D0_all_df,
    "KEGG",
    "D7_vs_D0"
  )
  
  pathway_D7_D0_Reactome <- prepare_pathway_summary(
    reactome_D7_D0_all_df,
    "Reactome",
    "D7_vs_D0"
  )
  
  pathway_D7_D3_GO <- prepare_pathway_summary(
    ego_D7_D3_all_df,
    "GO-BP",
    "D7_vs_D3"
  )
  
  pathway_D7_D3_KEGG <- prepare_pathway_summary(
    ekegg_D7_D3_all_df,
    "KEGG",
    "D7_vs_D3"
  )
  
  pathway_D7_D3_Reactome <- prepare_pathway_summary(
    reactome_D7_D3_all_df,
    "Reactome",
    "D7_vs_D3"
  )
  
  longitudinal_pathway_summary <- rbind(
    pathway_D3_D0_GO,
    pathway_D3_D0_KEGG,
    pathway_D3_D0_Reactome,
    pathway_D7_D0_GO,
    pathway_D7_D0_KEGG,
    pathway_D7_D0_Reactome,
    pathway_D7_D3_GO,
    pathway_D7_D3_KEGG,
    pathway_D7_D3_Reactome
  )
  
  longitudinal_pathway_summary <- longitudinal_pathway_summary[
    order(
      longitudinal_pathway_summary$Comparison,
      longitudinal_pathway_summary$Database,
      longitudinal_pathway_summary$AdjustedPValue
    ),
  ]
  
  rownames(longitudinal_pathway_summary) <- NULL
  
  write.csv(
    longitudinal_pathway_summary,
    file.path(
      output_dir,
      "COVID_Longitudinal_Pathway_Summary_All_DEGs.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    longitudinal_pathway_summary,
    file.path(
      output_dir,
      "COVID_Longitudinal_Pathway_Summary_All_DEGs.rds"
    )
  )
  
  cat("Longitudinal pathway summary created and saved successfully.\n")
  cat("Total pathway records:", nrow(longitudinal_pathway_summary), "\n")
  cat("Unique pathways:", length(unique(longitudinal_pathway_summary$Pathway)), "\n")
  
  cat("\nPathway records by comparison and database:\n")
  print(
    as.data.frame(
      table(
        longitudinal_pathway_summary$Comparison,
        longitudinal_pathway_summary$Database
      )
    )
  )
  
  cat("\nTop 10 pathways overall:\n")
  print(
    head(
      longitudinal_pathway_summary[
        order(longitudinal_pathway_summary$AdjustedPValue),
      ],
      10
    )
  )
  
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal"
  
  pathway_recurrence <- aggregate(
    Comparison ~ Database + Pathway,
    data = longitudinal_pathway_summary,
    FUN = function(x) paste(sort(unique(x)), collapse = "; ")
  )
  
  pathway_recurrence$N_Comparisons <- sapply(
    strsplit(pathway_recurrence$Comparison, "; "),
    length
  )
  
  pathway_recurrence$Recurrence <- ifelse(
    pathway_recurrence$N_Comparisons == 3,
    "All_3_Comparisons",
    ifelse(
      pathway_recurrence$N_Comparisons == 2,
      "2_Comparisons",
      "1_Comparison"
    )
  )
  
  pathway_recurrence <- pathway_recurrence[
    order(
      pathway_recurrence$Database,
      -pathway_recurrence$N_Comparisons,
      pathway_recurrence$Pathway
    ),
  ]
  
  rownames(pathway_recurrence) <- NULL
  
  write.csv(
    pathway_recurrence,
    file.path(
      output_dir,
      "COVID_Longitudinal_Pathway_Recurrence.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    pathway_recurrence,
    file.path(
      output_dir,
      "COVID_Longitudinal_Pathway_Recurrence.rds"
    )
  )
  
  cat("Pathway recurrence analysis completed and saved successfully.\n")
  cat("Unique pathway-database combinations:", nrow(pathway_recurrence), "\n")
  
  cat("\nRecurrence summary:\n")
  print(
    table(
      pathway_recurrence$Database,
      pathway_recurrence$Recurrence
    )
  )
  
  cat("\nPathways present in all 3 comparisons:\n")
  
  all_three <- pathway_recurrence[
    pathway_recurrence$N_Comparisons == 3,
  ]
  
  print(
    all_three[, c(
      "Database",
      "Pathway",
      "Comparison",
      "N_Comparisons"
    )]
  )  

  proteomics_file <- olink_files[1]
  
  proteomics_raw <- readxl::read_excel(proteomics_file)
  
  cat("Proteomics file loaded successfully.\n")
  cat("Rows:", nrow(proteomics_raw), "\n")
  cat("Columns:", ncol(proteomics_raw), "\n")
  
  cat("\nColumn names:\n")
  print(colnames(proteomics_raw))
  
  cat("\nFirst 6 rows:\n")
  print(head(proteomics_raw))
  
  cat("\nData types:\n")
  print(sapply(proteomics_raw, class)) 
  
  cat("=== Proteomics structural QC ===\n")
  
  protein_cols <- setdiff(colnames(proteomics_raw), c("Public ID", "Day"))
  
  cat("Number of protein assays:", length(protein_cols), "\n")
  cat("Number of samples:", nrow(proteomics_raw), "\n")
  
  cat("\nPublic ID unique:", length(unique(proteomics_raw$`Public ID`)) == nrow(proteomics_raw), "\n")
  cat("Missing Public IDs:", sum(is.na(proteomics_raw$`Public ID`) | proteomics_raw$`Public ID` == ""), "\n")
  
  cat("\nDay distribution:\n")
  print(table(proteomics_raw$Day, useNA = "ifany"))
  
  protein_classes <- sapply(proteomics_raw[protein_cols], class)
  
  cat("\nProtein columns by class:\n")
  print(table(protein_classes))
  
  non_numeric_proteins <- names(protein_classes)[protein_classes != "numeric"]
  
  cat("\nNon-numeric protein assays:\n")
  print(non_numeric_proteins)
  
  if (length(non_numeric_proteins) > 0) {
    cat("\nValues in non-numeric assays:\n")
    for (p in non_numeric_proteins) {
      cat("\n---", p, "---\n")
      print(head(unique(proteomics_raw[[p]]), 20))
    }
  }
  
  missing_counts <- sapply(proteomics_raw[protein_cols], function(x) {
    sum(is.na(x))
  })
  
  cat("\nMissing-value summary across protein assays:\n")
  cat("Assays with zero missing:", sum(missing_counts == 0), "\n")
  cat("Assays with any missing:", sum(missing_counts > 0), "\n")
  cat("Maximum missing per assay:", max(missing_counts), "\n")
  
  cat("\nSample-level missingness:\n")
  sample_missing <- rowSums(is.na(proteomics_raw[protein_cols]))
  print(summary(sample_missing))
  
  cat("\nConstant assays:\n")
  constant_assays <- protein_cols[sapply(proteomics_raw[protein_cols], function(x) {
    length(unique(x[!is.na(x)])) <= 1
  })]
  
  cat("Number of constant assays:", length(constant_assays), "\n")
  if (length(constant_assays) > 0) {
    print(constant_assays)
  }
  cat("=== Convert Olink assays to numeric ===\n")
  
  protein_numeric <- as.data.frame(
    lapply(proteomics_raw[protein_cols], function(x) {
      suppressWarnings(as.numeric(x))
    })
  )
  
  rownames(protein_numeric) <- proteomics_raw$`Public ID`
  
  conversion_na <- sapply(protein_numeric, function(x) {
    sum(is.na(x))
  })
  
  cat("Number of assays:", ncol(protein_numeric), "\n")
  cat("Number of samples:", nrow(protein_numeric), "\n")
  cat("Assays with conversion-generated NA:",
      sum(conversion_na > 0), "\n")
  cat("Maximum conversion-generated NA:",
      max(conversion_na), "\n")
  
  if (any(conversion_na > 0)) {
    cat("\nAssays affected by numeric conversion:\n")
    print(
      data.frame(
        Assay = names(conversion_na)[conversion_na > 0],
        N_NA = conversion_na[conversion_na > 0]
      )
    )
  }
  
  cat("\nNumeric structure check:\n")
  print(table(sapply(protein_numeric, class)))
  
  cat("\nOverall numeric value summary:\n")
  print(summary(as.vector(as.matrix(protein_numeric))))
  
  cat("\nQuantiles:\n")
  print(
    quantile(
      as.matrix(protein_numeric),
      probs = c(0, 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99, 1),
      na.rm = TRUE
    )
  )
  
  cat("=== Inspect non-numeric Olink values ===\n")
  
  problem_assays <- names(conversion_na)[conversion_na > 0]
  
  for (p in problem_assays) {
    original_values <- proteomics_raw[[p]]
    bad_idx <- which(
      is.na(protein_numeric[[p]]) &
        !is.na(original_values)
    )
    
    cat("\n---", p, "---\n")
    cat("Number of problematic values:", length(bad_idx), "\n")
    
    print(
      data.frame(
        Public_ID = proteomics_raw$`Public ID`[bad_idx],
        Day = proteomics_raw$Day[bad_idx],
        Original_Value = original_values[bad_idx],
        stringsAsFactors = FALSE
      )
    )
  }  
  
  cat("=== Final Proteomics Missingness QC ===\n")
  
  assay_missing_pct <- colSums(is.na(protein_numeric)) /
    nrow(protein_numeric) * 100
  
  sample_missing_pct <- rowSums(is.na(protein_numeric)) /
    ncol(protein_numeric) * 100
  
  proteomics_assay_qc <- data.frame(
    Assay = colnames(protein_numeric),
    Missing_N = colSums(is.na(protein_numeric)),
    Missing_Percent = assay_missing_pct,
    Mean = colMeans(protein_numeric, na.rm = TRUE),
    SD = apply(protein_numeric, 2, sd, na.rm = TRUE),
    Median = apply(protein_numeric, 2, median, na.rm = TRUE),
    Min = apply(protein_numeric, 2, min, na.rm = TRUE),
    Max = apply(protein_numeric, 2, max, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  proteomics_sample_qc <- data.frame(
    Public_ID = proteomics_raw$`Public ID`,
    Day = proteomics_raw$Day,
    Missing_N = rowSums(is.na(protein_numeric)),
    Missing_Percent = sample_missing_pct,
    stringsAsFactors = FALSE
  )
  
  cat("Total protein measurements:",
      nrow(protein_numeric) * ncol(protein_numeric), "\n")
  
  cat("Total missing measurements:",
      sum(is.na(protein_numeric)), "\n")
  
  cat("Overall missing percentage:",
      round(mean(is.na(protein_numeric)) * 100, 6), "%\n")
  
  cat("\nAssays with missing values:\n")
  print(
    proteomics_assay_qc[
      proteomics_assay_qc$Missing_N > 0,
    ]
  )
  
  cat("\nSample missingness summary:\n")
  print(summary(proteomics_sample_qc$Missing_N))
  
  cat("\nSamples with any missing protein values:\n")
  print(
    proteomics_sample_qc[
      proteomics_sample_qc$Missing_N > 0,
    ]
  )
  
  cat("\nMaximum sample missing percentage:",
      max(proteomics_sample_qc$Missing_Percent), "%\n")
  
  cat("\nAssays passing >=95% completeness:",
      sum(proteomics_assay_qc$Missing_Percent <= 5), "/",
      nrow(proteomics_assay_qc), "\n")
  
  cat("\nSamples passing >=95% completeness:",
      sum(proteomics_sample_qc$Missing_Percent <= 5), "/",
      nrow(proteomics_sample_qc), "\n")  
  
  cat("=== Load Olink assay annotation ===\n")
  
  supplemental_file <- file.path(
    dirname(proteomics_file),
    "Supplemental-Table-2-Olink-Assays-NPX-values_v2.xlsx"
  )
  
  cat("File exists:", file.exists(supplemental_file), "\n")
  cat("File path:", supplemental_file, "\n")
  
  olink_annotation_raw <- readxl::read_excel(
    supplemental_file
  )
  
  cat("\nAnnotation dimensions:\n")
  cat(
    nrow(olink_annotation_raw),
    "rows x",
    ncol(olink_annotation_raw),
    "columns\n"
  )
  
  cat("\nColumn names:\n")
  print(colnames(olink_annotation_raw))
  
  cat("\nFirst rows:\n")
  print(
    head(
      olink_annotation_raw,
      10
    )
  ) 
  
  cat("=== Build Olink assay annotation ===\n")
  
  olink_annotation <- readxl::read_excel(
    supplemental_file,
    skip = 1
  )
  
  cat("Annotation dimensions:",
      nrow(olink_annotation), "rows x",
      ncol(olink_annotation), "columns\n")
  
  cat("\nAnnotation columns:\n")
  print(colnames(olink_annotation))
  
  cat("\nFirst 10 annotation rows:\n")
  print(head(olink_annotation, 10))
  
  cat("\nOlinkID unique:\n")
  cat(
    length(unique(olink_annotation$OlinkID)) ==
      nrow(olink_annotation),
    "\n"
  )
  
  cat("\nMissing OlinkIDs:\n")
  print(sum(is.na(olink_annotation$OlinkID) |
              olink_annotation$OlinkID == ""))
  
  cat("\nProteomics assays found in annotation:\n")
  
  assay_mapping_check <- data.frame(
    OlinkID = protein_cols,
    In_Annotation = protein_cols %in% olink_annotation$OlinkID,
    stringsAsFactors = FALSE
  )
  
  cat(
    "Mapped:",
    sum(assay_mapping_check$In_Annotation),
    "/",
    nrow(assay_mapping_check),
    "\n"
  )
  
  cat(
    "Unmapped:",
    sum(!assay_mapping_check$In_Annotation),
    "\n"
  )
  
  if (any(!assay_mapping_check$In_Annotation)) {
    cat("\nUnmapped OlinkIDs:\n")
    print(
      assay_mapping_check$OlinkID[
        !assay_mapping_check$In_Annotation
      ]
    )
  }
  
  cat("=== Build final Olink proteomics dataset ===\n")
  
  proteomics_annotated <- olink_annotation[
    match(protein_cols, olink_annotation$OlinkID),
  ]
  
  cat("Annotation rows matched:",
      nrow(proteomics_annotated), "\n")
  
  cat("Unique OlinkIDs:",
      length(unique(proteomics_annotated$OlinkID)), "\n")
  
  cat("Missing annotation rows:",
      sum(is.na(proteomics_annotated$OlinkID)), "\n")
  
  protein_matrix_final <- protein_numeric
  
  colnames(protein_matrix_final) <- proteomics_annotated$Assay
  
  cat("\nProtein matrix dimensions:\n")
  cat(
    nrow(protein_matrix_final),
    "samples x",
    ncol(protein_matrix_final),
    "proteins\n"
  )
  
  cat("\nProtein names unique:\n")
  cat(
    length(unique(colnames(protein_matrix_final))) ==
      ncol(protein_matrix_final),
    "\n"
  )
  
  proteomics_sample_metadata <- data.frame(
    Public_ID = proteomics_raw$`Public ID`,
    Day = proteomics_raw$Day,
    stringsAsFactors = FALSE
  )
  
  cat("\nSample metadata dimensions:\n")
  print(dim(proteomics_sample_metadata))
  
  cat("\nFirst 10 mapped proteins:\n")
  print(
    proteomics_annotated[
      1:10,
      c("OlinkID", "Assay", "UniProt", "MissingFreq", "LOD", "Panel")
    ]
  )
  
  cat("\nFirst 10 protein matrix columns:\n")
  print(colnames(protein_matrix_final)[1:10])
  
  cat("\nFinal missing values:\n")
  cat(sum(is.na(protein_matrix_final)), "\n")
  
  cat("=== Check protein assay naming and duplicates ===\n")
  
  assay_names <- proteomics_annotated$Assay
  uniprot_ids <- proteomics_annotated$UniProt
  
  cat("Total assays:", length(assay_names), "\n")
  cat("Unique assay names:", length(unique(assay_names)), "\n")
  cat("Duplicated assay names:",
      sum(duplicated(assay_names)), "\n")
  
  cat("\nNumber of assay names occurring more than once:\n")
  assay_frequency <- sort(
    table(assay_names),
    decreasing = TRUE
  )
  
  print(
    assay_frequency[assay_frequency > 1]
  )
  
  cat("\nUniProt IDs:\n")
  cat("Unique UniProt IDs:",
      length(unique(uniprot_ids)),
      "\n")
  
  cat("Duplicated UniProt IDs:",
      sum(duplicated(uniprot_ids)),
      "\n")
  
  cat("\nMissing UniProt IDs:",
      sum(is.na(uniprot_ids) | uniprot_ids == ""),
      "\n")
  
  cat("\nExamples of duplicated assay names:\n")
  
  duplicated_assays <- names(
    assay_frequency[assay_frequency > 1]
  )
  
  if (length(duplicated_assays) > 0) {
    print(
      proteomics_annotated[
        proteomics_annotated$Assay %in% duplicated_assays,
        c("OlinkID", "Assay", "UniProt", "Panel")
      ]
    )
  } 
  
  cat("=== Save finalized proteomics objects ===\n")
  
  proteomics_results_dir <- file.path(
    "C:/Users/ibrah/OneDrive/Documents/Results",
    "Proteomics"
  )
  
  dir.create(
    proteomics_results_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  protein_matrix_assay_id <- protein_numeric
  rownames(protein_matrix_assay_id) <- proteomics_raw$`Public ID`
  
  proteomics_assay_annotation <- as.data.frame(
    proteomics_annotated,
    stringsAsFactors = FALSE
  )
  
  proteomics_final <- list(
    protein_matrix = protein_matrix_assay_id,
    sample_metadata = proteomics_sample_metadata,
    assay_annotation = proteomics_assay_annotation
  )
  
  saveRDS(
    proteomics_final,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Final_1429_Assays.rds"
    )
  )
  
  saveRDS(
    protein_matrix_assay_id,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Matrix_784x1429.rds"
    )
  )
  
  saveRDS(
    proteomics_sample_metadata,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Sample_Metadata.rds"
    )
  )
  
  write.csv(
    proteomics_assay_annotation,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Assay_Annotation.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    proteomics_assay_qc,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Assay_QC.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    proteomics_sample_qc,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Sample_QC.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    proteomics_assay_qc,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Assay_QC.rds"
    )
  )
  
  saveRDS(
    proteomics_sample_qc,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Sample_QC.rds"
    )
  )
  
  cat("\nSaved files:\n")
  print(
    list.files(
      proteomics_results_dir,
      pattern = "Olink_Proteomics",
      full.names = FALSE
    )
  )
  
  cat("\nFinal matrix dimensions:",
      nrow(protein_matrix_assay_id),
      "samples x",
      ncol(protein_matrix_assay_id),
      "assays\n")
  
  cat("Final missing values:",
      sum(is.na(protein_matrix_assay_id)),
      "\n")
  
  cat("Final assay annotation rows:",
      nrow(proteomics_assay_annotation),
      "\n")
  
  cat("=== Locate and load Clinical Metadata ===\n")
  
  clinical_files <- list.files(
    "C:/Users/ibrah/OneDrive/Desktop/GSE212041_MultiOmics_Project",
    pattern = "^Clinical Metadata\\.xlsx$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  cat("Clinical metadata files found:", length(clinical_files), "\n")
  print(clinical_files)
  
  clinical_file <- clinical_files[1]
  
  cat("\nFile exists:", file.exists(clinical_file), "\n")
  
  clinical_raw <- readxl::read_excel(
    clinical_file
  )
  
  cat("\nClinical metadata dimensions:\n")
  cat(
    nrow(clinical_raw),
    "rows x",
    ncol(clinical_raw),
    "columns\n"
  )
  
  cat("\nColumn names:\n")
  print(colnames(clinical_raw))
  
  cat("\nFirst 10 rows:\n")
  print(head(clinical_raw, 10))
  cat("=== Match proteomics samples to clinical metadata ===\n")
  
  clinical_patient_ids <- clinical_raw$`Public ID`
  
  proteomics_public_ids <- as.character(
    proteomics_raw$`Public ID`
  )
  
  clinical_patient_ids_chr <- as.character(
    clinical_patient_ids
  )
  
  proteomics_patient_ids <- sub(
    "_.*$",
    "",
    proteomics_public_ids
  )
  
  cat("Clinical patients:", length(clinical_patient_ids_chr), "\n")
  cat("Unique clinical patients:",
      length(unique(clinical_patient_ids_chr)), "\n")
  
  cat("\nProteomics samples:", length(proteomics_public_ids), "\n")
  cat("Unique proteomics sample IDs:",
      length(unique(proteomics_public_ids)), "\n")
  
  cat("\nProteomics patient IDs extracted:",
      length(unique(proteomics_patient_ids)), "\n")
  
  cat("\nProteomics patients found in clinical metadata:\n")
  
  patient_match <- proteomics_patient_ids %in%
    clinical_patient_ids_chr
  
  cat(
    "Matched samples:",
    sum(patient_match),
    "/",
    length(patient_match),
    "\n"
  )
  
  cat(
    "Unmatched samples:",
    sum(!patient_match),
    "\n"
  )
  
  cat(
    "Matched unique patients:",
    length(unique(
      proteomics_patient_ids[patient_match]
    )),
    "\n"
  )
  
  cat("\nUnmatched proteomics patient IDs:\n")
  print(
    sort(
      unique(
        proteomics_patient_ids[!patient_match]
      )
    )
  )
  
  cat("\nClinical COVID distribution:\n")
  print(
    table(
      clinical_raw$COVID,
      useNA = "ifany"
    )
  )
  
  cat("\nClinical draw availability:\n")
  print(
    sapply(
      clinical_raw[
        ,
        c(
          "D0_draw",
          "D3_draw",
          "D7_draw",
          "DE_draw"
        )
      ],
      function(x) table(x, useNA = "ifany")
    )
  )
  
  cat("\nProteomics sample Day distribution:\n")
  print(
    table(
      proteomics_raw$Day,
      useNA = "ifany"
    )
  )
  
  cat("\nProteomics samples per patient:\n")
  samples_per_patient <- table(
    proteomics_patient_ids
  )
  
  print(
    summary(
      as.numeric(samples_per_patient)
    )
  )
  
  cat("\nPatients with repeated proteomics samples:\n")
  cat(
    sum(samples_per_patient > 1),
    "\n"
  )
  
  cat("\nMaximum proteomics samples per patient:\n")
  cat(
    max(samples_per_patient),
    "\n"
  ) 
  
  cat("=== Clinical–Proteomics patient and timepoint consistency QC ===\n")
  
  clinical_ids <- as.character(clinical_raw$`Public ID`)
  proteomics_ids <- unique(proteomics_patient_ids)
  
  clinical_only <- setdiff(
    clinical_ids,
    proteomics_ids
  )
  
  proteomics_only <- setdiff(
    proteomics_ids,
    clinical_ids
  )
  
  cat("Clinical-only patients:",
      length(clinical_only), "\n")
  print(clinical_only)
  
  cat("\nProteomics-only patients:",
      length(proteomics_only), "\n")
  print(proteomics_only)
  
  cat("\nProteomics sample counts by patient and Day:\n")
  
  patient_day_counts <- aggregate(
    proteomics_public_ids,
    by = list(
      Patient_ID = proteomics_patient_ids,
      Day = as.character(proteomics_raw$Day)
    ),
    FUN = length
  )
  
  colnames(patient_day_counts)[3] <- "N_Samples"
  
  print(
    patient_day_counts[
      patient_day_counts$N_Samples > 1,
    ]
  )
  
  cat("\nChecking clinical draw flags against proteomics Day...\n")
  
  clinical_lookup <- clinical_raw
  clinical_lookup$Patient_ID <- as.character(
    clinical_lookup$`Public ID`
  )
  
  proteomics_clinical <- merge(
    data.frame(
      Public_ID = proteomics_public_ids,
      Patient_ID = proteomics_patient_ids,
      Day = as.character(proteomics_raw$Day),
      stringsAsFactors = FALSE
    ),
    clinical_lookup,
    by = "Patient_ID",
    all.x = TRUE,
    sort = FALSE
  )
  
  proteomics_clinical$Expected_Draw <- NA_integer_
  
  proteomics_clinical$Expected_Draw[
    proteomics_clinical$Day == "0"
  ] <- proteomics_clinical$D0_draw[
    proteomics_clinical$Day == "0"
  ]
  
  proteomics_clinical$Expected_Draw[
    proteomics_clinical$Day == "3"
  ] <- proteomics_clinical$D3_draw[
    proteomics_clinical$Day == "3"
  ]
  
  proteomics_clinical$Expected_Draw[
    proteomics_clinical$Day == "7"
  ] <- proteomics_clinical$D7_draw[
    proteomics_clinical$Day == "7"
  ]
  
  proteomics_clinical$Expected_Draw[
    proteomics_clinical$Day == "E"
  ] <- proteomics_clinical$DE_draw[
    proteomics_clinical$Day == "E"
  ]
  
  cat("\nExpected draw flag distribution:\n")
  print(
    table(
      proteomics_clinical$Expected_Draw,
      useNA = "ifany"
    )
  )
  
  draw_mismatch <- proteomics_clinical[
    !is.na(proteomics_clinical$Expected_Draw) &
      proteomics_clinical$Expected_Draw != 1,
  ]
  
  cat("\nProteomics samples whose clinical draw flag is not 1:\n")
  cat(nrow(draw_mismatch), "\n")
  
  if (nrow(draw_mismatch) > 0) {
    print(
      draw_mismatch[
        ,
        c(
          "Public_ID",
          "Patient_ID",
          "Day",
          "Expected_Draw"
        )
      ]
    )
  }
  
  cat("\nClinical variables missing after merge:\n")
  print(
    colSums(
      is.na(
        proteomics_clinical[
          ,
          c(
            "COVID",
            "Acuity max",
            "D0_draw",
            "D3_draw",
            "D7_draw",
            "DE_draw"
          )
        ]
      )
    )
  )  
  
  cat("=== Proteomics–Clinical master metadata QC ===\n")
  
  master_metadata <- proteomics_clinical
  
  master_metadata$Timepoint <- dplyr::case_when(
    master_metadata$Day == "0" ~ "D0",
    master_metadata$Day == "3" ~ "D3",
    master_metadata$Day == "7" ~ "D7",
    master_metadata$Day == "E" ~ "DE",
    TRUE ~ NA_character_
  )
  
  master_metadata$Timepoint <- factor(
    master_metadata$Timepoint,
    levels = c("D0", "D3", "D7", "DE")
  )
  
  master_metadata$COVID <- as.integer(master_metadata$COVID)
  master_metadata$Acuity_max <- as.numeric(master_metadata$`Acuity max`)
  
  cat("\n=== Samples by COVID and Timepoint ===\n")
  print(
    table(
      master_metadata$COVID,
      master_metadata$Timepoint,
      useNA = "ifany"
    )
  )
  
  cat("\n=== Samples by Timepoint ===\n")
  print(table(master_metadata$Timepoint, useNA = "ifany"))
  
  cat("\n=== COVID+ samples by Timepoint ===\n")
  print(
    table(
      master_metadata$Timepoint[
        master_metadata$COVID == 1
      ],
      useNA = "ifany"
    )
  )
  
  cat("\n=== The 13 draw-flag discrepancies ===\n")
  
  print(
    master_metadata[
      !is.na(master_metadata$Expected_Draw) &
        master_metadata$Expected_Draw != 1,
      c(
        "Public_ID",
        "Patient_ID",
        "Day",
        "Timepoint",
        "COVID",
        "Acuity_max",
        "D0_draw",
        "D3_draw",
        "D7_draw",
        "DE_draw"
      )
    ]
  )
  
  cat("\n=== Duplicate patient/timepoint records ===\n")
  
  duplicate_keys <- duplicated(
    master_metadata[
      ,
      c("Patient_ID", "Timepoint")
    ]
  ) |
    duplicated(
      master_metadata[
        ,
        c("Patient_ID", "Timepoint")
      ],
      fromLast = TRUE
    )
  
  duplicates_detail <- master_metadata[
    duplicate_keys,
    c(
      "Public_ID",
      "Patient_ID",
      "Day",
      "Timepoint",
      "COVID",
      "Acuity_max"
    )
  ]
  
  print(duplicates_detail)
  
  cat("\nNumber of duplicate patient-timepoint samples:",
      nrow(duplicates_detail), "\n")
  
  cat("\n=== All Public IDs for duplicate patients ===\n")
  
  duplicate_patients <- unique(
    duplicates_detail$Patient_ID
  )
  
  for (p in duplicate_patients) {
    
    cat("\nPatient", p, ":\n")
    
    print(
      master_metadata[
        master_metadata$Patient_ID == p,
        c(
          "Public_ID",
          "Patient_ID",
          "Day",
          "Timepoint",
          "COVID",
          "Acuity_max"
        )
      ]
    )
  }
  
  cat("\n=== Patient-level COVID consistency ===\n")
  
  covid_consistency <- aggregate(
    COVID ~ Patient_ID,
    data = master_metadata,
    FUN = function(x) length(unique(x))
  )
  
  print(
    table(
      covid_consistency$COVID
    )
  )
  
  cat("\n=== Patient-level Acuity consistency ===\n")
  
  acuity_consistency <- aggregate(
    Acuity_max ~ Patient_ID,
    data = master_metadata,
    FUN = function(x) length(unique(x))
  )
  
  print(
    table(
      acuity_consistency$Acuity_max
    )
  )
  
  cat("\n=== Master metadata dimensions ===\n")
  cat(
    "Rows:",
    nrow(master_metadata),
    "\n"
  )
  
  cat(
    "Columns:",
    ncol(master_metadata),
    "\n"
  )
  
  cat(
    "Unique patients:",
    length(unique(master_metadata$Patient_ID)),
    "\n"
  ) 
  
  cat("=== Duplicate Proteomics Sample QC ===\n")
  
  duplicate_sample_ids <- duplicates_detail$Public_ID
  
  duplicate_values <- protein_matrix_assay_id[
    duplicate_sample_ids,
    ,
    drop = FALSE
  ]
  
  cat("Duplicate samples found:",
      nrow(duplicate_values), "\n")
  
  cat("Proteins available:",
      ncol(duplicate_values), "\n\n")
  
  duplicate_pairs <- list(
    Patient_59_D3 = c("59_D3", "59_D3.1"),
    Patient_321_D7 = c("321_D7.1", "321_D7.2"),
    Patient_344_DE = c("344_DE.1", "344_DE.2")
  )
  
  duplicate_qc <- data.frame(
    Patient = character(),
    Timepoint = character(),
    Sample_1 = character(),
    Sample_2 = character(),
    Pearson_Correlation = numeric(),
    Spearman_Correlation = numeric(),
    Mean_Absolute_Difference = numeric(),
    Median_Absolute_Difference = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (pair_name in names(duplicate_pairs)) {
    
    ids <- duplicate_pairs[[pair_name]]
    
    x <- as.numeric(
      protein_matrix_assay_id[ids[1], ]
    )
    
    y <- as.numeric(
      protein_matrix_assay_id[ids[2], ]
    )
    
    valid <- is.finite(x) & is.finite(y)
    
    pearson_cor <- cor(
      x[valid],
      y[valid],
      method = "pearson"
    )
    
    spearman_cor <- cor(
      x[valid],
      y[valid],
      method = "spearman"
    )
    
    abs_diff <- abs(
      x[valid] - y[valid]
    )
    
    patient_id <- sub(
      "_.*$",
      "",
      ids[1]
    )
    
    timepoint <- sub(
      "^[^.]+_",
      "",
      ids[1]
    )
    
    duplicate_qc <- rbind(
      duplicate_qc,
      data.frame(
        Patient = patient_id,
        Timepoint = timepoint,
        Sample_1 = ids[1],
        Sample_2 = ids[2],
        Pearson_Correlation = pearson_cor,
        Spearman_Correlation = spearman_cor,
        Mean_Absolute_Difference = mean(abs_diff),
        Median_Absolute_Difference = median(abs_diff),
        stringsAsFactors = FALSE
      )
    )
  }
  
  print(duplicate_qc)
  
  cat("\n=== Duplicate QC interpretation ranges ===\n")
  cat("Pearson/Spearman close to 1 = highly similar samples\n")
  cat("Lower correlation + large absolute differences = potentially distinct samples\n")
  
  cat("\n=== Missing values in duplicate pairs ===\n")
  
  for (pair_name in names(duplicate_pairs)) {
    
    ids <- duplicate_pairs[[pair_name]]
    
    cat("\n", pair_name, "\n", sep = "")
    
    for (id in ids) {
      
      cat(
        id,
        ":",
        sum(
          !is.finite(
            protein_matrix_assay_id[id, ]
          )
        ),
        "missing/non-finite values\n"
      )
    }
  } 
  
  cat("=== Detailed duplicate protein-level QC ===\n")
  
  duplicate_detail_qc <- data.frame(
    Patient = character(),
    Timepoint = character(),
    Sample_1 = character(),
    Sample_2 = character(),
    N_Valid_Proteins = integer(),
    Median_Absolute_Difference = numeric(),
    Mean_Absolute_Difference = numeric(),
    P95_Absolute_Difference = numeric(),
    Max_Absolute_Difference = numeric(),
    Correlation = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (pair_name in names(duplicate_pairs)) {
    
    ids <- duplicate_pairs[[pair_name]]
    
    x <- as.numeric(
      unlist(
        protein_matrix_assay_id[ids[1], ],
        use.names = FALSE
      )
    )
    
    y <- as.numeric(
      unlist(
        protein_matrix_assay_id[ids[2], ],
        use.names = FALSE
      )
    )
    
    valid <- is.finite(x) & is.finite(y)
    
    abs_diff <- abs(x[valid] - y[valid])
    
    patient_id <- sub(
      "_.*$",
      "",
      ids[1]
    )
    
    timepoint <- sub(
      "^[^.]+_",
      "",
      ids[1]
    )
    
    duplicate_detail_qc <- rbind(
      duplicate_detail_qc,
      data.frame(
        Patient = patient_id,
        Timepoint = timepoint,
        Sample_1 = ids[1],
        Sample_2 = ids[2],
        N_Valid_Proteins = sum(valid),
        Median_Absolute_Difference = median(abs_diff),
        Mean_Absolute_Difference = mean(abs_diff),
        P95_Absolute_Difference = as.numeric(
          quantile(abs_diff, 0.95)
        ),
        Max_Absolute_Difference = max(abs_diff),
        Correlation = cor(
          x[valid],
          y[valid],
          method = "pearson"
        ),
        stringsAsFactors = FALSE
      )
    )
  }
  
  print(duplicate_detail_qc)
  
  cat("\n=== Missing/non-finite values in duplicate samples ===\n")
  
  for (pair_name in names(duplicate_pairs)) {
    
    ids <- duplicate_pairs[[pair_name]]
    
    for (id in ids) {
      
      values <- as.numeric(
        unlist(
          protein_matrix_assay_id[id, ],
          use.names = FALSE
        )
      )
      
      cat(
        id,
        ":",
        sum(!is.finite(values)),
        "missing/non-finite values\n"
      )
    }
  }
  
  cat("\n=== Overall duplicate QC conclusion ===\n")
  
  cat(
    "All duplicate pairs have Pearson correlation >= ",
    round(min(duplicate_detail_qc$Correlation), 3),
    "\n",
    sep = ""
  )
  
  cat(
    "Highest mean absolute difference: ",
    round(
      max(duplicate_detail_qc$Mean_Absolute_Difference),
      3
    ),
    "\n",
    sep = ""
  )
  cat("=== Building final proteomics master dataset ===\n")
  
  protein_matrix_raw <- protein_matrix_assay_id
  
  master_metadata_final <- master_metadata
  
  master_metadata_final$Patient_ID <- as.character(
    master_metadata_final$Patient_ID
  )
  
  master_metadata_final$Public_ID <- as.character(
    master_metadata_final$Public_ID
  )
  
  master_metadata_final$Timepoint <- factor(
    master_metadata_final$Timepoint,
    levels = c("D0", "D3", "D7", "DE")
  )
  
  master_metadata_final$Duplicate_Patient_Timepoint <- duplicate_keys
  
  master_metadata_final$Draw_Metadata_Discrepancy <- (
    !is.na(master_metadata_final$Expected_Draw) &
      master_metadata_final$Expected_Draw != 1
  )
  
  cat("\nRaw matrix:\n")
  cat(
    "Samples:",
    nrow(protein_matrix_raw),
    "\n"
  )
  cat(
    "Assays:",
    ncol(protein_matrix_raw),
    "\n"
  )
  
  cat("\nCreating one-row-per-patient-timepoint matrix...\n")
  
  unique_keys <- paste(
    master_metadata_final$Patient_ID,
    master_metadata_final$Timepoint,
    sep = "_"
  )
  
  unique_key_order <- unique(unique_keys)
  
  protein_matrix_collapsed <- matrix(
    NA_real_,
    nrow = length(unique_key_order),
    ncol = ncol(protein_matrix_raw)
  )
  
  colnames(protein_matrix_collapsed) <- colnames(
    protein_matrix_raw
  )
  
  rownames(protein_matrix_collapsed) <- unique_key_order
  
  for (i in seq_along(unique_key_order)) {
    
    key <- unique_key_order[i]
    
    sample_ids <- master_metadata_final$Public_ID[
      unique_keys == key
    ]
    
    values <- protein_matrix_raw[
      sample_ids,
      ,
      drop = FALSE
    ]
    
    if (nrow(values) == 1) {
      
      protein_matrix_collapsed[i, ] <- as.numeric(
        values[1, ]
      )
      
    } else {
      
      protein_matrix_collapsed[i, ] <- colMeans(
        values,
        na.rm = TRUE
      )
    }
  }
  
  protein_matrix_collapsed <- as.data.frame(
    protein_matrix_collapsed,
    check.names = FALSE
  )
  
  cat("\nCollapsed matrix dimensions:\n")
  cat(
    "Patient-timepoints:",
    nrow(protein_matrix_collapsed),
    "\n"
  )
  cat(
    "Assays:",
    ncol(protein_matrix_collapsed),
    "\n"
  )
  
  cat("\nChecking for duplicated patient-timepoints:\n")
  
  collapsed_keys <- rownames(
    protein_matrix_collapsed
  )
  
  cat(
    "Unique patient-timepoints:",
    length(unique(collapsed_keys)),
    "\n"
  )
  
  cat(
    "Duplicated patient-timepoints:",
    sum(duplicated(collapsed_keys)),
    "\n"
  )
  
  cat("\nMissing values after collapsing:\n")
  
  collapsed_missing <- sum(
    is.na(
      as.matrix(
        protein_matrix_collapsed
      )
    )
  )
  
  cat(
    "Total missing values:",
    collapsed_missing,
    "\n"
  )
  
  cat("\nPatients represented:\n")
  cat(
    length(
      unique(
        master_metadata_final$Patient_ID
      )
    ),
    "\n"
  )
  
  cat("\nTimepoint distribution after collapsing:\n")
  print(
    table(
      master_metadata_final$Timepoint
    )
  )
  
  cat("\n=== Saving final proteomics master objects ===\n")
  
  saveRDS(
    master_metadata_final,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Clinical_Master_Metadata_784_samples.rds"
    )
  )
  
  write.csv(
    master_metadata_final,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Clinical_Master_Metadata_784_samples.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    protein_matrix_collapsed,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Matrix_Patient_Timepoint_Collapsed.rds"
    )
  )
  
  saveRDS(
    list(
      raw_matrix = protein_matrix_raw,
      collapsed_matrix = protein_matrix_collapsed,
      metadata = master_metadata_final,
      assay_annotation = proteomics_assay_annotation
    ),
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Final_Integrated_Master_Object.rds"
    )
  )
  
  cat("\nSaved successfully.\n")  
  cat("=== Proteomics PCA QC ===\n")
  
  protein_matrix_pca <- as.matrix(
    protein_matrix_collapsed
  )
  
  storage.mode(protein_matrix_pca) <- "numeric"
  
  cat("Matrix dimensions:\n")
  cat(
    "Samples:",
    nrow(protein_matrix_pca),
    "\n"
  )
  cat(
    "Assays:",
    ncol(protein_matrix_pca),
    "\n"
  )
  
  cat("\nMissing values:\n")
  cat(
    sum(is.na(protein_matrix_pca)),
    "\n"
  )
  
  protein_matrix_pca_imputed <- protein_matrix_pca
  
  for (j in seq_len(ncol(protein_matrix_pca_imputed))) {
    
    missing_j <- is.na(
      protein_matrix_pca_imputed[, j]
    )
    
    if (any(missing_j)) {
      
      protein_matrix_pca_imputed[
        missing_j,
        j
      ] <- median(
        protein_matrix_pca_imputed[
          !missing_j,
          j
        ],
        na.rm = TRUE
      )
    }
  }
  
  cat(
    "\nMissing values after median imputation:",
    sum(is.na(protein_matrix_pca_imputed)),
    "\n"
  )
  
  protein_variance <- apply(
    protein_matrix_pca_imputed,
    2,
    var
  )
  
  cat("\nProtein variance summary:\n")
  print(
    summary(protein_variance)
  )
  
  top_n <- min(
    500,
    length(protein_variance)
  )
  
  top_proteins <- names(
    sort(
      protein_variance,
      decreasing = TRUE
    )
  )[seq_len(top_n)]
  
  cat(
    "\nTop variable proteins used for PCA:",
    top_n,
    "\n"
  )
  
  pca_proteomics <- prcomp(
    protein_matrix_pca_imputed[
      ,
      top_proteins,
      drop = FALSE
    ],
    center = TRUE,
    scale. = TRUE
  )
  
  cat("\nPCA completed.\n")
  
  variance_explained <- (
    pca_proteomics$sdev^2
  ) / sum(
    pca_proteomics$sdev^2
  ) * 100
  
  cat("\nVariance explained:\n")
  print(
    round(
      variance_explained[1:10],
      3
    )
  )
  
  cat(
    "\nPC1 + PC2:",
    round(
      sum(variance_explained[1:2]),
      3
    ),
    "%\n"
  )
  
  pca_proteomics_scores <- as.data.frame(
    pca_proteomics$x[, 1:5]
  )
  
  pca_proteomics_scores$Patient_ID <- rownames(
    protein_matrix_collapsed
  )
  
  pca_proteomics_scores$Patient_ID <- sub(
    "_.*$",
    "",
    pca_proteomics_scores$Patient_ID
  )
  
  pca_proteomics_scores$Timepoint <- sub(
    "^[^_]+_",
    "",
    rownames(
      protein_matrix_collapsed
    )
  )
  
  metadata_lookup <- master_metadata_final[
    !duplicated(
      master_metadata_final$Public_ID
    ),
  ]
  
  metadata_lookup$Patient_ID <- as.character(
    metadata_lookup$Patient_ID
  )
  
  pca_proteomics_scores$COVID <- metadata_lookup$COVID[
    match(
      rownames(protein_matrix_collapsed),
      metadata_lookup$Public_ID
    )
  ]
  
  pca_proteomics_scores$Acuity_max <- metadata_lookup$Acuity_max[
    match(
      rownames(protein_matrix_collapsed),
      metadata_lookup$Public_ID
    )
  ]
  
  cat("\nPCA score dimensions:\n")
  print(
    dim(pca_proteomics_scores)
  )
  
  cat("\nTimepoint distribution in PCA scores:\n")
  print(
    table(
      pca_proteomics_scores$Timepoint,
      useNA = "ifany"
    )
  )
  
  cat("\nCOVID distribution in PCA scores:\n")
  print(
    table(
      pca_proteomics_scores$COVID,
      useNA = "ifany"
    )
  )
  
  cat("\nSaving PCA objects...\n")
  
  saveRDS(
    pca_proteomics,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_Object.rds"
    )
  )
  
  write.csv(
    pca_proteomics_scores,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_Scores.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    pca_proteomics_scores,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_Scores.rds"
    )
  )
  
  write.csv(
    data.frame(
      PC = seq_along(variance_explained),
      Variance_Explained_Percent = variance_explained
    ),
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_Variance_Explained.csv"
    ),
    row.names = FALSE
  )
  
  cat("\nPCA QC objects saved successfully.\n")
  
  cat("=== Correcting proteomics PCA metadata mapping ===\n")
  
  pca_key <- rownames(
    protein_matrix_collapsed
  )
  
  pca_patient <- sub(
    "_.*$",
    "",
    pca_key
  )
  
  pca_timepoint <- sub(
    "^[^_]+_",
    "",
    pca_key
  )
  
  pca_metadata <- master_metadata_final
  
  pca_metadata$Patient_ID <- as.character(
    pca_metadata$Patient_ID
  )
  
  pca_metadata$Timepoint <- as.character(
    pca_metadata$Timepoint
  )
  
  pca_metadata$key <- paste(
    pca_metadata$Patient_ID,
    pca_metadata$Timepoint,
    sep = "_"
  )
  
  pca_metadata_unique <- pca_metadata[
    !duplicated(pca_metadata$key),
  ]
  
  pca_proteomics_scores$Patient_ID <- pca_patient
  pca_proteomics_scores$Timepoint <- pca_timepoint
  
  pca_proteomics_scores$COVID <- pca_metadata_unique$COVID[
    match(
      pca_key,
      pca_metadata_unique$key
    )
  ]
  
  pca_proteomics_scores$Acuity_max <- pca_metadata_unique$Acuity_max[
    match(
      pca_key,
      pca_metadata_unique$key
    )
  ]
  
  pca_proteomics_scores$Public_ID <- pca_key
  
  cat("\n=== Validation ===\n")
  
  cat(
    "Total PCA records:",
    nrow(pca_proteomics_scores),
    "\n"
  )
  
  cat(
    "Missing Patient IDs:",
    sum(is.na(pca_proteomics_scores$Patient_ID)),
    "\n"
  )
  
  cat(
    "Missing Timepoints:",
    sum(is.na(pca_proteomics_scores$Timepoint)),
    "\n"
  )
  
  cat(
    "Missing COVID:",
    sum(is.na(pca_proteomics_scores$COVID)),
    "\n"
  )
  
  cat(
    "Missing Acuity:",
    sum(is.na(pca_proteomics_scores$Acuity_max)),
    "\n"
  )
  
  cat("\nTimepoint distribution:\n")
  print(
    table(
      pca_proteomics_scores$Timepoint,
      useNA = "ifany"
    )
  )
  
  cat("\nCOVID distribution:\n")
  print(
    table(
      pca_proteomics_scores$COVID,
      useNA = "ifany"
    )
  )
  
  cat("\nAcuity distribution:\n")
  print(
    table(
      pca_proteomics_scores$Acuity_max,
      useNA = "ifany"
    )
  )
  
  cat("\n=== Checking unique patient-timepoint keys ===\n")
  
  cat(
    "PCA keys:",
    length(unique(pca_key)),
    "\n"
  )
  
  cat(
    "Metadata keys:",
    length(unique(pca_metadata_unique$key)),
    "\n"
  )
  
  cat(
    "Unmatched PCA keys:",
    sum(
      !pca_key %in%
        pca_metadata_unique$key
    ),
    "\n"
  )
  
  cat("\nSaving corrected PCA scores...\n")
  
  saveRDS(
    pca_proteomics_scores,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_Scores.rds"
    )
  )
  
  write.csv(
    pca_proteomics_scores,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_Scores.csv"
    ),
    row.names = FALSE
  )
  
  cat("Corrected PCA scores saved successfully.\n") 

  
  
  cat("=== Proteomics PCA visualization and outlier screening ===\n")
  
  library(ggplot2)
  
  pca_plot_data <- pca_proteomics_scores
  
  pca_plot_data$Timepoint <- factor(
    pca_plot_data$Timepoint,
    levels = c("D0", "D3", "D7", "DE")
  )
  
  pca_plot_data$COVID_Group <- factor(
    pca_plot_data$COVID,
    levels = c(0, 1),
    labels = c("COVID-negative", "COVID-positive")
  )
  
  pca_plot_data$Acuity_max <- factor(
    pca_plot_data$Acuity_max,
    levels = c(1, 2, 3, 4, 5)
  )
  
  cat("\n=== PCA range ===\n")
  cat(
    "PC1:",
    round(min(pca_plot_data$PC1), 2),
    "to",
    round(max(pca_plot_data$PC1), 2),
    "\n"
  )
  
  cat(
    "PC2:",
    round(min(pca_plot_data$PC2), 2),
    "to",
    round(max(pca_plot_data$PC2), 2),
    "\n"
  )
  
  cat("\n=== Mean PC scores by Timepoint ===\n")
  print(
    aggregate(
      cbind(PC1, PC2) ~ Timepoint,
      data = pca_plot_data,
      FUN = mean
    )
  )
  
  cat("\n=== Mean PC scores by COVID group ===\n")
  print(
    aggregate(
      cbind(PC1, PC2) ~ COVID_Group,
      data = pca_plot_data,
      FUN = mean
    )
  )
  
  pca_plot_timepoint <- ggplot(
    pca_plot_data,
    aes(
      x = PC1,
      y = PC2,
      color = Timepoint
    )
  ) +
    geom_point(
      alpha = 0.65,
      size = 2
    ) +
    theme_classic() +
    labs(
      title = "Olink Proteomics PCA",
      subtitle = "Top 500 variable assays",
      x = paste0(
        "PC1 (",
        round(variance_explained[1], 2),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(variance_explained[2], 2),
        "%)"
      )
    )
  
  print(pca_plot_timepoint)
  
  ggsave(
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_PC1_PC2_by_Timepoint.png"
    ),
    pca_plot_timepoint,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  pca_plot_covid <- ggplot(
    pca_plot_data,
    aes(
      x = PC1,
      y = PC2,
      color = COVID_Group
    )
  ) +
    geom_point(
      alpha = 0.65,
      size = 2
    ) +
    theme_classic() +
    labs(
      title = "Olink Proteomics PCA",
      subtitle = "Colored by COVID status",
      x = paste0(
        "PC1 (",
        round(variance_explained[1], 2),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(variance_explained[2], 2),
        "%)"
      )
    )
  
  print(pca_plot_covid)
  
  ggsave(
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_PC1_PC2_by_COVID.png"
    ),
    pca_plot_covid,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  cat("\n=== Mahalanobis-style PCA distance screening ===\n")
  
  pca_coordinates <- pca_plot_data[
    ,
    c("PC1", "PC2", "PC3", "PC4", "PC5")
  ]
  
  pca_center <- colMeans(
    pca_coordinates
  )
  
  pca_cov <- cov(
    pca_coordinates
  )
  
  pca_distance <- mahalanobis(
    pca_coordinates,
    center = pca_center,
    cov = pca_cov
  )
  
  pca_plot_data$PCA_Mahalanobis_Distance <- pca_distance
  
  outlier_threshold <- qchisq(
    0.999,
    df = 5
  )
  
  pca_outliers <- pca_plot_data[
    pca_plot_data$PCA_Mahalanobis_Distance >
      outlier_threshold,
    ,
    drop = FALSE
  ]
  
  cat(
    "Outlier threshold (99.9%):",
    round(outlier_threshold, 3),
    "\n"
  )
  
  cat(
    "Potential PCA outliers:",
    nrow(pca_outliers),
    "\n"
  )
  
  if (nrow(pca_outliers) > 0) {
    print(
      pca_outliers[
        order(
          -pca_outliers$PCA_Mahalanobis_Distance
        ),
        c(
          "Public_ID",
          "Patient_ID",
          "Timepoint",
          "COVID_Group",
          "Acuity_max",
          "PCA_Mahalanobis_Distance"
        )
      ]
    )
  }
  
  cat("\nSaving PCA outlier screening...\n")
  
  write.csv(
    pca_plot_data,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_Scores_with_Outlier_Distance.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    pca_outliers,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_Potential_Outliers.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    pca_outliers,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_PCA_Potential_Outliers.rds"
    )
  )
  
  cat("\nPCA visualization and outlier screening saved successfully.\n")    
  
  cat("=== Proteomics sample-to-sample correlation QC ===\n")
  
  protein_qc_matrix <- as.matrix(protein_matrix_collapsed)
  
  storage.mode(protein_qc_matrix) <- "numeric"
  
  cat(
    "Matrix:",
    nrow(protein_qc_matrix),
    "samples x",
    ncol(protein_qc_matrix),
    "assays\n"
  )
  
  cat(
    "Missing values before imputation:",
    sum(is.na(protein_qc_matrix)),
    "\n"
  )
  
  protein_qc_variance <- apply(
    protein_qc_matrix,
    2,
    var,
    na.rm = TRUE
  )
  
  top_qc_n <- min(500, length(protein_qc_variance))
  
  top_qc_assays <- names(
    sort(
      protein_qc_variance,
      decreasing = TRUE
    )[seq_len(top_qc_n)]
  )
  
  protein_qc_top <- protein_qc_matrix[
    ,
    top_qc_assays,
    drop = FALSE
  ]
  
  for (j in seq_len(ncol(protein_qc_top))) {
    missing_idx <- is.na(protein_qc_top[, j])
    
    if (any(missing_idx)) {
      protein_qc_top[missing_idx, j] <- median(
        protein_qc_top[, j],
        na.rm = TRUE
      )
    }
  }
  
  cat(
    "Top variable assays used:",
    ncol(protein_qc_top),
    "\n"
  )
  
  cat(
    "Missing values after imputation:",
    sum(is.na(protein_qc_top)),
    "\n"
  )
  
  sample_correlation <- cor(
    t(protein_qc_top),
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  
  cat(
    "Correlation matrix:",
    nrow(sample_correlation),
    "x",
    ncol(sample_correlation),
    "\n"
  )
  
  median_sample_correlation <- apply(
    sample_correlation,
    1,
    function(x) {
      median(
        x[x < 0.999999],
        na.rm = TRUE
      )
    }
  )
  
  mean_sample_correlation <- apply(
    sample_correlation,
    1,
    function(x) {
      mean(
        x[x < 0.999999],
        na.rm = TRUE
      )
    }
  )
  
  sample_correlation_qc <- data.frame(
    Patient_Timepoint = rownames(protein_qc_top),
    Median_Correlation = median_sample_correlation,
    Mean_Correlation = mean_sample_correlation,
    stringsAsFactors = FALSE
  )
  
  sample_correlation_qc$PCA_Mahalanobis_Distance <-
    pca_plot_data$PCA_Mahalanobis_Distance[
      match(
        sample_correlation_qc$Patient_Timepoint,
        pca_plot_data$Public_ID
      )
    ]
  
  sample_correlation_qc <- sample_correlation_qc[
    order(
      sample_correlation_qc$Median_Correlation
    ),
    ,
    drop = FALSE
  ]
  
  cat("\n=== Correlation summary ===\n")
  
  print(
    summary(
      sample_correlation_qc$Median_Correlation
    )
  )
  
  cat("\n=== Lowest-correlation samples ===\n")
  
  print(
    head(
      sample_correlation_qc,
      20
    )
  )
  
  correlation_threshold <- quantile(
    sample_correlation_qc$Median_Correlation,
    0.01,
    na.rm = TRUE
  )
  
  correlation_flags <- sample_correlation_qc[
    sample_correlation_qc$Median_Correlation <
      correlation_threshold,
    ,
    drop = FALSE
  ]
  
  cat(
    "\n1st percentile correlation threshold:",
    round(correlation_threshold, 4),
    "\n"
  )
  
  cat(
    "Low-correlation samples:",
    nrow(correlation_flags),
    "\n"
  )
  
  if (nrow(correlation_flags) > 0) {
    print(correlation_flags)
  }
  
  correlation_plot <- ggplot(
    sample_correlation_qc,
    aes(
      x = Median_Correlation
    )
  ) +
    geom_histogram(
      bins = 40
    ) +
    theme_classic() +
    labs(
      title = "Proteomics Sample-to-Sample Correlation",
      subtitle = "Based on top 500 variable assays",
      x = "Median Pearson correlation",
      y = "Number of samples"
    )
  
  print(correlation_plot)
  
  ggsave(
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Sample_Correlation_Distribution.png"
    ),
    correlation_plot,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  write.csv(
    sample_correlation_qc,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Sample_Correlation_QC.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    correlation_flags,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Low_Correlation_Samples.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    sample_correlation,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Sample_Correlation_Matrix.rds"
    )
  )
  
  saveRDS(
    sample_correlation_qc,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Sample_Correlation_QC.rds"
    )
  )
  
  cat(
    "\nSample correlation QC saved successfully.\n"
  ) 

  cat("=== Cross-checking PCA and correlation outliers ===\n")
  
  qc_crosscheck <- sample_correlation_qc
  
  qc_crosscheck$PCA_Outlier <- qc_crosscheck$PCA_Mahalanobis_Distance >
    outlier_threshold
  
  qc_crosscheck$Patient_ID <- sub(
    "_.*$",
    "",
    qc_crosscheck$Patient_Timepoint
  )
  
  qc_crosscheck$Timepoint <- sub(
    "^[^_]+_",
    "",
    qc_crosscheck$Patient_Timepoint
  )
  
  qc_crosscheck$Correlation_Flag <- qc_crosscheck$Median_Correlation <
    correlation_threshold
  
  qc_crosscheck$Any_QC_Flag <- 
    qc_crosscheck$PCA_Outlier |
    qc_crosscheck$Correlation_Flag
  
  cat("\n=== QC flag summary ===\n")
  
  cat(
    "PCA outliers:",
    sum(qc_crosscheck$PCA_Outlier),
    "\n"
  )
  
  cat(
    "Low-correlation samples:",
    sum(qc_crosscheck$Correlation_Flag),
    "\n"
  )
  
  cat(
    "Samples flagged by either method:",
    sum(qc_crosscheck$Any_QC_Flag),
    "\n"
  )
  
  cat(
    "Samples flagged by BOTH methods:",
    sum(
      qc_crosscheck$PCA_Outlier &
        qc_crosscheck$Correlation_Flag
    ),
    "\n"
  )
  
  cat("\n=== Samples flagged by BOTH methods ===\n")
  
  both_flags <- qc_crosscheck[
    qc_crosscheck$PCA_Outlier &
      qc_crosscheck$Correlation_Flag,
    ,
    drop = FALSE
  ]
  
  print(
    both_flags[
      ,
      c(
        "Patient_Timepoint",
        "Patient_ID",
        "Timepoint",
        "Median_Correlation",
        "Mean_Correlation",
        "PCA_Mahalanobis_Distance"
      )
    ]
  )
  
  cat("\n=== All samples flagged by either method ===\n")
  
  all_flags <- qc_crosscheck[
    qc_crosscheck$Any_QC_Flag,
    ,
    drop = FALSE
  ]
  
  print(
    all_flags[
      order(
        -all_flags$PCA_Mahalanobis_Distance
      ),
      c(
        "Patient_Timepoint",
        "Patient_ID",
        "Timepoint",
        "Median_Correlation",
        "Mean_Correlation",
        "PCA_Mahalanobis_Distance",
        "PCA_Outlier",
        "Correlation_Flag"
      )
    ]
  )
  
  cat("\n=== Checking repeated timepoints for flagged patients ===\n")
  
  flagged_patients <- unique(
    all_flags$Patient_ID
  )
  
  flagged_patient_samples <- qc_crosscheck[
    qc_crosscheck$Patient_ID %in%
      flagged_patients,
    ,
    drop = FALSE
  ]
  
  print(
    flagged_patient_samples[
      order(
        flagged_patient_samples$Patient_ID,
        flagged_patient_samples$Timepoint
      ),
      c(
        "Patient_Timepoint",
        "Patient_ID",
        "Timepoint",
        "Median_Correlation",
        "PCA_Mahalanobis_Distance",
        "PCA_Outlier",
        "Correlation_Flag"
      )
    ]
  )
  
  cat("\n=== Within-patient correlation check ===\n")
  
  within_patient_results <- list()
  
  for (patient in flagged_patients) {
    
    patient_samples <- rownames(protein_qc_top)[
      sub(
        "_.*$",
        "",
        rownames(protein_qc_top)
      ) == patient
    ]
    
    if (length(patient_samples) >= 2) {
      
      patient_cor <- sample_correlation[
        patient_samples,
        patient_samples,
        drop = FALSE
      ]
      
      cor_values <- patient_cor[
        upper.tri(patient_cor)
      ]
      
      within_patient_results[[patient]] <- data.frame(
        Patient_ID = patient,
        N_Timepoints = length(patient_samples),
        Mean_Within_Patient_Correlation = mean(
          cor_values,
          na.rm = TRUE
        ),
        Minimum_Within_Patient_Correlation = min(
          cor_values,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  
  within_patient_qc <- do.call(
    rbind,
    within_patient_results
  )
  
  rownames(within_patient_qc) <- NULL
  
  cat("\nWithin-patient correlation results:\n")
  
  print(
    within_patient_qc[
      order(
        within_patient_qc$Mean_Within_Patient_Correlation
      ),
      ,
      drop = FALSE
    ]
  )
  
  cat("\n=== Final QC table ===\n")
  
  qc_final_flags <- merge(
    all_flags[
      ,
      c(
        "Patient_ID",
        "Timepoint",
        "Patient_Timepoint",
        "Median_Correlation",
        "Mean_Correlation",
        "PCA_Mahalanobis_Distance",
        "PCA_Outlier",
        "Correlation_Flag"
      )
    ],
    within_patient_qc,
    by = "Patient_ID",
    all.x = TRUE
  )
  
  qc_final_flags <- qc_final_flags[
    order(
      -qc_final_flags$PCA_Mahalanobis_Distance
    ),
    ,
    drop = FALSE
  ]
  
  print(qc_final_flags)
  
  write.csv(
    qc_final_flags,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Final_QC_Flagged_Samples_Crosscheck.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    qc_final_flags,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Final_QC_Flagged_Samples_Crosscheck.rds"
    )
  )
  
  saveRDS(
    within_patient_qc,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Within_Patient_Correlation_QC.rds"
    )
  )
  
  cat("\nCross-check QC saved successfully.\n")  
  
  
  cat("=== Detailed investigation of high-priority proteomics samples ===\n")
  
  investigation_ids <- c(
    "101_D7",
    "223_D3",
    "243_D0",
    "360_D3",
    "154_D0"
  )
  
  investigation_ids <- investigation_ids[
    investigation_ids %in%
      rownames(protein_matrix_collapsed)
  ]
  
  cat(
    "Samples selected for detailed investigation:",
    paste(investigation_ids, collapse = ", "),
    "\n"
  )
  
  protein_investigation <- protein_matrix_collapsed[
    investigation_ids,
    ,
    drop = FALSE
  ]
  
  sample_profile_qc <- data.frame(
    Patient_Timepoint = investigation_ids,
    Median = apply(
      protein_investigation,
      1,
      median,
      na.rm = TRUE
    ),
    Mean = apply(
      protein_investigation,
      1,
      mean,
      na.rm = TRUE
    ),
    SD = apply(
      protein_investigation,
      1,
      sd,
      na.rm = TRUE
    ),
    IQR = apply(
      protein_investigation,
      1,
      IQR,
      na.rm = TRUE
    ),
    Minimum = apply(
      protein_investigation,
      1,
      min,
      na.rm = TRUE
    ),
    Maximum = apply(
      protein_investigation,
      1,
      max,
      na.rm = TRUE
    ),
    Missing = apply(
      protein_investigation,
      1,
      function(x) sum(is.na(x))
    ),
    stringsAsFactors = FALSE
  )
  
  cat("\n=== Protein distribution profiles ===\n")
  print(sample_profile_qc)
  
  cat("\n=== Comparison of Patient 101 timepoints ===\n")
  
  patient_101_ids <- c(
    "101_D0",
    "101_D3",
    "101_D7"
  )
  
  patient_101_ids <- patient_101_ids[
    patient_101_ids %in%
      rownames(protein_matrix_collapsed)
  ]
  
  patient_101_matrix <- protein_matrix_collapsed[
    patient_101_ids,
    ,
    drop = FALSE
  ]
  
  patient_101_summary <- data.frame(
    Patient_Timepoint = patient_101_ids,
    Median = apply(
      patient_101_matrix,
      1,
      median,
      na.rm = TRUE
    ),
    Mean = apply(
      patient_101_matrix,
      1,
      mean,
      na.rm = TRUE
    ),
    SD = apply(
      patient_101_matrix,
      1,
      sd,
      na.rm = TRUE
    ),
    IQR = apply(
      patient_101_matrix,
      1,
      IQR,
      na.rm = TRUE
    ),
    Missing = apply(
      patient_101_matrix,
      1,
      function(x) sum(is.na(x))
    ),
    stringsAsFactors = FALSE
  )
  
  print(patient_101_summary)
  
  cat("\n=== Patient 101 pairwise correlations ===\n")
  
  print(
    round(
      cor(
        t(patient_101_matrix),
        method = "pearson",
        use = "pairwise.complete.obs"
      ),
      4
    )
  )
  
  cat("\n=== Patient 101 absolute protein-level differences ===\n")
  
  if (all(c("101_D0", "101_D7") %in% patient_101_ids)) {
    
    diff_101_D7_D0 <- abs(
      patient_101_matrix["101_D7", ] -
        patient_101_matrix["101_D0", ]
    )
    
    cat(
      "101_D7 vs 101_D0 median absolute difference:",
      round(
        median(
          diff_101_D7_D0,
          na.rm = TRUE
        ),
        4
      ),
      "\n"
    )
    
    cat(
      "101_D7 vs 101_D0 mean absolute difference:",
      round(
        mean(
          diff_101_D7_D0,
          na.rm = TRUE
        ),
        4
      ),
      "\n"
    )
    
    cat(
      "101_D7 vs 101_D0 95th percentile absolute difference:",
      round(
        quantile(
          diff_101_D7_D0,
          0.95,
          na.rm = TRUE
        ),
        4
      ),
      "\n"
    )
    
    cat(
      "101_D7 vs 101_D0 maximum absolute difference:",
      round(
        max(
          diff_101_D7_D0,
          na.rm = TRUE
        ),
        4
      ),
      "\n"
    )
  }
  
  cat("\n=== Comparing flagged samples with global distribution ===\n")
  
  global_profile <- data.frame(
    Patient_Timepoint = rownames(protein_matrix_collapsed),
    Median = apply(
      protein_matrix_collapsed,
      1,
      median,
      na.rm = TRUE
    ),
    Mean = apply(
      protein_matrix_collapsed,
      1,
      mean,
      na.rm = TRUE
    ),
    SD = apply(
      protein_matrix_collapsed,
      1,
      sd,
      na.rm = TRUE
    ),
    IQR = apply(
      protein_matrix_collapsed,
      1,
      IQR,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
  
  sample_profile_qc$Median_Z <- (
    sample_profile_qc$Median -
      median(global_profile$Median, na.rm = TRUE)
  ) / mad(
    global_profile$Median,
    na.rm = TRUE
  )
  
  sample_profile_qc$Mean_Z <- (
    sample_profile_qc$Mean -
      median(global_profile$Mean, na.rm = TRUE)
  ) / mad(
    global_profile$Mean,
    na.rm = TRUE
  )
  
  sample_profile_qc$SD_Z <- (
    sample_profile_qc$SD -
      median(global_profile$SD, na.rm = TRUE)
  ) / mad(
    global_profile$SD,
    na.rm = TRUE
  )
  
  cat("\n=== Robust distribution scores ===\n")
  print(sample_profile_qc)
  
  write.csv(
    sample_profile_qc,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_High_Priority_Sample_Investigation.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    patient_101_summary,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Patient_101_Timepoint_Profile.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    sample_profile_qc,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_High_Priority_Sample_Investigation.rds"
    )
  )
  
  saveRDS(
    patient_101_summary,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Patient_101_Timepoint_Profile.rds"
    )
  )
  
  cat("\nDetailed sample investigation saved successfully.\n")
  
  
  cat("=== Patient 101 corrected protein-level difference analysis ===\n")
  
  patient_101_numeric <- as.matrix(
    patient_101_matrix
  )
  
  storage.mode(patient_101_numeric) <- "numeric"
  
  diff_101_D7_D0 <- abs(
    patient_101_numeric["101_D7", ] -
      patient_101_numeric["101_D0", ]
  )
  
  diff_101_D7_D3 <- abs(
    patient_101_numeric["101_D7", ] -
      patient_101_numeric["101_D3", ]
  )
  
  cat("\n=== 101_D7 vs 101_D0 ===\n")
  
  cat(
    "Median absolute difference:",
    round(
      median(diff_101_D7_D0, na.rm = TRUE),
      4
    ),
    "\n"
  )
  
  cat(
    "Mean absolute difference:",
    round(
      mean(diff_101_D7_D0, na.rm = TRUE),
      4
    ),
    "\n"
  )
  
  cat(
    "95th percentile absolute difference:",
    round(
      quantile(
        diff_101_D7_D0,
        0.95,
        na.rm = TRUE
      ),
      4
    ),
    "\n"
  )
  
  cat(
    "Maximum absolute difference:",
    round(
      max(diff_101_D7_D0, na.rm = TRUE),
      4
    ),
    "\n"
  )
  
  cat("\n=== 101_D7 vs 101_D3 ===\n")
  
  cat(
    "Median absolute difference:",
    round(
      median(diff_101_D7_D3, na.rm = TRUE),
      4
    ),
    "\n"
  )
  
  cat(
    "Mean absolute difference:",
    round(
      mean(diff_101_D7_D3, na.rm = TRUE),
      4
    ),
    "\n"
  )
  
  cat(
    "95th percentile absolute difference:",
    round(
      quantile(
        diff_101_D7_D3,
        0.95,
        na.rm = TRUE
      ),
      4
    ),
    "\n"
  )
  
  cat(
    "Maximum absolute difference:",
    round(
      max(diff_101_D7_D3, na.rm = TRUE),
      4
    ),
    "\n"
  )
  
  cat("\n=== Checking numeric validity ===\n")
  
  cat(
    "D7-D0 numeric:",
    is.numeric(diff_101_D7_D0),
    "\n"
  )
  
  cat(
    "D7-D3 numeric:",
    is.numeric(diff_101_D7_D3),
    "\n"
  )
  
  cat(
    "D7-D0 valid proteins:",
    sum(is.finite(diff_101_D7_D0)),
    "\n"
  )
  
  cat(
    "D7-D3 valid proteins:",
    sum(is.finite(diff_101_D7_D3)),
    "\n"
  )
  
  cat("\nSaving corrected Patient 101 investigation...\n")
  
  patient_101_difference_summary <- data.frame(
    Comparison = c(
      "101_D7_vs_101_D0",
      "101_D7_vs_101_D3"
    ),
    Median_Absolute_Difference = c(
      median(diff_101_D7_D0, na.rm = TRUE),
      median(diff_101_D7_D3, na.rm = TRUE)
    ),
    Mean_Absolute_Difference = c(
      mean(diff_101_D7_D0, na.rm = TRUE),
      mean(diff_101_D7_D3, na.rm = TRUE)
    ),
    P95_Absolute_Difference = c(
      quantile(
        diff_101_D7_D0,
        0.95,
        na.rm = TRUE
      ),
      quantile(
        diff_101_D7_D3,
        0.95,
        na.rm = TRUE
      )
    ),
    Maximum_Absolute_Difference = c(
      max(diff_101_D7_D0, na.rm = TRUE),
      max(diff_101_D7_D3, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    patient_101_difference_summary,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Patient_101_Absolute_Difference_Summary.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    patient_101_difference_summary,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Patient_101_Absolute_Difference_Summary.rds"
    )
  )
  
  cat(
    "\nCorrected Patient 101 analysis saved successfully.\n"
  )  
  
  cat("=== Final proteomics QC candidate classification ===\n")
  
  qc_candidates <- qc_crosscheck[
    qc_crosscheck$Any_QC_Flag,
    ,
    drop = FALSE
  ]
  
  qc_candidates$QC_Category <- "Investigate"
  
  qc_candidates$QC_Category[
    qc_candidates$PCA_Outlier &
      qc_candidates$Correlation_Flag
  ] <- "High-priority investigation"
  
  qc_candidates$QC_Category[
    qc_candidates$PCA_Outlier &
      !qc_candidates$Correlation_Flag
  ] <- "PCA-only candidate"
  
  qc_candidates$QC_Category[
    !qc_candidates$PCA_Outlier &
      qc_candidates$Correlation_Flag
  ] <- "Correlation-only candidate"
  
  qc_candidates <- qc_candidates[
    order(
      factor(
        qc_candidates$QC_Category,
        levels = c(
          "High-priority investigation",
          "PCA-only candidate",
          "Correlation-only candidate"
        )
      ),
      qc_candidates$Median_Correlation
    ),
    ,
    drop = FALSE
  ]
  
  cat("\n=== Candidate classification ===\n")
  
  print(
    qc_candidates[
      ,
      c(
        "Patient_Timepoint",
        "Patient_ID",
        "Timepoint",
        "Median_Correlation",
        "Mean_Correlation",
        "PCA_Mahalanobis_Distance",
        "PCA_Outlier",
        "Correlation_Flag",
        "QC_Category"
      )
    ]
  )
  
  cat("\n=== QC category counts ===\n")
  
  print(
    table(
      qc_candidates$QC_Category
    )
  )
  
  cat("\n=== High-priority samples ===\n")
  
  high_priority <- qc_candidates[
    qc_candidates$QC_Category ==
      "High-priority investigation",
    ,
    drop = FALSE
  ]
  
  print(high_priority)
  
  cat("\n=== Correlation-only candidates ===\n")
  
  correlation_only <- qc_candidates[
    qc_candidates$QC_Category ==
      "Correlation-only candidate",
    ,
    drop = FALSE
  ]
  
  print(correlation_only)
  
  cat("\n=== PCA-only candidates ===\n")
  
  pca_only <- qc_candidates[
    qc_candidates$QC_Category ==
      "PCA-only candidate",
    ,
    drop = FALSE
  ]
  
  print(pca_only)
  
  write.csv(
    qc_candidates,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Final_QC_Candidate_Classification.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    qc_candidates,
    file.path(
      proteomics_results_dir,
      "Olink_Proteomics_Final_QC_Candidate_Classification.rds"
    )
  )
  
  cat("\nFinal QC candidate classification saved successfully.\n")  

  
  cat("=== Proteomics Differential Analysis: D0 COVID-positive vs COVID-negative ===\n")
  
  library(limma)
  
  cat("\n=== Preparing D0 dataset ===\n")
  
  d0_keys <- rownames(protein_matrix_collapsed)
  
  d0_metadata <- master_metadata_final[
    match(
      d0_keys,
      paste(
        master_metadata_final$Patient_ID,
        master_metadata_final$Timepoint,
        sep = "_"
      )
    ),
    ,
    drop = FALSE
  ]
  
  d0_metadata <- d0_metadata[
    as.character(d0_metadata$Timepoint) == "D0",
    ,
    drop = FALSE
  ]
  
  d0_keys_final <- paste(
    d0_metadata$Patient_ID,
    d0_metadata$Timepoint,
    sep = "_"
  )
  
  d0_matrix <- protein_matrix_collapsed[
    d0_keys_final,
    ,
    drop = FALSE
  ]
  
  rownames(d0_metadata) <- d0_keys_final
  
  cat(
    "D0 samples:",
    nrow(d0_matrix),
    "\n"
  )
  
  cat(
    "D0 assays:",
    ncol(d0_matrix),
    "\n"
  )
  
  cat(
    "Missing values:",
    sum(is.na(d0_matrix)),
    "\n"
  )
  
  cat("\n=== COVID group distribution ===\n")
  
  d0_metadata$COVID_Group <- ifelse(
    as.character(d0_metadata$COVID) == "1",
    "COVID-positive",
    "COVID-negative"
  )
  
  print(
    table(
      d0_metadata$COVID_Group,
      useNA = "ifany"
    )
  )
  
  cat("\n=== QC flags in D0 dataset ===\n")
  
  d0_metadata$PCA_Outlier <- qc_crosscheck$PCA_Outlier[
    match(
      rownames(d0_metadata),
      qc_crosscheck$Patient_Timepoint
    )
  ]
  
  d0_metadata$Correlation_Flag <- qc_crosscheck$Correlation_Flag[
    match(
      rownames(d0_metadata),
      qc_crosscheck$Patient_Timepoint
    )
  ]
  
  d0_metadata$Any_QC_Flag <- 
    d0_metadata$PCA_Outlier |
    d0_metadata$Correlation_Flag
  
  print(
    table(
      d0_metadata$Any_QC_Flag,
      useNA = "ifany"
    )
  )
  
  cat("\n=== Building limma model ===\n")
  
  covid_factor <- factor(
    d0_metadata$COVID_Group,
    levels = c(
      "COVID-negative",
      "COVID-positive"
    )
  )
  
  design_d0 <- model.matrix(
    ~ covid_factor
  )
  
  colnames(design_d0) <- c(
    "Intercept",
    "COVID_positive_vs_negative"
  )
  
  cat(
    "Design dimensions:",
    nrow(design_d0),
    "x",
    ncol(design_d0),
    "\n"
  )
  
  cat(
    "Design rank:",
    qr(design_d0)$rank,
    "\n"
  )
  
  print(
    colnames(design_d0)
  )
  
  cat("\n=== Fitting limma model ===\n")
  
  protein_matrix_d0_limma <- t(
    d0_matrix
  )
  
  fit_d0 <- lmFit(
    protein_matrix_d0_limma,
    design_d0
  )
  
  fit_d0 <- eBayes(
    fit_d0,
    robust = TRUE
  )
  
  results_d0 <- topTable(
    fit_d0,
    coef = "COVID_positive_vs_negative",
    number = Inf,
    sort.by = "P"
  )
  
  results_d0$OlinkID <- rownames(results_d0)
  
  rownames(results_d0) <- NULL
  
  cat("\n=== Differential analysis completed ===\n")
  
  cat(
    "Assays tested:",
    nrow(results_d0),
    "\n"
  )
  
  cat(
    "Finite P values:",
    sum(is.finite(results_d0$P.Value)),
    "\n"
  )
  
  cat(
    "FDR < 0.05:",
    sum(
      results_d0$adj.P.Val < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "FDR < 0.05 and |logFC| >= 1:",
    sum(
      results_d0$adj.P.Val < 0.05 &
        abs(results_d0$logFC) >= 1,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat("\n=== Top 20 differential proteins ===\n")
  
  print(
    results_d0[
      ,
      c(
        "OlinkID",
        "logFC",
        "AveExpr",
        "t",
        "P.Value",
        "adj.P.Val"
      )
    ][
      seq_len(
        min(20, nrow(results_d0))
      ),
      ,
      drop = FALSE
    ]
  )
  
  cat("\n=== Adding Olink assay annotation ===\n")
  
  assay_annotation_d0 <- proteomics_assay_annotation
  
  results_d0_annotated <- merge(
    results_d0,
    assay_annotation_d0,
    by = "OlinkID",
    all.x = TRUE,
    sort = FALSE
  )
  
  results_d0_annotated <- results_d0_annotated[
    match(
      results_d0$OlinkID,
      results_d0_annotated$OlinkID
    ),
    ,
    drop = FALSE
  ]
  
  rownames(results_d0_annotated) <- NULL
  
  cat(
    "Annotation rows:",
    nrow(results_d0_annotated),
    "\n"
  )
  
  cat(
    "Missing Assay annotation:",
    sum(
      is.na(results_d0_annotated$Assay)
    ),
    "\n"
  )
  
  cat(
    "Missing UniProt annotation:",
    sum(
      is.na(results_d0_annotated$UniProt)
    ),
    "\n"
  )
  
  cat("\n=== Creating significant/up/down tables ===\n")
  
  significant_d0 <- results_d0_annotated[
    results_d0_annotated$adj.P.Val < 0.05 &
      abs(results_d0_annotated$logFC) >= 1,
    ,
    drop = FALSE
  ]
  
  upregulated_d0 <- significant_d0[
    significant_d0$logFC > 0,
    ,
    drop = FALSE
  ]
  
  downregulated_d0 <- significant_d0[
    significant_d0$logFC < 0,
    ,
    drop = FALSE
  ]
  
  cat(
    "Significant proteins:",
    nrow(significant_d0),
    "\n"
  )
  
  cat(
    "Upregulated:",
    nrow(upregulated_d0),
    "\n"
  )
  
  cat(
    "Downregulated:",
    nrow(downregulated_d0),
    "\n"
  )
  
  cat("\n=== Saving D0 differential analysis ===\n")
  
  d0_de_dir <- file.path(
    proteomics_results_dir,
    "Differential_Expression"
  )
  
  dir.create(
    d0_de_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write.csv(
    results_d0_annotated,
    file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_All_Assays.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    significant_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_Significant.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    upregulated_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_Upregulated.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    downregulated_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_Downregulated.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    results_d0_annotated,
    file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_All_Assays.rds"
    )
  )
  
  saveRDS(
    significant_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_Significant.rds"
    )
  )
  
  saveRDS(
    upregulated_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_Upregulated.rds"
    )
  )
  
  saveRDS(
    downregulated_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_Downregulated.rds"
    )
  )
  
  saveRDS(
    fit_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_limma_Fit.rds"
    )
  )
  
  saveRDS(
    d0_metadata,
    file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_Metadata.rds"
    )
  )
  
  cat(
    "\nD0 proteomics differential analysis saved successfully.\n"
  )  
  
  cat("=== D0 Proteomics Visualization ===\n")
  
  library(ggplot2)
  library(pheatmap)
  
  cat("\n=== Preparing volcano plot ===\n")
  
  volcano_data <- results_d0_annotated
  
  volcano_data$Significance <- "Not significant"
  
  volcano_data$Significance[
    volcano_data$adj.P.Val < 0.05 &
      abs(volcano_data$logFC) >= 1
  ] <- "Significant"
  
  cat(
    "Total assays:",
    nrow(volcano_data),
    "\n"
  )
  
  cat(
    "Significant assays:",
    sum(
      volcano_data$Significance == "Significant"
    ),
    "\n"
  )
  
  volcano_data$MinusLog10FDR <- -log10(
    volcano_data$adj.P.Val
  )
  
  cat("\n=== Creating volcano plot ===\n")
  
  volcano_plot_d0 <- ggplot(
    volcano_data,
    aes(
      x = logFC,
      y = MinusLog10FDR
    )
  ) +
    geom_point(
      aes(
        shape = Significance
      ),
      alpha = 0.65,
      size = 2
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    theme_classic() +
    labs(
      title = "D0 Olink Proteomics",
      subtitle = "COVID-positive vs COVID-negative",
      x = "logFC",
      y = "-log10(FDR)",
      shape = "Classification"
    )
  
  print(volcano_plot_d0)
  
  ggsave(
    filename = file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_Volcano.png"
    ),
    plot = volcano_plot_d0,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  cat("Volcano plot saved.\n")
  
  cat("\n=== Preparing significant protein heatmap ===\n")
  
  heatmap_ids <- significant_d0$OlinkID
  
  heatmap_matrix <- d0_matrix[
    ,
    heatmap_ids,
    drop = FALSE
  ]
  
  heatmap_matrix <- t(
    heatmap_matrix
  )
  
  cat(
    "Heatmap proteins:",
    nrow(heatmap_matrix),
    "\n"
  )
  
  cat(
    "Heatmap samples:",
    ncol(heatmap_matrix),
    "\n"
  )
  
  cat(
    "Missing values before heatmap processing:",
    sum(is.na(heatmap_matrix)),
    "\n"
  )
  
  heatmap_matrix_imputed <- heatmap_matrix
  
  for (i in seq_len(nrow(heatmap_matrix_imputed))) {
    
    missing_idx <- is.na(
      heatmap_matrix_imputed[i, ]
    )
    
    if (any(missing_idx)) {
      
      heatmap_matrix_imputed[i, missing_idx] <- median(
        heatmap_matrix_imputed[i, ],
        na.rm = TRUE
      )
      
    }
  }
  
  cat(
    "Missing values after imputation:",
    sum(is.na(heatmap_matrix_imputed)),
    "\n"
  )
  
  cat("\n=== Scaling heatmap data ===\n")
  
  heatmap_matrix_scaled <- t(
    scale(
      t(heatmap_matrix_imputed)
    )
  )
  
  cat(
    "Missing values after scaling:",
    sum(is.na(heatmap_matrix_scaled)),
    "\n"
  )
  
  heatmap_annotation <- data.frame(
    COVID = d0_metadata$COVID_Group
  )
  
  rownames(
    heatmap_annotation
  ) <- rownames(d0_metadata)
  
  cat("\n=== Creating heatmap ===\n")
  
  pheatmap(
    heatmap_matrix_scaled,
    annotation_col = heatmap_annotation,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 7,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "D0 Significant Olink Proteins"
  )
  
  cat("\n=== Saving heatmap ===\n")
  
  png(
    filename = file.path(
      d0_de_dir,
      "Olink_D0_COVID_Positive_vs_Negative_Significant_Protein_Heatmap.png"
    ),
    width = 2600,
    height = 2200,
    res = 300
  )
  
  pheatmap(
    heatmap_matrix_scaled,
    annotation_col = heatmap_annotation,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 7,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "D0 Significant Olink Proteins"
  )
  
  dev.off()
  
  cat("Heatmap saved.\n")
  
  cat("\n=== Saving visualization data ===\n")
  
  write.csv(
    volcano_data,
    file.path(
      d0_de_dir,
      "Olink_D0_Volcano_Data.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    heatmap_matrix_scaled,
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_Protein_Heatmap_Matrix.csv"
    ),
    row.names = TRUE
  )
  
  saveRDS(
    heatmap_matrix_scaled,
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_Protein_Heatmap_Matrix.rds"
    )
  )
  
  cat("\n=== D0 visualization completed successfully ===\n")
  
  cat(
    "Volcano PNG: Olink_D0_COVID_Positive_vs_Negative_Volcano.png\n"
  )
  
  cat(
    "Heatmap PNG: Olink_D0_COVID_Positive_vs_Negative_Significant_Protein_Heatmap.png\n"
  )
  
  cat(
    "Visualization data saved in:",
    d0_de_dir,
    "\n"
  )
  
  
  cat("=== Recreating D0 proteomics heatmap ===\n")
  
  library(pheatmap)
  library(grid)
  
  heatmap_file <- file.path(
    d0_de_dir,
    "Olink_D0_COVID_Positive_vs_Negative_Significant_Protein_Heatmap_v2.png"
  )
  
  cat(
    "Heatmap matrix:",
    nrow(heatmap_matrix_scaled),
    "proteins x",
    ncol(heatmap_matrix_scaled),
    "samples\n"
  )
  
  cat(
    "Missing values:",
    sum(is.na(heatmap_matrix_scaled)),
    "\n"
  )
  
  heatmap_obj <- pheatmap(
    heatmap_matrix_scaled,
    annotation_col = heatmap_annotation,
    show_colnames = FALSE,
    show_rownames = TRUE,
    fontsize_row = 7,
    fontsize_col = 6,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    silent = TRUE,
    main = "D0 Significant Olink Proteins"
  )
  
  png(
    filename = heatmap_file,
    width = 3000,
    height = 2400,
    res = 300
  )
  
  grid::grid.newpage()
  grid::grid.draw(heatmap_obj$gtable)
  
  dev.off()
  
  cat(
    "Heatmap saved:",
    file.exists(heatmap_file),
    "\n"
  )
  
  cat(
    "Heatmap file size:",
    file.info(heatmap_file)$size,
    "bytes\n"
  )
  
  cat(
    "File:",
    heatmap_file,
    "\n"
  ) 
  
  cat("=== Preparing D0 Proteomics IDs for Enrichment ===\n")
  
  d0_sig <- significant_d0_annotated
  
  cat(
    "Significant assays:",
    nrow(d0_sig),
    "\n"
  )
  
  cat(
    "Unique Olink IDs:",
    length(unique(d0_sig$OlinkID)),
    "\n"
  )
  
  cat(
    "Unique UniProt IDs:",
    length(unique(d0_sig$UniProt)),
    "\n"
  )
  
  cat("\n=== Mapping UniProt to Entrez ===\n")
  
  uniprot_ids <- unique(
    d0_sig$UniProt[
      !is.na(d0_sig$UniProt) &
        d0_sig$UniProt != ""
    ]
  )
  
  uniprot_mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = uniprot_ids,
    keytype = "UNIPROT",
    columns = c(
      "UNIPROT",
      "ENTREZID",
      "SYMBOL",
      "GENENAME"
    )
  )
  
  uniprot_mapping <- uniprot_mapping[
    !is.na(uniprot_mapping$ENTREZID),
  ]
  
  uniprot_mapping <- uniprot_mapping[
    !duplicated(uniprot_mapping$UNIPROT),
  ]
  
  cat(
    "Unique UniProt IDs submitted:",
    length(uniprot_ids),
    "\n"
  )
  
  cat(
    "UniProt IDs mapped to Entrez:",
    length(unique(uniprot_mapping$UNIPROT)),
    "\n"
  )
  
  cat(
    "Unique Entrez IDs:",
    length(unique(uniprot_mapping$ENTREZID)),
    "\n"
  )
  
  cat(
    "Mapping percentage:",
    round(
      100 *
        length(unique(uniprot_mapping$UNIPROT)) /
        length(uniprot_ids),
      2
    ),
    "%\n"
  )
  
  cat("\n=== Mapping preview ===\n")
  
  print(
    head(
      uniprot_mapping,
      15
    )
  )
  
  cat("\n=== Saving mapping ===\n")
  
  write.csv(
    uniprot_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_UniProt_to_Entrez_Mapping.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    uniprot_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_UniProt_to_Entrez_Mapping.rds"
    )
  )
  
  cat("\n=== ID mapping completed ===\n") 
  
  cat("=== Loading human gene annotation database ===\n")
  
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install(
      "org.Hs.eg.db",
      ask = FALSE,
      update = FALSE
    )
  }
  
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
  cat(
    "org.Hs.eg.db loaded successfully.\n"
  )
  
  cat(
    "Annotation package version:",
    as.character(
      packageVersion("org.Hs.eg.db")
    ),
    "\n"
  )
  
  cat(
    "AnnotationDbi loaded:",
    requireNamespace("AnnotationDbi", quietly = TRUE),
    "\n"
  ) 
  
  cat("=== D0 Proteomics: UniProt to Entrez Mapping ===\n")
  
  uniprot_ids <- unique(
    d0_sig$UniProt[
      !is.na(d0_sig$UniProt) &
        d0_sig$UniProt != ""
    ]
  )
  
  cat(
    "Unique UniProt IDs submitted:",
    length(uniprot_ids),
    "\n"
  )
  
  uniprot_mapping_raw <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = uniprot_ids,
    keytype = "UNIPROT",
    columns = c(
      "UNIPROT",
      "ENTREZID",
      "SYMBOL",
      "GENENAME"
    )
  )
  
  uniprot_mapping_raw <- uniprot_mapping_raw[
    !is.na(uniprot_mapping_raw$ENTREZID),
  ]
  
  uniprot_mapping <- uniprot_mapping_raw[
    !duplicated(uniprot_mapping_raw$UNIPROT),
  ]
  
  mapped_uniprot <- length(
    unique(uniprot_mapping$UNIPROT)
  )
  
  unique_entrez <- length(
    unique(uniprot_mapping$ENTREZID)
  )
  
  mapping_percentage <- round(
    100 * mapped_uniprot / length(uniprot_ids),
    2
  )
  
  cat(
    "UniProt IDs mapped to Entrez:",
    mapped_uniprot,
    "\n"
  )
  
  cat(
    "Unique Entrez IDs:",
    unique_entrez,
    "\n"
  )
  
  cat(
    "Mapping percentage:",
    mapping_percentage,
    "%\n"
  )
  
  unmapped_uniprot <- setdiff(
    uniprot_ids,
    uniprot_mapping$UNIPROT
  )
  
  cat(
    "Unmapped UniProt IDs:",
    length(unmapped_uniprot),
    "\n"
  )
  
  if (length(unmapped_uniprot) > 0) {
    print(unmapped_uniprot)
  }
  
  cat("\n=== Mapping preview ===\n")
  
  print(uniprot_mapping)
  
  cat("\n=== Saving mapping results ===\n")
  
  write.csv(
    uniprot_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_UniProt_to_Entrez_Mapping.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      UniProt = unmapped_uniprot
    ),
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_Unmapped_UniProt.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    uniprot_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_UniProt_to_Entrez_Mapping.rds"
    )
  )
  
  cat("\n=== UniProt to Entrez mapping completed ===\n")
  
  cat("=== D0 Proteomics: Background Universe Mapping ===\n")
  
  all_uniprot_ids <- unique(
    proteomics_assay_annotation$UniProt[
      !is.na(proteomics_assay_annotation$UniProt) &
        proteomics_assay_annotation$UniProt != ""
    ]
  )
  
  cat(
    "Unique UniProt IDs in all 1429 assays:",
    length(all_uniprot_ids),
    "\n"
  )
  
  background_mapping_raw <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = all_uniprot_ids,
    keytype = "UNIPROT",
    columns = c(
      "UNIPROT",
      "ENTREZID",
      "SYMBOL",
      "GENENAME"
    )
  )
  
  background_mapping_raw <- background_mapping_raw[
    !is.na(background_mapping_raw$ENTREZID),
  ]
  
  background_mapping <- background_mapping_raw[
    !duplicated(background_mapping_raw$UNIPROT),
  ]
  
  background_entrez <- unique(
    background_mapping$ENTREZID
  )
  
  cat(
    "UniProt IDs mapped:",
    length(unique(background_mapping$UNIPROT)),
    "\n"
  )
  
  cat(
    "Unique Entrez IDs in background:",
    length(background_entrez),
    "\n"
  )
  
  cat(
    "Background mapping percentage:",
    round(
      100 * length(unique(background_mapping$UNIPROT)) /
        length(all_uniprot_ids),
      2
    ),
    "%\n"
  )
  
  cat(
    "Unmapped UniProt IDs:",
    length(
      setdiff(
        all_uniprot_ids,
        background_mapping$UNIPROT
      )
    ),
    "\n"
  )
  
  write.csv(
    background_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Background_UniProt_to_Entrez_Mapping.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    background_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Background_UniProt_to_Entrez_Mapping.rds"
    )
  )
  
  write.csv(
    data.frame(
      EntrezID = background_entrez
    ),
    file.path(
      d0_de_dir,
      "Olink_D0_Background_Entrez_Universe.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    background_entrez,
    file.path(
      d0_de_dir,
      "Olink_D0_Background_Entrez_Universe.rds"
    )
  )
  
  cat(
    "\n=== Background universe mapping completed ===\n"
  )  
  
  
  cat("=== D0 Proteomics: GO Biological Process Enrichment ===\n")
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    BiocManager::install(
      "clusterProfiler",
      ask = FALSE,
      update = FALSE
    )
  }
  
  library(clusterProfiler)
  
  significant_entrez <- unique(
    uniprot_mapping$ENTREZID
  )
  
  background_entrez <- unique(
    background_entrez
  )
  
  cat(
    "Significant Entrez IDs:",
    length(significant_entrez),
    "\n"
  )
  
  cat(
    "Background Entrez IDs:",
    length(background_entrez),
    "\n"
  )
  
  ego_d0 <- enrichGO(
    gene = significant_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  ego_d0_df <- as.data.frame(ego_d0)
  
  cat(
    "Significant GO-BP pathways:",
    nrow(ego_d0_df),
    "\n"
  )
  
  if (nrow(ego_d0_df) > 0) {
    cat("\n=== Top GO-BP pathways ===\n")
    
    print(
      ego_d0_df[
        order(ego_d0_df$p.adjust),
        c(
          "ID",
          "Description",
          "GeneRatio",
          "BgRatio",
          "Count",
          "pvalue",
          "p.adjust",
          "qvalue",
          "geneID"
        )
      ][1:min(20, nrow(ego_d0_df)), ]
    )
  }
  
  write.csv(
    ego_d0_df,
    file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    ego_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment.rds"
    )
  )
  
  cat(
    "\nSaved:",
    file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment.csv"
    ),
    "\n"
  )
  
  cat(
    "=== GO Biological Process enrichment completed ===\n"
  )  

  cat("=== D0 Proteomics: KEGG Enrichment ===\n")
  
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install(
      "org.Hs.eg.db",
      ask = FALSE,
      update = FALSE
    )
  }
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  ekegg_d0 <- enrichKEGG(
    gene = significant_entrez,
    universe = background_entrez,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
  
  ekegg_d0_df <- as.data.frame(ekegg_d0)
  
  cat(
    "Significant KEGG pathways:",
    nrow(ekegg_d0_df),
    "\n"
  )
  
  if (nrow(ekegg_d0_df) > 0) {
    cat("\n=== Top KEGG pathways ===\n")
    
    print(
      ekegg_d0_df[
        order(ekegg_d0_df$p.adjust),
        c(
          "ID",
          "Description",
          "GeneRatio",
          "BgRatio",
          "Count",
          "pvalue",
          "p.adjust",
          "qvalue",
          "geneID"
        )
      ][1:min(20, nrow(ekegg_d0_df)), ]
    )
  }
  
  write.csv(
    ekegg_d0_df,
    file.path(
      d0_de_dir,
      "Olink_D0_KEGG_Enrichment.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    ekegg_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_KEGG_Enrichment.rds"
    )
  )
  
  cat(
    "\nSaved:",
    file.path(
      d0_de_dir,
      "Olink_D0_KEGG_Enrichment.csv"
    ),
    "\n"
  )
  
  cat("=== KEGG enrichment completed ===\n") 
  
  cat("=== D0 Proteomics: Reactome Enrichment ===\n")
  
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    BiocManager::install(
      "ReactomePA",
      ask = FALSE,
      update = FALSE
    )
  }
  
  library(ReactomePA)
  
  ereactome_d0 <- enrichPathway(
    gene = significant_entrez,
    universe = background_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  ereactome_d0_df <- as.data.frame(ereactome_d0)
  
  cat(
    "Significant Reactome pathways:",
    nrow(ereactome_d0_df),
    "\n"
  )
  
  if (nrow(ereactome_d0_df) > 0) {
    cat("\n=== Top Reactome pathways ===\n")
    
    print(
      ereactome_d0_df[
        order(ereactome_d0_df$p.adjust),
        c(
          "ID",
          "Description",
          "GeneRatio",
          "BgRatio",
          "Count",
          "pvalue",
          "p.adjust",
          "qvalue",
          "geneID"
        )
      ][1:min(20, nrow(ereactome_d0_df)), ]
    )
  }
  
  write.csv(
    ereactome_d0_df,
    file.path(
      d0_de_dir,
      "Olink_D0_Reactome_Enrichment.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    ereactome_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_Reactome_Enrichment.rds"
    )
  )
  
  cat(
    "\nSaved:",
    file.path(
      d0_de_dir,
      "Olink_D0_Reactome_Enrichment.csv"
    ),
    "\n"
  )
  
  cat("=== Reactome enrichment completed ===\n")
  
  
  cat("=== D0 Proteomics: GO-BP Visualization ===\n")
  
  library(ggplot2)
  
  go_plot_df <- ego_d0_df
  
  go_plot_df <- go_plot_df[
    order(go_plot_df$p.adjust),
  ]
  
  go_plot_df$Description <- factor(
    go_plot_df$Description,
    levels = rev(go_plot_df$Description)
  )
  
  go_plot_df$minus_log10_padj <- -log10(
    go_plot_df$p.adjust
  )
  
  go_plot <- ggplot(
    go_plot_df,
    aes(
      x = minus_log10_padj,
      y = Description,
      size = Count
    )
  ) +
    geom_point() +
    labs(
      title = "D0 Proteomics: GO Biological Process Enrichment",
      x = "-log10 adjusted P-value",
      y = NULL,
      size = "Protein count"
    ) +
    theme_minimal(base_size = 12)
  
  print(go_plot)
  
  ggsave(
    filename = file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment_DotPlot.png"
    ),
    plot = go_plot,
    width = 9,
    height = 5,
    dpi = 300
  )
  
  write.csv(
    go_plot_df,
    file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment_Plot_Data.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    go_plot,
    file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment_DotPlot.rds"
    )
  )
  
  cat(
    "Plot saved:",
    file.exists(
      file.path(
        d0_de_dir,
        "Olink_D0_GO_BP_Enrichment_DotPlot.png"
      )
    ),
    "\n"
  )
  
  cat("=== GO-BP visualization completed ===\n") 
  
  
  cat("=== D0 Proteomics: Final Significant Protein Annotation ===\n")
  
  final_d0_protein_table <- d0_sig
  
  final_d0_protein_table$ENTREZID <- uniprot_mapping$ENTREZID[
    match(
      final_d0_protein_table$UniProt,
      uniprot_mapping$UNIPROT
    )
  ]
  
  final_d0_protein_table$SYMBOL <- uniprot_mapping$SYMBOL[
    match(
      final_d0_protein_table$UniProt,
      uniprot_mapping$UNIPROT
    )
  ]
  
  final_d0_protein_table$GENENAME <- uniprot_mapping$GENENAME[
    match(
      final_d0_protein_table$UniProt,
      uniprot_mapping$UNIPROT
    )
  ]
  
  final_d0_protein_table$Direction <- ifelse(
    final_d0_protein_table$logFC > 0,
    "Upregulated",
    "Downregulated"
  )
  
  chemokine_genes <- c(
    "CXCL11",
    "CCL16",
    "CCL7",
    "TFF2",
    "CXCL10",
    "CCL24",
    "CCL8"
  )
  
  final_d0_protein_table$GO_Chemokine_Enrichment <- 
    final_d0_protein_table$SYMBOL %in% chemokine_genes
  
  final_d0_protein_table <- final_d0_protein_table[
    order(
      final_d0_protein_table$adj.P.Val,
      decreasing = FALSE
    ),
  ]
  
  cat(
    "Total significant proteins:",
    nrow(final_d0_protein_table),
    "\n"
  )
  
  cat(
    "Upregulated:",
    sum(final_d0_protein_table$Direction == "Upregulated"),
    "\n"
  )
  
  cat(
    "Downregulated:",
    sum(final_d0_protein_table$Direction == "Downregulated"),
    "\n"
  )
  
  cat(
    "Proteins contributing to GO chemokine enrichment:",
    sum(final_d0_protein_table$GO_Chemokine_Enrichment),
    "\n"
  )
  
  cat("\n=== Significant proteins ===\n")
  
  print(
    final_d0_protein_table[
      ,
      c(
        "OlinkID",
        "Assay",
        "UniProt",
        "ENTREZID",
        "SYMBOL",
        "logFC",
        "P.Value",
        "adj.P.Val",
        "Direction",
        "GO_Chemokine_Enrichment"
      )
    ]
  )
  
  cat("\n=== GO chemokine proteins ===\n")
  
  print(
    final_d0_protein_table[
      final_d0_protein_table$GO_Chemokine_Enrichment,
      c(
        "OlinkID",
        "Assay",
        "UniProt",
        "ENTREZID",
        "SYMBOL",
        "logFC",
        "adj.P.Val",
        "Direction"
      )
    ]
  )
  
  write.csv(
    final_d0_protein_table,
    file.path(
      d0_de_dir,
      "Olink_D0_Final_42_Significant_Proteins_Annotated.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    final_d0_protein_table,
    file.path(
      d0_de_dir,
      "Olink_D0_Final_42_Significant_Proteins_Annotated.rds"
    )
  )
  
  write.csv(
    final_d0_protein_table[
      final_d0_protein_table$GO_Chemokine_Enrichment,
    ],
    file.path(
      d0_de_dir,
      "Olink_D0_GO_Chemokine_Leading_Proteins.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    final_d0_protein_table[
      final_d0_protein_table$GO_Chemokine_Enrichment,
    ],
    file.path(
      d0_de_dir,
      "Olink_D0_GO_Chemokine_Leading_Proteins.rds"
    )
  )
  
  cat("\n=== Final protein annotation completed ===\n")  
  
  cat("=== D0 Proteomics: Chemokine-Enrichment Protein Plot ===\n")
  
  chemokine_plot_df <- final_d0_protein_table[
    final_d0_protein_table$GO_Chemokine_Enrichment,
    c(
      "SYMBOL",
      "Assay",
      "logFC",
      "adj.P.Val",
      "Direction"
    )
  ]
  
  chemokine_plot_df <- chemokine_plot_df[
    order(chemokine_plot_df$logFC),
  ]
  
  chemokine_plot_df$SYMBOL <- factor(
    chemokine_plot_df$SYMBOL,
    levels = chemokine_plot_df$SYMBOL
  )
  
  chemokine_plot <- ggplot(
    chemokine_plot_df,
    aes(
      x = SYMBOL,
      y = logFC
    )
  ) +
    geom_col() +
    geom_hline(
      yintercept = 0,
      linewidth = 0.5
    ) +
    labs(
      title = "D0 Proteomics: Proteins Contributing to Chemokine-Enriched GO Terms",
      x = "Protein",
      y = "logFC (COVID+ vs COVID−)"
    ) +
    theme_minimal(base_size = 12)
  
  print(chemokine_plot)
  
  ggsave(
    filename = file.path(
      d0_de_dir,
      "Olink_D0_Chemokine_Enrichment_Proteins_logFC.png"
    ),
    plot = chemokine_plot,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  write.csv(
    chemokine_plot_df,
    file.path(
      d0_de_dir,
      "Olink_D0_Chemokine_Enrichment_Proteins_Plot_Data.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    chemokine_plot,
    file.path(
      d0_de_dir,
      "Olink_D0_Chemokine_Enrichment_Proteins_logFC.rds"
    )
  )
  
  cat(
    "Plot saved:",
    file.exists(
      file.path(
        d0_de_dir,
        "Olink_D0_Chemokine_Enrichment_Proteins_logFC.png"
      )
    ),
    "\n"
  )
  
  cat(
    "Proteins plotted:",
    nrow(chemokine_plot_df),
    "\n"
  )
  
  cat("=== Chemokine-enrichment protein visualization completed ===\n")
  
  
  cat("=== MOFA2 Preparation: RNA–Proteomics Sample Matching QC ===\n")
  
  # Load RNA longitudinal cohort
  rna_cohort <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/Longitudinal/COVID_Positive_Longitudinal_Cohort_635_samples.rds"
  )
  
  # Load collapsed proteomics matrix
  protein_collapsed <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/Olink_Proteomics_Matrix_Patient_Timepoint_Collapsed.rds"
  )
  
  # RNA keys
  rna_patient <- as.character(rna_cohort$Patient_ID)
  rna_timepoint <- as.character(rna_cohort$Timepoint)
  
  rna_keys <- paste(
    rna_patient,
    rna_timepoint,
    sep = "_"
  )
  
  # Proteomics keys
  protein_keys <- rownames(protein_collapsed)
  
  # Check uniqueness
  cat(
    "RNA samples:",
    length(rna_keys),
    "\n"
  )
  
  cat(
    "RNA unique patient-timepoint keys:",
    length(unique(rna_keys)),
    "\n"
  )
  
  cat(
    "Proteomics patient-timepoint records:",
    length(protein_keys),
    "\n"
  )
  
  cat(
    "Proteomics unique patient-timepoint keys:",
    length(unique(protein_keys)),
    "\n"
  )
  
  # Common keys
  common_keys <- intersect(
    rna_keys,
    protein_keys
  )
  
  cat(
    "Common RNA–Proteomics patient-timepoints:",
    length(common_keys),
    "\n"
  )
  
  cat(
    "RNA-only patient-timepoints:",
    length(setdiff(rna_keys, protein_keys)),
    "\n"
  )
  
  cat(
    "Proteomics-only patient-timepoints:",
    length(setdiff(protein_keys, rna_keys)),
    "\n"
  )
  
  # Timepoint distribution in common dataset
  common_rna <- rna_cohort[
    rna_keys %in% common_keys,
    ,
    drop = FALSE
  ]
  
  cat("\n=== Common dataset by timepoint ===\n")
  print(table(common_rna$Timepoint))
  
  # Patient count
  common_patients <- unique(
    as.character(common_rna$Patient_ID)
  )
  
  cat(
    "\nCommon patients:",
    length(common_patients),
    "\n"
  )
  
  # Expected MOFA key preview
  cat("\n=== First common keys ===\n")
  print(head(common_keys, 20))
  
  cat("\n=== RNA–Proteomics matching QC completed ===\n")  
  
  cat("=== MOFA2 Preparation: Matching QC Details ===\n")
  
  # Identify RNA-only and Proteomics-only keys
  rna_only_keys <- setdiff(
    rna_keys,
    protein_keys
  )
  
  protein_only_keys <- setdiff(
    protein_keys,
    rna_keys
  )
  
  cat("\nRNA-only patient-timepoints:\n")
  print(rna_only_keys)
  
  cat(
    "\nNumber of RNA-only:",
    length(rna_only_keys),
    "\n"
  )
  
  cat("\nProteomics-only timepoint distribution:\n")
  
  protein_only_patient <- sub(
    "_.*$",
    "",
    protein_only_keys
  )
  
  protein_only_timepoint <- sub(
    "^[^_]+_",
    "",
    protein_only_keys
  )
  
  print(
    table(
      protein_only_timepoint
    )
  )
  
  cat(
    "\nProteomics-only unique patients:",
    length(unique(protein_only_patient)),
    "\n"
  )
  
  cat("\nCommon patient-timepoint distribution:\n")
  print(
    table(
      sub(
        "^[^_]+_",
        "",
        common_keys
      )
    )
  )
  
  # Save matching QC
  mofa_matching_qc <- data.frame(
    Metric = c(
      "RNA patient-timepoints",
      "Proteomics patient-timepoints",
      "Common patient-timepoints",
      "RNA-only patient-timepoints",
      "Proteomics-only patient-timepoints",
      "Common patients"
    ),
    Value = c(
      length(unique(rna_keys)),
      length(unique(protein_keys)),
      length(common_keys),
      length(rna_only_keys),
      length(protein_only_keys),
      length(common_patients)
    )
  )
  
  write.csv(
    mofa_matching_qc,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2_RNA_Proteomics_Matching_QC.csv",
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      Patient_Timepoint = rna_only_keys
    ),
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2_RNA_Only_Patient_Timepoints.csv",
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      Patient_Timepoint = protein_only_keys
    ),
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2_Proteomics_Only_Patient_Timepoints.csv",
    row.names = FALSE
  )
  
  cat(
    "\nMatching QC saved successfully.\n"
  )
  
  cat("=== Matching QC completed ===\n")  

  cat("=== MOFA2 Preparation: Building Matched Metadata ===\n")
  
  # Keep only RNA observations that have a matching proteomics patient-timepoint
  mofa_metadata <- rna_cohort[
    rna_keys %in% common_keys,
    ,
    drop = FALSE
  ]
  
  # Build stable MOFA key
  mofa_metadata$Patient_ID <- as.character(
    mofa_metadata$Patient_ID
  )
  
  mofa_metadata$Timepoint <- as.character(
    mofa_metadata$Timepoint
  )
  
  mofa_metadata$MOFA_Key <- paste(
    mofa_metadata$Patient_ID,
    mofa_metadata$Timepoint,
    sep = "_"
  )
  
  # Create the final metadata table using columns that actually exist
  mofa_metadata_final <- data.frame(
    MOFA_Key = mofa_metadata$MOFA_Key,
    Patient_ID = mofa_metadata$Patient_ID,
    Timepoint = mofa_metadata$Timepoint,
    Acuity_max = mofa_metadata$Acuity_max_GEO,
    Patient_category = mofa_metadata$patient_category,
    COVID_status = mofa_metadata$covid_status,
    stringsAsFactors = FALSE
  )
  
  # Sort by patient and timepoint
  mofa_metadata_final <- mofa_metadata_final[
    order(
      as.numeric(mofa_metadata_final$Patient_ID),
      factor(
        mofa_metadata_final$Timepoint,
        levels = c("D0", "D3", "D7")
      )
    ),
    ,
    drop = FALSE
  ]
  
  rownames(mofa_metadata_final) <- mofa_metadata_final$MOFA_Key
  
  cat(
    "MOFA matched observations:",
    nrow(mofa_metadata_final),
    "\n"
  )
  
  cat(
    "Unique MOFA keys:",
    length(unique(mofa_metadata_final$MOFA_Key)),
    "\n"
  )
  
  cat(
    "Unique patients:",
    length(unique(mofa_metadata_final$Patient_ID)),
    "\n"
  )
  
  cat("\n=== Timepoint distribution ===\n")
  print(
    table(
      factor(
        mofa_metadata_final$Timepoint,
        levels = c("D0", "D3", "D7")
      )
    )
  )
  
  cat("\n=== COVID status ===\n")
  print(
    table(
      mofa_metadata_final$COVID_status,
      useNA = "ifany"
    )
  )
  
  cat("\n=== Acuity max ===\n")
  print(
    table(
      mofa_metadata_final$Acuity_max,
      useNA = "ifany"
    )
  )
  
  cat(
    "\nMissing MOFA keys:",
    sum(
      is.na(mofa_metadata_final$MOFA_Key) |
        mofa_metadata_final$MOFA_Key == ""
    ),
    "\n"
  )
  
  cat(
    "Duplicated MOFA keys:",
    sum(
      duplicated(mofa_metadata_final$MOFA_Key)
    ),
    "\n"
  )
  
  # Save
  mofa_metadata_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  if (!dir.exists(mofa_metadata_dir)) {
    dir.create(
      mofa_metadata_dir,
      recursive = TRUE
    )
  }
  
  write.csv(
    mofa_metadata_final,
    file.path(
      mofa_metadata_dir,
      "MOFA2_Matched_Metadata_631_Observations.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    mofa_metadata_final,
    file.path(
      mofa_metadata_dir,
      "MOFA2_Matched_Metadata_631_Observations.rds"
    )
  )
  
  cat(
    "\nMetadata saved successfully.\n"
  )
  
  cat("=== Matched MOFA metadata completed ===\n") 
  
  
  cat("=== MOFA2 Preparation: RNA Matrix Matching ===\n")
  
  # Load the longitudinal RNA-seq dataset
  dds_mofa <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Differential_Expression/dds_primary_filtered.rds"
  )
  
  # RNA count matrix
  rna_counts_all <- counts(
    dds_mofa,
    normalized = FALSE
  )
  
  cat(
    "RNA genes:",
    nrow(rna_counts_all),
    "\n"
  )
  
  cat(
    "RNA samples:",
    ncol(rna_counts_all),
    "\n"
  )
  
  # Build RNA sample keys
  rna_sample_keys <- colnames(rna_counts_all)
  
  cat(
    "Unique RNA sample keys:",
    length(unique(rna_sample_keys)),
    "\n"
  )
  
  # Match RNA samples to the MOFA metadata
  rna_match <- match(
    mofa_metadata_final$MOFA_Key,
    rna_sample_keys
  )
  
  cat(
    "MOFA observations matched to RNA:",
    sum(!is.na(rna_match)),
    "\n"
  )
  
  cat(
    "MOFA observations NOT matched to RNA:",
    sum(is.na(rna_match)),
    "\n"
  )
  
  # Extract matched RNA counts
  rna_counts_mofa <- rna_counts_all[
    ,
    rna_match[!is.na(rna_match)],
    drop = FALSE
  ]
  
  # Reorder metadata to exactly match the RNA matrix
  mofa_metadata_rna <- mofa_metadata_final[
    !is.na(rna_match),
    ,
    drop = FALSE
  ]
  
  # Check exact alignment
  cat(
    "RNA matrix dimensions after matching:",
    paste(dim(rna_counts_mofa), collapse = " x "),
    "\n"
  )
  
  cat(
    "Metadata rows after matching:",
    nrow(mofa_metadata_rna),
    "\n"
  )
  
  cat(
    "RNA/metadata sample order identical:",
    identical(
      colnames(rna_counts_mofa),
      mofa_metadata_rna$MOFA_Key
    ),
    "\n"
  )
  
  cat(
    "Missing RNA values:",
    sum(is.na(rna_counts_mofa)),
    "\n"
  )
  
  cat(
    "Zero-count entries:",
    sum(rna_counts_mofa == 0),
    "\n"
  )
  
  cat("\n=== First RNA sample keys ===\n")
  print(
    head(colnames(rna_counts_mofa), 10)
  )
  
  # Save matched RNA matrix
  write.csv(
    rna_counts_mofa,
    file.path(
      mofa_metadata_dir,
      "MOFA2_RNA_Counts_Matched_631_Observations.csv"
    )
  )
  
  saveRDS(
    rna_counts_mofa,
    file.path(
      mofa_metadata_dir,
      "MOFA2_RNA_Counts_Matched_631_Observations.rds"
    )
  )
  
  saveRDS(
    mofa_metadata_rna,
    file.path(
      mofa_metadata_dir,
      "MOFA2_RNA_Metadata_Matched_631_Observations.rds"
    )
  )
  
  cat("\nRNA matching and QC completed successfully.\n")
  cat("=== RNA matrix preparation completed ===\n") 
  
  cat("=== Loading DESeq2 for RNA matrix extraction ===\n")
  
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    BiocManager::install(
      "DESeq2",
      ask = FALSE,
      update = FALSE
    )
  }
  
  library(DESeq2)
  
  rna_counts_all <- counts(
    dds_mofa,
    normalized = FALSE
  )
  
  cat(
    "RNA genes:",
    nrow(rna_counts_all),
    "\n"
  )
  
  cat(
    "RNA samples:",
    ncol(rna_counts_all),
    "\n"
  )
  
  cat(
    "RNA matrix extracted successfully:",
    is.matrix(rna_counts_all),
    "\n"
  )  
  
  # Build the RNA sample key from the cohort metadata
  rna_key <- paste(
    rna_cohort$Patient_ID,
    rna_cohort$Timepoint,
    sep = "_"
  )
  
  # Match the 631 MOFA observations to the 773 RNA samples
  rna_key_match <- match(
    mofa_metadata_final$Patient_Timepoint,
    rna_key
  )
  
  cat(
    "Matched RNA observations:",
    sum(!is.na(rna_key_match)),
    "\n"
  )
  
  cat(
    "Missing RNA observations:",
    sum(is.na(rna_key_match)),
    "\n"
  )
  
  # Extract the corresponding RNA counts
  rna_counts_mofa <- rna_counts_all[, rna_key_match, drop = FALSE]
  
  # Rename columns to the common MOFA patient-timepoint keys
  colnames(rna_counts_mofa) <- mofa_metadata_final$Patient_Timepoint
  
  cat(
    "RNA MOFA matrix:",
    nrow(rna_counts_mofa),
    "genes x",
    ncol(rna_counts_mofa),
    "observations\n"
  )
  
  cat(
    "Unique RNA MOFA observation IDs:",
    length(unique(colnames(rna_counts_mofa))),
    "\n"
  )
  
  cat(
    "Any duplicated observation IDs:",
    anyDuplicated(colnames(rna_counts_mofa)) > 0,
    "\n"
  )
  
  cat("=== MOFA metadata keys ===\n")
  print(head(mofa_metadata_final$Patient_Timepoint, 10))
  
  cat("\n=== RNA cohort keys ===\n")
  print(head(rna_key, 10))
  
  cat("\n=== RNA counts column names ===\n")
  print(head(colnames(rna_counts_all), 10))
  
  cat("\n=== Lengths ===\n")
  cat("MOFA keys:", length(mofa_metadata_final$Patient_Timepoint), "\n")
  cat("RNA cohort keys:", length(rna_key), "\n")
  cat("RNA count columns:", ncol(rna_counts_all), "\n")
  
  cat("\n=== Example RNA cohort metadata ===\n")
  print(
    rna_cohort[
      1:min(10, nrow(rna_cohort)),
      c("Patient_ID", "Timepoint", "rna_sample_id", "sample_id"),
      drop = FALSE
    ]
  )
  
  cat("\n=== Data types ===\n")
  cat("MOFA key class:", class(mofa_metadata_final$Patient_Timepoint), "\n")
  cat("RNA key class:", class(rna_key), "\n")  
  
  # Match RNA samples to the 631 MOFA observations
  rna_key_match <- match(
    mofa_metadata_final$MOFA_Key,
    colnames(rna_counts_all)
  )
  
  cat(
    "Matched RNA observations:",
    sum(!is.na(rna_key_match)),
    "\n"
  )
  
  cat(
    "Missing RNA observations:",
    sum(is.na(rna_key_match)),
    "\n"
  )
  
  # Extract RNA counts for the matched 631 observations
  rna_counts_mofa <- rna_counts_all[, rna_key_match, drop = FALSE]
  
  # Use the MOFA keys as the observation names
  colnames(rna_counts_mofa) <- mofa_metadata_final$MOFA_Key
  
  cat(
    "RNA MOFA matrix:",
    nrow(rna_counts_mofa),
    "genes x",
    ncol(rna_counts_mofa),
    "observations\n"
  )
  
  cat(
    "Unique RNA MOFA observation IDs:",
    length(unique(colnames(rna_counts_mofa))),
    "\n"
  )
  
  cat(
    "Any duplicated observation IDs:",
    anyDuplicated(colnames(rna_counts_mofa)) > 0,
    "\n"
  ) 
  
  
  # Check the protein matrix structure
  cat("Protein matrix dimensions:\n")
  print(dim(protein_collapsed))
  
  cat("\nFirst protein observation IDs:\n")
  print(head(rownames(protein_collapsed), 10))
  
  cat("\nUnique protein observation IDs:\n")
  cat(length(unique(rownames(protein_collapsed))), "\n")
  
  # Match protein observations to the 631 MOFA observations
  protein_key_match <- match(
    mofa_metadata_final$MOFA_Key,
    rownames(protein_collapsed)
  )
  
  cat("\nMatched protein observations:",
      sum(!is.na(protein_key_match)),
      "\n")
  
  cat("Missing protein observations:",
      sum(is.na(protein_key_match)),
      "\n")
  
  # Extract matched protein matrix
  protein_matrix_mofa <- protein_collapsed[
    protein_key_match,
    ,
    drop = FALSE
  ]
  
  # Set observation names
  rownames(protein_matrix_mofa) <- mofa_metadata_final$MOFA_Key
  
  cat("\nProtein MOFA matrix:",
      nrow(protein_matrix_mofa),
      "observations x",
      ncol(protein_matrix_mofa),
      "proteins\n")
  
  cat("Unique protein observation IDs:",
      length(unique(rownames(protein_matrix_mofa))),
      "\n")
  
  cat("Any duplicated observation IDs:",
      anyDuplicated(rownames(protein_matrix_mofa)) > 0,
      "\n")
  
  
  # Final alignment validation
  rna_ids <- colnames(rna_counts_mofa)
  protein_ids <- rownames(protein_matrix_mofa)
  metadata_ids <- mofa_metadata_final$MOFA_Key
  
  cat("RNA IDs == Protein IDs:",
      identical(rna_ids, protein_ids),
      "\n")
  
  cat("RNA IDs == Metadata IDs:",
      identical(rna_ids, metadata_ids),
      "\n")
  
  cat("Protein IDs == Metadata IDs:",
      identical(protein_ids, metadata_ids),
      "\n")
  
  cat("\nRNA dimensions:",
      nrow(rna_counts_mofa),
      "genes x",
      ncol(rna_counts_mofa),
      "observations\n")
  
  cat("Protein dimensions:",
      nrow(protein_matrix_mofa),
      "observations x",
      ncol(protein_matrix_mofa),
      "proteins\n")
  
  cat("Metadata dimensions:",
      nrow(mofa_metadata_final),
      "observations x",
      ncol(mofa_metadata_final),
      "variables\n")
  
  # Save aligned MOFA preparation objects
  mofa_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  dir.create(
    mofa_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(
    rna_counts_mofa,
    file.path(mofa_dir, "MOFA2_RNA_Counts_Aligned_631.rds")
  )
  
  saveRDS(
    protein_matrix_mofa,
    file.path(mofa_dir, "MOFA2_Protein_Matrix_Aligned_631x1429.rds")
  )
  
  saveRDS(
    mofa_metadata_final,
    file.path(mofa_dir, "MOFA2_Metadata_Aligned_631.rds")
  )
  
  cat("\nAligned MOFA preparation objects saved successfully.\n")
  
  # Load aligned RNA matrix
  rna_counts_mofa <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_RNA_Counts_Aligned_631.rds"
  )
  
  cat("RNA matrix dimensions:",
      nrow(rna_counts_mofa),
      "genes x",
      ncol(rna_counts_mofa),
      "observations\n")
  
  cat("Total RNA counts:",
      sum(rna_counts_mofa),
      "\n")
  
  cat("Zero-count entries:",
      sum(rna_counts_mofa == 0),
      "\n")
  
  cat("Genes with zero counts in all observations:",
      sum(rowSums(rna_counts_mofa) == 0),
      "\n")
  
  cat("Genes detected in at least 10% of observations:",
      sum(rowSums(rna_counts_mofa > 0) >= ceiling(0.10 * ncol(rna_counts_mofa))),
      "\n")
  
  cat("Genes detected in at least 50% of observations:",
      sum(rowSums(rna_counts_mofa > 0) >= ceiling(0.50 * ncol(rna_counts_mofa))),
      "\n")  
  
  # RNA library-size normalization and log2-CPM transformation
  
  library(edgeR)
  
  rna_dge_mofa <- DGEList(
    counts = rna_counts_mofa
  )
  
  rna_dge_mofa <- calcNormFactors(
    rna_dge_mofa,
    method = "TMM"
  )
  
  rna_logCPM_mofa <- cpm(
    rna_dge_mofa,
    log = TRUE,
    prior.count = 1
  )
  
  cat("RNA log2-CPM dimensions:",
      nrow(rna_logCPM_mofa),
      "genes x",
      ncol(rna_logCPM_mofa),
      "observations\n")
  
  cat("Any NA values:",
      anyNA(rna_logCPM_mofa),
      "\n")
  
  cat("Any infinite values:",
      any(!is.finite(rna_logCPM_mofa)),
      "\n")
  
  gene_variance <- apply(
    rna_logCPM_mofa,
    1,
    var
  )
  
  cat("\nRNA variance summary:\n")
  print(summary(gene_variance))
  
  cat("\nGenes with variance > 0.5:",
      sum(gene_variance > 0.5),
      "\n")
  
  cat("Genes with variance > 1:",
      sum(gene_variance > 1),
      "\n")
  
  cat("Genes with variance > 2:",
      sum(gene_variance > 2),
      "\n")  
  
  # Select the most variable RNA features for MOFA2
  
  top_rna_n <- 2000
  
  rna_variance_order <- order(
    gene_variance,
    decreasing = TRUE
  )
  
  top_rna_genes <- names(gene_variance)[
    rna_variance_order[seq_len(top_rna_n)]
  ]
  
  rna_mofa <- rna_logCPM_mofa[
    top_rna_genes,
    ,
    drop = FALSE
  ]
  
  cat("Selected RNA features:",
      nrow(rna_mofa),
      "\n")
  
  cat("RNA observations:",
      ncol(rna_mofa),
      "\n")
  
  cat("Minimum selected gene variance:",
      min(gene_variance[top_rna_genes]),
      "\n")
  
  cat("Median selected gene variance:",
      median(gene_variance[top_rna_genes]),
      "\n")
  
  cat("Maximum selected gene variance:",
      max(gene_variance[top_rna_genes]),
      "\n")
  
  cat("Any NA values:",
      anyNA(rna_mofa),
      "\n")
  
  cat("Any infinite values:",
      any(!is.finite(rna_mofa)),
      "\n")
  
  # Save RNA MOFA view
  mofa_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  saveRDS(
    rna_mofa,
    file.path(mofa_dir, "MOFA2_RNA_View_Top2000_Variable_Genes.rds")
  )
  
  write.csv(
    data.frame(
      Gene = top_rna_genes,
      Variance = gene_variance[top_rna_genes]
    ),
    file.path(mofa_dir, "MOFA2_RNA_Top2000_Variable_Genes.csv"),
    row.names = FALSE
  )
  
  cat("\nRNA MOFA view saved successfully.\n") 
  
  # Load aligned protein matrix
  protein_matrix_mofa <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Protein_Matrix_Aligned_631x1429.rds"
  )
  
  cat("Protein matrix dimensions:",
      nrow(protein_matrix_mofa),
      "observations x",
      ncol(protein_matrix_mofa),
      "proteins\n")
  
  cat("Missing values:",
      sum(is.na(protein_matrix_mofa)),
      "\n")
  
  cat("Missing percentage:",
      round(
        100 * sum(is.na(protein_matrix_mofa)) /
          length(protein_matrix_mofa),
        6
      ),
      "%\n")
  
  protein_variance <- apply(
    protein_matrix_mofa,
    2,
    var,
    na.rm = TRUE
  )
  
  cat("\nProtein variance summary:\n")
  print(summary(protein_variance))
  
  cat("\nProteins with variance > 0.1:",
      sum(protein_variance > 0.1),
      "\n")
  
  cat("Proteins with variance > 0.25:",
      sum(protein_variance > 0.25),
      "\n")
  
  cat("Proteins with variance > 0.5:",
      sum(protein_variance > 0.5),
      "\n")
  
  cat("Proteins with variance > 1:",
      sum(protein_variance > 1),
      "\n")  
  
  # Select the most variable protein features for MOFA2
  
  top_protein_n <- 1000
  
  protein_variance_order <- order(
    protein_variance,
    decreasing = TRUE
  )
  
  top_proteins <- names(protein_variance)[
    protein_variance_order[seq_len(top_protein_n)]
  ]
  
  protein_mofa <- protein_matrix_mofa[
    ,
    top_proteins,
    drop = FALSE
  ]
  
  # Median imputation for the MOFA view only
  protein_imputation_values <- apply(
    protein_mofa,
    2,
    median,
    na.rm = TRUE
  )
  
  for (j in seq_len(ncol(protein_mofa))) {
    missing_idx <- is.na(protein_mofa[, j])
    
    if (any(missing_idx)) {
      protein_mofa[missing_idx, j] <-
        protein_imputation_values[j]
    }
  }
  
  cat("Selected protein features:",
      ncol(protein_mofa),
      "\n")
  
  cat("Protein observations:",
      nrow(protein_mofa),
      "\n")
  
  cat("Minimum selected protein variance:",
      min(protein_variance[top_proteins]),
      "\n")
  
  cat("Median selected protein variance:",
      median(protein_variance[top_proteins]),
      "\n")
  
  cat("Maximum selected protein variance:",
      max(protein_variance[top_proteins]),
      "\n")
  
  cat("Missing values after imputation:",
      sum(is.na(protein_mofa)),
      "\n")
  
  cat("Any infinite values:",
      any(!is.finite(as.matrix(protein_mofa))),
      "\n")
  
  # Save protein MOFA view
  mofa_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  saveRDS(
    protein_mofa,
    file.path(mofa_dir, "MOFA2_Protein_View_Top1000_Variable_Proteins.rds")
  )
  
  write.csv(
    data.frame(
      Protein = top_proteins,
      Variance = protein_variance[top_proteins]
    ),
    file.path(mofa_dir, "MOFA2_Protein_Top1000_Variable_Proteins.csv"),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      Protein = top_proteins,
      Median_Imputation_Value = protein_imputation_values[top_proteins]
    ),
    file.path(mofa_dir, "MOFA2_Protein_Imputation_Values_Top1000.csv"),
    row.names = FALSE
  )
  
  cat("\nProtein MOFA view saved successfully.\n")  
  
  
  # Standardize RNA and protein views for MOFA2
  
  rna_mofa_scaled <- t(
    scale(t(rna_mofa))
  )
  
  protein_mofa_scaled <- t(
    scale(protein_mofa)
  )
  
  cat("RNA scaled dimensions:",
      nrow(rna_mofa_scaled),
      "features x",
      ncol(rna_mofa_scaled),
      "observations\n")
  
  cat("Protein scaled dimensions:",
      nrow(protein_mofa_scaled),
      "observations x",
      ncol(protein_mofa_scaled),
      "\n")
  
  cat("\nRNA mean range:\n")
  print(range(
    rowMeans(rna_mofa_scaled)
  ))
  
  cat("RNA SD range:\n")
  print(range(
    apply(rna_mofa_scaled, 1, sd)
  ))
  
  cat("\nProtein mean range:\n")
  print(range(
    colMeans(protein_mofa_scaled)
  ))
  
  cat("Protein SD range:\n")
  print(range(
    apply(protein_mofa_scaled, 2, sd)
  ))
  
  cat("\nRNA NA:",
      sum(is.na(rna_mofa_scaled)),
      "\n")
  
  cat("Protein NA:",
      sum(is.na(protein_mofa_scaled)),
      "\n")
  
  cat("RNA infinite:",
      sum(!is.finite(rna_mofa_scaled)),
      "\n")
  
  cat("Protein infinite:",
      sum(!is.finite(protein_mofa_scaled)),
      "\n")  
  
 
  # Correct protein feature-wise standardization
  
  protein_mofa_scaled <- scale(
    protein_mofa,
    center = TRUE,
    scale = TRUE
  )
  
  protein_mofa_scaled <- as.matrix(
    protein_mofa_scaled
  )
  
  protein_feature_means <- colMeans(
    protein_mofa_scaled
  )
  
  protein_feature_sds <- apply(
    protein_mofa_scaled,
    2,
    sd
  )
  
  cat("Protein scaled dimensions:",
      nrow(protein_mofa_scaled),
      "observations x",
      ncol(protein_mofa_scaled),
      "proteins\n")
  
  cat("\nProtein feature mean range:\n")
  print(range(protein_feature_means))
  
  cat("Protein feature SD range:\n")
  print(range(protein_feature_sds))
  
  cat("\nProtein features with mean deviation > 1e-10:",
      sum(abs(protein_feature_means) > 1e-10),
      "\n")
  
  cat("Protein features with SD deviation > 1e-10:",
      sum(abs(protein_feature_sds - 1) > 1e-10),
      "\n")
  
  cat("\nProtein NA:",
      sum(is.na(protein_mofa_scaled)),
      "\n")
  
  cat("Protein infinite:",
      sum(!is.finite(protein_mofa_scaled)),
      "\n")  
  
  # Final MOFA2 input validation
  
  rna_ids_final <- colnames(rna_mofa_scaled)
  protein_ids_final <- rownames(protein_mofa_scaled)
  metadata_ids_final <- mofa_metadata_final$MOFA_Key
  
  cat("RNA IDs:", length(rna_ids_final), "\n")
  cat("Protein IDs:", length(protein_ids_final), "\n")
  cat("Metadata IDs:", length(metadata_ids_final), "\n")
  
  cat("\nRNA vs Protein IDs identical:",
      identical(rna_ids_final, protein_ids_final),
      "\n")
  
  cat("RNA vs Metadata IDs identical:",
      identical(rna_ids_final, metadata_ids_final),
      "\n")
  
  cat("Protein vs Metadata IDs identical:",
      identical(protein_ids_final, metadata_ids_final),
      "\n")
  
  cat("\nRNA dimensions:",
      nrow(rna_mofa_scaled),
      "features x",
      ncol(rna_mofa_scaled),
      "observations\n")
  
  cat("Protein dimensions:",
      nrow(protein_mofa_scaled),
      "observations x",
      ncol(protein_mofa_scaled),
      "features\n")
  
  cat("\nRNA NA:", sum(is.na(rna_mofa_scaled)), "\n")
  cat("Protein NA:", sum(is.na(protein_mofa_scaled)), "\n")
  
  cat("RNA infinite:",
      sum(!is.finite(rna_mofa_scaled)),
      "\n")
  
  cat("Protein infinite:",
      sum(!is.finite(protein_mofa_scaled)),
      "\n")
  
  cat("\nMetadata rows:", nrow(mofa_metadata_final), "\n")
  cat("Metadata duplicate IDs:",
      sum(duplicated(metadata_ids_final)),
      "\n")
  
  # Check MOFA2 availability
  
  cat("Checking MOFA2 package...\n")
  
  cat("MOFA2 installed:",
      requireNamespace("MOFA2", quietly = TRUE),
      "\n")
  
  if (requireNamespace("MOFA2", quietly = TRUE)) {
    library(MOFA2)
    cat("MOFA2 version:",
        as.character(packageVersion("MOFA2")),
        "\n")
  }  
  
  # Create MOFA2 object
  
  mofa_views <- list(
    RNA = t(rna_mofa_scaled),
    Protein = protein_mofa_scaled
  )
  
  cat("RNA view:\n")
  cat("  Observations:", nrow(mofa_views$RNA), "\n")
  cat("  Features:", ncol(mofa_views$RNA), "\n")
  
  cat("\nProtein view:\n")
  cat("  Observations:", nrow(mofa_views$Protein), "\n")
  cat("  Features:", ncol(mofa_views$Protein), "\n")
  
  cat("\nIDs identical between views:",
      identical(
        rownames(mofa_views$RNA),
        rownames(mofa_views$Protein)
      ),
      "\n")
  
  mofa_object <- create_mofa(mofa_views)
  
  cat("\nMOFA object created successfully.\n")
  print(mofa_object)  
  
  # Prepare MOFA2 views in the required orientation
  
  mofa_views <- list(
    RNA = rna_mofa_scaled,
    Protein = t(protein_mofa_scaled)
  )
  
  cat("RNA view:\n")
  cat("  Features:", nrow(mofa_views$RNA), "\n")
  cat("  Samples:", ncol(mofa_views$RNA), "\n")
  
  cat("\nProtein view:\n")
  cat("  Features:", nrow(mofa_views$Protein), "\n")
  cat("  Samples:", ncol(mofa_views$Protein), "\n")
  
  cat("\nRNA sample IDs identical to Protein:",
      identical(
        colnames(mofa_views$RNA),
        colnames(mofa_views$Protein)
      ),
      "\n")
  
  cat("\nRNA sample IDs identical to Metadata:",
      identical(
        colnames(mofa_views$RNA),
        mofa_metadata_final$MOFA_Key
      ),
      "\n")  
  # Create the MOFA2 object
  
  mofa_object <- create_mofa(mofa_views)
  
  cat("MOFA2 object created successfully.\n\n")
  
  print(mofa_object)  
  # Save the initial MOFA2 object before training
  
  mofa2_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  if (!dir.exists(mofa2_dir)) {
    dir.create(mofa2_dir, recursive = TRUE)
  }
  
  saveRDS(
    mofa_object,
    file.path(
      mofa2_dir,
      "MOFA2_Untrained_Object_RNA2000_Protein1000_631.rds"
    )
  )
  
  cat("MOFA2 untrained object saved successfully.\n")
  
  cat("File:\n",
      file.path(
        mofa2_dir,
        "MOFA2_Untrained_Object_RNA2000_Protein1000_631.rds"
      ),
      "\n")
  
  cat("Object status:\n")
  print(mofa_object)  
  
  # Configure MOFA2 model options
  
  mofa_data_options <- get_default_data_options(mofa_object)
  
  mofa_model_options <- get_default_model_options(mofa_object)
  
  mofa_train_options <- get_default_training_options(mofa_object)
  
  cat("DATA OPTIONS:\n")
  print(mofa_data_options)
  
  cat("\nMODEL OPTIONS:\n")
  print(mofa_model_options)
  
  cat("\nTRAINING OPTIONS:\n")
  print(mofa_train_options) 
  
  # Configure MOFA2 training
  
  mofa_model_options$num_factors <- 15
  mofa_model_options$ard_weights <- TRUE
  mofa_model_options$ard_factors <- FALSE
  mofa_model_options$spikeslab_weights <- FALSE
  mofa_model_options$spikeslab_factors <- FALSE
  
  mofa_train_options$maxiter <- 1000
  mofa_train_options$convergence_mode <- "fast"
  mofa_train_options$verbose <- TRUE
  mofa_train_options$seed <- 42
  mofa_train_options$stochastic <- FALSE
  mofa_train_options$gpu_mode <- FALSE
  mofa_train_options$save_interrupted <- TRUE
  mofa_train_options$outfile <- file.path(
    mofa2_dir,
    "MOFA2_RNA_Protein_631_training.hdf5"
  )
  
  cat("Configured MODEL OPTIONS:\n")
  print(mofa_model_options)
  
  cat("\nConfigured TRAINING OPTIONS:\n")
  print(mofa_train_options)  
  
  # Prepare MOFA2 object for training
  
  mofa_object <- prepare_mofa(
    object = mofa_object,
    data_options = mofa_data_options,
    model_options = mofa_model_options,
    training_options = mofa_train_options
  )
  
  cat("MOFA2 object prepared successfully.\n\n")
  print(mofa_object)  
  
  # Save the prepared MOFA2 object before training
  
  prepared_mofa_path <- file.path(
    mofa2_dir,
    "MOFA2_Prepared_RNA2000_Protein1000_631.rds"
  )
  
  saveRDS(
    mofa_object,
    prepared_mofa_path
  )
  
  cat("Prepared MOFA2 object saved successfully.\n")
  cat("File:\n", prepared_mofa_path, "\n")
  cat("\nCurrent object:\n")
  print(mofa_object)  
  
  # Train the MOFA2 model
  
  set.seed(42)
  
  mofa_object <- run_mofa(
    mofa_object,
    use_basilisk = FALSE
  )
  
  cat("\nMOFA2 training completed.\n\n")
  print(mofa_object)  
  
  # Check Basilisk availability for MOFA2
  
  cat("Checking basilisk package...\n")
  
  cat("Basilisk installed:",
      requireNamespace("basilisk", quietly = TRUE),
      "\n")
  
  if (requireNamespace("basilisk", quietly = TRUE)) {
    cat("Basilisk version:",
        as.character(packageVersion("basilisk")),
        "\n")
  } 
  # Diagnose Basilisk configuration
  
  library(basilisk)
  
  cat("R version:\n")
  print(R.version.string)
  
  cat("\nBasilisk version:\n")
  print(packageVersion("basilisk"))
  
  cat("\nBasilisk installation path:\n")
  print(system.file(package = "basilisk"))
  
  cat("\nBasilisk cache directory:\n")
  print(basilisk.utils::getBasiliskDir())
  
  cat("\nR_LIBS_USER:\n")
  print(Sys.getenv("R_LIBS_USER"))
  
  cat("\nTemporary directory:\n")
  print(tempdir()) 
  
  # Check the MOFA2 Basilisk environment
  
  cat("MOFA2 package path:\n")
  print(system.file(package = "MOFA2"))
  
  cat("\nMOFA2 Basilisk environment files:\n")
  
  mofa_files <- list.files(
    system.file(package = "MOFA2"),
    recursive = TRUE,
    full.names = TRUE
  )
  
  print(
    mofa_files[
      grepl(
        "basilisk|environment|mofapy2",
        mofa_files,
        ignore.case = TRUE
      )
    ]
  )
  
  cat("\nBasilisk package dependencies:\n")
  print(
    packageDescription("MOFA2")$Imports
  ) 
  
  cat("Basilisk default Python version:\n")
  
  default_python <- get(
    "defaultPythonVersion",
    envir = basilisk_ns
  )
  
  print(default_python)
  
  cat("\nBasilisk internal Python version:\n")
  
  python_version_internal <- get(
    ".python_version",
    envir = basilisk_ns
  )
  
  print(python_version_internal)
  
  cat("\nBasilisk listPythonVersion object:\n")
  
  list_python <- get(
    "listPythonVersion",
    envir = basilisk_ns
  )
  
  print(list_python)
  
  cat("\nObject classes:\n")
  
  cat("defaultPythonVersion: ")
  print(class(default_python))
  
  cat(".python_version: ")
  print(class(python_version_internal))
  
  cat("listPythonVersion: ")
  print(class(list_python)) 
  
  cat("MOFA object class:\n")
  print(class(mofa_object))
  
  cat("\nTraining file exists:\n")
  print(file.exists(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_RNA_Protein_631_training.hdf5"
  ))
  
  cat("\nNumber of factors:\n")
  print(mofa_object@dimensions$K) 
  
  trained_mofa <- MOFA2::load_model(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_RNA_Protein_631_training.hdf5"
  )
  
  cat("Class:\n")
  print(class(trained_mofa))
  
  cat("\nNumber of factors:\n")
  print(trained_mofa@dimensions$K)
  
  cat("\nNumber of samples:\n")
  print(trained_mofa@dimensions$N)
  
  cat("\nNumber of views:\n")
  print(trained_mofa@dimensions$M) 
  library(MOFA2)
  
  mofa_object <- trained_mofa
  
  variance_explained <- MOFA2::plot_variance_explained(
    mofa_object,
    x = "view",
    y = "factor",
    plot_total = TRUE,
    return_data = TRUE
  )
  
  print(variance_explained)  
  
  
  variance_data <- MOFA2::calculate_variance_explained(
    mofa_object
  )
  
  print(variance_data)
  
  library(ggplot2)
  
  r2_plot_data <- reshape2::melt(
    r2_factor,
    id.vars = "Factor",
    variable.name = "View",
    value.name = "Variance_Explained"
  )
  
  r2_plot_data$Factor <- factor(
    r2_plot_data$Factor,
    levels = paste0("Factor", 1:15)
  )
  
  p_mofa_variance <- ggplot(
    r2_plot_data,
    aes(x = Factor, y = Variance_Explained, fill = View)
  ) +
    geom_col(position = "dodge") +
    labs(
      title = "MOFA2 Variance Explained by Factor and View",
      x = "MOFA2 factor",
      y = "Variance explained (%)",
      fill = "Omics view"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(p_mofa_variance)
  
  ggsave(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Variance_Explained_Per_Factor.png",
    p_mofa_variance,
    width = 11,
    height = 6,
    dpi = 300
  )
  
  write.csv(
    r2_plot_data,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Variance_Explained_Per_Factor_PlotData.csv",
    row.names = FALSE
  )
  
  saveRDS(
    p_mofa_variance,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Variance_Explained_Per_Factor_Plot.rds"
  ) 
  
  factor_scores <- MOFA2::get_factors(
    mofa_object,
    factors = "all",
    groups = "all",
    as.data.frame = TRUE
  )
  
  cat("Factor score data:\n")
  print(dim(factor_scores))
  print(head(factor_scores))
  
  cat("\nColumns:\n")
  print(names(factor_scores))  
  
  
  factor_scores$MOFA_Key <- factor_scores$sample
  
  mofa_meta_scores <- mofa_metadata_final
  
  factor_scores_annotated <- merge(
    factor_scores,
    mofa_meta_scores,
    by.x = "MOFA_Key",
    by.y = "MOFA_Key",
    all.x = TRUE,
    sort = FALSE
  )
  
  cat("Annotated factor-score dimensions:\n")
  print(dim(factor_scores_annotated))
  
  cat("\nMissing Timepoint:\n")
  print(sum(is.na(factor_scores_annotated$Timepoint)))
  
  cat("\nMissing Acuity:\n")
  print(sum(is.na(factor_scores_annotated$Acuity_max)))
  
  cat("\nMissing Patient ID:\n")
  print(sum(is.na(factor_scores_annotated$Patient_ID)))
  
  cat("\nFactor counts:\n")
  print(table(factor_scores_annotated$factor))
  
  cat("\nTimepoint counts:\n")
  print(table(factor_scores_annotated$Timepoint)) 
  
  
  mofa_factor_scores_annotated <- factor_scores_annotated
  
  saveRDS(
    mofa_factor_scores_annotated,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_Annotated.rds"
  )
  
  write.csv(
    mofa_factor_scores_annotated,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_Annotated.csv",
    row.names = FALSE
  )
  
  cat("Saved MOFA2 annotated factor scores.\n")
  cat("Rows:", nrow(mofa_factor_scores_annotated), "\n")
  cat("Columns:", ncol(mofa_factor_scores_annotated), "\n") 
  
  mofa_timepoint_summary <- aggregate(
    value ~ factor + Timepoint,
    data = mofa_factor_scores_annotated,
    FUN = mean
  )
  
  mofa_timepoint_sd <- aggregate(
    value ~ factor + Timepoint,
    data = mofa_factor_scores_annotated,
    FUN = sd
  )
  
  mofa_timepoint_n <- aggregate(
    value ~ factor + Timepoint,
    data = mofa_factor_scores_annotated,
    FUN = length
  )
  
  names(mofa_timepoint_summary)[3] <- "Mean"
  names(mofa_timepoint_sd)[3] <- "SD"
  names(mofa_timepoint_n)[3] <- "N"
  
  mofa_timepoint_summary$SD <- mofa_timepoint_sd$SD
  mofa_timepoint_summary$N <- mofa_timepoint_n$N
  
  mofa_timepoint_summary$SE <- 
    mofa_timepoint_summary$SD / sqrt(mofa_timepoint_summary$N)
  
  mofa_timepoint_summary$Timepoint <- factor(
    mofa_timepoint_summary$Timepoint,
    levels = c("D0", "D3", "D7")
  )
  
  cat("Timepoint summary dimensions:\n")
  print(dim(mofa_timepoint_summary))
  
  cat("\nFirst rows:\n")
  print(head(mofa_timepoint_summary, 12))
  
  cat("\nObservations per factor/timepoint:\n")
  print(table(
    mofa_timepoint_summary$factor,
    mofa_timepoint_summary$Timepoint
  ))  
  
  n_check <- xtabs(
    N ~ factor + Timepoint,
    data = mofa_timepoint_summary
  )
  
  cat("N per factor/timepoint:\n")
  print(n_check)
  
  cat("\nUnique N values:\n")
  print(unique(mofa_timepoint_summary$N))
  
  cat("\nExpected sample counts:\n")
  cat("D0 =", length(unique(mofa_factor_scores_annotated$sample[
    mofa_factor_scores_annotated$Timepoint == "D0"
  ])), "\n")
  
  cat("D3 =", length(unique(mofa_factor_scores_annotated$sample[
    mofa_factor_scores_annotated$Timepoint == "D3"
  ])), "\n")
  
  cat("D7 =", length(unique(mofa_factor_scores_annotated$sample[
    mofa_factor_scores_annotated$Timepoint == "D7"
  ])), "\n")
  
  saveRDS(
    mofa_timepoint_summary,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_vs_Timepoint_Summary.rds"
  )
  
  write.csv(
    mofa_timepoint_summary,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_vs_Timepoint_Summary.csv",
    row.names = FALSE
  )
  
  cat("\nTimepoint summary saved successfully.\n")  

  
  library(ggplot2)
  
  mofa_timepoint_plot <- ggplot(
    mofa_timepoint_summary,
    aes(x = Timepoint, y = Mean, group = 1)
  ) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(
      aes(
        ymin = Mean - SE,
        ymax = Mean + SE
      ),
      width = 0.12
    ) +
    facet_wrap(
      ~ factor,
      scales = "free_y",
      ncol = 3
    ) +
    labs(
      title = "MOFA2 Factor Scores Across Longitudinal Timepoints",
      x = "Timepoint",
      y = "Mean Factor Score ± SE"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
  
  print(mofa_timepoint_plot)
  
  ggsave(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_vs_Timepoint.png",
    mofa_timepoint_plot,
    width = 12,
    height = 15,
    dpi = 300
  )
  
  write.csv(
    mofa_timepoint_summary,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_vs_Timepoint_PlotData.csv",
    row.names = FALSE
  )
  
  saveRDS(
    mofa_timepoint_plot,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_vs_Timepoint_Plot.rds"
  )
  
  cat("MOFA2 Factor Scores vs Timepoint plot saved successfully.\n")   
  
  mofa_acuity_summary <- aggregate(
    value ~ factor + Acuity_max,
    data = mofa_factor_scores_annotated,
    FUN = mean
  )
  
  mofa_acuity_sd <- aggregate(
    value ~ factor + Acuity_max,
    data = mofa_factor_scores_annotated,
    FUN = sd
  )
  
  mofa_acuity_n <- aggregate(
    value ~ factor + Acuity_max,
    data = mofa_factor_scores_annotated,
    FUN = length
  )
  
  names(mofa_acuity_summary)[3] <- "Mean"
  names(mofa_acuity_sd)[3] <- "SD"
  names(mofa_acuity_n)[3] <- "N"
  
  mofa_acuity_summary$SD <- mofa_acuity_sd$SD
  mofa_acuity_summary$N <- mofa_acuity_n$N
  
  mofa_acuity_summary$SE <-
    mofa_acuity_summary$SD / sqrt(mofa_acuity_summary$N)
  
  mofa_acuity_summary$Acuity_max <- factor(
    mofa_acuity_summary$Acuity_max,
    levels = sort(unique(mofa_acuity_summary$Acuity_max))
  )
  
  cat("Acuity summary dimensions:\n")
  print(dim(mofa_acuity_summary))
  
  cat("\nFirst rows:\n")
  print(head(mofa_acuity_summary, 15))
  
  cat("\nN per factor/acutity:\n")
  print(xtabs(N ~ factor + Acuity_max, data = mofa_acuity_summary))
  
  saveRDS(
    mofa_acuity_summary,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_vs_Acuity_Summary.rds"
  )
  
  write.csv(
    mofa_acuity_summary,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_vs_Acuity_Summary.csv",
    row.names = FALSE
  )
  
  cat("\nAcuity summary saved successfully.\n")
  
  library(ggplot2)
  
  mofa_acuity_plot <- ggplot(
    mofa_acuity_summary,
    aes(x = Acuity_max, y = Mean, group = 1)
  ) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(
      aes(
        ymin = Mean - SE,
        ymax = Mean + SE
      ),
      width = 0.12
    ) +
    facet_wrap(
      ~ factor,
      scales = "free_y",
      ncol = 3
    ) +
    labs(
      title = "MOFA2 Factor Scores Across Disease Acuity",
      x = "Acuity",
      y = "Mean Factor Score ± SE"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
  
  print(mofa_acuity_plot)
  
  ggsave(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_vs_Acuity.png",
    mofa_acuity_plot,
    width = 12,
    height = 15,
    dpi = 300
  )
  
  write.csv(
    mofa_acuity_summary,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_vs_Acuity_PlotData.csv",
    row.names = FALSE
  )
  
  saveRDS(
    mofa_acuity_plot,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Scores_vs_Acuity_Plot.rds"
  )
  
  cat("MOFA2 Factor Scores vs Acuity plot saved successfully.\n")
  
  trajectory_data <- mofa_factor_scores_annotated[
    mofa_factor_scores_annotated$Patient_ID %in%
      unique(
        mofa_factor_scores_annotated$Patient_ID[
          ave(
            mofa_factor_scores_annotated$Timepoint,
            mofa_factor_scores_annotated$Patient_ID,
            FUN = function(x) length(unique(x))
          ) >= 2
        ]
      ),
  ]
  
  trajectory_data$Timepoint <- factor(
    trajectory_data$Timepoint,
    levels = c("D0", "D3", "D7")
  )
  
  patient_timepoint_counts <- aggregate(
    Timepoint ~ Patient_ID,
    data = trajectory_data,
    FUN = function(x) length(unique(x))
  )
  
  cat("Trajectory dataset dimensions:\n")
  print(dim(trajectory_data))
  
  cat("\nUnique patients:\n")
  print(length(unique(trajectory_data$Patient_ID)))
  
  cat("\nPatients by number of timepoints:\n")
  print(table(patient_timepoint_counts$Timepoint))
  
  cat("\nSamples by timepoint:\n")
  print(table(trajectory_data$Timepoint))
  
  cat("\nFactors:\n")
  print(table(trajectory_data$factor))
  
  cat("\nDuplicate Patient × Timepoint × Factor records:\n")
  print(
    sum(
      duplicated(
        trajectory_data[
          , c("Patient_ID", "Timepoint", "factor")
        ]
      )
    )
  )
  
  saveRDS(
    trajectory_data,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Trajectory_Data.rds"
  )
  
  write.csv(
    trajectory_data,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Trajectory_Data.csv",
    row.names = FALSE
  )
  
  cat("\nTrajectory data saved successfully.\n")  
  library(ggplot2)
  
  mofa_trajectory_plot <- ggplot(
    trajectory_data,
    aes(
      x = Timepoint,
      y = value,
      group = Patient_ID
    )
  ) +
    geom_line(
      alpha = 0.12,
      linewidth = 0.3
    ) +
    geom_point(
      alpha = 0.12,
      size = 0.5
    ) +
    stat_summary(
      aes(group = 1),
      fun = mean,
      geom = "line",
      linewidth = 1
    ) +
    stat_summary(
      aes(group = 1),
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.12
    ) +
    stat_summary(
      aes(group = 1),
      fun = mean,
      geom = "point",
      size = 2
    ) +
    facet_wrap(
      ~ factor,
      scales = "free_y",
      ncol = 3
    ) +
    labs(
      title = "Within-Patient MOFA2 Factor Trajectories",
      subtitle = "Individual trajectories with mean ± SE",
      x = "Timepoint",
      y = "MOFA2 Factor Score"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
  
  print(mofa_trajectory_plot)
  
  ggsave(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Within_Patient_Factor_Trajectories.png",
    mofa_trajectory_plot,
    width = 12,
    height = 15,
    dpi = 300
  )
  
  write.csv(
    trajectory_data,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Within_Patient_Factor_Trajectories_PlotData.csv",
    row.names = FALSE
  )
  
  saveRDS(
    mofa_trajectory_plot,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Within_Patient_Factor_Trajectories_Plot.rds"
  )
  
  cat("Within-patient MOFA2 trajectory plot saved successfully.\n")  
  
  library(dplyr)
  
  paired_factor_test <- function(data, factor_name, tp1, tp2) {
    
    x <- data %>%
      filter(
        factor == factor_name,
        Timepoint %in% c(tp1, tp2)
      ) %>%
      select(Patient_ID, Timepoint, value) %>%
      tidyr::pivot_wider(
        names_from = Timepoint,
        values_from = value
      ) %>%
      filter(
        !is.na(.data[[tp1]]),
        !is.na(.data[[tp2]])
      )
    
    if (nrow(x) < 3) {
      return(
        data.frame(
          factor = factor_name,
          Timepoint_1 = tp1,
          Timepoint_2 = tp2,
          N = nrow(x),
          Mean_Difference = NA_real_,
          SD_Difference = NA_real_,
          t_statistic = NA_real_,
          p_value = NA_real_
        )
      )
    }
    
    d <- x[[tp2]] - x[[tp1]]
    test <- t.test(
      x[[tp2]],
      x[[tp1]],
      paired = TRUE
    )
    
    data.frame(
      factor = factor_name,
      Timepoint_1 = tp1,
      Timepoint_2 = tp2,
      N = length(d),
      Mean_Difference = mean(d),
      SD_Difference = sd(d),
      t_statistic = unname(test$statistic),
      p_value = test$p.value
    )
  }
  
  factor_names <- unique(as.character(trajectory_data$factor))
  
  trajectory_stats <- bind_rows(
    lapply(
      factor_names,
      function(f) paired_factor_test(
        trajectory_data,
        f,
        "D0",
        "D3"
      )
    ),
    lapply(
      factor_names,
      function(f) paired_factor_test(
        trajectory_data,
        f,
        "D0",
        "D7"
      )
    ),
    lapply(
      factor_names,
      function(f) paired_factor_test(
        trajectory_data,
        f,
        "D3",
        "D7"
      )
    )
  )
  
  trajectory_stats$FDR <- NA_real_
  
  for (comparison in unique(
    paste(
      trajectory_stats$Timepoint_1,
      trajectory_stats$Timepoint_2,
      sep = "_"
    )
  )) {
    
    idx <- paste(
      trajectory_stats$Timepoint_1,
      trajectory_stats$Timepoint_2,
      sep = "_"
    ) == comparison
    
    trajectory_stats$FDR[idx] <- p.adjust(
      trajectory_stats$p_value[idx],
      method = "BH"
    )
  }
  
  trajectory_stats$Significant_FDR05 <- trajectory_stats$FDR < 0.05
  
  trajectory_stats$Direction <- ifelse(
    trajectory_stats$Mean_Difference > 0,
    "Increased",
    ifelse(
      trajectory_stats$Mean_Difference < 0,
      "Decreased",
      "No change"
    )
  )
  
  cat("Longitudinal factor statistics:\n")
  print(dim(trajectory_stats))
  
  cat("\nD3 vs D0:\n")
  print(
    trajectory_stats %>%
      filter(
        Timepoint_1 == "D0",
        Timepoint_2 == "D3"
      ) %>%
      arrange(FDR)
  )
  
  cat("\nD7 vs D0:\n")
  print(
    trajectory_stats %>%
      filter(
        Timepoint_1 == "D0",
        Timepoint_2 == "D7"
      ) %>%
      arrange(FDR)
  )
  
  cat("\nD7 vs D3:\n")
  print(
    trajectory_stats %>%
      filter(
        Timepoint_1 == "D3",
        Timepoint_2 == "D7"
      ) %>%
      arrange(FDR)
  )
  
  saveRDS(
    trajectory_stats,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Longitudinal_Factor_Statistics.rds"
  )
  
  write.csv(
    trajectory_stats,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Longitudinal_Factor_Statistics.csv",
    row.names = FALSE
  )
  
  cat("\nLongitudinal factor statistics saved successfully.\n")  
  
  library(dplyr)
  
  factor_clinical_correlation <- bind_rows(
    
    lapply(
      c("Overall", "D0", "D3", "D7"),
      function(tp) {
        
        dat <- if (tp == "Overall") {
          mofa_factor_scores_annotated
        } else {
          mofa_factor_scores_annotated %>%
            filter(Timepoint == tp)
        }
        
        results <- lapply(
          unique(as.character(dat$factor)),
          function(f) {
            
            x <- dat %>%
              filter(factor == f) %>%
              select(value, Acuity_max) %>%
              filter(
                !is.na(value),
                !is.na(Acuity_max)
              )
            
            test <- suppressWarnings(
              cor.test(
                x$value,
                as.numeric(as.character(x$Acuity_max)),
                method = "spearman",
                exact = FALSE
              )
            )
            
            data.frame(
              Timepoint = tp,
              factor = f,
              N = nrow(x),
              Spearman_rho = unname(test$estimate),
              p_value = test$p.value
            )
          }
        )
        
        bind_rows(results)
      }
    )
  )
  
  factor_clinical_correlation <- factor_clinical_correlation %>%
    group_by(Timepoint) %>%
    mutate(
      FDR = p.adjust(p_value, method = "BH"),
      Significant_FDR05 = FDR < 0.05,
      Direction = case_when(
        Spearman_rho > 0 ~ "Positive",
        Spearman_rho < 0 ~ "Negative",
        TRUE ~ "No correlation"
      )
    ) %>%
    ungroup()
  
  cat("Factor-clinical correlation dimensions:\n")
  print(dim(factor_clinical_correlation))
  
  cat("\nOverall correlations:\n")
  print(
    factor_clinical_correlation %>%
      filter(Timepoint == "Overall") %>%
      arrange(FDR)
  )
  
  cat("\nD0 correlations:\n")
  print(
    factor_clinical_correlation %>%
      filter(Timepoint == "D0") %>%
      arrange(FDR)
  )
  
  cat("\nD3 correlations:\n")
  print(
    factor_clinical_correlation %>%
      filter(Timepoint == "D3") %>%
      arrange(FDR)
  )
  
  cat("\nD7 correlations:\n")
  print(
    factor_clinical_correlation %>%
      filter(Timepoint == "D7") %>%
      arrange(FDR)
  )
  
  cat("\nSignificant correlations by timepoint:\n")
  print(
    factor_clinical_correlation %>%
      group_by(Timepoint) %>%
      summarise(
        Significant = sum(Significant_FDR05),
        Total = n(),
        .groups = "drop"
      )
  )
  
  saveRDS(
    factor_clinical_correlation,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Acuity_Spearman_Correlations.rds"
  )
  
  write.csv(
    factor_clinical_correlation,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Factor_Acuity_Spearman_Correlations.csv",
    row.names = FALSE
  )
  
  cat("\nFactor-clinical correlation results saved successfully.\n")  

  cat("Extracting MOFA2 RNA feature weights...\n")
  
  rna_weights <- MOFA2::get_weights(
    mofa_object,
    views = "RNA",
    factors = "all",
    scale = FALSE,
    as.data.frame = TRUE
  )
  
  cat("\nRNA weights dimensions:\n")
  print(dim(rna_weights))
  
  cat("\nRNA weights columns:\n")
  print(names(rna_weights))
  
  cat("\nFirst rows:\n")
  print(head(rna_weights))
  
  cat("\nNumber of unique features:\n")
  print(length(unique(rna_weights$feature)))
  
  cat("\nNumber of factors:\n")
  print(length(unique(rna_weights$factor)))
  
  cat("\nMissing weights:\n")
  print(sum(is.na(rna_weights$value)))
  
  cat("\nRNA feature weights extracted successfully.\n")  
  
  protein_weights <- MOFA2::get_weights(
    mofa_object,
    views = "Protein",
    factors = "all",
    as.data.frame = TRUE
  )
  
  cat("Protein weights dimensions:\n")
  print(dim(protein_weights))
  
  cat("\nProtein weights columns:\n")
  print(names(protein_weights))
  
  cat("\nFirst rows:\n")
  print(head(protein_weights))
  
  cat("\nNumber of unique features:\n")
  print(length(unique(protein_weights$feature)))
  
  cat("\nNumber of factors:\n")
  print(length(unique(protein_weights$factor)))
  
  cat("\nMissing weights:\n")
  print(sum(is.na(protein_weights$value)))
  
  saveRDS(
    protein_weights,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Protein_Feature_Weights.rds"
  )
  
  write.csv(
    protein_weights,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Protein_Feature_Weights.csv",
    row.names = FALSE
  )
  
  cat("\nProtein feature weights extracted and saved successfully.\n") 
  
  library(dplyr)
  
  rna_top_features <- rna_weights %>%
    mutate(
      abs_weight = abs(value)
    ) %>%
    group_by(factor) %>%
    arrange(desc(abs_weight), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    mutate(
      rank = row_number()
    ) %>%
    ungroup() %>%
    select(
      factor,
      rank,
      feature,
      value,
      abs_weight,
      view
    )
  
  cat("Top RNA features dimensions:\n")
  print(dim(rna_top_features))
  
  cat("\nTop RNA features per factor:\n")
  print(
    rna_top_features %>%
      group_by(factor) %>%
      summarise(
        n_features = n(),
        top_feature = feature[1],
        top_abs_weight = abs_weight[1],
        .groups = "drop"
      )
  )
  
  cat("\nFirst 30 RNA top features:\n")
  print(head(rna_top_features, 30))
  
  saveRDS(
    rna_top_features,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_RNA_Features_Per_Factor.rds"
  )
  
  write.csv(
    rna_top_features,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_RNA_Features_Per_Factor.csv",
    row.names = FALSE
  )
  
  cat("\nTop RNA features saved successfully.\n")  
  
  protein_top_features <- protein_weights %>%
    mutate(
      abs_weight = abs(value)
    ) %>%
    group_by(factor) %>%
    arrange(desc(abs_weight), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    mutate(
      rank = row_number()
    ) %>%
    ungroup() %>%
    select(
      factor,
      rank,
      feature,
      value,
      abs_weight,
      view
    )
  
  cat("Top Protein features dimensions:\n")
  print(dim(protein_top_features))
  
  cat("\nTop Protein features per factor:\n")
  print(
    protein_top_features %>%
      group_by(factor) %>%
      summarise(
        n_features = n(),
        top_feature = feature[1],
        top_abs_weight = abs_weight[1],
        .groups = "drop"
      )
  )
  
  cat("\nFirst 30 Protein top features:\n")
  print(head(protein_top_features, 30))
  
  saveRDS(
    protein_top_features,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_Protein_Features_Per_Factor.rds"
  )
  
  write.csv(
    protein_top_features,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_Protein_Features_Per_Factor.csv",
    row.names = FALSE
  )
  
  cat("\nTop Protein features saved successfully.\n")  
  
  library(dplyr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(readr)
  
  cat("Annotating top RNA features...\n")
  
  rna_top_annotated <- rna_top_features %>%
    mutate(
      feature = as.character(feature)
    )
  
  rna_gene_symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique(rna_top_annotated$feature),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  rna_top_annotated <- rna_top_annotated %>%
    mutate(
      gene_symbol = unname(rna_gene_symbols[feature]),
      gene_symbol = ifelse(
        is.na(gene_symbol) | gene_symbol == "",
        feature,
        gene_symbol
      )
    ) %>%
    select(
      factor,
      rank,
      feature,
      gene_symbol,
      value,
      abs_weight,
      view
    )
  
  cat("\nRNA annotation summary:\n")
  cat("Total RNA features:", nrow(rna_top_annotated), "\n")
  cat(
    "Mapped to gene symbol:",
    sum(rna_top_annotated$gene_symbol != rna_top_annotated$feature),
    "\n"
  )
  cat(
    "Unique gene symbols:",
    length(unique(rna_top_annotated$gene_symbol)),
    "\n"
  )
  
  cat("\nFirst 20 annotated RNA features:\n")
  print(head(rna_top_annotated, 20))
  
  
  cat("\n\nAnnotating top Protein features...\n")
  
  protein_annotation <- read.csv(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/Olink_Proteomics_Assay_Annotation.csv",
    stringsAsFactors = FALSE
  )
  
  protein_top_annotated <- protein_top_features %>%
    mutate(
      feature = as.character(feature)
    ) %>%
    left_join(
      protein_annotation %>%
        select(
          OlinkID,
          Assay,
          UniProt,
          Panel,
          Panel_Version
        ),
      by = c("feature" = "OlinkID")
    ) %>%
    mutate(
      Assay = ifelse(
        is.na(Assay) | Assay == "",
        feature,
        Assay
      )
    ) %>%
    select(
      factor,
      rank,
      feature,
      Assay,
      UniProt,
      Panel,
      Panel_Version,
      value,
      abs_weight,
      view
    )
  
  cat("\nProtein annotation summary:\n")
  cat("Total Protein features:", nrow(protein_top_annotated), "\n")
  cat(
    "Mapped to Assay:",
    sum(protein_top_annotated$Assay != protein_top_annotated$feature),
    "\n"
  )
  cat(
    "Unique Assays:",
    length(unique(protein_top_annotated$Assay)),
    "\n"
  )
  cat(
    "Missing UniProt:",
    sum(is.na(protein_top_annotated$UniProt) | protein_top_annotated$UniProt == ""),
    "\n"
  )
  
  cat("\nFirst 20 annotated Protein features:\n")
  print(head(protein_top_annotated, 20))
  
  
  saveRDS(
    rna_top_annotated,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_RNA_Features_Per_Factor_Annotated.rds"
  )
  
  write.csv(
    rna_top_annotated,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_RNA_Features_Per_Factor_Annotated.csv",
    row.names = FALSE
  )
  
  saveRDS(
    protein_top_annotated,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_Protein_Features_Per_Factor_Annotated.rds"
  )
  
  write.csv(
    protein_top_annotated,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_Protein_Features_Per_Factor_Annotated.csv",
    row.names = FALSE
  )
  
  cat("\nAnnotated RNA and Protein feature tables saved successfully.\n")
  
  cat("Re-running feature annotation with explicit dplyr::select...\n")
  
  rna_top_annotated <- rna_top_features %>%
    dplyr::mutate(
      feature = as.character(feature)
    )
  
  rna_gene_symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique(rna_top_annotated$feature),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  rna_top_annotated <- rna_top_annotated %>%
    dplyr::mutate(
      gene_symbol = unname(rna_gene_symbols[feature]),
      gene_symbol = ifelse(
        is.na(gene_symbol) | gene_symbol == "",
        feature,
        gene_symbol
      )
    ) %>%
    dplyr::select(
      factor,
      rank,
      feature,
      gene_symbol,
      value,
      abs_weight,
      view
    )
  
  cat("\nRNA annotation:\n")
  cat("Total features:", nrow(rna_top_annotated), "\n")
  cat(
    "Mapped gene symbols:",
    sum(rna_top_annotated$gene_symbol != rna_top_annotated$feature),
    "\n"
  )
  
  cat("\nFirst 20 annotated RNA features:\n")
  print(head(rna_top_annotated, 20))
  
  
  protein_annotation <- read.csv(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/Olink_Proteomics_Assay_Annotation.csv",
    stringsAsFactors = FALSE
  )
  
  protein_top_annotated <- protein_top_features %>%
    dplyr::mutate(
      feature = as.character(feature)
    ) %>%
    dplyr::left_join(
      protein_annotation %>%
        dplyr::select(
          OlinkID,
          Assay,
          UniProt,
          Panel,
          Panel_Version
        ),
      by = c("feature" = "OlinkID")
    ) %>%
    dplyr::mutate(
      Assay = ifelse(
        is.na(Assay) | Assay == "",
        feature,
        Assay
      )
    ) %>%
    dplyr::select(
      factor,
      rank,
      feature,
      Assay,
      UniProt,
      Panel,
      Panel_Version,
      value,
      abs_weight,
      view
    )
  
  cat("\nProtein annotation:\n")
  cat("Total features:", nrow(protein_top_annotated), "\n")
  cat(
    "Mapped to Assay:",
    sum(protein_top_annotated$Assay != protein_top_annotated$feature),
    "\n"
  )
  cat(
    "Missing UniProt:",
    sum(
      is.na(protein_top_annotated$UniProt) |
        protein_top_annotated$UniProt == ""
    ),
    "\n"
  )
  
  cat("\nFirst 20 annotated Protein features:\n")
  print(head(protein_top_annotated, 20))
  
  
  saveRDS(
    rna_top_annotated,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_RNA_Features_Per_Factor_Annotated.rds"
  )
  
  write.csv(
    rna_top_annotated,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_RNA_Features_Per_Factor_Annotated.csv",
    row.names = FALSE
  )
  
  saveRDS(
    protein_top_annotated,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_Protein_Features_Per_Factor_Annotated.rds"
  )
  
  write.csv(
    protein_top_annotated,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Top20_Protein_Features_Per_Factor_Annotated.csv",
    row.names = FALSE
  )
  
  cat("\nAnnotated RNA and Protein feature tables saved successfully.\n")  
  
  
  library(ggplot2)
  library(dplyr)
  
  mofa_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  rna_top_annotated <- readRDS(
    file.path(
      mofa_dir,
      "MOFA2_Top20_RNA_Features_Per_Factor_Annotated.rds"
    )
  )
  
  rna_heatmap_data <- rna_top_annotated %>%
    dplyr::mutate(
      feature_label = ifelse(
        is.na(gene_symbol) | gene_symbol == "",
        feature,
        gene_symbol
      ),
      factor = factor(
        factor,
        levels = paste0("Factor", 1:15)
      )
    ) %>%
    dplyr::arrange(factor, rank)
  
  cat("\nRNA heatmap data:\n")
  cat("Rows:", nrow(rna_heatmap_data), "\n")
  cat("Factors:", dplyr::n_distinct(rna_heatmap_data$factor), "\n")
  cat("Features per factor:\n")
  print(table(rna_heatmap_data$factor))
  
  p_rna_weights <- ggplot(
    rna_heatmap_data,
    aes(
      x = factor,
      y = reorder(feature_label, abs_weight),
      fill = value
    )
  ) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient2(
      midpoint = 0,
      name = "MOFA weight"
    ) +
    labs(
      title = "Top 20 RNA Features per MOFA2 Factor",
      subtitle = "Features ranked by absolute feature weight",
      x = "MOFA2 Factor",
      y = "RNA feature"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(
        size = 7
      ),
      panel.grid = element_blank()
    )
  
  print(p_rna_weights)
  
  ggsave(
    file.path(
      mofa_dir,
      "MOFA2_RNA_Feature_Weights_Top20_Per_Factor_Heatmap.png"
    ),
    p_rna_weights,
    width = 12,
    height = 14,
    dpi = 300
  )
  
  saveRDS(
    p_rna_weights,
    file.path(
      mofa_dir,
      "MOFA2_RNA_Feature_Weights_Top20_Per_Factor_Heatmap.rds"
    )
  )
  
  write.csv(
    rna_heatmap_data,
    file.path(
      mofa_dir,
      "MOFA2_RNA_Feature_Weights_Top20_Per_Factor_Heatmap_Data.csv"
    ),
    row.names = FALSE
  )
  
  cat(
    "\nRNA feature-weight heatmap saved successfully.\n"
    
    
    library(ggplot2)
    library(dplyr)
    
    mofa_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
    
    protein_top_annotated <- readRDS(
      file.path(
        mofa_dir,
        "MOFA2_Top20_Protein_Features_Per_Factor_Annotated.rds"
      )
    )
    
    protein_heatmap_data <- protein_top_annotated %>%
      dplyr::mutate(
        feature_label = ifelse(
          is.na(Assay) | Assay == "",
          feature,
          Assay
        ),
        factor = factor(
          factor,
          levels = paste0("Factor", 1:15)
        )
      ) %>%
      dplyr::arrange(factor, rank)
    
    cat("\nProtein heatmap data:\n")
    cat("Rows:", nrow(protein_heatmap_data), "\n")
    cat("Factors:", dplyr::n_distinct(protein_heatmap_data$factor), "\n")
    cat("Features per factor:\n")
    print(table(protein_heatmap_data$factor))
    
    p_protein_weights <- ggplot(
      protein_heatmap_data,
      aes(
        x = factor,
        y = reorder(feature_label, abs_weight),
        fill = value
      )
    ) +
      geom_tile(color = "white", linewidth = 0.15) +
      scale_fill_gradient2(
        midpoint = 0,
        name = "MOFA weight"
      ) +
      labs(
        title = "Top 20 Protein Features per MOFA2 Factor",
        subtitle = "Features ranked by absolute feature weight",
        x = "MOFA2 Factor",
        y = "Protein assay"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(
          size = 7
        ),
        panel.grid = element_blank()
      )
    
    print(p_protein_weights)
    
    ggsave(
      file.path(
        mofa_dir,
        "MOFA2_Protein_Feature_Weights_Top20_Per_Factor_Heatmap.png"
      ),
      p_protein_weights,
      width = 12,
      height = 14,
      dpi = 300
    )
    
    saveRDS(
      p_protein_weights,
      file.path(
        mofa_dir,
        "MOFA2_Protein_Feature_Weights_Top20_Per_Factor_Heatmap.rds"
      )
    )
    
    write.csv(
      protein_heatmap_data,
      file.path(
        mofa_dir,
        "MOFA2_Protein_Feature_Weights_Top20_Per_Factor_Heatmap_Data.csv"
      ),
      row.names = FALSE
    )
    
    cat(
      "\nProtein feature-weight heatmap saved successfully.\n"
    )    
  )  
  
  mofa_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  files_to_check <- c(
    "MOFA2_RNA_Feature_Weights_Top20_Per_Factor_Heatmap.png",
    "MOFA2_RNA_Feature_Weights_Top20_Per_Factor_Heatmap.rds",
    "MOFA2_RNA_Feature_Weights_Top20_Per_Factor_Heatmap_Data.csv",
    "MOFA2_Protein_Feature_Weights_Top20_Per_Factor_Heatmap.png",
    "MOFA2_Protein_Feature_Weights_Top20_Per_Factor_Heatmap.rds",
    "MOFA2_Protein_Feature_Weights_Top20_Per_Factor_Heatmap_Data.csv"
  )
  
  qc_files <- data.frame(
    file = files_to_check,
    exists = file.exists(file.path(mofa_dir, files_to_check)),
    size_bytes = ifelse(
      file.exists(file.path(mofa_dir, files_to_check)),
      file.info(file.path(mofa_dir, files_to_check))$size,
      NA
    )
  )
  
  print(qc_files)
  
  rna_heatmap_check <- readRDS(
    file.path(
      mofa_dir,
      "MOFA2_RNA_Feature_Weights_Top20_Per_Factor_Heatmap.rds"
    )
  )
  
  protein_heatmap_check <- readRDS(
    file.path(
      mofa_dir,
      "MOFA2_Protein_Feature_Weights_Top20_Per_Factor_Heatmap.rds"
    )
  )
  
  cat("\nRNA heatmap object class:", class(rna_heatmap_check), "\n")
  cat("Protein heatmap object class:", class(protein_heatmap_check), "\n")
  
  cat(
    "\nRNA data rows:",
    nrow(read.csv(
      file.path(
        mofa_dir,
        "MOFA2_RNA_Feature_Weights_Top20_Per_Factor_Heatmap_Data.csv"
      )
    )),
    "\n"
  )
  
  cat(
    "Protein data rows:",
    nrow(read.csv(
      file.path(
        mofa_dir,
        "MOFA2_Protein_Feature_Weights_Top20_Per_Factor_Heatmap_Data.csv"
      )
    )),
    "\n"
  )
  
  cat("\nFeature-weight heatmap QC completed.\n")  
  
  
  library(MOFA2)
  
  mofa_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  mofa_file <- file.path(
    mofa_dir,
    "MOFA2_RNA_Protein_631_training.hdf5"
  )
  
  cat("\nChecking trained MOFA2 model:\n")
  cat("Model exists:", file.exists(mofa_file), "\n")
  
  mofa_object <- MOFA2::load_model(mofa_file)
  
  cat("\nMOFA2 model loaded successfully.\n")
  cat("Number of factors:", ncol(MOFA2::get_factors(
    mofa_object,
    factors = "all",
    groups = "all",
    as.data.frame = TRUE
  )) %>% { length(unique(MOFA2::get_factors(
    mofa_object,
    factors = "all",
    groups = "all",
    as.data.frame = TRUE
  )$factor)) }, "\n")
  
  variance_data <- MOFA2::calculate_variance_explained(
    mofa_object
  )
  
  cat("\nVariance explained object:\n")
  print(names(variance_data))
  
  cat("\nTotal variance explained:\n")
  print(variance_data$r2_total)
  
  cat("\nPer-factor variance explained:\n")
  print(variance_data$r2_per_factor)
  
  saveRDS(
    variance_data,
    file.path(
      mofa_dir,
      "MOFA2_Variance_Explained.rds"
    )
  )
  
  if (!is.null(variance_data$r2_per_factor)) {
    r2_factor <- as.data.frame(
      variance_data$r2_per_factor
    )
    
    r2_factor$Factor <- rownames(r2_factor)
    
    r2_factor <- r2_factor %>%
      dplyr::relocate(Factor)
    
    r2_factor$Total_R2 <- rowSums(
      r2_factor[, setdiff(
        names(r2_factor),
        c("Factor", "Total_R2")
      ), drop = FALSE],
      na.rm = TRUE
    )
    
    write.csv(
      r2_factor,
      file.path(
        mofa_dir,
        "MOFA2_Variance_Explained_Per_Factor.csv"
      ),
      row.names = FALSE
    )
    
    cat(
      "\nPer-factor variance table saved:",
      nrow(r2_factor),
      "factors\n"
    )
  }
  
  if (!is.null(variance_data$r2_total)) {
    r2_total <- as.data.frame(
      variance_data$r2_total
    )
    
    r2_total$View <- rownames(r2_total)
    
    r2_total <- r2_total %>%
      dplyr::relocate(View)
    
    write.csv(
      r2_total,
      file.path(
        mofa_dir,
        "MOFA2_Variance_Explained_Total.csv"
      ),
      row.names = FALSE
    )
    
    cat(
      "Total variance table saved:",
      nrow(r2_total),
      "views\n"
    )
  }
  
  cat("\nVariance explained files recreated successfully.\n") 
  
  
  library(MOFA2)
  
  model_metrics_table <- data.frame(
    metric = c(
      "Number_of_views",
      "Number_of_samples",
      "Number_of_factors",
      "Total_factor_score_rows"
    ),
    value = c(
      length(mofa_object@data),
      ncol(mofa_object@data[[1]]),
      length(unique(model_metrics$factors$factor)),
      nrow(model_metrics$factors)
    )
  )
  
  print(model_metrics_table)
  
  write.csv(
    model_metrics_table,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Model_Metrics.csv",
    row.names = FALSE
  )
  
  saveRDS(
    model_metrics_table,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Model_Metrics.rds"
  )
  
  cat("\nFiles saved successfully:\n")
  cat("MOFA2_Model_Metrics.csv\n")
  cat("MOFA2_Model_Metrics.rds\n")  
  
  library(MOFA2)
  
  MOFAobject_tutorial <- create_mofa(mofa_views)
  
  print(MOFAobject_tutorial)
  
  saveRDS(
    MOFAobject_tutorial,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_DataOverview_Object.rds"
  )  
  plot_data_overview(MOFAobject_tutorial) 
  library(MOFA2)
  
  png(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Data_Overview.png",
    width = 1800,
    height = 1200,
    res = 180
  )
  
  plot_data_overview(MOFAobject_tutorial)
  
  dev.off() 
  
  library(MOFA2)
  
  mofa_metadata <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Metadata_Aligned_631.rds"
  )
  
  mofa_metadata <- as.data.frame(mofa_metadata)
  
  # Check current columns
  cat("Current metadata columns:\n")
  print(colnames(mofa_metadata))
  
  # Rename MOFA_Key to sample
  colnames(mofa_metadata)[colnames(mofa_metadata) == "MOFA_Key"] <- "sample"
  
  cat("\nUpdated metadata columns:\n")
  print(colnames(mofa_metadata))
  
  # Validate sample IDs
  cat("\nNumber of samples:", nrow(mofa_metadata), "\n")
  cat("Unique sample IDs:", length(unique(mofa_metadata$sample)), "\n")
  cat("Missing sample IDs:", sum(is.na(mofa_metadata$sample)), "\n")
  
  # Add metadata to MOFA model
  samples_metadata(MOFAobject_tutorial) <- mofa_metadata
  
  cat("\nMetadata successfully added to MOFA object.\n")
  
  # Display first rows
  print(head(samples_metadata(MOFAobject_tutorial), 3))
  
  # Save updated object and metadata
  saveRDS(
    MOFAobject_tutorial,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Object_With_Metadata.rds"
  )
  
  saveRDS(
    mofa_metadata,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Samples_Metadata.rds"
  )
  
  cat("\nFiles saved successfully.\n")  
  
  
  library(MOFA2)
  
  # Restore the MOFA object with metadata
  MOFAobject_tutorial <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Object_With_Metadata.rds"
  )
  
  # Check MOFA dimensions
  cat("MOFA dimensions:\n")
  print(get_dimensions(MOFAobject_tutorial))
  
  # Overview of the data
  p_overview <- plot_data_overview(MOFAobject_tutorial)
  
  print(p_overview)
  
  # Save overview plot
  ggsave(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Data_Overview.png",
    p_overview,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  saveRDS(
    p_overview,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Data_Overview.rds"
  )
  
  cat("\nData overview plot saved successfully.\n")  
  
  library(MOFA2)
  
  # Load the trained MOFA model
  MOFAobject_tutorial <- MOFA2::load_model(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_RNA_Protein_631_training.hdf5"
  )
  
  cat("Trained MOFA model dimensions:\n")
  print(get_dimensions(MOFAobject_tutorial))
  
  # Load the validated metadata
  mofa_metadata <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Samples_Metadata.rds"
  )
  
  mofa_metadata <- as.data.frame(mofa_metadata)
  
  # Add metadata to the trained model
  samples_metadata(MOFAobject_tutorial) <- mofa_metadata
  
  cat("\nMetadata added to trained model.\n")
  
  # Check factors
  cat("\nNumber of factors (K):\n")
  print(get_dimensions(MOFAobject_tutorial)$K)
  
  # Check model class
  cat("\nMOFA object class:\n")
  print(class(MOFAobject_tutorial))
  
  # Save the correct downstream-analysis object
  saveRDS(
    MOFAobject_tutorial,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Trained_Model_With_Metadata.rds"
  )
  
  cat("\nTrained MOFA model with metadata saved successfully.\n") 
  
  
  library(MOFA2)
  
  MOFAobject_tutorial <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Trained_Model_With_Metadata.rds"
  )
  
  # Calculate variance explained
  variance_data <- get_variance_explained(MOFAobject_tutorial)
  
  # Total variance explained by each view
  r2_total <- variance_data$r2_total[[1]]
  
  cat("Total variance explained:\n")
  print(r2_total)
  
  # Variance explained by each factor and view
  r2_per_factor <- variance_data$r2_per_factor[[1]]
  
  cat("\nVariance explained per factor:\n")
  print(r2_per_factor)
  
  # Save the complete variance decomposition object
  saveRDS(
    variance_data,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Variance_Decomposition.rds"
  )
  
  # Save total variance explained
  write.csv(
    as.data.frame(r2_total),
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Variance_Explained_Total.csv"
  )
  
  # Save per-factor variance explained
  write.csv(
    as.data.frame(r2_per_factor),
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Variance_Explained_Per_Factor.csv"
  )
  
  cat("\nVariance decomposition calculated and saved successfully.\n") 
  
  
  library(MOFA2)
  library(ggplot2)
  
  MOFAobject_tutorial <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Trained_Model_With_Metadata.rds"
  )
  
  # Official MOFA2 variance-explained plot
  p_variance_factors <- plot_variance_explained(
    MOFAobject_tutorial,
    x = "view",
    y = "factor"
  )
  
  print(p_variance_factors)
  
  # Save official MOFA2 plot
  ggsave(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Variance_Explained_Per_Factor_Official.png",
    p_variance_factors,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  saveRDS(
    p_variance_factors,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Variance_Explained_Per_Factor_Official.rds"
  )
  
  cat("\nOfficial MOFA2 per-factor variance plot saved successfully.\n")  
  
  
  library(MOFA2)
  library(ggplot2)
  
  MOFAobject_tutorial <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Trained_Model_With_Metadata.rds"
  )
  
  # Official MOFA2 factor visualization
  p_factor_1_3 <- plot_factor(
    MOFAobject_tutorial,
    factors = c(1, 2, 3),
    color_by = "Timepoint",
    shape_by = "Acuity_max",
    dot_size = 3,
    dodge = TRUE,
    legend = TRUE,
    add_violin = TRUE,
    violin_alpha = 0.25
  )
  
  print(p_factor_1_3)
  
  # Save plot
  ggsave(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Factor_1_3_Timepoint_Acuity.png",
    p_factor_1_3,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  saveRDS(
    p_factor_1_3,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Factor_1_3_Timepoint_Acuity.rds"
  )
  
  cat("\nFactor 1-3 visualization saved successfully.\n") 
  
  library(MOFA2)
  library(ggplot2)
  
  MOFAobject_tutorial <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Trained_Model_With_Metadata.rds"
  )
  
  # Official MOFA2 factor combination visualization
  p_factor_combination <- plot_factors(
    MOFAobject_tutorial,
    factors = c(1, 2, 3),
    color_by = "Timepoint"
  )
  
  print(p_factor_combination)
  
  # Save plot
  ggsave(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Factor_Combination_F1_F2_F3.png",
    p_factor_combination,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  saveRDS(
    p_factor_combination,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Factor_Combination_F1_F2_F3.rds"
  )
  
  cat("\nFactor combination visualization saved successfully.\n")  
  
  
  library(MOFA2)
  
  MOFAobject_tutorial <- readRDS(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Trained_Model_With_Metadata.rds"
  )
  
  p_data_overview <- plot_data_overview(MOFAobject_tutorial)
  
  print(p_data_overview)
  
  ggsave(
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Data_Overview.png",
    plot = p_data_overview,
    width = 10,
    height = 7,
    dpi = 300
  )
  
  saveRDS(
    p_data_overview,
    "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2/MOFA2_Tutorial_Data_Overview.rds"
  )  
  
  library(MOFA2)
  library(ggplot2)
  
  mofa_factor_batches <- list(
    Factor_1_5 = 1:5,
    Factor_6_10 = 6:10,
    Factor_11_15 = 11:15
  )
  
  for (batch_name in names(mofa_factor_batches)) {
    
    factor_set <- mofa_factor_batches[[batch_name]]
    
    p <- plot_factor(
      MOFAobject_tutorial,
      factors = factor_set,
      color_by = "Timepoint",
      shape_by = "Acuity_max",
      dot_size = 3,
      dodge = TRUE,
      legend = TRUE,
      add_violin = TRUE,
      violin_alpha = 0.25
    )
    
    print(p)
    
    ggsave(
      file.path(
        "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2",
        paste0("MOFA2_Tutorial_Factor_Visualization_", batch_name, ".png")
      ),
      plot = p,
      width = 14,
      height = 10,
      dpi = 300
    )
    
    saveRDS(
      p,
      file.path(
        "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2",
        paste0("MOFA2_Tutorial_Factor_Visualization_", batch_name, ".rds")
      )
    )
  }
  
  cat("All 15 factor visualization plots were generated and saved successfully.\n")  
  
  
  library(MOFA2)
  library(ggplot2)
  
  factor_combination_batches <- list(
    Factor_1_5 = 1:5,
    Factor_6_10 = 6:10,
    Factor_11_15 = 11:15
  )
  
  for (batch_name in names(factor_combination_batches)) {
    
    factor_set <- factor_combination_batches[[batch_name]]
    
    p <- plot_factors(
      MOFAobject_tutorial,
      factors = factor_set,
      color_by = "Timepoint"
    )
    
    print(p)
    
    ggsave(
      file.path(
        "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2",
        paste0("MOFA2_Tutorial_Factor_Combinations_", batch_name, ".png")
      ),
      plot = p,
      width = 14,
      height = 12,
      dpi = 300
    )
    
    saveRDS(
      p,
      file.path(
        "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2",
        paste0("MOFA2_Tutorial_Factor_Combinations_", batch_name, ".rds")
      )
    )
  }
  
  cat("All 15 factors were included in the factor-combination analysis.\n") 
  
  library(MOFA2)
  library(ggplot2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  for (f in 1:15) {
    
    p <- plot_factor(
      MOFAobject_tutorial,
      factor = f,
      color_by = "Timepoint",
      shape_by = "Acuity_max",
      dot_size = 3,
      dodge = TRUE,
      legend = TRUE,
      add_violin = TRUE,
      violin_alpha = 0.25
    )
    
    print(p)
    
    ggsave(
      file.path(
        output_dir,
        paste0("MOFA2_Tutorial_Factor_", f, "_One_at_a_Time.png")
      ),
      plot = p,
      width = 8,
      height = 7,
      dpi = 300
    )
    
    saveRDS(
      p,
      file.path(
        output_dir,
        paste0("MOFA2_Tutorial_Factor_", f, "_One_at_a_Time.rds")
      )
    )
  }
  
  cat("All 15 factors were visualized individually and saved successfully.\n")  
  
  library(MOFA2)
  library(ggplot2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  views <- c("RNA", "Protein")
  
  for (v in views) {
    
    for (f in 1:15) {
      
      p_weights <- plot_weights(
        MOFAobject_tutorial,
        view = v,
        factor = f,
        nfeatures = 10,
        scale = TRUE,
        abs = FALSE
      )
      
      print(p_weights)
      
      ggsave(
        file.path(
          output_dir,
          paste0(
            "MOFA2_Tutorial_Plot_Weights_",
            v,
            "_Factor_",
            f,
            ".png"
          )
        ),
        plot = p_weights,
        width = 9,
        height = 7,
        dpi = 300
      )
      
      saveRDS(
        p_weights,
        file.path(
          output_dir,
          paste0(
            "MOFA2_Tutorial_Plot_Weights_",
            v,
            "_Factor_",
            f,
            ".rds"
          )
        )
      )
      
      p_top <- plot_top_weights(
        MOFAobject_tutorial,
        view = v,
        factor = f,
        nfeatures = 10
      )
      
      print(p_top)
      
      ggsave(
        file.path(
          output_dir,
          paste0(
            "MOFA2_Tutorial_Plot_Top_Weights_",
            v,
            "_Factor_",
            f,
            ".png"
          )
        ),
        plot = p_top,
        width = 9,
        height = 7,
        dpi = 300
      )
      
      saveRDS(
        p_top,
        file.path(
          output_dir,
          paste0(
            "MOFA2_Tutorial_Plot_Top_Weights_",
            v,
            "_Factor_",
            f,
            ".rds"
          )
        )
      )
    }
  }
  
  cat("Feature-weight analysis completed for all 15 factors across RNA and Protein views.\n")
  cat("Both plot_weights() and plot_top_weights() were generated and saved.\n") 
  
  
  library(MOFA2)
  library(ggplot2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  for (v in c("RNA", "Protein")) {
    
    for (f in 1:15) {
      
      p_heatmap <- plot_data_heatmap(
        MOFAobject_tutorial,
        view = v,
        factor = f,
        features = 20,
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_rownames = TRUE,
        show_colnames = FALSE
      )
      
      print(p_heatmap)
      
      ggsave(
        file.path(
          output_dir,
          paste0(
            "MOFA2_Tutorial_Z_Input_Heatmap_",
            v,
            "_Factor_",
            f,
            ".png"
          )
        ),
        plot = p_heatmap,
        width = 10,
        height = 8,
        dpi = 300
      )
      
      saveRDS(
        p_heatmap,
        file.path(
          output_dir,
          paste0(
            "MOFA2_Tutorial_Z_Input_Heatmap_",
            v,
            "_Factor_",
            f,
            ".rds"
          )
        )
      )
    }
  }
  
  cat("Z-input data heatmaps completed for all 15 factors across RNA and Protein views.\n")  
  library(MOFA2)
  library(ggplot2)
  library(ggpubr)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  for (v in c("RNA", "Protein")) {
    
    for (f in 1:15) {
      
      p_scatter <- plot_data_scatter(
        MOFAobject_tutorial,
        view = v,
        factor = f,
        features = 5,
        add_lm = TRUE,
        color_by = "Timepoint"
      )
      
      print(p_scatter)
      
      ggsave(
        file.path(
          output_dir,
          paste0(
            "MOFA2_Tutorial_Scatter_",
            v,
            "_Factor_",
            f,
            ".png"
          )
        ),
        plot = p_scatter,
        width = 12,
        height = 9,
        dpi = 300
      )
      
      saveRDS(
        p_scatter,
        file.path(
          output_dir,
          paste0(
            "MOFA2_Tutorial_Scatter_",
            v,
            "_Factor_",
            f,
            ".rds"
          )
        )
      )
    }
  }
  
  cat("Scatter plots completed for all 15 factors across RNA and Protein views.\n") 
  
  
  
  library(MOFA2)
  library(ggplot2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  set.seed(42)
  
  MOFAobject_tutorial <- run_umap(MOFAobject_tutorial)
  MOFAobject_tutorial <- run_tsne(MOFAobject_tutorial)
  
  p_umap <- plot_dimred(
    MOFAobject_tutorial,
    method = "UMAP",
    color_by = "Timepoint",
    dot_size = 5
  )
  
  print(p_umap)
  
  ggsave(
    file.path(output_dir, "MOFA2_Tutorial_UMAP_Timepoint.png"),
    plot = p_umap,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  saveRDS(
    p_umap,
    file.path(output_dir, "MOFA2_Tutorial_UMAP_Timepoint.rds")
  )
  
  p_tsne <- plot_dimred(
    MOFAobject_tutorial,
    method = "TSNE",
    color_by = "Timepoint",
    dot_size = 5
  )
  
  print(p_tsne)
  
  ggsave(
    file.path(output_dir, "MOFA2_Tutorial_TSNE_Timepoint.png"),
    plot = p_tsne,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  saveRDS(
    p_tsne,
    file.path(output_dir, "MOFA2_Tutorial_TSNE_Timepoint.rds")
  )
  
  saveRDS(
    MOFAobject_tutorial,
    file.path(output_dir, "MOFA2_Tutorial_Model_With_UMAP_TSNE.rds")
  )
  
  cat("UMAP and t-SNE dimensionality-reduction analyses completed successfully.\n")  
  
  library(MOFA2)
  library(MOFAdata)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  data("reactomeGS", package = "MOFAdata")
  
  cat("Reactome gene-set matrix dimensions:\n")
  print(dim(reactomeGS))
  
  cat("\nFirst feature sets:\n")
  print(head(rownames(reactomeGS)))
  
  cat("\nFirst features:\n")
  print(head(colnames(reactomeGS)))
  
  rna_features <- features_names(MOFAobject_tutorial)[["RNA"]]
  
  common_rna_features <- intersect(rna_features, colnames(reactomeGS))
  
  cat("\nRNA features in MOFA model:", length(rna_features), "\n")
  cat("Reactome features:", ncol(reactomeGS), "\n")
  cat("Common RNA features:", length(common_rna_features), "\n")
  
  reactomeGS_mofa <- reactomeGS[, common_rna_features, drop = FALSE]
  
  cat("\nFinal Reactome matrix for MOFA:\n")
  print(dim(reactomeGS_mofa))
  
  saveRDS(
    reactomeGS_mofa,
    file.path(output_dir, "MOFA2_Reactome_GeneSets_Matched_to_RNA.rds")
  )
  
  write.csv(
    reactomeGS_mofa,
    file.path(output_dir, "MOFA2_Reactome_GeneSets_Matched_to_RNA.csv"),
    quote = FALSE
  )
  
  cat("\nReactome gene-set preparation completed successfully.\n")
  
  library(MOFA2)
  
  cat("run_enrichment() arguments:\n")
  print(args(run_enrichment))
  
  cat("\nplot_enrichment() arguments:\n")
  print(args(plot_enrichment))
  
  cat("\nplot_enrichment_heatmap() arguments:\n")
  print(args(plot_enrichment_heatmap))


  library(MOFA2)
  
  cat("RNA feature IDs:\n")
  rna_features <- features_names(MOFAobject_tutorial)[["RNA"]]
  print(length(rna_features))
  print(head(rna_features))
  
  cat("\nProtein feature IDs:\n")
  protein_features <- features_names(MOFAobject_tutorial)[["Protein"]]
  print(length(protein_features))
  print(head(protein_features))
  
  cat("\nCurrent MOFA views:\n")
  print(names(MOFAobject_tutorial@data))
  
  library(MOFA2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  rna_features <- features_names(MOFAobject_tutorial)[["RNA"]]
  
  gene_symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = rna_features,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  gene_symbols <- unique(gene_symbols)
  
  cat("MOFA RNA features:", length(rna_features), "\n")
  cat("Mapped unique gene symbols:", length(gene_symbols), "\n")
  
  go_bp <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_symbols,
    keytype = "SYMBOL",
    columns = c("GO", "ONTOLOGY")
  )
  
  go_bp <- go_bp[
    !is.na(go_bp$GO) &
      go_bp$ONTOLOGY == "BP",
    c("GO", "SYMBOL")
  ]
  
  go_bp <- unique(go_bp)
  
  feature_sets_go_bp <- split(go_bp$SYMBOL, go_bp$GO)
  
  feature_sets_go_bp <- feature_sets_go_bp[
    lengths(feature_sets_go_bp) >= 10
  ]
  
  cat("GO-BP gene sets with >=10 genes:", length(feature_sets_go_bp), "\n")
  cat("Genes represented in GO-BP sets:", length(unique(unlist(feature_sets_go_bp))), "\n")
  
  saveRDS(
    feature_sets_go_bp,
    file.path(output_dir, "MOFA2_GO_BP_Feature_Sets.rds")
  )
  
  cat("\nGO-BP feature-set preparation completed successfully.\n")  

  library(MOFA2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  feature_sets_go_bp <- readRDS(
    file.path(output_dir, "MOFA2_GO_BP_Feature_Sets.rds")
  )
  
  cat("Running MOFA2 GO-BP enrichment for all 15 factors...\n")
  
  mofa_go_bp_enrichment <- run_enrichment(
    object = MOFAobject_tutorial,
    view = "RNA",
    feature.sets = feature_sets_go_bp,
    factors = 1:15,
    set.statistic = "mean.diff",
    statistical.test = "parametric",
    sign = "all",
    min.size = 10,
    p.adj.method = "BH",
    alpha = 0.1,
    verbose = TRUE
  )
  
  cat("\nEnrichment completed successfully.\n")
  cat("Result dimensions:\n")
  print(dim(mofa_go_bp_enrichment))
  
  cat("\nResult columns:\n")
  print(colnames(mofa_go_bp_enrichment))
  
  saveRDS(
    mofa_go_bp_enrichment,
    file.path(output_dir, "MOFA2_GO_BP_Enrichment_All_15_Factors.rds")
  )
  
  write.csv(
    mofa_go_bp_enrichment,
    file.path(output_dir, "MOFA2_GO_BP_Enrichment_All_15_Factors.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("\nMOFA2 GO-BP enrichment results saved successfully.\n")  
  
  library(MOFA2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  feature_sets_go_bp <- readRDS(
    file.path(output_dir, "MOFA2_GO_BP_Feature_Sets.rds")
  )
  
  rna_features <- features_names(MOFAobject_tutorial)[["RNA"]]
  
  rna_symbol_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = rna_features,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  rna_symbol_map <- rna_symbol_map[!is.na(rna_symbol_map)]
  
  go_bp_binary <- matrix(
    0L,
    nrow = length(feature_sets_go_bp),
    ncol = length(rna_features),
    dimnames = list(
      names(feature_sets_go_bp),
      rna_features
    )
  )
  
  for (i in seq_along(feature_sets_go_bp)) {
    genes_in_set <- feature_sets_go_bp[[i]]
    
    matching_features <- names(rna_symbol_map)[
      rna_symbol_map %in% genes_in_set
    ]
    
    if (length(matching_features) > 0) {
      go_bp_binary[i, matching_features] <- 1L
    }
  }
  
  cat("GO-BP binary matrix dimensions:\n")
  print(dim(go_bp_binary))
  
  cat("\nNumber of gene sets:\n")
  print(nrow(go_bp_binary))
  
  cat("\nNumber of MOFA RNA features:\n")
  print(ncol(go_bp_binary))
  
  cat("\nGene-set sizes in MOFA features:\n")
  print(summary(rowSums(go_bp_binary)))
  
  cat("\nFeatures represented in at least one gene set:\n")
  print(sum(colSums(go_bp_binary) > 0))
  
  saveRDS(
    go_bp_binary,
    file.path(output_dir, "MOFA2_GO_BP_Feature_Sets_Binary_Matrix.rds")
  )
  
  write.csv(
    go_bp_binary,
    file.path(output_dir, "MOFA2_GO_BP_Feature_Sets_Binary_Matrix.csv"),
    quote = FALSE
  )
  
  cat("\nGO-BP binary feature-set matrix created and saved successfully.\n")  
  
  library(MOFA2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  go_bp_binary <- readRDS(
    file.path(output_dir, "MOFA2_GO_BP_Feature_Sets_Binary_Matrix.rds")
  )
  
  cat("Running MOFA2 GO-BP enrichment for all 15 factors...\n")
  
  mofa_go_bp_enrichment <- run_enrichment(
    object = MOFAobject_tutorial,
    view = "RNA",
    feature.sets = go_bp_binary,
    factors = 1:15,
    set.statistic = "mean.diff",
    statistical.test = "parametric",
    sign = "all",
    min.size = 10,
    p.adj.method = "BH",
    alpha = 0.1,
    verbose = TRUE
  )
  
  cat("\nEnrichment completed successfully.\n")
  
  cat("\nResult dimensions:\n")
  print(dim(mofa_go_bp_enrichment))
  
  cat("\nResult columns:\n")
  print(colnames(mofa_go_bp_enrichment))
  
  cat("\nSignificant results (FDR < 0.1):\n")
  print(sum(mofa_go_bp_enrichment$pval.adj < 0.1, na.rm = TRUE))
  
  saveRDS(
    mofa_go_bp_enrichment,
    file.path(output_dir, "MOFA2_GO_BP_Enrichment_All_15_Factors.rds")
  )
  
  write.csv(
    mofa_go_bp_enrichment,
    file.path(output_dir, "MOFA2_GO_BP_Enrichment_All_15_Factors.csv"),
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("\nMOFA2 GO-BP enrichment results saved successfully.\n")
  
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  pval_matrix <- mofa_go_bp_enrichment$pval
  padj_matrix <- mofa_go_bp_enrichment$pval.adj
  
  cat("P-value matrix dimensions:\n")
  print(dim(pval_matrix))
  
  cat("\nAdjusted p-value matrix dimensions:\n")
  print(dim(padj_matrix))
  
  cat("\nNumber of significant factor-pathway pairs (FDR < 0.1):\n")
  print(sum(padj_matrix < 0.1, na.rm = TRUE))
  
  cat("\nSignificant pathways per factor:\n")
  print(colSums(padj_matrix < 0.1, na.rm = TRUE))
  
  pathway_ids <- rownames(padj_matrix)
  factor_names <- colnames(padj_matrix)
  
  enrichment_table <- do.call(
    rbind,
    lapply(seq_along(factor_names), function(i) {
      data.frame(
        Factor = factor_names[i],
        GO_ID = pathway_ids,
        P_value = pval_matrix[, i],
        FDR = padj_matrix[, i],
        stringsAsFactors = FALSE
      )
    })
  )
  
  enrichment_table <- enrichment_table[
    order(enrichment_table$Factor, enrichment_table$FDR),
  ]
  
  cat("\nFinal enrichment table dimensions:\n")
  print(dim(enrichment_table))
  
  cat("\nTop significant results:\n")
  print(
    head(
      enrichment_table[
        enrichment_table$FDR < 0.1,
      ],
      20
    )
  )
  
  saveRDS(
    enrichment_table,
    file.path(
      output_dir,
      "MOFA2_GO_BP_Enrichment_All_15_Factors_Tidy.rds"
    )
  )
  
  write.csv(
    enrichment_table,
    file.path(
      output_dir,
      "MOFA2_GO_BP_Enrichment_All_15_Factors_Tidy.csv"
    ),
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("\nTidy MOFA2 enrichment results saved successfully.\n")  
  
  library(MOFA2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  p_enrichment_F1 <- plot_enrichment(
    enrichment.results = mofa_go_bp_enrichment,
    factor = 1,
    alpha = 0.1,
    max.pathways = 25,
    text_size = 1,
    dot_size = 5
  )
  
  print(p_enrichment_F1)
  
  ggsave(
    filename = file.path(
      output_dir,
      "MOFA2_GO_BP_Enrichment_Factor1_Official.png"
    ),
    plot = p_enrichment_F1,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  saveRDS(
    p_enrichment_F1,
    file.path(
      output_dir,
      "MOFA2_GO_BP_Enrichment_Factor1_Official.rds"
    )
  )
  
  cat("\nFactor 1 enrichment plot saved successfully.\n")  
  
  library(MOFA2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  for (f in 1:15) {
    
    p <- plot_enrichment(
      enrichment.results = mofa_go_bp_enrichment,
      factor = f,
      alpha = 0.1,
      max.pathways = 25,
      text_size = 1,
      dot_size = 5
    )
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0("MOFA2_GO_BP_Enrichment_Factor", f, "_Official.png")
      ),
      plot = p,
      width = 10,
      height = 8,
      dpi = 300
    )
    
    saveRDS(
      p,
      file.path(
        output_dir,
        paste0("MOFA2_GO_BP_Enrichment_Factor", f, "_Official.rds")
      )
    )
  }
  
  cat("Official MOFA2 enrichment plots completed for all 15 factors.\n")
  cat("15 PNG files and 15 RDS files were saved successfully.\n") 
  
  library(MOFA2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  p_enrichment_heatmap <- plot_enrichment_heatmap(
    enrichment.results = mofa_go_bp_enrichment,
    alpha = 0.1,
    cap = 1e-50,
    log_scale = TRUE
  )
  
  print(p_enrichment_heatmap)
  
  ggsave(
    filename = file.path(
      output_dir,
      "MOFA2_GO_BP_Enrichment_Heatmap_Official.png"
    ),
    plot = p_enrichment_heatmap,
    width = 12,
    height = 10,
    dpi = 300
  )
  
  saveRDS(
    p_enrichment_heatmap,
    file.path(
      output_dir,
      "MOFA2_GO_BP_Enrichment_Heatmap_Official.rds"
    )
  )
  
  cat("\nOfficial MOFA2 enrichment heatmap saved successfully.\n")  
  
  library(MOFA2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  p_factor_combination <- plot_factors(
    mofa_object,
    factors = 1:15
  )
  
  print(p_factor_combination)
  
  ggsave(
    filename = file.path(
      output_dir,
      "MOFA2_Factor_Combination_15_Factors_Official.png"
    ),
    plot = p_factor_combination,
    width = 14,
    height = 10,
    dpi = 300
  )
  
  saveRDS(
    p_factor_combination,
    file.path(
      output_dir,
      "MOFA2_Factor_Combination_15_Factors_Official.rds"
    )
  )
  
  cat("Official MOFA2 factor-combination plot saved successfully.\n") 
  
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  cat("Plot class:\n")
  print(class(p_z_heatmap_f1_rna))
  
  ggsave(
    filename = file.path(
      output_dir,
      "MOFA2_Z_Input_Data_Heatmap_Factor1_RNA_Official.png"
    ),
    plot = p_z_heatmap_f1_rna,
    width = 14,
    height = 12,
    dpi = 300
  )
  
  saveRDS(
    p_z_heatmap_f1_rna,
    file.path(
      output_dir,
      "MOFA2_Z_Input_Data_Heatmap_Factor1_RNA_Official.rds"
    )
  )
  
  cat("\nOfficial MOFA2 Z-input heatmap for Factor 1 (RNA) saved successfully.\n")  
  
  library(ggplot2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  factor5_data <- factor_scores_long[
    factor_scores_long$factor == "Factor5",
    c("sample", "value")
  ]
  
  metadata_for_merge <- mofa_metadata_final
  metadata_for_merge$sample <- rownames(metadata_for_merge)
  
  factor5_data <- merge(
    factor5_data,
    metadata_for_merge[, c("sample", "Acuity_max")],
    by = "sample"
  )
  
  factor5_data$Acuity_max <- as.numeric(factor5_data$Acuity_max)
  
  cat("Factor 5 observations:", nrow(factor5_data), "\n")
  cat("Missing Acuity:", sum(is.na(factor5_data$Acuity_max)), "\n")
  
  p_factor5_acuity <- ggplot(
    factor5_data,
    aes(x = Acuity_max, y = value)
  ) +
    geom_point(alpha = 0.5) +
    geom_smooth(
      method = "lm",
      se = TRUE
    ) +
    theme_classic() +
    labs(
      title = "MOFA2 Factor 5 vs Clinical Acuity",
      x = "Maximum acuity",
      y = "Factor 5 score"
    )
  
  print(p_factor5_acuity)
  
  ggsave(
    filename = file.path(
      output_dir,
      "MOFA2_Factor5_vs_Acuity_Scatter.png"
    ),
    plot = p_factor5_acuity,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  saveRDS(
    p_factor5_acuity,
    file.path(
      output_dir,
      "MOFA2_Factor5_vs_Acuity_Scatter.rds"
    )
  )
  
  write.csv(
    factor5_data,
    file.path(
      output_dir,
      "MOFA2_Factor5_vs_Acuity_Scatter_Data.csv"
    ),
    row.names = FALSE
  )
  
  cat("MOFA2 Factor 5 vs Acuity scatter plot saved successfully.\n") 
  
  library(MOFA2)
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  variance_data <- calculate_variance_explained(mofa_object)
  
  r2_factor <- variance_data$r2_per_factor$group1
  
  factor_summary <- data.frame(
    Factor = paste0("Factor", seq_len(nrow(r2_factor))),
    RNA_Variance = r2_factor[, "RNA"],
    Protein_Variance = r2_factor[, "Protein"]
  )
  
  gsea_pval <- mofa_go_bp_enrichment$pval
  gsea_fdr <- mofa_go_bp_enrichment$pval.adj
  
  gsea_sig_counts <- data.frame(
    Factor = colnames(gsea_fdr),
    Significant_GO_BP_Pathways = colSums(
      gsea_fdr < 0.1,
      na.rm = TRUE
    )
  )
  
  factor_summary <- merge(
    factor_summary,
    gsea_sig_counts,
    by = "Factor",
    all.x = TRUE
  )
  
  factor_summary$Significant_GO_BP_Pathways[
    is.na(factor_summary$Significant_GO_BP_Pathways)
  ] <- 0
  
  factor_summary$Total_Variance <-
    factor_summary$RNA_Variance +
    factor_summary$Protein_Variance
  
  factor_summary <- factor_summary[
    order(-factor_summary$Total_Variance),
  ]
  
  rownames(factor_summary) <- NULL
  
  print(factor_summary)
  
  write.csv(
    factor_summary,
    file.path(
      output_dir,
      "MOFA2_Final_Factor_Summary.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    factor_summary,
    file.path(
      output_dir,
      "MOFA2_Final_Factor_Summary.rds"
    )
  )
  
  cat("\nFinal MOFA2 factor summary saved successfully.\n")
  cat("Number of factors:", nrow(factor_summary), "\n")
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  factor_summary <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Final_Factor_Summary.rds"
    )
  )
  
  acuity_corr <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Factor_Acuity_Spearman_Correlations.rds"
    )
  )
  
  trajectory_stats <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Longitudinal_Factor_Statistics.rds"
    )
  )
  
  factor_summary$Acuity_Rho <- NA_real_
  factor_summary$Acuity_FDR <- NA_real_
  factor_summary$Acuity_Significant <- FALSE
  
  for (i in seq_len(nrow(factor_summary))) {
    
    current_factor <- factor_summary$Factor[i]
    
    hit <- acuity_corr[
      acuity_corr$factor == current_factor &
        acuity_corr$Timepoint == "Overall",
      ,
      drop = FALSE
    ]
    
    if (nrow(hit) == 1) {
      factor_summary$Acuity_Rho[i] <- hit$Spearman_rho
      factor_summary$Acuity_FDR[i] <- hit$FDR
      factor_summary$Acuity_Significant[i] <- hit$FDR < 0.05
    }
  }
  
  factor_summary$Longitudinal_Significant <- FALSE
  
  for (i in seq_len(nrow(factor_summary))) {
    
    current_factor <- factor_summary$Factor[i]
    
    hit <- trajectory_stats[
      trajectory_stats$factor == current_factor,
      ,
      drop = FALSE
    ]
    
    if (nrow(hit) > 0) {
      factor_summary$Longitudinal_Significant[i] <-
        any(hit$FDR < 0.05, na.rm = TRUE)
    }
  }
  
  factor_summary$Acuity_Abs_Rho <- abs(
    factor_summary$Acuity_Rho
  )
  
  factor_summary$Evidence_Count <-
    (factor_summary$Total_Variance >= 5) +
    factor_summary$Acuity_Significant +
    factor_summary$Longitudinal_Significant +
    (factor_summary$Significant_GO_BP_Pathways >= 20)
  
  factor_summary <- factor_summary[
    order(
      -factor_summary$Evidence_Count,
      -factor_summary$Acuity_Abs_Rho,
      -factor_summary$Total_Variance
    ),
  ]
  
  rownames(factor_summary) <- NULL
  
  cat("\nFINAL MOFA2 FACTOR EVIDENCE TABLE\n")
  print(factor_summary)
  
  write.csv(
    factor_summary,
    file.path(
      output_dir,
      "MOFA2_Final_Key_Factor_Evidence_Table.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    factor_summary,
    file.path(
      output_dir,
      "MOFA2_Final_Key_Factor_Evidence_Table.rds"
    )
  )
  
  key_factors <- subset(
    factor_summary,
    Evidence_Count >= 3
  )
  
  cat("\nKEY FACTORS IDENTIFIED\n")
  print(
    key_factors[, c(
      "Factor",
      "RNA_Variance",
      "Protein_Variance",
      "Total_Variance",
      "Significant_GO_BP_Pathways",
      "Acuity_Rho",
      "Acuity_FDR",
      "Longitudinal_Significant",
      "Evidence_Count"
    )]
  )
  
  write.csv(
    key_factors,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Final.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    key_factors,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Final.rds"
    )
  )
  
  cat("\nFinal key-factor selection saved successfully.\n")
  cat("Number of key factors:", nrow(key_factors), "\n") 
  
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  factor_summary <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Final_Factor_Summary.rds"
    )
  )
  
  acuity_corr <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Factor_Acuity_Spearman_Correlations.rds"
    )
  )
  
  trajectory_stats <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Longitudinal_Factor_Statistics.rds"
    )
  )
  
  factor_summary$Acuity_Rho <- NA_real_
  factor_summary$Acuity_FDR <- NA_real_
  factor_summary$Acuity_Significant <- FALSE
  
  for (i in seq_len(nrow(factor_summary))) {
    
    current_factor <- factor_summary$Factor[i]
    
    hit <- acuity_corr[
      acuity_corr$factor == current_factor &
        acuity_corr$Timepoint == "Overall",
      ,
      drop = FALSE
    ]
    
    if (nrow(hit) == 1) {
      factor_summary$Acuity_Rho[i] <- hit$Spearman_rho
      factor_summary$Acuity_FDR[i] <- hit$FDR
      factor_summary$Acuity_Significant[i] <- hit$FDR < 0.05
    }
  }
  
  factor_summary$Longitudinal_Significant <- FALSE
  
  for (i in seq_len(nrow(factor_summary))) {
    
    current_factor <- factor_summary$Factor[i]
    
    hit <- trajectory_stats[
      trajectory_stats$factor == current_factor,
      ,
      drop = FALSE
    ]
    
    if (nrow(hit) > 0) {
      factor_summary$Longitudinal_Significant[i] <-
        any(hit$FDR < 0.05, na.rm = TRUE)
    }
  }
  
  factor_summary$Acuity_Abs_Rho <- abs(
    factor_summary$Acuity_Rho
  )
  
  factor_summary$Evidence_Count <-
    (factor_summary$Total_Variance >= 5) +
    factor_summary$Acuity_Significant +
    factor_summary$Longitudinal_Significant +
    (factor_summary$Significant_GO_BP_Pathways >= 20)
  
  factor_summary <- factor_summary[
    order(
      -factor_summary$Evidence_Count,
      -factor_summary$Acuity_Abs_Rho,
      -factor_summary$Total_Variance
    ),
  ]
  
  rownames(factor_summary) <- NULL
  
  cat("\nFINAL MOFA2 FACTOR EVIDENCE TABLE\n")
  print(factor_summary)
  
  write.csv(
    factor_summary,
    file.path(
      output_dir,
      "MOFA2_Final_Key_Factor_Evidence_Table.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    factor_summary,
    file.path(
      output_dir,
      "MOFA2_Final_Key_Factor_Evidence_Table.rds"
    )
  )
  
  key_factors <- subset(
    factor_summary,
    Evidence_Count >= 3
  )
  
  cat("\nKEY FACTORS IDENTIFIED\n")
  print(
    key_factors[, c(
      "Factor",
      "RNA_Variance",
      "Protein_Variance",
      "Total_Variance",
      "Significant_GO_BP_Pathways",
      "Acuity_Rho",
      "Acuity_FDR",
      "Longitudinal_Significant",
      "Evidence_Count"
    )]
  )
  
  write.csv(
    key_factors,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Final.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    key_factors,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Final.rds"
    )
  )
  
  cat("\nFinal key-factor selection saved successfully.\n")
  cat("Number of key factors:", nrow(key_factors), "\n")  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  key_factors <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Final.rds"
    )
  )
  
  rna_top <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Top20_RNA_Features_Per_Factor_Annotated.rds"
    )
  )
  
  protein_top <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Top20_Protein_Features_Per_Factor_Annotated.rds"
    )
  )
  
  key_factor_names <- key_factors$Factor
  
  rna_key <- rna_top[
    rna_top$factor %in% key_factor_names,
    ,
    drop = FALSE
  ]
  
  protein_key <- protein_top[
    protein_top$factor %in% key_factor_names,
    ,
    drop = FALSE
  ]
  
  rna_key <- rna_key[
    order(rna_key$factor, -rna_key$abs_weight),
  ]
  
  protein_key <- protein_key[
    order(protein_key$factor, -protein_key$abs_weight),
  ]
  
  cat("\nRNA features for key factors:\n")
  print(dim(rna_key))
  print(table(rna_key$factor))
  
  cat("\nTop RNA features:\n")
  print(
    rna_key[
      ,
      intersect(
        c("factor", "rank", "feature", "symbol", "value", "abs_weight"),
        colnames(rna_key)
      )
    ]
  )
  
  cat("\nProtein features for key factors:\n")
  print(dim(protein_key))
  print(table(protein_key$factor))
  
  cat("\nTop protein features:\n")
  print(
    protein_key[
      ,
      intersect(
        c(
          "factor",
          "rank",
          "feature",
          "Assay",
          "UniProt",
          "Panel",
          "value",
          "abs_weight"
        ),
        colnames(protein_key)
      )
    ]
  )
  
  write.csv(
    rna_key,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Top20_RNA_Features.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    protein_key,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Top20_Protein_Features.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    rna_key,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Top20_RNA_Features.rds"
    )
  )
  
  saveRDS(
    protein_key,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Top20_Protein_Features.rds"
    )
  )
  
  cat("\nKey-factor feature extraction completed successfully.\n")
  cat("RNA rows:", nrow(rna_key), "\n")
  cat("Protein rows:", nrow(protein_key), "\n")  
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  key_factors <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Final.rds"
    )
  )
  
  gsea_tidy <- readRDS(
    file.path(
      output_dir,
      "MOFA2_GO_BP_Enrichment_All_15_Factors_Tidy.rds"
    )
  )
  
  key_factor_names <- as.character(key_factors$Factor)
  
  key_gsea <- gsea_tidy[
    as.character(gsea_tidy$Factor) %in% key_factor_names,
    ,
    drop = FALSE
  ]
  
  key_gsea <- key_gsea[
    order(
      as.character(key_gsea$Factor),
      key_gsea$FDR,
      key_gsea$P_value
    ),
  ]
  
  significant_key_gsea <- key_gsea[
    key_gsea$FDR < 0.10,
    ,
    drop = FALSE
  ]
  
  cat("\nKey-factor GO-BP enrichment:\n")
  print(dim(significant_key_gsea))
  
  cat("\nSignificant pathways per key factor:\n")
  print(table(significant_key_gsea$Factor))
  
  cat("\nTop 10 pathways for each key factor:\n")
  
  for (f in key_factor_names) {
    
    cat("\n====================\n")
    cat(f, "\n")
    cat("====================\n")
    
    tmp <- significant_key_gsea[
      as.character(significant_key_gsea$Factor) == f,
      ,
      drop = FALSE
    ]
    
    tmp <- head(tmp, 10)
    
    print(tmp)
  }
  
  write.csv(
    significant_key_gsea,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_GO_BP_Enrichment.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    significant_key_gsea,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_GO_BP_Enrichment.rds"
    )
  )
  
  cat("\nKey-factor GO-BP interpretation table saved successfully.\n")  

  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  key_factors <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Final.rds"
    )
  )
  
  key_gsea <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Key_Factors_GO_BP_Enrichment.rds"
    )
  )
  
  key_factor_names <- as.character(key_factors$Factor)
  
  go_ids <- unique(key_gsea$GO_ID)
  
  go_annotation <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = go_ids,
    keytype = "GO",
    columns = c("GO", "ONTOLOGY")
  )
  
  go_annotation <- unique(
    go_annotation[
      go_annotation$ONTOLOGY == "BP",
      c("GO"),
      drop = FALSE
    ]
  )
  
  go_annotation$GO_Term <- vapply(
    go_annotation$GO,
    function(go_id) {
      tryCatch(
        AnnotationDbi::Term(
          AnnotationDbi::GOID(
            org.Hs.eg.db,
            keys = go_id,
            keytype = "GO"
          )
        ),
        error = function(e) NA_character_
      )
    },
    character(1)
  )
  
  go_annotation <- unique(
    go_annotation[
      !is.na(go_annotation$GO_Term),
      c("GO", "GO_Term"),
      drop = FALSE
    ]
  )
  
  key_gsea_annotated <- merge(
    key_gsea,
    go_annotation,
    by.x = "GO_ID",
    by.y = "GO",
    all.x = TRUE
  )
  
  key_gsea_annotated <- key_gsea_annotated[
    order(
      as.character(key_gsea_annotated$Factor),
      key_gsea_annotated$FDR
    ),
  ]
  
  cat("\nAnnotated GO terms:", sum(!is.na(key_gsea_annotated$GO_Term)), "\n")
  cat("Total enrichment records:", nrow(key_gsea_annotated), "\n")
  
  print(
    head(
      key_gsea_annotated[
        ,
        c("Factor", "GO_ID", "GO_Term", "P_value", "FDR")
      ],
      20
    ),
    row.names = FALSE
  ) 
  
  if (!requireNamespace("GO.db", quietly = TRUE)) {
    stop("Package 'GO.db' is not installed.")
  }
  
  library(GO.db)
  library(AnnotationDbi)
  
  go_ids <- unique(key_gsea$GO_ID)
  
  go_terms <- AnnotationDbi::select(
    GO.db,
    keys = go_ids,
    keytype = "GOID",
    columns = c("GOID", "TERM", "ONTOLOGY")
  )
  
  go_terms <- unique(go_terms)
  
  go_terms <- go_terms[
    go_terms$ONTOLOGY == "BP",
    c("GOID", "TERM"),
    drop = FALSE
  ]
  
  key_gsea_annotated <- merge(
    key_gsea,
    go_terms,
    by.x = "GO_ID",
    by.y = "GOID",
    all.x = TRUE
  )
  
  key_gsea_annotated <- key_gsea_annotated[
    order(
      as.character(key_gsea_annotated$Factor),
      key_gsea_annotated$FDR
    ),
  ]
  
  cat("\nAnnotated GO terms:", sum(!is.na(key_gsea_annotated$TERM)), "\n")
  cat("Total enrichment records:", nrow(key_gsea_annotated), "\n\n")
  
  print(
    head(
      key_gsea_annotated[
        ,
        c("Factor", "GO_ID", "TERM", "P_value", "FDR")
      ],
      20
    ),
    row.names = FALSE
  ) 
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  key_factor_names <- as.character(key_factors$Factor)
  
  top_pathways <- do.call(
    rbind,
    lapply(
      key_factor_names,
      function(f) {
        
        tmp <- key_gsea_annotated[
          as.character(key_gsea_annotated$Factor) == f,
          ,
          drop = FALSE
        ]
        
        tmp <- head(tmp[order(tmp$FDR), ], 10)
        
        data.frame(
          Factor = f,
          GO_ID = tmp$GO_ID,
          GO_Term = tmp$TERM,
          P_value = tmp$P_value,
          FDR = tmp$FDR,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  write.csv(
    key_gsea_annotated,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_GO_BP_Enrichment_Annotated.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    key_gsea_annotated,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_GO_BP_Enrichment_Annotated.rds"
    )
  )
  
  write.csv(
    top_pathways,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Top10_GO_BP_Pathways_Annotated.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    top_pathways,
    file.path(
      output_dir,
      "MOFA2_Key_Factors_Top10_GO_BP_Pathways_Annotated.rds"
    )
  )
  
  cat("\nFinal MOFA2 GO-BP annotation saved successfully.\n")
  cat("Annotated records:", nrow(key_gsea_annotated), "\n")
  cat("Top pathways table:", nrow(top_pathways), "rows\n")
  cat("Factors covered:", paste(key_factor_names, collapse = ", "), "\n")  
  
  output_dir <- "C:/Users/ibrah/OneDrive/Documents/Results/Proteomics/MOFA2"
  
  factor_evidence <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Final_Key_Factor_Evidence_Table.rds"
    )
  )
  
  key_gsea_annotated <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Key_Factors_GO_BP_Enrichment_Annotated.rds"
    )
  )
  
  trajectory_stats <- readRDS(
    file.path(
      output_dir,
      "MOFA2_Longitudinal_Factor_Statistics.rds"
    )
  )
  
  key_factors <- factor_evidence[
    factor_evidence$Evidence_Count >= 3,
    ,
    drop = FALSE
  ]
  
  cat("\n===== FINAL KEY FACTORS =====\n\n")
  
  print(
    key_factors[
      ,
      c(
        "Factor",
        "RNA_Variance",
        "Protein_Variance",
        "Significant_GO_BP_Pathways",
        "Total_Variance",
        "Acuity_Rho",
        "Acuity_FDR",
        "Acuity_Significant",
        "Longitudinal_Significant",
        "Evidence_Count"
      )
    ],
    row.names = FALSE
  )
  
  cat("\n===== TOP 3 BIOLOGICAL PATHWAYS =====\n")
  
  for (f in as.character(key_factors$Factor)) {
    
    tmp <- key_gsea_annotated[
      as.character(key_gsea_annotated$Factor) == f,
      ,
      drop = FALSE
    ]
    
    tmp <- tmp[order(tmp$FDR), ]
    tmp <- head(tmp, 3)
    
    cat("\n", f, ":\n", sep = "")
    
    print(
      tmp[
        ,
        c("GO_ID", "TERM", "P_value", "FDR")
      ],
      row.names = FALSE
    )
  }
  
  cat("\n===== LONGITUDINAL EVIDENCE =====\n")
  
  long_sig <- trajectory_stats[
    trajectory_stats$Significant_FDR05 == TRUE &
      as.character(trajectory_stats$factor) %in%
      as.character(key_factors$Factor),
    ,
    drop = FALSE
  ]
  
  print(
    long_sig[
      ,
      c(
        "factor",
        "Timepoint_1",
        "Timepoint_2",
        "N",
        "Mean_Difference",
        "FDR",
        "Direction"
      )
    ],
    row.names = FALSE
  )
  
  cat("\n===== FINAL SUMMARY =====\n")
  cat("Key factors:", nrow(key_factors), "\n")
  
  
  
  cat("MOFA GO-BP annotated records:", nrow(key_gsea_annotated), "\n")
  cat("Longitudinal significant comparisons:", nrow(long_sig), "\n") 