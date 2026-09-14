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

# ============================================================
# STEP 1 - DOWNLOAD GEO METADATA
# Project: GSE212041 Multi-Omics COVID-19 Analysis
# ============================================================

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
