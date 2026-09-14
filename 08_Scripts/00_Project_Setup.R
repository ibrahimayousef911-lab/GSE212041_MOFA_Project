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

# Install packages only when missing.
cran_packages <- c("ggplot2", "pheatmap", "dplyr", "tidyr", "readxl", "igraph", "reticulate")
bioc_packages <- c(
  "GEOquery", "DESeq2", "edgeR", "limma", "clusterProfiler",
  "org.Hs.eg.db", "AnnotationDbi", "ReactomePA", "enrichplot",
  "GO.db", "DOSE", "fgsea", "STRINGdb", "MOFA2"
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# Load commonly used packages.
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(pheatmap)
})

# Project identifiers.
geo_id <- "GSE212041"

# Main project directories.
project_dir <- getwd()
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(results_dir, "figures")
tables_dir <- file.path(results_dir, "tables")
supplementary_dir <- file.path(results_dir, "supplementary")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplementary_dir, recursive = TRUE, showWarnings = FALSE)

sessionInfo()
