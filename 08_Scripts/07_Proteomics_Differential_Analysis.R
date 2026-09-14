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
  
