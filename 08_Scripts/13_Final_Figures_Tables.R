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
