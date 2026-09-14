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
