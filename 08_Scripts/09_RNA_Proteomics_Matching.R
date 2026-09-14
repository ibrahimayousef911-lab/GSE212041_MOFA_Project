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
