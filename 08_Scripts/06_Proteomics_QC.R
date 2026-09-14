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

  
