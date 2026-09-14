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

