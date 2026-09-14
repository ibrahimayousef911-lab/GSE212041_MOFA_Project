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
  
