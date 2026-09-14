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
