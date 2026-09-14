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

  cat("=== Preparing D0 Proteomics IDs for Enrichment ===\n")
  
  d0_sig <- significant_d0_annotated
  
  cat(
    "Significant assays:",
    nrow(d0_sig),
    "\n"
  )
  
  cat(
    "Unique Olink IDs:",
    length(unique(d0_sig$OlinkID)),
    "\n"
  )
  
  cat(
    "Unique UniProt IDs:",
    length(unique(d0_sig$UniProt)),
    "\n"
  )
  
  cat("\n=== Mapping UniProt to Entrez ===\n")
  
  uniprot_ids <- unique(
    d0_sig$UniProt[
      !is.na(d0_sig$UniProt) &
        d0_sig$UniProt != ""
    ]
  )
  
  uniprot_mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = uniprot_ids,
    keytype = "UNIPROT",
    columns = c(
      "UNIPROT",
      "ENTREZID",
      "SYMBOL",
      "GENENAME"
    )
  )
  
  uniprot_mapping <- uniprot_mapping[
    !is.na(uniprot_mapping$ENTREZID),
  ]
  
  uniprot_mapping <- uniprot_mapping[
    !duplicated(uniprot_mapping$UNIPROT),
  ]
  
  cat(
    "Unique UniProt IDs submitted:",
    length(uniprot_ids),
    "\n"
  )
  
  cat(
    "UniProt IDs mapped to Entrez:",
    length(unique(uniprot_mapping$UNIPROT)),
    "\n"
  )
  
  cat(
    "Unique Entrez IDs:",
    length(unique(uniprot_mapping$ENTREZID)),
    "\n"
  )
  
  cat(
    "Mapping percentage:",
    round(
      100 *
        length(unique(uniprot_mapping$UNIPROT)) /
        length(uniprot_ids),
      2
    ),
    "%\n"
  )
  
  cat("\n=== Mapping preview ===\n")
  
  print(
    head(
      uniprot_mapping,
      15
    )
  )
  
  cat("\n=== Saving mapping ===\n")
  
  write.csv(
    uniprot_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_UniProt_to_Entrez_Mapping.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    uniprot_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_UniProt_to_Entrez_Mapping.rds"
    )
  )
  
  cat("\n=== ID mapping completed ===\n") 
  
  cat("=== Loading human gene annotation database ===\n")
  
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install(
      "org.Hs.eg.db",
      ask = FALSE,
      update = FALSE
    )
  }
  
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
  cat(
    "org.Hs.eg.db loaded successfully.\n"
  )
  
  cat(
    "Annotation package version:",
    as.character(
      packageVersion("org.Hs.eg.db")
    ),
    "\n"
  )
  
  cat(
    "AnnotationDbi loaded:",
    requireNamespace("AnnotationDbi", quietly = TRUE),
    "\n"
  ) 
  
  cat("=== D0 Proteomics: UniProt to Entrez Mapping ===\n")
  
  uniprot_ids <- unique(
    d0_sig$UniProt[
      !is.na(d0_sig$UniProt) &
        d0_sig$UniProt != ""
    ]
  )
  
  cat(
    "Unique UniProt IDs submitted:",
    length(uniprot_ids),
    "\n"
  )
  
  uniprot_mapping_raw <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = uniprot_ids,
    keytype = "UNIPROT",
    columns = c(
      "UNIPROT",
      "ENTREZID",
      "SYMBOL",
      "GENENAME"
    )
  )
  
  uniprot_mapping_raw <- uniprot_mapping_raw[
    !is.na(uniprot_mapping_raw$ENTREZID),
  ]
  
  uniprot_mapping <- uniprot_mapping_raw[
    !duplicated(uniprot_mapping_raw$UNIPROT),
  ]
  
  mapped_uniprot <- length(
    unique(uniprot_mapping$UNIPROT)
  )
  
  unique_entrez <- length(
    unique(uniprot_mapping$ENTREZID)
  )
  
  mapping_percentage <- round(
    100 * mapped_uniprot / length(uniprot_ids),
    2
  )
  
  cat(
    "UniProt IDs mapped to Entrez:",
    mapped_uniprot,
    "\n"
  )
  
  cat(
    "Unique Entrez IDs:",
    unique_entrez,
    "\n"
  )
  
  cat(
    "Mapping percentage:",
    mapping_percentage,
    "%\n"
  )
  
  unmapped_uniprot <- setdiff(
    uniprot_ids,
    uniprot_mapping$UNIPROT
  )
  
  cat(
    "Unmapped UniProt IDs:",
    length(unmapped_uniprot),
    "\n"
  )
  
  if (length(unmapped_uniprot) > 0) {
    print(unmapped_uniprot)
  }
  
  cat("\n=== Mapping preview ===\n")
  
  print(uniprot_mapping)
  
  cat("\n=== Saving mapping results ===\n")
  
  write.csv(
    uniprot_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_UniProt_to_Entrez_Mapping.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      UniProt = unmapped_uniprot
    ),
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_Unmapped_UniProt.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    uniprot_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Significant_UniProt_to_Entrez_Mapping.rds"
    )
  )
  
  cat("\n=== UniProt to Entrez mapping completed ===\n")
  
  cat("=== D0 Proteomics: Background Universe Mapping ===\n")
  
  all_uniprot_ids <- unique(
    proteomics_assay_annotation$UniProt[
      !is.na(proteomics_assay_annotation$UniProt) &
        proteomics_assay_annotation$UniProt != ""
    ]
  )
  
  cat(
    "Unique UniProt IDs in all 1429 assays:",
    length(all_uniprot_ids),
    "\n"
  )
  
  background_mapping_raw <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = all_uniprot_ids,
    keytype = "UNIPROT",
    columns = c(
      "UNIPROT",
      "ENTREZID",
      "SYMBOL",
      "GENENAME"
    )
  )
  
  background_mapping_raw <- background_mapping_raw[
    !is.na(background_mapping_raw$ENTREZID),
  ]
  
  background_mapping <- background_mapping_raw[
    !duplicated(background_mapping_raw$UNIPROT),
  ]
  
  background_entrez <- unique(
    background_mapping$ENTREZID
  )
  
  cat(
    "UniProt IDs mapped:",
    length(unique(background_mapping$UNIPROT)),
    "\n"
  )
  
  cat(
    "Unique Entrez IDs in background:",
    length(background_entrez),
    "\n"
  )
  
  cat(
    "Background mapping percentage:",
    round(
      100 * length(unique(background_mapping$UNIPROT)) /
        length(all_uniprot_ids),
      2
    ),
    "%\n"
  )
  
  cat(
    "Unmapped UniProt IDs:",
    length(
      setdiff(
        all_uniprot_ids,
        background_mapping$UNIPROT
      )
    ),
    "\n"
  )
  
  write.csv(
    background_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Background_UniProt_to_Entrez_Mapping.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    background_mapping,
    file.path(
      d0_de_dir,
      "Olink_D0_Background_UniProt_to_Entrez_Mapping.rds"
    )
  )
  
  write.csv(
    data.frame(
      EntrezID = background_entrez
    ),
    file.path(
      d0_de_dir,
      "Olink_D0_Background_Entrez_Universe.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    background_entrez,
    file.path(
      d0_de_dir,
      "Olink_D0_Background_Entrez_Universe.rds"
    )
  )
  
  cat(
    "\n=== Background universe mapping completed ===\n"
  )  
  
  
  cat("=== D0 Proteomics: GO Biological Process Enrichment ===\n")
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    BiocManager::install(
      "clusterProfiler",
      ask = FALSE,
      update = FALSE
    )
  }
  
  library(clusterProfiler)
  
  significant_entrez <- unique(
    uniprot_mapping$ENTREZID
  )
  
  background_entrez <- unique(
    background_entrez
  )
  
  cat(
    "Significant Entrez IDs:",
    length(significant_entrez),
    "\n"
  )
  
  cat(
    "Background Entrez IDs:",
    length(background_entrez),
    "\n"
  )
  
  ego_d0 <- enrichGO(
    gene = significant_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  ego_d0_df <- as.data.frame(ego_d0)
  
  cat(
    "Significant GO-BP pathways:",
    nrow(ego_d0_df),
    "\n"
  )
  
  if (nrow(ego_d0_df) > 0) {
    cat("\n=== Top GO-BP pathways ===\n")
    
    print(
      ego_d0_df[
        order(ego_d0_df$p.adjust),
        c(
          "ID",
          "Description",
          "GeneRatio",
          "BgRatio",
          "Count",
          "pvalue",
          "p.adjust",
          "qvalue",
          "geneID"
        )
      ][1:min(20, nrow(ego_d0_df)), ]
    )
  }
  
  write.csv(
    ego_d0_df,
    file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    ego_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment.rds"
    )
  )
  
  cat(
    "\nSaved:",
    file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment.csv"
    ),
    "\n"
  )
  
  cat(
    "=== GO Biological Process enrichment completed ===\n"
  )  

  cat("=== D0 Proteomics: KEGG Enrichment ===\n")
  
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install(
      "org.Hs.eg.db",
      ask = FALSE,
      update = FALSE
    )
  }
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  ekegg_d0 <- enrichKEGG(
    gene = significant_entrez,
    universe = background_entrez,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
  
  ekegg_d0_df <- as.data.frame(ekegg_d0)
  
  cat(
    "Significant KEGG pathways:",
    nrow(ekegg_d0_df),
    "\n"
  )
  
  if (nrow(ekegg_d0_df) > 0) {
    cat("\n=== Top KEGG pathways ===\n")
    
    print(
      ekegg_d0_df[
        order(ekegg_d0_df$p.adjust),
        c(
          "ID",
          "Description",
          "GeneRatio",
          "BgRatio",
          "Count",
          "pvalue",
          "p.adjust",
          "qvalue",
          "geneID"
        )
      ][1:min(20, nrow(ekegg_d0_df)), ]
    )
  }
  
  write.csv(
    ekegg_d0_df,
    file.path(
      d0_de_dir,
      "Olink_D0_KEGG_Enrichment.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    ekegg_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_KEGG_Enrichment.rds"
    )
  )
  
  cat(
    "\nSaved:",
    file.path(
      d0_de_dir,
      "Olink_D0_KEGG_Enrichment.csv"
    ),
    "\n"
  )
  
  cat("=== KEGG enrichment completed ===\n") 
  
  cat("=== D0 Proteomics: Reactome Enrichment ===\n")
  
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    BiocManager::install(
      "ReactomePA",
      ask = FALSE,
      update = FALSE
    )
  }
  
  library(ReactomePA)
  
  ereactome_d0 <- enrichPathway(
    gene = significant_entrez,
    universe = background_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  
  ereactome_d0_df <- as.data.frame(ereactome_d0)
  
  cat(
    "Significant Reactome pathways:",
    nrow(ereactome_d0_df),
    "\n"
  )
  
  if (nrow(ereactome_d0_df) > 0) {
    cat("\n=== Top Reactome pathways ===\n")
    
    print(
      ereactome_d0_df[
        order(ereactome_d0_df$p.adjust),
        c(
          "ID",
          "Description",
          "GeneRatio",
          "BgRatio",
          "Count",
          "pvalue",
          "p.adjust",
          "qvalue",
          "geneID"
        )
      ][1:min(20, nrow(ereactome_d0_df)), ]
    )
  }
  
  write.csv(
    ereactome_d0_df,
    file.path(
      d0_de_dir,
      "Olink_D0_Reactome_Enrichment.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    ereactome_d0,
    file.path(
      d0_de_dir,
      "Olink_D0_Reactome_Enrichment.rds"
    )
  )
  
  cat(
    "\nSaved:",
    file.path(
      d0_de_dir,
      "Olink_D0_Reactome_Enrichment.csv"
    ),
    "\n"
  )
  
  cat("=== Reactome enrichment completed ===\n")
  
  
  cat("=== D0 Proteomics: GO-BP Visualization ===\n")
  
  library(ggplot2)
  
  go_plot_df <- ego_d0_df
  
  go_plot_df <- go_plot_df[
    order(go_plot_df$p.adjust),
  ]
  
  go_plot_df$Description <- factor(
    go_plot_df$Description,
    levels = rev(go_plot_df$Description)
  )
  
  go_plot_df$minus_log10_padj <- -log10(
    go_plot_df$p.adjust
  )
  
  go_plot <- ggplot(
    go_plot_df,
    aes(
      x = minus_log10_padj,
      y = Description,
      size = Count
    )
  ) +
    geom_point() +
    labs(
      title = "D0 Proteomics: GO Biological Process Enrichment",
      x = "-log10 adjusted P-value",
      y = NULL,
      size = "Protein count"
    ) +
    theme_minimal(base_size = 12)
  
  print(go_plot)
  
  ggsave(
    filename = file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment_DotPlot.png"
    ),
    plot = go_plot,
    width = 9,
    height = 5,
    dpi = 300
  )
  
  write.csv(
    go_plot_df,
    file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment_Plot_Data.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    go_plot,
    file.path(
      d0_de_dir,
      "Olink_D0_GO_BP_Enrichment_DotPlot.rds"
    )
  )
  
  cat(
    "Plot saved:",
    file.exists(
      file.path(
        d0_de_dir,
        "Olink_D0_GO_BP_Enrichment_DotPlot.png"
      )
    ),
    "\n"
  )
  
  cat("=== GO-BP visualization completed ===\n") 
  
  
  cat("=== D0 Proteomics: Final Significant Protein Annotation ===\n")
  
  final_d0_protein_table <- d0_sig
  
  final_d0_protein_table$ENTREZID <- uniprot_mapping$ENTREZID[
    match(
      final_d0_protein_table$UniProt,
      uniprot_mapping$UNIPROT
    )
  ]
  
  final_d0_protein_table$SYMBOL <- uniprot_mapping$SYMBOL[
    match(
      final_d0_protein_table$UniProt,
      uniprot_mapping$UNIPROT
    )
  ]
  
  final_d0_protein_table$GENENAME <- uniprot_mapping$GENENAME[
    match(
      final_d0_protein_table$UniProt,
      uniprot_mapping$UNIPROT
    )
  ]
  
  final_d0_protein_table$Direction <- ifelse(
    final_d0_protein_table$logFC > 0,
    "Upregulated",
    "Downregulated"
  )
  
  chemokine_genes <- c(
    "CXCL11",
    "CCL16",
    "CCL7",
    "TFF2",
    "CXCL10",
    "CCL24",
    "CCL8"
  )
  
  final_d0_protein_table$GO_Chemokine_Enrichment <- 
    final_d0_protein_table$SYMBOL %in% chemokine_genes
  
  final_d0_protein_table <- final_d0_protein_table[
    order(
      final_d0_protein_table$adj.P.Val,
      decreasing = FALSE
    ),
  ]
  
  cat(
    "Total significant proteins:",
    nrow(final_d0_protein_table),
    "\n"
  )
  
  cat(
    "Upregulated:",
    sum(final_d0_protein_table$Direction == "Upregulated"),
    "\n"
  )
  
  cat(
    "Downregulated:",
    sum(final_d0_protein_table$Direction == "Downregulated"),
    "\n"
  )
  
  cat(
    "Proteins contributing to GO chemokine enrichment:",
    sum(final_d0_protein_table$GO_Chemokine_Enrichment),
    "\n"
  )
  
  cat("\n=== Significant proteins ===\n")
  
  print(
    final_d0_protein_table[
      ,
      c(
        "OlinkID",
        "Assay",
        "UniProt",
        "ENTREZID",
        "SYMBOL",
        "logFC",
        "P.Value",
        "adj.P.Val",
        "Direction",
        "GO_Chemokine_Enrichment"
      )
    ]
  )
  
  cat("\n=== GO chemokine proteins ===\n")
  
  print(
    final_d0_protein_table[
      final_d0_protein_table$GO_Chemokine_Enrichment,
      c(
        "OlinkID",
        "Assay",
        "UniProt",
        "ENTREZID",
        "SYMBOL",
        "logFC",
        "adj.P.Val",
        "Direction"
      )
    ]
  )
  
  write.csv(
    final_d0_protein_table,
    file.path(
      d0_de_dir,
      "Olink_D0_Final_42_Significant_Proteins_Annotated.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    final_d0_protein_table,
    file.path(
      d0_de_dir,
      "Olink_D0_Final_42_Significant_Proteins_Annotated.rds"
    )
  )
  
  write.csv(
    final_d0_protein_table[
      final_d0_protein_table$GO_Chemokine_Enrichment,
    ],
    file.path(
      d0_de_dir,
      "Olink_D0_GO_Chemokine_Leading_Proteins.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    final_d0_protein_table[
      final_d0_protein_table$GO_Chemokine_Enrichment,
    ],
    file.path(
      d0_de_dir,
      "Olink_D0_GO_Chemokine_Leading_Proteins.rds"
    )
  )
  
  cat("\n=== Final protein annotation completed ===\n")  
  
  cat("=== D0 Proteomics: Chemokine-Enrichment Protein Plot ===\n")
  
  chemokine_plot_df <- final_d0_protein_table[
    final_d0_protein_table$GO_Chemokine_Enrichment,
    c(
      "SYMBOL",
      "Assay",
      "logFC",
      "adj.P.Val",
      "Direction"
    )
  ]
  
  chemokine_plot_df <- chemokine_plot_df[
    order(chemokine_plot_df$logFC),
  ]
  
  chemokine_plot_df$SYMBOL <- factor(
    chemokine_plot_df$SYMBOL,
    levels = chemokine_plot_df$SYMBOL
  )
  
  chemokine_plot <- ggplot(
    chemokine_plot_df,
    aes(
      x = SYMBOL,
      y = logFC
    )
  ) +
    geom_col() +
    geom_hline(
      yintercept = 0,
      linewidth = 0.5
    ) +
    labs(
      title = "D0 Proteomics: Proteins Contributing to Chemokine-Enriched GO Terms",
      x = "Protein",
      y = "logFC (COVID+ vs COVID−)"
    ) +
    theme_minimal(base_size = 12)
  
  print(chemokine_plot)
  
  ggsave(
    filename = file.path(
      d0_de_dir,
      "Olink_D0_Chemokine_Enrichment_Proteins_logFC.png"
    ),
    plot = chemokine_plot,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  write.csv(
    chemokine_plot_df,
    file.path(
      d0_de_dir,
      "Olink_D0_Chemokine_Enrichment_Proteins_Plot_Data.csv"
    ),
    row.names = FALSE
  )
  
  saveRDS(
    chemokine_plot,
    file.path(
      d0_de_dir,
      "Olink_D0_Chemokine_Enrichment_Proteins_logFC.rds"
    )
  )
  
  cat(
    "Plot saved:",
    file.exists(
      file.path(
        d0_de_dir,
        "Olink_D0_Chemokine_Enrichment_Proteins_logFC.png"
      )
    ),
    "\n"
  )
  
  cat(
    "Proteins plotted:",
    nrow(chemokine_plot_df),
    "\n"
  )
  
  cat("=== Chemokine-enrichment protein visualization completed ===\n")
  
  
