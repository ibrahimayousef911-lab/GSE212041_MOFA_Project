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

# STEP 41-A - Restore Previous DE Results
# ============================================================

res_covid_vs_symptomatic <- readRDS(
  "Results/Differential_Expression/COVID_vs_Symptomatic_all_genes.rds"
)

sig_degs <- readRDS(
  "Results/Differential_Expression/COVID_vs_Symptomatic_significant_DEGs.rds"
)

upregulated_genes <- read.csv(
  "Results/Differential_Expression/COVID_vs_Symptomatic_upregulated.csv",
  stringsAsFactors = FALSE
)

downregulated_genes <- read.csv(
  "Results/Differential_Expression/COVID_vs_Symptomatic_downregulated.csv",
  stringsAsFactors = FALSE
)

cat("\n===== DE RESULTS RESTORED =====\n")

cat(
  "Complete DESeq2 results:",
  nrow(res_covid_vs_symptomatic),
  "\n"
)

cat(
  "Significant DEGs:",
  nrow(sig_degs),
  "\n"
)

cat(
  "Upregulated:",
  nrow(upregulated_genes),
  "\n"
)

cat(
  "Downregulated:",
  nrow(downregulated_genes),
  "\n"
)
# ============================================================
# STEP 41-B - Restore GSEA Ranked Gene Lists
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Create Ensembl-ranked list
# ------------------------------------------------------------

gsea_df <- as.data.frame(res_covid_vs_symptomatic)

gsea_df <- gsea_df[
  !is.na(gsea_df$stat) &
    !is.na(gsea_df$gene_id),
]

gsea_df <- gsea_df[
  !duplicated(gsea_df$gene_id),
]

gene_ranking <- gsea_df$stat

names(gene_ranking) <- gsea_df$gene_id

gene_ranking <- sort(
  gene_ranking,
  decreasing = TRUE
)

# ------------------------------------------------------------
# Convert Ensembl -> Entrez
# ------------------------------------------------------------

gene_map_gsea <- bitr(
  names(gene_ranking),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Remove duplicate Entrez IDs
gene_map_gsea <- gene_map_gsea[
  !duplicated(gene_map_gsea$ENTREZID),
]

# Match statistics
gene_map_gsea$stat <- gene_ranking[
  match(
    gene_map_gsea$ENSEMBL,
    names(gene_ranking)
  )
]

# Create Entrez-ranked vector
gene_ranking_entrez <- gene_map_gsea$stat

names(gene_ranking_entrez) <- gene_map_gsea$ENTREZID

gene_ranking_entrez <- sort(
  gene_ranking_entrez,
  decreasing = TRUE
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

cat("\n===== GSEA RANKING RESTORED =====\n")

cat(
  "Ensembl ranked genes:",
  length(gene_ranking),
  "\n"
)

cat(
  "Entrez ranked genes:",
  length(gene_ranking_entrez),
  "\n"
)

cat(
  "Highest statistic:",
  max(gene_ranking_entrez),
  "\n"
)

cat(
  "Lowest statistic:",
  min(gene_ranking_entrez),
  "\n"
)

cat(
  "NA values:",
  sum(is.na(gene_ranking_entrez)),
  "\n"
)

cat(
  "Duplicated Entrez IDs:",
  sum(duplicated(names(gene_ranking_entrez))),
  "\n"
)
# ============================================================
# STEP 41 - Leading-edge Genes
# ============================================================

# ------------------------------------------------------------
# Select top Reactome pathways
# ------------------------------------------------------------

selected_reactome <- c(
  "R-HSA-69481",
  "R-HSA-69620",
  "R-HSA-69306"
)

reactome_selected <- gsea_reactome_df[
  gsea_reactome_df$ID %in% selected_reactome,
]

reactome_selected <- reactome_selected[
  order(reactome_selected$p.adjust),
]

# ------------------------------------------------------------
# Display pathway statistics
# ------------------------------------------------------------

cat("\n===== SELECTED REACTOME PATHWAYS =====\n\n")

print(
  reactome_selected[
    ,
    c(
      "ID",
      "Description",
      "setSize",
      "enrichmentScore",
      "NES",
      "pvalue",
      "p.adjust"
    )
  ]
)

# ------------------------------------------------------------
# Extract core enrichment genes
# ------------------------------------------------------------

leading_edge_list <- list()

for (i in seq_len(nrow(reactome_selected))) {
  
  pathway_id <- reactome_selected$ID[i]
  
  pathway_name <- reactome_selected$Description[i]
  
  genes <- unlist(
    strsplit(
      reactome_selected$core_enrichment[i],
      "/"
    )
  )
  
  leading_edge_list[[pathway_id]] <- genes
  
  cat("\n============================================\n")
  cat("Pathway:", pathway_name, "\n")
  cat("Pathway ID:", pathway_id, "\n")
  cat("Leading-edge genes:", length(genes), "\n")
  cat("============================================\n")
  
  print(genes)
}

# ============================================================
# STEP 42 — Convert Leading-edge Entrez IDs to Gene Symbols
# ============================================================

library(org.Hs.eg.db)
library(AnnotationDbi)

# ------------------------------------------------------------
# Combine all leading-edge genes
# ------------------------------------------------------------

leading_edge_df <- data.frame(
  EntrezID = unique(unlist(leading_edge_list)),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Map Entrez IDs -> Gene Symbols
# ------------------------------------------------------------

gene_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = leading_edge_df$EntrezID,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "ENTREZID"
)

gene_mapping <- gene_mapping[
  !duplicated(gene_mapping$ENTREZID),
]

# ------------------------------------------------------------
# Merge mapping
# ------------------------------------------------------------

leading_edge_df <- merge(
  leading_edge_df,
  gene_mapping,
  by.x = "EntrezID",
  by.y = "ENTREZID",
  all.x = TRUE
)

# ------------------------------------------------------------
# Add pathway membership
# ------------------------------------------------------------

leading_edge_df$G2M_Checkpoints <- leading_edge_df$EntrezID %in%
  leading_edge_list[["R-HSA-69481"]]

leading_edge_df$Cell_Cycle_Checkpoints <- leading_edge_df$EntrezID %in%
  leading_edge_list[["R-HSA-69620"]]

leading_edge_df$DNA_Replication <- leading_edge_df$EntrezID %in%
  leading_edge_list[["R-HSA-69306"]]

# Number of pathways containing each gene
leading_edge_df$N_Pathways <- rowSums(
  leading_edge_df[
    ,
    c(
      "G2M_Checkpoints",
      "Cell_Cycle_Checkpoints",
      "DNA_Replication"
    )
  ]
)

# ------------------------------------------------------------
# Sort: genes shared by all pathways first
# ------------------------------------------------------------

leading_edge_df <- leading_edge_df[
  order(
    -leading_edge_df$N_Pathways,
    leading_edge_df$SYMBOL
  ),
]

# ------------------------------------------------------------
# Display summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 42 — LEADING-EDGE GENE SUMMARY\n")
cat("============================================\n\n")

cat("Unique leading-edge Entrez IDs:",
    nrow(leading_edge_df), "\n")

cat(
  "Mapped Gene Symbols:",
  sum(!is.na(leading_edge_df$SYMBOL)),
  "\n"
)

cat("\nGenes shared by ALL 3 pathways:\n")

shared_all <- leading_edge_df$SYMBOL[
  leading_edge_df$N_Pathways == 3 &
    !is.na(leading_edge_df$SYMBOL)
]

print(shared_all)

cat("\nNumber shared by all 3:",
    length(shared_all), "\n")

cat("\nTop leading-edge genes by pathway overlap:\n")

print(
  leading_edge_df[
    ,
    c(
      "EntrezID",
      "SYMBOL",
      "G2M_Checkpoints",
      "Cell_Cycle_Checkpoints",
      "DNA_Replication",
      "N_Pathways"
    )
  ],
  row.names = FALSE
)

# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

dir.create(
  "Results/Enrichment/GSEA/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  leading_edge_df,
  "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.csv",
  row.names = FALSE
)

write.csv(
  data.frame(GeneSymbol = shared_all),
  "Results/Enrichment/GSEA/Shared_Leading_Edge_Genes_Top3.csv",
  row.names = FALSE
)

saveRDS(
  leading_edge_df,
  "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.rds"
)

cat("\n============================================\n")
cat("FILES SAVED SUCCESSFULLY\n")
cat("============================================\n")

cat(
  "Leading-edge table: ",
  file.exists(
    "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.csv"
  ),
  "\n"
)

cat(
  "Shared genes table: ",
  file.exists(
    "Results/Enrichment/GSEA/Shared_Leading_Edge_Genes_Top3.csv"
  ),
  "\n"
)

# ============================================================
# STEP 43 — Core Leading-edge Genes + DESeq2 Statistics
# ============================================================

core_genes <- shared_all

# ------------------------------------------------------------
# Extract DESeq2 results
# ------------------------------------------------------------

core_de <- res_covid_vs_symptomatic[
  res_covid_vs_symptomatic$SYMBOL %in% core_genes |
    rownames(res_covid_vs_symptomatic) %in% core_genes,
]

# ------------------------------------------------------------
# If SYMBOL column is unavailable, map rownames
# ------------------------------------------------------------

if (nrow(core_de) == 0) {
  
  core_entrez <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = core_genes,
    columns = c("SYMBOL"),
    keytype = "SYMBOL"
  )
  
  core_de <- res_covid_vs_symptomatic[
    rownames(res_covid_vs_symptomatic) %in%
      core_entrez$SYMBOL,
  ]
}

# ------------------------------------------------------------
# Add Gene Symbol from rownames when necessary
# ------------------------------------------------------------

if (!"SYMBOL" %in% colnames(core_de)) {
  core_de$SYMBOL <- rownames(core_de)
}

# ------------------------------------------------------------
# Select useful columns
# ------------------------------------------------------------

core_de <- core_de[
  ,
  intersect(
    c(
      "SYMBOL",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    ),
    colnames(core_de)
  ),
  drop = FALSE
]

# ------------------------------------------------------------
# Remove duplicates
# ------------------------------------------------------------

core_de <- core_de[
  !duplicated(core_de$SYMBOL),
]

# ------------------------------------------------------------
# Add regulation status
# ------------------------------------------------------------

core_de$Regulation <- ifelse(
  !is.na(core_de$padj) &
    core_de$padj < 0.05 &
    core_de$log2FoldChange > 1,
  "Upregulated",
  ifelse(
    !is.na(core_de$padj) &
      core_de$padj < 0.05 &
      core_de$log2FoldChange < -1,
    "Downregulated",
    "Not significant"
  )
)

# ------------------------------------------------------------
# Sort by log2FC
# ------------------------------------------------------------

core_de <- core_de[
  order(
    -core_de$log2FoldChange
  ),
]

# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 43 — CORE LEADING-EDGE + DESEQ2\n")
cat("============================================\n\n")

cat(
  "Core genes:",
  length(core_genes),
  "\n"
)

cat(
  "Genes found in DESeq2:",
  nrow(core_de),
  "\n"
)

cat("\nRegulation summary:\n")
print(table(core_de$Regulation))

cat("\nCore genes with DESeq2 statistics:\n\n")

print(
  core_de,
  row.names = FALSE
)

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  core_de,
  "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.csv",
  row.names = FALSE
)

saveRDS(
  core_de,
  "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.rds"
)

cat("\n============================================\n")
cat("FILES SAVED\n")
cat("============================================\n")

cat(
  "CSV:",
  file.exists(
    "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.csv"
  ),
  "\n"
)

cat(
  "RDS:",
  file.exists(
    "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.rds"
  ),
  "\n"
)
# ============================================================
# STEP 43 FIX — Map DESeq2 Ensembl IDs to Gene Symbols
# ============================================================

library(org.Hs.eg.db)
library(AnnotationDbi)

# ------------------------------------------------------------
# Check DESeq2 structure
# ------------------------------------------------------------

cat("\n============================================\n")
cat("DESeq2 RESULT STRUCTURE\n")
cat("============================================\n\n")

cat("Number of DESeq2 genes:",
    nrow(res_covid_vs_symptomatic), "\n\n")

cat("DESeq2 columns:\n")
print(colnames(res_covid_vs_symptomatic))

cat("\nFirst 5 rownames:\n")
print(head(rownames(res_covid_vs_symptomatic)))

# ------------------------------------------------------------
# Create Ensembl ID column
# ------------------------------------------------------------

de_ids <- rownames(res_covid_vs_symptomatic)

# Remove version numbers if present
# Example: ENSG00000123456.12 -> ENSG00000123456

de_ensembl <- sub(
  "\\..*$",
  "",
  de_ids
)

# ------------------------------------------------------------
# Map Ensembl -> Gene Symbol
# ------------------------------------------------------------

de_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(de_ensembl),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
)

# Keep one mapping per Ensembl ID
de_mapping <- de_mapping[
  !duplicated(de_mapping$ENSEMBL),
]

# ------------------------------------------------------------
# Add Ensembl + Symbol
# ------------------------------------------------------------

de_table <- as.data.frame(res_covid_vs_symptomatic)

de_table$ENSEMBL <- de_ensembl

de_table$SYMBOL <- de_mapping$SYMBOL[
  match(
    de_table$ENSEMBL,
    de_mapping$ENSEMBL
  )
]

# ------------------------------------------------------------
# Extract the 54 core genes
# ------------------------------------------------------------

core_de <- de_table[
  !is.na(de_table$SYMBOL) &
    de_table$SYMBOL %in% core_genes,
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Remove duplicate symbols
# ------------------------------------------------------------

core_de <- core_de[
  !duplicated(core_de$SYMBOL),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Add regulation status
# ------------------------------------------------------------

core_de$Regulation <- ifelse(
  !is.na(core_de$padj) &
    core_de$padj < 0.05 &
    core_de$log2FoldChange > 1,
  "Upregulated",
  ifelse(
    !is.na(core_de$padj) &
      core_de$padj < 0.05 &
      core_de$log2FoldChange < -1,
    "Downregulated",
    "Not significant"
  )
)

# ------------------------------------------------------------
# Sort by log2FC
# ------------------------------------------------------------

core_de <- core_de[
  order(
    -core_de$log2FoldChange
  ),
]

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 43 FIX — CORE GENE RESULTS\n")
cat("============================================\n\n")

cat(
  "Core genes from GSEA:",
  length(core_genes),
  "\n"
)

cat(
  "Core genes found in DESeq2:",
  nrow(core_de),
  "\n"
)

cat("\nRegulation summary:\n")
print(table(core_de$Regulation))

cat("\nCore genes + DESeq2 statistics:\n\n")

print(
  core_de[
    ,
    c(
      "ENSEMBL",
      "SYMBOL",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj",
      "Regulation"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

# ------------------------------------------------------------
# Save corrected results
# ------------------------------------------------------------

write.csv(
  core_de,
  "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.csv",
  row.names = FALSE
)

saveRDS(
  core_de,
  "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.rds"
)

cat("\n============================================\n")
cat("FILES SAVED\n")
cat("============================================\n")

cat(
  "CSV:",
  file.exists(
    "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.csv"
  ),
  "\n"
)

cat(
  "RDS:",
  file.exists(
    "Results/Enrichment/GSEA/Core_Leading_Edge_Genes_DESeq2.rds"
  ),
  "\n"
)
# ============================================================
# STEP 44 — Core Cell-Cycle Signature
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# Select significant core genes
# ------------------------------------------------------------

core_sig <- core_de[
  core_de$Regulation == "Upregulated",
  ,
  drop = FALSE
]

# Sort by log2FC
core_sig <- core_sig[
  order(core_sig$log2FoldChange),
  ,
  drop = FALSE
]

# ------------------------------------------------------------
# Print
# ------------------------------------------------------------

cat("\n============================================\n")
cat("STEP 44 — CORE CELL-CYCLE SIGNATURE\n")
cat("============================================\n\n")

cat(
  "Significant core genes:",
  nrow(core_sig),
  "\n\n"
)

print(
  core_sig[
    ,
    c(
      "SYMBOL",
      "baseMean",
      "log2FoldChange",
      "pvalue",
      "padj"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

# ------------------------------------------------------------
# Save table
# ------------------------------------------------------------

write.csv(
  core_sig,
  "Results/Enrichment/GSEA/Core_Significant_Genes.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Create plot
# ------------------------------------------------------------

p_core_signature <- ggplot(
  core_sig,
  aes(
    x = log2FoldChange,
    y = SYMBOL
  )
) +
  geom_point(
    size = 4
  ) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  labs(
    title = "Core Cell-Cycle / DNA-Replication Signature",
    subtitle = "Significant leading-edge genes from top Reactome GSEA pathways",
    x = "log2 Fold Change",
    y = "Gene"
  ) +
  theme_minimal(base_size = 13)

print(p_core_signature)

# ------------------------------------------------------------
# Save plot
# ------------------------------------------------------------

ggsave(
  "Results/Enrichment/GSEA/Figures/Core_Cell_Cycle_Signature.png",
  p_core_signature,
  width = 9,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

cat("\n============================================\n")
cat("FILES SAVED\n")
cat("============================================\n")

cat(
  "Table:",
  file.exists(
    "Results/Enrichment/GSEA/Core_Significant_Genes.csv"
  ),
  "\n"
)

cat(
  "Plot:",
  file.exists(
    "Results/Enrichment/GSEA/Figures/Core_Cell_Cycle_Signature.png"
  ),
  "\n"
)
# ============================================================
# STEP 45 — Core Gene–Pathway Network
# ============================================================

cat("\n============================================================\n")
cat("STEP 45 — CORE GENE–PATHWAY NETWORK\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load required objects
# ------------------------------------------------------------

core_sig <- read.csv(
  "Results/Enrichment/GSEA/Core_Significant_Genes.csv",
  stringsAsFactors = FALSE
)

leading_edge <- read.csv(
  "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Check column names
# ------------------------------------------------------------

cat("Core signature columns:\n")
print(colnames(core_sig))

cat("\nLeading-edge columns:\n")
print(colnames(leading_edge))

# ------------------------------------------------------------
# 3. Keep only the 12 significant core genes
# ------------------------------------------------------------

core_genes <- unique(core_sig$SYMBOL)

cat("\nNumber of core genes:", length(core_genes), "\n")
print(core_genes)

# ------------------------------------------------------------
# 4. Identify pathway column
# ------------------------------------------------------------

pathway_col <- intersect(
  c("Description", "Pathway", "pathway", "ID"),
  colnames(leading_edge)
)

gene_col <- intersect(
  c("SYMBOL", "gene_symbol", "Gene", "gene"),
  colnames(leading_edge)
)

cat("\nDetected pathway column:", pathway_col[1], "\n")
cat("Detected gene column:", gene_col[1], "\n")

# ------------------------------------------------------------
# 5. Build gene–pathway network
# ------------------------------------------------------------

network_df <- leading_edge[
  leading_edge[[gene_col[1]]] %in% core_genes,
  c(pathway_col[1], gene_col[1])
]

colnames(network_df) <- c(
  "Pathway",
  "Gene"
)

network_df <- unique(network_df)

# ------------------------------------------------------------
# 6. Add gene regulation information
# ------------------------------------------------------------

network_df <- merge(
  network_df,
  core_sig[, c("SYMBOL", "log2FoldChange", "padj", "Regulation")],
  by.x = "Gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

network_df <- network_df[
  order(network_df$Pathway, -network_df$log2FoldChange),
]

# ------------------------------------------------------------
# 7. Print network
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE–PATHWAY NETWORK\n")
cat("============================================================\n\n")

print(network_df)

# ------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------

cat("\nNumber of network edges:", nrow(network_df), "\n")
cat("Number of genes:", length(unique(network_df$Gene)), "\n")
cat("Number of pathways:", length(unique(network_df$Pathway)), "\n")

# ------------------------------------------------------------
# 9. Save network table
# ------------------------------------------------------------

dir.create(
  "Results/Network",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  network_df,
  "Results/Network/Core_Gene_Pathway_Network.csv",
  row.names = FALSE
)

saveRDS(
  network_df,
  "Results/Network/Core_Gene_Pathway_Network.rds"
)

# ------------------------------------------------------------
# 10. Validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n")

cat(
  "CSV:",
  file.exists(
    "Results/Network/Core_Gene_Pathway_Network.csv"
  ),
  "\n"
)

cat(
  "RDS:",
  file.exists(
    "Results/Network/Core_Gene_Pathway_Network.rds"
  ),
  "\n"
)

cat("\nSTEP 45 completed successfully.\n")
# ============================================================
# STEP 45 — FIX
# Core Gene–Pathway Network
# ============================================================

cat("\n============================================================\n")
cat("STEP 45 — CORE GENE–PATHWAY NETWORK\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

core_sig <- read.csv(
  "Results/Enrichment/GSEA/Core_Significant_Genes.csv",
  stringsAsFactors = FALSE
)

leading_edge <- read.csv(
  "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 2. Core genes
# ------------------------------------------------------------

core_genes <- unique(core_sig$SYMBOL)

cat("Number of core genes:", length(core_genes), "\n")
print(core_genes)

# ------------------------------------------------------------
# 3. Reactome pathway columns
# ------------------------------------------------------------

pathway_cols <- c(
  "G2M_Checkpoints",
  "Cell_Cycle_Checkpoints",
  "DNA_Replication"
)

cat("\nPathway columns:\n")
print(pathway_cols)

# Check that all columns exist
cat("\nAll pathway columns exist:",
    all(pathway_cols %in% colnames(leading_edge)),
    "\n")

# ------------------------------------------------------------
# 4. Create gene–pathway edges
# ------------------------------------------------------------

network_list <- lapply(pathway_cols, function(pw) {
  
  x <- leading_edge[
    leading_edge$SYMBOL %in% core_genes &
      !is.na(leading_edge[[pw]]) &
      leading_edge[[pw]] != FALSE &
      leading_edge[[pw]] != 0,
    c("SYMBOL", pw)
  ]
  
  if (nrow(x) == 0) {
    return(NULL)
  }
  
  data.frame(
    Gene = x$SYMBOL,
    Pathway = pw,
    stringsAsFactors = FALSE
  )
})

network_df <- do.call(rbind, network_list)

# Remove duplicates
network_df <- unique(network_df)

# ------------------------------------------------------------
# 5. Add DESeq2 information
# ------------------------------------------------------------

network_df <- merge(
  network_df,
  core_sig[, c(
    "SYMBOL",
    "log2FoldChange",
    "padj",
    "Regulation"
  )],
  by.x = "Gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

# Sort
network_df <- network_df[
  order(network_df$Pathway, -network_df$log2FoldChange),
]

# ------------------------------------------------------------
# 6. Clean pathway names
# ------------------------------------------------------------

network_df$Pathway <- factor(
  network_df$Pathway,
  levels = pathway_cols
)

network_df$Pathway <- as.character(network_df$Pathway)

network_df$Pathway[
  network_df$Pathway == "G2M_Checkpoints"
] <- "G2/M Checkpoints"

network_df$Pathway[
  network_df$Pathway == "Cell_Cycle_Checkpoints"
] <- "Cell Cycle Checkpoints"

network_df$Pathway[
  network_df$Pathway == "DNA_Replication"
] <- "DNA Replication"

# ------------------------------------------------------------
# 7. Print network
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GENE–PATHWAY NETWORK\n")
cat("============================================================\n\n")

print(network_df)

# ------------------------------------------------------------
# 8. Summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("NETWORK SUMMARY\n")
cat("============================================================\n\n")

cat("Number of edges:",
    nrow(network_df),
    "\n")

cat("Number of unique genes:",
    length(unique(network_df$Gene)),
    "\n")

cat("Number of pathways:",
    length(unique(network_df$Pathway)),
    "\n")

cat("\nEdges per pathway:\n")
print(table(network_df$Pathway))

# ------------------------------------------------------------
# 9. Save results
# ------------------------------------------------------------

dir.create(
  "Results/Network",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  network_df,
  "Results/Network/Core_Gene_Pathway_Network.csv",
  row.names = FALSE
)

saveRDS(
  network_df,
  "Results/Network/Core_Gene_Pathway_Network.rds"
)

# ------------------------------------------------------------
# 10. Validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n")

cat(
  "CSV:",
  file.exists(
    "Results/Network/Core_Gene_Pathway_Network.csv"
  ),
  "\n"
)

cat(
  "RDS:",
  file.exists(
    "Results/Network/Core_Gene_Pathway_Network.rds"
  ),
  "\n"
)

cat("\nSTEP 45 FIX completed successfully.\n")


# ============================================================
# STEP 46 — PREPARE CYTOSCAPE NETWORK FILES
# ============================================================

cat("\n============================================================\n")
cat("STEP 46 — CYTOSCAPE NETWORK PREPARATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load network
# ------------------------------------------------------------

network_df <- readRDS(
  "Results/Network/Core_Gene_Pathway_Network.rds"
)

# ------------------------------------------------------------
# 2. Create EDGE table
# ------------------------------------------------------------

edges <- network_df[, c(
  "Gene",
  "Pathway"
)]

colnames(edges) <- c(
  "source",
  "target"
)

edges$interaction <- "core_gene_in_pathway"

# ------------------------------------------------------------
# 3. Create GENE nodes
# ------------------------------------------------------------

gene_nodes <- unique(
  network_df[, c(
    "Gene",
    "log2FoldChange",
    "padj",
    "Regulation"
  )]
)

gene_nodes$NodeType <- "Gene"

colnames(gene_nodes)[1] <- "Node"

# ------------------------------------------------------------
# 4. Create PATHWAY nodes
# ------------------------------------------------------------

pathway_nodes <- data.frame(
  Node = unique(network_df$Pathway),
  log2FoldChange = NA,
  padj = NA,
  Regulation = NA,
  NodeType = "Pathway",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 5. Combine nodes
# ------------------------------------------------------------

nodes <- rbind(
  gene_nodes,
  pathway_nodes
)

# Remove duplicate nodes if any
nodes <- unique(nodes)

# ------------------------------------------------------------
# 6. Save Cytoscape files
# ------------------------------------------------------------

dir.create(
  "Results/Network/Cytoscape",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  edges,
  "Results/Network/Cytoscape/Core_Network_Edges.csv",
  row.names = FALSE
)

write.csv(
  nodes,
  "Results/Network/Cytoscape/Core_Network_Nodes.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n\n")

cat("Number of edges:", nrow(edges), "\n")
cat("Number of nodes:", nrow(nodes), "\n")

cat("\nNode types:\n")
print(table(nodes$NodeType))

cat("\nEdges file:",
    file.exists(
      "Results/Network/Cytoscape/Core_Network_Edges.csv"
    ),
    "\n")

cat("Nodes file:",
    file.exists(
      "Results/Network/Cytoscape/Core_Network_Nodes.csv"
    ),
    "\n")

cat("\nSTEP 46 completed successfully.\n")

# ============================================================
# STEP 47 — FIX
# CORE GENE–PATHWAY NETWORK VISUALIZATION
# ============================================================

cat("\n============================================================\n")
cat("STEP 47 — NETWORK VISUALIZATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(igraph)
library(ggraph)
library(ggplot2)

# ------------------------------------------------------------
# 2. Load network files
# ------------------------------------------------------------

edges <- read.csv(
  "Results/Network/Cytoscape/Core_Network_Edges.csv",
  stringsAsFactors = FALSE
)

nodes <- read.csv(
  "Results/Network/Cytoscape/Core_Network_Nodes.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 3. Create TRUE bipartite node type
# ------------------------------------------------------------

nodes$type <- nodes$NodeType == "Pathway"

# ------------------------------------------------------------
# 4. Create igraph object
# ------------------------------------------------------------

graph <- graph_from_data_frame(
  d = edges[, c("source", "target")],
  vertices = nodes,
  directed = FALSE
)

# ------------------------------------------------------------
# 5. Add readable node labels
# ------------------------------------------------------------

V(graph)$type_label <- ifelse(
  V(graph)$type,
  "Pathway",
  "Gene"
)

# ------------------------------------------------------------
# 6. Add node sizes
# ------------------------------------------------------------

logfc_values <- nodes$log2FoldChange[
  match(V(graph)$name, nodes$Node)
]

V(graph)$node_size <- ifelse(
  V(graph)$type_label == "Gene",
  4 + abs(logfc_values) * 4,
  9
)

# ------------------------------------------------------------
# 7. Create network plot
# ------------------------------------------------------------

network_plot <- ggraph(
  graph,
  layout = "bipartite"
) +
  geom_edge_link(
    linewidth = 0.7,
    alpha = 0.6
  ) +
  geom_node_point(
    aes(
      size = node_size,
      shape = type_label
    )
  ) +
  geom_node_text(
    aes(
      label = name
    ),
    repel = TRUE,
    size = 4
  ) +
  scale_shape_manual(
    values = c(
      Gene = 16,
      Pathway = 17
    )
  ) +
  scale_size_identity() +
  labs(
    title = "Core Gene–Pathway Network",
    subtitle = "12 significant core genes connected to top Reactome pathways",
    shape = "Node type"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 12
    ),
    legend.position = "bottom"
  )

# ------------------------------------------------------------
# 8. Display
# ------------------------------------------------------------

print(network_plot)

# ------------------------------------------------------------
# 9. Save figure
# ------------------------------------------------------------

dir.create(
  "Results/Network/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "Results/Network/Figures/Core_Gene_Pathway_Network.png",
  network_plot,
  width = 12,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n\n")

cat("Nodes:", vcount(graph), "\n")
cat("Edges:", ecount(graph), "\n")

cat(
  "Gene nodes:",
  sum(V(graph)$type_label == "Gene"),
  "\n"
)

cat(
  "Pathway nodes:",
  sum(V(graph)$type_label == "Pathway"),
  "\n"
)

cat(
  "Plot:",
  file.exists(
    "Results/Network/Figures/Core_Gene_Pathway_Network.png"
  ),
  "\n"
)

cat("\nSTEP 47 FIX completed successfully.\n")
# ============================================================
# STEP 48 — IMPROVED CORE NETWORK FIGURE
# ============================================================

cat("\n============================================================\n")
cat("STEP 48 — IMPROVED NETWORK FIGURE\n")
cat("============================================================\n\n")

library(igraph)
library(ggraph)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load network
# ------------------------------------------------------------

network_df <- readRDS(
  "Results/Network/Core_Gene_Pathway_Network.rds"
)

# ------------------------------------------------------------
# 2. Create edges
# ------------------------------------------------------------

edges <- network_df[, c("Gene", "Pathway")]
colnames(edges) <- c("from", "to")

# ------------------------------------------------------------
# 3. Create nodes
# ------------------------------------------------------------

gene_nodes <- unique(
  network_df[, c(
    "Gene",
    "log2FoldChange",
    "padj",
    "Regulation"
  )]
)

colnames(gene_nodes)[1] <- "name"
gene_nodes$NodeType <- "Gene"

pathway_nodes <- data.frame(
  name = unique(network_df$Pathway),
  log2FoldChange = NA,
  padj = NA,
  Regulation = NA,
  NodeType = "Pathway",
  stringsAsFactors = FALSE
)

nodes <- rbind(
  gene_nodes,
  pathway_nodes
)

# ------------------------------------------------------------
# 4. Create graph
# ------------------------------------------------------------

graph <- graph_from_data_frame(
  edges,
  vertices = nodes,
  directed = FALSE
)

# TRUE = pathway
V(graph)$type <- V(graph)$NodeType == "Pathway"

# ------------------------------------------------------------
# 5. Node size
# ------------------------------------------------------------

V(graph)$size_value <- ifelse(
  V(graph)$NodeType == "Gene",
  5 + abs(V(graph)$log2FoldChange) * 3,
  10
)

# ------------------------------------------------------------
# 6. Improved plot
# ------------------------------------------------------------

network_plot_improved <- ggraph(
  graph,
  layout = "bipartite"
) +
  
  geom_edge_link(
    linewidth = 0.5,
    alpha = 0.35
  ) +
  
  geom_node_point(
    aes(
      size = size_value,
      shape = NodeType
    )
  ) +
  
  geom_node_text(
    aes(
      label = name
    ),
    repel = TRUE,
    size = 3.8,
    fontface = ifelse(
      V(graph)$NodeType == "Pathway",
      "bold",
      "plain"
    )
  ) +
  
  scale_size_identity() +
  
  scale_shape_manual(
    values = c(
      Gene = 16,
      Pathway = 17
    )
  ) +
  
  labs(
    title = "Core Gene–Pathway Network",
    subtitle = paste(
      "12 significant core genes across 3 Reactome pathways"
    ),
    shape = "Node type"
  ) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      size = 20,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5
    ),
    legend.position = "bottom",
    plot.margin = margin(
      20, 40, 20, 40
    )
  )

# ------------------------------------------------------------
# 7. Display
# ------------------------------------------------------------

print(network_plot_improved)

# ------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------

ggsave(
  "Results/Network/Figures/Core_Gene_Pathway_Network_Improved.png",
  network_plot_improved,
  width = 14,
  height = 10,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VALIDATION\n")
cat("============================================================\n\n")

cat("Nodes:", vcount(graph), "\n")
cat("Edges:", ecount(graph), "\n")

cat(
  "Improved plot:",
  file.exists(
    "Results/Network/Figures/Core_Gene_Pathway_Network_Improved.png"
  ),
  "\n"
)

cat("\nSTEP 48 completed successfully.\n")

# ============================================================
# STEP 49 — FINAL CYTOSCAPE FILE VALIDATION
# ============================================================

cat("\n============================================================\n")
cat("STEP 49 — CYTOSCAPE FILE VALIDATION\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Load files
# ------------------------------------------------------------

edges_cyto <- read.csv(
  "Results/Network/Cytoscape/Core_Network_Edges.csv",
  stringsAsFactors = FALSE
)

nodes_cyto <- read.csv(
  "Results/Network/Cytoscape/Core_Network_Nodes.csv",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Basic dimensions
# ------------------------------------------------------------

cat("Edges:", nrow(edges_cyto), "\n")
cat("Nodes:", nrow(nodes_cyto), "\n")

# ------------------------------------------------------------
# Node types
# ------------------------------------------------------------

cat("\nNode types:\n")
print(table(nodes_cyto$NodeType))

# ------------------------------------------------------------
# Check genes
# ------------------------------------------------------------

gene_nodes <- nodes_cyto[
  nodes_cyto$NodeType == "Gene",
]

cat("\nGene nodes:", nrow(gene_nodes), "\n")

cat("\nCore genes:\n")
print(gene_nodes$Node)

# ------------------------------------------------------------
# Check pathways
# ------------------------------------------------------------

pathway_nodes <- nodes_cyto[
  nodes_cyto$NodeType == "Pathway",
]

cat("\nPathway nodes:", nrow(pathway_nodes), "\n")

cat("\nPathways:\n")
print(pathway_nodes$Node)

# ------------------------------------------------------------
# Check edge integrity
# ------------------------------------------------------------

cat("\nAll edge sources are valid:",
    all(edges_cyto$source %in% nodes_cyto$Node),
    "\n")

cat("All edge targets are valid:",
    all(edges_cyto$target %in% nodes_cyto$Node),
    "\n")

# ------------------------------------------------------------
# Check duplicate edges
# ------------------------------------------------------------

cat(
  "Duplicate edges:",
  sum(duplicated(edges_cyto[, c("source", "target")])),
  "\n"
)

# ------------------------------------------------------------
# Final validation
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL VALIDATION\n")
cat("============================================================\n\n")

cat(
  "Edges file:",
  file.exists(
    "Results/Network/Cytoscape/Core_Network_Edges.csv"
  ),
  "\n"
)

cat(
  "Nodes file:",
  file.exists(
    "Results/Network/Cytoscape/Core_Network_Nodes.csv"
  ),
  "\n"
)

cat("\nSTEP 49 completed successfully.\n")
# STEP 52 — Prepare Core Genes for STRING

string_genes <- core_genes$SYMBOL

cat("Genes prepared for STRING:", length(string_genes), "\n")
print(string_genes)

# STEP 53 — STRING/PPI analysis

if (!requireNamespace("STRINGdb", quietly = TRUE)) {
  BiocManager::install("STRINGdb", ask = FALSE, update = FALSE)
}

library(STRINGdb)

string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400,
  input_directory = ""
)

string_mapped <- string_db$map(
  data.frame(gene = string_genes),
  "gene",
  removeUnmappedRows = TRUE
)

cat("Input genes:", length(string_genes), "\n")
cat("Mapped genes:", nrow(string_mapped), "\n")

print(string_mapped[, c("gene", "STRING_id")])

# STEP 54 — Extract STRING interactions

string_ids <- string_mapped$STRING_id

ppi <- string_db$get_interactions(string_ids)

cat("Number of PPI interactions:", nrow(ppi), "\n")
print(head(ppi))


# STEP 55 — Save STRING PPI network

dir.create("Results/Network/STRING", recursive = TRUE, showWarnings = FALSE)

# Add DESeq2 information to STRING nodes
string_nodes <- string_mapped

string_nodes <- merge(
  string_nodes,
  core_genes[, c("SYMBOL", "log2FoldChange", "padj", "Regulation")],
  by.x = "gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

# Save interactions
write.csv(
  ppi,
  "Results/Network/STRING/Core_Genes_STRING_PPI.csv",
  row.names = FALSE
)

# Save nodes
write.csv(
  string_nodes,
  "Results/Network/STRING/Core_Genes_STRING_Nodes.csv",
  row.names = FALSE
)

cat("PPI edges saved:", file.exists(
  "Results/Network/STRING/Core_Genes_STRING_PPI.csv"
), "\n")

cat("PPI nodes saved:", file.exists(
  "Results/Network/STRING/Core_Genes_STRING_Nodes.csv"
), "\n")

cat("Edges:", nrow(ppi), "\n")
cat("Nodes:", nrow(string_nodes), "\n")

# STEP 56 — Calculate PPI degree

library(igraph)

ppi_graph <- graph_from_data_frame(
  ppi[, c("from_name", "to_name", "combined_score")],
  directed = FALSE,
  vertices = string_nodes[, c("gene", "log2FoldChange", "padj", "Regulation")]
)

V(ppi_graph)$degree <- degree(ppi_graph)

hub_genes <- data.frame(
  Gene = V(ppi_graph)$name,
  Degree = V(ppi_graph)$degree
)

hub_genes <- hub_genes[order(-hub_genes$Degree), ]

print(hub_genes)
# STEP 57 — Save hub gene analysis

write.csv(
  hub_genes,
  "Results/Network/STRING/PPI_Hub_Genes.csv",
  row.names = FALSE
)

cat(
  "Hub genes file saved:",
  file.exists("Results/Network/STRING/PPI_Hub_Genes.csv"),
  "\n"
)

cat("Top hub gene:", hub_genes$Gene[1], "\n")
cat("Top degree:", hub_genes$Degree[1], "\n")

# STEP 58 — Prepare STRING edges for Cytoscape

cyto_string_edges <- ppi[, c(
  "from_name",
  "to_name",
  "combined_score"
)]

colnames(cyto_string_edges) <- c(
  "source",
  "target",
  "combined_score"
)

write.csv(
  cyto_string_edges,
  "Results/Network/STRING/Core_Genes_STRING_Cytoscape_Edges.csv",
  row.names = FALSE
)

cat(
  "Cytoscape STRING edges saved:",
  file.exists(
    "Results/Network/STRING/Core_Genes_STRING_Cytoscape_Edges.csv"
  ),
  "\n"
)

cat("Edges:", nrow(cyto_string_edges), "\n")
print(head(cyto_string_edges))

# STEP 59 — Prepare STRING nodes for Cytoscape

cyto_string_nodes <- merge(
  hub_genes,
  core_genes[, c("SYMBOL", "log2FoldChange", "padj", "Regulation")],
  by.x = "Gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

cyto_string_nodes$NodeType <- "Gene"

colnames(cyto_string_nodes)[1] <- "shared.name"

cyto_string_nodes <- cyto_string_nodes[, c(
  "shared.name",
  "Degree",
  "log2FoldChange",
  "padj",
  "Regulation",
  "NodeType"
)]

write.csv(
  cyto_string_nodes,
  "Results/Network/STRING/Core_Genes_STRING_Cytoscape_Nodes.csv",
  row.names = FALSE
)

cat(
  "Cytoscape STRING nodes saved:",
  file.exists(
    "Results/Network/STRING/Core_Genes_STRING_Cytoscape_Nodes.csv"
  ),
  "\n"
)

cat("Nodes:", nrow(cyto_string_nodes), "\n")
print(cyto_string_nodes)

# STEP 60 — Final STRING network validation

edge_genes <- unique(c(
  cyto_string_edges$source,
  cyto_string_edges$target
))

node_genes <- cyto_string_nodes$shared.name

cat("Total nodes:", length(node_genes), "\n")
cat("Total edges:", nrow(cyto_string_edges), "\n")
cat("Unique genes in edges:", length(edge_genes), "\n")

cat(
  "All edge genes exist in nodes:",
  all(edge_genes %in% node_genes),
  "\n"
)

cat(
  "Duplicate edges:",
  sum(duplicated(
    apply(
      cyto_string_edges[, c("source", "target")],
      1,
      function(x) paste(sort(x), collapse = "_")
    )
  )),
  "\n"
)
# STEP 61 — Extended PPI network centrality

V(ppi_graph)$degree <- degree(ppi_graph)
V(ppi_graph)$betweenness <- betweenness(
  ppi_graph,
  directed = FALSE,
  normalized = TRUE
)
V(ppi_graph)$closeness <- closeness(
  ppi_graph,
  normalized = TRUE
)

ppi_centrality <- data.frame(
  Gene = V(ppi_graph)$name,
  Degree = V(ppi_graph)$degree,
  Betweenness = V(ppi_graph)$betweenness,
  Closeness = V(ppi_graph)$closeness
)

ppi_centrality <- ppi_centrality[
  order(-ppi_centrality$Degree,
        -ppi_centrality$Betweenness),
]

print(ppi_centrality)

# STEP 62 — Save PPI centrality analysis

ppi_centrality <- merge(
  ppi_centrality,
  core_genes[, c("SYMBOL", "log2FoldChange", "padj", "Regulation")],
  by.x = "Gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

write.csv(
  ppi_centrality,
  "Results/Network/STRING/PPI_Centrality_Analysis.csv",
  row.names = FALSE
)

cat(
  "Centrality file saved:",
  file.exists(
    "Results/Network/STRING/PPI_Centrality_Analysis.csv"
  ),
  "\n"
)

cat("Genes analyzed:", nrow(ppi_centrality), "\n")

# STEP 63 — STRING PPI Network Figure

set.seed(123)

# Add edge weights
E(ppi_graph)$weight <- E(ppi_graph)$combined_score / 1000

# Node size based on degree
node_size <- 8 + (V(ppi_graph)$degree * 2)

# Node labels
node_labels <- V(ppi_graph)$name

png(
  "Results/Network/STRING/Core_Genes_STRING_PPI_Network.png",
  width = 2400,
  height = 2000,
  res = 300
)

plot(
  ppi_graph,
  layout = layout_with_fr(ppi_graph),
  vertex.size = node_size,
  vertex.label = node_labels,
  vertex.label.cex = 0.9,
  vertex.label.color = "black",
  vertex.frame.color = "black",
  edge.width = 1 + (E(ppi_graph)$weight * 3),
  edge.color = "grey70",
  main = "STRING PPI Network of Core Genes"
)

dev.off()

cat(
  "PPI network figure saved:",
  file.exists(
    "Results/Network/STRING/Core_Genes_STRING_PPI_Network.png"
  ),
  "\n"
)
# STEP 64 — Link PPI hubs to GSEA core pathways

leading_edge <- read.csv(
  "Results/Enrichment/GSEA/Leading_Edge_Genes_Reactome_Top3.csv",
  stringsAsFactors = FALSE
)

# Keep only the 12 significant core genes
hub_core <- ppi_centrality[
  ppi_centrality$Gene %in% core_genes$SYMBOL,
]

# Pathway membership
pathway_cols <- c(
  "G2M_Checkpoints",
  "Cell_Cycle_Checkpoints",
  "DNA_Replication"
)

hub_core_pathways <- leading_edge[
  leading_edge$SYMBOL %in% hub_core$Gene,
  c("SYMBOL", pathway_cols)
]

# Merge centrality with pathway membership
hub_core_pathways <- merge(
  hub_core,
  hub_core_pathways,
  by.x = "Gene",
  by.y = "SYMBOL",
  all.x = TRUE
)

# Save integrated table
write.csv(
  hub_core_pathways,
  "Results/Network/STRING/Hub_Genes_GSEA_Integration.csv",
  row.names = FALSE
)

cat(
  "Integrated file saved:",
  file.exists(
    "Results/Network/STRING/Hub_Genes_GSEA_Integration.csv"
  ),
  "\n"
)

cat("Genes integrated:", nrow(hub_core_pathways), "\n")

print(hub_core_pathways)

# STEP 65 — Final integrated core signature

final_core_signature <- hub_core_pathways

final_core_signature$N_Core_Pathways <- rowSums(
  final_core_signature[, c(
    "G2M_Checkpoints",
    "Cell_Cycle_Checkpoints",
    "DNA_Replication"
  )]
)

final_core_signature <- final_core_signature[
  order(
    -final_core_signature$Degree,
    -final_core_signature$log2FoldChange
  ),
]

write.csv(
  final_core_signature,
  "Results/Network/STRING/Final_Integrated_Core_Signature.csv",
  row.names = FALSE
)

cat(
  "Final signature saved:",
  file.exists(
    "Results/Network/STRING/Final_Integrated_Core_Signature.csv"
  ),
  "\n"
)

cat(
  "Genes:",
  nrow(final_core_signature),
  "\n"
)

cat(
  "Genes shared across all 3 pathways:",
  sum(final_core_signature$N_Core_Pathways == 3),
  "\n"
)

print(
  final_core_signature[, c(
    "Gene",
    "Degree",
    "Betweenness",
    "log2FoldChange",
    "padj",
    "N_Core_Pathways"
  )]
)
# STEP 66 — Final Integrated Core Signature Figure

library(ggplot2)

p_final_signature <- ggplot(
  final_core_signature,
  aes(
    x = log2FoldChange,
    y = Degree,
    size = Betweenness,
    label = Gene
  )
) +
  geom_point() +
  geom_text(
    vjust = -0.8,
    size = 3.5
  ) +
  scale_size_continuous(
    range = c(4, 10)
  ) +
  labs(
    title = "Integrated Core Gene Signature",
    subtitle = "Expression change and STRING network centrality",
    x = "log2 Fold Change",
    y = "STRING Degree",
    size = "Betweenness"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10)
  )

ggsave(
  "Results/Network/STRING/Final_Integrated_Core_Signature.png",
  plot = p_final_signature,
  width = 9,
  height = 7,
  dpi = 300
)

cat(
  "Final integrated figure saved:",
  file.exists(
    "Results/Network/STRING/Final_Integrated_Core_Signature.png"
  ),
  "\n"
)
# STEP 67 — Final Master Core Gene Table

master_core_table <- final_core_signature[, c(
  "Gene",
  "log2FoldChange",
  "padj",
  "Degree",
  "Betweenness",
  "Closeness",
  "N_Core_Pathways"
)]

master_core_table <- master_core_table[
  order(-master_core_table$Degree,
        -master_core_table$log2FoldChange),
]

write.csv(
  master_core_table,
  "Results/Network/STRING/MASTER_Core_Gene_Table.csv",
  row.names = FALSE
)

cat(
  "Master table saved:",
  file.exists(
    "Results/Network/STRING/MASTER_Core_Gene_Table.csv"
  ),
  "\n"
)

cat("Genes:", nrow(master_core_table), "\n")

print(master_core_table)

# STEP 68 — Inspect current project metadata

cat("Metadata dimensions:", dim(metadata_all), "\n\n")

cat("Metadata columns:\n")
print(names(metadata_all))

cat("\nFirst 6 rows:\n")
print(head(metadata_all))
# STEP 69 — Extract metadata from current DESeq2 object

metadata_primary <- as.data.frame(colData(dds_primary))

cat("Metadata dimensions:", dim(metadata_primary), "\n\n")

cat("Metadata columns:\n")
print(names(metadata_primary))

cat("\nFirst 6 rows:\n")
print(head(metadata_primary))
# STEP 70 — Inspect clinical metadata values

cat("Patient category:\n")
print(table(metadata_primary$patient_category, useNA = "ifany"))

cat("\nCOVID status:\n")
print(table(metadata_primary$covid_status, useNA = "ifany"))

cat("\nTime point:\n")
print(table(metadata_primary$time_point, useNA = "ifany"))

cat("\nAcuity:\n")
print(table(metadata_primary$acuity, useNA = "ifany"))

cat("\nCell type:\n")
print(table(metadata_primary$cell_type, useNA = "ifany"))
# STEP 71 — Extract and validate patient IDs

metadata_primary$patient_id <- sub("_.*$", "", rownames(metadata_primary))

cat("Number of unique patients:", 
    length(unique(metadata_primary$patient_id)), "\n\n")

cat("Samples per patient — first 20 patients:\n")
print(head(table(metadata_primary$patient_id), 20))

cat("\nExample longitudinal patients:\n")
print(
  head(
    metadata_primary[
      metadata_primary$patient_id %in%
        names(sort(table(metadata_primary$patient_id), decreasing = TRUE))[1:10],
      c("patient_id", "patient_category", "time_point", "acuity")
    ],
    30
  )
)
# STEP 72 — Longitudinal coverage per patient

patient_time_table <- table(
  metadata_primary$patient_id,
  metadata_primary$time_point
)

cat("Patient × Time-point matrix dimensions:\n")
print(dim(patient_time_table))

cat("\nNumber of patients with each time point:\n")
print(colSums(patient_time_table > 0))

cat("\nNumber of patients with each longitudinal pattern:\n")

patient_patterns <- apply(
  patient_time_table,
  1,
  function(x) paste(names(x)[x > 0], collapse = " + ")
)

print(sort(table(patient_patterns), decreasing = TRUE))

# STEP 73 — Check whether acuity is constant within each patient

acuity_by_patient <- tapply(
  metadata_primary$acuity,
  metadata_primary$patient_id,
  function(x) length(unique(x))
)

cat("Patients with constant acuity:",
    sum(acuity_by_patient == 1), "\n")

cat("Patients with changing acuity:",
    sum(acuity_by_patient > 1), "\n")

cat("\nPatients with changing acuity:\n")

changing_acuity_patients <- names(acuity_by_patient[acuity_by_patient > 1])

print(
  metadata_primary[
    metadata_primary$patient_id %in% changing_acuity_patients,
    c("patient_id", "patient_category", "time_point", "acuity")
  ]
)
# STEP 74 — Define complete longitudinal cohort

complete_longitudinal_patients <- names(
  which(
    rowSums(patient_time_table[, c("D0", "D3", "D7")] > 0) == 3
  )
)

metadata_longitudinal <- metadata_primary[
  metadata_primary$patient_id %in% complete_longitudinal_patients,
]

cat("Complete longitudinal patients:",
    length(complete_longitudinal_patients), "\n")

cat("Samples in longitudinal cohort:",
    nrow(metadata_longitudinal), "\n")

cat("\nTime-point distribution:\n")
print(table(metadata_longitudinal$time_point))

cat("\nPatient category distribution:\n")
print(table(
  metadata_longitudinal$patient_id,
  metadata_longitudinal$patient_category
))
# STEP 75 — Check duplicate samples within patient × time point

patient_time_counts <- as.data.frame(
  table(
    metadata_longitudinal$patient_id,
    metadata_longitudinal$time_point
  )
)

names(patient_time_counts) <- c(
  "patient_id",
  "time_point",
  "n_samples"
)

duplicates <- patient_time_counts[
  patient_time_counts$n_samples > 1,
]

cat("Patient × time-point combinations with >1 sample:",
    nrow(duplicates), "\n\n")

cat("Duplicate combinations:\n")
print(duplicates)

cat("\nCorresponding metadata:\n")

if (nrow(duplicates) > 0) {
  
  duplicate_ids <- paste(
    duplicates$patient_id,
    duplicates$time_point,
    sep = "_"
  )
  
  print(
    metadata_longitudinal[
      paste(
        metadata_longitudinal$patient_id,
        metadata_longitudinal$time_point,
        sep = "_"
      ) %in% duplicate_ids,
      c(
        "patient_id",
        "sample_id",
        "geo_accession",
        "time_point",
        "acuity",
        "rna_sample_id"
      )
    ]
  )
}

# STEP 76 — Compare duplicate D7 samples

duplicate_d7 <- c("154_D7A", "154_D7B", "43_D7A", "43_D7B")

duplicate_d7_counts <- counts(dds_primary, normalized = TRUE)[
  ,
  intersect(duplicate_d7, colnames(dds_primary))
]

cat("Duplicate D7 samples found:\n")
print(colnames(duplicate_d7_counts))

cat("\nSample correlations:\n")
print(round(cor(duplicate_d7_counts, method = "spearman"), 3))

# STEP 77 — Inspect metadata for duplicate D7 samples

duplicate_samples <- c(
  "154_D7A",
  "154_D7B",
  "43_D7A",
  "43_D7B"
)

print(
  metadata_primary[
    rownames(metadata_primary) %in% duplicate_samples,
    ,
    drop = FALSE
  ]
)
# ============================================================
# STEP 73 — Validate completed analysis before moving forward
# ============================================================

required_files <- c(
  "Results/Differential_Expression/COVID_vs_Symptomatic_all_genes.rds",
  "Results/Differential_Expression/COVID_vs_Symptomatic_significant_DEGs.csv",
  "Results/Enrichment/GSEA/GO_BP_GSEA.rds",
  "Results/Enrichment/GSEA/Reactome_GSEA.rds",
  "Results/Enrichment/GSEA/Core_Significant_Genes.csv",
  "Results/Network/STRING/Core_Genes_STRING_PPI.csv",
  "Results/Network/STRING/PPI_Centrality_Analysis.csv",
  "Results/Network/STRING/Final_Integrated_Core_Signature.csv",
  "Results/Network/STRING/MASTER_Core_Gene_Table.csv"
)

validation <- data.frame(
  File = required_files,
  Exists = file.exists(required_files)
)

print(validation)

cat("\nAll required files available:",
    all(validation$Exists), "\n")

# ============================================================
# STEP 74 — Inspect final integrated core signature
# ============================================================

core_signature <- read.csv(
  "Results/Network/STRING/Final_Integrated_Core_Signature.csv",
  stringsAsFactors = FALSE
)

cat("Number of core genes:", nrow(core_signature), "\n\n")

print(core_signature)

cat("\nRegulation:\n")
print(table(core_signature$Regulation))

cat("\nShared Reactome pathways:\n")
print(table(core_signature$N_Core_Pathways))

cat("\nDegree summary:\n")
print(summary(core_signature$Degree))

# ============================================================
# STEP 75 — Rank Core Genes
# ============================================================

core_ranked <- core_signature

# Rank by network degree first, then statistical significance
core_ranked <- core_ranked[
  order(
    -core_ranked$Degree,
    core_ranked$padj,
    -core_ranked$log2FoldChange
  ),
]

# Add rank
core_ranked$Rank <- seq_len(nrow(core_ranked))

# Reorder columns
core_ranked <- core_ranked[, c(
  "Rank",
  "Gene",
  "Degree",
  "Betweenness",
  "Closeness",
  "log2FoldChange",
  "padj",
  "Regulation",
  "N_Core_Pathways"
)]

print(core_ranked)

# Save ranked table
write.csv(
  core_ranked,
  "Results/Network/STRING/Ranked_Core_Genes.csv",
  row.names = FALSE
)

saveRDS(
  core_ranked,
  "Results/Network/STRING/Ranked_Core_Genes.rds"
)

cat("\nRanked core genes saved successfully.\n")

# ============================================================
# STEP 76 — Core Gene log2FC Visualization
# ============================================================

library(ggplot2)

core_plot <- ggplot(
  core_ranked,
  aes(
    x = reorder(Gene, log2FoldChange),
    y = log2FoldChange
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Core Gene Signature",
    subtitle = "COVID+ vs Symptomatic Controls",
    x = "Core Gene",
    y = "log2 Fold Change"
  ) +
  theme_minimal(base_size = 12)

print(core_plot)

ggsave(
  "Results/Network/STRING/Core_Gene_log2FC_Barplot.png",
  plot = core_plot,
  width = 8,
  height = 6,
  dpi = 300
)

cat(
  "\nPlot saved:",
  file.exists(
    "Results/Network/STRING/Core_Gene_log2FC_Barplot.png"
  ),
  "\n"
)

# ============================================================
# STEP 77 — Final Core Signature Annotation
# ============================================================

core_annotated <- core_ranked

core_annotated$Functional_Group <- NA_character_

core_annotated$Functional_Group[
  core_annotated$Gene %in% c("CCNA1", "CCNA2")
] <- "Cell cycle regulation"

core_annotated$Functional_Group[
  core_annotated$Gene %in% c("CDC6", "CDC45", "MCM2", "MCM4", "MCM10", "ORC1")
] <- "DNA replication"

core_annotated$Functional_Group[
  core_annotated$Gene %in% c("H2BC5", "H2BC7", "H2BC9", "H2BC17")
] <- "Chromatin / Histone"

print(core_annotated[, c(
  "Rank",
  "Gene",
  "Functional_Group",
  "Degree",
  "log2FoldChange",
  "padj",
  "N_Core_Pathways"
)])

# Save final annotated table
write.csv(
  core_annotated,
  "Results/Network/STRING/Final_Core_Signature_Annotated.csv",
  row.names = FALSE
)

saveRDS(
  core_annotated,
  "Results/Network/STRING/Final_Core_Signature_Annotated.rds"
)

cat(
  "\nAnnotated core signature saved:",
  file.exists(
