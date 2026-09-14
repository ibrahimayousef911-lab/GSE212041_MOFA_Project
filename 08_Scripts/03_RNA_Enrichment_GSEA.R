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

# STEP 19 - Functional Enrichment Analysis
# ============================================================

# ------------------------------------------------------------
# STEP 19A - Gene ID Conversion
# ------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Clean Ensembl IDs
# ------------------------------------------------------------

all_gene_ids <- sig_degs$gene_id
up_gene_ids <- upregulated_genes$gene_id
down_gene_ids <- downregulated_genes$gene_id

# Remove Ensembl version numbers if present
all_gene_ids <- sub("\\..*$", "", all_gene_ids)
up_gene_ids <- sub("\\..*$", "", up_gene_ids)
down_gene_ids <- sub("\\..*$", "", down_gene_ids)

# Remove duplicates
all_gene_ids <- unique(all_gene_ids)
up_gene_ids <- unique(up_gene_ids)
down_gene_ids <- unique(down_gene_ids)

# ------------------------------------------------------------
# Convert Ensembl IDs → Entrez IDs
# ------------------------------------------------------------

all_id_map <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_id_map <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_id_map <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# ------------------------------------------------------------
# Check conversion
# ------------------------------------------------------------

cat("Total significant DEGs:", length(all_gene_ids), "\n")
cat("DEGs mapped to Entrez:", nrow(all_id_map), "\n")

cat("Upregulated genes:", length(up_gene_ids), "\n")
cat("Upregulated mapped:", nrow(up_id_map), "\n")

cat("Downregulated genes:", length(down_gene_ids), "\n")
cat("Downregulated mapped:", nrow(down_id_map), "\n")

# Show first mappings
head(all_id_map)

# ============================================================
# STEP 19A - Check Enrichment Packages
# ============================================================

library(BiocManager)

BiocManager::valid()

requireNamespace("GO.db", quietly = TRUE)
requireNamespace("org.Hs.eg.db", quietly = TRUE)
requireNamespace("clusterProfiler", quietly = TRUE)
requireNamespace("DOSE", quietly = TRUE)
requireNamespace("enrichplot", quietly = TRUE)
# ============================================================
# STEP 19B - Install Enrichment Analysis Packages
# ============================================================

BiocManager::install(
  "clusterProfiler",
  ask = FALSE,
  update = FALSE
)

BiocManager::install(
  "DOSE",
  ask = FALSE,
  update = FALSE
)

BiocManager::install(
  "enrichplot",
  ask = FALSE,
  update = FALSE
)
# ============================================================
# STEP 19C - Ensembl to Entrez ID Mapping
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# Get Ensembl IDs
all_gene_ids <- unique(sig_degs$gene_id)
up_gene_ids <- unique(upregulated_genes$gene_id)
down_gene_ids <- unique(downregulated_genes$gene_id)

# Remove Ensembl version numbers if present
all_gene_ids <- sub("\\..*$", "", all_gene_ids)
up_gene_ids <- sub("\\..*$", "", up_gene_ids)
down_gene_ids <- sub("\\..*$", "", down_gene_ids)

# Convert Ensembl → Entrez
all_id_map <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_id_map <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_id_map <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Create unique Entrez gene lists
all_entrez <- unique(all_id_map$ENTREZID)
up_entrez <- unique(up_id_map$ENTREZID)
down_entrez <- unique(down_id_map$ENTREZID)

# Check mapping
cat("Significant DEGs:", length(all_gene_ids), "\n")
cat("Mapped Entrez IDs:", length(all_entrez), "\n\n")

cat("Upregulated genes:", length(up_gene_ids), "\n")
cat("Mapped Up Entrez IDs:", length(up_entrez), "\n\n")

cat("Downregulated genes:", length(down_gene_ids), "\n")
cat("Mapped Down Entrez IDs:", length(down_entrez), "\n")

# Preview mappings
head(all_id_map)
# ============================================================
# STEP 19D - GO Biological Process Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Create background gene list
# ------------------------------------------------------------

background_ensembl <- unique(res_valid$gene_id)

background_ensembl <- sub(
  "\\..*$",
  "",
  background_ensembl
)

background_map <- bitr(
  background_ensembl,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_entrez <- unique(
  background_map$ENTREZID
)

cat(
  "Background genes tested:",
  length(background_ensembl),
  "\n"
)

cat(
  "Background Entrez IDs:",
  length(background_entrez),
  "\n"
)

# ------------------------------------------------------------
# GO-BP enrichment - All significant DEGs
# ------------------------------------------------------------

ego_all <- enrichGO(
  gene = all_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP enrichment - Upregulated genes
# ------------------------------------------------------------

ego_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP enrichment - Downregulated genes
# ------------------------------------------------------------

ego_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# Check results
# ------------------------------------------------------------

cat(
  "GO-BP terms - All DEGs:",
  nrow(as.data.frame(ego_all)),
  "\n"
)

cat(
  "GO-BP terms - Upregulated:",
  nrow(as.data.frame(ego_up)),
  "\n"
)

cat(
  "GO-BP terms - Downregulated:",
  nrow(as.data.frame(ego_down)),
  "\n"
)

# ============================================================
# STEP 19D - Prepare GO Background
# ============================================================

background_ensembl <- unique(
  res_covid_vs_symptomatic$gene_id
)

background_ensembl <- sub(
  "\\..*$",
  "",
  background_ensembl
)

background_map <- bitr(
  background_ensembl,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_entrez <- unique(
  background_map$ENTREZID
)

cat(
  "Background Ensembl genes:",
  length(background_ensembl),
  "\n"
)

cat(
  "Background Entrez genes:",
  length(background_entrez),
  "\n"
)
# ============================================================
# STEP 19E - GO Biological Process Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# GO-BP: All Significant DEGs
# ------------------------------------------------------------

ego_all <- enrichGO(
  gene = all_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP: Upregulated Genes
# ------------------------------------------------------------

ego_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP: Downregulated Genes
# ------------------------------------------------------------

ego_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# Check Results
# ------------------------------------------------------------

cat(
  "GO-BP terms - All DEGs:",
  nrow(as.data.frame(ego_all)),
  "\n"
)

cat(
  "GO-BP terms - Upregulated:",
  nrow(as.data.frame(ego_up)),
  "\n"
)

cat(
  "GO-BP terms - Downregulated:",
  nrow(as.data.frame(ego_down)),
  "\n"
)
# ============================================================
# STEP 19E - Recreate Entrez Gene Lists
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# Significant DEG lists
all_gene_ids <- unique(sig_degs$gene_id)
up_gene_ids <- unique(upregulated_genes$gene_id)
down_gene_ids <- unique(downregulated_genes$gene_id)

# Remove Ensembl version numbers
all_gene_ids <- sub("\\..*$", "", all_gene_ids)
up_gene_ids <- sub("\\..*$", "", up_gene_ids)
down_gene_ids <- sub("\\..*$", "", down_gene_ids)

# Ensembl -> Entrez
all_id_map <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_id_map <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_id_map <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Unique Entrez IDs
all_entrez <- unique(all_id_map$ENTREZID)
up_entrez <- unique(up_id_map$ENTREZID)
down_entrez <- unique(down_id_map$ENTREZID)

# Check
cat("All DEGs:", length(all_gene_ids), "\n")
cat("Mapped all:", length(all_entrez), "\n")

cat("Upregulated:", length(up_gene_ids), "\n")
cat("Mapped up:", length(up_entrez), "\n")

cat("Downregulated:", length(down_gene_ids), "\n")
cat("Mapped down:", length(down_entrez), "\n")

# ============================================================
# STEP 19D - Recreate Entrez IDs for Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# Ensembl IDs
all_gene_ids <- unique(sig_degs$gene_id)
up_gene_ids <- unique(upregulated_genes$gene_id)
down_gene_ids <- unique(downregulated_genes$gene_id)

# Remove Ensembl version numbers
all_gene_ids <- sub("\\..*$", "", all_gene_ids)
up_gene_ids <- sub("\\..*$", "", up_gene_ids)
down_gene_ids <- sub("\\..*$", "", down_gene_ids)

# Ensembl -> Entrez
all_id_map <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_id_map <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_id_map <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Unique Entrez IDs
all_entrez <- unique(all_id_map$ENTREZID)
up_entrez <- unique(up_id_map$ENTREZID)
down_entrez <- unique(down_id_map$ENTREZID)

# Check
cat("All DEGs:", length(all_gene_ids), "\n")
cat("Mapped all:", length(all_entrez), "\n")
cat("Upregulated:", length(up_gene_ids), "\n")
cat("Mapped up:", length(up_entrez), "\n")
cat("Downregulated:", length(down_gene_ids), "\n")
cat("Mapped down:", length(down_entrez), "\n")
# ============================================================
# STEP 19E - GO Biological Process Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# GO-BP: All Significant DEGs
# ------------------------------------------------------------

ego_all <- enrichGO(
  gene = all_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP: Upregulated Genes
# ------------------------------------------------------------

ego_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# GO-BP: Downregulated Genes
# ------------------------------------------------------------

ego_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# ------------------------------------------------------------
# Check Results
# ------------------------------------------------------------

cat(
  "GO-BP terms - All DEGs:",
  nrow(as.data.frame(ego_all)),
  "\n"
)

cat(
  "GO-BP terms - Upregulated:",
  nrow(as.data.frame(ego_up)),
  "\n"
)

cat(
  "GO-BP terms - Downregulated:",
  nrow(as.data.frame(ego_down)),
  "\n"
)
# ============================================================
# STEP 19F - Inspect and Save GO-BP Results
# ============================================================

# Create GO results directory
if (!dir.exists("Results/Enrichment/GO")) {
  dir.create("Results/Enrichment/GO", recursive = TRUE)
}

# Convert enrichment results to data frames
go_all_df <- as.data.frame(ego_all)
go_up_df <- as.data.frame(ego_up)
go_down_df <- as.data.frame(ego_down)

# Check number of enriched terms
cat("GO-BP - All DEGs:", nrow(go_all_df), "\n")
cat("GO-BP - Upregulated:", nrow(go_up_df), "\n")
cat("GO-BP - Downregulated:", nrow(go_down_df), "\n")

# Show top 10 terms
cat("\n===== TOP 10 GO-BP: UPREGULATED =====\n")
print(
  go_up_df[
    order(go_up_df$p.adjust),
    c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "Count")
  ][1:min(10, nrow(go_up_df)), ]
)

cat("\n===== TOP 10 GO-BP: DOWNREGULATED =====\n")
print(
  go_down_df[
    order(go_down_df$p.adjust),
    c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "Count")
  ][1:min(10, nrow(go_down_df)), ]
)

# Save tables
write.csv(
  go_all_df,
  "Results/Enrichment/GO/GO_BP_all_DEGs.csv",
  row.names = FALSE
)

write.csv(
  go_up_df,
  "Results/Enrichment/GO/GO_BP_upregulated.csv",
  row.names = FALSE
)

write.csv(
  go_down_df,
  "Results/Enrichment/GO/GO_BP_downregulated.csv",
  row.names = FALSE
)

# Save RDS objects
saveRDS(
  ego_all,
  "Results/Enrichment/GO/GO_BP_all_DEGs.rds"
)

saveRDS(
  ego_up,
  "Results/Enrichment/GO/GO_BP_upregulated.rds"
)

saveRDS(
  ego_down,
  "Results/Enrichment/GO/GO_BP_downregulated.rds"
)

# Verify files
cat("\n===== FILE CHECK =====\n")

cat(
  "All GO:",
  file.exists("Results/Enrichment/GO/GO_BP_all_DEGs.csv"),
  "\n"
)

cat(
  "Up GO:",
  file.exists("Results/Enrichment/GO/GO_BP_upregulated.csv"),
  "\n"
)

cat(
  "Down GO:",
  file.exists("Results/Enrichment/GO/GO_BP_downregulated.csv"),
  "\n"
)
# ============================================================
# STEP 19G - GO-BP Visualization
# ============================================================

library(enrichplot)
library(ggplot2)

# Create GO figures directory
if (!dir.exists("Results/Enrichment/GO/Figures")) {
  dir.create(
    "Results/Enrichment/GO/Figures",
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# UPREGULATED GO-BP
# ------------------------------------------------------------

p_go_up <- dotplot(
  ego_up,
  showCategory = 20
) +
  ggtitle("GO Biological Process - Upregulated Genes")

print(p_go_up)

ggsave(
  "Results/Enrichment/GO/Figures/GO_BP_Upregulated_Dotplot.png",
  p_go_up,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# DOWNREGULATED GO-BP
# ------------------------------------------------------------

p_go_down <- dotplot(
  ego_down,
  showCategory = 14
) +
  ggtitle("GO Biological Process - Downregulated Genes")

print(p_go_down)

ggsave(
  "Results/Enrichment/GO/Figures/GO_BP_Downregulated_Dotplot.png",
  p_go_down,
  width = 10,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# ALL DEGs GO-BP
# ------------------------------------------------------------

p_go_all <- dotplot(
  ego_all,
  showCategory = 20
) +
  ggtitle("GO Biological Process - All Significant DEGs")

print(p_go_all)

ggsave(
  "Results/Enrichment/GO/Figures/GO_BP_All_DEGs_Dotplot.png",
  p_go_all,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Verify figures
# ------------------------------------------------------------

cat(
  "Upregulated figure:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_BP_Upregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "Downregulated figure:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_BP_Downregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "All DEGs figure:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_BP_All_DEGs_Dotplot.png"
  ),
  "\n"
)

# ============================================================
# STEP 19 - KEGG Pathway Enrichment Analysis
# ============================================================

library(clusterProfiler)
library(ggplot2)

# ------------------------------------------------------------
# STEP 19A - Create KEGG results directory
# ------------------------------------------------------------

if (!dir.exists("Results/Enrichment/KEGG")) {
  dir.create(
    "Results/Enrichment/KEGG",
    recursive = TRUE
  )
}

if (!dir.exists("Results/Enrichment/KEGG/Figures")) {
  dir.create(
    "Results/Enrichment/KEGG/Figures",
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# STEP 19B - Check input gene lists
# ------------------------------------------------------------

cat("All significant Entrez genes:",
    length(all_entrez), "\n")

cat("Upregulated Entrez genes:",
    length(up_entrez), "\n")

cat("Downregulated Entrez genes:",
    length(down_entrez), "\n")

cat("Background Entrez genes:",
    length(background_entrez), "\n")


# ------------------------------------------------------------
# STEP 19C - KEGG enrichment: ALL DEGs
# ------------------------------------------------------------

ekegg_all <- enrichKEGG(
  gene = all_entrez,
  universe = background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500
)

# ------------------------------------------------------------
# STEP 19D - KEGG enrichment: UPREGULATED
# ------------------------------------------------------------

ekegg_up <- enrichKEGG(
  gene = up_entrez,
  universe = background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 500
)

# ------------------------------------------------------------
# STEP 19E - KEGG enrichment: DOWNREGULATED
# ------------------------------------------------------------

ekegg_down <- enrichKEGG(
  gene = down_entrez,
  universe = background_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 500
)

# ------------------------------------------------------------
# STEP 19F - Convert to data frames
# ------------------------------------------------------------

kegg_all_df <- as.data.frame(ekegg_all)
kegg_up_df <- as.data.frame(ekegg_up)
kegg_down_df <- as.data.frame(ekegg_down)

# ------------------------------------------------------------
# STEP 19G - Check number of pathways
# ------------------------------------------------------------

cat("\n===== KEGG RESULTS =====\n")

cat(
  "KEGG pathways - All DEGs:",
  nrow(kegg_all_df),
  "\n"
)

cat(
  "KEGG pathways - Upregulated:",
  nrow(kegg_up_df),
  "\n"
)

cat(
  "KEGG pathways - Downregulated:",
  nrow(kegg_down_df),
  "\n"
)

# ------------------------------------------------------------
# STEP 19H - Show top pathways
# ------------------------------------------------------------

if (nrow(kegg_all_df) > 0) {
  
  cat("\n===== TOP KEGG - ALL DEGs =====\n")
  
  print(
    kegg_all_df[
      order(kegg_all_df$p.adjust),
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ][1:min(15, nrow(kegg_all_df)), ]
  )
}


if (nrow(kegg_up_df) > 0) {
  
  cat("\n===== TOP KEGG - UPREGULATED =====\n")
  
  print(
    kegg_up_df[
      order(kegg_up_df$p.adjust),
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ][1:min(15, nrow(kegg_up_df)), ]
  )
}


if (nrow(kegg_down_df) > 0) {
  
  cat("\n===== TOP KEGG - DOWNREGULATED =====\n")
  
  print(
    kegg_down_df[
      order(kegg_down_df$p.adjust),
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ][1:min(15, nrow(kegg_down_df)), ]
  )
}

# ------------------------------------------------------------
# STEP 19I - Save KEGG tables
# ------------------------------------------------------------

write.csv(
  kegg_all_df,
  "Results/Enrichment/KEGG/KEGG_all_DEGs.csv",
  row.names = FALSE
)

write.csv(
  kegg_up_df,
  "Results/Enrichment/KEGG/KEGG_upregulated.csv",
  row.names = FALSE
)

write.csv(
  kegg_down_df,
  "Results/Enrichment/KEGG/KEGG_downregulated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# STEP 19J - Save RDS objects
# ------------------------------------------------------------

saveRDS(
  ekegg_all,
  "Results/Enrichment/KEGG/KEGG_all_DEGs.rds"
)

saveRDS(
  ekegg_up,
  "Results/Enrichment/KEGG/KEGG_upregulated.rds"
)

saveRDS(
  ekegg_down,
  "Results/Enrichment/KEGG/KEGG_downregulated.rds"
)

# ------------------------------------------------------------
# STEP 19K - Verify files
# ------------------------------------------------------------

cat("\n===== KEGG FILE CHECK =====\n")

cat(
  "All KEGG:",
  file.exists(
    "Results/Enrichment/KEGG/KEGG_all_DEGs.csv"
  ),
  "\n"
)

cat(
  "Up KEGG:",
  file.exists(
    "Results/Enrichment/KEGG/KEGG_upregulated.csv"
  ),
  "\n"
)

cat(
  "Down KEGG:",
  file.exists(
    "Results/Enrichment/KEGG/KEGG_downregulated.csv"
  ),
  "\n"
)
# ============================================================
# STEP 20 - KEGG Enrichment Figures
# ============================================================

library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# STEP 20A - Create Figures directory
# ------------------------------------------------------------

if (!dir.exists("Results/Enrichment/KEGG/Figures")) {
  dir.create(
    "Results/Enrichment/KEGG/Figures",
    recursive = TRUE
  )
}

# ------------------------------------------------------------
# STEP 20B - KEGG All DEGs
# ------------------------------------------------------------

p_kegg_all <- dotplot(
  ekegg_all,
  showCategory = min(20, nrow(as.data.frame(ekegg_all)))
) +
  ggtitle("KEGG Pathway Enrichment - All Significant DEGs")

print(p_kegg_all)

ggsave(
  "Results/Enrichment/KEGG/Figures/KEGG_All_DEGs_Dotplot.png",
  p_kegg_all,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# STEP 20C - KEGG Upregulated
# ------------------------------------------------------------

if (nrow(as.data.frame(ekegg_up)) > 0) {
  
  p_kegg_up <- dotplot(
    ekegg_up,
    showCategory = min(20, nrow(as.data.frame(ekegg_up)))
  ) +
    ggtitle("KEGG Pathway Enrichment - Upregulated Genes")
  
  print(p_kegg_up)
  
  ggsave(
    "Results/Enrichment/KEGG/Figures/KEGG_Upregulated_Dotplot.png",
    p_kegg_up,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ------------------------------------------------------------
# STEP 20D - KEGG Downregulated
# ------------------------------------------------------------

if (nrow(as.data.frame(ekegg_down)) > 0) {
  
  p_kegg_down <- dotplot(
    ekegg_down,
    showCategory = min(20, nrow(as.data.frame(ekegg_down)))
  ) +
    ggtitle("KEGG Pathway Enrichment - Downregulated Genes")
  
  print(p_kegg_down)
  
  ggsave(
    "Results/Enrichment/KEGG/Figures/KEGG_Downregulated_Dotplot.png",
    p_kegg_down,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ------------------------------------------------------------
# STEP 20E - Verify Figures
# ------------------------------------------------------------

cat("\n===== KEGG FIGURE CHECK =====\n")

cat(
  "All KEGG figure:",
  file.exists(
    "Results/Enrichment/KEGG/Figures/KEGG_All_DEGs_Dotplot.png"
  ),
  "\n"
)

cat(
  "Up KEGG figure:",
  file.exists(
    "Results/Enrichment/KEGG/Figures/KEGG_Upregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "Down KEGG figure:",
  file.exists(
    "Results/Enrichment/KEGG/Figures/KEGG_Downregulated_Dotplot.png"
  ),
  "\n"
)
# ============================================================
# STEP 21A - Check ReactomePA
# ============================================================

cat(
  "ReactomePA:",
  requireNamespace("ReactomePA", quietly = TRUE),
  "\n"
)
# ============================================================
# STEP 21A - Install ReactomePA
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("ReactomePA", ask = FALSE, update = FALSE)

# ============================================================
# STEP 22 - Restore Gene IDs and Background
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Define gene IDs
# ------------------------------------------------------------

all_gene_ids <- unique(sig_degs$gene_id)
up_gene_ids <- unique(upregulated_genes$gene_id)
down_gene_ids <- unique(downregulated_genes$gene_id)

# ------------------------------------------------------------
# Convert Ensembl -> Entrez
# ------------------------------------------------------------

all_id_map <- bitr(
  all_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

up_id_map <- bitr(
  up_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

down_id_map <- bitr(
  down_gene_ids,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# ------------------------------------------------------------
# Unique Entrez IDs
# ------------------------------------------------------------

all_entrez <- unique(all_id_map$ENTREZID)
up_entrez <- unique(up_id_map$ENTREZID)
down_entrez <- unique(down_id_map$ENTREZID)

# ------------------------------------------------------------
# Background
# ------------------------------------------------------------

background_ensembl <- unique(res_covid_vs_symptomatic$gene_id)

background_id_map <- bitr(
  background_ensembl,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_entrez <- unique(
  background_id_map$ENTREZID
)

# ------------------------------------------------------------
# Check
# ------------------------------------------------------------

cat("\n===== GENE ID CHECK =====\n")

cat(
  "Significant DEGs:",
  length(all_gene_ids),
  "\n"
)

cat(
  "Mapped all Entrez:",
  length(all_entrez),
  "\n"
)

cat(
  "Upregulated:",
  length(up_gene_ids),
  "\n"
)

cat(
  "Mapped up Entrez:",
  length(up_entrez),
  "\n"
)

cat(
  "Downregulated:",
  length(down_gene_ids),
  "\n"
)

cat(
  "Mapped down Entrez:",
  length(down_entrez),
  "\n"
)

cat(
  "Background Ensembl:",
  length(background_ensembl),
  "\n"
)

cat(
  "Background Entrez:",
  length(background_entrez),
  "\n"
)
# ============================================================
# STEP 23 - Install ReactomePA
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  "ReactomePA",
  ask = FALSE,
  update = FALSE
)
# ============================================================
# STEP 23A - Install missing Reactome dependency
# ============================================================

BiocManager::install(
  "reactome.db",
  ask = FALSE,
  update = FALSE,
  force = TRUE
)
# ============================================================
# STEP 24 - Reactome Pathway Enrichment
# ============================================================

library(ReactomePA)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# 1. Reactome enrichment - All significant DEGs
# ------------------------------------------------------------

reactome_all <- enrichPathway(
  gene = all_entrez,
  universe = background_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 2. Reactome enrichment - Upregulated genes
# ------------------------------------------------------------

reactome_up <- enrichPathway(
  gene = up_entrez,
  universe = background_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 3. Reactome enrichment - Downregulated genes
# ------------------------------------------------------------

reactome_down <- enrichPathway(
  gene = down_entrez,
  universe = background_entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 4. Check results
# ------------------------------------------------------------

cat("\n===== REACTOME RESULTS =====\n")
cat("Reactome All pathways:", nrow(as.data.frame(reactome_all)), "\n")
cat("Reactome Up pathways:", nrow(as.data.frame(reactome_up)), "\n")
cat("Reactome Down pathways:", nrow(as.data.frame(reactome_down)), "\n")

# ============================================================
# STEP 25 - Save Reactome Results
# ============================================================

# Create output directories
dir.create(
  "Results/Enrichment/Reactome",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/Enrichment/Reactome/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Save tables
# ------------------------------------------------------------

write.csv(
  as.data.frame(reactome_all),
  "Results/Enrichment/Reactome/Reactome_all_DEGs.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(reactome_up),
  "Results/Enrichment/Reactome/Reactome_upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(reactome_down),
  "Results/Enrichment/Reactome/Reactome_downregulated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Save RDS objects
# ------------------------------------------------------------

saveRDS(
  reactome_all,
  "Results/Enrichment/Reactome/Reactome_all_DEGs.rds"
)

saveRDS(
  reactome_up,
  "Results/Enrichment/Reactome/Reactome_upregulated.rds"
)

saveRDS(
  reactome_down,
  "Results/Enrichment/Reactome/Reactome_downregulated.rds"
)

# ------------------------------------------------------------
# Verify files
# ------------------------------------------------------------

cat("\n===== REACTOME FILE CHECK =====\n")

cat(
  "All CSV:",
  file.exists("Results/Enrichment/Reactome/Reactome_all_DEGs.csv"),
  "\n"
)

cat(
  "Up CSV:",
  file.exists("Results/Enrichment/Reactome/Reactome_upregulated.csv"),
  "\n"
)

cat(
  "Down CSV:",
  file.exists("Results/Enrichment/Reactome/Reactome_downregulated.csv"),
  "\n"
)

cat(
  "All RDS:",
  file.exists("Results/Enrichment/Reactome/Reactome_all_DEGs.rds"),
  "\n"
)

cat(
  "Up RDS:",
  file.exists("Results/Enrichment/Reactome/Reactome_upregulated.rds"),
  "\n"
)

cat(
  "Down RDS:",
  file.exists("Results/Enrichment/Reactome/Reactome_downregulated.rds"),
  "\n"
)

# ============================================================
# STEP 26 - Reactome Dotplots
# ============================================================

library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# 1. All significant DEGs
# ------------------------------------------------------------

reactome_all_plot <- dotplot(
  reactome_all,
  showCategory = 20
) +
  ggtitle("Reactome Pathway Enrichment - All Significant DEGs")

ggsave(
  "Results/Enrichment/Reactome/Figures/Reactome_All_DEGs_Dotplot.png",
  reactome_all_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 2. Upregulated genes
# ------------------------------------------------------------

reactome_up_plot <- dotplot(
  reactome_up,
  showCategory = 20
) +
  ggtitle("Reactome Pathway Enrichment - Upregulated Genes")

ggsave(
  "Results/Enrichment/Reactome/Figures/Reactome_Upregulated_Dotplot.png",
  reactome_up_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 3. Downregulated genes
# ------------------------------------------------------------

if (nrow(as.data.frame(reactome_down)) > 0) {
  
  reactome_down_plot <- dotplot(
    reactome_down,
    showCategory = 20
  ) +
    ggtitle("Reactome Pathway Enrichment - Downregulated Genes")
  
  ggsave(
    "Results/Enrichment/Reactome/Figures/Reactome_Downregulated_Dotplot.png",
    reactome_down_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  
  cat(
    "\nNo significant Reactome pathways for downregulated genes.\n"
  )
}

# ------------------------------------------------------------
# 4. Verify figures
# ------------------------------------------------------------

cat("\n===== REACTOME FIGURE CHECK =====\n")

cat(
  "All:",
  file.exists(
    "Results/Enrichment/Reactome/Figures/Reactome_All_DEGs_Dotplot.png"
  ),
  "\n"
)

cat(
  "Up:",
  file.exists(
    "Results/Enrichment/Reactome/Figures/Reactome_Upregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "Down:",
  file.exists(
    "Results/Enrichment/Reactome/Figures/Reactome_Downregulated_Dotplot.png"
  ),
  "\n"
)
# ============================================================
# STEP 27 - Top Reactome Pathways
# ============================================================

# Top 20 Reactome pathways - All DEGs
top_reactome_all <- as.data.frame(reactome_all) |>
  dplyr::arrange(p.adjust) |>
  dplyr::select(
    ID,
    Description,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    geneID
  ) |>
  head(20)

# Top 20 Reactome pathways - Upregulated
top_reactome_up <- as.data.frame(reactome_up) |>
  dplyr::arrange(p.adjust) |>
  dplyr::select(
    ID,
    Description,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    qvalue,
    Count,
    geneID
  ) |>
  head(20)

# Print results
cat("\n===== TOP REACTOME - ALL DEGs =====\n")
print(top_reactome_all)

cat("\n===== TOP REACTOME - UPREGULATED =====\n")
print(top_reactome_up)
# ============================================================
# STEP 28 - Save Top Reactome Pathways
# ============================================================

write.csv(
  top_reactome_all,
  "Results/Enrichment/Reactome/Top20_Reactome_All_DEGs.csv",
  row.names = FALSE
)

write.csv(
  top_reactome_up,
  "Results/Enrichment/Reactome/Top20_Reactome_Upregulated.csv",
  row.names = FALSE
)

cat("\n===== TOP REACTOME TABLE CHECK =====\n")

cat(
  "All:",
  file.exists(
    "Results/Enrichment/Reactome/Top20_Reactome_All_DEGs.csv"
  ),
  "\n"
)

cat(
  "Up:",
  file.exists(
    "Results/Enrichment/Reactome/Top20_Reactome_Upregulated.csv"
  ),
  "\n"
)
# ============================================================
# STEP 29 - GO Molecular Function Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# 1. All significant DEGs
# ------------------------------------------------------------

ego_mf_all <- enrichGO(
  gene = all_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 2. Upregulated genes
# ------------------------------------------------------------

ego_mf_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 3. Downregulated genes
# ------------------------------------------------------------

ego_mf_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 4. Check results
# ------------------------------------------------------------

cat("\n===== GO MOLECULAR FUNCTION RESULTS =====\n")

cat(
  "GO-MF All pathways:",
  nrow(as.data.frame(ego_mf_all)),
  "\n"
)

cat(
  "GO-MF Up pathways:",
  nrow(as.data.frame(ego_mf_up)),
  "\n"
)

cat(
  "GO-MF Down pathways:",
  nrow(as.data.frame(ego_mf_down)),
  "\n"
)
# ============================================================
# STEP 30 - Save GO Molecular Function Results
# ============================================================

dir.create(
  "Results/Enrichment/GO/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Save CSV tables
# ------------------------------------------------------------

write.csv(
  as.data.frame(ego_mf_all),
  "Results/Enrichment/GO/GO_MF_all_DEGs.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_mf_up),
  "Results/Enrichment/GO/GO_MF_upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_mf_down),
  "Results/Enrichment/GO/GO_MF_downregulated.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Save RDS objects
# ------------------------------------------------------------

saveRDS(
  ego_mf_all,
  "Results/Enrichment/GO/GO_MF_all_DEGs.rds"
)

saveRDS(
  ego_mf_up,
  "Results/Enrichment/GO/GO_MF_upregulated.rds"
)

saveRDS(
  ego_mf_down,
  "Results/Enrichment/GO/GO_MF_downregulated.rds"
)

# ------------------------------------------------------------
# Verify
# ------------------------------------------------------------

cat("\n===== GO-MF FILE CHECK =====\n")

cat(
  "All CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_all_DEGs.csv"
  ),
  "\n"
)

cat(
  "Up CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_upregulated.csv"
  ),
  "\n"
)

cat(
  "Down CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_downregulated.csv"
  ),
  "\n"
)

cat(
  "All RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_all_DEGs.rds"
  ),
  "\n"
)

cat(
  "Up RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_upregulated.rds"
  ),
  "\n"
)

cat(
  "Down RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_MF_downregulated.rds"
  ),
  "\n"
)
# ============================================================
# STEP 31 - GO Molecular Function Dotplots
# ============================================================

library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# All significant DEGs
# ------------------------------------------------------------

go_mf_all_plot <- dotplot(
  ego_mf_all,
  showCategory = 20
) +
  ggtitle("GO Molecular Function - All Significant DEGs")

ggsave(
  "Results/Enrichment/GO/Figures/GO_MF_All_DEGs_Dotplot.png",
  go_mf_all_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Upregulated
# ------------------------------------------------------------

go_mf_up_plot <- dotplot(
  ego_mf_up,
  showCategory = 20
) +
  ggtitle("GO Molecular Function - Upregulated Genes")

ggsave(
  "Results/Enrichment/GO/Figures/GO_MF_Upregulated_Dotplot.png",
  go_mf_up_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Downregulated
# ------------------------------------------------------------

if (nrow(as.data.frame(ego_mf_down)) > 0) {
  
  go_mf_down_plot <- dotplot(
    ego_mf_down,
    showCategory = 20
  ) +
    ggtitle("GO Molecular Function - Downregulated Genes")
  
  ggsave(
    "Results/Enrichment/GO/Figures/GO_MF_Downregulated_Dotplot.png",
    go_mf_down_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  
  cat(
    "\nNo significant GO-MF terms for downregulated genes.\n"
  )
}

# ------------------------------------------------------------
# Verify figures
# ------------------------------------------------------------

cat("\n===== GO-MF FIGURE CHECK =====\n")

cat(
  "All:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_MF_All_DEGs_Dotplot.png"
  ),
  "\n"
)

cat(
  "Up:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_MF_Upregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "Down:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_MF_Downregulated_Dotplot.png"
  ),
  "\n"
)

# ============================================================
# STEP 32 - GO Cellular Component Enrichment
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# 1. All significant DEGs
# ------------------------------------------------------------

ego_cc_all <- enrichGO(
  gene = all_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 2. Upregulated genes
# ------------------------------------------------------------

ego_cc_up <- enrichGO(
  gene = up_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 3. Downregulated genes
# ------------------------------------------------------------

ego_cc_down <- enrichGO(
  gene = down_entrez,
  universe = background_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# ------------------------------------------------------------
# 4. Check results
# ------------------------------------------------------------

cat("\n===== GO CELLULAR COMPONENT RESULTS =====\n")

cat(
  "GO-CC All pathways:",
  nrow(as.data.frame(ego_cc_all)),
  "\n"
)

cat(
  "GO-CC Up pathways:",
  nrow(as.data.frame(ego_cc_up)),
  "\n"
)

cat(
  "GO-CC Down pathways:",
  nrow(as.data.frame(ego_cc_down)),
  "\n"
)
# ============================================================
# STEP 33 - Save GO Cellular Component Results
# ============================================================

# Save CSV tables

write.csv(
  as.data.frame(ego_cc_all),
  "Results/Enrichment/GO/GO_CC_all_DEGs.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_cc_up),
  "Results/Enrichment/GO/GO_CC_upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_cc_down),
  "Results/Enrichment/GO/GO_CC_downregulated.csv",
  row.names = FALSE
)

# Save RDS objects

saveRDS(
  ego_cc_all,
  "Results/Enrichment/GO/GO_CC_all_DEGs.rds"
)

saveRDS(
  ego_cc_up,
  "Results/Enrichment/GO/GO_CC_upregulated.rds"
)

saveRDS(
  ego_cc_down,
  "Results/Enrichment/GO/GO_CC_downregulated.rds"
)

# Verify files

cat("\n===== GO-CC FILE CHECK =====\n")

cat(
  "All CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_all_DEGs.csv"
  ),
  "\n"
)

cat(
  "Up CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_upregulated.csv"
  ),
  "\n"
)

cat(
  "Down CSV:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_downregulated.csv"
  ),
  "\n"
)

cat(
  "All RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_all_DEGs.rds"
  ),
  "\n"
)

cat(
  "Up RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_upregulated.rds"
  ),
  "\n"
)

cat(
  "Down RDS:",
  file.exists(
    "Results/Enrichment/GO/GO_CC_downregulated.rds"
  ),
  "\n"
)
# ============================================================
# STEP 34 - GO Cellular Component Dotplots
# ============================================================

library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# All significant DEGs
# ------------------------------------------------------------

go_cc_all_plot <- dotplot(
  ego_cc_all,
  showCategory = 20
) +
  ggtitle("GO Cellular Component - All Significant DEGs")

ggsave(
  "Results/Enrichment/GO/Figures/GO_CC_All_DEGs_Dotplot.png",
  go_cc_all_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Upregulated
# ------------------------------------------------------------

go_cc_up_plot <- dotplot(
  ego_cc_up,
  showCategory = 20
) +
  ggtitle("GO Cellular Component - Upregulated Genes")

ggsave(
  "Results/Enrichment/GO/Figures/GO_CC_Upregulated_Dotplot.png",
  go_cc_up_plot,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Downregulated
# ------------------------------------------------------------

if (nrow(as.data.frame(ego_cc_down)) > 0) {
  
  go_cc_down_plot <- dotplot(
    ego_cc_down,
    showCategory = 20
  ) +
    ggtitle("GO Cellular Component - Downregulated Genes")
  
  ggsave(
    "Results/Enrichment/GO/Figures/GO_CC_Downregulated_Dotplot.png",
    go_cc_down_plot,
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  
  cat(
    "\nNo significant GO-CC terms for downregulated genes.\n"
  )
}

# ------------------------------------------------------------
# Verify figures
# ------------------------------------------------------------

cat("\n===== GO-CC FIGURE CHECK =====\n")

cat(
  "All:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_CC_All_DEGs_Dotplot.png"
  ),
  "\n"
)

cat(
  "Up:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_CC_Upregulated_Dotplot.png"
  ),
  "\n"
)

cat(
  "Down:",
  file.exists(
    "Results/Enrichment/GO/Figures/GO_CC_Downregulated_Dotplot.png"
  ),
  "\n"
)

# ============================================================
# STEP 35 - Prepare Ranked Gene List for GSEA
# ============================================================

# Make sure DESeq2 results are available
head(res_covid_vs_symptomatic)

# ------------------------------------------------------------
# Create ranked statistics
# ------------------------------------------------------------

gsea_df <- as.data.frame(res_covid_vs_symptomatic)

# Keep genes with valid statistics
gsea_df <- gsea_df[
  !is.na(gsea_df$stat) &
    !is.na(gsea_df$gene_id),
]

# Remove duplicated Ensembl IDs
gsea_df <- gsea_df[
  !duplicated(gsea_df$gene_id),
]

# Create ranked vector
gene_ranking <- gsea_df$stat

names(gene_ranking) <- gsea_df$gene_id

# Sort decreasing
gene_ranking <- sort(
  gene_ranking,
  decreasing = TRUE
)

# ------------------------------------------------------------
# Check ranking
# ------------------------------------------------------------

cat("\n===== GSEA RANKING CHECK =====\n")

cat(
  "Total ranked genes:",
  length(gene_ranking),
  "\n"
)

cat(
  "Highest statistic:",
  max(gene_ranking),
  "\n"
)

cat(
  "Lowest statistic:",
  min(gene_ranking),
  "\n"
)

cat(
  "NA values:",
  sum(is.na(gene_ranking)),
  "\n"
)

cat("\nTop 10 genes:\n")
print(head(gene_ranking, 10))

cat("\nBottom 10 genes:\n")
print(tail(gene_ranking, 10))
# ============================================================
# STEP 36 - Convert GSEA Ranking: Ensembl -> Entrez
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ------------------------------------------------------------
# Convert Ensembl IDs to Entrez IDs
# ------------------------------------------------------------

gene_map_gsea <- bitr(
  names(gene_ranking),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# ------------------------------------------------------------
# Remove duplicated mappings
# Keep one Ensembl gene per Entrez ID
# ------------------------------------------------------------

gene_map_gsea <- gene_map_gsea[
  !duplicated(gene_map_gsea$ENTREZID),
]

# ------------------------------------------------------------
# Match ranking statistics to mapped genes
# ------------------------------------------------------------

gene_map_gsea$stat <- gene_ranking[
  match(
    gene_map_gsea$ENSEMBL,
    names(gene_ranking)
  )
]

# ------------------------------------------------------------
# Create Entrez-ranked vector
# ------------------------------------------------------------

gene_ranking_entrez <- gene_map_gsea$stat

names(gene_ranking_entrez) <- gene_map_gsea$ENTREZID

gene_ranking_entrez <- sort(
  gene_ranking_entrez,
  decreasing = TRUE
)

# ------------------------------------------------------------
# Check
# ------------------------------------------------------------

cat("\n===== ENTREZ GSEA RANKING CHECK =====\n")

cat(
  "Original Ensembl genes:",
  length(gene_ranking),
  "\n"
)

cat(
  "Mapped Entrez genes:",
  length(gene_ranking_entrez),
  "\n"
)

cat(
  "Mapping percentage:",
  round(
    100 * length(gene_ranking_entrez) /
      length(gene_ranking),
    2
  ),
  "%\n"
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

cat("\nTop 10 Entrez genes:\n")
print(head(gene_ranking_entrez, 10))

cat("\nBottom 10 Entrez genes:\n")
print(tail(gene_ranking_entrez, 10))

# ============================================================
# STEP 37 - GSEA: GO Biological Process
# ============================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# Run GO-BP GSEA
# ------------------------------------------------------------

gsea_go_bp <- gseGO(
  geneList = gene_ranking_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE
)

# ------------------------------------------------------------
# Check result
# ------------------------------------------------------------

cat("\n===== GO-BP GSEA RESULT =====\n")

cat(
  "Number of enriched pathways:",
  nrow(as.data.frame(gsea_go_bp)),
  "\n"
)

if (nrow(as.data.frame(gsea_go_bp)) > 0) {
  
  cat("\nTop 10 enriched GO-BP pathways:\n")
  
  print(
    as.data.frame(gsea_go_bp)[
      1:min(10, nrow(as.data.frame(gsea_go_bp))),
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
  
} else {
  
  cat("\nNo significant GO-BP pathways detected.\n")
}

# ============================================================
# STEP 38 - Save GO-BP GSEA Results + Figures
# ============================================================

# Create output directories
dir.create(
  "Results/Enrichment/GSEA",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "Results/Enrichment/GSEA/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Convert result to data frame
# ------------------------------------------------------------

gsea_go_bp_df <- as.data.frame(gsea_go_bp)

# ------------------------------------------------------------
# Save complete results
# ------------------------------------------------------------

write.csv(
  gsea_go_bp_df,
  "Results/Enrichment/GSEA/GO_BP_GSEA.csv",
  row.names = FALSE
)

saveRDS(
  gsea_go_bp,
  "Results/Enrichment/GSEA/GO_BP_GSEA.rds"
)

# ------------------------------------------------------------
# Save top 20 pathways
# ------------------------------------------------------------

top20_gsea_go_bp <- gsea_go_bp_df[
  order(gsea_go_bp_df$p.adjust),
]

top20_gsea_go_bp <- head(
  top20_gsea_go_bp,
  20
)

write.csv(
  top20_gsea_go_bp,
  "Results/Enrichment/GSEA/Top20_GO_BP_GSEA.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Dotplot
# ------------------------------------------------------------

p_gsea_go_bp <- dotplot(
  gsea_go_bp,
  showCategory = 20,
  title = "GO Biological Process GSEA"
)

ggsave(
  "Results/Enrichment/GSEA/Figures/GO_BP_GSEA_Dotplot.png",
  p_gsea_go_bp,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Enrichment plots for selected top pathways
# ------------------------------------------------------------

top_pathways <- gsea_go_bp_df$ID[
  order(gsea_go_bp_df$p.adjust)
]

top_pathways <- head(
  top_pathways,
  3
)

for (i in seq_along(top_pathways)) {
  
  pathway_id <- top_pathways[i]
  
  p <- gseaplot2(
    gsea_go_bp,
    geneSetID = pathway_id,
    title = gsea_go_bp_df$Description[
      match(pathway_id, gsea_go_bp_df$ID)
    ]
  )
  
  ggsave(
    paste0(
      "Results/Enrichment/GSEA/Figures/",
      "GO_BP_GSEA_",
      i,
      "_",
      pathway_id,
      ".png"
    ),
    p,
    width = 10,
    height = 7,
    dpi = 300
  )
}

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

cat("\n===== GSEA SAVE CHECK =====\n")

cat(
  "GO-BP GSEA CSV exists:",
  file.exists(
    "Results/Enrichment/GSEA/GO_BP_GSEA.csv"
  ),
  "\n"
)

cat(
  "GO-BP GSEA RDS exists:",
  file.exists(
    "Results/Enrichment/GSEA/GO_BP_GSEA.rds"
  ),
  "\n"
)

cat(
  "GO-BP GSEA Dotplot exists:",
  file.exists(
    "Results/Enrichment/GSEA/Figures/GO_BP_GSEA_Dotplot.png"
  ),
  "\n"
)

cat(
  "Number of GSEA pathways saved:",
  nrow(gsea_go_bp_df),
  "\n"
)

# ============================================================
# STEP 39 - Reactome GSEA
# ============================================================

library(ReactomePA)
library(enrichplot)
library(ggplot2)

# ------------------------------------------------------------
# Run Reactome GSEA
# ------------------------------------------------------------

gsea_reactome <- gsePathway(
  geneList = gene_ranking_entrez,
  organism = "human",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE
)

# ------------------------------------------------------------
# Check result
# ------------------------------------------------------------

gsea_reactome_df <- as.data.frame(gsea_reactome)

cat("\n===== REACTOME GSEA RESULT =====\n")

cat(
  "Number of enriched Reactome pathways:",
  nrow(gsea_reactome_df),
  "\n"
)

if (nrow(gsea_reactome_df) > 0) {
  
  cat("\nTop 10 Reactome pathways:\n")
  
  print(
    gsea_reactome_df[
      1:min(10, nrow(gsea_reactome_df)),
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
  
} else {
  
  cat("\nNo significant Reactome pathways detected.\n")
}

# ============================================================
# STEP 40 - Save Reactome GSEA Results + Figures
# ============================================================

# ------------------------------------------------------------
# Save complete Reactome GSEA results
# ------------------------------------------------------------

write.csv(
  gsea_reactome_df,
  "Results/Enrichment/GSEA/Reactome_GSEA.csv",
  row.names = FALSE
)

saveRDS(
  gsea_reactome,
  "Results/Enrichment/GSEA/Reactome_GSEA.rds"
)

# ------------------------------------------------------------
# Top 20 Reactome pathways
# ------------------------------------------------------------

top20_gsea_reactome <- gsea_reactome_df[
  order(gsea_reactome_df$p.adjust),
]

top20_gsea_reactome <- head(
  top20_gsea_reactome,
  20
)

write.csv(
  top20_gsea_reactome,
  "Results/Enrichment/GSEA/Top20_Reactome_GSEA.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Reactome GSEA Dotplot
# ------------------------------------------------------------

p_gsea_reactome <- dotplot(
  gsea_reactome,
  showCategory = 20,
  title = "Reactome GSEA"
)

ggsave(
  "Results/Enrichment/GSEA/Figures/Reactome_GSEA_Dotplot.png",
  p_gsea_reactome,
  width = 10,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

cat("\n===== REACTOME GSEA SAVE CHECK =====\n")

cat(
  "Reactome GSEA CSV exists:",
  file.exists(
    "Results/Enrichment/GSEA/Reactome_GSEA.csv"
  ),
  "\n"
)

cat(
  "Reactome GSEA RDS exists:",
  file.exists(
    "Results/Enrichment/GSEA/Reactome_GSEA.rds"
  ),
  "\n"
)

cat(
  "Reactome GSEA Dotplot exists:",
  file.exists(
    "Results/Enrichment/GSEA/Figures/Reactome_GSEA_Dotplot.png"
  ),
  "\n"
)

cat(
  "Number of Reactome pathways saved:",
  nrow(gsea_reactome_df),
  "\n"
)

# ============================================================
