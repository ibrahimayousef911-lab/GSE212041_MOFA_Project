# GSE212041 Longitudinal Multi-Omics Analysis

![R](https://img.shields.io/badge/R-4.6.1-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![MOFA2](https://img.shields.io/badge/MOFA2-multi--omics-orange)

## Integrating Neutrophil Transcriptomics and Plasma Proteomics Using MOFA2

Dataset: NCBI GEO GSE212041  
Analysis environment: Egypt  
Author: Ibrahim Ashraf
R version: 4.6.1

## Overview

This repository contains a reproducible computational re-analysis of the public GSE212041 dataset, integrating longitudinal neutrophil transcriptomics with plasma Olink proteomics using MOFA2.

The workflow combines RNA-seq quality control, differential expression, functional enrichment, GSEA, longitudinal analysis, core gene/PPI analysis, proteomic differential analysis, and multi-omics factor analysis.

## Data Availability

The raw RNA-seq data (Count and TPM) are available on GEO under accession number GSE212041.

- Link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212041

> Note: Due to file size limitations, raw data files are not included in this repository. Please download them directly from GEO.

## Main Results

### Transcriptomics

The primary RNA-seq analysis evaluated 20,044 genes and identified 458 significant DEGs between COVID-19-positive and symptomatic COVID-19-negative samples:

- 406 upregulated
- 52 downregulated

Enrichment and GSEA highlighted antiviral/interferon, innate and adaptive immune, inflammatory, cell-cycle, and DNA-replication programs.

### Core Gene Signature and PPI

Leading-edge analysis identified 54 shared core genes across recurrent Reactome cell-cycle pathways. Twelve genes were significantly differentially expressed and formed the final core signature.

The signature was further evaluated using STRING protein–protein interaction analysis and network centrality measures, highlighting connected cell-cycle, DNA-replication, and chromatin-associated genes.

### Longitudinal Transcriptomics

A COVID-19-positive longitudinal cohort contained 635 samples from 304 patients across D0, D3, and D7. 115 patients had complete measurements at all three timepoints.

Within-patient analysis identified:

- 747 genes in D3 vs D0
- 2,355 genes in D7 vs D0
- 44 genes in D7 vs D3

The longitudinal results support substantial remodeling of antiviral, interferon, innate-immune, inflammatory, and neutrophil-associated programs during hospitalization.

### Plasma Proteomics

Olink plasma proteomics contained 784 samples and 1,429 protein assays. Quality control showed very low overall missingness.

The primary D0 comparison identified 42 proteins meeting the predefined stringent differential-abundance criteria. Functional analysis highlighted chemokine-related and immune-signaling processes.

### RNA–Protein Integration

RNA and protein measurements were matched by patient and timepoint, yielding 631 common patient-timepoint observations from 303 patients.

The MOFA2 input contained:

- 2 molecular views
- 2,000 variable RNA features
- 1,000 variable protein features
- 15 latent factors

### MOFA2

The trained MOFA2 model explained 48.42% of RNA variance and 62.50% of protein variance.

Downstream analysis included variance decomposition, factor visualization, factor combinations, longitudinal factor trajectories, clinical associations, feature weights, RNA/protein feature-weight heatmaps, factor-associated GO biological processes, and covariation analysis.

Six factors showed the strongest convergent biological, clinical, longitudinal, and multi-omic evidence. Factor 5 showed the strongest positive association with disease acuity, Factor 1 showed a strong negative association, and Factor 4 showed the strongest longitudinal increase.

###Key Result
- 20,044 RNA-seq genes analyzed
- 458 significant COVID-related DEGs
- 12-gene core molecular signature
- 635 longitudinal COVID-positive RNA-seq samples
- 1,429 plasma protein assays
- 42 significant D0 plasma proteins
- 631 matched RNA–protein patient-timepoint observations
- 15 MOFA2 latent factors
- 6 key MOFA2 factors with convergent biological, clinical, and longitudinal evidence

## Biological Interpretation

The integrated results support a dynamic molecular landscape in hospitalized COVID-19 characterized by:
- antiviral and interferon-associated responses
- innate and adaptive immune activation
- inflammatory and TNF-related signaling
- neutrophil-associated molecular remodeling
- cell-cycle and DNA-replication programs
- coordinated longitudinal changes across RNA and protein measurements

MOFA2 provides a unified representation of coordinated variation across the transcriptomic and proteomic layers.

## Repository Structure

GSE212041_MOFA2_Project/
├── README.md
├── LICENSE
├── .gitignore
├── .gitattributes
├── 01_Project_Overview/
│   ├── METHODS.md
│   └── workflow_diagram.png
├── 02_Data/
│   ├── RNA-Seq/
│   ├── Proteomics/
│   └── README.md
├── 03_Metadata/
│   ├── GSE212041_metadata_all_781_samples.csv
│   └── GSE212041_metadata_clean.csv
├── 04_Analysis/
├── 05_Main_Figures/
├── 06_Supplementary_Figures/
├── 07_Network/
├── 08_Scripts/
│   ├── 00_Project_Setup.R
│   ├── 01_RNA_Preprocessing_QC.R
│   ├── 02_RNA_Differential_Expression.R
│   ├── 03_RNA_Enrichment_GSEA.R
│   ├── 04_Core_Signature_PPI.R
│   ├── 05_Longitudinal_RNA_Analysis.R
│   ├── 06_Proteomics_QC.R
│   ├── 07_Proteomics_Differential_Analysis.R
│   ├── 08_Proteomics_Enrichment.R
│   ├── 09_RNA_Proteomics_Matching.R
│   ├── 10_MOFA2_Preparation.R
│   ├── 11_MOFA2_Training.R
│   ├── 12_MOFA2_Downstream_Analysis.R
│   └── 13_Final_Figures_Tables.R
└── 09_References/
## Software

The analysis uses R 4.6.1 with CRAN/Bioconductor packages including DESeq2, edgeR, limma, GEOquery, clusterProfiler, org.Hs.eg.db, ReactomePA, enrichplot, STRINGdb, igraph, ggplot2, pheatmap, dplyr, tidyr, readxl, MOFA2, and reticulate.

MOFA2 model training uses a dedicated Python environment containing the required mofapy2 dependencies.

## Reproducibility

The organized scripts separate the workflow into focused analytical modules so that preprocessing, differential analysis, enrichment, longitudinal analysis, proteomics, integration, MOFA2 analysis, and final reporting can be maintained independently.

Raw/private input files are not included in the repository. See 02_Data/README.md for data organization and reproducibility instructions.
