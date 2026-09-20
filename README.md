# GSE212041 Longitudinal Multi-Omics Analysis

![R](https://img.shields.io/badge/R-4.6.1-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![MOFA2](https://img.shields.io/badge/MOFA2-multi--omics-orange)

## Integrating Neutrophil Transcriptomics and Plasma Proteomics Using MOFA2

**Author:** Ibrahim Ashraf  
**Dataset:** NCBI GEO GSE212041  
**R Version:** 4.6.1  
**Project Type:** Longitudinal Multi-Omics Computational Re-analysis

---

## Overview

This project presents a reproducible computational re-analysis of publicly available longitudinal molecular data associated with GSE212041, integrating neutrophil transcriptomics with plasma Olink proteomics using MOFA2.

The main goal was to move from individual genes and proteins toward coordinated molecular programs that can be studied across molecular features, biological pathways, longitudinal timepoints, and clinical severity.

The workflow combines RNA-seq analysis, differential expression, functional enrichment, GSEA, core gene and network analysis, longitudinal analysis, plasma proteomics, metadata integration, and multi-omics factor analysis.

---

## Scientific Question

Can longitudinal integration of neutrophil transcriptomics and plasma proteomics reveal coordinated molecular programs associated with COVID-19 severity and temporal molecular remodeling?

Rather than focusing only on individual biomarkers, the project investigates broader molecular patterns that change over time and show relationships with clinical severity.

---

## Data Sources & Provenance

This project integrates data from **two distinct sources**:

### 1. Transcriptomic Data (RNA-seq)
- **Source:** NCBI GEO
- **Accession:** [GSE212041](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212041)
- **Data:** Gene expression counts, sample metadata, clinical annotations.

### 2. Plasma Olink Proteomics & Clinical Metadata
- **Source:** Filbin et al. (Supplementary Tables)
- **Paper Title:** *"Plasma proteomics reveals tissue-specific cell death and mediators of cell-cell interactions in severe COVID-19 patients"*
- **Data Used:**
  - Olink Proteomics (Supplementary Table 2)
  - Clinical Metadata (Supplementary Table 1)
  - Olink Assay Annotation
  - Additional Supplementary Tables (3–7)

### 3. Data Integration
- RNA-seq and Proteomics datasets were matched at the **patient-timepoint** level.
- Final integrated dataset: **631 patient-timepoint observations** from **303 patients** (D0 = 298, D3 = 206, D7 = 127).

> **Note:** Raw RNA-seq data files are not included due to GitHub file size limitations. Please download them directly from GEO.

---

## Statistical Design

The project applies **two complementary statistical designs**:

### 1. Cross-Sectional Analysis
- Comparison: **COVID-positive vs COVID-symptomatic controls**
- Performed **separately per timepoint** where relevant.
- Avoids confounding between COVID status and timepoint.

### 2. Longitudinal Analysis (Within-Patient)
- Paired comparisons: **D3 vs D0**, **D7 vs D0**, **D7 vs D3**
- Complete-case cohort: **115 patients** with D0/D3/D7 measurements
- Within-patient trajectory analysis, PCA trajectory, Friedman tests, and pathway recurrence.
- Leverages repeated measures from the same patient.

---

## Workflow

![Workflow](01_Project_Overview/workflow_diagram.png)

R Basics → Data Manipulation → Data Visualization → Statistics
→ RNA-seq QC → Differential Expression → Enrichment → GSEA
→ Core Gene Signature → PPI / Network Analysis
→ Longitudinal Transcriptomics → Plasma Proteomics
→ Proteomic QC → Proteomic Differential Analysis → Proteomic Enrichment
→ Metadata Integration → RNA–Protein Matching
→ MOFA2 Integration → Factor Analysis → Feature-Weight Analysis
→ Clinical Associations → Longitudinal Factor Trajectories
→ Biological Interpretation

`

---

## 1. RNA-seq AnalysiCOVID-positive vs COVID-symptomatic controls**20,044 genes**.

The main comparison was: **COVID-positive vs COVID-symptomatic controls**.
The differential-expression analysis identified:
- 458 significant DEGs
- 406 upregulated genes
- 52 downregulated genes

Results were explored using PCA, Volcano plots, Heatmaps, and statistical summaries.

---

## 2. Functional Enrichment & GSEA

Differentially expressed genes were investigated using GO (BP/MF/CC), KEGG, and Reactome. Ranked gene statistics were analyzed using GSEA.

Main biological programs:
- Antiviral responses
- Interferon signaling
- Innate immune responses
- Inflammatory processes
- Cell-cycle regulation
- DNA replication
- Chromatin-related processes

Strong Reactome GSEA signals: G2/M Checkpoints, Cell Cycle Checkpoints, DNA Replication.

---

## 3. Core Gene Signature & Network Analysis

Leading-edge gene analysis produced:
- 159 unique leading-edge genes
- 54 genes shared across the three major pathways
- 12 significantly upregulated core genes:
  CCNA2, ORC1, CCNA1, CDC6, CDC45, MCM4, MCM2, MCM10, H2BC5, H2BC17, H2BC9, H2BC7

These represent:
- Cell Cycle: CCNA1, CCNA2
- DNA Replication: CDC6, CDC45, MCM2, MCM4, MCM10, ORC1
- Histone / Chromatin: H2BC5, H2BC7, H2BC9, H2BC17

The 12-gene signature was analyzed using STRINGdb PPI analysis, producing 43 STRING interactions with degree, betweenness, and closeness centrality.

---

## 4. Longitudinal Transcriptomic Analysis

The COVID-positive longitudinal cohort contained:
- 635 samples
- 304 patients
- D0, D3, and D7 measurements
- 115 patients with complete D0/D3/D7 measurements

Main within-patient comparisons:

| Comparison | Significant genes |
|------------|-------------------|
| D3 vs D0   | 747               |
| D7 vs D0   | 2,355             |
| D7 vs D3   | 44                |

The largest transcriptomic remodeling occurred between D0 and D7.

Longitudinal enrichment highlighted antiviral, interferon, innate immunity, inflammatory, neutrophil-related, and disease-associated programs.

---

## 5. Plasma Proteomics Analysis

The Olink dataset contained:
- 784 samples
- 1,429 protein assays

Proteomic QC included: numeric validation, assay completeness, missing-value assessment, sample-level QC, UniProt mapping, clinical metadata matching, and duplicate patient-timepoint assessment. Overall missingness was extremely low.

---

## 6. Proteomic Differential Analysis

Baseline D0 comparison: COVID-positive vs COVID-negative.

Across 1,429 protein assays, the analysis identified 42 significant proteins, including immune- and inflammatory-associated proteins:
CXCL10, CXCL11, CCL7, CCL8, CCL16, CCL24

---

## 7. Proteomic Functional Enrichment

The strongest enrichment was associated with chemokine-mediated signaling. Seven proteins contributed:
CXCL10, CCL7, CCL8, CXCL11, CCL16, CCL24, TFF2

---

## 8. RNA–Protein Metadata Integration

The transcriptomic and proteomic datasets were integrated using clinical and sample metadata (Patient ID, Sample ID, GEO accession, COVID status, patient category, timepoint, acuity, cell type).

Final matched dataset: 631 patient-timepoint observations from 303 patients (D0 = 298, D3 = 206, D7 = 127).

---

## 9. Multi-Omics Integration Using MOFA2

MOFA2 (Multi-Omics Factor Analysis) integrated two molecular views:
- RNA: 2,000 highly variable genes
- Protein: 1,000 highly variable proteins

Trained on 631 matched patient-timepoint observations, generating 15 latent molecular factors.

---

## 10. MOFA2 Variance Decomposition

| Molecular View | Variance Explained |
|----------------|-------------------|
| RNA            | 48.42%            |
| Protein        | 62.50%            |

---

## 11. MOFA2 Downstream Analysis

Analyses included: factor visualization, factor combinations, factor scores, factor trajectories, feature weights, RNA/protein contribution analysis, clinical associations, GO-BP enrichment, heatmaps, and scatter plots.

---

## 12. Key MOFA2 Factors

Six factors showed convergent evidence:
- Factor 5: strong positive association with clinical acuity (Spearman ρ ≈ +0.654). Enrichment: adaptive immune response, immune system process.
- Factor 1: strong negative association with clinical acuity (Spearman ρ ≈ −0.500). Enrichment: innate immune response, defense response to virus.
- Factor 4: strongest longitudinal increase (D3−D0 ≈ +1.41, D7−D0 ≈ +2.19, D7−D3 ≈ +0.72).
- Factor 2: strong RNA-level contribution and longitudinal behavior.
- Factor 3: longitudinal changes with innate immunity, viral response, TNF-related processes.
- Factor 6: protein-associated variation with LPS and bacterial defense response.

---

## 13. MOFA2 Feature Weights

Feature weights were used to identify the molecular measurements contributing most strongly to each factor. Top RNA and protein features were extracted and visualized. Dedicated RNA and protein feature-weight heatmaps are included in the repository.

---

## 14. MOFA2 GO-BP Enrichment

MOFA2 factor weights were analyzed using GO Biological Process gene sets. Across the six key factors, the analysis identified 204 significant factor–pathway associations with immune, antiviral, inflammatory, and disease-related themes.

---

## 15. Clinical Associations & Longitudinal Factor Trajectories

MOFA2 factor scores were evaluated against clinical acuity and across D0 → D3 → D7. Significant longitudinal changes were observed among key factors, shifting interpretation from individual genes/proteins toward coordinated multi-omic programs.

---

## 16. Biological Interpretation

The combined analyses highlighted:
- Antiviral and Interferon Responses
- Innate and Adaptive Immunity
- Inflammatory Signaling (including TNF-related)
- Neutrophil Molecular Remodeling across timepoints
- Cell Cycle and DNA Replication (core transcriptomic signature)
- RNA–Protein Coordination via MOFA2 latent programs

---

## Main Results

| Analysis | Result |
|----------|--------|
| Genes analyzed | 20,044 |
| Significant primary DEGs | 458 |
| Upregulated DEGs | 406 |
| Downregulated DEGs | 52 |
| Shared leading-edge genes | 54 |
| Core significant genes | 12 |
| STRING interactions | 43 |
| Longitudinal RNA samples | 635 |
| Longitudinal patients | 304 |
| Complete D0/D3/D7 patients | 115 |
| D3 vs D0 significant genes | 747 |
| D7 vs D0 significant genes | 2,355 |
| D7 vs D3 significant genes | 44 |
| Proteomics samples | 784 |
| Olink assays | 1,429 |
| Significant D0 proteins | 42 |
| Matched RNA–protein observations | 631 |
| Matched patients | 303 |
| MOFA2 RNA features | 2,000 |
| MOFA2 protein features | 1,000 |
| MOFA2 factors | 15 |
| Key MOFA2 factors | 6 |
| RNA variance explained | 48.42% |
| Protein variance explained | 62.50% |

---

## Repository Structure

| Folder | Contents |
|--------|----------|
| 01_Project_Overview/ | Workflow diagram, METHODS.md |
| 02_Data/ | RNA-seq and Proteomics data folders |
| 03_Metadata/ | Clinical and sample metadata (CSV) |
| 04_Analysis/ | Analysis results (Differential Expression, Enrichment, Network, Proteomics) |
| 05_Main_Figures/ | Main figures organized by analysis type (RNA, Proteomics, MOFA2, Network, PCA) |
| 06_Supplementary_Figures/ | Supplementary figures |
| 07_Network/ | Cytoscape and STRING network files |
| 08_Scripts/ | Modular R scripts (00 to 13) |
| 09_References/ | Bibliography and reference papers |
| README.md | Project documentation |
| LICENSE | MIT License |
| .gitignore | Files excluded from the repository |

---

## Analysis Scripts

The workflow is organized into modular R scripts in 08_Scripts/:

| Script | Purpose |
|--------|---------|
| 00_Project_Setup.R | Libraries and project setup |
| 01_RNA_Preprocessing_QC.R | RNA-seq QC and preprocessing |
| 02_RNA_Differential_Expression.R | Primary DE analysis (COVID vs symptomatic) |
| 03_RNA_Enrichment_GSEA.R | Functional enrichment and GSEA |
| 04_Core_Signature_PPI.R | Core gene signature and PPI network |
| 05_Longitudinal_RNA_Analysis.R | Longitudinal transcriptomic analysis |
| 06_Proteomics_QC.R | Olink proteomics QC |
| 07_Proteomics_Differential_Analysis.R | Proteomic DE analysis |
| 08_Proteomics_Enrichment.R | Proteomic enrichment |
| 09_RNA_Proteomics_Matching.R | Patient-timepoint matching |
| 10_MOFA2_Preparation.R | MOFA2 input preparation |
| 11_MOFA2_Training.R | MOFA2 model training |
| 12_MOFA2_Downstream_Analysis.R | Factor analysis and interpretation |
| 13_Final_Figures_Tables.R | Final figures and tables |

---

## Software & Methods

Developed primarily in R 4.6.1.

Major packages: GEOquery, DESeq2, edgeR, limma, ggplot2, pheatmap, dplyr, readxl, clusterProfiler, ReactomePA, org.Hs.eg.db, AnnotationDbi, STRINGdb, enrichplot, MOFA2, reticulate, Python / mofapy2.

Workflow: RNA-seq → Differential Expression → Enrichment → GSEA → Network Analysis → Longitudinal Analysis → Proteomics → Metadata Integration → Multi-Omics Integration → MOFA2 → Biological Interpretation.

---

## Project Motivation

The main purpose of this project was to transform a self-learning journey in R into a complete scientific workflow using real biological data and real clinical/sample metadata.

The project gradually connected R programming, data manipulation, visualization, statistics, RNA-seq analysis, proteomics, longitudinal analysis, network analysis, multi-omics integration, and biological interpretation — moving beyond tutorials toward an end-to-end computational biology project.

---

## Conclusion

This project integrates longitudinal neutrophil transcriptomics and plasma proteomics to investigate molecular changes associated with COVID-19 and clinical severity.

The workflow moves from individual genes and proteins to pathways, networks, longitudinal changes, and finally integrated latent molecular programs using MOFA2.

The combined results highlight coordinated molecular patterns involving antiviral and interferon responses, innate and adaptive immunity, inflammatory signaling, neutrophil remodeling, cell-cycle and DNA-replication programs, RNA–protein coordination, and longitudinal molecular changes associated with clinical severity.

This is my first project — but definitely not my last.

---

##  Full Repository

[GSE212041_MOFA2_Project](https://github.com/ibrahimayousef911-lab/GSE212041_MOFA2_Project)

Feedback, suggestions, scientific discussion, and constructive criticism are welcome.

From learning R to analyzing real biological data.  
From writing code to building a complete project.
`

