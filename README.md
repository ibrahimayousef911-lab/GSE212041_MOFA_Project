# GSE212041 Longitudinal Multi-Omics Analysis

## Integrating Neutrophil Transcriptomics and Plasma Proteomics Using MOFA2

**Author:** Ibrahim Ashraf  
**Dataset:** NCBI GEO GSE212041  
**R Version:** 4.6.1  
**Project Type:** Longitudinal Multi-Omics Computational Re-analysis

---

## Overview

This project presents a reproducible computational re-analysis of publicly available longitudinal molecular data associated with **GSE212041**, integrating neutrophil transcriptomics with plasma Olink proteomics using **MOFA2**.

The main goal was to move from individual genes and proteins toward **coordinated molecular programs** that can be studied across molecular features, biological pathways, longitudinal timepoints, and clinical severity.

The workflow combines RNA-seq analysis, differential expression, functional enrichment, GSEA, core gene and network analysis, longitudinal analysis, plasma proteomics, metadata integration, and multi-omics factor analysis.

---

## Scientific Question

> **Can longitudinal integration of neutrophil transcriptomics and plasma proteomics reveal coordinated molecular programs associated with COVID-19 severity and temporal molecular remodeling?**

Rather than focusing only on individual biomarkers, the project investigates broader molecular patterns that change over time and show relationships with clinical severity.

---

## Data Sources

### Transcriptomic Data

The RNA-seq component was obtained from the public **NCBI GEO dataset GSE212041**.

The transcriptomic data and associated sample metadata provided information including:

- Gene expression counts
- GEO sample accession
- Patient/sample identifiers
- COVID-19 status
- Patient category
- Timepoint
- Clinical acuity
- Cell-type information

The primary transcriptomic analysis contained:

**20,044 genes**

---

### Plasma Proteomics Data

The plasma proteomics component was obtained from the **Olink plasma proteomics data associated with the GSE212041 study**.

The proteomics dataset contained:

- **784 samples**
- **1,429 Olink protein assays**

The corresponding Olink assay annotation was used to connect assay identifiers with:

- Assay names
- UniProt identifiers
- Olink panel information
- Panel version information

**OlinkID** was retained as the primary assay identifier in order to preserve assay-level information and avoid incorrectly collapsing measurements from different panels.

The proteomics data were independently quality-controlled and then integrated with the clinical/sample metadata.

---

### Clinical and Sample Metadata

Clinical and sample metadata were used to connect the transcriptomic and proteomic measurements.

The metadata included:

- Patient ID
- Sample ID
- GEO accession
- COVID-19 status
- Patient category
- Timepoint
- Clinical acuity
- Cell type
- Sample/draw information

This metadata layer was essential for identifying repeated measurements from the same patients and constructing the longitudinal and multi-omics cohorts.

---

## Data Integration Strategy

The project integrates the three main data layers:

```text
GSE212041 RNA-seq
        +
Olink Plasma Proteomics
        +
Clinical / Sample Metadata
        ↓
Patient-Level & Timepoint Matching
        ↓
Integrated RNA–Protein Dataset
        ↓
MOFA2 Multi-Omics Analysis
The RNA and protein datasets were matched using a patient-timepoint identifier.

The final matched multi-omics dataset contained:

631 patient-timepoint observations

from:

303 patients

with:
D0 = 298
D3 = 206
D7 = 127
The original RNA and proteomics datasets were retained separately for their individual analyses, while the matched observations were used for multi-omics integration.

##Analysis Workflow
R Basics
   ↓
Data Manipulation
   ↓
Data Visualization
   ↓
Statistics
   ↓
RNA-seq Quality Control
   ↓
Differential Expression
   ↓
Functional Enrichment
   ↓
GSEA
   ↓
Core Gene Signature
   ↓
PPI / Network Analysis
   ↓
Longitudinal Transcriptomics
   ↓
Plasma Proteomics
   ↓
Proteomic Quality Control
   ↓
Proteomic Differential Analysis
   ↓
Proteomic Enrichment
   ↓
Metadata Integration
   ↓
RNA–Protein Matching
   ↓
MOFA2 Multi-Omics Integration
   ↓
Factor Analysis
   ↓
Feature-Weight Analysis
   ↓
Clinical Associations
   ↓
Longitudinal Factor Trajectories
   ↓
Biological Interpretation

1. RNA-seq Analysis

The transcriptomic analysis began with sample and metadata validation, RNA-seq quality control, and exploratory analysis.

The primary transcriptomic dataset contained:

20,044 genes

The main comparison was:

COVID-positive vs COVID-symptomatic controls

The differential-expression analysis identified:

458 significant DEGs
406 upregulated genes
52 downregulated genes

The transcriptomic results were explored using:

PCA
Volcano plots
Differential-expression heatmaps
Statistical summaries
2. Functional Enrichment & GSEA

Differentially expressed genes were investigated using:

Gene Ontology Biological Process
Gene Ontology Molecular Function
Gene Ontology Cellular Component
KEGG
Reactome

Ranked gene-level statistics were also analyzed using Gene Set Enrichment Analysis (GSEA).

The main biological programs identified across the transcriptomic analyses included:

Antiviral responses
Interferon signaling
Innate immune responses
Inflammatory processes
Cell-cycle regulation
DNA replication
Chromatin-related processes

Strong Reactome GSEA signals included:

G2/M Checkpoints
Cell Cycle Checkpoints
DNA Replication
3. Core Gene Signature & Network Analysis

Leading-edge genes from the major recurring Reactome pathways were compared to identify a focused molecular signature.

This produced:

159 unique leading-edge genes

with:

54 genes shared across the three major pathways.

Among these, 12 genes were significantly upregulated:
CCNA2
ORC1
CCNA1
CDC6
CDC45
MCM4
MCM2
MCM10
H2BC5
H2BC17
H2BC9
H2BC7

The signature mainly represented:

Cell Cycle
CCNA1
CCNA2
DNA Replication
CDC6
CDC45
MCM2
MCM4
MCM10
ORC1
Histone / Chromatin
H2BC5
H2BC7
H2BC9
H2BC17

The 12-gene signature was further investigated using STRINGdb and protein–protein interaction analysis.

The resulting network contained:

43 STRING interactions

with degree, betweenness, and closeness centrality used to characterize network topology.

4. Longitudinal Transcriptomic Analysis

The transcriptomic analysis was extended from cross-sectional comparisons to patient-level longitudinal analysis.

The COVID-positive longitudinal cohort contained:

635 samples
304 patients
D0, D3, and D7 measurements
115 patients with complete D0/D3/D7 measurements

The main within-patient comparisons were:

Comparison	Significant genes
D3 vs D0	747
D7 vs D0	2,355
D7 vs D3	44

The largest transcriptomic remodeling occurred between D0 and D7.

Longitudinal enrichment highlighted recurring biological programs involving antiviral responses, interferon signaling, innate immunity, inflammatory signaling, neutrophil-related processes, and disease-associated molecular pathways.

5. Plasma Proteomics Analysis

The plasma proteomics component provided a complementary molecular layer to the transcriptomic analysis.

The Olink dataset contained:

784 samples

and:

1,429 protein assays

Proteomic quality control included:

Numeric validation
Assay completeness
Missing-value assessment
Sample-level QC
Assay annotation
UniProt mapping
Clinical metadata matching
Duplicate patient-timepoint assessment

The overall missingness was extremely low.

The Olink assays were retained using their unique OlinkID identifiers to preserve assay-level information.

6. Proteomic Differential Analysis

Baseline plasma proteomic differences were analyzed at D0 between:

COVID-positive vs COVID-negative samples

Across the:

1,429 protein assays

the analysis identified:

42 significant proteins

under the project-specific significance and effect-size criteria.

The significant proteins included immune- and inflammatory-associated proteins such as:
CXCL10
CXCL11
CCL7
CCL8
CCL16
CCL24

7. Proteomic Functional Enrichment

The significant proteomic features were mapped to biological processes.

The strongest enrichment was associated with:

chemokine-mediated signaling and cellular responses to chemokines.

Seven proteins contributed to the enriched GO biological processes:
CXCL10
CCL7
CCL8
CXCL11
CCL16
CCL24
TFF2
These were interpreted as proteins contributing to the enriched biological processes rather than treating every contributing protein as a chemokine.


8. RNA–Protein Metadata Integration

The transcriptomic and proteomic datasets were integrated using the available clinical and sample metadata.

The matching process used:

Patient ID
Sample ID
GEO accession
COVID status
Patient category
Timepoint
Acuity
Cell type

The two molecular layers were matched using a patient-timepoint identifier.

The final matched dataset contained:

631 patient-timepoint observations

from:

303 patients

with:
D0 = 298
D3 = 206
D7 = 127
The RNA and protein observations were aligned using identical identifiers before multi-omics modeling.

9. Multi-Omics Integration Using MOFA2

The central integration analysis was performed using:

MOFA2 — Multi-Omics Factor Analysis

Two molecular views were integrated:
RNA
2,000 highly variable genes

        +

Protein
1,000 highly variable proteins

        ↓

MOFA2

        ↓

15 latent factors
The model was trained using:

631 matched patient-timepoint observations

and generated:

15 latent molecular factors

representing coordinated variation across the transcriptomic and proteomic layers.

10. MOFA2 Variance Decomposition

The trained model explained:

Molecular View	Variance Explained
RNA	48.42%
Protein	62.50%

Variance decomposition was examined at both the global and individual-factor levels to understand the contribution of different latent factors to each molecular view.

11. MOFA2 Downstream Analysis

The trained MOFA2 model was investigated using downstream analyses including:

Factor visualization
Factor combinations
Factor scores
Factor trajectories
Feature weights
RNA/protein contribution analysis
Clinical associations
GO biological-process enrichment
Heatmap visualization
Scatter plots

This connected latent factors with their molecular features, biological pathways, clinical associations, and longitudinal behavior.

12. Key MOFA2 Factors

Six factors showed convergent evidence across variance contribution, longitudinal behavior, clinical association, and/or functional enrichment:
Factor 5
Factor 1
Factor 6
Factor 3
Factor 4
Factor 2
These factors represent different components of the integrated molecular structure rather than a simple ranking.

Factor 5

Factor 5 showed a strong positive association with clinical acuity:

Spearman ρ ≈ +0.654

Its biological enrichment included:

Adaptive immune response
Immune system process
Immune response
Factor 1

Factor 1 showed a strong negative association with clinical acuity:

Spearman ρ ≈ −0.500

Its biological enrichment included:

Innate immune response
Defense response to virus
Response to virus
Negative regulation of viral genome replication
Factor 4

Factor 4 showed the strongest longitudinal increase among the key factors:
D3 − D0 ≈ +1.41
D7 − D0 ≈ +2.19
D7 − D3 ≈ +0.72
This factor therefore captured a strong temporal component of the integrated molecular response.

Other Key Factors

Factor 2 showed strong RNA-level contribution and significant longitudinal behavior.

Factor 3 showed longitudinal changes and enrichment involving innate immunity, viral response, and TNF-related processes.

Factor 6 showed protein-associated variation, clinical association, and enrichment involving cellular responses to lipopolysaccharide and bacterial defense.

13. MOFA2 Feature Weights

Feature weights were used to identify the molecular measurements contributing most strongly to each factor.

For every factor, the top RNA and protein features were extracted and visualized.

The interpretation follows:
MOFA2 Factor
     ↓
Top RNA Features
     +
Top Protein Features
     ↓
Biological Pathways
     ↓
Clinical / Longitudinal Pattern
The repository contains dedicated RNA and protein feature-weight heatmaps.

14. MOFA2 GO-BP Enrichment

MOFA2 factor weights were analyzed using GO Biological Process gene sets.

This analysis is distinct from the earlier DEG enrichment.

The earlier enrichment asks:

Which pathways are represented among differentially expressed genes?

The MOFA2 enrichment asks:

Which biological processes are represented among the molecular features contributing strongly to a latent factor?

Across the six key factors, the analysis identified:

204 significant factor–pathway associations

with biological themes involving immune, antiviral, inflammatory, and other disease-related processes.

15. Clinical Associations & Longitudinal Factor Trajectories

MOFA2 factor scores were evaluated against clinical acuity and across longitudinal timepoints.

The main trajectory framework was:
D0 → D3 → D7
Significant longitudinal changes were observed among the key factors.

This analysis shifts the interpretation from individual genes and proteins toward coordinated multi-omic molecular programs changing over time.

16. Biological Interpretation

The combined transcriptomic, proteomic, longitudinal, and MOFA2 analyses highlighted several major biological themes.

Antiviral and Interferon Responses

Strong antiviral and interferon-associated programs were identified across transcriptomic and integrated analyses.

Innate and Adaptive Immunity

Multiple molecular layers showed coordinated immune-related processes involving innate and adaptive responses.

Inflammatory Signaling

Inflammatory and TNF-associated programs were repeatedly observed across the molecular analyses.

Neutrophil Molecular Remodeling

Longitudinal analyses demonstrated substantial changes in circulating neutrophil molecular states across disease timepoints.

Cell Cycle and DNA Replication

The core transcriptomic signature highlighted coordinated cell-cycle, DNA-replication, and chromatin-related programs.

RNA–Protein Coordination

MOFA2 provided a framework for identifying latent molecular programs representing coordinated variation across transcriptomic and plasma protein measurements.

Main Results
Analysis	Result
Genes analyzed	20,044
Significant primary DEGs	458
Upregulated DEGs	406
Downregulated DEGs	52
Shared leading-edge genes	54
Core significant genes	12
STRING interactions	43
Longitudinal RNA samples	635
Longitudinal patients	304
Complete D0/D3/D7 patients	115
D3 vs D0 significant genes	747
D7 vs D0 significant genes	2,355
D7 vs D3 significant genes	44
Proteomics samples	784
Olink assays	1,429
Significant D0 proteins	42
Matched RNA–protein observations	631
Matched patients	303
MOFA2 RNA features	2,000
MOFA2 protein features	1,000
MOFA2 factors	15
Key MOFA2 factors	6
RNA variance explained	48.42%
Protein variance explained	62.50%
Main Figures

The repository contains selected representative figures covering the major stages of the workflow:

RNA-seq PCA
Differential-expression volcano plot
Reactome GSEA
Core gene/PPI network
Proteomic differential-expression heatmap
MOFA2 variance explained
MOFA2 factor combination
MOFA2 factor trajectories
RNA feature-weight heatmap
Protein feature-weight heatmap
MOFA2 GO-BP enrichment
Factor–acuity associations

Additional QC and supporting figures are available in the supplementary figures directory.

##Repository Structure
GSE212041_MOFA_Project/
│
├── 01_Project_Overview/
├── 02_Data/
├── 03_Metadata/
├── 04_Results/
├── 05_Main_Figures/
├── 06_Supplementary_Figures/
├── 07_Network/
├── 08_Scripts/
│
├── README.md
├── LICENSE
└── .gitignore

##Analysis Scripts
08_Scripts/
│
├── 00_Project_Setup.R
├── 01_RNA_Preprocessing_QC.R
├── 02_RNA_Differential_Expression.R
├── 03_RNA_Enrichment_GSEA.R
├── 04_Core_Signature_PPI.R
├── 05_Longitudinal_RNA_Analysis.R
├── 06_Proteomics_QC.R
├── 07_Proteomics_Differential_Analysis.R
├── 08_Proteomics_Enrichment.R
├── 09_RNA_Proteomics_Matching.R
├── 10_MOFA2_Preparation.R
├── 11_MOFA2_Training.R
├── 12_MOFA2_Downstream_Analysis.R
└── 13_Final_Figures_Tables.R

Software & Methods

The project was developed primarily in R 4.6.1.

Major packages and tools include:

GEOquery
DESeq2
edgeR
limma
ggplot2
pheatmap
dplyr
readxl
clusterProfiler
ReactomePA
org.Hs.eg.db
AnnotationDbi
STRINGdb
enrichplot
MOFA2
reticulate
Python / mofapy2

The workflow combines:

RNA-seq → Differential Expression → Enrichment → GSEA → Network Analysis → Longitudinal Analysis → Proteomics → Metadata Integration → Multi-Omics Integration → MOFA2 → Biological Interpretation

Project Motivation

The main purpose of this project was to transform a self-learning journey in R into a complete scientific workflow using real biological data and real clinical/sample metadata.

The project gradually connected:

R programming
Data manipulation
Data visualization
Statistics
RNA-seq analysis
Proteomics
Longitudinal analysis
Network analysis
Multi-omics integration
Biological interpretation

Rather than keeping the learning process limited to tutorials or isolated code examples, the goal was to build an end-to-end computational biology project.

Conclusion

This project integrates longitudinal neutrophil transcriptomics and plasma proteomics to investigate molecular changes associated with COVID-19 and clinical severity.

The workflow moves from individual genes and proteins to pathways, networks, longitudinal changes, and finally integrated latent molecular programs using MOFA2.

The combined results highlight coordinated molecular patterns involving:

Antiviral and interferon responses
Innate and adaptive immunity
Inflammatory signaling
Neutrophil remodeling
Cell-cycle and DNA-replication programs
RNA–protein coordination
Longitudinal molecular changes associated with clinical severity

The project demonstrates how multiple biological data layers can be integrated into a reproducible computational framework for multi-omic interpretation.

From Learning R to Building a Multi-Omics Project
R
↓
Data Analysis
↓
RNA-seq
↓
Proteomics
↓
Longitudinal Analysis
↓
Network Analysis
↓
Multi-Omics Integration
↓
MOFA2
↓
Biological Interpretation
This is my first project — but definitely not my last.

Full Repository

GSE212041_MOFA_Project — GitHub Repository

Feedback, suggestions, scientific discussion, and constructive criticism are welcome.

From learning R to analyzing real biological data.
From writing code to building a complete project.
