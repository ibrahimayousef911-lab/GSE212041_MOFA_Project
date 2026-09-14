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
