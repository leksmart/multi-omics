# Multi-Omics Integration Pipeline: SNV, RNA-Seq & DNA Methylation
## Overvew
This repository contains a reproducible and modular bioinformatics workflow for integrating single nucleotide variants (SNVs), RNA sequencing (RNA-Seq), and DNA methylation (DNAm) data. The pipeline is designed to facilitate multi-omics analysis in cancer and other complex diseases, enabling discovery of driver mutations, regulatory alterations, and prognostic biomarkers.
The workflow supports data preprocessing, quality control, normalization, feature selection, integrative analysis, and visualization, and is built using Nextflow, R, and Python.

## Key Features
### Preprocessing

SNV annotation (using ANNOVAR or VEP)

RNA-Seq normalization (e.g., TPM, DESeq2)

DNA methylation filtering and β-value transformation

### Data Integration

Sample-matched omics alignment

Common identifier harmonization (e.g., HGNC, Ensembl)

Matrix fusion and dimensionality reduction (e.g., PCA, MOFA, or SNF)

### Analysis Modules

Mutation impact on gene expression and methylation

Differential expression and differential methylation analysis

Multi-omics clustering and subtype discovery

Prognostic biomarker identification (CoxPH, LASSO)

Visualization

Heatmaps, volcano plots, oncoprints

