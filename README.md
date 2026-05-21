# README

## Overview

This repository contains analysis scripts used for scRNA-seq and scATAC-seq analyses of control and Dnd1 conditional knockout (cKO) embryonic testes.

The analyses were performed mainly using Seurat, ArchR, and Monocle3.

---

## Repository Structure

````text
.
├── 1_scRNA dataset process and annotation_Dnd1cKO.R     # Complehensive scRNA-seq preprocessing pipeline
├── 2_GastrulationData.R        # MouseGastrulationData preprocessing
├── 3_scATAC_dataset_process.R # Comprehensive scATAC-seq preprocessing pipeline
├── Fig1_S1.R                  # Figure 1 and Supplementary Figure 1 analyses
├── Fig2_S2.R                  # Figure 2 and Supplementary Figure 2 analyses
├── Fig3_S3.R                  # Figure 3 and Supplementary Figure 3 analyses
├── Fig4_S4.R                  # Figure 4 and Supplementary Figure 4 analyses
├── Fig5_S5.R                  # Figure 5 and Supplementary Figure 5 analyses
├── Fig6_S6.R                  # Figure 6 and Supplementary Figure 6 analyses
└── README.md
````

---

## Requirements

### Main R packages

* Seurat
* ArchR
* monocle3
* Signac
* MouseGastrulationData
* ggplot2
* dplyr
* ComplexHeatmap
* velocyto.R


# Analysis Description

## Figure 1 and Supplementary Figure 1

Script:

```text
Fig1_S1.R
```

Main analyses:

* UMAP visualization of all annotated cell populations
* Marker gene FeaturePlot visualization
* Cell type annotation
* Differentially expressed gene (DEG) analysis
* Heatmap generation using top marker genes

## Figure 2 and Supplementary Figure 2

Script:

```text
Fig2_S2.R
```

Main analyses:

* Stage-specific UMAP visualization of Dnd1-cKO PGCs
* Cell number summaries across developmental stages
* DotPlot analyses of germ cell and pluripotency markers
* Comparative marker expression analyses across E14.5–E17.5

## Mouse Gastrulation Dataset Preprocessing

Script:

```text
2_GastrulationData.R
```

Main analyses:

* Download and preprocessing of MouseGastrulationData datasets
* Conversion from SingleCellExperiment objects to Seurat objects
* Removal of doublets and stripped cells
* Selection of primitive streak and epiblast populations
* Preparation of gastrulation reference datasets for integration analyses

## Figure 3 and Supplementary Figure 3

Script:

```text
Fig3_S3.R
```

Main analyses:

* Integration of current scRNA-seq datasets with published reference datasets
* Cross-study comparative analysis
* UMAP visualization after integration
* Cluster composition analysis
* Relative proportion analysis across developmental stages and culture conditions


## Figure 4 and Supplementary Figure 4

Script:

```text
Fig4_S4.R
```

Main analyses:

* Marker gene analysis in E16.5 cKO germ cells
* DotPlot visualization of Car4-related genes
* Two-way ANOVA
* Tukey–Kramer post hoc testing
* Boxplot visualization of quantitative experimental results


## scATAC-seq Dataset Processing Pipeline

Script:

```text
3_scATAC_dataset_process.R
```

Main analyses:

* Creation of ArchR Arrow files
* Doublet detection and filtering
* Iterative LSI dimensional reduction
* UMAP visualization and clustering
* Harmony batch correction
* RNA–ATAC integration
* Peak calling using MACS2
* Motif enrichment analysis
* PGC subclustering

Main analyses:

* ArchR preprocessing
* Clustering and UMAP analysis
* RNA–ATAC integration
* Peak calling
* Motif enrichment analysis

## Figure 5 and Supplementary Figure 5

Script:

```text
Fig5_S5.R
```

Main analyses:

* scRNA-seq analysis of E13.5–E15.5 PGCs
* Integration of published E13.5 PGC datasets with current datasets
* RNA velocity analysis using velocyto.R
* Cluster composition analysis
* DotPlot visualization of pluripotency and PGC markers

Reference dataset:

* Nguyen D.H. et al., 2021


## Figure 6 and Supplementary Figure 6

Script:

```text
Fig6_S6.R
```

Main analyses:

* scATAC-seq trajectory analysis using Monocle3
* Pseudotime analysis of PGC populations
* RNA–ATAC integration
* Gene score visualization
* Motif enrichment heatmap generation
* Cell-type annotation of scATAC-seq clusters
* Quality control metrics for fragment counts and TSS enrichment