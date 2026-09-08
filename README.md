# Single-Cell-RNA-seq-Analysis-Using-R
This repository contains the computational workflow used for the analysis of  scRNA-seq data using R.  The analysis was done using NCBI geo dataset (GSE213625) and the analysis includes quality control, normalization, dimensionality reduction, cell-type annotation, marker-gene identification, differential expression analysis, and visualization.
The primary objective of this analysis was to identify distinct cellular populations and characterize their transcriptional profiles based on cluster-specific marker genes.

# Analysis Workflow

Raw scRNA-seq data
        │
Quality Control
        │
Normalization
        │
Highly Variable Gene Identification
        │
PCA
        │
Nearest-neighbor Graph
        │
Clustering
        │
UMAP Visualization
        │
Cell-type Annotation
        │
Marker Gene Identification
        │
Differential Expression Analysis
        │
Visualization & Biological Interpretation

# Packages required 
Seurat
ggplot2
dplyr
patchwork
SingleR
celldex
RColorBrewer

# Data
The raw sequencing data are not included in this repository because of their size.
Raw scRNA-seq data are available from:
Database: GEO GSE213625

# Input files
barcode.tsv
feature.tsv
matrix.mtx

# Output files
1. QC plots- Assessment of cell quality
2. PCA- Dimensionality reduction
3. UMAP- Visualization of cell populations
4. Cluster assignments-	Unsupervised cell grouping
5. Cell-type annotations-	Biological identification of clusters
6. Marker genes-	Cluster-specific genes
7. Differential expression-	Comparison between populations
8. Volcano plots-	Visualization of differential expression
