# Single-cell ATAC-seq Analysis with Motif Discovery

This R script performs comprehensive single-cell ATAC-seq analysis using Signac and Seurat, focusing on motif analysis and chromatin accessibility patterns.

## Required Libraries

```r
# Core single-cell analysis
library(Signac)
library(Seurat)

# Genomic annotation and reference
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# Clustering and statistical analysis
library(ConsensusClusterPlus)
library(aricode)
library(corrplot)

# Data manipulation and processing
library(data.table)
library(pryr)
library(Matrix)
library(dplyr)

# Genomic data handling
library(GenomicRanges)
library(rtracklayer)
library(GenomicInteractions)

# Parallel processing
library(parallel)
library(foreach)
library(doParallel)

# Additional genomic tools
library(Rsamtools)
library(AnnotationHub)
library(stringr)
library(vroom)
library(GenomeInfoDb)

# Visualization
library(ggplot2)
library(Gviz)
library(InteractionSet)
library(patchwork)

# Motif analysis
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)

# Machine learning
library(glmnet)          # For LASSO/Ridge regression
library(randomForest)    # For random forest
library(caret)          # For model training and validation
```

## Setup and Data Loading

```r
# Set working directory
setwd("./script")

# Define data paths
data_read <- "../data/"
analysis_save <- "../analysis/"
data_save <- "../data/"

# Load genomic annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

# Optional: Remove gap regions
annotation <- annotation[!(annotation$type %in% "gap")]

# Load preprocessed Seurat object with clustering
seuratObj <- readRDS(file = paste0(data_read, "seuratObj_with_clustering.rds"))
```

## Motif Analysis Pipeline

### 1. Motif Database Setup

```r
# Set ATAC assay as default
DefaultAssay(seuratObj) <- "ATAC"

# Load position frequency matrices from JASPAR2020
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(
    collection = "CORE", 
    tax_group = "vertebrates",
    species = 9606,  # Homo sapiens
    all_versions = FALSE
  )
)

# Add motifs to Seurat object
seuratObj <- AddMotifs(
  object = seuratObj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
```

### 2. Differential Accessibility Analysis

```r
# Find differentially accessible peaks between cell types
da_peaks <- FindMarkers(
  object = seuratObj,
  ident.1 = 0,        # Cell type/cluster 0
  ident.2 = 1,        # Cell type/cluster 1
  only.pos = TRUE,    # Only upregulated peaks
  test.use = "LR",    # Logistic regression test
  latent.vars = "nCount_ATAC"  # Control for sequencing depth
)

# Filter for significant peaks
top_da_peaks <- rownames(da_peaks[da_peaks$p_val_adj < 0.05, ])
```

### 3. Motif Enrichment Analysis

```r
# Find enriched motifs in differential peaks
enriched_motifs <- FindMotifs(
  object = seuratObj,
  features = top_da_peaks
)

# Visualize top enriched motifs
MotifPlot(
  object = seuratObj,
  motifs = head(rownames(enriched_motifs))
)
```

### 4. ChromVAR Motif Activity Analysis

```r
# Run chromVAR to compute motif activities
seuratObj <- RunChromVAR(
  object = seuratObj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Switch to chromVAR assay for differential analysis
DefaultAssay(seuratObj) <- 'chromvar'

# Find differential motif activities
differential.activity <- FindMarkers(
  object = seuratObj,
  ident.1 = 0,
  ident.2 = 1,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

# Plot motifs with differential activity
MotifPlot(
  object = seuratObj,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)
```

## Analysis Overview

This pipeline performs:

1. **Data Setup**: Loads genomic annotations and preprocessed single-cell ATAC-seq data
2. **Motif Database Integration**: Incorporates JASPAR2020 motif database for human transcription factors
3. **Differential Peak Analysis**: Identifies peaks with different accessibility between cell types
4. **Motif Enrichment**: Determines which transcription factor motifs are enriched in differential peaks
5. **Activity Scoring**: Uses chromVAR to compute per-cell motif activity scores
6. **Differential Activity**: Identifies motifs with significantly different activities between cell types

## Key Features

- **Comprehensive motif analysis** using state-of-the-art tools (Signac, chromVAR)
- **Statistical rigor** with multiple testing correction and proper controls
- **Visualization capabilities** for motif logos and activity patterns
- **Scalable processing** with parallel computing support
- **Integration ready** for downstream machine learning analyses

## Workflow Summary

### Step-by-Step Process

1. **Environment Setup**
   - Load required R packages for single-cell genomics analysis
   - Set up file paths and working directories

2. **Data Preparation**
   - Load genomic annotations from Ensembl
   - Import preprocessed Seurat object with cell clustering information
   - Optional filtering of genomic gap regions

3. **Motif Database Configuration**
   - Access JASPAR2020 transcription factor motif database
   - Filter for human (vertebrate) core motifs
   - Integrate motif information into Seurat object

4. **Differential Accessibility Testing**
   - Compare chromatin accessibility between cell clusters
   - Apply statistical testing with proper controls
   - Filter for significantly differential peaks

5. **Motif Enrichment Analysis**
   - Test for over-representation of motifs in differential peaks
   - Generate motif logos and visualization plots

6. **Activity-Based Analysis**
   - Compute per-cell motif activity scores using chromVAR
   - Identify motifs with differential activity patterns
   - Create comprehensive visualization outputs

## Output Files and Results

The analysis generates several key outputs:

- **Differential peak lists**: Genomic regions with cell-type-specific accessibility
- **Enriched motif tables**: Transcription factor motifs over-represented in peaks
- **Activity score matrices**: Per-cell quantification of motif activities
- **Visualization plots**: Motif logos, activity heatmaps, and statistical summaries

## Applications

This analysis framework enables:

- **Cell type characterization** through chromatin accessibility patterns
- **Transcriptional regulator identification** driving cell fate decisions
- **Regulatory network reconstruction** using motif-peak associations
- **Comparative genomics** across different cell types or conditions
- **Integration with RNA-seq** data for multi-modal analysis

## Technical Notes

### Performance Considerations
- Uses sparse matrix representations for memory efficiency
- Supports parallel processing for computationally intensive steps
- Optimized for large-scale single-cell datasets

### Statistical Methods
- Logistic regression for differential accessibility testing
- Hypergeometric testing for motif enrichment
- Bias-corrected deviation scores for activity quantification
- Multiple testing correction using Benjamini-Hochberg procedure

### Quality Control
- Sequencing depth normalization through latent variables
- Gap region exclusion to avoid technical artifacts
- Strict significance thresholds (p_val_adj < 0.05)

This comprehensive pipeline provides a robust foundation for understanding transcriptional regulation in single-cell ATAC-seq data through motif-based analysis approaches.