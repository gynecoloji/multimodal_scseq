# Multimodal scRNA-seq and scATAC-seq Analysis Pipeline

This pipeline performs integrated analysis of single-cell RNA-seq and ATAC-seq data using Seurat and Signac packages in R.

## Table of Contents
1. [Package Loading](#package-loading)
2. [Data Loading and Setup](#data-loading-and-setup)
3. [Data Normalization](#data-normalization)
4. [Neighbor Graph Construction](#neighbor-graph-construction)
5. [Clustering Analysis](#clustering-analysis)
6. [UMAP Visualization](#umap-visualization)
7. [Clustering Evaluation](#clustering-evaluation)
8. [Data Export](#data-export)

## Package Loading

```r
# Load required packages
library(Signac)           # scATAC-seq analysis
library(Seurat)           # Single-cell analysis framework
library(EnsDb.Hsapiens.v86)
library(hdf5r)            # HDF5 file handling
library(BSgenome.Hsapiens.UCSC.hg38)
library(ConsensusClusterPlus)
library(aricode)          # Clustering evaluation metrics
library(corrplot)         # Correlation plots
library(data.table)
library(Matrix)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(AnnotationHub)
library(stringr)
library(vroom)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(cowplot)
```

## Data Loading and Setup

```r
# Set working directories
setwd("./script/")
data_read = "../data/"
analysis_save = "../analysis/"
data_save = "../data/"

# Load pre-filtered Seurat object
seuratObj = readRDS(paste0(data_read, "seuratObj_filtered.rds"))
```

## Data Normalization

### RNA Data Normalization
```r
# Normalize RNA data using SCTransform
DefaultAssay(seuratObj) <- "RNA"
seuratObj <- SCTransform(seuratObj)
seuratObj <- RunPCA(seuratObj)
```

### ATAC Data Normalization
```r
# Normalize ATAC data using TF-IDF and LSI
DefaultAssay(seuratObj) <- "ATAC"
seuratObj <- FindTopFeatures(seuratObj, min.cutoff = 5)
seuratObj <- RunTFIDF(seuratObj)
seuratObj <- RunSVD(seuratObj)
```

## Neighbor Graph Construction

### Joint Multi-modal Neighbors
```r
# Calculate weighted nearest neighbors from both modalities
DefaultAssay(seuratObj) <- "SCT"
seuratObj <- FindMultiModalNeighbors(
  object = seuratObj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = c("SCT.weight", "ATAC.weight"),
  verbose = TRUE
)
```

### Individual Modality Neighbors

#### SCT/RNA Neighbors
```r
DefaultAssay(seuratObj) <- "SCT"

# Create neighbor object for return
seuratObj <- FindNeighbors(
  object = seuratObj,
  reduction = "pca",
  dims = 1:50,
  k.param = 20,
  compute.SNN = FALSE,
  graph.name = c("sct.nn"),
  verbose = TRUE,
  return.neighbor = TRUE
)

# Create SNN graph for clustering
seuratObj <- FindNeighbors(
  object = seuratObj,
  reduction = "pca",
  dims = 1:50,
  k.param = 20,
  compute.SNN = TRUE,
  graph.name = c("SCT_nn", "SCT_snn"),
  verbose = TRUE,
  return.neighbor = FALSE
)
```

#### ATAC Neighbors
```r
DefaultAssay(seuratObj) <- "ATAC"

# Create neighbor object for return
seuratObj <- FindNeighbors(
  object = seuratObj,
  reduction = "lsi",
  dims = 2:40,  # Skip first LSI component (depth-related)
  k.param = 20,
  compute.SNN = FALSE,
  graph.name = c("atac.nn"),
  verbose = TRUE,
  return.neighbor = TRUE
)

# Create SNN graph for clustering
seuratObj <- FindNeighbors(
  object = seuratObj,
  reduction = "lsi",
  dims = 2:40,
  k.param = 20,
  compute.SNN = TRUE,
  graph.name = c("ATAC_nn", "ATAC_snn"),
  verbose = TRUE,
  return.neighbor = FALSE
)
```

## Clustering Analysis

### RNA-only Clustering
```r
seuratObj <- FindClusters(
  seuratObj,
  graph.name = "SCT_snn",
  resolution = 0.5,
  cluster.name = "SCT_clusters"
)
```

### ATAC-only Clustering
```r
seuratObj <- FindClusters(
  seuratObj,
  graph.name = "ATAC_snn",
  resolution = 0.5,
  cluster.name = "ATAC_clusters"
)
```

### Joint RNA+ATAC Clustering
```r
seuratObj <- FindClusters(
  object = seuratObj,
  graph.name = "wsnn",        # Use weighted graph
  algorithm = 3,              # SLM algorithm
  resolution = 0.5,
  verbose = TRUE,
  cluster.name = "weighted_clusters"
)
```

## UMAP Visualization

```r
# Joint UMAP using weighted neighbors
seuratObj <- RunUMAP(
  object = seuratObj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  verbose = TRUE
)

# SCT-based UMAP
seuratObj <- RunUMAP(
  object = seuratObj,
  nn.name = "sct.nn",
  reduction.name = "sct.umap",
  reduction.key = "sctUMAP_",
  verbose = TRUE
)

# ATAC-based UMAP
seuratObj <- RunUMAP(
  object = seuratObj,
  nn.name = "atac.nn",
  reduction.name = "atac.umap",
  reduction.key = "atacUMAP_",
  verbose = TRUE
)
```

### Visualization of Clustering Results
```r
DefaultAssay(seuratObj) <- "SCT"

# Create comparison plots
p1 = DimPlot(seuratObj, label = TRUE, repel = TRUE,
             reduction = "wnn.umap",
             group.by = "weighted_clusters") + NoLegend()
p2 = DimPlot(seuratObj, label = TRUE, repel = TRUE,
             reduction = "sct.umap",
             group.by = "SCT_clusters") + NoLegend()
p3 = DimPlot(seuratObj, label = TRUE, repel = TRUE, 
             reduction = "atac.umap",
             group.by = "ATAC_clusters") + NoLegend()

# Save combined plot
graphics.off()
pdf(str_c(analysis_save, "UMAP_summary.pdf"), width = 20, height = 6)
plot_grid(p1, p2, p3, ncol = 3)
dev.off()
```

## Clustering Evaluation

### Pairwise Clustering Comparison
```r
# Define clustering methods to compare
clustering_methods <- c("weighted_clusters", 
                        "SCT_clusters",
                        "ATAC_clusters")

cluster_matrix <- seuratObj@meta.data[, clustering_methods]
cluster_list <- as.list(cluster_matrix)
```

### Custom Function for Pairwise Metrics
```r
pairwise_clustering_metrics <- function(clustering_list) {
  
  method_names <- names(clustering_list)
  n_methods <- length(clustering_list)
  
  # Initialize matrices
  ari_matrix <- matrix(1, nrow = n_methods, ncol = n_methods,
                       dimnames = list(method_names, method_names))
  nmi_matrix <- matrix(1, nrow = n_methods, ncol = n_methods,
                       dimnames = list(method_names, method_names))
  
  # Calculate pairwise metrics
  cat("Calculating pairwise metrics...\n")
  
  for (i in 1:(n_methods-1)) {
    for (j in (i+1):n_methods) {
      cat(sprintf("Comparing %s vs %s\n", method_names[i], method_names[j]))
      
      ari_val <- ARI(clustering_list[[i]], clustering_list[[j]])
      nmi_val <- NMI(clustering_list[[i]], clustering_list[[j]])
      
      ari_matrix[i, j] <- ari_val
      ari_matrix[j, i] <- ari_val
      nmi_matrix[i, j] <- nmi_val
      nmi_matrix[j, i] <- nmi_val
    }
  }
  
  # Calculate summary statistics
  ari_upper <- ari_matrix[upper.tri(ari_matrix)]
  nmi_upper <- nmi_matrix[upper.tri(nmi_matrix)]
  
  summary_stats <- data.frame(
    Metric = c("ARI", "NMI"),
    Mean = c(mean(ari_upper), mean(nmi_upper)),
    Median = c(median(ari_upper), median(nmi_upper)),
    Min = c(min(ari_upper), min(nmi_upper)),
    Max = c(max(ari_upper), max(nmi_upper)),
    SD = c(sd(ari_upper), sd(nmi_upper))
  )
  
  # Create pairwise comparison table
  comparison_table <- data.frame(
    Method1 = character(),
    Method2 = character(),
    ARI = numeric(),
    NMI = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(n_methods-1)) {
    for (j in (i+1):n_methods) {
      comparison_table <- rbind(comparison_table, data.frame(
        Method1 = method_names[i],
        Method2 = method_names[j],
        ARI = round(ari_matrix[i, j], 4),
        NMI = round(nmi_matrix[i, j], 4)
      ))
    }
  }
  
  return(list(
    ari_matrix = ari_matrix,
    nmi_matrix = nmi_matrix,
    summary_stats = summary_stats,
    pairwise_table = comparison_table
  ))
}
```

### Visualization Function
```r
plot_pairwise_metrics <- function(results) {
  
  par(mfrow = c(1, 2))
  
  # ARI Heatmap
  corrplot(results$ari_matrix, 
           method = "color", 
           type = "full",
           title = "Adjusted Rand Index",
           mar = c(0,0,2,0),
           tl.cex = 0.8,
           cl.cex = 0.8,
           addCoef.col = "black",
           number.cex = 0.7,
           col = colorRampPalette(c("blue", "white", "red"))(200))
  
  # NMI Heatmap
  corrplot(results$nmi_matrix, 
           method = "color", 
           type = "full",
           title = "Normalized Mutual Information",
           mar = c(0,0,2,0),
           tl.cex = 0.8,
           cl.cex = 0.8,
           addCoef.col = "black",
           number.cex = 0.7,
           col = colorRampPalette(c("blue", "white", "red"))(200))
  
  par(mfrow = c(1, 1))
}
```

### Execute Evaluation
```r
# Run pairwise comparison
pairwise_results <- pairwise_clustering_metrics(cluster_list)

# Save evaluation plots
graphics.off()
pdf(str_c(analysis_save, "pairwise_clustering_metrics.pdf"),
    height = 10,
    width = 20)
plot_pairwise_metrics(pairwise_results)
dev.off()
```

## Data Export

### Save R Object
```r
# Save the complete Seurat object
saveRDS(seuratObj, file = paste0(data_save, "seuratObj_with_clustering.rds"))
```

### Export to H5 Format for Python

#### Data Extraction
```r
# Create directory for Python export
dir.create(paste0(data_save, "for_python"), showWarnings = FALSE)
python_save_rnaseq = paste0(data_save, "for_python/rnaseq_normalized.h5")
python_save_atacseq = paste0(data_save, "for_python/atacseq_normalized.h5")

# Extract matrices and metadata
rna_counts <- GetAssayData(seuratObj, assay = "RNA", slot = "counts")
rna_data <- GetAssayData(seuratObj, assay = "SCT", slot = "data")
atac_counts <- GetAssayData(seuratObj, assay = "ATAC", slot = "counts")  
atac_data <- GetAssayData(seuratObj, assay = "ATAC", slot = "data")

# Combine metadata with embeddings
metadata <- seuratObj@meta.data
metadata <- cbind(metadata, Embeddings(seuratObj))
metadata <- cbind(metadata, Embeddings(seuratObj, reduction = "wnn.umap"))
metadata <- cbind(metadata, Embeddings(seuratObj, reduction = "sct.umap"))
metadata <- cbind(metadata, Embeddings(seuratObj, reduction = "atac.umap"))

# Convert factors to characters
for(col in colnames(metadata)) {
  if(is.factor(metadata[[col]])) {
    metadata[[col]] <- as.character(metadata[[col]])
  }
}

# Convert to dense matrices
rna_counts <- as.matrix(rna_counts)
rna_data <- as.matrix(rna_data)
atac_counts <- as.matrix(atac_counts)
atac_data <- as.matrix(atac_data)
gc()  # Garbage collection
```

#### RNA-seq H5 Export
```r
h5file <- H5File$new(python_save_rnaseq, mode = "w")

h5file[["rna_counts"]] <- rna_counts
h5file[["rna_normalized"]] <- rna_data
h5file[["rna_gene_names"]] <- rownames(rna_counts)
h5file[["cell_names"]] <- colnames(rna_counts)

metadata_group <- h5file$create_group("metadata")
for (col in colnames(metadata)) {
  metadata_group[[col]] <- metadata[[col]]
}

h5file$close_all()
```

#### ATAC-seq H5 Export
```r
h5file <- H5File$new(python_save_atacseq, mode = "w")

h5file[["atac_counts"]] <- atac_counts
h5file[["atac_normalized"]] <- atac_data
h5file[["atac_feature_names"]] <- rownames(atac_counts)
h5file[["cell_names"]] <- colnames(atac_counts)

metadata_group <- h5file$create_group("metadata")
for (col in colnames(metadata)) {
  metadata_group[[col]] <- metadata[[col]]
}

h5file$close_all()

cat("Data saved to rnaseq/atacseq_normalized.h5 files\n")
cat("Matrices saved: rna_counts, rna_normalized, atac_counts, atac_normalized\n")
cat("Metadata saved in metadata/ group\n")
```

## Key Features

- **Multi-modal Integration**: Combines RNA-seq and ATAC-seq data using weighted nearest neighbors
- **Multiple Clustering Approaches**: Compares RNA-only, ATAC-only, and joint clustering
- **Comprehensive Evaluation**: Uses ARI and NMI metrics to assess clustering consistency
- **Flexible Export**: Saves data in both R (RDS) and Python-compatible (H5) formats
- **Visualization**: Creates UMAP plots and correlation heatmaps for comparison

## Output Files

- `seuratObj_with_clustering.rds`: Complete R object with all analyses
- `UMAP_summary.pdf`: Side-by-side UMAP comparisons
- `pairwise_clustering_metrics.pdf`: Clustering evaluation heatmaps
- `rnaseq_normalized.h5`: RNA data for Python analysis
- `atacseq_normalized.h5`: ATAC data for Python analysis