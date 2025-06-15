# load required packages ---------
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(hdf5r)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ConsensusClusterPlus)
library(aricode)
library(corrplot)
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

# load data ---------
setwd("./script/")
data_read = "../data/"
analysis_save = "../analysis/"
data_save = "../data/"

seuratObj = readRDS(paste0(data_read, "seuratObj_filtered.rds"))

# normalization ----------
DefaultAssay(seuratObj) <- "RNA"
seuratObj <- SCTransform(seuratObj)
seuratObj <- RunPCA(seuratObj)

DefaultAssay(seuratObj) <- "ATAC"
seuratObj <- FindTopFeatures(seuratObj, min.cutoff = 5)
seuratObj <- RunTFIDF(seuratObj)
seuratObj <- RunSVD(seuratObj)

# weight calculations ----------
## joint weight calculations  ---------
# (output wknn, wsnn and neighbor object (weight.nn))
DefaultAssay(seuratObj) <- "SCT"
seuratObj <- FindMultiModalNeighbors(
  object = seuratObj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = c("SCT.weight", "ATAC.weight"),
  verbose = TRUE
)

## invididual weight calculations  ---------
### use return.neighbor = TRUE to get the neighbor object but you should set
### compute.SNN = FALSE to avoid computing SNN graph at the same time
### please pay attention to nomenclature for sct.nn/atac.nn 
### (You cannot use the same name as you use in graph object)

### SCT -----------
DefaultAssay(seuratObj) <- "SCT"  # or "RNA"

seuratObj <- FindNeighbors(
  object = seuratObj,
  reduction = "pca",           # Use PCA reduction from RNA data
  dims = 1:50,                 # Use first 50 PCs
  k.param = 20,                # Number of nearest neighbors (default: 20)
  compute.SNN = FALSE,          # Compute shared nearest neighbor graph
  graph.name = c("sct.nn"),  # Names for NN and SNN graphs
  verbose = TRUE,
  return.neighbor = TRUE
)

seuratObj <- FindNeighbors(
  object = seuratObj,
  reduction = "pca",           # Use PCA reduction from RNA data
  dims = 1:50,                 # Use first 50 PCs
  k.param = 20,                # Number of nearest neighbors (default: 20)
  compute.SNN = TRUE,          # Compute shared nearest neighbor graph
  graph.name = c("SCT_nn", "SCT_snn"),  # Names for NN and SNN graphs
  verbose = TRUE,
  return.neighbor = FALSE
)

### ATAC -----------
DefaultAssay(seuratObj) <- "ATAC"

seuratObj <- FindNeighbors(
  object = seuratObj,
  reduction = "lsi",           # Use LSI reduction from ATAC data
  dims = 2:40,                 # Skip first LSI component (depth-related)
  k.param = 20,                # Number of nearest neighbors
  compute.SNN = FALSE,          # Compute shared nearest neighbor graph
  graph.name = c("atac.nn"),  # Names for graphs
  verbose = TRUE,
  return.neighbor = TRUE
)

seuratObj <- FindNeighbors(
  object = seuratObj,
  reduction = "lsi",           # Use LSI reduction from ATAC data
  dims = 2:40,                 # Skip first LSI component (depth-related)
  k.param = 20,                # Number of nearest neighbors
  compute.SNN = TRUE,          # Compute shared nearest neighbor graph
  graph.name = c("ATAC_nn", "ATAC_snn"),  # Names for graphs
  verbose = TRUE,
  return.neighbor = FALSE
)

# clustering ------------
## based on RNA data only ------
seuratObj <- FindClusters(
  seuratObj,
  graph.name = "SCT_snn",
  resolution = 0.5,
  cluster.name = "SCT_clusters"
)

## based on ATAC data only ------
seuratObj <- FindClusters(
  seuratObj,
  graph.name = "ATAC_snn",
  resolution = 0.5,
  cluster.name = "ATAC_clusters"
)

## based on joint RNA and ATAC data -------
seuratObj <- FindClusters(
  object = seuratObj,
  graph.name = "wsnn",        # Use weighted graph
  algorithm = 3,              # SLM algorithm (recommended for large datasets)
  resolution = 0.5,           # Adjust based on desired granularity
  verbose = TRUE,
  cluster.name = "weighted_clusters"
)

## build a joint UMAP visualization ----------
# Using weighted k-nearest neighbors (most common for multimodal)
seuratObj <- RunUMAP(
  object = seuratObj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  verbose = TRUE
)

# SCT-based UMAPs (SCT-transformed data)
seuratObj <- RunUMAP(
  object = seuratObj,
  nn.name = "sct.nn",
  reduction.name = "sct.umap",
  reduction.key = "sctUMAP_",
  verbose = TRUE
)

# ATAC-based UMAPs
seuratObj <- RunUMAP(
  object = seuratObj,
  nn.name = "atac.nn",
  reduction.name = "atac.umap",
  reduction.key = "atacUMAP_",
  verbose = TRUE
)

## visualization of clustering results ---------
DefaultAssay(seuratObj) <- "SCT"
p1 = DimPlot(seuratObj, label = TRUE, repel = TRUE,
             reduction = "wnn.umap",
             group.by = "weighted_clusters") + NoLegend()
p2 = DimPlot(seuratObj, label = TRUE, repel = TRUE,
             reduction = "sct.umap",
             group.by = "SCT_clusters") + NoLegend()
p3 = DimPlot(seuratObj, label = TRUE, repel = TRUE, 
             reduction = "atac.umap",
             group.by = "ATAC_clusters") + NoLegend()

graphics.off()
pdf(str_c(analysis_save, "UMAP_summary.pdf"), width = 20, height = 6)
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

# clustering evaluation based on different clustering usage ---------
# Create consensus clustering from multiple methods
clustering_methods <- c("weighted_clusters", 
                        "SCT_clusters",
                        "ATAC_clusters")

# Calculate consensus matrix
cluster_matrix <- seuratObj@meta.data[, clustering_methods]
cluster_list <- as.list(cluster_matrix)

# Pairwise Clustering Comparison Metrics
# Compare consistency between multiple clustering methods using ARI and NMI

# Function to calculate pairwise clustering metrics
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

# Function to visualize pairwise metrics
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

pairwise_results <- pairwise_clustering_metrics(cluster_list)

graphics.off()
pdf(str_c(analysis_save, "pairwise_clustering_metrics.pdf"),
    height = 10,
    width = 20)
plot_pairwise_metrics(pairwise_results)
dev.off()

# Save results ----------
saveRDS(seuratObj, file = paste0(data_save, "seuratObj_with_clustering.rds"))

# Save h5ad ---------
dir.create(paste0(data_save, "for_python"), showWarnings = FALSE)
python_save_rnaseq = paste0(data_save, "for_python/rnaseq_normalized.h5")
python_save_atacseq = paste0(data_save, "for_python/atacseq_normalized.h5")

## Extract the data --------
rna_counts <- GetAssayData(seuratObj, assay = "RNA", slot = "counts")
rna_data <- GetAssayData(seuratObj, assay = "SCT", slot = "data")
atac_counts <- GetAssayData(seuratObj, assay = "ATAC", slot = "counts")  
atac_data <- GetAssayData(seuratObj, assay = "ATAC", slot = "data")
metadata <- seuratObj@meta.data
metadata <- cbind(metadata, Embeddings(seuratObj))
metadata <- cbind(metadata, Embeddings(seuratObj, reduction = "wnn.umap"))
metadata <- cbind(metadata, Embeddings(seuratObj, reduction = "sct.umap"))
metadata <- cbind(metadata, Embeddings(seuratObj, reduction = "atac.umap"))
for(col in colnames(metadata)) {
  if(is.factor(metadata[[col]])) {
    metadata[[col]] <- as.character(metadata[[col]])
  }
}

rna_counts <- as.matrix(rna_counts)
rna_data <- as.matrix(rna_data)
atac_counts <- as.matrix(atac_counts)
atac_data <- as.matrix(atac_data)
gc()
gc()

## h5ad file saving --------
### rnaseq -------
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

### atacseq -------
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
