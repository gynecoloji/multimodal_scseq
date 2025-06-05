# Joint scATAC-seq&scRNAseq analysis with Seurat and Signac
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ConsensusClusterPlus)
library(aricode)
library(corrplot)

data_read = "../data/"
analysis_save = "../analysis/"
data_save = "../data/"

# load data ---------
counts <- Read10X_h5(str_c(data_read,"pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"))
fragpath <- str_c(data_read,"pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")

# get gene annotations for hg38 -------
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA and atac data ---------
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(pbmc) <- "ATAC"

peaks.keep <- seqnames(granges(pbmc)) %in% str_c("chr", c(1:22, "X"))
pbmc[["ATAC"]] <- subset(pbmc[["ATAC"]],
                         features = rownames(pbmc[["ATAC"]])[as.vector(peaks.keep)])


pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)
DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)


# filtering using data from both atac and rna datasets-------------
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

# gene expression & atac data analysis (normalization) ----------
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

DefaultAssay(pbmc) <- "ATAC"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)

# weight calculations based on joint atac and rna seq----------
## joint weight calculations (output wknn, wsnn and neighbor object (weight.nn)) ---------
DefaultAssay(pbmc) <- "SCT"
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = c("SCT.weight", "ATAC.weight"),
  verbose = TRUE
)


## weight calculations based on atac and rna seq respectively  ---------
### use return.neighbor = TRUE to get the neighbor object but you should set
### compute.SNN = FALSE to avoid computing SNN graph at the same time
### please pay attention to nomenclature for sct.nn/atac.nn (You cannot use the same name as you use in graph object)

### SCT -----------

DefaultAssay(pbmc) <- "SCT"  # or "RNA"


pbmc <- FindNeighbors(
  object = pbmc,
  reduction = "pca",           # Use PCA reduction from RNA data
  dims = 1:50,                 # Use first 50 PCs
  k.param = 20,                # Number of nearest neighbors (default: 20)
  compute.SNN = FALSE,          # Compute shared nearest neighbor graph
  graph.name = c("sct.nn"),  # Names for NN and SNN graphs
  verbose = TRUE,
  return.neighbor = TRUE
)

pbmc <- FindNeighbors(
  object = pbmc,
  reduction = "pca",           # Use PCA reduction from RNA data
  dims = 1:50,                 # Use first 50 PCs
  k.param = 20,                # Number of nearest neighbors (default: 20)
  compute.SNN = TRUE,          # Compute shared nearest neighbor graph
  graph.name = c("SCT_nn", "SCT_snn"),  # Names for NN and SNN graphs
  verbose = TRUE,
  return.neighbor = FALSE
)

### ATAC -----------
DefaultAssay(pbmc) <- "ATAC"

pbmc <- FindNeighbors(
  object = pbmc,
  reduction = "lsi",           # Use LSI reduction from ATAC data
  dims = 2:40,                 # Skip first LSI component (depth-related)
  k.param = 20,                # Number of nearest neighbors
  compute.SNN = FALSE,          # Compute shared nearest neighbor graph
  graph.name = c("atac.nn"),  # Names for graphs
  verbose = TRUE,
  return.neighbor = TRUE
)

pbmc <- FindNeighbors(
  object = pbmc,
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
pbmc <- FindClusters(
  pbmc,
  graph.name = "SCT_snn",
  resolution = 0.5,
  cluster.name = "SCT_clusters"
)


## based on ATAC data only ------
pbmc <- FindClusters(
  pbmc,
  graph.name = "ATAC_snn",
  resolution = 0.5,
  cluster.name = "ATAC_clusters"
)

## based on joint RNA and ATAC data -------
pbmc <- FindClusters(
  object = pbmc,
  graph.name = "wsnn",        # Use weighted graph
  algorithm = 3,              # SLM algorithm (recommended for large datasets)
  resolution = 0.5,           # Adjust based on desired granularity
  verbose = TRUE,
  cluster.name = "weighted_clusters"
)


## build a joint UMAP visualization ----------
# Using weighted k-nearest neighbors (most common for multimodal)
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_",
  verbose = TRUE
)


# SCT-based UMAPs (SCT-transformed data)
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "sct.nn",
  reduction.name = "sct.umap",
  reduction.key = "sctUMAP_",
  verbose = TRUE
)

# ATAC-based UMAPs
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "atac.nn",
  reduction.name = "atac.umap",
  reduction.key = "atacUMAP_",
  verbose = TRUE
)




# visualization of clustering results ---------
DefaultAssay(pbmc) <- "SCT"
DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "wnn.umap") + NoLegend()
DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "sct.umap") + NoLegend()
DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "atac.umap") + NoLegend()


# clustering evaluation based on different clustering usage ---------
# Create consensus clustering from multiple methods
clustering_methods <- c("weighted_clusters", 
                        "SCT_clusters",
                        "ATAC_clusters")

# Calculate consensus matrix
cluster_matrix <- pbmc@meta.data[, clustering_methods]
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
plot_pairwise_metrics(pairwise_results)


# linking peaks to genes -----------
DefaultAssay(pbmc) <- "ATAC"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = c("LYZ", "MS4A1")
)

Idents(pbmc) <- "weighted_clusters"
idents.plot <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)

pbmc <- SortIdents(pbmc)

p1 <- CoveragePlot(
  object = pbmc,
  region = "MS4A1",
  features = "MS4A1",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)

p2 <- CoveragePlot(
  object = pbmc,
  region = "LYZ",
  features = "LYZ",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)

# save the object
saveRDS(pbmc, file = str_c(data_save, "pbmc_multiomic_seurat.rds"))







