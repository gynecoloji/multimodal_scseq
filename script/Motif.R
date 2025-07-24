library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ConsensusClusterPlus)
library(aricode)
library(corrplot)
library(data.table)
library(pryr)
library(Matrix)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
library(foreach)
library(doParallel)
library(Rsamtools)
library(AnnotationHub)
library(stringr)
library(vroom)
library(GenomeInfoDb)
library(ggplot2)
library(Gviz)
library(InteractionSet)
library(patchwork)
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)

library(glmnet)          # For LASSO/Ridge regression
library(randomForest)    # For random forest
library(caret)          # For model training and validation
library(GenomicRanges)
library(GenomicInteractions)# For genomic data handling
library(Matrix)         # For sparse matrices
library(parallel)       # For parallel processing

setwd("./script")

# load data and dataset --------
data_read <- "../data/"
analysis_save <- "../analysis/"
data_save <- "../data/"

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))
## exclusion of gap regions is optional and if so run below one line code
annotation <- annotation[!(annotation$type %in% "gap")]

seuratObj <- readRDS(file = paste0(data_read, "seuratObj_with_clustering.rds")) # nolint: object_name_linter.

## motif analysis using Signac ----------
# load pfm into seurat object
DefaultAssay(seuratObj) = "ATAC" # nolint: object_name_linter.
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates",
              species = 9606,
              all_versions = FALSE)
)

seuratObj <- AddMotifs(
  object = seuratObj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# Find differentially accessible peaks between cell types

da_peaks <- FindMarkers(
  object = seuratObj,
  ident.1 = 0,
  ident.2 = 1,
  only.pos = TRUE,
  test.use = "LR",
  latent.vars = "nCount_ATAC"
)

# Get top differential peaks
top_da_peaks <- rownames(da_peaks[da_peaks$p_val_adj < 0.05, ])

# Find enriched motifs in differential peaks
enriched_motifs <- FindMotifs(
  object = seuratObj,
  features = top_da_peaks
)

# plot motif plots
MotifPlot(
  object = seuratObj,
  motifs = head(rownames(enriched_motifs))
)


# motif activity analysis --------------
seuratObj <- RunChromVAR(
  object = seuratObj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# differential motif analysis -------------
DefaultAssay(seuratObj) <- 'chromvar'
differential.activity <- FindMarkers(
  object = seuratObj,
  ident.1 = 0,
  ident.2 = 1,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = seuratObj,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)