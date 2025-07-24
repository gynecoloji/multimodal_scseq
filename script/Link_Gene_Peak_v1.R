# load required packages ---------
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
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
library(Gviz)
library(InteractionSet)
library(patchwork)

library(glmnet)          # For LASSO/Ridge regression
library(randomForest)    # For random forest
library(caret)          # For model training and validation
library(GenomicRanges)
library(GenomicInteractions)# For genomic data handling
library(Matrix)         # For sparse matrices
library(parallel)       # For parallel processing


# load data and dataset --------
data_read = "../data/"
analysis_save = "../analysis/"
data_save = "../data/"


annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
## exclusion of gap regions is optional and if so run below one line code
## annotation <- annotation[!(annotation$type %in% "gap")]

seuratObj <- readRDS(file = paste0(data_read, "seuratObj_with_clustering.rds"))




## for multiple genes -----------
perform_peak_gene_regression <- function(seuratObj,
                                         gene_list,
                                         annotation,
                                         window_size = 3000,  # 500kb window
                                         min_peaks = 2) {
  regression_results <- list()
  gr_peak = granges(seuratObj[['ATAC']])
  gr_gene = annotation[annotation$gene_name %in% rownames(seuratObj[['SCT']]), ]
  
  atac_matrix = seuratObj[['ATAC']]@data
  rna_matrix = seuratObj[['SCT']]@data
  
  
  for (gene in gene_list) {
    gr_select_gene = gr_gene[gr_gene$gene_name == gene, ]
    overlaps = overlapsAny(gr_peak, gr_select_gene, maxgap = window_size)
    
    
    if (sum(overlaps) < min_peaks) {
      message(paste("# of peaks found within", window_size, "bp of gene is less than specified values", gene, ". Skipping."))
      next
    }
    
    gr_select_peak_matrix = atac_matrix[overlaps, ]
    gr_select_gene_matrix = rna_matrix[rownames(seuratObj[['SCT']]) == gene, ]
    peak_names = rownames(gr_select_peak_matrix)
    

    y = as.numeric(gr_select_gene_matrix)
    X = t(as.matrix(gr_select_peak_matrix))
    df <- data.frame(y = y, X, check.names = FALSE)
    lm_fit <- lm(y ~ ., data = df)
    
    
    regression_results[[gene]] <- list(
      gene = gene,
      n_peaks = ncol(X),
      peak_names = peak_names,
      r_squared = summary(lm_fit)$r.squared,
      adj_r_squared = summary(lm_fit)$adj.r.squared,
      model = lm_fit
    )
    
  }
  
  return(regression_results)
}

linear_results = perform_peak_gene_regression(
  seuratObj = seuratObj,
  gene_list = VariableFeatures(seuratObj[['SCT']]),
  annotation = annotation,
  window_size = 3000,  # 3kb window
  min_peaks = 2
)

saveRDS(linear_results, file = paste0(data_save, "LPG_linear_results.rds"))

interpret_peak_gene_results <- function(regression_results, 
                                        r_squared_threshold = 0.3,
                                        p_value_threshold = 0.05) {
  
  summary_df <- data.frame(
    gene = character(),
    n_peaks = numeric(),
    r_squared = numeric(),
    adj_r_squared = numeric(),
    f_p_value = numeric(),
    model_significant = logical(),
    strong_relationship = logical(),
    n_significant_peaks = numeric(),
    top_peak = character(),
    top_peak_effect = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (gene_name in names(regression_results)) {
    result <- regression_results[[gene_name]]
    fstatistic <- summary(result$model)$fstatistic
    p <- pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=F)
    df_coef <- as.data.frame(summary(result$model)$coefficients)[-1,]# Exclude intercept
    
    # Overall model assessment
    model_significant <- p < p_value_threshold
    strong_relationship <- result$adj_r_squared > r_squared_threshold
    
    # Count significant peaks
    sig_peaks <- sum(df_coef[,4] < p_value_threshold, na.rm = TRUE)  
    
    # Identify most important peaks
    peak_effects <- abs(df_coef[,1]) 
    top_peak_idx <- which.max(peak_effects)
    top_peak <- result$peak_names[top_peak_idx]
    top_peak_effect <- df_coef[,1][top_peak_idx]
    
    summary_df <- rbind(summary_df, data.frame(
      gene = gene_name,
      n_peaks = result$n_peaks,
      r_squared = result$r_squared,
      adj_r_squared = result$adj_r_squared,
      f_p_value = unname(p),
      model_significant = model_significant,
      strong_relationship = strong_relationship,
      n_significant_peaks = sig_peaks,
      top_peak = top_peak,
      top_peak_effect = top_peak_effect
    ))
  }
  
  return(summary_df)
}

summary_results <- interpret_peak_gene_results(linear_results, 
                                                r_squared_threshold = 0.3, 
                                                p_value_threshold = 0.05)



# regression visualization analysis for single gene --------

## select sig. peak-gene-links ------
summary_results_select <- summary_results %>%
  filter(strong_relationship == TRUE) %>%
  arrange(desc(adj_r_squared))
gene <- summary_results_select$gene[1]

summary_results_select <- summary_results_select[1,]

## prepare granges ------
gr_peaks = granges(seuratObj[['ATAC']])
names(gr_peaks) <- rownames(seuratObj[['ATAC']])

gr_genes = annotation[annotation$gene_name %in% rownames(seuratObj[['SCT']]), ]
gr_select_gene = gr_genes[gr_genes$gene_name == gene, ]
gen <- 'hg38'
chr <- as.character(unique(seqnames(gr_select_gene)))

result = linear_results[[gene]]

### gr_select_peak --------
gr_select_peak = gr_peaks[names(gr_peaks) %in% result$peak_names]


### gene model for exons ------------
geneModels = gr_select_gene
geneModels$exon <- names(geneModels)
geneModels <- as.data.frame(geneModels, row.names = NULL)
colnames(geneModels) <- c("chromosome", "start", "end", "width", "strand",
                          "transcript", "symbol", "gene", "feature", "type",
                          "exon")
geneModels = geneModels[,c("chromosome", "start", "end", "width", "strand", 
                           "feature", "gene", "exon", "transcript", "symbol",
                           "type")]
geneModels = geneModels[geneModels$type %in% c('exon'), ]
geneModels = geneModels[,-ncol(geneModels)]



### gene model for cds ------------
geneModels_cds = gr_select_gene
geneModels_cds$exon <- names(geneModels_cds)
geneModels_cds <- as.data.frame(geneModels_cds, row.names = NULL)
colnames(geneModels_cds) <- c("chromosome", "start", "end", "width", "strand",
                          "transcript", "symbol", "gene", "feature", "type",
                          "exon")
geneModels_cds = geneModels_cds[,c("chromosome", "start", "end", "width", "strand", 
                           "feature", "gene", "exon", "transcript", "symbol",
                           "type")]
geneModels_cds = geneModels_cds[geneModels_cds$type %in% c('cds'), ]
geneModels_cds = geneModels_cds[,-ncol(geneModels_cds)]


### gene model for gap ------------
geneModels_gap = gr_select_gene
geneModels_gap$exon <- names(geneModels_gap)
geneModels_gap <- as.data.frame(geneModels_gap, row.names = NULL)
colnames(geneModels_gap) <- c("chromosome", "start", "end", "width", "strand",
                              "transcript", "symbol", "gene", "feature", "type",
                              "exon")
geneModels_gap = geneModels_gap[,c("chromosome", "start", "end", "width", "strand", 
                                   "feature", "gene", "exon", "transcript", "symbol",
                                   "type")]
geneModels_gap = geneModels_gap[geneModels_gap$type %in% c('gap'), ]
geneModels_gap = geneModels_gap[,-ncol(geneModels_gap)]



## prepare quantitative data --------
atac_peak_matrix = seuratObj[['ATAC']]@data[names(gr_select_peak),]
atac_peak_matrix = as.data.frame(Matrix::rowMeans(atac_peak_matrix))
colnames(atac_peak_matrix) = "peak_score"
atac_peak_coef = as.data.frame(summary(result$model)$coefficients)[-1,]

identical(rownames(atac_peak_matrix), gsub("`","",rownames(atac_peak_coef)))
atac_peak_matrix$estimate = atac_peak_coef$Estimate
atac_peak_matrix$sig = atac_peak_coef[,4] < 0.05
atac_peak_matrix$chr = sapply(str_split(rownames(atac_peak_matrix),  "[:-]"), `[`, 1)
atac_peak_matrix$start = sapply(str_split(rownames(atac_peak_matrix),  "[:-]"), `[`, 2)
atac_peak_matrix$end = sapply(str_split(rownames(atac_peak_matrix),  "[:-]"), `[`, 3)






## prepare interaction data ---------
pg_interaction = findOverlapPairs(gr_select_peak, 
                                  gr_select_gene, 
                                  maxgap = 3000)
peak_interaction = pg_interaction@first
gene_interaction = pg_interaction@second

mcols(peak_interaction)$strength = atac_peak_matrix[match(names(peak_interaction), 
                                                      rownames(atac_peak_matrix)),][,'estimate']
mcols(peak_interaction)$strength = atac_peak_matrix[match(names(peak_interaction), 
                                                          rownames(atac_peak_matrix)),][,'estimate']


mcols(peak_interaction)$midpoint = (start(peak_interaction) + end(peak_interaction)) / 2
mcols(gene_interaction)$midpoint = (start(gene_interaction) + end(gene_interaction)) / 2


link_num = length(pg_interaction@second)


# Create anchor points for interactions
# Anchor 1: Peak positions
anchor1 <- GRanges(
  seqnames = rep(chr, link_num),
  ranges = IRanges(
    start = mcols(peak_interaction)$midpoint,
    end = mcols(peak_interaction)$midpoint + 1  # Make valid range
  )
)

# Anchor 2: Gene positions  
anchor2 <- GRanges(
  seqnames = rep(chr, link_num),
  ranges = IRanges(
    start = mcols(gene_interaction)$midpoint,
    end = mcols(gene_interaction)$midpoint + 1  # Make valid range
  )
)

# Create GenomicInteractions object
gi_object <- GenomicInteractions(
  anchor1 = anchor1,
  anchor2 = anchor2,
  counts = round(mcols(peak_interaction)$strength * 100)
)

# Add interaction metadata
gi_object$strength <- mcols(peak_interaction)$strength
gi_object$from_feature <- names(peak_interaction)
gi_object$to_feature <- rep(gene, link_num)
gi_object$interaction_type <- rep("peak_gene", link_num)


# selective interaction sets
gi_object_sig = gi_object[mcols(gi_object)$counts == max(mcols(gi_object)$counts)]



## prepare frag coverage ---------
fragment_file <- "../data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
tbx <- TabixFile(fragment_file)

calculate_bin_heights <- function(coverage_rle, regions_gr, bin_size = 100) {
  
  bin_results <- list()
  
  for(i in seq_along(regions_gr)) {
    region <- regions_gr[i]
    chr_name <- as.character(seqnames(region))
    region_start <- start(region)
    region_end <- end(region)
    
    # Extract coverage for this region
    if(chr_name %in% names(coverage_rle)) {
      region_coverage <- coverage_rle[[chr_name]][region_start:region_end]
      coverage_vector <- as.numeric(region_coverage)
    } else {
      coverage_vector <- rep(0, width(region))
    }
    
    # Calculate bins
    n_bins <- ceiling(length(coverage_vector) / bin_size)
    bin_heights <- numeric(n_bins)
    bin_starts <- numeric(n_bins)
    bin_ends <- numeric(n_bins)
    
    for(j in 1:n_bins) {
      start_idx <- (j-1) * bin_size + 1
      end_idx <- min(j * bin_size, length(coverage_vector))
      
      # Calculate bin height (sum, mean, or max)
      bin_heights[j] <- sum(coverage_vector[start_idx:end_idx])
      
      # Calculate genomic coordinates
      bin_starts[j] <- region_start + start_idx - 1
      bin_ends[j] <- region_start + end_idx - 1
    }
    
    bin_results[[i]] <- data.frame(
      chr = chr_name,
      bin_start = bin_starts,
      bin_end = bin_ends,
      bin_center = (bin_starts + bin_ends) / 2,
      height = bin_heights,
      region_id = paste0(chr_name, ":", region_start, "-", region_end),
      stringsAsFactors = FALSE
    )
  }
  
  return(do.call(rbind, bin_results))
}

parse_fragments <- function(fragment_lines) {
  if(length(fragment_lines) == 0) return(NULL)
  
  # Split lines and create data frame
  fragment_data <- do.call(rbind, strsplit(fragment_lines, "\t"))
  colnames(fragment_data) <- c("chr", "start", "end", "barcode", "count")
  
  # Convert to appropriate types
  fragment_df <- data.frame(
    chr = fragment_data[,1],
    start = as.numeric(fragment_data[,2]),
    end = as.numeric(fragment_data[,3]),
    barcode = fragment_data[,4],
    count = as.numeric(fragment_data[,5]),
    stringsAsFactors = FALSE
  )
  
  return(fragment_df)
}

# Define regions to extract
regions <- GRanges(
  seqnames = chr,
  ranges = IRanges(
    start = min(start(atrack)),
    end = max(end(atrack)))
  )

fragments_list <- scanTabix(tbx, param = regions)


# Parse extracted fragments
### for all data -------
extracted_data <- lapply(fragments_list, parse_fragments)
extracted_data <- extracted_data[[1]]
extracted_data <- extracted_data[extracted_data$barcode %in% colnames(seuratObj),]

fragments <- GRanges(
  seqnames = extracted_data$chr,
  ranges = IRanges(start = extracted_data$start,
                   end = extracted_data$end)
)
coverage_rle <- coverage(fragments)


bin_heights <- calculate_bin_heights(
  coverage_rle = coverage_rle,
  regions_gr = regions,
  bin_size = 50
)

### for grouped cells ---------
frag_summary <- CountFragments(fragments = fragment_file, 
                               cells = colnames(seuratObj))
identical(colnames(seuratObj), frag_summary$CB)

frag_summary$group = seuratObj$weighted_clusters
scaling_factor_ncell_ls = as.vector(table(frag_summary$group))
names(scaling_factor_ncell_ls) = levels(frag_summary$group)
scaling_factor_depth_ls = frag_summary[,c("frequency_count", "group")] %>%
  group_by(group) %>%
  summarise(total_count = sum(frequency_count), .groups = "drop")
scaling_factor_depth_ls = scaling_factor_depth_ls$total_count
names(scaling_factor_depth_ls) = levels(frag_summary$group)


bin_heights_by_group = function(extracted_data,
                                cell_barcodes,
                                regions,
                                fragments_list,
                                scaling_factor){
  
  
  extracted_data <- lapply(fragments_list, parse_fragments)
  extracted_data <- extracted_data[[1]]
  extracted_data <- extracted_data[extracted_data$barcode %in% cell_barcodes,]
  
  fragments <- GRanges(
    seqnames = extracted_data$chr,
    ranges = IRanges(start = extracted_data$start,
                     end = extracted_data$end)
  )
  coverage_rle <- coverage(fragments)
  
  
  bin_heights <- calculate_bin_heights(
    coverage_rle = coverage_rle,
    regions_gr = regions,
    bin_size = 50
  )
  
  bin_heights$height = bin_heights$height/scaling_factor*10000
  return(bin_heights)
}

bin_heights_depth = data.frame(chr = character(),
                               bin_start = numeric(),
                               bin_end = numeric(),
                               bin_center = numeric(),
                               height = numeric(),
                               region_id = character(),
                               group = factor())
for (i in c(0,1,2)) {
  cell_barcodes = WhichCells(seuratObj, idents = i)
  scaling_factor = scaling_factor_depth_ls[names(scaling_factor_depth_ls) %in% i]
  scaling_factor = unname(scaling_factor)
  tmp = bin_heights_by_group(extracted_data,
                             cell_barcodes,
                             regions,
                             fragments_list,
                             scaling_factor)
  tmp$group = i
  bin_heights_depth = rbind(bin_heights_depth, tmp)
}





## Combine together to plot --------
### plot single ftrack ------------
Idtrack <- IdeogramTrack(genome = gen, chromosome = chr)
gtrack <- GenomeAxisTrack()
ptrack <- AnnotationTrack(gr_select_peak, name = "peaks")
ftrack <- DataTrack(data = bin_heights$height,
                    start = bin_heights$bin_start,
                    end = bin_heights$bin_end, 
                    chromosome = bin_heights$chr, 
                    genome = gen, 
                    name = "frag\ncounts",
                    type = 'histogram',
                    col.histogram = 'pink')
grtrack <- GeneRegionTrack(geneModels, name = 'Exon', genome = gen, chromosome = chr)
cdstrack <- GeneRegionTrack(geneModels_cds, name = 'CDS', genome = gen, chromosome = chr)
gaptrack <- GeneRegionTrack(geneModels_gap, name = 'gap', genome = gen, chromosome = chr)
pstrack <- DataTrack(data = atac_peak_matrix$peak_score,
                    start = atac_peak_matrix$start,
                    end = atac_peak_matrix$end, 
                    chromosome = atac_peak_matrix$chr, 
                    genome = gen, 
                    name = "peak score", type = "smooth",
                    col = 'blue')
itrack <- InteractionTrack(gi_object, name = "Peak-Gene Interactions")
pitrack <- InteractionTrack(gi_object_sig, name = "Prominent\nPeak-Gene Interactions")

displayPars(pitrack) <- list(
  col.interactions = "blue",           # Interaction line color
  col.anchors.line = "black",         # Anchor border color  
  col.anchors.fill = "lightblue",     # Anchor fill color
  lwd = 10,               # Line width
  alpha.interactions = 0.8            # Transparency
)




plotTracks(list(Idtrack, gtrack, ptrack, ftrack, grtrack, cdstrack, pstrack, itrack, pitrack),
           sizes = c(1,2,2,2,2,2,3,4,4),
           from = min(start(ptrack))-5000,
           to = max(end(ptrack)) + 5000)


### plot multiple ftracks ------------
colors_group = c("green", "lightblue", "pink")
names(colors_group) = c(0, 1, 2)
ylim_range = c(0, max(bin_heights_depth$height))

for (i in c(0, 1, 2)) {
  tmp <- bin_heights_depth[bin_heights_depth$group %in% i,]
  tmp <- DataTrack(
    data = tmp$height,
    start = tmp$bin_start,
    end = tmp$bin_end,
    chromosome = tmp$chr,
    genome = gen,
    name = as.character(i),
    type = 'histogram',
    col.histogram = unname(colors_group[names(colors_group) %in% i]),
    ylim = ylim_range
  )
  
  assign(str_c("ftrack",i), tmp)
  
}

graphics.off()
pdf(file = str_c(analysis_save, "Gviz_tracks_", gene,".pdf"),
    height = 15,
    width = 10)
plotTracks(list(Idtrack, gtrack, ptrack, ftrack0, ftrack1, ftrack2, grtrack, cdstrack, pstrack, itrack, pitrack),
           sizes = c(1, # Idtrack
                     2, # gtrack
                     2, # ptrack for peak regions
                     2,2,2, # ftracks for normalized fragment counts
                     3, # grtrack for stadard gene models using exons
                     3, # grtrack for gene models using cds
                     4, # pstrack for peak scores for peak regions
                     3, # itrack for interactions
                     3), # pitrack for prominent interactions
           from = min(start(ptrack))-5000,
           to = max(end(ptrack)) + 5000)
dev.off()







