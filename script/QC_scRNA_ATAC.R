# Load Required Packages ----------
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
library(patchwork)
library(DoubletFinder)
library(scDblFinder)
library(patchwork)  # Make sure to load patchwork for combining plots

# QC for scATACseq ----------
## Setup Directories and Load Data ----------
setwd("./script/")
data_read <- "../data/"
analysis_save <- "../analysis/"
data_save <- "../data/"
counts <- Read10X_h5(str_c(data_read,
                           "filtered_feature_bc_matrix.h5"))
fragpath <- str_c(data_read,
                  "atac_fragments.tsv.gz")

### Separate Counts by Sample IDs ----------
# separate counts based on sample IDs (optional)
sample_name <- "Experimental"
counts$Peaks <- counts$Peaks[, grepl("-2", colnames(counts$Peaks))]
counts$`Gene Expression` <- counts$`Gene Expression`[, grepl("-2", colnames(counts$`Gene Expression`))] # nolint: line_length_linter.

## Prepare Repetitive Elements Annotation ----------
ah <- AnnotationHub()
query(ah, c("RepeatMasker", "Homo sapiens"))
repeats_hg38 <- ah[["AH99003"]]
repeats_hg38 <- repeats_hg38[!grepl("\\?", repeats_hg38$repClass), ]
repeats_hg38 <- repeats_hg38[!grepl("Unknown", repeats_hg38$repClass), ]

### Filter Chromosomes ----------
# because non-uniform nomenclature for chromosomes except for
# chr1-22, chrX, chrY and chrM, we just discard them.
repeats_hg38 <- repeats_hg38[seqnames(repeats_hg38) %in% c(paste0("chr", 1:22), "chrX", "chrY", "chrM"), ] 
repeats_hg38 <- keepSeqlevels(repeats_hg38, c(paste0("chr", 1:22), "chrX", "chrY", "chrM"), pruning.mode = "coarse") 

## Parse Peak Coordinates and Setup ----------
peak_names <- rownames(counts$Peaks)

### Define Coordinate Parser Function ----------
parse_coordinates <- function(coord_strings) {
  # Split by ":" and "-"
  parts <- str_split(coord_strings, "[:-]")

  # Extract components
  chr <- sapply(parts, `[`, 1)
  start <- as.numeric(sapply(parts, `[`, 2))
  end <- as.numeric(sapply(parts, `[`, 3))

  # Create GRanges
  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = start, end = end) # nolint: object_usage_linter.
  )

  return(gr) # nolint: return_linter.
}

### Create Peak GRanges Object ----------
peak_gr <- parse_coordinates(peak_names)
seqlevels(peak_gr)[seqlevels(peak_gr) == "chrMT"] <- "chrM"
peak_chrs <- sapply(strsplit(peak_names, ":"), `[`, 1)

## Peak Cell-Level QC ----------

### Standard Chromosome Analysis ----------
# Define chromosomes to keep for research
chromosomes_to_keep <- c(paste0("chr", 1:22), "chrX")

# Create logical vector for peaks to keep
peaks_keep <- peak_chrs %in% chromosomes_to_keep

peak_std_chr_ratio <- Matrix::colSums(keep_peaks <- counts$Peaks[peaks_keep, ]) / Matrix::colSums(counts$Peaks) # nolint: line_length_linter.
df_peak_qc_cell_lvl <- as.data.frame(peak_std_chr_ratio)
colnames(df_peak_qc_cell_lvl) <- "std_chr"

### Repetitive Elements Analysis ----------

#### Define Repetitive Element Analysis Functions ----------
calculate_repetitive_element_ratios_by_category = function(peak_gr, # nolint: object_length_linter.
                                                           repeats_hg38,
                                                           peak_matrix,
                                                           Category="all"){
  # Classify repeats by type
  high_confidence_removes <- c("Simple_repeat", "Low_complexity", "Satellite",
                               "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA")

  moderate_removes <- c("LINE", "SINE", "LTR")  # More controversial

  keep_types <- c("DNA", "RC")  # Often contain regulatory elements

  # select specific categories of repetitive elements
  repeats_hg38 <- switch(Category,
    'all' = repeats_hg38,
    'high_confidence_removes' = repeats_hg38[repeats_hg38$repClass %in% high_confidence_removes,],
    'moderate_removes' = repeats_hg38[repeats_hg38$repClass %in% moderate_removes,],
    'keep_types' = repeats_hg38[repeats_hg38$repClass %in% keep_types,],
    stop("Invalid category specified.")
  )
  
  # Create logical vector for peaks to keep
  peaks_keep <- overlapsAny(peak_gr, repeats_hg38)
  
  # Filter the peak matrix
  filtered_peaks <- peak_matrix[peaks_keep, ]
  
  peak_ratio = Matrix::colSums(filtered_peaks)/Matrix::colSums(peak_matrix)
  
  return(peak_ratio)
}

calculate_individual_repeat_type_ratios = function(peak_gr, 
                               repeats_hg38, 
                               peak_matrix){
  # all types
  categories <- c("Simple_repeat", "Low_complexity", "Satellite",
                  "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA",
                  "LINE", "SINE", "LTR", "DNA", "RC")

  df = data.frame(row.names = colnames(peak_matrix))
  for (i in categories) {
    # Select specific category of repetitive elements
    repeats_subset <- repeats_hg38[repeats_hg38$repClass == i, ]

    # Create logical vector for peaks to keep
    peaks_keep <- overlapsAny(peak_gr, repeats_subset)

    # Filter the peak matrix
    filtered_peaks <- peak_matrix[peaks_keep, ]

    # Calculate ratio
    peak_ratio = Matrix::colSums(filtered_peaks)/Matrix::colSums(peak_matrix)

    df[paste0(i,'_rep')] <- peak_ratio
  }
  
  return(df)
}

#### Calculate Repetitive Element Ratios by Category ----------
for (i in c('all', 'high_confidence_removes', 'moderate_removes', 'keep_types')) {
  peak_ratio = calculate_repetitive_element_ratios_by_category(
    peak_gr = peak_gr,
    repeats_hg38 = repeats_hg38,
    peak_matrix = counts$Peaks,
    Category = i
  )

  df_peak_qc_cell_lvl[paste0(i, '_rep')] <- unname(peak_ratio)
}

#### Calculate Individual Repeat Type Ratios ----------
tmp = calculate_individual_repeat_type_ratios(
  peak_gr = peak_gr,
  repeats_hg38 = repeats_hg38,
  peak_matrix = counts$Peaks
)

df_peak_qc_cell_lvl <- cbind(df_peak_qc_cell_lvl, tmp)

### Mitochondrial Chromosome Analysis ----------
peak_chrmt_ratio = Matrix::colSums(counts$Peaks[peak_chrs %in% c("chrM"), ]) / Matrix::colSums(counts$Peaks)
tmp = as.data.frame(peak_chrmt_ratio)
colnames(tmp) <- 'chrmt'
df_peak_qc_cell_lvl = cbind(df_peak_qc_cell_lvl, tmp)

### Blacklisted Regions Analysis ----------
peaks_keep <- overlapsAny(peak_gr, blacklist_hg38_unified)
filtered_peaks <- counts$Peaks[peaks_keep, ]
peak_bl_ratio = Matrix::colSums(filtered_peaks)/Matrix::colSums(counts$Peaks)
tmp = as.data.frame(peak_bl_ratio)
colnames(tmp) <- 'blacklist'
df_peak_qc_cell_lvl = cbind(df_peak_qc_cell_lvl,tmp)
### Save Results ----------
write.csv(
  df_peak_qc_cell_lvl,
  file = str_c(data_save, "cell_level_peak_qc.csv"),
  row.names = TRUE
)

## Peak Peak-Level QC ----------
df_peak_qc = data.frame(row.names = 'peak_ratio')
df_peak_qc_peak_lvl = data.frame(row.names = rownames(counts$Peaks))

### Standard Chromosome Analysis ----------
chromosomes_to_keep = c(paste0("chr", 1:22), "chrX")

peaks_keep <- seqnames(peak_gr) %in% chromosomes_to_keep
peak_std_chr_ratio = sum(peaks_keep)/length(peaks_keep)
df_peak_qc['std_chr'] <- peak_std_chr_ratio
df_peak_qc_peak_lvl['std_chr'] <- as.vector(peaks_keep)

### Repetitive Elements Analysis ----------

#### Define Repetitive Element Analysis Functions ----------
get_repeat_categories <- function() {
  list(
    high_confidence_removes = c("Simple_repeat", "Low_complexity", "Satellite", 
                                "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA"),
    moderate_removes = c("LINE", "SINE", "LTR"),
    keep_types = c("DNA", "RC"),
    all_individual = c("Simple_repeat", "Low_complexity", "Satellite", 
                       "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA",
                       "LINE", "SINE", "LTR", "DNA", "RC")
  )
}

# Filter repeats by category
filter_repeats_by_category <- function(repeats_hg38, category) {
  categories <- get_repeat_categories()
  
  if (category == "all") {
    return(repeats_hg38)
  }
  
  if (!category %in% names(categories)) {
    stop("Invalid category. Choose from: ", paste(names(categories), collapse = ", "))
  }
  
  return(repeats_hg38[repeats_hg38$repClass %in% categories[[category]], ])
}

# Calculate overlap ratios for grouped categories
calculate_category_ratios <- function(peak_gr, repeats_hg38) {
  categories <- c("all", "high_confidence_removes", "moderate_removes", "keep_types")
  
  ratios <- sapply(categories, function(cat) {
    filtered_repeats <- filter_repeats_by_category(repeats_hg38, cat)
    peaks_overlap <- overlapsAny(peak_gr, filtered_repeats)
    sum(peaks_overlap) / length(peaks_overlap)
  })

  # Convert to data frame with descriptive names
  df <- data.frame(t(ratios))
  colnames(df) <- paste0(categories, "_rep")
  rownames(df) <- "peak_ratio"

  return(df)
}

# Calculate overlap ratios for individual repeat types
calculate_individual_ratios <- function(peak_gr, repeats_hg38) {
  categories <- get_repeat_categories()$all_individual
  
  ratios <- sapply(categories, function(rep_type) {
    repeats_subset <- repeats_hg38[repeats_hg38$repClass == rep_type, ]
    peaks_overlap <- overlapsAny(peak_gr, repeats_subset)
    sum(peaks_overlap) / length(peaks_overlap)
  })
  
  # Convert to data frame with descriptive names
  df <- data.frame(t(ratios))
  colnames(df) <- paste0(categories, "_rep")
  rownames(df) <- "peak_ratio"
  
  return(df)
}

# Add binary overlap indicators for grouped categories
add_category_indicators <- function(peak_gr, repeats_hg38, df_peaks) {
  categories <- c("all", "high_confidence_removes", "moderate_removes", "keep_types")
  
  for (cat in categories) {
    filtered_repeats <- filter_repeats_by_category(repeats_hg38, cat)
    peaks_overlap <- overlapsAny(peak_gr, filtered_repeats)
    df_peaks[[paste0(cat, "_rep")]] <- as.vector(peaks_overlap)
  }
  
  return(df_peaks)
}

# Add binary overlap indicators for individual repeat types
add_individual_indicators <- function(peak_gr, repeats_hg38, df_peaks) {
  categories <- get_repeat_categories()$all_individual
  
  for (rep_type in categories) {
    repeats_subset <- repeats_hg38[repeats_hg38$repClass == rep_type, ]
    peaks_overlap <- overlapsAny(peak_gr, repeats_subset)
    df_peaks[[paste0(rep_type, "_rep")]] <- as.vector(peaks_overlap)
  }
  
  return(df_peaks)
}

# Main analysis function that combines everything
analyze_repetitive_elements <- function(peak_gr, repeats_hg38, peak_matrix = NULL, 
                                        df_peak_annotations = NULL, 
                                        include_ratios = TRUE, 
                                        include_indicators = TRUE) {
  
  results <- list()
  
  # Calculate summary ratios if requested
  if (include_ratios) {
    cat("Calculating category ratios...\n")
    category_ratios <- calculate_category_ratios(peak_gr, repeats_hg38)
    
    cat("Calculating individual repeat type ratios...\n")
    individual_ratios <- calculate_individual_ratios(peak_gr, repeats_hg38)
    
    # Combine ratio results
    results$summary_ratios <- cbind(category_ratios, individual_ratios)
  }
  
  # Add binary indicators if requested and data frame provided
  if (include_indicators && !is.null(df_peak_annotations)) {
    cat("Adding category indicators to peak annotations...\n")
    df_with_categories <- add_category_indicators(peak_gr, repeats_hg38, df_peak_annotations)
    
    cat("Adding individual repeat type indicators...\n")
    df_with_all <- add_individual_indicators(peak_gr, repeats_hg38, df_with_categories)
    
    results$annotated_peaks <- df_with_all
  }
  
  # Add summary statistics
  results$summary_stats <- list(
    total_peaks = length(peak_gr),
    total_repeats = length(repeats_hg38),
    repeat_classes = table(repeats_hg38$repClass)
  )
  
  return(results)
}

#### Run Repetitive Elements Analysis ----------
results <- analyze_repetitive_elements(
  peak_gr = peak_gr,
  repeats_hg38 = repeats_hg38,
  df_peak_annotations = df_peak_qc_peak_lvl
)

df_peak_qc_peak_lvl = results$annotated_peaks
df_peak_qc = cbind(df_peak_qc,results$summary_ratios)

### Mitochondrial Chromosome Analysis ----------
peaks_keep <- seqnames(peak_gr) %in% c("chrM")
peak_chrmt_ratio = sum(peaks_keep)/length(peaks_keep)
df_peak_qc['chrmt'] <- peak_chrmt_ratio
df_peak_qc_peak_lvl['chrmt'] <- as.vector(peaks_keep)

### Blacklisted Regions Analysis ----------
peaks_keep <- overlapsAny(peak_gr, blacklist_hg38_unified)
peak_bl_ratio = sum(peaks_keep)/length(peaks_keep)
df_peak_qc['blacklist'] <- peak_bl_ratio
df_peak_qc_peak_lvl['blacklist'] <- as.vector(peaks_keep)
df_peak_qc <- as.data.frame(t(df_peak_qc))

### Peak Expression Statistics ----------
df_peak_qc_peak_lvl['peak_cell_ratio'] <- unname(Matrix::rowSums(counts$Peaks > 0)/ncol(counts$Peaks))
df_peak_qc_peak_lvl['peak_cell'] <- unname(Matrix::rowSums(counts$Peaks > 0))
df_peak_qc_peak_lvl['peak_reads'] <- unname(Matrix::rowSums(counts$Peaks))
df_peak_qc_peak_lvl['peak_width'] <- width(peak_gr)

### Save Results ----------
write.csv(
  df_peak_qc,
  file = str_c(data_save, "peak_level_peak_qc_summary.csv"),
  row.names = TRUE
)

write.csv(
  df_peak_qc_peak_lvl,
  file = str_c(data_save, "peak_level_peak_qc.csv"),
  row.names = TRUE
)
## Setup Gene Annotations and Seurat Object ----------

### Get Gene Annotations for hg38 ----------
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

### Create Seurat Object with RNA and ATAC Data ----------
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

### Clean Up Memory ----------
rm(list = c("counts", "df_peak_qc", "df_peak_qc_cell_lvl", "df_peak_qc_peak_lvl", "keep_peaks", "results", "tmp", "filtered_peaks"))
gc()
gc()

## Fragment Data Processing ----------

### Get Fragment Information ----------
DefaultAssay(pbmc) <- "ATAC"
frag_summary <- CountFragments(fragments = fragpath, 
                               cells = colnames(pbmc))

cell_barcodes <- colnames(pbmc)

### Define Fragment Processing Functions ----------

#### Main GRanges Creation Function ----------
create_granges_from_tsv <- function(file_path,
                                    cell_barcodes,
                                    chunk_size = NULL,
                                    verbose = TRUE) {
  
  if (verbose) {
    cat("Creating GRanges from:", basename(file_path), "\n")
    cat("Target cell barcodes:", length(cell_barcodes), "\n")
  }
  
  # Determine if we need chunked reading
  if (is.null(chunk_size)) {
    # Try to estimate if chunked reading is needed
    file_size_mb <- file.size(file_path) / (1024^2)
    
    if (file_size_mb > 1000) {  # > 1GB, use chunked reading
      chunk_size <- 100000
      if (verbose) cat("Large file detected, using chunked reading\n")
    }
  }
  
  if (!is.null(chunk_size)) {
    # Use chunked reading for large files
    gr <- create_granges_chunked(file_path, cell_barcodes, chunk_size, verbose)
  } else {
    # Direct reading for smaller files
    gr <- create_granges_direct(file_path, cell_barcodes, verbose)
  }
  
  return(gr)
}

#### Direct File Reading Function ----------
create_granges_direct <- function(file_path, cell_barcodes, verbose = TRUE) {
  
  if (verbose) cat("Reading file directly...\n")
  
  # Read only first 4 columns
  data <- vroom(
    file_path,
    delim = "\t",
    col_names = FALSE,
    col_select = 1:4,
    col_types = cols(
      X1 = col_character(),  # chromosome
      X2 = col_integer(),    # start
      X3 = col_integer(),    # end  
      X4 = col_character()   # cell_barcode
    ),
    altrep = FALSE,
    show_col_types = FALSE,
    progress = verbose
  )
  
  # Filter by cell barcodes
  filtered_data <- data[data$X4 %in% cell_barcodes, ]
  
  # Convert to GRanges
  if (nrow(filtered_data) == 0) {
    warning("No rows match the provided cell barcodes!")
    return(GRanges())
  }
  
  gr <- GRanges(
    seqnames = filtered_data$X1,
    ranges = IRanges(start = filtered_data$X2, end = filtered_data$X3),
    cell_barcode = filtered_data$X4
  )
  
  if (verbose) {
    cat("Created GRanges with", length(gr), "ranges\n")
    cat("Unique chromosomes:", length(unique(seqnames(gr))), "\n")
    cat("Unique cell barcodes:", length(unique(gr$cell_barcode)), "\n")
  }
  
  return(gr)
}

#### Chunked File Reading Function ----------
create_granges_chunked <- function(file_path, cell_barcodes, chunk_size, verbose = TRUE) {
  
  if (verbose) cat("Using chunked reading with chunk size:", chunk_size, "\n")
  
  # Use the helper function to read in chunks
  filtered_chunks <- read_file_in_chunks(file_path, cell_barcodes, chunk_size, verbose)
  
  if (verbose) {
    cat("Chunked reading complete\n")
  }
  
  # Combine filtered chunks
  if (length(filtered_chunks) == 0) {
    warning("No rows match the provided cell barcodes!")
    return(GRanges())
  }
  
  combined_data <- bind_rows(filtered_chunks)
  
  # Convert to GRanges
  gr <- GRanges(
    seqnames = combined_data$X1,
    ranges = IRanges(start = combined_data$X2, end = combined_data$X3),
    cell_barcode = combined_data$X4
  )
  
  if (verbose) {
    cat("Created GRanges with", length(gr), "ranges\n")
    cat("Unique chromosomes:", length(unique(seqnames(gr))), "\n")
    cat("Unique cell barcodes:", length(unique(gr$cell_barcode)), "\n")
  }
  
  return(gr)
}

#### Chunk Processing Helper Function ----------
read_file_in_chunks <- function(file_path, cell_barcodes, chunk_size, verbose) {
  
  # Open file connection
  con <- file(file_path, "r")
  on.exit(close(con))
  
  filtered_chunks <- list()
  chunk_number <- 1
  
  # Define column types
  col_spec <- cols(
    X1 = col_character(),  # chromosome
    X2 = col_integer(),    # start
    X3 = col_integer(),    # end
    X4 = col_character()   # cell_barcode
  )
  
  # Read chunks
  repeat {
    chunk_lines <- readLines(con, n = chunk_size)
    
    if (length(chunk_lines) == 0) {
      break
    }
    
    # Parse the chunk
    chunk_df <- vroom(
      I(chunk_lines),  # Read from character vector
      delim = "\t",
      col_names = FALSE,
      col_select = 1:4,
      col_types = col_spec,
      altrep = FALSE,
      show_col_types = FALSE
    )
    
    # Filter by cell barcodes
    chunk_filtered <- chunk_df[chunk_df$X4 %in% cell_barcodes, ]
    
    if (nrow(chunk_filtered) > 0) {
      filtered_chunks[[chunk_number]] <- chunk_filtered
    }
    
    if (verbose && chunk_number %% 10 == 0) {
      cat("Processed", chunk_number, "chunks\n")
    }
    
    chunk_number <- chunk_number + 1
  }
  
  return(filtered_chunks)
}

### Create Fragment GRanges Object ----------
peaks_gr <- create_granges_from_tsv(
  file_path = fragpath,
  cell_barcodes = cell_barcodes,
  verbose = TRUE,
  chunk_size = 1000000  # Adjust chunk size as needed
)

gc()
gc()

## Fragment-Level QC ----------

### Initialize QC Data Frame ----------
identical(frag_summary$CB, cell_barcodes)
df_frag_qc_cell_lvl <- data.frame(row.names = cell_barcodes,
                                  fragments = frag_summary$frequency_count)
tmp <- peaks_gr$cell_barcode

### Fragments Within Peaks Analysis ----------
peak_frag <- overlapsAny(peaks_gr, granges(pbmc[["ATAC"]]))
peak_frag <- factor(peak_frag, levels = c(TRUE, FALSE))
peak_frag <- as.data.frame(table(peak_frag, tmp))
peak_frag <- peak_frag[peak_frag$peak_frag == TRUE,]
peak_frag <- peak_frag[match(cell_barcodes, peak_frag$tmp),]
gc()
gc()
df_frag_qc_cell_lvl$peak_frag <- peak_frag$Freq

### Fragments Within Blacklisted Regions Analysis ----------
bl_frag <- overlapsAny(peaks_gr, blacklist_hg38_unified)
bl_frag <- factor(bl_frag, levels = c(TRUE, FALSE))
bl_frag <- as.data.frame(table(bl_frag, tmp))
bl_frag <- bl_frag[bl_frag$bl_frag == TRUE,]
bl_frag <- bl_frag[match(cell_barcodes, bl_frag$tmp),]
gc()
gc()
df_frag_qc_cell_lvl$bl_frag <- bl_frag$Freq

### Fragments Within Standard Chromosomes Analysis ----------
chromosomes_to_keep <- c(paste0("chr", 1:22), "chrX")
std_chr_frag <- seqnames(peaks_gr) %in% chromosomes_to_keep
std_chr_frag <- as.vector(std_chr_frag)
std_chr_frag <- factor(std_chr_frag, levels = c(TRUE, FALSE))
std_chr_frag <- as.data.frame(table(std_chr_frag, tmp))
std_chr_frag <- std_chr_frag[std_chr_frag$std_chr_frag == TRUE,]
std_chr_frag <- std_chr_frag[match(cell_barcodes, std_chr_frag$tmp),]
gc()
gc()
df_frag_qc_cell_lvl$std_chr_frag <- std_chr_frag$Freq

### Fragments Within Mitochondrial Chromosome Analysis ----------
chromosomes_to_keep <- c("chrMT") # fragment information directly derived from cellranger output fragment.tsv.gz
chrmt_frag <- seqnames(peaks_gr) %in% chromosomes_to_keep
chrmt_frag <- as.vector(chrmt_frag)
chrmt_frag <- factor(chrmt_frag, levels = c(TRUE, FALSE))
chrmt_frag <- as.data.frame(table(chrmt_frag, tmp))
chrmt_frag <- chrmt_frag[chrmt_frag$chrmt_frag == TRUE,]
chrmt_frag <- chrmt_frag[match(cell_barcodes, chrmt_frag$tmp),]
gc()
gc()
df_frag_qc_cell_lvl$chrmt_frag <- chrmt_frag$Freq

### Fragments Within Repetitive Elements Analysis ----------

#### Define Repetitive Elements Analysis Function ----------
calculate_repetitive_element_frag <- function(peaks_gr,
                                              peaks_gr_barcodes,
                                              repeats_hg38,
                                              cell_barcodes,
                                              chr_change = TRUE,
                                              category){
  # pay attention to diff. between UCSC and Ensembl chr nomenclature
  if (chr_change) {
    seqlevels(repeats_hg38)[seqlevels(repeats_hg38) == "chrM"] <- "chrMT"
  }
  
  # Classify repeats by type
  high_confidence_removes <- c("Simple_repeat", "Low_complexity", "Satellite", 
                               "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA")
  
  moderate_removes <- c("LINE", "SINE", "LTR")  # More controversial
  
  keep_types <- c("DNA", "RC")  # Often contain regulatory elements
  
  # choose test_gr
  if (category == "all"){
    test_gr <- repeats_hg38
  } else if (category %in% c("high_confidence_removes", "moderate_removes", "keep_types")) {
    test_gr <- switch(category,
                     "high_confidence_removes" = repeats_hg38[repeats_hg38$repClass %in% high_confidence_removes,],
                     "moderate_removes" = repeats_hg38[repeats_hg38$repClass %in% moderate_removes,],
                     "keep_types" = repeats_hg38[repeats_hg38$repClass %in% keep_types,])
    
  } else if (category %in% c("Simple_repeat", "Low_complexity", "Satellite", 
                       "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA",
                       "LINE", "SINE", "LTR", "DNA", "RC")) {
    test_gr <- repeats_hg38[repeats_hg38$repClass == category,]
  } else {
    stop("Invalid category specified.")
  }
  
  # Select specific category of repetitive elements
  peak_frag <- overlapsAny(peaks_gr, test_gr)
  peak_frag <- factor(peak_frag, levels = c(TRUE, FALSE))
  peak_frag <- as.data.frame(table(peak_frag, peaks_gr_barcodes))
  peak_frag <- peak_frag[peak_frag$peak_frag == TRUE,]
  peak_frag <- peak_frag[match(cell_barcodes, peak_frag$peaks_gr_barcodes),]
  gc()
  gc()
  
  return(peak_frag)
}

#### Define Repetitive Element Categories ----------
# Classify repeats by type
high_confidence_removes <- c("Simple_repeat", "Low_complexity", "Satellite", 
                             "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA")

moderate_removes <- c("LINE", "SINE", "LTR")  # More controversial

keep_types <- c("DNA", "RC")  # Often contain regulatory elements

#### Calculate Repetitive Element Fragment Overlaps ----------
for (i in c("all", 
           "high_confidence_removes", 
           "moderate_removes", 
           "keep_types",
           "Simple_repeat", "Low_complexity", "Satellite", 
           "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA",
           "LINE", "SINE", "LTR", "DNA", "RC")) {

  rep_frag <- calculate_repetitive_element_frag(
    peaks_gr = peaks_gr,
    peaks_gr_barcodes = tmp,
    repeats_hg38 = repeats_hg38,
    cell_barcodes = cell_barcodes,
    chr_change = TRUE,
    category = i
  )

  df_frag_qc_cell_lvl[paste0(i, "_frag_rep")] <- rep_frag$Freq
}

### Nucleosome Signal and TSS Enrichment Analysis ----------
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)
df_frag_qc_cell_lvl$TSS.enrichment <- pbmc$TSS.enrichment
df_frag_qc_cell_lvl$nucleosome_signal <- pbmc$nucleosome_signal
df_frag_qc_cell_lvl$TSS.percentile <- pbmc$TSS.percentile
df_frag_qc_cell_lvl$nucleosome_percentile <- pbmc$nucleosome_percentile
df_frag_qc_cell_lvl$nCount_ATAC <- pbmc$nCount_ATAC
df_frag_qc_cell_lvl$nFeature_ATAC <- pbmc$nFeature_ATAC

### Calculate Fragment Ratios ----------
for (i in colnames(df_frag_qc_cell_lvl)[!(colnames(df_frag_qc_cell_lvl) %in% c("TSS.enrichment", "nucleosome_signal", "TSS.percentile", "nucleosome_percentile", "fragments", "nCount_ATAC", "nFeature_ATAC"))]){
  df_frag_qc_cell_lvl[[str_c(i, "_ratio")]] <- df_frag_qc_cell_lvl[[i]]/df_frag_qc_cell_lvl[["fragments"]]
}

### Save Results ----------
write.csv(
  df_frag_qc_cell_lvl,
  file = str_c(data_save, "cell_level_frag_qc.csv"),
  row.names = TRUE
)
## Doublets QC for atacseq -----------------
counts <- Read10X_h5(str_c(data_read,
                           "filtered_feature_bc_matrix.h5"))
# separate counts based on sample IDs (optional)
counts$Peaks <- counts$Peaks[, grepl("-2", colnames(counts$Peaks))]
counts$`Gene Expression` <- counts$`Gene Expression`[, grepl("-2", colnames(counts$`Gene Expression`))]

sce <- scDblFinder(SingleCellExperiment(list(counts = counts$Peaks)),
                   clusters = NULL,
                   aggregateFeatures = TRUE, nfeatures = 25,
                   processing = "normFeatures")
df_doublet_qc <- data.frame(row.names = colnames(counts$Peaks),
                            DF_atac_score = sce$scDblFinder.score,
                            DF_atac_classification = sce$scDblFinder.class)

write.csv(
  df_doublet_qc,
  file = str_c(data_save, "cell_level_doublet_qc.csv")
)

## QC Metrics Visualization ----------

### Clear Environment and Setup Paths ----------
rm(list = ls())
gc()

data_read <- "../data/"
analysis_save <- "../analysis/"
data_save <- "../data/"
sample_name <- "Ctrl"

### Load QC Data Files ----------
df_peak_qc_cell_lvl <- read.csv(
  str_c(data_read, "cell_level_peak_qc.csv"),
  row.name = 1
)

df_frag_qc_cell_lvl <- read.csv(
  str_c(data_read, "cell_level_frag_qc.csv"),
  row.name = 1
)

df_DF_qc_cell_lvl <- read.csv(
  str_c(data_save, "cell_level_doublet_qc.csv"),
  row.name = 1
)

### Combine QC Data ----------
df_atac_qc_cell_lvl <- cbind(df_peak_qc_cell_lvl, df_frag_qc_cell_lvl, df_DF_qc_cell_lvl)
df_atac_qc_cell_lvl$orig.ident <- sample_name 

### Define Visualization Themes ----------

#### Ridge Plot Theme ----------
ridgeplot_clean_theme <- function() {
  theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.line.x = element_line(color = "black", size = 0.5),
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "none"
    )
}

#### Bar Plot Theme ----------
barplot_clean_theme <- function() {
  theme_bw() +
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", 
                               arrow = arrow(length = unit(0.25, "cm"), type = "closed")),
      text = element_text(size = 20)
    )
}

### Cell-Level QC ----------

#### Doublet Detection Visualization ----------
df_atac_qc_cell_lvl %>%
  ggplot(aes(x=DF_atac_classification, fill=DF_atac_classification)) +
  geom_bar(stat = "count") +
  theme_classic() +
  xlab("Duplication rate") +
  NoLegend()
ggsave(file.path(analysis_save, paste0("QC_DF_atac.pdf")), width = 4, height = 7)

#### Fragments per Cell Distribution ----------
ggplot(df_atac_qc_cell_lvl, 
       aes(x=fragments,fill=orig.ident))+
  geom_density()+
  scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),  # original x values
    labels = c("0", "1", "10", "100", "1000", "10000")
  )+
  geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle("Distribution of fragments per cell")+
  xlab("Fragments per cell")+
  ridgeplot_clean_theme()
ggsave(
  str_c(analysis_save, "fragments_per_cell.pdf"),
  width = 8, height = 6
)

#### Peak Reads per Cell Distribution ----------
ggplot(df_atac_qc_cell_lvl, 
       aes(x=nCount_ATAC,fill=orig.ident))+
  geom_density()+
    scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),  # original x values
    labels = c("0", "1", "10", "100", "1000", "10000")
  )+
  geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 1)+
  geom_vline(xintercept = 100000, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle("Distribution of # of peak reads per cell")+
  xlab("Fragments per cell")+
  ridgeplot_clean_theme()
ggsave(
  str_c(analysis_save, "peak_reads_per_cell.pdf"),
  width = 8, height = 6
)

#### Peaks per Cell Distribution ----------
ggplot(df_atac_qc_cell_lvl, 
       aes(x=nFeature_ATAC,fill=orig.ident))+
  geom_density()+
    scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),  # original x values
    labels = c("0", "1", "10", "100", "1000", "10000")
  )+
  geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle("Distribution of # of peaks per cell")+
  xlab("Fragments per cell")+
  ridgeplot_clean_theme()
ggsave(
  str_c(analysis_save, "peak_per_cell.pdf"),
  width = 8, height = 6
)

#### Frag/Peak/Chr QC Metrics ----------
qc_metrics <- c("std_chr", 
               "chrmt", 
               "blacklist",
               "peak_frag", 
               "bl_frag",
               "std_chr_frag",
               "chrmt_frag", 
               "peak_frag_ratio",
               "bl_frag_ratio", 
               "std_chr_frag_ratio",
               "chrmt_frag_ratio",
               "nucleosome_signal", 
               "TSS.enrichment",
               "nucleosome_percentile",
               "TSS.percentile")

names(qc_metrics) <- c("Ratio of peaks within standard chromosomes",
                      "Ratio of peaks within chrMT",
                      "Ratio of peaks within blacklisted regions",
                      "# of fragments within peaks",
                      "# of fragments within blacklisted regions",
                      "# of fragments within standard chromosomes",
                      "# of fragments within chrMT",
                      "Ratio of fragments within peaks",
                      "Ratio of fragments within blacklisted regions",
                      "Ratio of fragments within standard chromosomes",
                      "Ratio of fragments within chrMT",
                      "Nucleosome signal",
                      "TSS enrichment",
                      "Nucleosome percentile",
                      "TSS percentile")

graphics.off()
for (i in unname(qc_metrics)) {
  cut_off <- switch(i,
                   "std_chr" = 0.95, 
                   "chrmt" = 0.05, 
                   "blacklist" = 0.01,
                   "peak_frag"= 500, 
                   "bl_frag"= NULL,
                   "std_chr_frag" = NULL,
                   "chrmt_frag" = NULL, 
                   "peak_frag_ratio" = 0.4,
                   "bl_frag_ratio" = 0.01, 
                   "std_chr_frag_ratio" = 0.95,
                   "chrmt_frag_ratio" = 0.05,
                   "nucleosome_signal" = 2, 
                   "TSS.enrichment" = 1,
                   "nucleosome_percentile" = 0.7,
                   "TSS.percentile" = 0.7
                   )
  
  linecolor <- switch(i,
                     "std_chr" = "red", 
                     "chrmt" = "red", 
                     "blacklist" = "red",
                     "peak_frag" = "blue", 
                     "bl_frag" = "",
                     "std_chr_frag" = "",
                     "chrmt_frag" = "", 
                     "peak_frag_ratio" = "blue",
                     "bl_frag_ratio" = "blue", 
                     "std_chr_frag_ratio" = "blue",
                     "chrmt_frag_ratio" = "blue",
                     "nucleosome_signal" = "red", 
                     "TSS.enrichment" = "red",
                     "nucleosome_percentile" = "blue",
                     "TSS.percentile" = "blue"
                     )
  
  if (length(unique(df_atac_qc_cell_lvl[[i]])) == 1) {
    print(paste("Skipping", i, "as all values are same."))
    next
  }

  ggplot(df_atac_qc_cell_lvl, 
         aes_string(x=i,fill="orig.ident"))+
    geom_density()+
    geom_vline(xintercept = cut_off, linetype = "dashed", color = linecolor, linewidth = 1)+
    NoLegend()+
    ggtitle(str_c("Distribution of data"))+
    xlab(names(qc_metrics)[qc_metrics == i])+
    ridgeplot_clean_theme()

  ggsave(
    str_c(analysis_save, i, ".pdf"),
    width = 8, height = 6
  )
}

#### Repetitive Elements QC ----------

qc_genomic_metrics <- c("all repetitive elements",
                               "repetitive elements(microsatellites, tandem repeats, etc.)",
                               "repetitive elements(LINEs, SINEs and LTRs)",
                               "repetitive elements(DNA transponsons and Rolling circle elements)",
                               "microsatellites",
                               "low complexity regions",
                               "satellite DNA",
                               "tRNA genes",
                               "rRNA genes",
                               "scRNA genes",
                               "snRNA genes",
                               "srpRNA genes",
                               "LINEs",
                               "SINEs",
                               "LTRs",
                               "DNA transposons",
                               "Rolling circle elements")

##### Peak Ratio ---------
names(qc_genomic_metrics) <- c("all_rep", 
                       "high_confidence_removes_rep", 
                       "moderate_removes_rep", 
                       "keep_types_rep",
                       "Simple_repeat_rep",
                       "Low_complexity_rep",
                       "Satellite_rep",
                       "tRNA_rep",
                       "rRNA_rep",
                       "scRNA_rep",
                       "snRNA_rep",
                       "srpRNA_rep",
                       "LINE_rep",
                       "SINE_rep",
                       "LTR_rep",
                       "DNA_rep",
                       "RC_rep")

df <- df_atac_qc_cell_lvl[,names(qc_genomic_metrics)]
df <- df[, !sapply(df, function(col) length(unique(col)) == 1)]
df <- tidyr::pivot_longer(df, 
                        cols = everything(), 
                        names_to = "QC_metric_rep", 
                        values_to = "Ratio")
df$QC_metric_rep <- factor(df$QC_metric_rep, 
                           levels = names(qc_genomic_metrics))

graphics.off()
ggplot(df, 
       aes_string(x="Ratio",
                  fill="QC_metric_rep"))+
  geom_density()+
  facet_wrap(~ QC_metric_rep, 
             labeller = labeller(QC_metric_rep = qc_genomic_metrics), 
             scales = "free_y", ncol=4) +
  ggtitle("Distribution of ratios of peaks within genomic repetitive elements")+
  xlab("")+
  ridgeplot_clean_theme()

ggsave(
  str_c(analysis_save, "Raio_peaks_of_genomic_rep_summary_cell_lvl.pdf"),
  width = 20, height = 8
)

##### Frag # ---------
names(qc_genomic_metrics) <- c("all_frag_rep", 
                              "high_confidence_removes_frag_rep", 
                              "moderate_removes_frag_rep", 
                              "keep_types_frag_rep",
                              "Simple_repeat_frag_rep",
                              "Low_complexity_frag_rep",
                              "Satellite_frag_rep",
                              "tRNA_frag_rep",
                              "rRNA_frag_rep",
                              "scRNA_frag_rep",
                              "snRNA_frag_rep",
                              "srpRNA_frag_rep",
                              "LINE_frag_rep",
                              "SINE_frag_rep",
                              "LTR_frag_rep",
                              "DNA_frag_rep",
                              "RC_frag_rep")

df <- df_atac_qc_cell_lvl[, names(qc_genomic_metrics)]
df <- df[, !sapply(df, function(col) length(unique(col)) == 1)]
df <- tidyr::pivot_longer(df, 
                         cols = everything(), 
                         names_to = "QC_metric_rep", 
                         values_to = "Ratio")
df$QC_metric_rep <- factor(df$QC_metric_rep, 
                           levels = names(qc_genomic_metrics))

graphics.off()
ggplot(df,
       aes_string(x="Ratio",
                  fill="QC_metric_rep"))+
  geom_density()+
  scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    labels = c("0", "1", "10", "100", "1000", "10000")
  )+
  facet_wrap(~ QC_metric_rep, 
             labeller = labeller(QC_metric_rep = qc_genomic_metrics), 
             scales = "free_y", ncol=4) +
  
  ggtitle("Distribution of # of fragments within genomic repetitive elements")+
  xlab("")+
  ridgeplot_clean_theme()

ggsave(
  str_c(analysis_save, "Number_frag_of_genomic_rep_summary_cell_lvl.pdf"),
  width = 20, height = 8
)

##### Frag Ratio ---------
names(qc_genomic_metrics) <- c("all_frag_rep_ratio", 
                              "high_confidence_removes_frag_rep_ratio", 
                              "moderate_removes_frag_rep_ratio", 
                              "keep_types_frag_rep_ratio",
                              "Simple_repeat_frag_rep_ratio",
                              "Low_complexity_frag_rep_ratio",
                              "Satellite_frag_rep_ratio",
                              "tRNA_frag_rep_ratio",
                              "rRNA_frag_rep_ratio",
                              "scRNA_frag_rep_ratio",
                              "snRNA_frag_rep_ratio",
                              "srpRNA_frag_rep_ratio",
                              "LINE_frag_rep_ratio",
                              "SINE_frag_rep_ratio",
                              "LTR_frag_rep_ratio",
                              "DNA_frag_rep_ratio",
                              "RC_frag_rep_ratio")

df <- df_atac_qc_cell_lvl[, names(qc_genomic_metrics)]
df <- df[, !sapply(df, function(col) length(unique(col)) == 1)]
df <- tidyr::pivot_longer(df, 
                         cols = everything(), 
                         names_to = "QC_metric_rep", 
                         values_to = "Ratio")
df$QC_metric_rep <- factor(df$QC_metric_rep, 
                           levels = names(qc_genomic_metrics))

graphics.off()
ggplot(df, 
       aes_string(x="Ratio",
                  fill="QC_metric_rep"))+
  geom_density()+
  facet_wrap(~ QC_metric_rep, 
             labeller = labeller(QC_metric_rep = qc_genomic_metrics), 
             scales = "free_y", ncol=4) +
  
  ggtitle("Distribution of # of fragments within genomic repetitive elements")+
  xlab("")+
  ridgeplot_clean_theme()

ggsave(
  str_c(analysis_save, "Ratio_frag_of_genomic_rep_summary_cell_lvl.pdf"),
  width = 20, height = 8
)
### Peak-Level QC ----------

df_peak_qc <- read.csv(str_c(data_save, "peak_level_peak_qc_summary.csv"), row.names = 1)
df_peak_qc$peak_cat <- rownames(df_peak_qc)

peak_qc <- c("Ratio of peaks within standard chromosomes",
                   "Ratio of peaks within blacklisted regions",
                   "Ratio of peaks within chrMT",
                   "Ratio of peaks within all repetitive elements",
                   "Ratio of peaks within repetitive elements(microsatellites, tandem repeats, etc.)",
                   "Ratio of peaks within repetitive elements(LINEs, SINEs and LTRs)",
                   "Ratio of peaks within other repetitive elements(DNA transponsons and Rolling circle elements)",
                   "Ratio of peaks within microsatellites",
                   "Ratio of peaks within low complexity regions",
                   "Ratio of peaks within satellite DNA",
                   "Ratio of peaks within tRNA genes",
                   "Ratio of peaks within rRNA genes",
                   "Ratio of peaks within scRNA genes",
                   "Ratio of peaks within snRNA genes",
                   "Ratio of peaks within srpRNA genes",
                   "Ratio of peaks within LINEs",
                   "Ratio of peaks within SINEs",
                   "Ratio of peaks within LTRs",
                   "Ratio of peaks within DNA transposons",
                   "Ratio of peaks within Rolling circle elements")

names(peak_qc) <- c("std_chr", 
            "blacklist", 
            "chrmt",
            "all_rep",
            "high_confidence_removes_rep",
            "moderate_removes_rep", 
            "keep_types_rep",
            "Simple_repeat_rep",
            "Low_complexity_rep",
            "Satellite_rep",
            "tRNA_rep",
            "rRNA_rep",
            "scRNA_rep",
            "snRNA_rep",
            "srpRNA_rep",
            "LINE_rep",
            "SINE_rep",
            "LTR_rep",
            "DNA_rep",
            "RC_rep")

df_peak_qc$peak_cat <- factor(df_peak_qc$peak_cat, 
                             levels = names(peak_qc))

ggplot(df_peak_qc, 
       aes_string(x="peak_cat", y="peak_ratio", fill = "peak_cat"))+
  geom_bar(stat = "identity")+
  ggtitle("Ratio of peaks within different regions")+
  xlab("")+
  ylab("Ratio")+
  scale_fill_discrete(labels=peak_qc)+
  scale_y_continuous(limits = c(0,1))+
  barplot_clean_theme()+
  theme(axis.text.x = element_blank())
ggsave(
  str_c(analysis_save, "Ratio_peaks_of_genomic_regions_summary_peak_lvl.pdf"),
  width = 20, height = 6
)

#### Load Detailed Peak-Level Data ----------
df_peak_qc_peak_lvl <- read.csv(str_c(data_save, "peak_level_peak_qc.csv"), row.names = 1)
df_peak_qc_peak_lvl$orig.ident <- sample_name 

#### Cell Expression Ratio per Peak ----------
ggplot(df_peak_qc_peak_lvl, 
       aes_string(x="peak_cell_ratio",fill="orig.ident"))+
  geom_density()+
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle(str_c("Distribution of ratios of cells expressing corresponding peaks"))+
  xlab("Ratio of cells per peak")+
  ridgeplot_clean_theme()
ggsave(str_c(analysis_save, "cell_per_peak_ratio_peak_lvl.pdf"),
       width = 8, height = 6)

#### # of Cells per Peak ----------
ggplot(df_peak_qc_peak_lvl, 
       aes_string(x="peak_cell",fill="orig.ident"))+
  geom_density()+
  geom_vline(xintercept = 5, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle(str_c("Distribution of # of cells expressing corresponding peaks"))+
  xlab("# of cells per peak")+
  ridgeplot_clean_theme()
ggsave(str_c(analysis_save, "cell_per_peak_peak_lvl.pdf"),
       width = 8, height = 6)

#### # of Reads per Peak ----------
ggplot(df_peak_qc_peak_lvl, 
       aes_string(x="peak_reads",fill="orig.ident"))+
  geom_density()+
  geom_vline(xintercept = 100, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle(str_c("Distribution of # of reads for peaks"))+
  xlab("# of reads for peaks")+
  scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    labels = c("0", "1", "10", "100", "1000", "10000")
  )+
  ridgeplot_clean_theme()
ggsave(str_c(analysis_save, "peak_reads_peak_lvl.pdf"),
       width = 8, height = 6)

#### Peak Width Distribution ----------
ggplot(df_peak_qc_peak_lvl, 
       aes_string(x="peak_width",fill="orig.ident"))+
  geom_density()+
  geom_vline(xintercept = 200, linetype = "dashed", color = "blue", linewidth = 1)+
  NoLegend()+
  ggtitle(str_c("Distribution of widths of peaks"))+
  xlab("peak widths")+
  scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    labels = c("0", "1", "10", "100", "1000", "10000")
  )+
  ridgeplot_clean_theme()
ggsave(str_c(analysis_save, "peak_width_peak_lvl.pdf"),
       width = 8, height = 6)

## Filter cells and peaks based on QC metrics --------------
df_atac_qc_summary = data.frame(row.names = unique(df_atac_qc_cell_lvl$orig.ident),
                                nCells_before = nrow(df_atac_qc_cell_lvl),
                                nCells_after = NA,
                                nPeaks = nrow(df_peak_qc_peak_lvl),
                                nPeaks_after = NA)
keep_atacseq = df_atac_qc_cell_lvl$fragments > 1000 & 
  df_atac_qc_cell_lvl$nCount_ATAC > 1000 & 
  df_atac_qc_cell_lvl$nCount_ATAC < 100000 &
  df_atac_qc_cell_lvl$nFeature_ATAC > 1000 & 
  df_atac_qc_cell_lvl$peak_frag_ratio > 0.4 & 
  df_atac_qc_cell_lvl$bl_frag_ratio < 0.01 & 
  df_atac_qc_cell_lvl$std_chr_frag_ratio > 0.95 & 
  df_atac_qc_cell_lvl$chrmt_frag_ratio < 0.05 & 
  df_atac_qc_cell_lvl$std_chr > 0.95 & 
  df_atac_qc_cell_lvl$blacklist < 0.01 & 
  df_atac_qc_cell_lvl$chrmt < 0.05 &
  df_atac_qc_cell_lvl$nucleosome_signal < 2 &
  df_atac_qc_cell_lvl$TSS.enrichment > 1 &
  df_atac_qc_cell_lvl$DF_atac_classification == "singlet"

df_cell_atac_qc = data.frame(
  row.names = rownames(df_atac_qc_cell_lvl),
  keep_atacseq_cell = keep_atacseq
)

keep_atacseq_peak = df_peak_qc_peak_lvl$peak_cell_ratio > 0.01 & 
  df_peak_qc_peak_lvl$peak_cell > 5 & 
  df_peak_qc_peak_lvl$peak_reads > 100 & 
  df_peak_qc_peak_lvl$std_chr == TRUE & 
  df_peak_qc_peak_lvl$blacklist == FALSE & 
  df_peak_qc_peak_lvl$chrmt == FALSE

df_peak_atac_qc = data.frame(
  row.names = rownames(df_peak_qc_peak_lvl),
  keep_atacseq_peak = keep_atacseq_peak
)

df_atac_qc_summary['nCells_after'] <- sum(keep_atacseq)
df_atac_qc_summary['nPeaks_after'] <- sum(keep_atacseq_peak)

## Save Results ----------
write.csv(
  df_atac_qc_summary,
  file = str_c(data_save, "df_atac_qc_summary.csv"),
  row.names = TRUE
)

write.csv(
  df_cell_atac_qc,
  file = str_c(data_save, "df_cell_atac_keep.csv"),
  row.names = TRUE
)

write.csv(
  df_peak_atac_qc,
  file = str_c(data_save, "df_peak_atac_keep.csv"),
  row.names = TRUE
)

# QC for scRNAseq ----------
## Setup and Data Loading ----------
### Clear Environment and Setup Paths ----------
rm(list = ls())
gc()

data_read <- "../data/"
analysis_save <- "../analysis/"
data_save <- "../data/"
sample_name <- "Ctrl"

### Load 10X Data ----------
counts <- Read10X_h5(str_c(data_read,
                           "filtered_feature_bc_matrix.h5"))

### Separate Counts by Sample IDs ----------
# separate counts based on sample IDs (optional)
counts$Peaks <- counts$Peaks[, grepl("-2", colnames(counts$Peaks))]
counts$`Gene Expression` <- counts$`Gene Expression`[, grepl("-2", colnames(counts$`Gene Expression`))]

### Create Seurat Object ----------
seurat_obj <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA",
  min.cells = 0,
  min.features = 0
)

### Setup Gene Annotations ----------
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

## Cell-Level QC ----------

### Calculate Basic QC Metrics ----------
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA) # nolint: line_length_linter.
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-") # nolint: line_length_linter.
seurat_obj$mitoRatio <- seurat_obj$mitoRatio / 100

### Prepare Data Frame for Visualization ----------
df <- seurat_obj@meta.data
df <- df %>% dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)

### Cell-Level QC Visualizations ----------

#### UMI Count Distribution ----------
p1 <- ggplot(df, aes(color=orig.ident, x=nUMI, fill=orig.ident)) +
  geom_density() +
    scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    labels = c("0", "1", "10", "100", "1000", "10000")
  ) +
  theme_classic() +
  xlab("# of UMIs") +
  NoLegend()+
  geom_vline(xintercept = 500)
ggsave(file.path(analysis_save, paste0("QC_Cell_Density_UMI.pdf")),
       plot = p1, width = 8, height = 7)

#### Gene Count Distribution ----------
p2 <- df %>% 
  ggplot(aes(color=orig.ident, x=nGene, fill=orig.ident)) + 
  geom_density() + 
  theme_classic() +
    scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    labels = c("0", "1", "10", "100", "1000", "10000")
  ) +
  xlab("# of genes")+
  NoLegend()+
  geom_vline(xintercept = 300)
ggsave(file.path(analysis_save, paste0("QC_Cell_Density_genes.pdf")),
       plot = p2, width = 8, height = 7)

#### Mitochondrial Ratio Distribution ----------
p3 <- df %>% 
  ggplot(aes(color=orig.ident, x=mitoRatio, fill=orig.ident)) + 
  geom_density() +
  theme_classic() +
  xlab("Mitochondrial ratio") +
  NoLegend()+
  geom_vline(xintercept = 0.25)
ggsave(file.path(analysis_save, paste0("QC_Cell_Density_mitoRatio.pdf")),
       plot = p3, width = 8, height = 7)

#### UMI vs Gene Correlation Plot ----------
p4 <- df %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point(aes(size=log10GenesPerUMI)) + 
  scale_colour_gradient(low = "gray90", high = "black") +
  geom_smooth(method = "lm") +
  scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    labels = c("0", "1", "10", "100", "1000", "10000")
  ) +
  scale_y_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    labels = c("0", "1", "10", "100", "1000", "10000")
  ) +
  theme_classic() +
  xlab("# of UMIs") +
  ylab("# of genes") +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 300)
ggsave(file.path(analysis_save, paste0("QC_Cell_correlations.pdf")),
       plot = p4, width = 8, height = 7)

#### Gene Complexity Distribution ----------
p5 <- df %>%
  ggplot(aes(x=log10GenesPerUMI, color=orig.ident, fill=orig.ident)) +
  geom_density() +
  theme_classic() +
  xlab("log10(Genes per UMI)") +
  NoLegend()+
  geom_vline(xintercept = 0.8)
ggsave(file.path(analysis_save, paste0("QC_Cell_complexity.pdf")),
       plot = p5, width = 8, height = 7)

### Doublet Detection Analysis ----------

#### Perform Standard Processing for Doublet Detection ----------
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj)

#### Optimize DoubletFinder Parameters ----------
sweep.res <- paramSweep(seurat_obj, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

best_pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
best_pK <- as.numeric(as.character(best_pK[[1]]))

#### Estimate Multiplet Rate ----------
seu_clusters <- seurat_obj@meta.data$seurat_clusters # use the clusters as the user-defined cell types
homotypic.prop <- modelHomotypic(seu_clusters)

multiplet_rate <- NULL
if(is.null(multiplet_rate)){
  print("multiplet_rate not provided....... estimating multiplet rate from cells in dataset")
  
  # 10X multiplet rates table
  #https://rpubs.com/kenneditodd/doublet_finder_example
  multiplet_rates_10x <- data.frame("Multiplet_rate"= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                    "Loaded_cells" = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                    "Recovered_cells" = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
  
  print(multiplet_rates_10x)
  
  multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(seurat_obj@meta.data)) %>% 
    dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
    dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells

  print(paste("Setting multiplet rate to", multiplet_rate))
}

#### Run DoubletFinder ----------
nExp.poi <- round(multiplet_rate * nrow(seurat_obj@meta.data)) # multiply by number of cells to get the number of expected multiplets
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets

seurat_obj <- doubletFinder(seu = seurat_obj,
                        PCs = 1:10,
                        pK = best_pK,
                        nExp = nExp.poi.adj)
df$DF <- seurat_obj@meta.data[, grepl("DF.classifications", colnames(seurat_obj@meta.data))]

#### Visualize Doublet Detection Results ----------
p6 <- df %>%
  ggplot(aes(x=DF, fill=DF)) +
  geom_bar(stat = "count") +
  theme_classic() +
  xlab("Duplication rate") +
  NoLegend()
ggsave(file.path(analysis_save, paste0("QC_DF.pdf")),
       plot = p6, width = 4, height = 7)

### Save Results ----------
write.csv(
  df,
  file = str_c(data_save, "cell_level_rna_qc.csv"),
  row.names = TRUE
)

## Gene-Level QC ----------

### Calculate Gene-Level QC Metrics ----------
counts <- GetAssayData(object = seurat_obj, layer = "counts")
nonzero <- counts > 0

df_gene_qc <- data.frame(
  row.names = rownames(counts),
  nCells = Matrix::rowSums(nonzero),
  nCells_ratio = Matrix::rowSums(nonzero) / ncol(counts),
  nCounts = Matrix::rowSums(counts)
)

### Annotate Genes by Genomic Location ----------
tmp <- annotation[seqnames(annotation) == "chrMT"]
df_gene_qc$mito <- rownames(df_gene_qc) %in% tmp$gene_name
tmp <- annotation[seqnames(annotation) %in% c(paste0("chr", 1:22), "chrX")]
df_gene_qc$std_chr <- rownames(df_gene_qc) %in% tmp$gene_name
df_gene_qc$orig.ident <- sample_name

### Generate Gene-Level QC Visualizations ----------

#### # of Cells Expressing Each Gene ----------
p1 <- df_gene_qc %>%
  ggplot(aes(x=nCells, color=orig.ident, fill=orig.ident)) +
  geom_density() +
  theme_classic() +
  scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    labels = c("0", "1", "10", "100", "1000", "10000")
  ) +
  NoLegend()+
  xlab("# of cells expressing gene") +
  geom_vline(xintercept = 10)
ggsave(file.path(analysis_save, paste0("QC_gene_nCells.pdf")),
       plot = p1, width = 8, height = 7)

#### Ratio of Cells Expressing Each Gene ----------
p2 <- df_gene_qc %>%
  ggplot(aes(x=nCells_ratio, color=orig.ident, fill=orig.ident)) +
  geom_density() +
  theme_classic() +
  scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    labels = c("0", "1", "10", "100", "1000", "10000")
  ) +
  NoLegend()+
  xlab("Ratio of cells expressing gene") +
  geom_vline(xintercept = 0.01)
ggsave(file.path(analysis_save, paste0("QC_gene_nCells_ratio.pdf")),
       plot = p2, width = 8, height = 7)

#### Total Counts per Gene ----------
p3 <- df_gene_qc %>%
  ggplot(aes(x=nCounts, color=orig.ident, fill=orig.ident)) +
  geom_density() +
  theme_classic() +
  scale_x_continuous(
    trans = scales::log1p_trans(),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    labels = c("0", "1", "10", "100", "1000", "10000")
  ) +
  NoLegend()+
  xlab("# of counts for gene") +
  geom_vline(xintercept = 10)
ggsave(file.path(analysis_save, paste0("QC_gene_nCounts.pdf")),
       plot = p3, width = 8, height = 7)

## Filter Cells and Genes Based on QC Metrics ----------

### Initialize QC Summary ----------
df_rna_qc_summary <- data.frame(row.names = sample_name,
                                nCells_before = nrow(df),
                                nCells_after = NA,
                                nGenes = nrow(df_gene_qc),
                                nGenes_after = NA)

### Define Cell Filtering Criteria ----------
keep_rnaseq <- df$nUMI > 500 &
  df$nGene > 300 &
  df$mitoRatio < 0.25 &
  df$log10GenesPerUMI > 0.8 &
  df$DF == "Singlet"

### Create Cell Filtering Data Frame ----------
df_cell_rna_qc <- data.frame(row.names = rownames(df),
                            keep_rnaseq_cell = keep_rnaseq)

### Define Gene Filtering Criteria ----------
df_gene_rna_qc <- data.frame(
  row.names = rownames(df_gene_qc),
  keep_rnaseq_gene = df_gene_qc$nCells_ratio > 0.01 & 
    df_gene_qc$nCells > 10 & 
    df_gene_qc$nCounts > 10 & 
    df_gene_qc$std_chr == TRUE & 
    df_gene_qc$mito == FALSE
)

### Update QC Summary ----------
df_rna_qc_summary["nCells_after"] <- sum(keep_rnaseq)
df_rna_qc_summary["nGenes_after"] <- sum(df_gene_rna_qc$keep_rnaseq_gene)

## Save QC Results ----------

### Save QC Summary ----------
write.csv(
  df_rna_qc_summary,
  file = str_c(data_save, "df_rna_qc_summary.csv"),
  row.names = TRUE
)

### Save Cell Filtering Results ----------
write.csv(
  df_cell_rna_qc,
  file = str_c(data_save, "df_cell_rna_keep.csv"),
  row.names = TRUE
)

### Save Gene Filtering Results ----------
write.csv(
  df_gene_rna_qc,
  file = str_c(data_save, "df_gene_rna_keep.csv"),
  row.names = TRUE
)

# Integrated filtering based on scRNAseq and scATACseq ----------
rm(list = ls())
gc()

data_read <- "../data/"
analysis_save <- "../analysis/"
data_save <- "../data/"
sample_name <- "Ctrl"
counts <- Read10X_h5(str_c(data_read,
                           "filtered_feature_bc_matrix.h5"))

# separate counts based on sample IDs (optional)
counts$Peaks <- counts$Peaks[, grepl("-2", colnames(counts$Peaks))]
counts$`Gene Expression` <- counts$`Gene Expression`[, grepl("-2", colnames(counts$`Gene Expression`))]

df_cell_rna_keep = read.csv(
  str_c(data_save, "df_cell_rna_keep.csv"),
  row.names = 1
)

df_cell_atac_keep = read.csv(
  str_c(data_save, "df_cell_atac_keep.csv"),
  row.names = 1
)

df_gene_rna_keep = read.csv(
  str_c(data_save, "df_gene_rna_keep.csv"),
  row.names = 1
)

df_peak_atac_keep = read.csv(
  str_c(data_save, "df_peak_atac_keep.csv"),
  row.names = 1
)

## check rownames and colnames are identical ---------
identical(
  rownames(df_cell_rna_keep), 
  rownames(df_cell_atac_keep)
)

identical(
  rownames(df_cell_rna_keep), 
  colnames(counts$`Gene Expression`)
)

identical(
  rownames(df_peak_atac_keep), 
  rownames(counts$`Peaks`)
)

identical(
  rownames(df_gene_rna_keep), 
  rownames(counts$`Gene Expression`)
)

## filter cells based on scRNAseq and scATACseq QC metrics ---------
df_integrated_qc_summary = data.frame(
  row.names = "Ctrl",
  nCells = nrow(df_cell_rna_keep),
  nCells_after = NA,
  nPeaks = nrow(df_peak_atac_keep),
  nPeaks_after = NA,
  nGenes = nrow(df_gene_rna_keep),
  nGenes_after = NA
)

keep_cells = df_cell_rna_keep$keep_rnaseq_cell & 
  df_cell_atac_keep$keep_atacseq_cell

keep_genes = df_gene_rna_keep$keep_rnaseq_gene

keep_peaks = df_peak_atac_keep$keep_atacseq_peak

keep_gene_expr = counts$`Gene Expression`[keep_genes, keep_cells]
keep_peak_expr = counts$`Peaks`[keep_peaks, keep_cells]

df_integrated_qc_summary['nCells_after'] <- sum(keep_cells)
df_integrated_qc_summary['nPeaks_after'] <- sum(keep_peaks)
df_integrated_qc_summary['nGenes_after'] <- sum(keep_genes)

write.csv(
  df_integrated_qc_summary,
  file = str_c(data_save, "df_integrated_qc_summary.csv"),
  row.names = TRUE
)


## create a filtered seurat object ---------
seuratObj <- CreateSeuratObject(
  counts = keep_gene_expr,
  assay = "RNA",
  min.cells = 0,
  min.features = 0
)

seuratObj$orig.ident <- sample_name

fragpath <- str_c(data_read,"atac_fragments.tsv.gz")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
seuratObj[["ATAC"]] <- CreateChromatinAssay(
  counts = keep_peak_expr,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

saveRDS(
  seuratObj,
  file = str_c(data_save, "seuratObj_filtered.rds")
)

# Integrated QC Visualization ----------

## Define Visualization Theme ----------
barplot_clean_theme <- function() {
  theme_bw() +
    theme(
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", 
                               arrow = arrow(length = unit(0.25, "cm"), type = "closed")),
      text = element_text(size = 20)
    )
}

## Data Preparation for Visualization ----------

### Transform and Structure QC Summary Data ----------
df_integrated_qc_summary <- as.data.frame(t(df_integrated_qc_summary))
df_integrated_qc_summary$metric <- rownames(df_integrated_qc_summary)
df_integrated_qc_summary$metric <- factor(df_integrated_qc_summary$metric, 
                                          levels = c("nCells", "nCells_after", 
                                                     "nPeaks", "nPeaks_after", 
                                                     "nGenes", "nGenes_after"))

### Add Category Labels ----------
df_integrated_qc_summary$cat <- c("Cells",
                                  "Cells",
                                  "Peaks",
                                  "Peaks",
                                  "Genes",
                                  "Genes")

## Generate Category-Specific Visualizations ----------

### Define Common Legend Labels ----------
legend_labels <- c( 
  "nCells" = "# of cells before filtering", 
  "nCells_after" = "# of cells after filtering", 
  "nPeaks" = "# of peaks before filtering", 
  "nPeaks_after" = "# of peaks after filtering", 
  "nGenes" = "# of genes before filtering", 
  "nGenes_after" = "# of genes after filtering" 
)

### Create Cells Category Plot ----------
p1 <- df_integrated_qc_summary %>%
  filter(cat == "Cells") %>%
  ggplot(aes_string(x="metric", y="Ctrl", fill = "metric")) +
  geom_bar(stat = "identity") +
  ggtitle("Cells") +
  xlab("") +
  ylab("") +
  scale_fill_discrete(labels = legend_labels) +
  barplot_clean_theme() +
  theme(axis.text.x = element_blank()) +
  scale_x_discrete(drop = TRUE)

### Create Genes Category Plot ----------
p2 <- df_integrated_qc_summary %>%
  filter(cat == "Genes") %>%
  ggplot(aes_string(x="metric", y="Ctrl", fill = "metric")) +
  geom_bar(stat = "identity") +
  ggtitle("Genes") +
  xlab("") +
  ylab("") +
  scale_fill_discrete(labels = legend_labels) +
  barplot_clean_theme() +
  theme(axis.text.x = element_blank()) +
  scale_x_discrete(drop = TRUE)

### Create Peaks Category Plot ----------
p3 <- df_integrated_qc_summary %>%
  filter(cat == "Peaks") %>%
  ggplot(aes_string(x="metric", y="Ctrl", fill = "metric")) +
  geom_bar(stat = "identity") +
  ggtitle("Peaks") +
  xlab("") +
  ylab("") +
  scale_fill_discrete(labels = legend_labels) +
  barplot_clean_theme() +
  theme(axis.text.x = element_blank()) +
  scale_x_discrete(drop = TRUE)

## Combine and Save Visualization ----------
combined_plot <- p1 + p2 + p3 + 
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(title = "Summary of QC metrics after integrated filtering")

ggsave(
  str_c(analysis_save, "integrated_qc_summary.pdf"),
  plot = combined_plot,
  width = 12, height = 6
)
