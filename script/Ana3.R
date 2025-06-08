# Detailed Joint scATAC-seq&scRNAseq analysis with Seurat and Signac
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
library(ggplot2)



# Load the data ------------
data_read = "../data/"
analysis_save = "../analysis/"
data_save = "../data/"
counts <- Read10X_h5(str_c(data_read,"pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"))
fragpath <- str_c(data_read,"pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")

# QC for scATACseq ----------
# QC using peaks -----------
peak_names <- rownames(counts$Peaks)

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
    ranges = IRanges(start = start, end = end)
  )
  
  return(gr)
}

peak_gr <- parse_coordinates(peak_names)
peak_chrs <- sapply(strsplit(peak_names, ":"), `[`, 1)

## cell level ---------
### check peaks within standard chr --------------
# Define chromosomes to keep
chromosomes_to_keep <- c(paste0("chr", 1:22), "chrX")

# Create logical vector for peaks to keep
peaks_keep <- peak_chrs %in% chromosomes_to_keep

peak_std_chr_ratio = Matrix::colSums(keep_peaks <- counts$Peaks[peaks_keep, ])/Matrix::colSums(counts$Peaks)
df_peak_qc = as.data.frame(peak_std_chr_ratio)

### check peaks within genomic repetitive elements (filtering is controversial)--------------
ah <- AnnotationHub()
query(ah, c("RepeatMasker", "Homo sapiens"))
repeats_hg38 <- ah[["AH99003"]]
repeats_hg38 <- repeats_hg38[!grepl("\\?", repeats_hg38$repClass),]

summarize_rep_func1 = function(peak_gr, 
                               repeats_hg38, 
                               peak_matrix, 
                               Category='all'){
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

summarize_rep_func2 = function(peak_gr, 
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


peak_ratio_all <- summarize_rep_func1(
  peak_gr = peak_gr,
  repeats_hg38 = repeats_hg38,
  peak_matrix = counts$Peaks,
  Category = 'all'
)

peak_ratio_high_confidence_removes <- summarize_rep_func1(
  peak_gr = peak_gr,
  repeats_hg38 = repeats_hg38,
  peak_matrix = counts$Peaks,
  Category = 'high_confidence_removes'
)

peak_ratio_moderate_removes <- summarize_rep_func1(
  peak_gr = peak_gr,
  repeats_hg38 = repeats_hg38,
  peak_matrix = counts$Peaks,
  Category = 'moderate_removes'
)

peak_ratio_keep_types <- summarize_rep_func1(
  peak_gr = peak_gr,
  repeats_hg38 = repeats_hg38,
  peak_matrix = counts$Peaks,
  Category = 'keep_types'
)

df_rep <- rbind(
  peak_ratio_all,
  peak_ratio_high_confidence_removes,
  peak_ratio_moderate_removes,
  peak_ratio_keep_types
)
rownames(df_rep) = c(
  "all_rep",
  "high_confidence_removes_rep",
  "moderate_removes_rep",
  "keep_types_rep"
)
df_rep = as.data.frame(t(df_rep))
df_peak_qc <- cbind(df_peak_qc, df_rep)

tmp = summarize_rep_func2(
  peak_gr = peak_gr,
  repeats_hg38 = repeats_hg38,
  peak_matrix = counts$Peaks
)

df_peak_qc <- cbind(df_peak_qc, tmp)


### check peaks within MT--------------
peak_chrmt_ratio = Matrix::colSums(counts$Peaks[peak_chrs %in% c("chrMT"), ])/Matrix::colSums(counts$Peaks)
df_peak_qc = cbind(df_peak_qc,as.data.frame(peak_chrmt_ratio))


### check peaks within blacklisted regions --------------
peaks_keep <- overlapsAny(peak_gr, blacklist_hg38_unified)
filtered_peaks <- counts$Peaks[peaks_keep, ]
peak_ratio = Matrix::colSums(filtered_peaks)/Matrix::colSums(counts$Peaks)
df_peak_qc = cbind(df_peak_qc,as.data.frame(peak_ratio))
colnames(df_peak_qc)[ncol(df_peak_qc)] <- 'peak_black_ratio'

### save peak QC results for cells -----------
write.csv(
  df_peak_qc,
  file = str_c(data_save, "cell_level_peak_qc.csv"),
  row.names = TRUE
)



## peak level ---------
df_peak_qc = data.frame(row.names = 'peak_ratio')
df_peak_qc_peak_name = data.frame(row.names = rownames(counts$Peaks))

### check peaks within standard chr --------------
chromosomes_to_keep = c(paste0("chr", 1:22), "chrX")

peaks_keep <- seqnames(peak_gr) %in% chromosomes_to_keep
peak_ratio = sum(peaks_keep)/length(peaks_keep)
df_peak_qc['peak_std_chr_ratio'] <- peak_ratio
df_peak_qc_peak_name['peak_std_chr'] <- as.vector(peaks_keep)


### check peaks within genomic repetitive elements --------------
summarize_rep_func3 = function(peak_gr, 
                               repeats_hg38, 
                               peak_matrix, 
                               Category='all'){
  
  df = data.frame(row.names = 'peak_raio')
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
  
  
  peak_ratio = sum(peaks_keep)/length(peaks_keep)
  
  df[paste0(i,'_rep')] <- peak_ratio
  
  return(df)
}

summarize_rep_func4 = function(peak_gr, 
                               repeats_hg38, 
                               peak_matrix){
  # all types
  categories <- c("Simple_repeat", "Low_complexity", "Satellite", 
                  "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA",
                  "LINE", "SINE", "LTR", "DNA", "RC")
  
  df = data.frame(row.names = 'peak_raio')
  for (i in categories) {
    # Select specific category of repetitive elements
    repeats_subset <- repeats_hg38[repeats_hg38$repClass == i, ]
    
    
    # Create logical vector for peaks to keep
    peaks_keep <- overlapsAny(peak_gr, repeats_subset)
    
    peak_ratio = sum(peaks_keep)/length(peaks_keep)
    
    df[paste0(i,'_rep')] <- peak_ratio
  }
  
  
  return(df)
}

summarize_rep_func5 = function(peak_gr, 
                               repeats_hg38, 
                               peak_matrix,
                               df_peak_qc_peak_name,
                               Category='all'){
  
  df = df_peak_qc_peak_name
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
  df[paste0(i,'_rep')] <- as.vector(peaks_keep)
  
  return(df)
}

summarize_rep_func6 = function(peak_gr, 
                               repeats_hg38,
                               df_peak_qc_peak_name,
                               peak_matrix){
  # all types
  categories <- c("Simple_repeat", "Low_complexity", "Satellite", 
                  "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA",
                  "LINE", "SINE", "LTR", "DNA", "RC")
  
  df = df_peak_qc_peak_name
  for (i in categories) {
    # Select specific category of repetitive elements
    repeats_subset <- repeats_hg38[repeats_hg38$repClass == i, ]
    
    
    # Create logical vector for peaks to keep
    peaks_keep <- overlapsAny(peak_gr, repeats_subset)
    
    df[paste0(i,'_rep')] <- as.vector(peaks_keep)
  }
  
  
  return(df)
}

for (i in c('all', 'high_confidence_removes', 'moderate_removes', 'keep_types')) {
  tmp <- summarize_rep_func3(
    peak_gr = peak_gr,
    repeats_hg38 = repeats_hg38,
    peak_matrix = counts$Peaks,
    Category = i
  )
  
  df_peak_qc <- cbind(df_peak_qc, tmp)
  
}


tmp <- summarize_rep_func4(
  peak_gr = peak_gr,
  repeats_hg38 = repeats_hg38,
  peak_matrix = counts$Peaks
)

df_peak_qc <- cbind(df_peak_qc, tmp)


for (i in c('all', 'high_confidence_removes', 'moderate_removes', 'keep_types')) {
  df_peak_qc_peak_name <- summarize_rep_func5(
    peak_gr = peak_gr,
    repeats_hg38 = repeats_hg38,
    peak_matrix = counts$Peaks,
    df_peak_qc_peak_name = df_peak_qc_peak_name,
    Category = i
  )

}


df_peak_qc_peak_name <- summarize_rep_func6(
  peak_gr = peak_gr,
  repeats_hg38 = repeats_hg38,
  df_peak_qc_peak_name = df_peak_qc_peak_name,
  peak_matrix = counts$Peaks
)








### check peaks within MT --------------
peaks_keep <- seqnames(peak_gr) %in% 'chrMT'
peak_ratio = sum(peaks_keep)/length(peaks_keep)
df_peak_qc['peak_chrmt_ratio'] <- peak_ratio
df_peak_qc_peak_name['peak_chrMT'] <- as.vector(peaks_keep)

### check peaks within blacklisted regions --------------
peaks_keep <- overlapsAny(peak_gr, blacklist_hg38_unified)
peak_ratio = sum(peaks_keep)/length(peaks_keep)
df_peak_qc['peak_black_ratio'] <- peak_ratio
df_peak_qc_peak_name['peak_black_ratio'] <- as.vector(peaks_keep)

### check # of cells expressing corresponding peaks --------------
peak_cell_ratio = Matrix::rowSums(counts$Peaks > 0)/ncol(counts$Peaks)
df_peak_qc_peak_name['peak_cell_ratio'] <- unname(peak_cell_ratio)
df_peak_qc_peak_name['peak_cell'] <- unname(Matrix::rowSums(counts$Peaks > 0))
df_peak_qc_peak_name['peak_reads'] <- unname(Matrix::rowSums(counts$Peaks))
df_peak_qc_peak_name['peak_width'] <- width(peak_gr)

### save peak QC results for peaks ---------
df_peak_qc <- as.data.frame(t(df_peak_qc))
write.csv(
  df_peak_qc,
  file = str_c(data_save, "peak_level_peak_qc.csv"),
  row.names = TRUE
)

write.csv(
  df_peak_qc_peak_name,
  file = str_c(data_save, "peak_level_peak_qc_peak_name.csv"),
  row.names = TRUE
)





# QC using fragments--------------
# Comprehensive scATAC-seq Data Filtering in R (without Seurat/Signac)
# Implements all standard QC criteria for single-cell ATAC-seq

## get gene annotations for hg38 -------
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

ah <- AnnotationHub()
query(ah, c("RepeatMasker", "Homo sapiens"))
repeats_hg38 <- ah[["AH99003"]]
repeats_hg38 <- repeats_hg38[!grepl("\\?", repeats_hg38$repClass),]

# Classify repeats by type
high_confidence_removes <- c("Simple_repeat", "Low_complexity", "Satellite", 
                             "tRNA", "rRNA", "scRNA", "snRNA", "srpRNA")

moderate_removes <- c("LINE", "SINE", "LTR")  # More controversial

keep_types <- c("DNA", "RC")  # Often contain regulatory elements

## create a Seurat object containing the RNA and atac data ---------
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

rm(list = 'counts')
gc()
gc()




## get frag information --------
frag_summary <- CountFragments(fragments = fragpath, 
                               cells = colnames(pbmc))


cell_barcodes <- colnames(pbmc)


data <- vroom("../data/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz", 
              delim = "\t", col_names = FALSE)
data <- data[data$X4 %in% cell_barcodes, ]

tmp = data$X4
gc()
gc()

data <- GRanges(
  seqnames = data$X1,
  ranges = IRanges(start = data$X2, end = data$X3))
gc()
gc()



### overview of fragments information -------
DefaultAssay(pbmc) <- "ATAC"
frag_summary <- CountFragments(fragments = fragpath, 
                       cells = colnames(pbmc))

gc()
gc()


### check fragments within peaks --------------
df_peak_frag = overlapsAny(data, granges(pbmc[['ATAC']]))
df_peak_frag = as.data.frame(table(df_peak_frag, tmp))
df_peak_frag = df_peak_frag[df_peak_frag$df_peak_frag == TRUE,]
df_peak_frag = df_peak_frag[match(cell_barcodes, df_peak_frag$tmp),]
gc()
gc()

### check fragments within blacklisted regions ----------
df_bl_frag = overlapsAny(data, blacklist_hg38_unified)
df_bl_frag = as.data.frame(table(df_bl_frag, tmp))
df_bl_frag = df_bl_frag[df_bl_frag$df_bl_frag == TRUE,]
df_bl_frag = df_bl_frag[match(cell_barcodes, df_bl_frag$tmp),]
gc()
gc()


### check fragments within standard chromosomes ----------
chromosomes_to_keep <- c(paste0("chr", 1:22), "chrX")
df_chr_frag = seqnames(data) %in% chromosomes_to_keep
df_chr_frag = as.vector(df_chr_frag)
df_chr_frag = as.data.frame(table(df_chr_frag, tmp))
df_chr_frag = df_chr_frag[df_chr_frag$df_chr_frag == TRUE,]
df_chr_frag = df_chr_frag[match(cell_barcodes, df_chr_frag$tmp),]
gc()
gc()

### check fragments within chrMT ----------
chromosomes_to_keep <- c("chrMT")
df_chrMT_frag = seqnames(data) %in% chromosomes_to_keep
df_chrMT_frag = as.vector(df_chrMT_frag)
df_chrMT_frag = as.data.frame(table(df_chrMT_frag, tmp))
df_chrMT_frag = df_chrMT_frag[df_chrMT_frag$df_chrMT_frag == TRUE,]
df_chrMT_frag = df_chrMT_frag[match(cell_barcodes, df_chrMT_frag$tmp),] # if all output is NA, that means all values are 0
gc()
gc()

write.csv(df_chrMT_frag, str_c(data_save, "df_chrMT_frag.csv"), row.names = FALSE)
write.csv(df_chr_frag, str_c(data_save, "df_chr_frag.csv"), row.names = FALSE)
write.csv(df_peak_frag, str_c(data_save, "df_peak_frag.csv"), row.names = FALSE)
write.csv(df_bl_frag, str_c(data_save, "df_bl_frag.csv"), row.names = FALSE)

### check fragments within genomic reptitive elements (optioonal) ----------
### very computationally intensive


### check nucleosome signal and TSS enrichment -------
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)


# Integrate QC metrics into metadata --------------
# df_chrMT_frag = read.csv(str_c(data_save, "df_chrMT_frag.csv"))
# df_chr_frag = read.csv(str_c(data_save, "df_chr_frag.csv"))
# df_peak_frag = read.csv(str_c(data_save, "df_peak_frag.csv"))
# df_bl_frag = read.csv(str_c(data_save, "df_bl_frag.csv"))
# df_cell_peak = read.csv(str_c(data_save, "cell_level_peak_qc.csv"))

pbmc$fragments <- frag_summary$frequency_count
pbmc$peak_frag_ratio <- df_peak_frag$Freq/pbmc$fragments
pbmc$blacklist_frag_ratio <- df_bl_frag$Freq/pbmc$fragments
pbmc$standard_chr_frag_ratio <- df_chr_frag$Freq/pbmc$fragments
pbmc$chrMT_frag_ratio <- 0/pbmc$fragments
pbmc@meta.data <- cbind(pbmc@meta.data, df_cell_peak[,-1])

df_atac_qc_cell_lvl = pbmc@meta.data[,!grepl('(nCount_RNA)|(nFeature_RNA)', colnames(pbmc@meta.data))]

## visualization of QC metrics ----------
### cell level ----------
ggplot(df_atac_qc_cell_lvl, 
       aes(x=fragments,fill=orig.ident))+
  geom_density()+
  scale_x_log10()+
  geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle("Distribution of fragments per cell")+
  xlab("Fragments per cell")+
  ridgeplot_clean_theme()
ggsave(
  str_c(analysis_save, "fragments_per_cell.pdf"),
  width = 8, height = 6
)

ggplot(df_atac_qc_cell_lvl, 
       aes(x=nCount_ATAC,fill=orig.ident))+
  geom_density()+
  scale_x_log10()+
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


ggplot(df_atac_qc_cell_lvl, 
       aes(x=nFeature_ATAC,fill=orig.ident))+
  geom_density()+
  scale_x_log10()+
  geom_vline(xintercept = 1000, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle("Distribution of # of peaks per cell")+
  xlab("Fragments per cell")+
  ridgeplot_clean_theme()
ggsave(
  str_c(analysis_save, "peak_per_cell.pdf"),
  width = 8, height = 6
)


qc_metrics = c("peak_frag_ratio", 
               "blacklist_frag_ratio", 
               "standard_chr_frag_ratio", 
               "chrMT_frag_ratio",
               "peak_std_chr_ratio", 
               "peak_black_ratio", 
               "peak_chrmt_ratio", 
               "high_confidence_removes_rep",
               "moderate_removes_rep", 
               "keep_types_rep",
               "nucleosome_signal", 
               "TSS.enrichment")

names(qc_metrics) = c("Rario of fragments within peaks",
                      "Ratio of blacklisted fragments",
                      "Ratio of fragments within standard chromosome (chr1-22, chrX)",
                      "Ratio of fragments within chrMT ",
                      "Ratio of peaks within standard chromosomes",
                      "Ratio of peaks within blacklisted regions",
                      "Ratio of peaks within chrMT",
                      "Ratio of peaks within repetitive elements(microsatellites, tandem repeats, etc.)",
                      "Ratio of peaks within repetitive elements(LINEs, SINEs and LTRs)",
                      "Ratio of peaks within other repetitive elements(DNA transponsons and Rolling circle elements)",
                      "Nucleosome signal", 
                      "TSS enrichment")


graphics.off()
for (i in unname(qc_metrics)) {
  cut_off = switch(i,
                   "peak_frag_ratio" = 0.4,
                   "blacklist_frag_ratio" = 0.05,
                   "standard_chr_frag_ratio" = 0.95,
                   "chrMT_frag_ratio" = 0.05,
                   "peak_std_chr_ratio" = 0.95,
                   "peak_black_ratio" = 0.05,
                   "peak_chrmt_ratio" = 0.05,
                   "high_confidence_removes_rep" = 0,
                   "moderate_removes_rep" = 0,
                   "keep_types_rep" = 0,
                   "nucleosome_signal" = 2,
                   "TSS.enrichment" = 1
                   )
  
  linecolor = switch(i,
                     "peak_frag_ratio" = "blue",
                     "blacklist_frag_ratio" = "red",
                     "standard_chr_frag_ratio" = "blue",
                     "chrMT_frag_ratio" = "red",
                     "peak_std_chr_ratio" = "blue",
                     "peak_black_ratio" = "red",
                     "peak_chrmt_ratio" = "red",
                     "high_confidence_removes_rep" = "blue",
                     "moderate_removes_rep" = "blue",
                     "keep_types_rep" = "blue",
                     "nucleosome_signal" = "red",
                     "TSS.enrichment" = "red"
                     )
  
  ggplot(df_atac_qc_cell_lvl, 
         aes_string(x=i,fill='orig.ident'))+
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


qc_genomic_metrics = c("all repetitive elements",
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

names(qc_genomic_metrics) = c("all_rep", 
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

df = df_atac_qc_cell_lvl[,names(qc_genomic_metrics)]
df = tidyr::pivot_longer(df, 
                        cols = everything(), 
                        names_to = "QC_metric_rep", 
                        values_to = "Ratio")
df$QC_metric_rep <- factor(df$QC_metric_rep, 
                           levels = names(qc_genomic_metrics))



graphics.off()
ggplot(df, 
       aes_string(x="Ratio",
                  fill='QC_metric_rep'))+
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


### peak level ----------
# df_peak_qc = read.csv(str_c(data_save, "peak_level_peak_qc.csv"), row.names = 1)
# df_peak_qc$peak_cat = rownames(df_peak_qc)

peak_qc = c("Ratio of peaks within standard chromosomes",
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

names(peak_qc) = c("peak_std_chr_ratio", 
            "peak_black_ratio", 
            "peak_chrmt_ratio",
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


df_peak_qc$peak_cat = factor(df_peak_qc$peak_cat, 
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

df_peak_qc_peak_name$orig.ident = unique(pbmc@meta.data$orig.ident)
ggplot(df_peak_qc_peak_name, 
       aes_string(x="peak_cell_ratio",fill='orig.ident'))+
  geom_density()+
  geom_vline(xintercept = 0.01, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle(str_c("Distribution of ratios of cells expressing corresponding peaks"))+
  xlab("Ratio of cells per peak")+
  ridgeplot_clean_theme()
ggsave(str_c(analysis_save, "cell_per_peak_ratio_peak_lvl.pdf"),
       width = 8, height = 6)

ggplot(df_peak_qc_peak_name, 
       aes_string(x="peak_cell",fill='orig.ident'))+
  geom_density()+
  geom_vline(xintercept = 5, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle(str_c("Distribution of # of cells expressing corresponding peaks"))+
  xlab("# of cells per peak")+
  ridgeplot_clean_theme()
ggsave(str_c(analysis_save, "cell_per_peak_peak_lvl.pdf"),
       width = 8, height = 6)

ggplot(df_peak_qc_peak_name, 
       aes_string(x="peak_reads",fill='orig.ident'))+
  geom_density()+
  geom_vline(xintercept = 100, linetype = "dashed", color = "red", linewidth = 1)+
  NoLegend()+
  ggtitle(str_c("Distribution of # of reads for peaks"))+
  xlab("# of reads for peaks")+
  scale_x_log10()+
  ridgeplot_clean_theme()
ggsave(str_c(analysis_save, "peak_reads_peak_lvl.pdf"),
       width = 8, height = 6)

ggplot(df_peak_qc_peak_name, 
       aes_string(x="peak_width",fill='orig.ident'))+
  geom_density()+
  geom_vline(xintercept = 200, linetype = "dashed", color = "blue", linewidth = 1)+
  NoLegend()+
  ggtitle(str_c("Distribution of widths of peaks"))+
  xlab("peak widths")+
  ridgeplot_clean_theme()+
  scale_x_log10()
ggsave(str_c(analysis_save, "peak_width_peak_lvl.pdf"),
       width = 8, height = 6)



## filter cells and peaks based on QC metrics --------------
df_atac_qc_summary = data.frame(row.names = unique(pbmc$orig.ident),
                                nCells_before = ncol(pbmc),
                                nCells_after = NA,
                                nPeaks = nrow(pbmc[['ATAC']]),
                                nPeaks_after = NA)


keep_atacseq = df_atac_qc_cell_lvl$fragments > 1000 & 
  df_atac_qc_cell_lvl$nCount_ATAC > 1000 & 
  df_atac_qc_cell_lvl$nCount_ATAC < 100000 &
  df_atac_qc_cell_lvl$nFeature_ATAC > 1000 & 
  df_atac_qc_cell_lvl$peak_frag_ratio > 0.4 & 
  df_atac_qc_cell_lvl$blacklist_frag_ratio < 0.05 & 
  df_atac_qc_cell_lvl$standard_chr_frag_ratio > 0.95 & 
  df_atac_qc_cell_lvl$chrMT_frag_ratio < 0.05 & 
  df_atac_qc_cell_lvl$peak_std_chr_ratio > 0.95 & 
  df_atac_qc_cell_lvl$peak_black_ratio < 0.05 & 
  df_atac_qc_cell_lvl$peak_chrmt_ratio < 0.05 &
  pbmc$nucleosome_signal < 2 &
  pbmc$TSS.enrichment > 1

df_cell_atac_qc = data.frame(
  row.names = colnames(pbmc),
  keep_atacseq_cell = keep_atacseq
)


keep_atacseq_peak = df_peak_qc_peak_name$peak_cell_ratio > 0.01 & 
  df_peak_qc_peak_name$peak_cell > 5 & 
  df_peak_qc_peak_name$peak_reads > 100 & 
  df_peak_qc_peak_name$peak_std_chr == TRUE & 
  df_peak_qc_peak_name$peak_black_ratio == FALSE & 
  df_peak_qc_peak_name$peak_chrMT == FALSE

df_peak_atac_qc = data.frame(
  row.names = rownames(pbmc),
  keep_atacseq_peak = keep_atacseq_peak
)

df_atac_qc_summary['nCells_after'] <- sum(keep_atacseq)
df_atac_qc_summary['nPeaks_after'] <- sum(keep_atacseq_peak)

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

## scRNA-seq ----------
DefaultAssay(pbmc) <- "RNA"











