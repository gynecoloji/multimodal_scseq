# load required packages ---------
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


atac_matrix <- seuratObj[["ATAC"]]@data
rna_matrix <- seuratObj[["SCT"]]@data[rownames(seuratObj[["SCT"]]) %in% VariableFeatures(seuratObj[["SCT"]]), ] # nolint: line_length_linter.

atac_matrix <- t(atac_matrix) # now: cell*feature amtrix
rna_matrix <- t(rna_matrix)

dim(atac_matrix)
dim(rna_matrix)

# Ridge regression -------------
# Feature selection function for your large dataset
select_relevant_peaks <- function(atac_data, rna_data, target_gene,
                                  max_peaks = 5000, method = "correlation") {

  y <- rna_data[, target_gene]
  
  if(method == "correlation") {
    # Calculate correlation for all peaks (memory efficient)
    cat("Calculating correlations for", ncol(atac_data), "peaks...\n")
    
    # Use chunked processing for memory efficiency
    chunk_size <- 10000
    n_chunks <- ceiling(ncol(atac_data) / chunk_size)
    
    peak_cors <- numeric(ncol(atac_data))
    
    for(i in 1:n_chunks) {
      start_idx <- (i-1) * chunk_size + 1
      end_idx <- min(i * chunk_size, ncol(atac_data))
      
      chunk_cors <- apply(atac_data[, start_idx:end_idx], 2, 
                         function(peak) cor(peak, y, use = "complete.obs"))
      peak_cors[start_idx:end_idx] <- chunk_cors
      
      if(i %% 10 == 0) cat("Processed chunk", i, "of", n_chunks, "\n")
    }
    
    # Select top peaks by absolute correlation
    top_indices <- order(abs(peak_cors), decreasing = TRUE)[1:min(max_peaks, length(peak_cors))]
    
    cat("Selected", length(top_indices), "peaks with correlation range:", 
        round(range(abs(peak_cors[top_indices])), 3), "\n")
    
  }
  
  return(list(
    selected_atac = atac_data[, top_indices],
    selected_indices = top_indices,
    selection_scores = peak_cors[top_indices]
  ))
}

# Comprehensive single gene analysis
analyze_single_gene <- function(gene_name, rna_data, atac_data,
                                max_peaks = 3000, cv_folds = 5) {

  cat("\n=== Analyzing", gene_name, "===\n")

  # Check if gene exists
  if (!(gene_name %in% colnames(rna_data))) {
    stop("Gene not found in data")
  }
  
  # Get target gene expression
  y <- rna_data[, gene_name]
  
  # Feature selection
  cat("Step 1: Selecting relevant peaks...\n")
  peak_selection <- select_relevant_peaks(atac_data, rna_data, gene_name, 
                                        max_peaks = max_peaks)
  X <- peak_selection$selected_atac


  cat("Step 2: Running Ridge regression...\n")
  # Ridge regression with cross-validation
  cv_ridge <- cv.glmnet(X, y, alpha = 0, nfolds = cv_folds, 
                       type.measure = "mse", standardize = TRUE)
  
  # Make predictions
  predictions <- predict(cv_ridge, s = "lambda.min", newx = X)
  
  # Calculate performance metrics
  correlation <- cor(y, predictions)
  if (is.na(correlation)) {
    correlation <- 0
    r_squared <- 0
  } else {
    r_squared <- correlation^2
  }
  rmse <- sqrt(mean((y - predictions)^2))

  # Get important peaks
  coefficients <- coef(cv_ridge, s = "lambda.min")[-1]  # Remove intercept
  peak_importance <- data.frame(
    Peak_Index = peak_selection$selected_indices,
    Peak_Name = colnames(X),
    Coefficient = as.numeric(coefficients),
    Abs_Coefficient = abs(as.numeric(coefficients)),
    Initial_Correlation = peak_selection$selection_scores
  )
  peak_importance <- peak_importance[order(peak_importance$Abs_Coefficient, decreasing = TRUE), ]

  cat("Performance: R =", round(correlation, 3), 
      "R² =", round(r_squared, 3), 
      "RMSE =", round(rmse, 3), "\n")

  return(list(
    gene = gene_name,
    model = cv_ridge,
    predictions = as.numeric(predictions),
    actual = y,
    performance = list(
      correlation = correlation,
      r_squared = r_squared,
      rmse = rmse,
      n_cells = length(y),
      n_peaks_used = ncol(X),
      lambda_min = cv_ridge$lambda.min
    ),
    peak_importance = peak_importance,
    selected_peaks = colnames(X)
  ))
}

# Define performance evaluation framework
evaluate_prediction_quality <- function(correlation, r_squared, rmse, gene_name,
                                        actual_values, n_cells, n_peaks) {

  # Performance categories for single-cell RNA-ATAC relationships
  quality_assessment <- list()

  # 1. Correlation-based assessment
  if (correlation >= 0.7) {
    corr_quality <- "Excellent"
  } else if (correlation >= 0.5) {
    corr_quality <- "Good"
  } else if (correlation >= 0.3) {
    corr_quality <- "Moderate"
  } else if (correlation >= 0.1) {
    corr_quality <- "Weak"
  } else {
    corr_quality <- "Poor"
  }

  # 2. R² interpretation
  if (r_squared >= 0.5) {
    r2_quality <- "High explanatory power"
  } else if (r_squared >= 0.25) {
    r2_quality <- "Moderate explanatory power"
  } else if (r_squared >= 0.1) {
    r2_quality <- "Low explanatory power"
  } else {
    r2_quality <- "Very low explanatory power"
  }

  # 3. RMSE relative to gene expression range
  gene_range <- max(actual_values) - min(actual_values)
  relative_rmse <- rmse / gene_range

  if (relative_rmse <= 0.1) {
    rmse_quality <- "Very precise"
  } else if (relative_rmse <= 0.2) {
    rmse_quality <- "Precise"
  } else if (relative_rmse <= 0.3) {
    rmse_quality <- "Moderately precise"
  } else {
    rmse_quality <- "Imprecise"
  }

  # 4. Overall assessment
  # Convert to numerical scores
  corr_score <- case_when(
    correlation >= 0.7 ~ 5,
    correlation >= 0.5 ~ 4,
    correlation >= 0.3 ~ 3,
    correlation >= 0.1 ~ 2,
    TRUE ~ 1
  )

  r2_score <- case_when(
    r_squared >= 0.5 ~ 5,
    r_squared >= 0.25 ~ 4,
    r_squared >= 0.1 ~ 3,
    r_squared >= 0.05 ~ 2,
    TRUE ~ 1
  )

  rmse_score <- case_when(
    relative_rmse <= 0.1 ~ 5,
    relative_rmse <= 0.2 ~ 4,
    relative_rmse <= 0.3 ~ 3,
    relative_rmse <= 0.4 ~ 2,
    TRUE ~ 1
  )

  overall_score <- (corr_score + r2_score + rmse_score) / 3

  if (overall_score >= 4.5) {
    overall_quality <- "Excellent prediction"
  } else if (overall_score >= 3.5) {
    overall_quality <- "Good prediction"
  } else if (overall_score >= 2.5) {
    overall_quality <- "Moderate prediction"
  } else {
    overall_quality <- "Poor prediction"
  }

  return(list(
    gene = gene_name,
    correlation = correlation,
    correlation_quality = corr_quality,
    r_squared = r_squared,
    r2_quality = r2_quality,
    rmse = rmse,
    relative_rmse = relative_rmse,
    rmse_quality = rmse_quality,
    overall_score = overall_score,
    overall_quality = overall_quality,
    sample_size = n_cells,
    features_used = n_peaks
  ))
}

# Apply to your results
assess_gene_quality <- function(single_result) {
  evaluation <- evaluate_prediction_quality(
    correlation = single_result$performance$correlation,
    r_squared = single_result$performance$r_squared,
    rmse = single_result$performance$rmse,
    gene_name = single_result$gene,
    actual_values = single_result$actual,
    n_cells = single_result$performance$n_cells,
    n_peaks = single_result$performance$n_peaks_used
  )
  
  # Print detailed assessment
  cat("\n=== PREDICTION QUALITY ASSESSMENT ===\n")
  cat("Gene:", evaluation$gene, "\n")
  cat("Correlation:", round(evaluation$correlation, 3), "(", evaluation$correlation_quality, ")\n")
  cat("R²:", round(evaluation$r_squared, 3), "(", evaluation$r2_quality, ")\n")
  cat("Relative RMSE:", round(evaluation$relative_rmse, 3), "(", evaluation$rmse_quality, ")\n")
  cat("Overall Score:", round(evaluation$overall_score, 2), "/5\n")
  cat("Assessment:", evaluation$overall_quality, "\n")
  
  return(evaluation)
}

## analyze single gene using all peaks ----------
test_gene <- "HES4"  # Pick first gene for testing

single_result <- analyze_single_gene(test_gene, rna_data = rna_matrix,
                                     atac_data = atac_matrix,
                                     max_peaks = 3000)

gene_assessment <- assess_gene_quality(single_result)

peak_importance_df <- single_result$peak_importance



## analyze multiple genes with distance into consideration ----------
## match genes to peaks based on regions -----

single_result_ls <- list()
gene_assessment_ls <- list()
peak_gr <- granges(seuratObj[["ATAC"]])
gap_dist <- 500000
min_peaks <- 4

for (i in colnames(rna_matrix)){
  gene_gr <- annotation[annotation$gene_name %in% i]

  peak_gr_indices <- overlapsAny(peak_gr, gene_gr, maxgap = gap_dist)
  atac_matrix_select <- atac_matrix[, peak_gr_indices]

  if (sum(peak_gr_indices) <= min_peaks) next

  single_result <- analyze_single_gene(i, rna_data = rna_matrix,
                                       atac_data = atac_matrix_select,
                                       max_peaks = 3000)

  single_result_ls[[i]] <- single_result
  gene_assessment <- assess_gene_quality(single_result)
  gene_assessment_ls[[i]] <- gene_assessment

}

# Set up parallel backend
n_cores <- parallel::detectCores() - 4
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Shared data
peak_gr <- granges(seuratObj[["ATAC"]])
gap_dist <- 500000
min_peaks <- 4

clusterExport(cl,
              varlist = c("peak_gr", "gap_dist", "min_peaks"),
              envir = environment())

# Run foreach over genes
results <- foreach(gene_name = colnames(rna_matrix),
                   .packages = c("GenomicRanges", "Matrix", "glmnet", "Seurat", "dplyr")) %dopar% { # nolint: line_length_linter.
  gene_gr <- annotation[annotation$gene_name %in% gene_name]
  peak_gr_indices <- overlapsAny(peak_gr, gene_gr, maxgap = gap_dist)

  if (sum(peak_gr_indices) <= min_peaks) return(NULL)

  atac_matrix_select <- atac_matrix[, peak_gr_indices, drop = FALSE]

  single_result <- analyze_single_gene(gene_name, rna_data = rna_matrix,
                                       atac_data = atac_matrix_select,
                                       max_peaks = 3000)

  gene_assessment <- assess_gene_quality(single_result)

  list(gene = gene_name,
       result = single_result,
       assessment = gene_assessment)
}

stopCluster(cl)

results$gene_name
assessment_df <- do.call(rbind, lapply(results, function(res) {
  if (is.null(res)) return(NULL)
  data.frame(
    gene = res$gene,
    correlation = as.numeric(res$assessment$correlation),
    correlation_quality = res$assessment$correlation_quality,
    r_squared = as.numeric(res$assessment$r_squared),
    r2_quality = res$assessment$r2_quality,
    rmse = as.numeric(res$assessment$rmse),
    relative_rmse = as.numeric(res$assessment$relative_rmse),
    rmse_quality = res$assessment$rmse_quality,
    overall_score = as.numeric(res$assessment$overall_score),
    overall_quality = res$assessment$overall_quality,
    sample_size = as.integer(res$assessment$sample_size),
    features_used = as.integer(res$assessment$features_used),
    stringsAsFactors = FALSE
  )
}))

saveRDS(results, paste0(data_save, "results_Ridge.rds"))
write.csv(assessment_df, paste0(data_save, "assessment_df_Ridge.csv"))


select_relevant_peaks_parallel <- function(atac_data, rna_data, target_gene,
                                           max_peaks = 5000,
                                           method = "correlation",
                                           n_cores = NULL,
                                           chunk_size = 500) {

  if (is.null(n_cores)) n_cores <- detectCores() - 1

  y <- rna_data[, target_gene]

  if (method == "correlation") {
    cat("Calculating correlations for", ncol(atac_data), "peaks using", n_cores, "cores...\n") # nolint: line_length_linter.


    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    clusterExport(cl, varlist = c("atac_data", "y"), envir = environment())

    n_chunks <- ceiling(ncol(atac_data) / chunk_size)
    chunk_indices <- split(seq_len(ncol(atac_data)),
                           rep(1:n_chunks, each = chunk_size,
                               length.out = ncol(atac_data)))

    peak_cors <- foreach(chunk_idx = chunk_indices,
                         .combine = "c") %dopar% {
      library(Matrix)
      # Extract chunk data
      chunk_data <- atac_data[, chunk_idx, drop = FALSE]
      chunk_cors <- sapply(1:ncol(chunk_data),
                           function(i) {cor(chunk_data[,i], y)})
      names(chunk_cors) <- colnames(chunk_data)
      chunk_cors
    }

    stopCluster(cl)

    # Select top peaks by absolute correlation
    top_indices <- order(abs(peak_cors), decreasing = TRUE)[1:min(max_peaks, nrow(peak_cors))] # nolint: line_length_linter.

    cat("Selected", length(top_indices), "peaks with absolute correlation range:",
        round(range(abs(peak_cors[top_indices])), 3), "\n")
  }

  return(list(
    selected_atac = atac_data[, top_indices],
    selected_indices = top_indices,
    selection_scores = peak_cors[top_indices]
  ))
}

test_gene <- "HES4"
res <- select_relevant_peaks_parallel(atac_matrix,
                                      rna_matrix,
                                      test_gene,
                                      n_cores = 12)

# Lasso regression -------------
# Feature selection function for your large dataset
select_relevant_peaks <- function(atac_data, rna_data, target_gene,
                                  max_peaks = 5000, method = "correlation") {

  y <- rna_data[, target_gene]

  if (method == "correlation") {
    # Calculate correlation for all peaks (memory efficient)
    cat("Calculating correlations for", ncol(atac_data), "peaks...\n")

    # Use chunked processing for memory efficiency
    chunk_size <- 10000
    n_chunks <- ceiling(ncol(atac_data) / chunk_size)

    peak_cors <- numeric(ncol(atac_data))

    for (i in 1:n_chunks) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, ncol(atac_data))

      chunk_cors <- apply(atac_data[, start_idx:end_idx], 2,
                          function(peak) cor(peak, y, use = "complete.obs"))
      peak_cors[start_idx:end_idx] <- chunk_cors

      if (i %% 10 == 0) cat("Processed chunk", i, "of", n_chunks, "\n")
    }

    # Select top peaks by absolute correlation
    top_indices <- order(abs(peak_cors), decreasing = TRUE)[1:min(max_peaks, length(peak_cors))] # nolint: seq_linter

    cat("Selected", length(top_indices), "peaks with correlation range:",
        round(range(abs(peak_cors[top_indices])), 3), "\n")

  }

  return(list(
    selected_atac = atac_data[, top_indices],
    selected_indices = top_indices,
    selection_scores = peak_cors[top_indices]
  ))
}

# Comprehensive single gene analysis
analyze_single_gene <- function(gene_name, rna_data, atac_data,
                                max_peaks = 3000, cv_folds = 5) {

  cat("\n=== Analyzing", gene_name, "===\n")

  # Check if gene exists
  if (!(gene_name %in% colnames(rna_data))) {
    stop("Gene not found in data")
  }

  # Get target gene expression
  y <- rna_data[, gene_name]

  # Feature selection
  cat("Step 1: Selecting relevant peaks...\n")
  peak_selection <- select_relevant_peaks(atac_data, rna_data, gene_name,
                                          max_peaks = max_peaks)
  X <- peak_selection$selected_atac  # nolint


  cat("Step 2: Running Ridge regression...\n")
  # Ridge regression with cross-validation
  cv_lasso <- cv.glmnet(X, y, alpha = 1, nfolds = cv_folds,
                        type.measure = "mse", standardize = TRUE)

  # Make predictions
  predictions <- predict(cv_lasso, s = "lambda.min", newx = X)

  # Calculate performance metrics
  correlation <- cor(y, predictions)
  if (is.na(correlation)) {
    correlation <- 0
    r_squared <- 0
  } else {
    r_squared <- correlation^2
  }
  rmse <- sqrt(mean((y - predictions)^2))

  # Get important peaks
  coefficients <- coef(cv_lasso, s = "lambda.min")[-1]  # Remove intercept
  peak_importance <- data.frame(
    Peak_Index = peak_selection$selected_indices,
    Peak_Name = colnames(X),
    Coefficient = as.numeric(coefficients),
    Abs_Coefficient = abs(as.numeric(coefficients)),
    Initial_Correlation = peak_selection$selection_scores
  )
  peak_importance <- peak_importance[order(peak_importance$Abs_Coefficient, decreasing = TRUE), ] # nolint

  cat("Performance: R =", round(correlation, 3),
      "R² =", round(r_squared, 3),
      "RMSE =", round(rmse, 3), "\n")

  return(list(
    gene = gene_name,
    model = cv_lasso,
    predictions = as.numeric(predictions),
    actual = y,
    performance = list(
      correlation = correlation,
      r_squared = r_squared,
      rmse = rmse,
      n_cells = length(y),
      n_peaks_used = ncol(X),
      lambda_min = cv_lasso$lambda.min
    ),
    peak_importance = peak_importance,
    selected_peaks = colnames(X)
  ))
}

# Define performance evaluation framework
evaluate_prediction_quality <- function(correlation, r_squared, rmse, gene_name,
                                        actual_values, n_cells, n_peaks) {
  # 1. Correlation-based assessment
  if (correlation >= 0.7) {
    corr_quality <- "Excellent"
  } else if (correlation >= 0.5) {
    corr_quality <- "Good"
  } else if (correlation >= 0.3) {
    corr_quality <- "Moderate"
  } else if (correlation >= 0.1) {
    corr_quality <- "Weak"
  } else {
    corr_quality <- "Poor"
  }

  # 2. R² interpretation
  if (r_squared >= 0.5) {
    r2_quality <- "High explanatory power"
  } else if (r_squared >= 0.25) {
    r2_quality <- "Moderate explanatory power"
  } else if (r_squared >= 0.1) {
    r2_quality <- "Low explanatory power"
  } else {
    r2_quality <- "Very low explanatory power"
  }

  # 3. RMSE relative to gene expression range
  gene_range <- max(actual_values) - min(actual_values)
  relative_rmse <- rmse / gene_range

  if (relative_rmse <= 0.1) {
    rmse_quality <- "Very precise"
  } else if (relative_rmse <= 0.2) {
    rmse_quality <- "Precise"
  } else if (relative_rmse <= 0.3) {
    rmse_quality <- "Moderately precise"
  } else {
    rmse_quality <- "Imprecise"
  }

  # 4. Overall assessment
  # Convert to numerical scores
  corr_score <- case_when(
    correlation >= 0.7 ~ 5,
    correlation >= 0.5 ~ 4,
    correlation >= 0.3 ~ 3,
    correlation >= 0.1 ~ 2,
    TRUE ~ 1
  )

  r2_score <- case_when(
    r_squared >= 0.5 ~ 5,
    r_squared >= 0.25 ~ 4,
    r_squared >= 0.1 ~ 3,
    r_squared >= 0.05 ~ 2,
    TRUE ~ 1
  )

  rmse_score <- case_when(
    relative_rmse <= 0.1 ~ 5,
    relative_rmse <= 0.2 ~ 4,
    relative_rmse <= 0.3 ~ 3,
    relative_rmse <= 0.4 ~ 2,
    TRUE ~ 1
  )

  overall_score <- (corr_score + r2_score + rmse_score) / 3

  if (overall_score >= 4.5) {
    overall_quality <- "Excellent prediction"
  } else if (overall_score >= 3.5) {
    overall_quality <- "Good prediction"
  } else if (overall_score >= 2.5) {
    overall_quality <- "Moderate prediction"
  } else {
    overall_quality <- "Poor prediction"
  }

  return(list( # nolint: return_linter.
    gene = gene_name,
    correlation = correlation,
    correlation_quality = corr_quality,
    r_squared = r_squared,
    r2_quality = r2_quality,
    rmse = rmse,
    relative_rmse = relative_rmse,
    rmse_quality = rmse_quality,
    overall_score = overall_score,
    overall_quality = overall_quality,
    sample_size = n_cells,
    features_used = n_peaks
  ))
}

# Apply to your results
assess_gene_quality <- function(single_result) {
  evaluation <- evaluate_prediction_quality(
    correlation = single_result$performance$correlation,
    r_squared = single_result$performance$r_squared,
    rmse = single_result$performance$rmse,
    gene_name = single_result$gene,
    actual_values = single_result$actual,
    n_cells = single_result$performance$n_cells,
    n_peaks = single_result$performance$n_peaks_used
  )

  # Print detailed assessment
  cat("\n=== PREDICTION QUALITY ASSESSMENT ===\n")
  cat("Gene:", evaluation$gene, "\n")
  cat("Correlation:", round(evaluation$correlation, 3), "(", evaluation$correlation_quality, ")\n") # nolint
  cat("R²:", round(evaluation$r_squared, 3), "(", evaluation$r2_quality, ")\n")
  cat("Relative RMSE:", round(evaluation$relative_rmse, 3), "(", evaluation$rmse_quality, ")\n") # nolint
  cat("Overall Score:", round(evaluation$overall_score, 2), "/5\n")
  cat("Assessment:", evaluation$overall_quality, "\n")

  return(evaluation) # nolint: return_linter.
}

## analyze single gene using all peaks ----------
test_gene <- "HES4"  # Pick first gene for testing

single_result <- analyze_single_gene(test_gene, rna_data = rna_matrix,
                                     atac_data = atac_matrix,
                                     max_peaks = 3000)

gene_assessment <- assess_gene_quality(single_result)

peak_importance_df <- single_result$peak_importance



## analyze multiple genes with distance into consideration ----------
## match genes to peaks based on regions -----

single_result_ls <- list()
gene_assessment_ls <- list()
peak_gr <- granges(seuratObj[["ATAC"]])
gap_dist <- 500000
min_peaks <- 4

for (i in colnames(rna_matrix)){
  gene_gr <- annotation[annotation$gene_name %in% i]

  peak_gr_indices <- overlapsAny(peak_gr, gene_gr, maxgap = gap_dist)
  atac_matrix_select <- atac_matrix[, peak_gr_indices]

  if (sum(peak_gr_indices) <= min_peaks) next

  single_result <- analyze_single_gene(i, rna_data = rna_matrix,
                                       atac_data = atac_matrix_select,
                                       max_peaks = 3000)

  single_result_ls[[i]] <- single_result
  gene_assessment <- assess_gene_quality(single_result)
  gene_assessment_ls[[i]] <- gene_assessment

}

# Set up parallel backend
n_cores <- parallel::detectCores() - 4
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Shared data
peak_gr <- granges(seuratObj[["ATAC"]])
gap_dist <- 500000
min_peaks <- 4

clusterExport(cl,
              varlist = c("peak_gr", "gap_dist", "min_peaks"),
              envir = environment())

# Run foreach over genes
results <- foreach(gene_name = colnames(rna_matrix),
                   .packages = c("GenomicRanges", "Matrix", "glmnet", "Seurat", "dplyr")) %dopar% { # nolint: line_length_linter.
  gene_gr <- annotation[annotation$gene_name %in% gene_name]
  peak_gr_indices <- overlapsAny(peak_gr, gene_gr, maxgap = gap_dist)

  if (sum(peak_gr_indices) <= min_peaks) return(NULL)

  atac_matrix_select <- atac_matrix[, peak_gr_indices, drop = FALSE]

  single_result <- analyze_single_gene(gene_name, rna_data = rna_matrix,
                                       atac_data = atac_matrix_select,
                                       max_peaks = 3000)

  gene_assessment <- assess_gene_quality(single_result)

  list(gene = gene_name,
       result = single_result,
       assessment = gene_assessment)
}

stopCluster(cl)

results$gene_name
assessment_df <- do.call(rbind, lapply(results, function(res) {
  if (is.null(res)) return(NULL)
  data.frame(
    gene = res$gene,
    correlation = as.numeric(res$assessment$correlation),
    correlation_quality = res$assessment$correlation_quality,
    r_squared = as.numeric(res$assessment$r_squared),
    r2_quality = res$assessment$r2_quality,
    rmse = as.numeric(res$assessment$rmse),
    relative_rmse = as.numeric(res$assessment$relative_rmse),
    rmse_quality = res$assessment$rmse_quality,
    overall_score = as.numeric(res$assessment$overall_score),
    overall_quality = res$assessment$overall_quality,
    sample_size = as.integer(res$assessment$sample_size),
    features_used = as.integer(res$assessment$features_used),
    stringsAsFactors = FALSE
  )
}))

saveRDS(results, paste0(data_save, "results_Ridge.rds"))
write.csv(assessment_df, paste0(data_save, "assessment_df_Ridge.csv"))


# elastic ---------------
