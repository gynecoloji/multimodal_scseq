# Single-cell ATAC-seq Peak-Gene Analysis Guide

## Overview
This guide covers two main approaches for analyzing peak-gene relationships in single-cell ATAC-seq data:
1. **Linear regression approach** - Traditional statistical modeling
2. **Machine learning approach** - Ridge/LASSO regression with feature selection

## Table of Contents
1. [Environment Setup](#environment-setup)
2. [Data Loading and Preprocessing](#data-loading-and-preprocessing)
3. [Method 1: Linear Regression Analysis](#method-1-linear-regression-analysis)
4. [Method 2: Machine Learning Approach](#method-2-machine-learning-approach)
5. [Visualization with Gviz](#visualization-with-gviz)
6. [Python Implementation](#python-implementation)

---

## Environment Setup

### R Package Dependencies
```r
# Core single-cell packages
library(Signac)          # scATAC-seq analysis
library(Seurat)          # Single-cell analysis
library(EnsDb.Hsapiens.v86)  # Human genome annotation

# Genomics packages
library(GenomicRanges)
library(GenomicInteractions)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(Rsamtools)

# Machine learning packages
library(glmnet)          # LASSO/Ridge regression
library(randomForest)    # Random forest
library(caret)          # Model training and validation

# Visualization
library(Gviz)           # Genomic visualization
library(ggplot2)
library(patchwork)

# Data manipulation
library(Matrix)
library(dplyr)
library(data.table)

# Parallel processing
library(parallel)
library(foreach)
library(doParallel)
```

### Python Equivalent Setup
```python
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import Ridge, Lasso, ElasticNet
from sklearn.model_selection import cross_val_score, KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error
import scipy.sparse as sp
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')
```

---

## Data Loading and Preprocessing

### Step 1: Load Genomic Annotation and Data

**R Implementation:**
```r
# Set up paths
data_read <- "../data/"
analysis_save <- "../analysis/"
data_save <- "../data/"

# Load genome annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# Optional: exclude gap regions
annotation <- annotation[!(annotation$type %in% "gap")]

# Load preprocessed Seurat object
seuratObj <- readRDS(file = paste0(data_read, "seuratObj_with_clustering.rds"))
```

**Key Data Structures:**
- `seuratObj[['ATAC']]@data`: ATAC-seq peak accessibility matrix
- `seuratObj[['SCT']]@data`: RNA expression matrix (SCTransform normalized)
- `annotation`: Genomic features (genes, exons, etc.)

---

## Method 1: Linear Regression Analysis

### Concept
This method identifies peaks within a specified window around each gene and performs linear regression to model the relationship between peak accessibility and gene expression.

### Step 2A: Peak-Gene Regression Function

**R Implementation:**
```r
perform_peak_gene_regression <- function(seuratObj,
                                         gene_list,
                                         annotation,
                                         window_size = 3000,  # 3kb window
                                         min_peaks = 2) {
  regression_results <- list()
  
  # Extract genomic ranges and expression matrices
  gr_peak <- granges(seuratObj[['ATAC']])
  gr_gene <- annotation[annotation$gene_name %in% rownames(seuratObj[['SCT']]), ]
  atac_matrix <- seuratObj[['ATAC']]@data
  rna_matrix <- seuratObj[['SCT']]@data
  
  for (gene in gene_list) {
    # Find gene location
    gr_select_gene <- gr_gene[gr_gene$gene_name == gene, ]
    
    # Find overlapping peaks within window
    overlaps <- overlapsAny(gr_peak, gr_select_gene, maxgap = window_size)
    
    if (sum(overlaps) < min_peaks) {
      message(paste("Insufficient peaks for gene", gene, ". Skipping."))
      next
    }
    
    # Extract relevant data
    gr_select_peak_matrix <- atac_matrix[overlaps, ]
    gr_select_gene_matrix <- rna_matrix[rownames(seuratObj[['SCT']]) == gene, ]
    peak_names <- rownames(gr_select_peak_matrix)
    
    # Prepare regression data
    y <- as.numeric(gr_select_gene_matrix)
    X <- t(as.matrix(gr_select_peak_matrix))
    df <- data.frame(y = y, X, check.names = FALSE)
    
    # Fit linear model
    lm_fit <- lm(y ~ ., data = df)
    
    # Store results
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
```

### Step 2B: Run Analysis
```r
# Run regression analysis on variable features
linear_results <- perform_peak_gene_regression(
  seuratObj = seuratObj,
  gene_list = VariableFeatures(seuratObj[['SCT']]),
  annotation = annotation,
  window_size = 3000,
  min_peaks = 2
)

# Save results
saveRDS(linear_results, file = paste0(data_save, "LPG_linear_results.rds"))
```

---

## Method 2: Machine Learning Approach

### Concept
This method uses Ridge/LASSO regression with feature selection to predict gene expression from peak accessibility across the entire genome, with optional distance constraints.

### Step 3A: Feature Selection Function

**R Implementation:**
```r
select_relevant_peaks <- function(atac_data, rna_data, target_gene,
                                  max_peaks = 5000, method = "correlation") {
  y <- rna_data[, target_gene]
  
  if(method == "correlation") {
    cat("Calculating correlations for", ncol(atac_data), "peaks...\n")
    
    # Memory-efficient chunked processing
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
```

### Step 3B: Single Gene Analysis Function

**R Implementation:**
```r
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
```

### Step 3C: Performance Evaluation

**R Implementation:**
```r
evaluate_prediction_quality <- function(correlation, r_squared, rmse, gene_name,
                                        actual_values, n_cells, n_peaks) {
  # Correlation assessment
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
  
  # R² interpretation
  if (r_squared >= 0.5) {
    r2_quality <- "High explanatory power"
  } else if (r_squared >= 0.25) {
    r2_quality <- "Moderate explanatory power"
  } else if (r_squared >= 0.1) {
    r2_quality <- "Low explanatory power"
  } else {
    r2_quality <- "Very low explanatory power"
  }
  
  # RMSE relative to gene expression range
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
  
  # Overall scoring system
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
```

---

## Visualization with Gviz

### Step 4: Genomic Track Visualization

**Purpose:** Create comprehensive genomic visualizations showing:
- Peak locations and scores
- Gene models (exons, CDS)
- Fragment coverage by cell type
- Peak-gene interactions

**R Implementation:**
```r
# Select best gene for visualization
summary_results_select <- summary_results %>%
  filter(strong_relationship == TRUE) %>%
  arrange(desc(adj_r_squared))
gene <- summary_results_select$gene[1]

# Prepare genomic ranges
gr_peaks <- granges(seuratObj[['ATAC']])
names(gr_peaks) <- rownames(seuratObj[['ATAC']])
gr_genes <- annotation[annotation$gene_name %in% rownames(seuratObj[['SCT']]), ]
gr_select_gene <- gr_genes[gr_genes$gene_name == gene, ]

# Create tracks
Idtrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
gtrack <- GenomeAxisTrack()
ptrack <- AnnotationTrack(gr_select_peak, name = "peaks")
ftrack <- DataTrack(data = bin_heights$height,
                    start = bin_heights$bin_start,
                    end = bin_heights$bin_end, 
                    chromosome = bin_heights$chr, 
                    genome = "hg38", 
                    name = "frag\ncounts",
                    type = 'histogram',
                    col.histogram = 'pink')

# Plot comprehensive view
plotTracks(list(Idtrack, gtrack, ptrack, ftrack, grtrack, cdstrack, 
               pstrack, itrack, pitrack),
           sizes = c(1,2,2,2,2,2,3,4,4),
           from = min(start(ptrack))-5000,
           to = max(end(ptrack)) + 5000)
```

---

## Python Implementation

Now let me provide Python equivalents for the key analysis steps:

### Basic Setup and Data Loading
```python
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.linear_model import Ridge, Lasso
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

# Load data (assuming AnnData format)
adata = sc.read_h5ad("seuratObj_with_clustering.h5ad")

# Separate ATAC and RNA data
atac_data = adata.obsm['X_atac']  # Peak accessibility matrix
rna_data = adata.X  # Gene expression matrix
```

### Feature Selection Function
```python
def select_relevant_peaks(atac_data, rna_data, target_gene_idx, 
                         max_peaks=5000, method="correlation"):
    """
    Select peaks most correlated with target gene expression
    """
    y = rna_data[:, target_gene_idx]
    
    if method == "correlation":
        print(f"Calculating correlations for {atac_data.shape[1]} peaks...")
        
        # Calculate correlations in chunks for memory efficiency
        chunk_size = 10000
        n_chunks = int(np.ceil(atac_data.shape[1] / chunk_size))
        peak_cors = np.zeros(atac_data.shape[1])
        
        for i in range(n_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, atac_data.shape[1])
            
            chunk_data = atac_data[:, start_idx:end_idx]
            chunk_cors = [pearsonr(chunk_data[:, j], y)[0] 
                         for j in range(chunk_data.shape[1])]
            peak_cors[start_idx:end_idx] = chunk_cors
            
            if (i + 1) % 10 == 0:
                print(f"Processed chunk {i+1} of {n_chunks}")
        
        # Handle NaN values
        peak_cors = np.nan_to_num(peak_cors)
        
        # Select top peaks by absolute correlation
        top_indices = np.argsort(np.abs(peak_cors))[::-1][:max_peaks]
        
        print(f"Selected {len(top_indices)} peaks with correlation range: "
              f"{np.round(np.abs(peak_cors[top_indices]).min(), 3)} to "
              f"{np.round(np.abs(peak_cors[top_indices]).max(), 3)}")
    
    return {
        'selected_atac': atac_data[:, top_indices],
        'selected_indices': top_indices,
        'selection_scores': peak_cors[top_indices]
    }
```

### Single Gene Analysis Function
```python
def analyze_single_gene(gene_idx, gene_name, rna_data, atac_data, 
                       max_peaks=3000, cv_folds=5, alpha=0):
    """
    Analyze single gene using Ridge/LASSO regression
    alpha=0 for Ridge, alpha=1 for LASSO
    """
    print(f"\n=== Analyzing {gene_name} ===")
    
    # Get target gene expression
    y = rna_data[:, gene_idx]
    
    # Feature selection
    print("Step 1: Selecting relevant peaks...")
    peak_selection = select_relevant_peaks(atac_data, rna_data, gene_idx, 
                                         max_peaks=max_peaks)
    X = peak_selection['selected_atac']
    
    print("Step 2: Running regression...")
    # Standardize features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Ridge/LASSO regression
    if alpha == 0:
        model = Ridge(alpha=1.0)  # Ridge
    else:
        model = Lasso(alpha=0.1)  # LASSO
    
    # Cross-validation
    cv_scores = cross_val_score(model, X_scaled, y, cv=cv_folds, 
                               scoring='r2')
    
    # Fit final model
    model.fit(X_scaled, y)
    predictions = model.predict(X_scaled)
    
    # Calculate performance metrics
    correlation = pearsonr(y, predictions)[0]
    if np.isnan(correlation):
        correlation = 0
        r_squared = 0
    else:
        r_squared = correlation**2
    
    rmse = np.sqrt(mean_squared_error(y, predictions))
    
    # Peak importance
    coefficients = model.coef_
    peak_importance = pd.DataFrame({
        'Peak_Index': peak_selection['selected_indices'],
        'Coefficient': coefficients,
        'Abs_Coefficient': np.abs(coefficients),
        'Initial_Correlation': peak_selection['selection_scores']
    }).sort_values('Abs_Coefficient', ascending=False)
    
    print(f"Performance: R = {correlation:.3f}, "
          f"R² = {r_squared:.3f}, RMSE = {rmse:.3f}")
    
    return {
        'gene': gene_name,
        'model': model,
        'scaler': scaler,
        'predictions': predictions,
        'actual': y,
        'performance': {
            'correlation': correlation,
            'r_squared': r_squared,
            'rmse': rmse,
            'cv_scores': cv_scores,
            'n_cells': len(y),
            'n_peaks_used': X.shape[1]
        },
        'peak_importance': peak_importance
    }
```

### Performance Evaluation
```python
def evaluate_prediction_quality(result):
    """
    Evaluate prediction quality with interpretable categories
    """
    correlation = result['performance']['correlation']
    r_squared = result['performance']['r_squared']
    rmse = result['performance']['rmse']
    actual_values = result['actual']
    
    # Correlation assessment
    if correlation >= 0.7:
        corr_quality = "Excellent"
    elif correlation >= 0.5:
        corr_quality = "Good"
    elif correlation >= 0.3:
        corr_quality = "Moderate"
    elif correlation >= 0.1:
        corr_quality = "Weak"
    else:
        corr_quality = "Poor"
    
    # R² assessment
    if r_squared >= 0.5:
        r2_quality = "High explanatory power"
    elif r_squared >= 0.25:
        r2_quality = "Moderate explanatory power"
    elif r_squared >= 0.1:
        r2_quality = "Low explanatory power"
    else:
        r2_quality = "Very low explanatory power"
    
    # RMSE assessment
    gene_range = np.max(actual_values) - np.min(actual_values)
    relative_rmse = rmse / gene_range if gene_range > 0 else float('inf')
    
    if relative_rmse <= 0.1:
        rmse_quality = "Very precise"
    elif relative_rmse <= 0.2:
        rmse_quality = "Precise"
    elif relative_rmse <= 0.3:
        rmse_quality = "Moderately precise"
    else:
        rmse_quality = "Imprecise"
    
    # Overall scoring
    corr_score = 5 if correlation >= 0.7 else 4 if correlation >= 0.5 else 3 if correlation >= 0.3 else 2 if correlation >= 0.1 else 1
    r2_score = 5 if r_squared >= 0.5 else 4 if r_squared >= 0.25 else 3 if r_squared >= 0.1 else 2 if r_squared >= 0.05 else 1
    rmse_score = 5 if relative_rmse <= 0.1 else 4 if relative_rmse <= 0.2 else 3 if relative_rmse <= 0.3 else 2 if relative_rmse <= 0.4 else 1
    
    overall_score = (corr_score + r2_score + rmse_score) / 3
    
    if overall_score >= 4.5:
        overall_quality = "Excellent prediction"
    elif overall_score >= 3.5:
        overall_quality = "Good prediction"
    elif overall_score >= 2.5:
        overall_quality = "Moderate prediction"
    else:
        overall_quality = "Poor prediction"
    
    print(f"\n=== PREDICTION QUALITY ASSESSMENT ===")
    print(f"Gene: {result['gene']}")
    print(f"Correlation: {correlation:.3f} ({corr_quality})")
    print(f"R²: {r_squared:.3f} ({r2_quality})")
    print(f"Relative RMSE: {relative_rmse:.3f} ({rmse_quality})")
    print(f"Overall Score: {overall_score:.2f}/5")
    print(f"Assessment: {overall_quality}")
    
    return {
        'gene': result['gene'],
        'correlation': correlation,
        'correlation_quality': corr_quality,
        'r_squared': r_squared,
        'r2_quality': r2_quality,
        'rmse': rmse,
        'relative_rmse': relative_rmse,
        'rmse_quality': rmse_quality,
        'overall_score': overall_score,
        'overall_quality': overall_quality
    }
```

### Example Usage
```python
# Analyze a single gene
gene_idx = 0  # Index of target gene
gene_name = adata.var_names[gene_idx]

# Run Ridge regression
result = analyze_single_gene(gene_idx, gene_name, rna_data, atac_data, 
                           max_peaks=3000, alpha=0)

# Evaluate performance
assessment = evaluate_prediction_quality(result)

# View top important peaks
print(result['peak_importance'].head())
```

---

## Key Differences Between Methods

| Aspect | Linear Regression | Machine Learning |
|--------|------------------|-----------------|
| **Scope** | Local (gene vicinity) | Global (genome-wide) |
| **Window** | Fixed (e.g., 3kb) | Feature selection |
| **Regularization** | None | Ridge/LASSO |
| **Interpretability** | High | Moderate |
| **Scalability** | Good | Excellent |
| **Overfitting Risk** | Low | Moderate (with CV) |

## When to Use Each Method

- **Linear Regression**: When you want to understand local regulatory relationships and have biological hypotheses about cis-regulatory elements
- **Machine Learning**: When you want to maximize predictive power and discover global regulatory patterns

---

## Next Steps

This analysis provides a foundation for understanding peak-gene regulatory relationships. Consider extending with:

1. **Multi-gene analysis**: Scale to all genes in your dataset
2. **Comparative analysis**: Compare different cell types or conditions
3. **Pathway analysis**: Group genes by biological pathways
4. **Integration**: Combine with other omics data (ChIP-seq, Hi-C)
5. **Validation**: Experimental validation of predicted relationships