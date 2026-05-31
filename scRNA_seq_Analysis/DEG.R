#!/usr/bin/env Rscript

# ============================================
# Pseudobulk Differential Expression Analysis with limma
# From .h5ad to DEG results
# ============================================
#
# This script:
#   1. Converts .h5ad to Seurat object
#   2. Creates pseudobulk counts (by sample + cell type)
#   3. Runs limma-voom for each cell type (EX vs SED)
#   4. Adjusts for sex and BMI (if provided)
#   5. Saves DEG tables and summary plots
#
# Usage:
#   Modify the USER SETTINGS section below, then run:
#   Rscript pseudobulk_limma_DEG.R
#
# ============================================

# Load required packages
suppressPackageStartupMessages({
  library(sceasy)
  library(Seurat)
  library(limma)
  library(edgeR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reticulate)
})

# ============================================
# USER SETTINGS - MODIFY THESE!
# ============================================

# Input/output paths
input_h5ad <- "path/to/your/data.h5ad"      # CHANGE: path to your .h5ad file
output_dir <- "./limma_results/"            # CHANGE: output directory

# Metadata column names (must match your .h5ad)
sample_col <- "sample_id"                   # CHANGE: sample/patient ID column
celltype_col <- "cell_type_l3"              # CHANGE: cell type column
group_col <- "exercise_group"               # CHANGE: group column (control vs treatment)
control_group <- "SED"                 # CHANGE: control group name
treatment_group <- "EX"       # CHANGE: treatment group name
covariates <- c("sex", "BMI")               # CHANGE: covariates (e.g., c("sex", "age"))

# Analysis parameters
log2fc_threshold <- 0.25                    # |log2FC| cutoff for significance
p_adj_threshold <- 0.05                     # Adjusted p-value cutoff
min_cells_per_sample <- 5                   # Minimum cells per sample for pseudobulk
min_samples_per_group <- 2                  # Minimum samples per group for analysis

# ============================================
# DO NOT MODIFY BELOW THIS LINE
# ============================================

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
deg_dir <- file.path(output_dir, "all_DEGs")
sig_dir <- file.path(output_dir, "significant_DEGs")
plot_dir <- file.path(output_dir, "plots")
dir.create(deg_dir, showWarnings = FALSE)
dir.create(sig_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

cat("=", rep("=", 60), "\n")
cat("Pseudobulk DEG Analysis with limma\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=", rep("=", 60), "\n\n")

# 1. Convert .h5ad to Seurat
cat("Step 1: Converting .h5ad to Seurat...\n")
seurat_obj <- sceasy::convertFormat(
  obj = input_h5ad,
  from = "anndata",
  to = "seurat"
)
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Genes:", nrow(seurat_obj), "\n")

# 2. Extract data
cat("\nStep 2: Creating pseudobulk counts...\n")
counts <- GetAssayData(seurat_obj, slot = "counts")
metadata <- seurat_obj@meta.data

# Check required columns
required_cols <- c(sample_col, celltype_col, group_col)
missing_cols <- required_cols[!required_cols %in% colnames(metadata)]
if(length(missing_cols) > 0) {
  stop("Missing columns in metadata: ", paste(missing_cols, collapse=", "))
}

# Create pseudobulk IDs
metadata$pseudo_id <- paste0(metadata[[sample_col]], "#", metadata[[celltype_col]])

# Aggregate counts
unique_ids <- unique(metadata$pseudo_id)
pb_counts <- matrix(0, nrow = nrow(counts), ncol = length(unique_ids))
rownames(pb_counts) <- rownames(counts)
colnames(pb_counts) <- unique_ids

for(id in unique_ids) {
  cells <- rownames(metadata)[metadata$pseudo_id == id]
  if(length(cells) >= min_cells_per_sample) {
    pb_counts[, id] <- Matrix::rowSums(counts[, cells, drop = FALSE])
  }
}
pb_counts <- pb_counts[, colSums(pb_counts) > 0]

# Pseudobulk metadata
pb_meta <- data.frame(pseudo_id = colnames(pb_counts), stringsAsFactors = FALSE)
pb_meta$sample <- sapply(strsplit(pb_meta$pseudo_id, "#"), `[`, 1)
pb_meta$cell_type <- sapply(strsplit(pb_meta$pseudo_id, "#"), `[`, 2)

# Add group and covariates
sample_info <- metadata[!duplicated(metadata[[sample_col]]), 
                        c(sample_col, group_col, covariates)]
colnames(sample_info)[1] <- "sample"
colnames(sample_info)[2] <- "group"

pb_meta <- merge(pb_meta, sample_info, by = "sample", all.x = TRUE)
pb_meta <- pb_meta[!is.na(pb_meta$group), ]
pb_counts <- pb_counts[, pb_meta$pseudo_id]

cat("  Pseudobulk samples:", ncol(pb_counts), "\n")
cat("  Cell types:", paste(unique(pb_meta$cell_type), collapse=", "), "\n")

# 3. Run limma for each cell type
cat("\nStep 3: Running limma for each cell type...\n")

all_results <- list()
summary_df <- data.frame()

for(ct in unique(pb_meta$cell_type)) {
  
  cat("\n  Analyzing:", ct, "\n")
  
  sub_meta <- pb_meta[pb_meta$cell_type == ct, ]
  sub_counts <- pb_counts[, sub_meta$pseudo_id, drop = FALSE]
  
  n_ctl <- sum(sub_meta$group == control_group)
  n_trt <- sum(sub_meta$group == treatment_group)
  cat("    Control:", n_ctl, "Treatment:", n_trt, "\n")
  
  if(n_ctl < min_samples_per_group || n_trt < min_samples_per_group) {
    cat("    SKIP: insufficient samples\n")
    next
  }
  
  # Filter low-expression genes
  keep <- rowSums(sub_counts >= 10) >= min_samples_per_group
  sub_counts <- sub_counts[keep, ]
  cat("    Genes after filtering:", nrow(sub_counts), "\n")
  
  # Normalize
  dge <- DGEList(counts = sub_counts)
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Design matrix
  sub_meta$group <- factor(sub_meta$group, levels = c(control_group, treatment_group))
  
  design_vars <- "group"
  for(cov in covariates) {
    if(cov %in% colnames(sub_meta)) {
      if(is.character(sub_meta[[cov]])) {
        sub_meta[[cov]] <- factor(sub_meta[[cov]])
      }
      design_vars <- paste0(design_vars, " + ", cov)
    }
  }
  
  design <- model.matrix(as.formula(paste0("~ ", design_vars)), data = sub_meta)
  
  # Check design matrix
  if(qr(design)$rank < ncol(design)) {
    cat("    SKIP: singular design matrix\n")
    next
  }
  
  # Voom + limma
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  # Extract contrast
  contrast_name <- paste0("group", treatment_group)
  if(!contrast_name %in% colnames(design)) {
    cat("    SKIP: contrast not found\n")
    next
  }
  
  res <- topTable(fit, coef = contrast_name, number = Inf, adjust.method = "BH")
  
  if(nrow(res) == 0) {
    cat("    SKIP: no results\n")
    next
  }
  
  # Format results
  res$gene <- rownames(res)
  res$significant <- (abs(res$logFC) >= log2fc_threshold & res$adj.P.Val < p_adj_threshold)
  res <- res[, c("gene", "logFC", "P.Value", "adj.P.Val", "significant", "AveExpr", "t", "B")]
  colnames(res) <- c("gene", "log2FC", "p_value", "p_adj", "significant", "AveExpr", "t", "B")
  
  # Save results
  all_results[[ct]] <- res
  safe_name <- gsub("/", "_", ct)
  write.csv(res, file.path(deg_dir, paste0(safe_name, "_all_genes.csv")), row.names = FALSE)
  
  sig_res <- res[res$significant, ]
  if(nrow(sig_res) > 0) {
    write.csv(sig_res, file.path(sig_dir, paste0(safe_name, "_significant.csv")), row.names = FALSE)
  }
  
  summary_df <- rbind(summary_df, data.frame(
    Cell_Type = ct,
    Control_Samples = n_ctl,
    Treatment_Samples = n_trt,
    Genes_Tested = nrow(res),
    Significant_DEGs = sum(res$significant),
    Upregulated = sum(res$significant & res$log2FC > 0),
    Downregulated = sum(res$significant & res$log2FC < 0)
  ))
  
  cat("    Significant DEGs:", sum(res$significant), "\n")
}

# 4. Save summary
write.csv(summary_df, file.path(output_dir, "DEG_summary.csv"), row.names = FALSE)
cat("\n  Summary saved to:", file.path(output_dir, "DEG_summary.csv"), "\n")

# 5. Combine all significant DEGs
if(length(all_results) > 0) {
  combined <- do.call(rbind, lapply(names(all_results), function(ct) {
    df <- all_results[[ct]]
    df[df$significant, ]
  }))
  if(nrow(combined) > 0) {
    combined$Cell_Type <- rep(names(all_results), sapply(all_results, function(x) sum(x$significant)))
    write.csv(combined, file.path(output_dir, "combined_DEGs_all_celltypes.csv"), row.names = FALSE)
    cat("  Combined DEGs saved\n")
  }
}

# 6. Generate plot
if(nrow(summary_df) > 0) {
  p <- ggplot(summary_df, aes(x = reorder(Cell_Type, Significant_DEGs), y = Significant_DEGs)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Number of DEGs per Cell Type",
         x = "Cell Type", 
         y = paste0("Number of DEGs (|log2FC| >= ", log2fc_threshold, 
                    ", p_adj < ", p_adj_threshold, ")")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(plot_dir, "DEG_summary.pdf"), p, width = 10, height = max(6, nrow(summary_df)*0.3))
  cat("  Plot saved\n")
}

# 7. Done
cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nResults saved to:", output_dir, "\n")
cat("  ├── all_DEGs/\n")
cat("  ├── significant_DEGs/\n")
cat("  ├── plots/\n")
cat("  ├── DEG_summary.csv\n")
cat("  └── combined_DEGs_all_celltypes.csv\n")
cat("=", rep("=", 60), "\n")
