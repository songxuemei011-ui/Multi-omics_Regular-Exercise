#!/usr/bin/env Rscript

# ============================================
# DAR analysis for scATAC-seq data
# Step 1: Call peaks
# Step 2: Extract cell subsets, match cells, run DAR
# ============================================

library(ArchR)
library(GenomicRanges)
library(dplyr)
library(readr)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

# ============================================
# Step 1: Call peaks (RUN ONCE per project)
# ============================================

prepare_peak_project <- function(
    project_path,
    output_dir,
    group_col = "exercise_group",
    pathToMacs2 = "./anaconda3/envs/pbmc/bin/macs2"
) {
  
  cat(rep("=", 60), "\n")
  cat("Preparing project with peaks\n")
  cat(rep("=", 60), "\n")
  
  proj <- loadArchRProject(project_path)
  cat("Groups:\n")
  print(table(proj[[group_col]]))
  
  # Add coverages
  proj <- addGroupCoverages(proj, groupBy = group_col, force = TRUE)
  
  # Call peaks
  proj <- addReproduciblePeakSet(proj, groupBy = group_col, pathToMacs2 = pathToMacs2, force = TRUE)
  
  # Add peak matrix
  proj <- addPeakMatrix(proj, force = TRUE)
  
  # Save
  saveArchRProject(proj, outputDirectory = output_dir, overwrite = TRUE)
  
  cat("Project with peaks saved to:", output_dir, "\n")
  return(proj)
}

# ============================================
# Step 2: Extract subsets, match, and run DAR
# ============================================

run_DAR_analysis <- function(
    project_with_peaks_path,
    celltype_col,
    output_dir,
    group_col = "exercise_group"
) {
  
  cat(rep("=", 60), "\n")
  cat("Running DAR analysis\n")
  cat(rep("=", 60), "\n")
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load project (already has peaks)
  proj <- loadArchRProject(project_with_peaks_path)
  
  # Get cell subtypes
  cell_subtypes <- unique(na.omit(getCellColData(proj)[[celltype_col]]))
  cat("Cell subtypes found:", paste(cell_subtypes, collapse = ", "), "\n\n")
  
  set.seed(123)
  
  for (ct in cell_subtypes) {
    
    cat(rep("-", 50), "\n")
    cat("Processing:", ct, "\n")
    
    ct_dir <- file.path(output_dir, ct)
    dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Extract cells for this subtype
    idx <- which(getCellColData(proj)[[celltype_col]] == ct)
    if (length(idx) < 30) {
      cat("  Too few cells (", length(idx), "), skipping\n")
      next
    }
    
    proj_sub <- proj[proj$cellNames[idx], ]
    cat("  Total cells:", nCells(proj_sub), "\n")
    
    # Match cells by sex and BMI within this subtype
    cell_meta <- as.data.frame(getCellColData(proj_sub)) %>%
      mutate(
        cellName = rownames(.),
        BMI_group = cut(as.numeric(BMI),
                        breaks = c(0, 18.5, 25, 30, Inf),
                        labels = c("Underweight", "Normal", "Overweight", "Obese")),
        strata = paste(sex, BMI_group, sep = "_")
      ) %>%
      filter(!is.na(BMI_group))
    
    matched_cells <- cell_meta %>%
      group_by(strata) %>%
      mutate(
        n_SED = sum(.data[[group_col]] == "SED"),
        n_EX = sum(.data[[group_col]] == "EX"),
        n_keep = min(n_SED, n_EX)
      ) %>%
      filter(n_keep >= 10) %>%
      group_by(strata, .data[[group_col]]) %>%
      slice_sample(n = first(n_keep)) %>%
      ungroup() %>%
      pull(cellName)
    
    if (length(matched_cells) < 30) {
      cat("  Too few matched cells, skipping\n")
      next
    }
    
    proj_matched <- proj_sub[matched_cells, ]
    cat("  Matched cells:", nCells(proj_matched), "\n")
    
    # Check group balance
    group_counts <- table(getCellColData(proj_matched)[[group_col]])
    print(group_counts)
    
    if (min(group_counts) < 10) {
      cat("  Insufficient per group, skipping DAR\n")
      next
    }
    
    # DAR analysis
    peakMarkers <- getMarkerFeatures(
      ArchRProj = proj_matched,
      useMatrix = "PeakMatrix",
      groupBy = group_col,
      useGroups = "EX",
      bgdGroups = "SED",
      testMethod = "wilcoxon"
    )
    
    # Save all peaks
    all_markers <- getMarkers(peakMarkers, cutOff = "FDR <= 1 & abs(Log2FC) >= 0")
    if (!is.null(all_markers$EX) && nrow(all_markers$EX) > 0) {
      write.csv(as.data.frame(all_markers$EX),
                file.path(ct_dir, paste0(ct, "_all_peaks.csv")),
                row.names = FALSE)
      cat("  All peaks:", nrow(all_markers$EX), "\n")
    }
    
    # Save significant peaks
    sig_markers <- getMarkers(peakMarkers, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.25")
    if (!is.null(sig_markers$EX) && nrow(sig_markers$EX) > 0) {
      write.csv(as.data.frame(sig_markers$EX),
                file.path(ct_dir, paste0(ct, "_significant_peaks.csv")),
                row.names = FALSE)
      
      up <- getMarkers(peakMarkers, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")
      down <- getMarkers(peakMarkers, cutOff = "FDR <= 0.05 & Log2FC <= -0.25")
      
      n_up <- ifelse(is.null(up$EX), 0, nrow(up$EX))
      n_down <- ifelse(is.null(down$EX), 0, nrow(down$EX))
      cat("  Up:", n_up, "| Down:", n_down, "\n")
      
      # Save BED files
      if (n_up > 0) {
        write.table(
          data.frame(chrom = up$EX$seqnames,
                     start = up$EX$start - 1,
                     end = up$EX$end),
          file.path(ct_dir, "up_peaks.bed"),
          sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
        )
      }
      if (n_down > 0) {
        write.table(
          data.frame(chrom = down$EX$seqnames,
                     start = down$EX$start - 1,
                     end = down$EX$end),
          file.path(ct_dir, "down_peaks.bed"),
          sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
        )
      }
    } else {
      cat("  No significant peaks\n")
    }
    
    cat("  Done:", ct, "\n")
  }
  
  cat("\n", rep("=", 60), "\n")
  cat("DAR analysis complete!\n")
  cat("Results saved to:", output_dir, "\n")
  cat(rep("=", 60), "\n")
}

# ============================================
# Usage
# ============================================

# Step 1: Prepare project with peaks
# prepare_peak_project(
#   project_path = "/path/to/ArchRProject/",
#   output_dir = "/path/to/ArchRProject_with_peaks/",
#   pathToMacs2 = "/path/to/macs2"
# )

# Step 2: Run DAR analysis
# run_DAR_analysis(
#   project_with_peaks_path = "/path/to/ArchRProject_with_peaks/",
#   celltype_col = "celltype",
#   output_dir = "/path/to/DAR_results/"
# )
