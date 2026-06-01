#!/usr/bin/env Rscript

# ============================================
# GREAT enrichment analysis for DAR peaks
# ============================================

library(rGREAT)
library(dplyr)

run_GO_enrichment <- function(
    bed_file,
    output_dir,
    cell_type,
    peak_type = "up",  # "up", "down"
    species = "hg38",
    rule = "basalPlusExt",
    upstream = 5.0,      # kb, TSS upstream extension (default 5kb)
    downstream = 1.0,    # kb, TSS downstream extension (default 1kb)  
    span = 1000.0,       # kb, maximum distance to associate genes (default 1Mb)
    timeout_sec = 600    # seconds, HTTP timeout (increase if slow network)
) {
  
  # Print header
  cat(rep("=", 70), "\n")
  cat("GREAT Enrichment Analysis\n")
  cat("  Cell type: ", cell_type, "\n", sep = "")
  cat("  Peak type: ", peak_type, "\n", sep = "")
  cat("  BED file:  ", basename(bed_file), "\n", sep = "")
  cat(rep("=", 70), "\n\n")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ============================================
  # Step 1: Read and validate BED file
  # ============================================
  cat("[1/5] Reading BED file...\n")
  
  if (!file.exists(bed_file)) {
    stop("Error: BED file not found: ", bed_file)
  }
  
  bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
  cat("      File has", ncol(bed_data), "columns\n")
  
  # Extract first 3 columns as BED format
  if (ncol(bed_data) < 3) {
    stop("Error: BED file must have at least 3 columns")
  }
  
  great_bed <- bed_data[, 1:3]
  colnames(great_bed) <- c("chr", "start", "end")
  
  # Ensure start < end
  great_bed <- great_bed[great_bed$start < great_bed$end, ]
  
  # Remove any rows with NA
  great_bed <- na.omit(great_bed)
  
  cat("      Valid peaks:", nrow(great_bed), "\n")
  
  if (nrow(great_bed) == 0) {
    stop("Error: No valid peaks after filtering")
  }
  
  # Save formatted BED
  formatted_bed <- file.path(output_dir, paste0(cell_type, "_", peak_type, "_peaks.bed"))
  write.table(great_bed, formatted_bed,
              sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  cat("      BED saved:", basename(formatted_bed), "\n")
  
  # ============================================
  # Step 2: Submit GREAT job
  # ============================================
  cat("\n[2/5] Submitting to GREAT server...\n")
  cat("      Parameters:\n")
  cat("        - Species:", species, "\n")
  cat("        - Rule:", rule, "\n")
  cat("        - Upstream:", upstream, "kb\n")
  cat("        - Downstream:", downstream, "kb\n")
  cat("        - Max span:", span, "kb\n")
  
  # Set timeout
  options(timeout = timeout_sec)
  
  # Submit job with retry mechanism
  max_attempts <- 3
  job <- NULL
  
  for (attempt in 1:max_attempts) {
    tryCatch({
      cat("        Attempt", attempt, "...\n")
      job <- submitGreatJob(
        bed_file = formatted_bed,
        species = species,
        rule = rule,
        adv_upstream = upstream,
        adv_downstream = downstream,
        adv_span = span
      )
      break
    }, error = function(e) {
      cat("        Attempt", attempt, "failed:", e$message, "\n")
      if (attempt == max_attempts) {
        stop("Failed to submit GREAT job after ", max_attempts, " attempts")
      }
      Sys.sleep(10)  # Wait before retry
    })
  }
  
  cat("      Job submitted successfully\n")
  
  # ============================================
  # Step 3: Get enrichment tables
  # ============================================
  cat("\n[3/5] Retrieving enrichment results...\n")
  
  # Get tables (no extra parameters for compatibility)
  enrich_table <- getEnrichmentTables(job)
  
  # Display available tables
  cat("      Available tables:", paste(names(enrich_table), collapse = ", "), "\n")
  
  # Extract GO categories (rGREAT typically returns: MF, BP, CC in this order)
  mf_results <- NULL
  bp_results <- NULL
  cc_results <- NULL
  
  # Method 1: By name
  for (i in seq_along(enrich_table)) {
    tbl_name <- names(enrich_table)[i]
    if (grepl("Biological Process", tbl_name, ignore.case = TRUE)) {
      bp_results <- enrich_table[[i]]
      cat("      Found: GO Biological Process\n")
    } else if (grepl("Molecular Function", tbl_name, ignore.case = TRUE)) {
      mf_results <- enrich_table[[i]]
      cat("      Found: GO Molecular Function\n")
    } else if (grepl("Cellular Component", tbl_name, ignore.case = TRUE)) {
      cc_results <- enrich_table[[i]]
      cat("      Found: GO Cellular Component\n")
    }
  }
  
  # Method 2: Fallback by position (rGREAT typical order: MF, BP, CC)
  if (is.null(bp_results) && length(enrich_table) >= 2) {
    cat("      Using positional fallback (MF, BP, CC order)\n")
    mf_results <- enrich_table[[1]]
    bp_results <- enrich_table[[2]]
    if (length(enrich_table) >= 3) cc_results <- enrich_table[[3]]
  }
  
  # ============================================
  # Step 4: Save results
  # ============================================
  cat("\n[4/5] Saving results...\n")
  
  # Helper function to save results
  save_results <- function(results, category, suffix) {
    if (!is.null(results) && nrow(results) > 0) {
      # Sort by FDR Q-value (ascending)
      if ("BinomFdrQ" %in% colnames(results)) {
        results_sorted <- results[order(results$BinomFdrQ, decreasing = FALSE), ]
      } else if ("HyperFdrQ" %in% colnames(results)) {
        results_sorted <- results[order(results$HyperFdrQ, decreasing = FALSE), ]
      } else {
        results_sorted <- results
      }
      
      # Save CSV
      out_file <- file.path(output_dir, 
                           paste0(cell_type, "_", peak_type, "_", category, ".csv"))
      write.csv(results_sorted, out_file, row.names = FALSE)
      cat("      ", category, ": ", nrow(results_sorted), " terms saved\n", sep = "")
      
      # Print top 5
      if (category == "GO_BP" && nrow(results_sorted) > 0) {
        cat("\n      Top 5 Biological Process terms:\n")
        top_cols <- intersect(c("name", "BinomFdrQ", "BinomP", "HyperFdrQ"), 
                             colnames(results_sorted))
        print(head(results_sorted[, top_cols, drop = FALSE], 5))
        cat("\n")
      }
      
      return(results_sorted)
    } else {
      cat("      ", category, ": No results\n", sep = "")
      return(NULL)
    }
  }
  
  bp_saved <- save_results(bp_results, "GO_BP", "BP")
  mf_saved <- save_results(mf_results, "GO_MF", "MF")
  cc_saved <- save_results(cc_results, "GO_CC", "CC")
  
  # ============================================
  # Step 5: Get region-gene associations
  # ============================================
  cat("\n[5/5] Getting region-gene associations...\n")
  
  tryCatch({
    gene_assoc <- getRegionGeneAssociations(job)
    
    if (!is.null(gene_assoc)) {
      # Convert to data frame based on type
      if (is.data.frame(gene_assoc)) {
        gene_assoc_df <- gene_assoc
      } else if (is.list(gene_assoc)) {
        # Extract relevant information
        if (all(c("annotated_genes", "distances") %in% names(gene_assoc))) {
          gene_assoc_df <- data.frame(
            region = names(gene_assoc$annotated_genes),
            gene = unlist(gene_assoc$annotated_genes),
            distance = unlist(gene_assoc$distances),
            stringsAsFactors = FALSE
          )
        } else {
          gene_assoc_df <- as.data.frame(do.call(rbind, gene_assoc))
        }
      } else {
        gene_assoc_df <- as.data.frame(gene_assoc)
      }
      
      if (nrow(gene_assoc_df) > 0) {
        out_file <- file.path(output_dir, 
                             paste0(cell_type, "_", peak_type, "_gene_associations.csv"))
        write.csv(gene_assoc_df, out_file, row.names = FALSE)
        cat("      Saved:", nrow(gene_assoc_df), "region-gene associations\n")
      } else {
        cat("      No associations found\n")
      }
    } else {
      cat("      No associations returned\n")
    }
  }, error = function(e) {
    cat("      Warning: Could not retrieve associations:", e$message, "\n")
  })
  
  # ============================================
  # Summary
  # ============================================
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("✓ GREAT analysis complete!\n")
  cat("  Output directory:", output_dir, "\n")
  cat("  Files saved:\n")
  cat("    -", cell_type, "_", peak_type, "_GO_BP.csv\n", sep = "")
  cat("    -", cell_type, "_", peak_type, "_GO_MF.csv\n", sep = "")
  cat("    -", cell_type, "_", peak_type, "_GO_CC.csv\n", sep = "")
  cat("    -", cell_type, "_", peak_type, "_gene_associations.csv\n", sep = "")
  cat(rep("=", 70), "\n", sep = "")
  
  # Return results invisibly
  invisible(list(
    BP = bp_saved,
    MF = mf_saved,
    CC = cc_saved,
    job = job
  ))
}

# ============================================
# Example usage
# ============================================

# Analyze up-regulated peaks
run_GO_enrichment(
  bed_file = "/Subset_CD4/DAR_results/CD4_Tn/up_peaks.bed",
  output_dir = "/Subset_CD4/DAR_results/CD4_Tn/GO_results",
  cell_type = "Naive_CD4",
  peak_type = "up"
)

# Analyze down-regulated peaks
run_GO_enrichment(
  bed_file = "/Subset_CD4/DAR_results/CD4_Tn/down_peaks.bed",
  output_dir = "/Subset_CD4/DAR_results/CD4_Tn/GO_results",
  cell_type = "Naive_CD4",
  peak_type = "down"
)
