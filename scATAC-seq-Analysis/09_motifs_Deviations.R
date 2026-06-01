#!/usr/bin/env Rscript

# ============================================
# Motif analysis for specific cell subsets
# ============================================

library(ArchR)
library(ggplot2)
library(ggrepel)
library(dplyr)

addArchRThreads(threads = 1)
addArchRGenome("hg38")

# ============================================
# Helper function: check if not in
# ============================================

`%ni%` <- function(x, table) {
  match(x, table, nomatch = 0L) == 0L
}

# ============================================
# Main motif analysis function
# ============================================

run_motif_analysis <- function(
    project_path,
    cell_types,
    output_dir,
    group_col = "exercise_group",
    use_groups = "EX",
    bgd_groups = "SED"
) {
  
  cat(rep("=", 60), "\n")
  cat("Motif Analysis\n")
  cat(rep("=", 60), "\n\n")
  
  # Load project
  cat("Step 1: Loading project...\n")
  proj <- loadArchRProject(project_path)
  
  # Add motif annotations if not already present
  cat("Step 2: Adding motif annotations...\n")
  if ("Motif" %ni% names(proj@peakAnnotation)) {
    proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force = TRUE)
  }
  
  # Add background peaks and deviation matrix
  cat("Step 3: Adding background peaks and deviation matrix...\n")
  proj <- addBgdPeaks(proj, force = TRUE)
  proj <- addDeviationsMatrix(
    ArchRProj = proj,
    peakAnnotation = "Motif",
    force = TRUE
  )
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Process each cell type
  for (ct in cell_types) {
    
    cat("\n", rep("-", 50), "\n")
    cat("Processing:", ct, "\n")
    
    ct_dir <- file.path(output_dir, ct)
    dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Extract cells for this cell type
    idx <- which(proj[[paste0(ct, "_celltype")]] == ct)
    if (length(idx) < 30) {
      cat("  Too few cells (", length(idx), "), skipping\n")
      next
    }
    
    proj_sub <- proj[proj$cellNames[idx], ]
    cat("  Total cells:", length(proj_sub$cellNames), "\n")
    
    # ============================================
    # Match cells by sex
    # ============================================
    
    cat("  Matching cells by sex...\n")
    
    cell_meta <- as.data.frame(getCellColData(proj_sub)) %>%
      mutate(
        cellName = rownames(.),
        sex = as.character(sex)
      ) %>%
      filter(!is.na(sex))
    
    set.seed(123)
    matched_cells <- cell_meta %>%
      group_by(sex, .data[[group_col]]) %>%
      mutate(n_group = n()) %>%
      group_by(sex) %>%
      mutate(n_match = min(n_group)) %>%
      ungroup() %>%
      group_by(sex, .data[[group_col]]) %>%
      slice_sample(n = first(n_match)) %>%
      ungroup() %>%
      pull(cellName)
    
    if (length(matched_cells) >= 30) {
      proj_sub <- proj_sub[matched_cells, ]
      cat("  Cells after matching:", length(proj_sub$cellNames), "\n")
      
      # Show sex balance after matching
      meta_check <- as.data.frame(getCellColData(proj_sub))
      cat("  Sex balance after matching:\n")
      print(table(meta_check[[group_col]], meta_check$sex))
    } else {
      cat("  Using all cells (matching not possible)\n")
    }
    
    # ============================================
    # Get marker features for motifs
    # ============================================
    
    cat("  Computing motif deviations...\n")
    
    markerTest <- getMarkerFeatures(
      ArchRProj = proj_sub,
      useMatrix = "MotifMatrix",
      groupBy = group_col,
      testMethod = "wilcoxon",
      useSeqnames = "z",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = use_groups,
      bgdGroups = bgd_groups
    )
    
    # Get all markers
    markerList <- getMarkers(markerTest, cutOff = "FDR <= 1")
    df <- markerList[[use_groups]]
    
    if (is.null(df) || nrow(df) == 0) {
      cat("  No motifs found\n")
      next
    }
    
    # Convert to data.frame if needed
    if (!inherits(df, "data.frame")) {
      df <- as.data.frame(df)
    }
    
    # Remove suffix from motif names (e.g., "RUNX3_123" -> "RUNX3")
    df$motif_name <- gsub("_.*", "", df$name)
    
    # Calculate -log10 FDR
    df$neg_log10_FDR <- -log10(df$FDR)
    
    # Define significance groups
    mean_diff_thresh <- 0.05
    fdr_thresh <- 0.05
    
    df$direction <- with(df, ifelse(
      abs(MeanDiff) >= mean_diff_thresh & FDR < fdr_thresh,
      ifelse(MeanDiff > 0, "Upregulated (Active)", "Downregulated (Inactive)"),
      "Not Significant"
    ))
    
    # Save all results
    write.csv(df, file.path(ct_dir, paste0(ct, "_all_motifs.csv")), row.names = FALSE)
    
    # Save significant motifs
    sig_df <- df[df$direction != "Not Significant", ]
    if (nrow(sig_df) > 0) {
      write.csv(sig_df, file.path(ct_dir, paste0(ct, "_significant_motifs.csv")), row.names = FALSE)
    }
    
    # ============================================
    # Volcano plot
    # ============================================
    
    color_pal <- c(
      "Upregulated (Active)" = "#BCDF7A",
      "Downregulated (Inactive)" = "#7EBFC9",
      "Not Significant" = "grey70"
    )
    
    p <- ggplot(df, aes(x = MeanDiff, y = neg_log10_FDR, color = direction)) +
      geom_point(data = subset(df, direction == "Not Significant"), alpha = 0.4, size = 2.5) +
      geom_point(data = subset(df, direction != "Not Significant"), alpha = 0.6, size = 2.5) +
      scale_color_manual(values = color_pal) +
      labs(x = "Mean Difference", y = "-log10(FDR)", title = paste(ct, "- Motif Volcano Plot")) +
      theme_classic() +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    # Label top motifs
    top_up <- head(subset(df, direction == "Upregulated (Active)"), 10)
    top_down <- head(subset(df, direction == "Downregulated (Inactive)"), 10)
    
    if (nrow(top_up) > 0 || nrow(top_down) > 0) {
      p <- p + geom_text_repel(
        data = rbind(top_up, top_down),
        aes(label = motif_name),
        color = "black",
        size = 3.5,
        box.padding = 0.3,
        point.padding = 0.3,
        max.overlaps = Inf,
        segment.color = NA,
        min.segment.length = 0,
        force = 1,
        nudge_x = 0.1,
        nudge_y = 0.1
      )
    }
    
    # Save plot
    pdf(file.path(ct_dir, paste0(ct, "_motif_volcano.pdf")), width = 8, height = 6)
    print(p)
    dev.off()
    
    cat("  Motif results saved to:", ct_dir, "\n")
  }
  
  cat("\n", rep("=", 60), "\n")
  cat("Motif analysis complete!\n")
  cat("Results saved to:", output_dir, "\n")
  cat(rep("=", 60), "\n")
}

# ============================================
# Usage
# ============================================

# run_motif_analysis(
#   project_path = "/path/to/Subset_CD8/",
#   cell_types = c("CD8_CTL", "CD8_Tn","CD8_Tem","CD8_Tcm"),
#   output_dir = "/path/to/motif_deviations/",
#   use_groups = "EX",
#   bgd_groups = "SED"
# )





# ============================================
# Motif violin plot for specific cell subset
# Example: Mature_NK cells from NK project
# ============================================

library(ArchR)
library(ggplot2)
library(ggpubr)
library(gghalves)
library(tidyverse)
library(ggrepel)
library(dplyr)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

# ============================================
# Step 1: Load project
# ============================================

project_path <- "/path/to/Subset_NK"
proj <- loadArchRProject(project_path)

# ============================================
# Step 2: Extract cell subset
# ============================================

specific_cell_types <- c("Mature_NK")

idxSample <- which(proj$NK_celltype %in% specific_cell_types)
cellsSample <- proj$cellNames[idxSample]

proj_sub <- proj[cellsSample, ]

cat("Processing cell type:", paste(specific_cell_types, collapse = ", "), "\n")
cat("Cells in subset:", length(proj_sub$cellNames), "\n")

# ============================================
# Function: Match cells by sex
# ============================================

match_cells_by_sex <- function(
    archr_proj,
    group_col = "exercise_group",
    groups = c("SED", "EX"),
    seed = 123
) {
  
  cat("  Matching cells by sex...\n")
  
  cell_meta <- as.data.frame(getCellColData(archr_proj)) %>%
    mutate(
      cellName = rownames(.),
      sex = as.character(sex)
    ) %>%
    filter(!is.na(sex), !is.na(.data[[group_col]]))
  
  set.seed(seed)
  matched_cells <- cell_meta %>%
    group_by(sex, .data[[group_col]]) %>%
    mutate(n_group = n()) %>%
    group_by(sex) %>%
    mutate(n_match = min(n_group)) %>%
    ungroup() %>%
    group_by(sex, .data[[group_col]]) %>%
    slice_sample(n = first(n_match)) %>%
    ungroup() %>%
    pull(cellName)
  
  archr_proj_matched <- archr_proj[matched_cells, ]
  
  cat("    Cells before matching:", nrow(cell_meta), "\n")
  cat("    Cells after matching:", length(matched_cells), "\n")
  
  meta_check <- as.data.frame(getCellColData(archr_proj_matched))
  cat("    Sex balance after matching:\n")
  print(table(meta_check[[group_col]], meta_check$sex))
  
  return(archr_proj_matched)
}

# ============================================
# Function: Create motif violin plot
# ============================================

create_motif_violin_plot <- function(
    archr_proj,
    tf_list,
    group_col = "exercise_group",
    group_labels = c("SED", "EX"),
    colors = c("#7EBFC9", "#BCDF7A"),
    y_limit = c(-6, 7),
    match_sex = TRUE,
    return_data = FALSE
) {
  
  if (match_sex) {
    proj_plot <- match_cells_by_sex(archr_proj, group_col, group_labels)
  } else {
    proj_plot <- archr_proj
  }
  
  z_scores_se <- getMatrixFromProject(proj_plot, "MotifMatrix")
  z_matrix <- assay(z_scores_se, "z")
  
  plot_list <- lapply(tf_list, function(tf) {
    target_name <- grep(paste0("^", tf, "_"), rownames(z_matrix), value = TRUE)[1]
    if (is.na(target_name)) {
      cat("  Warning:", tf, "not found\n")
      return(NULL)
    }
    
    meta <- as.data.frame(getCellColData(proj_plot)) %>%
      rownames_to_column("cell_ID") %>%
      select(cell_ID, all_of(group_col))
    
    df <- data.frame(
      cell_ID = colnames(z_matrix),
      Values = as.numeric(z_matrix[target_name, ])
    ) %>%
      inner_join(meta, by = "cell_ID") %>%
      rename(group = !!sym(group_col)) %>%
      filter(group %in% group_labels) %>%
      mutate(Genes = tf)
    
    return(df)
  })
  
  data_long <- bind_rows(plot_list)
  data_long$group <- factor(data_long$group, levels = group_labels)
  data_long$Genes <- factor(data_long$Genes, levels = tf_list)
  
  star_height <- y_limit[2] - 1.0
  
  p <- ggplot(data_long, aes(x = Genes, y = Values, fill = group)) +
    geom_half_violin(
      data = data_long %>% filter(group == group_labels[1]),
      side = "l", alpha = 0.8, color = NA, trim = TRUE
    ) +
    geom_half_violin(
      data = data_long %>% filter(group == group_labels[2]),
      side = "r", alpha = 0.8, color = NA, trim = TRUE
    ) +
    stat_summary(
      fun.min = function(x) quantile(x, 0.25),
      fun.max = function(x) quantile(x, 0.75),
      geom = 'errorbar', color = 'black', width = 0.08, size = 0.3,
      position = position_dodge(width = 0.2)
    ) +
    stat_summary(
      fun = mean, geom = 'point', shape = 19, size = 1.0, color = "black",
      position = position_dodge(width = 0.2)
    ) +
    stat_compare_means(
      aes(group = group),
      label = "p.signif",
      method = "t.test",
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "ns")),
      label.y = star_height,
      size = 3.5
    ) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    labs(x = NULL, y = "Z-score") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.6)
    ) +
    scale_y_continuous(limits = y_limit, expand = c(0, 0))
  
  if (return_data) {
    return(list(plot = p, data = data_long))
  } else {
    return(p)
  }
}

# ============================================
# Step 3: Run the plot
# ============================================

my_tfs <- c("RUNX3", "ETS2", "ETS1", "EOMES", "STAT1", "TBX21")
my_colors <- c("#7EBFC9", "#BCDF7A")

p_final <- create_motif_violin_plot(
  archr_proj = proj_sub,
  tf_list = my_tfs,
  group_col = "exercise_group",
  group_labels = c("SED", "EX"),
  colors = my_colors,
  y_limit = c(-6, 7),
  match_sex = TRUE
)

print(p_final)

ggsave("Mature_NK_motif_violin_plot.pdf", p_final, width = 8, height = 6)

cat("Plot saved: motif_violin_plot.pdf\n")

