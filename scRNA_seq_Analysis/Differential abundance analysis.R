# ============================================
# Milo Differential Abundance Analysis
# Compare EX vs SED, adjusted for sex and BMI
# ============================================

library(miloR)
library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scales)

# ============================================
# USER SETTINGS - MODIFY THESE
# ============================================

INPUT_FILE <- "/data/work/Myeloid.rds"  # Change to your file
OUTPUT_DIR <- "./milo_results/"
SAMPLE_COL <- "sample"
CELL_TYPE_COL <- "cell_type_l3"

# ============================================
# DO NOT MODIFY BELOW
# ============================================

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=", rep("=", 60), "\n")
cat("Milo Differential Abundance Analysis\n")
cat("EX vs SED (adjusted for sex + BMI)\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=", rep("=", 60), "\n\n")

# 1. Load data
cat("Step 1: Loading data...\n")
sc <- readRDS(INPUT_FILE)
cat("  Cells:", ncol(sc), "\n")

# 2. Subset to SED and EX
cat("\nStep 2: Subsetting SED vs EX...\n")
Idents(sc) <- "exercise_group"
sc <- subset(sc, ident = c("SED", "EX"))
cat("  Cells after subset:", ncol(sc), "\n")
print(table(sc$exercise_group))

# 3. Convert to SingleCellExperiment
cat("\nStep 3: Converting to SingleCellExperiment...\n")
sce <- as.SingleCellExperiment(sc)

# 4. Create Milo object
cat("\nStep 4: Creating Milo object...\n")
milo <- miloR::Milo(sce)

# 5. Build graph
cat("\nStep 5: Building graph (k=30, d=50)...\n")
milo <- miloR::buildGraph(milo, k = 30, d = 50)

# 6. Make neighborhoods
cat("\nStep 6: Defining neighborhoods...\n")
milo <- makeNhoods(milo, prop = 0.2, k = 30, d = 50, refined = TRUE)
cat("  Number of neighborhoods:", nrow(nhoods(milo)), "\n")

# 7. Count cells per neighborhood
cat("\nStep 7: Counting cells...\n")
milo <- countCells(milo, 
                   meta.data = data.frame(colData(milo)),
                   samples = SAMPLE_COL)

# 8. Prepare design matrix
cat("\nStep 8: Preparing design matrix...\n")
traj_design <- data.frame(colData(milo))[, c("sample", "exercise_group", "sex", "BMI")]
traj_design$sample <- as.factor(traj_design$sample)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$sample
print(head(traj_design))

# 9. Calculate neighborhood distances
cat("\nStep 9: Calculating distances...\n")
milo <- calcNhoodDistance(milo, d = 50)

# 10. Run differential abundance test (adjusted for sex + BMI)
cat("\nStep 10: Running DA test (EX vs SED, adjusted for sex + BMI)...\n")
da_results <- testNhoods(milo, 
                         design = ~ sex + BMI + exercise_group,
                         design.df = traj_design)
da_results$SpatialFDR <- p.adjust(da_results$PValue, method = "BH")
cat("  Significant neighborhoods (FDR < 0.05):", sum(da_results$SpatialFDR < 0.05, na.rm = TRUE), "\n")

# 11. Annotate with cell types
cat("\nStep 11: Annotating neighborhoods...\n")
if(CELL_TYPE_COL %in% colnames(colData(milo))) {
  da_results <- annotateNhoods(milo, da_results, coldata_col = CELL_TYPE_COL)
} else {
  cat("  Warning:", CELL_TYPE_COL, "not found, skipping annotation\n")
}

# 12. Save results
cat("\nStep 12: Saving results...\n")
write.csv(da_results, file.path(OUTPUT_DIR, "DA_results.csv"), row.names = FALSE)
cat("  Saved:", file.path(OUTPUT_DIR, "DA_results.csv"), "\n")

# 13. Generate plot (using your original plot code)
cat("\nStep 13: Generating plot...\n")

# Convert to data.frame for plotting
da_results_df <- as.data.frame(da_results)

p <- plotDAbeeswarm(da_results_df, group.by = CELL_TYPE_COL, alpha = 1) +
  geom_boxplot(outlier.shape = NA) +
  scale_color_gradient2(midpoint = 0, low = "#1A4D8C", mid = "white", high = "#3D8E3D", space = "Lab") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "", y = "Log2 Fold Change (EX vs SED)") +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = 'black'))

print(p)
ggsave(file.path(OUTPUT_DIR, "DA_beeswarm.pdf"), p, width = 10, height = 6)
cat("  Saved:", file.path(OUTPUT_DIR, "DA_beeswarm.pdf"), "\n")

# 14. Summary
cat("\n", rep("=", 60), "\n")
cat("SUMMARY\n")
cat("=", rep("=", 60), "\n")
cat("Total neighborhoods:", nrow(da_results), "\n")
cat("Significant (FDR < 0.05):", sum(da_results$SpatialFDR < 0.05, na.rm = TRUE), "\n")

# Cell type breakdown
if(CELL_TYPE_COL %in% colnames(da_results)) {
  cat("\nSignificant by cell type:\n")
  sig_table <- da_results[da_results$SpatialFDR < 0.05, ] %>%
    group_by(!!sym(CELL_TYPE_COL)) %>%
    summarise(Count = n())
  print(sig_table)
}

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output saved to:", OUTPUT_DIR, "\n")
cat("=", rep("=", 60), "\n")
