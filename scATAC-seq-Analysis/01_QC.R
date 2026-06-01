# ============================================
# ArchR Analysis Pipeline for Single-cell ATAC-seq
# ============================================
# # Data source:
#   - Fragment files downloaded from CIMA database (https://db.cngb.org/trueblood/cima/resource)
#   - Sample list and metadata used in this study are provided in Supplementary Table 1
# This script performs:
#   1. Create Arrow files from fragment files
#   2. Doublet detection and filtering
#   3. Quality control
#   4. Create ArchR project with metadata
#
# Usage:
#   Modify USER SETTINGS below, then run the script
#
# ============================================

library(ArchR)
library(parallel)

# ============================================
# USER SETTINGS - MODIFY THESE
# ============================================

# Genome version
GENOME <- "hg38"

# Paths
FRAGMENT_DIR <- "data/fragments/"        # Directory containing .tsv.gz fragment files
OUTPUT_DIR <- "results/ArchR/"           # Output directory for Arrow files and project
METADATA_FILE <- "data/sample_metadata.csv"  # Sample metadata CSV

# QC thresholds
MIN_TSS <- 10
MIN_FRAGS <- 1000
MAX_FRAGS <- 100000
MAX_NUCLEOSOME_SIGNAL <- 4
MAX_BLACKLIST_RATIO <- 0.05

# Doublet detection parameters
DOUBLET_K <- 10
DOUBLET_KNN_METHOD <- "UMAP"

# Threads (adjust to your system)
THREADS <- 16

# Random seed
set.seed(1)

# ============================================
# DO NOT MODIFY BELOW THIS LINE
# ============================================

# Initialize ArchR
addArchRGenome(GENOME)
addArchRThreads(threads = THREADS)

# ============================================
# 1. Get fragment files
# ============================================

cat("=", rep("=", 60), "\n")
cat("ArchR Analysis Pipeline\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=", rep("=", 60), "\n\n")

cat("Step 1: Locating fragment files...\n")

files <- list.files(
    FRAGMENT_DIR, 
    full.names = TRUE, 
    recursive = TRUE, 
    pattern = "\\.tsv\\.gz$"
)

# Extract sample names from directory structure
sample_names <- basename(dirname(files))

sample_info <- data.frame(
    Sample = sample_names,
    FragmentFiles = files,
    stringsAsFactors = FALSE
)

cat("  Found", nrow(sample_info), "samples\n")
print(sample_info$Sample)

# ============================================
# 2. Create Arrow files
# ============================================

cat("\nStep 2: Creating Arrow files...\n")

arrow_dir <- file.path(OUTPUT_DIR, "ArrowFiles")
dir.create(arrow_dir, showWarnings = FALSE, recursive = TRUE)

ArrowFiles <- createArrowFiles(
    inputFiles = sample_info$FragmentFiles,
    sampleNames = sample_info$Sample,
    outputDir = arrow_dir,
    minTSS = MIN_TSS,
    minFrags = MIN_FRAGS,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    force = TRUE
)

# ============================================
# 3. Doublet detection
# ============================================

cat("\nStep 3: Doublet detection...\n")

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = DOUBLET_K,
    knnMethod = DOUBLET_KNN_METHOD,
    LSIMethod = 1
)

# ============================================
# 4. Create ArchR project
# ============================================

cat("\nStep 4: Creating ArchR project...\n")

proj_dir <- file.path(OUTPUT_DIR, "ArchRProject")
dir.create(proj_dir, showWarnings = FALSE, recursive = TRUE)

proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = proj_dir,
    copyArrows = TRUE
)

# Extract base sample names (remove technical suffixes)
proj$newSample <- sub("([A-Za-z0-9_]+)-[A-Za-z]-[0-9]+", "\\1", proj$Sample)
cat("  Samples:", paste(unique(proj$newSample), collapse=", "), "\n")

# ============================================
# 5. Add sample metadata
# ============================================

cat("\nStep 5: Adding sample metadata...\n")

if(file.exists(METADATA_FILE)) {
    sample_meta <- read.csv(METADATA_FILE)
    
    # Function to add metadata columns
    add_meta <- function(proj, meta_df, sample_col, value_col, new_name) {
        mapping <- setNames(meta_df[[value_col]], meta_df[[sample_col]])
        cell_values <- mapping[proj$newSample]
        proj <- addCellColData(
            ArchRProj = proj,
            data = cell_values,
            name = new_name,
            cells = proj$cellNames,
            force = TRUE
        )
        return(proj)
    }
    
    # Add columns (modify according to your CSV)
    if("exercise_group" %in% colnames(sample_meta)) {
        proj <- add_meta(proj, sample_meta, "Sample_name", "exercise_group", "exercise_group")
    }
    if("Age" %in% colnames(sample_meta)) {
        proj <- add_meta(proj, sample_meta, "Sample_name", "Age", "Age")
    }
    if("Sex" %in% colnames(sample_meta)) {
        proj <- add_meta(proj, sample_meta, "Sample_name", "Sex", "Sex")
    }
    if("BMI" %in% colnames(sample_meta)) {
        proj <- add_meta(proj, sample_meta, "Sample_name", "BMI", "BMI")
    }
    
    cat("  Metadata added successfully\n")
} else {
    cat("  Warning: Metadata file not found:", METADATA_FILE, "\n")
    cat("  Skipping metadata addition\n")
}

# ============================================
# 6. Quality control filtering
# ============================================

cat("\nStep 6: Quality control filtering...\n")

# Plot QC before filtering
p_tss <- plotGroups(proj, groupBy = "newSample", colorBy = "cellColData",
                    name = "TSSEnrichment", plotAs = "ridges")
p_frag <- plotGroups(proj, groupBy = "newSample", colorBy = "cellColData",
                     name = "log10(nFrags)", plotAs = "ridges")

plotPDF(p_tss, p_frag,
        name = "QC_before_filter",
        ArchRProj = proj,
        addDOC = FALSE, width = 8, height = 6)

# Filter doublets
proj <- filterDoublets(proj)
cat("  After doublet removal:", nrow(getCellColData(proj)), "cells\n")

# Calculate nucleosome signal
proj$NucleosomeSignal <- proj$nMonoFrags / proj$nDiFrags

# Apply QC filters
idxPass <- which(
    proj$TSSEnrichment >= MIN_TSS &
    proj$nFrags >= MIN_FRAGS &
    proj$nFrags <= MAX_FRAGS &
    proj$NucleosomeSignal < MAX_NUCLEOSOME_SIGNAL &
    proj$BlacklistRatio < MAX_BLACKLIST_RATIO
)

cellsPass <- proj$cellNames[idxPass]
cat("  After QC filtering:", length(cellsPass), "cells\n")

proj <- proj[cellsPass, ]

# Plot QC after filtering
p_tss2 <- plotGroups(proj, groupBy = "newSample", colorBy = "cellColData",
                     name = "TSSEnrichment", plotAs = "ridges")
p_frag2 <- plotGroups(proj, groupBy = "newSample", colorBy = "cellColData",
                      name = "log10(nFrags)", plotAs = "ridges")

plotPDF(p_tss2, p_frag2,
        name = "QC_after_filter",
        ArchRProj = proj,
        addDOC = FALSE, width = 8, height = 6)

# ============================================
# 7. Save project
# ============================================

cat("\nStep 7: Saving project...\n")

saveArchRProject(
    ArchRProj = proj,
    outputDirectory = proj_dir,
    overwrite = TRUE
)

# ============================================
# 8. Summary
# ============================================

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nFinal summary:\n")
cat("  - Total cells:", nrow(getCellColData(proj)), "\n")
cat("  - Samples:", length(unique(proj$newSample)), "\n")
cat("  - Output directory:", proj_dir, "\n")
cat("=", rep("=", 60), "\n")

# ============================================
# QC Visualization with ArchR
# ============================================

# Calculate nucleosome signal
proj$NucleosomeSignal <- proj$nMonoFrags / proj$nDiFrags

# Generate QC plots
p_tss_ridge <- plotGroups(proj, groupBy = "exercise_group", colorBy = "cellColData", 
                          name = "TSSEnrichment", plotAs = "ridges")

p_tss_violin <- plotGroups(proj, groupBy = "exercise_group", colorBy = "cellColData", 
                           name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

p_frags <- plotGroups(proj, groupBy = "exercise_group", colorBy = "cellColData", 
                      name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

p_nucleo <- plotGroups(proj, groupBy = "exercise_group", colorBy = "cellColData", 
                       name = "NucleosomeSignal", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

p_blacklist <- plotGroups(proj, groupBy = "exercise_group", colorBy = "cellColData", 
                          name = "BlacklistRatio", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)

# Save plots
plotPDF(p_tss_ridge, p_tss_violin, p_frags, p_nucleo, p_blacklist,
        name = "QC_Statistics.pdf",
        ArchRProj = proj,
        addDOC = FALSE,
        width = 10,
        height = 8)

cat("QC plots saved: QC_Statistics.pdf\n")







