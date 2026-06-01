#!/usr/bin/env Rscript

# ============================================
# Browser Track Visualization for CD8_CTL
# Genes: GZMB, TBX21
# ============================================

library(ArchR)
library(grid)
library(ggplot2)

addArchRThreads(threads = 8)
addArchRGenome("hg38")

# ============================================
# Step 1: Load project
# ============================================

project_path <- "/path/to/Subset_CD8_with_peaks/"
proj <- loadArchRProject(project_path)

cat("Loaded project with", nCells(proj), "cells\n")

# ============================================
# Step 2: Add co-accessibility (if not already done)
# ============================================

if (!"CoAccessibility" %in% names(proj@reducedDims)) {
  cat("Adding co-accessibility...\n")
  proj <- addCoAccessibility(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    force = TRUE
  )
}

# ============================================
# Step 3: Extract CD8_CTL subset
# ============================================

specific_cell_types <- c("CD8_CTL")
idxSample <- which(proj$CD8_celltype %in% specific_cell_types)
cellsSample <- proj$cellNames[idxSample]

proj_sub <- proj[cellsSample, ]

cat("Extracted CD8_CTL cells:", length(proj_sub$cellNames), "\n")

# ============================================
# Step 4: Plot browser tracks for GZMB and TBX21
# ============================================

markerGenes <- c("GZMB", "TBX21")

p <- plotBrowserTrack(
  ArchRProj = proj_sub,
  groupBy = "exercise_group",
  geneSymbol = markerGenes,
  upstream = 30000,
  downstream = 30000,
  pal = c("#7EBFC9", "#BCDF7A"),
  scCellsMax = 6000
)

# ============================================
# Step 5: Display and save plots
# ============================================

# Display GZMB
grid::grid.newpage()
grid::grid.draw(p$GZMB)

# Display TBX21
grid::grid.newpage()
grid::grid.draw(p$TBX21)

# Save all tracks as PDF
plotPDF(plotList = p,
        name = "CD8_CTL_Browser_Tracks_GZMB_TBX21.pdf",
        ArchRProj = proj_sub,
        addDOC = FALSE,
        width = 8,
        height = 6)

cat("Plots saved: CD8_CTL_Browser_Tracks_GZMB_TBX21.pdf\n")

# ============================================
# Optionally save individual plots
# ============================================

pdf("GZMB_browser_track.pdf", width = 8, height = 4)
grid::grid.newpage()
grid::grid.draw(p$GZMB)
dev.off()

pdf("TBX21_browser_track.pdf", width = 8, height = 4)
grid::grid.newpage()
grid::grid.draw(p$TBX21)
dev.off()

cat("Individual plots saved: GZMB_browser_track.pdf, TBX21_browser_track.pdf\n")
