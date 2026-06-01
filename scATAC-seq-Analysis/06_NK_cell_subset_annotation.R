#!/usr/bin/env Rscript

# ============================================
# NK cell subset analysis
# LSI + Harmony + Integration + Annotation
# ============================================

library(ArchR)
library(dbscan)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ComplexHeatmap)
library(circlize)

addArchRThreads(threads = 16)
addArchRGenome("hg38")

# Create output directories
dir.create("./results", showWarnings = FALSE)
dir.create("./umap-data", showWarnings = FALSE)

# ============================================
# 1. Load project and RNA reference
# ============================================

proj <- loadArchRProject("/scATAC/Subset_NK/")
seRNA <- readRDS("./scRNA/NK_drop.rds")

cat("=", rep("=", 60), "\n")
cat("NK Cell Subset Analysis\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=", rep("=", 60), "\n\n")

# ============================================
# 2. LSI dimensionality reduction
# ============================================

cat("Step 1: LSI dimensionality reduction...\n")

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSINK",
  iterations = 3,
  clusterParams = list(
    resolution = c(0.2, 0.4, 0.5),
    sampleCells = 50000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 2:30,
  force = TRUE
)

# ============================================
# 3. Harmony batch correction
# ============================================

cat("\nStep 2: Harmony batch correction...\n")

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSINK",
  name = "HarmonyNK",
  groupBy = "newSample",
  force = TRUE
)

# ============================================
# 4. Clustering
# ============================================

cat("\nStep 3: Clustering...\n")

proj <- addClusters(
  input = proj,
  reducedDims = "HarmonyNK",
  method = "Seurat",
  name = "ClustersNK",
  resolution = 0.5,
  force = TRUE,
  maxClusters = 12
)

cat("  Number of clusters:", length(unique(proj$ClustersNK)), "\n")
print(table(proj$ClustersNK))

# ============================================
# 5. UMAP visualization
# ============================================

cat("\nStep 4: UMAP visualization...\n")

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "HarmonyNK",
  name = "UMAPHarmony_NK",
  nNeighbors = 50,
  minDist = 0.4,
  metric = "cosine"
)

# ============================================
# 6. Unsupervised gene integration (no groupList)
# ============================================

cat("\nStep 5: Unsupervised gene integration...\n")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "HarmonyNK",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "celltype_l3",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  sampleCellsATAC = 30000,
  sampleCellsRNA = 30000,
  threads = 8,
  force = TRUE
)

# ============================================
# 7. NK cell type annotation
# ============================================

cat("\nStep 6: NK cell type annotation...\n")

original_labels <- getCellColData(proj, "predictedGroup")[, 1]

# Merge rules
new_labels <- ifelse(
  original_labels %in% c("Cycling_NK", "ILC2"),
  ifelse(
    original_labels == "Cycling_NK",
    "CD56_bright_NK",
    "Transitional_NK"
  ),
  as.character(original_labels)
)

proj$NK_celltype <- new_labels

cat("  NK cell type distribution:\n")
print(table(proj$NK_celltype))

# ============================================
# 8. DBSCAN outlier removal
# ============================================

cat("\nStep 7: DBSCAN outlier removal...\n")

umap_coords <- getEmbedding(proj, embedding = "UMAPHarmony_NK")

kNNdist <- dbscan::kNNdist(umap_coords, k = 10)
eps_auto <- quantile(kNNdist, 0.95)
cat("  Auto-selected eps:", eps_auto, "\n")

clusters <- dbscan(umap_coords, eps = eps_auto, minPts = 10)$cluster
proj$NK_DBSCAN <- paste0("D", clusters)

clusters_to_keep <- unique(proj$NK_DBSCAN[proj$NK_DBSCAN != "D0"])
cells_to_keep <- which(proj$NK_DBSCAN %in% clusters_to_keep)
proj_filtered <- proj[cells_to_keep, ]

cat("  Cells before filtering:", length(proj$cellNames), "\n")
cat("  Cells after filtering:", length(proj_filtered$cellNames), "\n")

# ============================================
# 9. Filtered UMAP visualization
# ============================================

cat("\nStep 8: Generating filtered UMAP...\n")

p_anno <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "cellColData",
  name = "NK_celltype",
  embedding = "UMAPHarmony_NK",
  size = 0.01,
  alpha = 0.6,
  plotAs = "points",
  labelMeans = TRUE,
  legendPosition = "right",
  title = "NK Cell Type UMAP"
)

plotPDF(p_anno, name = "NK_Celltype_UMAP.pdf",
        ArchRProj = proj_filtered, addDOC = FALSE, width = 14, height = 12)

# ============================================
# 10. Save UMAP coordinates
# ============================================

cat("\nStep 9: Saving UMAP coordinates...\n")

umap_coords_filtered <- getEmbedding(ArchRProj = proj_filtered,
                                      embedding = "UMAPHarmony_NK",
                                      returnDF = TRUE)
cell_annotations <- getCellColData(ArchRProj = proj_filtered, select = "NK_celltype")

umap_data <- data.frame(
  Cell_ID = rownames(umap_coords_filtered),
  UMAP_1 = umap_coords_filtered[, 1],
  UMAP_2 = umap_coords_filtered[, 2],
  CellType = cell_annotations$NK_celltype
)

write.csv(umap_data, "./umap-data/NK_umap.csv", row.names = FALSE)
cat("  UMAP coordinates saved: ./umap-data/NK_umap.csv\n")

# ============================================
# 11. Add imputation weights
# ============================================

proj_filtered <- addImputeWeights(proj_filtered)

# ============================================
# 12. NK marker genes UMAP
# ============================================

cat("\nStep 10: Plotting marker genes...\n")

markerGenes <- c(
  "NCAM1", "SELL", "IL2RB", "CD44", "KLRB1", "XCL1", "IFNG", "CCR7",
  "FCGR3A", "PRF1", "GZMB", "KLRD1", "NCR1", "KLRK1", "CD7", "TYROBP",
  "KLRG1", "B3GAT1", "KLRC1", "CD160", "LAG3", "HAVCR2", "CX3CR1",
  "CD27", "ITGAM", "TBX21", "ZEB2", "IKZF2", "CD244", "KLRF1"
)

all_genes <- getFeatures(proj_filtered, useMatrix = "GeneScoreMatrix")
existing_markerGenes <- markerGenes[markerGenes %in% all_genes]
missing_genes <- markerGenes[!markerGenes %in% all_genes]

if (length(missing_genes) > 0) {
  cat("  Missing genes:", paste(missing_genes, collapse = ", "), "\n")
}
cat("  Plotting", length(existing_markerGenes), "marker genes\n")

p <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "GeneScoreMatrix",
  name = existing_markerGenes,
  pal = c("grey", "#f32a1f"),
  embedding = "UMAPHarmony_NK",
  imputeWeights = getImputeWeights(proj_filtered)
)

plotPDF(plotList = p,
  name = "NK_Marker_Genes_UMAP.pdf",
  ArchRProj = proj_filtered,
  addDOC = FALSE,
  width = 5,
  height = 5)

# ============================================
# 13. NK subset marker gene heatmap
# ============================================

cat("\nStep 11: Generating NK subset marker heatmap...\n")

markersGS <- getMarkerFeatures(
  ArchRProj = proj_filtered,
  useMatrix = "GeneScoreMatrix",
  groupBy = "NK_celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmap_markerGenes <- c(
  "NCAM1", "SELL", "IL2RB", "KLRB1", "CCR7",  # CD56_bright_NK
  "FCGR3A", "PRF1", "GZMB", "KLRD1", "NCR1",  # Mature_NK
  "KLRG1", "B3GAT1", "CD160", "CX3CR1",       # Terminal_NK
  "CD27", "ITGAM", "ZEB2", "CD244"             # Transitional_NK
)

all_genes <- getFeatures(proj_filtered, useMatrix = "GeneScoreMatrix")
existing_heatmap_genes <- heatmap_markerGenes[heatmap_markerGenes %in% all_genes]
missing_genes <- heatmap_markerGenes[!heatmap_markerGenes %in% all_genes]

if (length(missing_genes) > 0) {
  cat("  Missing genes:", paste(missing_genes, collapse = ", "), "\n")
} else {
  cat("  All marker genes found!\n")
}

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")

cell_types <- names(markerList)
all_markers <- data.frame()

for (ct in cell_types) {
  markers <- markerList[[ct]]
  if (!is.null(markers) && nrow(markers) > 0) {
    markers_df <- data.frame(
      name = markers$name,
      Log2FC = markers$Log2FC,
      cell_type = ct
    )
    all_markers <- rbind(all_markers, markers_df)
  }
}

if (nrow(all_markers) > 0) {
  all_markers <- all_markers[order(all_markers$cell_type, -all_markers$Log2FC),]
  all_markers$rank <- 1:nrow(all_markers)

  gene_freq <- as.data.frame(table(all_markers$name))
  dup_genes <- as.character(gene_freq[gene_freq$Freq > 1, "Var1"])

  markerList_unique <- all_markers[!(all_markers$name %in% dup_genes),]
  markerList_dup <- all_markers[all_markers$name %in% dup_genes,]
  if (nrow(markerList_dup) > 0) {
    markerList_dup <- markerList_dup %>% group_by(name) %>% filter(Log2FC == max(Log2FC))
  }
  all_markers <- rbind(markerList_unique, markerList_dup)
  all_markers <- all_markers[order(all_markers$rank),]
}

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  labelMarkers = existing_heatmap_genes,
  binaryClusterRows = TRUE,
  limits = c(-2, 2),
  transpose = TRUE,
  returnMatrix = TRUE
)

if (!is.null(heatmapGS) && nrow(all_markers) > 0) {
  heatmapGS <- heatmapGS[, colnames(heatmapGS) %in% all_markers$name, drop = FALSE]
  heatmapGS <- t(heatmapGS)

  col <- paletteContinuous(set = "blueYellow")

  pdf("./results/NK_celltype_heatmap.pdf", width = 10, height = 6)

  p_heatmap <- ComplexHeatmap::pheatmap(heatmapGS,
    show_rownames = TRUE,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    col = col,
    name = "Z-Scores",
    fontsize_row = 8,
    fontsize_col = 10,
    main = "NK Cell Subset Marker Genes"
  )

  idx <- which(rownames(heatmapGS) %in% existing_heatmap_genes)
  if (length(idx) > 0) {
    print(p_heatmap + rowAnnotation(link = anno_mark(
      at = idx,
      labels = rownames(heatmapGS)[idx],
      labels_gp = gpar(fontsize = 8)
    )))
  } else {
    print(p_heatmap)
  }

  dev.off()
  cat("  Heatmap saved to: ./results/NK_celltype_heatmap.pdf\n")
} else {
  cat("  No marker genes found for heatmap\n")
}

# ============================================
# 14. Save project
# ============================================

cat("\nStep 12: Saving project...\n")

saveArchRProject(
  ArchRProj = proj_filtered,
  outputDirectory = "./ArchRProject_NK_annotated/",
  overwrite = TRUE
)

# ============================================
# 15. Final summary
# ============================================

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nFinal summary:\n")
cat("  - Total cells in filtered project:", nrow(getCellColData(proj_filtered)), "\n")
cat("  - NK cell types:", paste(unique(proj_filtered$NK_celltype), collapse=", "), "\n")
cat("\nOutput files:\n")
cat("  - ./umap-data/NK_umap.csv\n")
cat("  - NK_Marker_Genes_UMAP.pdf\n")
cat("  - NK_Celltype_UMAP.pdf\n")
cat("  - ./results/NK_celltype_heatmap.pdf\n")
cat("=", rep("=", 60), "\n")
