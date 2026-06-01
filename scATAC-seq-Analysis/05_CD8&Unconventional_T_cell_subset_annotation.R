#!/usr/bin/env Rscript

# ============================================
# CD8 & Unconventional T cell subset analysis
# LSI + Harmony + Integration (unsupervised + supervised)
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

proj <- loadArchRProject("/scATAC/Subset_CD8/")
seRNA <- readRDS("./data/CD8_drop.rds")

cat("=", rep("=", 60), "\n")
cat("CD8 & Unconventional T Cell Subset Analysis\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=", rep("=", 60), "\n\n")

# ============================================
# 2. LSI dimensionality reduction
# ============================================

cat("Step 1: LSI dimensionality reduction...\n")

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSICD8",
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
  reducedDims = "IterativeLSICD8",
  name = "HarmonyCD8",
  groupBy = "newSample",
  force = TRUE
)

# ============================================
# 4. Clustering
# ============================================

cat("\nStep 3: Clustering...\n")

proj <- addClusters(
  input = proj,
  reducedDims = "HarmonyCD8",
  method = "Seurat",
  name = "Clusters0.5",
  resolution = 0.5,
  force = TRUE,
  maxClusters = 12
)

cat("  Number of clusters:", length(unique(proj$Clusters0.5)), "\n")
print(table(proj$Clusters0.5))

# ============================================
# 5. UMAP visualization
# ============================================

cat("\nStep 4: UMAP visualization...\n")

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "HarmonyCD8",
  name = "UMAPHarmony",
  nNeighbors = 20,
  minDist = 0.5,
  metric = "cosine"
)

# ============================================
# 6. Unsupervised gene integration
# ============================================

cat("\nStep 5: Unsupervised gene integration...\n")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_Unsup",
  reducedDims = "IterativeLSICD8",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "celltype_l3",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  sampleCellsATAC = 50000,
  sampleCellsRNA = 50000,
  threads = 8,
  force = TRUE
)

# ============================================
# 7. Extract RNA subsets for supervised integration
# ============================================

cat("\nStep 6: Extracting RNA subsets...\n")

# Confident RNA subsets
CD8Naive <- colnames(seRNA)[grep('CD8_naive', seRNA$celltype_l3)]
MAIT <- colnames(seRNA)[grep('MAIT', seRNA$celltype_l3)]

# Patterns for other T cells
patterns <- c('CD8_CTL', 'CD8_Tem', 'NKT', 'gdT2', 'CD8_Tcm', 'DNT', 'gdT1', 'Cycling_T')

groupList <- SimpleList(
  CD8Naive = SimpleList(
    ATAC = proj$cellNames[proj$predictedGroup_Un %in% c('CD8_naive')],
    RNA = CD8Naive
  ),
  MAIT = SimpleList(
    ATAC = proj$cellNames[proj$predictedGroup_Un %in% c('MAIT')],
    RNA = MAIT
  ),
  Other = SimpleList(
    ATAC = proj$cellNames[proj$predictedGroup_Un %in% patterns],
    RNA = colnames(seRNA)[grep(paste(patterns, collapse = "|"), seRNA$celltype_l3)]
  )
)

cat("  Group list cell counts:\n")
print(lapply(groupList, function(x) {
  list(ATAC = length(x$ATAC), RNA = length(x$RNA))
}))

# ============================================
# 8. Supervised gene integration (constrained)
# ============================================

cat("\nStep 7: Supervised gene integration...\n")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_Sup",
  reducedDims = "IterativeLSICD8",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupList = groupList,
  groupRNA = "celltype_l3",
  nameCell = "predictedCell_Sup",
  nameGroup = "predictedGroup_Sup",
  nameScore = "predictedScore_Sup",
  sampleCellsATAC = 50000,
  sampleCellsRNA = 50000,
  threads = 8,
  force = TRUE
)

# ============================================
# 9. Cell type annotation with merge rules
# ============================================

cat("\nStep 8: Cell type annotation...\n")

original_labels <- getCellColData(proj, "predictedGroup_Un")[, 1]

# Merge rules
new_labels <- ifelse(
  original_labels %in% c("gdT1", "gdT2"),
  "gdT",
  ifelse(
    original_labels == "DNT",
    "CD8_Tcm",
    as.character(original_labels)
  )
)

proj$CD8_celltype <- new_labels

cat("  Cell type distribution:\n")
print(table(proj$CD8_celltype))

# ============================================
# 10. DBSCAN outlier removal
# ============================================

cat("\nStep 9: DBSCAN outlier removal...\n")

umap_coords <- getEmbedding(proj, embedding = "UMAPHarmony")

kNNdist <- dbscan::kNNdist(umap_coords, k = 10)
eps_auto <- quantile(kNNdist, 0.95)
cat("  Auto-selected eps:", eps_auto, "\n")

clusters <- dbscan(umap_coords, eps = eps_auto, minPts = 10)$cluster
proj$CD8_DBSCAN <- paste0("D", clusters)

clusters_to_keep <- unique(proj$CD8_DBSCAN[proj$CD8_DBSCAN != "D0"])
cells_to_keep <- which(proj$CD8_DBSCAN %in% clusters_to_keep)
proj_filtered <- proj[cells_to_keep, ]

cat("  Cells before filtering:", length(proj$cellNames), "\n")
cat("  Cells after filtering:", length(proj_filtered$cellNames), "\n")

# ============================================
# 11. Filtered UMAP visualization
# ============================================

cat("\nStep 10: Generating filtered UMAP...\n")

p_filtered <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "cellColData",
  name = "CD8_celltype",
  embedding = "UMAPHarmony",
  size = 0.01,
  alpha = 0.6,
  plotAs = "points",
  labelMeans = TRUE,
  legendPosition = "right",
  title = "CD8 & Unconventional T Cell UMAP"
)

plotPDF(p_filtered, name = "CD8_Celltype_UMAP.pdf",
        ArchRProj = proj_filtered, addDOC = FALSE, width = 14, height = 12)

# ============================================
# 12. Save UMAP coordinates
# ============================================

cat("\nStep 11: Saving UMAP coordinates...\n")

umap_coords_filtered <- getEmbedding(ArchRProj = proj_filtered,
                                      embedding = "UMAPHarmony",
                                      returnDF = TRUE)
cell_annotations <- getCellColData(ArchRProj = proj_filtered, select = "CD8_celltype")

umap_data <- data.frame(
  Cell_ID = rownames(umap_coords_filtered),
  UMAP_1 = umap_coords_filtered[, 1],
  UMAP_2 = umap_coords_filtered[, 2],
  CellType = cell_annotations$CD8_celltype
)

write.csv(umap_data, "./umap-data/CD8_umap.csv", row.names = FALSE)
cat("  UMAP coordinates saved: ./umap-data/CD8_umap.csv\n")

# ============================================
# 13. Add imputation weights
# ============================================

proj_filtered <- addImputeWeights(proj_filtered)

# ============================================
# 14. Marker genes UMAP
# ============================================

cat("\nStep 12: Plotting marker genes...\n")

markerGenes <- c(
  "CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "CX3CR1", "TBX21",
  "CCR7", "LEF1", "TCF7", "SELL", "CD27", "CD28",
  "IL2RB", "GZMK", "CD44",
  "CXCR3", "CD69", "RUNX3",
  "MKI67", "TOP2A", "PCNA", "STMN1",
  "TRDC", "TRGC1", "TRDV1", "KLRB1", "RORC",
  "SLC4A10", "RORA", "TRAV1-2",
  "CD161", "FCGR3A", "PLZF"
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
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(proj_filtered)
)

plotPDF(plotList = p,
  name = "CD8_Marker_Genes_UMAP.pdf",
  ArchRProj = proj_filtered,
  addDOC = FALSE,
  width = 5,
  height = 5)

# ============================================
# 15. CD8 subset marker gene heatmap
# ============================================

cat("\nStep 13: Generating CD8 subset marker heatmap...\n")

markersGS <- getMarkerFeatures(
  ArchRProj = proj_filtered,
  useMatrix = "GeneScoreMatrix",
  groupBy = "CD8_celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmap_markerGenes <- c(
  "CD8A", "CD8B", "GZMB", "PRF1", "NKG7",
  "CCR7", "LEF1", "TCF7", "SELL", "CD27",
  "IL2RB", "GZMK",
  "CXCR3", "CD69", "RUNX3",
  "MKI67", "TOP2A", "STMN1",
  "TRDC", "TRGC1", "KLRB1", "RORC",
  "SLC4A10", "RORA", "TRAV1-2",
  "CD161", "FCGR3A", "PLZF"
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

  pdf("./results/CD8_celltype_heatmap.pdf", width = 12, height = 8)

  p_heatmap <- ComplexHeatmap::pheatmap(heatmapGS,
    show_rownames = TRUE,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    col = col,
    name = "Z-Scores",
    fontsize_row = 8,
    fontsize_col = 10,
    main = "CD8 & Unconventional T Cell Marker Genes"
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
  cat("  Heatmap saved to: ./results/CD8_celltype_heatmap.pdf\n")
} else {
  cat("  No marker genes found for heatmap\n")
}

# ============================================
# 16. Save project
# ============================================

cat("\nStep 14: Saving project...\n")

saveArchRProject(
  ArchRProj = proj_filtered,
  outputDirectory = "./ArchRProject_CD8_annotated/",
  overwrite = TRUE
)

# ============================================
# 17. Final summary
# ============================================

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nFinal summary:\n")
cat("  - Total cells in filtered project:", nrow(getCellColData(proj_filtered)), "\n")
cat("  - Cell types:", paste(unique(proj_filtered$CD8_celltype), collapse=", "), "\n")
cat("\nOutput files:\n")
cat("  - ./umap-data/CD8_umap.csv\n")
cat("  - CD8_Marker_Genes_UMAP.pdf\n")
cat("  - CD8_Celltype_UMAP.pdf\n")
cat("  - ./results/CD8_celltype_heatmap.pdf\n")
cat("=", rep("=", 60), "\n")
