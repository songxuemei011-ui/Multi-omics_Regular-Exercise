#!/usr/bin/env Rscript

# ============================================
# Myeloid cell subset analysis: LSI + Harmony + Integration
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
# 1. Load project
# ============================================

proj <- loadArchRProject("./ArchRProject_Myeloid/")
seRNA <- readRDS("./data/Myeloid_drop.rds")

cat("=", rep("=", 60), "\n")
cat("Myeloid Cell Subset Analysis\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=", rep("=", 60), "\n\n")

# ============================================
# 2. LSI dimensionality reduction
# ============================================

cat("Step 1: LSI dimensionality reduction...\n")

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 3,
  clusterParams = list(
    resolution = c(0.2, 0.4, 0.5),
    sampleCells = 30000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 2:30,
  force = TRUE
)

# ============================================
# 3. Clustering
# ============================================

cat("\nStep 2: Clustering...\n")

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "MClusters",
  resolution = 0.6
)

cat("  Number of clusters:", length(unique(proj$MClusters)), "\n")
print(table(proj$MClusters))

# ============================================
# 4. Harmony batch correction
# ============================================

cat("\nStep 3: Harmony batch correction...\n")

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "HarmonyM",
  groupBy = "newSample"
)

# ============================================
# 5. UMAP visualization
# ============================================

cat("\nStep 4: UMAP visualization...\n")

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "HarmonyM",
  name = "UMAPHarmonyM",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine"
)

# ============================================
# 6. Extract RNA subsets for groupList
# ============================================

cat("\nStep 5: Extracting RNA subsets...\n")

rna_celltypes <- table(seRNA$celltype_l3)
cat("  Available RNA cell types:\n")
print(rna_celltypes)

# Myeloid subsets - modify based on your seRNA
rnaCMono <- rownames(seRNA@meta.data[seRNA$celltype_l3 == "cMono", ])
rnaNCMono <- rownames(seRNA@meta.data[seRNA$celltype_l3 == "ncMono", ])
rnaPDC <- rownames(seRNA@meta.data[seRNA$celltype_l3 == "pDC", ])

cat("  RNA subsets extracted:\n")
cat("    cMono:", length(rnaCMono), "cells\n")
cat("    ncMono:", length(rnaNCMono), "cells\n")
cat("    pDC:", length(rnaPDC), "cells\n")

# ============================================
# 7. Gene integration matrix (unsupervised)
# ============================================

cat("\nStep 6: Unsupervised gene integration...\n")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "celltype_l3",
  nameCell = "predictedCell_M",
  nameGroup = "predictedGroup_M",
  nameScore = "predictedScore_M",
  sampleCellsATAC = 30000,
  sampleCellsRNA = 30000,
  threads = 4
)

# ============================================
# 8. Group list for supervised integration
# ============================================

cat("\nStep 7: Creating group list for supervised integration...\n")

# Adjust cluster IDs based on your actual clustering results
# Run table(proj$MClusters) to see your cluster IDs
groupList <- SimpleList(
  cMono = SimpleList(
    ATAC = proj$cellNames[proj$MClusters %in% c('C1', 'C2')],
    RNA = rnaCMono
  ),
  ncMono = SimpleList(
    ATAC = proj$cellNames[proj$MClusters %in% c('C4', 'C5')],
    RNA = rnaNCMono
  ),
  pDC = SimpleList(
    ATAC = proj$cellNames[proj$MClusters %in% c('C8')],
    RNA = rnaPDC
  )

)

cat("  Group list created\n")

# ============================================
# 9. Supervised gene integration
# ============================================

cat("\nStep 8: Supervised gene integration...\n")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix_Sup",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupList = groupList,
  groupRNA = "celltype_l3",
  nameCell = "predictedCell_Sup",
  nameGroup = "predictedGroup_Sup",
  nameScore = "predictedScore_Sup",
  sampleCellsATAC = 30000,
  sampleCellsRNA = 30000,
  threads = 4
)

# ============================================
# 10. Confusion matrix
# ============================================

cat("\nStep 9: Generating confusion matrix...\n")

cM <- as.matrix(confusionMatrix(proj$MClusters, proj$predictedGroup_Sup))
preClust <- colnames(cM)[apply(cM, 1, which.max)]
confusion_result <- cbind(preClust, rownames(cM))
colnames(confusion_result) <- c("Predicted", "Cluster")
print(confusion_result)

# ============================================
# 11. Myeloid cell type annotation
# ============================================

cat("\nStep 10: Annotating Myeloid cell types...\n")

proj <- addCellColData(
  ArchRProj = proj,
  data = proj$predictedGroup_Sup,
  name = "Myeloid_celltype",
  cells = getCellNames(proj)
)

cat("  Myeloid cell type distribution:\n")
print(table(proj$Myeloid_celltype))

# ============================================
# 12. DBSCAN outlier removal
# ============================================

cat("\nStep 11: DBSCAN outlier removal...\n")

umap_coords <- getEmbedding(proj, embedding = "UMAPHarmonyM")

kNNdist <- dbscan::kNNdist(umap_coords, k = 10)
eps_auto <- quantile(kNNdist, 0.95)
cat("  Auto-selected eps:", eps_auto, "\n")

clusters <- dbscan(umap_coords, eps = eps_auto, minPts = 10)$cluster
proj$M_DBSCAN <- paste0("D", clusters)

clusters_to_keep <- unique(proj$M_DBSCAN[proj$M_DBSCAN != "D0"])
cells_to_keep <- which(proj$M_DBSCAN %in% clusters_to_keep)
proj_filtered <- proj[cells_to_keep, ]

cat("  Cells before filtering:", length(proj$cellNames), "\n")
cat("  Cells after filtering:", length(proj_filtered$cellNames), "\n")

# ============================================
# 13. Filtered UMAP visualization
# ============================================

cat("\nStep 12: Generating filtered UMAP...\n")

p_filtered <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "cellColData",
  name = "Myeloid_celltype",
  embedding = "UMAPHarmonyM",
  size = 0.01,
  alpha = 0.6,
  plotAs = "points",
  labelMeans = TRUE,
  legendPosition = "right",
  title = "Myeloid Cell Type UMAP"
)

plotPDF(p_filtered, name = "Myeloid_Celltype_UMAP.pdf",
        ArchRProj = proj_filtered, addDOC = FALSE, width = 14, height = 12)

# ============================================
# 14. Save UMAP coordinates
# ============================================

cat("\nStep 13: Saving UMAP coordinates...\n")

umap_coords_filtered <- getEmbedding(ArchRProj = proj_filtered,
                                      embedding = "UMAPHarmonyM",
                                      returnDF = TRUE)
cell_annotations <- getCellColData(ArchRProj = proj_filtered, select = "Myeloid_celltype")

umap_data <- data.frame(
  Cell_ID = rownames(umap_coords_filtered),
  UMAP_1 = umap_coords_filtered[, 1],
  UMAP_2 = umap_coords_filtered[, 2],
  CellType = cell_annotations$Myeloid_celltype
)

write.csv(umap_data, "./umap-data/Myeloid-umap.csv", row.names = FALSE)
cat("  UMAP coordinates saved: ./umap-data/Myeloid-umap.csv\n")

# ============================================
# 15. Myeloid marker genes UMAP
# ============================================

cat("\nStep 14: Plotting marker genes...\n")

markerGenes <- c(
  # Monocytes
  "CD14", "FCGR3A", "LYZ", "S100A8", "S100A9", "VCAN",
  # cDC
  "CD1C", "CLEC9A", "ITGAX", "IRF8", "CD74",
  # pDC
  "IRF7", "CLEC4C", "LILRA4", "IL3RA",
  # Common
  "PTPRC", "HLA-DRA"
)

all_genes <- getFeatures(proj_filtered, useMatrix = "GeneScoreMatrix")
existing_markerGenes <- markerGenes[markerGenes %in% all_genes]
missing_genes <- markerGenes[!markerGenes %in% all_genes]

if (length(missing_genes) > 0) {
  cat("  Missing genes:", paste(missing_genes, collapse = ", "), "\n")
}
cat("  Plotting", length(existing_markerGenes), "marker genes\n")

proj_filtered <- addImputeWeights(proj_filtered)

p <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "GeneScoreMatrix",
  name = existing_markerGenes,
  pal = c("grey", "#f32a1f"),
  embedding = "UMAPHarmonyM",
  imputeWeights = getImputeWeights(proj_filtered)
)

plotPDF(plotList = p,
  name = "Myeloid_Marker_Genes_UMAP.pdf",
  ArchRProj = proj_filtered,
  addDOC = FALSE,
  width = 5,
  height = 5)

# ============================================
# 16. Myeloid subset marker gene heatmap
# ============================================

cat("\nStep 15: Generating Myeloid subset marker heatmap...\n")

markersGS <- getMarkerFeatures(
  ArchRProj = proj_filtered,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Myeloid_celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmap_markerGenes <- c(
  # cMono
  "CD14", "S100A8", "S100A9", "VCAN", "LYZ",
  # intMono
  "CD14", "FCGR3A", "HLA-DRA",
  # ncMono
  "FCGR3A", "CX3CR1", "NCAM1",
  # cDC
  "CD1C", "ITGAX", "CLEC9A", "IRF8", "CD74",
  # pDC
  "IRF7", "CLEC4C", "LILRA4", "IL3RA"
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

  pdf("./results/Myeloid_celltype_heatmap.pdf", width = 10, height = 7)

  p_heatmap <- ComplexHeatmap::pheatmap(heatmapGS,
    show_rownames = TRUE,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    col = col,
    name = "Z-Scores",
    fontsize_row = 8,
    fontsize_col = 10,
    main = "Myeloid Cell Subset Marker Genes"
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
  cat("  Heatmap saved to: ./results/Myeloid_celltype_heatmap.pdf\n")
} else {
  cat("  No marker genes found for heatmap\n")
}

# ============================================
# 17. Save project
# ============================================

cat("\nStep 16: Saving project...\n")

saveArchRProject(
  ArchRProj = proj_filtered,
  outputDirectory = "./ArchRProject_Myeloid_annotated/",
  overwrite = TRUE
)

# ============================================
# 18. Final summary
# ============================================

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nFinal summary:\n")
cat("  - Total cells in filtered project:", nrow(getCellColData(proj_filtered)), "\n")
cat("  - Myeloid cell types:", paste(unique(proj_filtered$Myeloid_celltype), collapse=", "), "\n")
cat("\nOutput files:\n")
cat("  - ./umap-data/Myeloid-umap.csv\n")
cat("  - Myeloid_Marker_Genes_UMAP.pdf\n")
cat("  - Myeloid_Celltype_UMAP.pdf\n")
cat("  - ./results/Myeloid_celltype_heatmap.pdf\n")
cat("=", rep("=", 60), "\n")
