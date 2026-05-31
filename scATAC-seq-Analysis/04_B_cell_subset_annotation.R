#!/usr/bin/env Rscript

# ============================================
# B cell subset analysis: LSI + Harmony + Integration
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

proj <- loadArchRProject("/scATAC/Subset_B/")
seRNA <- readRDS("/scRNA/B-drop.rds")

cat("=", rep("=", 60), "\n")
cat("B Cell Subset Analysis\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("=", rep("=", 60), "\n\n")

# ============================================
# 2. LSI dimensionality reduction
# ============================================

cat("Step 1: LSI dimensionality reduction...\n")

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSIB",
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
# 3. Clustering
# ============================================

cat("\nStep 2: Clustering...\n")

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSIB",
  method = "Seurat",
  name = "BClusters0.6",
  resolution = 0.6
)

cat("  Number of clusters:", length(unique(proj$BClusters0.6)), "\n")

# ============================================
# 4. Harmony batch correction
# ============================================

cat("\nStep 3: Harmony batch correction...\n")

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "HarmonyB",
  groupBy = "newSample"
)

# ============================================
# 5. UMAP visualization
# ============================================

cat("\nStep 4: UMAP visualization...\n")

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "HarmonyB",
  name = "UMAPHarmonyB_50_0.4",
  nNeighbors = 50,
  minDist = 0.4,
  metric = "cosine"
)

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "IterativeLSIB",
  name = "UMAP40_0.35",
  nNeighbors = 40,
  minDist = 0.35,
  metric = "cosine"
)

# ============================================
# 6. Extract RNA subsets for groupList
# ============================================

cat("\nStep 5: Extracting RNA subsets...\n")

rna_celltypes <- table(seRNA$celltype_l3)
cat("  Available RNA cell types:\n")
print(rna_celltypes)

rnaNaiveB <- rownames(seRNA@meta.data[seRNA$celltype_l3 == "Naive_B", ])
rnaMemoryB <- rownames(seRNA@meta.data[seRNA$celltype_l3 == "Memory_B", ])
rnaSwitchedB <- rownames(seRNA@meta.data[seRNA$celltype_l3 == "Switched_memory_B", ])
rnaPlasma <- rownames(seRNA@meta.data[seRNA$celltype_l3 == "Plasma", ])

cat("  RNA subsets extracted:\n")
cat("    Naive_B:", length(rnaNaiveB), "cells\n")
cat("    Memory_B:", length(rnaMemoryB), "cells\n")
cat("    Switched_B:", length(rnaSwitchedB), "cells\n")
cat("    Plasma:", length(rnaPlasma), "cells\n")

# ============================================
# 7. Gene integration matrix (unsupervised)
# ============================================

cat("\nStep 6: Unsupervised gene integration...\n")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrixB",
  reducedDims = "IterativeLSIB",
  seRNA = seRNA,
  addToArrow = TRUE,
  force = TRUE,
  sampleCellsATAC = 50000,
  sampleCellsRNA = 50000,
  groupRNA = "celltype_l3",
  nameCell = "predictedCell_Un2",
  nameGroup = "predictedGroup_Un2",
  nameScore = "predictedScore_Un2",
  plotUMAP = TRUE,
  UMAPParams = list(n_neighbors = 35, min_dist = 0.4),
  nGenes = 3000,
  useImputation = TRUE,
  transferParams = list(),
  threads = 6
)

# ============================================
# 8. Group list for supervised integration
# ============================================

cat("\nStep 7: Creating group list for supervised integration...\n")

groupList <- SimpleList(
    Naive_B = SimpleList(
        ATAC = proj$cellNames[proj$BClusters0.6 %in% c('C5')],  
        RNA = rnaNaiveB
    ),
    Memory_B = SimpleList(
        ATAC = proj$cellNames[proj$BClusters0.6 %in% c('C8')],  
        RNA = rnaMemoryB 
    ),
    Switched_B = SimpleList(
        ATAC = proj$cellNames[proj$BClusters0.6 %in% c('C6','C7','C9','C10')],  
        RNA = rnaSwitchedB 
    ),
    Plasma = SimpleList(
        ATAC = proj$cellNames[proj$BClusters0.6 %in% c('C17')],  
        RNA = rnaPlasma
    )
)

cat("  Group list created\n")

# ============================================
# 9. Gene integration matrix (supervised)
# ============================================

cat("\nStep 8: Supervised gene integration...\n")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = TRUE,
  groupList = groupList,
  groupRNA = "celltype_l3",
  nameCell = "predictedCell_new",
  nameGroup = "predictedGroup_new",
  nameScore = "predictedScore_new",
  sampleCellsATAC = 50000,
  sampleCellsRNA = 50000,
  threads = 8
)

# ============================================
# 10. Confusion matrix
# ============================================

cat("\nStep 9: Generating confusion matrix...\n")

cM <- as.matrix(confusionMatrix(proj$BClusters0.6, proj$predictedGroup_new))
preClust <- colnames(cM)[apply(cM, 1, which.max)]
confusion_result <- cbind(preClust, rownames(cM))
colnames(confusion_result) <- c("Predicted", "Cluster")
print(confusion_result)

# ============================================
# 11. Prediction score visualization
# ============================================

cat("\nStep 10: Plotting prediction scores...\n")

p2 <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "predictedScore_new",
  embedding = "UMAPHarmonyB_50_0.4"
)

plotPDF(p2, name = "predictedScore_new.pdf", ArchRProj = proj, addDOC = FALSE, width = 14, height = 12)

# ============================================
# 12. B cell type annotation
# ============================================

cat("\nStep 11: Annotating B cell types...\n")

original_labels <- getCellColData(proj, "predictedGroup_Un2")[, 1]

new_labels <- ifelse(
  original_labels == "pre-Switched_memory_B",
  "Switched_memory_B",
  as.character(original_labels)
)

proj <- addCellColData(
  ArchRProj = proj,
  data = new_labels,
  name = "B_celltype",
  cells = getCellNames(proj)
)

cat("  B cell type distribution:\n")
print(table(proj$B_celltype))

# ============================================
# 13. DBSCAN outlier removal
# ============================================

cat("\nStep 12: DBSCAN outlier removal...\n")

umap_coords <- getEmbedding(proj, embedding = "UMAPHarmonyB_50_0.4")

pdf("kNNdist_plot.pdf", width = 6, height = 4)
dbscan::kNNdistplot(umap_coords, k = 10)
abline(h = 0.5, col = "red", lty = 2)
dev.off()

kNNdist <- dbscan::kNNdist(umap_coords, k = 10)
eps_auto <- quantile(kNNdist, 0.95)
cat("  Auto-selected eps:", eps_auto, "\n")

clusters <- dbscan(umap_coords, eps = eps_auto, minPts = 10)$cluster
proj$B_DBSCAN <- paste0("D", clusters)
cat("  DBSCAN cluster distribution:\n")
print(table(proj$B_DBSCAN))

clusters_to_keep <- unique(proj$B_DBSCAN[proj$B_DBSCAN != "D0"])
cells_to_keep <- which(proj$B_DBSCAN %in% clusters_to_keep)
proj_filtered <- proj[cells_to_keep, ]

cat("  Cells before filtering:", length(proj$cellNames), "\n")
cat("  Cells after filtering:", length(proj_filtered$cellNames), "\n")

# ============================================
# 14. Filtered UMAP visualization
# ============================================

cat("\nStep 13: Generating filtered UMAP...\n")

p_filtered <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "cellColData",
  name = "B_celltype",
  embedding = "UMAPHarmonyB_50_0.4",
  size = 0.01,
  alpha = 0.6,
  plotAs = "points",
  labelMeans = FALSE,
  legendPosition = "right",
  title = "B Cell Type UMAP",
  theme = theme_bw(),
  seed = 90
)

plotPDF(p_filtered, name = "B_predictedGroup_UMAPHarmony.pdf", 
        ArchRProj = proj_filtered, addDOC = FALSE, width = 14, height = 12)

# ============================================
# 15. Save UMAP coordinates
# ============================================

cat("\nStep 14: Saving UMAP coordinates...\n")

umap_coords_filtered <- getEmbedding(ArchRProj = proj_filtered, 
                                      embedding = "UMAPHarmonyB_50_0.4", 
                                      returnDF = TRUE)
cell_annotations <- getCellColData(ArchRProj = proj_filtered, select = "B_celltype")

umap_data <- data.frame(
  Cell_ID = rownames(umap_coords_filtered),
  UMAP_1 = umap_coords_filtered[, 1],
  UMAP_2 = umap_coords_filtered[, 2],
  CellType = cell_annotations$B_celltype
)

write.csv(umap_data, "./umap-data/B-umap.csv", row.names = FALSE)
cat("  UMAP coordinates saved: ./umap-data/B-umap.csv\n")

# ============================================
# 16. B cell marker genes UMAP
# ============================================

cat("\nStep 15: Plotting marker genes...\n")

markerGenes <- c(
  "CD19", "MS4A1", "NEIL1", "CD9", "CD24", "TCL1A", "SOX4",
  "CD27", "IGHD", "IL4R", "FCER2", "BACH2",
  "ITGAX", "FCRL5", "CD1C", "MKI67"
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
  embedding = "UMAPHarmonyB_50_0.4",
  imputeWeights = getImputeWeights(proj_filtered)
)

p2 <- lapply(p, function(x) {
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
})

plotPDF(plotList = p,
  name = "B_Marker_Genes_W_Imputation.pdf",
  ArchRProj = proj_filtered,
  addDOC = FALSE,
  width = 5,
  height = 5)

# ============================================
# 17. B cell subset marker gene heatmap
# ============================================

cat("\nStep 16: Generating B cell subset marker heatmap...\n")

markersGS <- getMarkerFeatures(
  ArchRProj = proj_filtered,
  useMatrix = "GeneScoreMatrix",
  groupBy = "B_celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmap_markerGenes <- c(
  "TCL1A", "IGHD", "CD24", "CD38",
  "CD27", "CD40", "CD80", "CD86",
  "IGHA1", "IGHG1", "AICDA",
  "IGHM",
  "CD11c", "FCRL3", "FCRL5", "CD21", "CD95",
  "CD5", "CD10", "MME",
  "XBP1", "IRF4", "PRDM1", "SDC1"
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
  
  pdf("./results/B_celltype_heatmap.pdf", width = 10, height = 7)
  
  p_heatmap <- ComplexHeatmap::pheatmap(heatmapGS,
    show_rownames = TRUE,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    col = col,
    name = "Z-Scores",
    fontsize_row = 8,
    fontsize_col = 10,
    main = "B Cell Subset Marker Genes"
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
  cat("  Heatmap saved to: ./results/B_celltype_heatmap.pdf\n")
} else {
  cat("  No marker genes found for heatmap\n")
}

# ============================================
# 18. Final summary
# ============================================

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETE!\n")
cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\nFinal summary:\n")
cat("  - Total cells in filtered project:", nrow(getCellColData(proj_filtered)), "\n")
cat("  - B cell types:", paste(unique(proj_filtered$B_celltype), collapse=", "), "\n")
cat("\nOutput files:\n")
cat("  - ./umap-data/B-umap.csv\n")
cat("  - B_Marker_Genes_W_Imputation.pdf\n")
cat("  - B_predictedGroup_UMAPHarmony.pdf\n")
cat("  - predictedScore_new.pdf\n")
cat("  - kNNdist_plot.pdf\n")
cat("  - ./results/B_celltype_heatmap.pdf\n")
cat("=", rep("=", 60), "\n")
