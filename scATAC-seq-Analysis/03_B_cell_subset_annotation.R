#!/usr/bin/env Rscript

# ============================================
# B cell subset analysis: LSI + Harmony + Integration
# ============================================

library(ArchR)
library(dbscan)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 16)
addArchRGenome("hg38")

# ============================================
# 1. Load project
# ============================================

proj <- loadArchRProject("/media/AnalysisTempDisk2/Songxuemei/Subset_B/")
seRNA <- readRDS("/media/AnalysisTempDisk2/Songxuemei/运动/B-drop.rds")

# ============================================
# 2. LSI dimensionality reduction
# ============================================

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

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSIB",
  method = "Seurat",
  name = "BClusters0.6",
  resolution = 0.6
)

# ============================================
# 4. Harmony batch correction
# ============================================

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "HarmonyB",
  groupBy = "newSample"
)

# ============================================
# 5. UMAP visualization
# ============================================

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
# 6. Gene integration matrix (unsupervised)
# ============================================

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
# 7. Group list for supervised integration
# ============================================

groupList <- SimpleList(
  ncMono = SimpleList(
    ATAC = proj$cellNames[proj$BClusters0.6 %in% c('C5')],
    RNA = ncMono
  ),
  pDC = SimpleList(
    ATAC = proj$cellNames[proj$BClusters0.6 %in% c('C8')],
    RNA = rnapDC
  ),
  Plasma = SimpleList(
    ATAC = proj$cellNames[proj$BClusters0.6 %in% c('C6', 'C7', 'C9', 'C10')],
    RNA = rnaPlasma
  ),
  Cmo = SimpleList(
    ATAC = proj$cellNames[proj$BClusters0.6 %in% c('C17')],
    RNA = rnaMo
  ),
  Macro = SimpleList(
    ATAC = proj$cellNames[proj$BClusters0.6 %in% c('C1', 'C2', 'C3', 'C4', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C18', 'C19')],
    RNA = rnamacro
  )
)

# ============================================
# 8. Gene integration matrix (supervised)
# ============================================

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
# 9. Confusion matrix
# ============================================

cM <- as.matrix(confusionMatrix(proj$BClusters0.6, proj$predictedGroup_Co))
preClust <- colnames(cM)[apply(cM, 1, which.max)]
cbind(preClust, rownames(cM))

# ============================================
# 10. Prediction score visualization
# ============================================

p2 <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "cellColData",
  name = "predictedScore_Co",
  embedding = "UMAPHarmonyB_50_0.4"
)
plotPDF(p2, name = "predictedScore_Co.pdf", ArchRProj = proj, addDOC = FALSE, width = 14, height = 12)

# ============================================
# 11. B cell type annotation
# ============================================

original_labels <- getCellColData(proj, "predictedGroup_UnB")[, 1]

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

table(proj$B_celltype)

# ============================================
# 12. DBSCAN outlier removal
# ============================================

umap_coords <- getEmbedding(proj, embedding = "UMAPHarmonyB_50_0.4")
dbscan::kNNdistplot(umap_coords, k = 10)

clusters <- dbscan(umap_coords, eps = 0.5, minPts = 10)$cluster
proj$B_DBSCAN <- paste0("D", clusters)
table(proj$B_DBSCAN)

clusters_to_keep <- c("D1", "D2")
cells_to_keep <- which(proj$B_DBSCAN %in% clusters_to_keep)
proj_filtered <- proj[cells_to_keep, ]

# ============================================
# 13. Filtered UMAP visualization
# ============================================

p_filtered <- plotEmbedding(
  ArchRProj = proj_filtered,
  colorBy = "cellColData",
  name = "B_celltype",
  embedding = "UMAPHarmonyB_50_0.4",
  size = 0.01,
  alpha = 0.6,
  palette = viridis::viridis(10),
  plotAs = "points",
  labelMeans = FALSE,
  legendPosition = "right",
  title = "B Cell Type UMAP",
  theme = theme_bw(),
  seed = 90
) + ggplot2::theme(legend.text = element_text(size = 12))

plotPDF(p_filtered, name = "B-predictedGroup_Un-UMAPHarmony.pdf", ArchRProj = proj_filtered, addDOC = FALSE, width = 14, height = 12)

# ============================================
# 14. Save UMAP coordinates
# ============================================

umap_coords <- getEmbedding(ArchRProj = proj_filtered, embedding = "UMAPHarmonyB_50_0.4", returnDF = TRUE)
cell_annotations <- getCellColData(ArchRProj = proj_filtered, select = "B_celltype")

umap_data <- data.frame(
  Cell_ID = rownames(umap_coords),
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  CellType = cell_annotations$B_celltype
)

write.csv(umap_data, "./umap-data/B-umap.csv", row.names = FALSE)

# ============================================
# 15. B cell marker genes UMAP
# ============================================

markerGenes <- c(
  "CD19", "MS4A1", "NEIL1", "CD9", "CD24", "TCL1A", "SOX4",
  "CD27", "IGHD", "IL4R", "FCER2", "BACH2",
  "ITGAX", "FCRL5", "CD1C", "MKI67"
)

all_genes <- getFeatures(proj, useMatrix = "GeneScoreMatrix")
existing_markerGenes <- markerGenes[markerGenes %in% all_genes]
missing_genes <- markerGenes[!markerGenes %in% all_genes]

if (length(missing_genes) > 0) {
  cat("Missing genes:", paste(missing_genes, collapse = ", "), "\n")
}

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
do.call(cowplot::plot_grid, c(list(ncol = 3), p2))

plotPDF(plotList = p,
  name = "1Plot-B-Marker-Genes-W-Imputation.pdf",
  ArchRProj = proj_filtered,
  addDOC = FALSE,
  width = 5,
  height = 5)

cat("Analysis complete!\n")
