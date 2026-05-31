#!/usr/bin/env Rscript

# ============================================
#  scATAC-seq: LSI + Harmony + Clustering
# ============================================

library(ArchR)
library(parallel)
library(dbscan)
set.seed(1)

# ============================================
# 1. Load ArchR project
# ============================================

proj <- loadArchRProject("./ArchRProject/")

# ============================================
# 2. LSI dimensionality reduction
# ============================================

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 4,
  clusterParams = list(
    resolution = c(0.1, 0.2, 0.4, 0.8),
    sampleCells = 30000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  force = TRUE
)

# ============================================
# 3. Harmony batch correction
# ============================================

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "newSample",
  force = TRUE
)

# ============================================
# 4. Clustering (major cell types)
# ============================================

proj <- addClusters(
  input = proj,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  maxClusters = 40,
  force = TRUE
)

cat("Clusters:\n")
print(table(proj$Clusters))

# ============================================
# 5. UMAP visualization
# ============================================

proj <- addUMAP(
  ArchRProj = proj,
  reducedDims = "Harmony",
  name = "UMAPHarmony",
  nNeighbors = 40,
  minDist = 0.35,
  metric = "cosine",
  force = TRUE
)

# Plot UMAP
p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "newSample", embedding = "UMAPHarmony")
p3 <- plotEmbedding(proj, colorBy = "cellColData", name = "exercise_group", embedding = "UMAPHarmony")

plotPDF(p1, p2, p3, name = "01_UMAP_Overview.pdf", ArchRProj = proj, addDOC = FALSE, width = 8, height = 7)

# ============================================
# 6. Remove outliers (DBSCAN)
# ============================================

umap_coords <- getEmbedding(proj, embedding = "UMAPHarmony")
db_result <- dbscan(umap_coords, eps = 0.5, minPts = 10)
proj$DBSCAN <- paste0("D", db_result$cluster)

cat("DBSCAN results:\n")
print(table(proj$DBSCAN))

# Keep only non-outlier cells
cells_keep <- proj$cellNames[proj$DBSCAN != "D0"]
proj <- proj[cells_keep, ]
cat("Cells after outlier removal:", nCells(proj), "\n")

# ============================================
# 7. Add imputation weights
# ============================================

proj <- addImputeWeights(proj)

# ============================================
# 8. Module scoring for cell type annotation
# ============================================

features <- list(
  B = c("MS4A1", "CD79A", "CD79B", "CD74"),
  CD4_T = c("CD4", "CCR7", "LEF1", "TCF7", "SELL"),
  CD8_T = c("CD8A", "CD8B", "GZMB", "PRF1"),
  NK = c("NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1"),
  Myeloid = c("CD14", "LYZ", "FCGR3A", "S100A8", "S100A9"),
  Unconv_T = c("TRGC1", "TRDV1", "SLC4A10","RORC")
)

proj <- addModuleScore(
  proj,
  useMatrix = "GeneScoreMatrix",
  name = "Module",
  features = features,
  force = TRUE
)

# Plot module scores
module_names <- paste0("Module.", names(features))
p_mod <- lapply(module_names, function(m) {
  plotEmbedding(proj, embedding = "UMAPHarmony",
                colorBy = "cellColData", name = m,
                imputeWeights = getImputeWeights(proj))
})
plotPDF(p_mod, name = "02_Module_Scores.pdf", ArchRProj = proj, addDOC = FALSE, width = 6, height = 5)

# ============================================
# 9. Marker genes
# ============================================

markerGenes <- c(
  "MS4A1", "CD79A", "CD79B",  # B
  "CD4", "CCR7", "LEF1", "TCF7",  # CD4_T
  "CD8A", "CD8B", "GZMB", "PRF1",  # CD8_T
  "NKG7", "GNLY", "KLRD1", "NCAM1",  # NK
  "CD14", "LYZ", "FCGR3A", "S100A8", "S100A9",  # Myeloid
  "TRGC1", "TRDV1", "SLC4A10","RORC"  # Unconv_T
)

all_genes <- getFeatures(proj, useMatrix = "GeneScoreMatrix")
markerGenes_ok <- markerGenes[markerGenes %in% all_genes]

p_genes <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "GeneScoreMatrix",
  name = markerGenes_ok,
  embedding = "UMAPHarmony",
  imputeWeights = getImputeWeights(proj)
)

plotPDF(p_genes, name = "03_Marker_Genes.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# ============================================
# 10. Cell type annotation
# ============================================

# Example annotation (modify according to your clusters)
cluster_annotation <- c(
  "C1" = "CD8_T", "C2" = "CD4_T", "C3" = "CD8_T",
  "C4" = "B", "C5" = "NK", "C6" = "Myeloid",
  "C7" = "CD4_T", "C8" = "CD8_T", "C9" = "Unconv_T",
  "C10" = "Myeloid", "C11" = "B", "C12" = "NK",
  "C13" = "CD8_T", "C14" = "CD4_T", "C15" = "Unconv_T"
)

proj$cell_type <- cluster_annotation[proj$Clusters]

p_anno <- plotEmbedding(
  proj, colorBy = "cellColData", name = "cell_type",
  embedding = "UMAPHarmony", labelMeans = TRUE
)
plotPDF(p_anno, name = "04_Cell_Type_Annotation.pdf", ArchRProj = proj, addDOC = FALSE, width = 8, height = 7)

# ============================================
# 11. Save project
# ============================================

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "./ArchRProject_annotated/",
  overwrite = TRUE
)

cat("Analysis complete!\n")
