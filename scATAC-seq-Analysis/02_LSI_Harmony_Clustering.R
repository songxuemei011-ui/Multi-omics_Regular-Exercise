#!/usr/bin/env Rscript

# ============================================
# scATAC-seq: LSI + Harmony + Clustering
# ============================================

library(ArchR)
library(parallel)
library(dbscan)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
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
  Unconv_T = c("TRGC1", "TRDV1", "SLC4A10", "RORC")
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
# 9. Cell type annotation
# ============================================

# Modify according to your clusters
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
plotPDF(p_anno, name = "03_Cell_Type_Annotation.pdf", ArchRProj = proj, addDOC = FALSE, width = 8, height = 7)

# ============================================
# 10. Marker genes UMAP
# ============================================

markerGenes <- c(
  "MS4A1", "CD79A", "CD79B",  # B
  "CD4", "CCR7", "LEF1", "TCF7",  # CD4_T
  "CD8A", "CD8B", "GZMB", "PRF1",  # CD8_T
  "NKG7", "GNLY", "KLRD1", "NCAM1",  # NK
  "CD14", "LYZ", "FCGR3A", "S100A8", "S100A9",  # Myeloid
  "TRGC1", "TRDV1", "SLC4A10", "RORC"  # Unconv_T
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

plotPDF(p_genes, name = "04_Marker_Genes_UMAP.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# ============================================
# 11. L1 marker gene heatmap
# ============================================

markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "cell_type",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Marker genes for heatmap
heatmap_markerGenes <- c(
  "MS4A1", "CD19", "CD79A",   # B
  "CD4", "IL7R", "FOXP3",     # CD4_T
  "CD8A", "GZMB", "CD69",     # CD8_T
  "CD14", "FCGR3A", "LYZ",    # Myeloid
  "NKG7", "GNLY", "KLRD1",    # NK
  "TRDC", "KLRB1", "FCER1G", "RORC", "RORA"   # Unconv_T
)

# Check which genes exist
all_genes <- getFeatures(proj, useMatrix = "GeneScoreMatrix")
existing_heatmap_genes <- heatmap_markerGenes[heatmap_markerGenes %in% all_genes]
missing_genes <- heatmap_markerGenes[!heatmap_markerGenes %in% all_genes]

if (length(missing_genes) > 0) {
  cat("Missing genes in heatmap:", paste(missing_genes, collapse = ", "), "\n")
} else {
  cat("All heatmap marker genes found!\n")
}

# Extract marker list
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")

# Combine all markers
cell_types <- names(markerList)
all_markers <- data.frame()

for (ct in cell_types) {
  markers <- markerList[[ct]]
  markers_df <- data.frame(
    name = markers$name,
    Log2FC = markers$Log2FC,
    cell_type = ct
  )
  all_markers <- rbind(all_markers, markers_df)
}

# Handle duplicate genes (keep highest Log2FC)
all_markers <- all_markers[order(all_markers$cell_type, -all_markers$Log2FC),]
all_markers$rank <- 1:nrow(all_markers)

gene_freq <- as.data.frame(table(all_markers$name))
dup_genes <- as.character(gene_freq[gene_freq$Freq > 1, "Var1"])

markerList_unique <- all_markers[!(all_markers$name %in% dup_genes),]
markerList_dup <- all_markers[all_markers$name %in% dup_genes,]
markerList_dup <- markerList_dup %>% group_by(name) %>% filter(Log2FC == max(Log2FC))
all_markers <- rbind(markerList_unique, markerList_dup)
all_markers <- all_markers[order(all_markers$rank),]

# Plot heatmap
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  labelMarkers = existing_heatmap_genes,
  binaryClusterRows = TRUE,
  limits = c(-2, 2),
  transpose = TRUE,
  returnMatrix = TRUE
)

# Filter and transpose
heatmapGS <- heatmapGS[, unique(all_markers$name)]
heatmapGS <- t(heatmapGS)

# Draw heatmap
col <- paletteContinuous(set = "blueYellow")

pdf("./results/celltype_L1_heatmap.pdf", width = 6, height = 8)
p <- ComplexHeatmap::pheatmap(heatmapGS,
  show_rownames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  col = col,
  name = paste0("Column Z-Scores\n", nrow(heatmapGS), " features\n", "GeneScoreMatrix")
)
p + rowAnnotation(link = anno_mark(
  at = which(rownames(heatmapGS) %in% existing_heatmap_genes),
  labels = rownames(heatmapGS)[which(rownames(heatmapGS) %in% existing_heatmap_genes)],
  labels_gp = gpar(fontsize = 10)
))
dev.off()

cat("Heatmap saved to: ./results/celltype_L1_heatmap.pdf\n")

# ============================================
# 12. Save project
# ============================================

saveArchRProject(
  ArchRProj = proj,
  outputDirectory = "./ArchRProject_annotated/",
  overwrite = TRUE
)

cat("Analysis complete!\n")
