#!/usr/bin/env Rscript

# ============================================
# CellChat Analysis for Cell-Cell Communication
# Compare SED vs EX groups
# ============================================

library(CellChat)
library(patchwork)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(dplyr)

set.seed(123)

# ============================================
# Step 1: Load data
# ============================================

subset_data <- readRDS("/media/AnalysisTempDisk2/Songxuemei/Miol/allcels.rds")

# ============================================
# Step 2: Balance samples by sex
# ============================================

cat("=", rep("=", 60), "\n")
cat("Step 1: Balancing samples by sex\n")
cat("=", rep("=", 60), "\n")

# Calculate minimum sample size per sex
sex_balance_info <- subset_data@meta.data %>%
  distinct(sample, group, sex) %>%
  group_by(sex) %>%
  summarise(
    min_in_sex = min(table(group)),
    .groups = "drop"
  )

cat("Original sample distribution:\n")
sample_counts <- subset_data@meta.data %>%
  distinct(sample, group, sex) %>%
  group_by(group, sex) %>%
  summarise(n_samples = n(), .groups = "drop")
print(sample_counts)

# Balance samples per sex
balanced_samples_list <- list()

for (current_sex in unique(subset_data@meta.data$sex)) {
  current_info <- sex_balance_info %>% filter(sex == current_sex)
  
  sex_samples <- subset_data@meta.data %>%
    filter(sex == current_sex) %>%
    distinct(sample, group)
  
  if (n_distinct(sex_samples$group) == 1 || 
      min(table(sex_samples$group)) == max(table(sex_samples$group))) {
    selected_samples <- sex_samples$sample
  } else {
    selected_samples <- sex_samples %>%
      group_by(group) %>%
      group_modify(~ {
        if (nrow(.x) > current_info$min_in_sex) {
          .x %>% sample_n(current_info$min_in_sex)
        } else {
          .x
        }
      }) %>%
      pull(sample)
  }
  
  balanced_samples_list[[current_sex]] <- selected_samples
}

balanced_samples <- unlist(balanced_samples_list)

# Create balanced Seurat object
seurat_balanced <- subset(subset_data, subset = sample %in% balanced_samples)

# Verify balance
final_sample_dist <- seurat_balanced@meta.data %>%
  distinct(sample, group, sex) %>%
  group_by(group, sex) %>%
  summarise(n_samples = n(), .groups = "drop")

cat("\nBalanced sample distribution:\n")
print(final_sample_dist)
cat("Total samples:", sum(final_sample_dist$n_samples), "\n")

# ============================================
# Step 3: Normalize data
# ============================================

cat("\n", rep("=", 60), "\n")
cat("Step 2: Normalizing data\n")
cat("=", rep("=", 60), "\n")

seurat_balanced <- NormalizeData(seurat_balanced)

# ============================================
# Step 4: Create CellChat objects for EX and SED
# ============================================

cat("\n", rep("=", 60), "\n")
cat("Step 3: Creating CellChat objects\n")
cat("=", rep("=", 60), "\n")

# EX group
data_ex <- GetAssayData(seurat_balanced, assay = "RNA", slot = "data")[, 
                seurat_balanced@meta.data$group == "EX"]
meta_ex <- seurat_balanced@meta.data[seurat_balanced@meta.data$group == "EX", ]
cellchat_ex <- createCellChat(object = data_ex, meta = meta_ex, group.by = "cell_type_l1")

# SED group
data_sed <- GetAssayData(seurat_balanced, assay = "RNA", slot = "data")[, 
                seurat_balanced@meta.data$group == "SED"]
meta_sed <- seurat_balanced@meta.data[seurat_balanced@meta.data$group == "SED", ]
cellchat_sed <- createCellChat(object = data_sed, meta = meta_sed, group.by = "cell_type_l1")

cat("EX group cells:", ncol(data_ex), "\n")
cat("SED group cells:", ncol(data_sed), "\n")

# ============================================
# Step 5: Set database
# ============================================

cat("\n", rep("=", 60), "\n")
cat("Step 4: Setting CellChat database\n")
cat("=", rep("=", 60), "\n")

CellChatDB <- CellChatDB.human  # Use CellChatDB.mouse for mouse data

cellchat_ex@DB <- CellChatDB
cellchat_sed@DB <- CellChatDB

showDatabaseCategory(CellChatDB)

# ============================================
# Step 6: Run CellChat for EX group
# ============================================

cat("\n", rep("=", 60), "\n")
cat("Step 5: Running CellChat for EX group\n")
cat("=", rep("=", 60), "\n")

cellchat_ex <- subsetData(cellchat_ex)
cellchat_ex <- identifyOverExpressedGenes(cellchat_ex)
cellchat_ex <- identifyOverExpressedInteractions(cellchat_ex)
cellchat_ex <- computeCommunProb(cellchat_ex, 
                                 type = "truncatedMean", 
                                 trim = 0.1,
                                 population.size = TRUE)
cellchat_ex <- filterCommunication(cellchat_ex, min.cells = 10)
cellchat_ex <- computeCommunProbPathway(cellchat_ex)
cellchat_ex <- aggregateNet(cellchat_ex)

cat("EX group analysis complete\n")

# ============================================
# Step 7: Run CellChat for SED group
# ============================================

cat("\n", rep("=", 60), "\n")
cat("Step 6: Running CellChat for SED group\n")
cat("=", rep("=", 60), "\n")

cellchat_sed <- subsetData(cellchat_sed)
cellchat_sed <- identifyOverExpressedGenes(cellchat_sed)
cellchat_sed <- identifyOverExpressedInteractions(cellchat_sed)
cellchat_sed <- computeCommunProb(cellchat_sed, 
                                 type = "truncatedMean", 
                                 trim = 0.1,
                                 population.size = TRUE)
cellchat_sed <- filterCommunication(cellchat_sed, min.cells = 10)
cellchat_sed <- computeCommunProbPathway(cellchat_sed)
cellchat_sed <- aggregateNet(cellchat_sed)

cat("SED group analysis complete\n")

# ============================================
# Step 8: Merge and compare
# ============================================

cat("\n", rep("=", 60), "\n")
cat("Step 7: Merging CellChat objects\n")
cat("=", rep("=", 60), "\n")

cellchat_list <- list(SED = cellchat_sed, EX = cellchat_ex)
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

cat("Merge complete\n")

# ============================================
# Step 9: Save objects
# ============================================

cat("\n", rep("=", 60), "\n")
cat("Step 8: Saving results\n")
cat("=", rep("=", 60), "\n")

saveRDS(cellchat_sed, "cellchat_SED.rds")
saveRDS(cellchat_ex, "cellchat_EX.rds")
saveRDS(cellchat_merged, "cellchat_merged.rds")

cat("Objects saved:\n")
cat("  - cellchat_SED.rds\n")
cat("  - cellchat_EX.rds\n")
cat("  - cellchat_merged.rds\n")

# ============================================
# Step 10: Basic visualization examples
# ============================================

cat("\n", rep("=", 60), "\n")
cat("Step 9: Generating basic plots\n")
cat("=", rep("=", 60), "\n")

# Number of interactions comparison
p1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2))
ggsave("compare_interactions_number.pdf", p1, width = 5, height = 4)

# Interaction strength comparison
p2 <- compareInteractions(cellchat_merged, measure = "weight", show.legend = FALSE, group = c(1,2))
ggsave("compare_interactions_weight.pdf", p2, width = 5, height = 4)

# Circle plot for SED
p3 <- netVisual_circle(cellchat_sed@net$count,
                       vertex.weight = as.numeric(table(cellchat_sed@idents)),
                       weight.scale = T,
                       label.edge = F,
                       title.name = "Number of interactions - SED")
ggsave("SED_circle_plot.pdf", p3, width = 8, height = 7)

# Circle plot for EX
p4 <- netVisual_circle(cellchat_ex@net$count,
                       vertex.weight = as.numeric(table(cellchat_ex@idents)),
                       weight.scale = T,
                       label.edge = F,
                       title.name = "Number of interactions - EX")
ggsave("EX_circle_plot.pdf", p4, width = 8, height = 7)

cat("Plots saved:\n")
cat("  - compare_interactions_number.pdf\n")
cat("  - compare_interactions_weight.pdf\n")
cat("  - SED_circle_plot.pdf\n")
cat("  - EX_circle_plot.pdf\n")

# ============================================
# Final summary
# ============================================

cat("\n", rep("=", 60), "\n")
cat("CellChat Analysis Complete!\n")
cat("=", rep("=", 60), "\n")
cat("\nOutput files:\n")
cat("  - cellchat_SED.rds\n")
cat("  - cellchat_EX.rds\n")
cat("  - cellchat_merged.rds\n")
cat("  - compare_interactions_number.pdf\n")
cat("  - compare_interactions_weight.pdf\n")
cat("  - SED_circle_plot.pdf\n")
cat("  - EX_circle_plot.pdf\n")

#!/usr/bin/env Rscript

# ============================================
# CellChat Visualization
# Heatmaps, comparison plots, and rank plots
# ============================================

library(CellChat)
library(ComplexHeatmap)
library(patchwork)
library(ggplot2)
library(dplyr)

# ============================================
# Load merged CellChat object
# ============================================

cellchat <- readRDS("cellchat_merged.rds")

cat("=", rep("=", 60), "\n")
cat("CellChat Visualization\n")
cat("=", rep("=", 60), "\n\n")

# ============================================
# 1. Pathway union for heatmaps
# ============================================

cat("Step 1: Creating pathway union\n")

# Get cellchat list from merged object
cellchat_list <- list(SED = cellchat@netP$pathways, EX = cellchat@netP$pathways)

# Union of pathways from both groups
pathway.union <- union(cellchat@netP$pathways, cellchat@netP$pathways)

cat("Number of pathways in union:", length(pathway.union), "\n")

# ============================================
# 2. Signaling role heatmaps (incoming and outgoing)
# ============================================

cat("\nStep 2: Generating signaling role heatmaps\n")

# Incoming signaling heatmap - SED
ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming",
  signaling = pathway.union,
  title = "Incoming - SED",
  width = 10,
  height = 17
)

# Outgoing signaling heatmap - SED
ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing",
  signaling = pathway.union,
  title = "Outgoing - SED",
  width = 10,
  height = 17
)

# Save combined heatmap
pdf("signaling_role_heatmap_SED.pdf", width = 12, height = 18)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))
dev.off()
cat("  Saved: signaling_role_heatmap_SED.pdf\n")

# For EX group
ht3 <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming",
  signaling = pathway.union,
  title = "Incoming - EX",
  width = 10,
  height = 17
)

ht4 <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing",
  signaling = pathway.union,
  title = "Outgoing - EX",
  width = 10,
  height = 17
)

pdf("signaling_role_heatmap_EX.pdf", width = 12, height = 18)
draw(ht3 + ht4, ht_gap = unit(1, "cm"))
dev.off()
cat("  Saved: signaling_role_heatmap_EX.pdf\n")

# ============================================
# 3. Compare interactions (number and weight)
# ============================================

cat("\nStep 3: Comparing interactions\n")

p1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
p2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")

pdf("compare_interactions.pdf", width = 8, height = 6)
print(p1 + p2)
dev.off()
cat("  Saved: compare_interactions.pdf\n")

# ============================================
# 4. Differential interaction networks
# ============================================

cat("\nStep 4: Differential interaction networks\n")

pdf("diff_interaction_networks.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()
cat("  Saved: diff_interaction_networks.pdf\n")

# ============================================
# 5. Heatmaps for interaction strength
# ============================================

cat("\nStep 5: Interaction strength heatmaps\n")

p3 <- netVisual_heatmap(cellchat)
p4 <- netVisual_heatmap(cellchat, measure = "weight")

pdf("interaction_heatmaps.pdf", width = 12, height = 8)
print(p3 + p4)
dev.off()
cat("  Saved: interaction_heatmaps.pdf\n")

# ============================================
# 6. Rank plots for signaling pathways
# ============================================

cat("\nStep 6: Rank plots for signaling pathways\n")

p5 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
p6 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)

pdf("rank_net_plots.pdf", width = 14, height = 8)
print(p5 + p6)
dev.off()
cat("  Saved: rank_net_plots.pdf\n")



# ============================================
# 7. Summary
# ============================================

cat("\n", rep("=", 60), "\n")
cat("Visualization Complete!\n")
cat("=", rep("=", 60), "\n")
cat("\nOutput files:\n")
cat("  - signaling_role_heatmap_SED.pdf\n")
cat("  - signaling_role_heatmap_EX.pdf\n")
cat("  - compare_interactions.pdf\n")
cat("  - diff_interaction_networks.pdf\n")
cat("  - interaction_heatmaps.pdf\n")
cat("  - rank_net_plots.pdf\n")



