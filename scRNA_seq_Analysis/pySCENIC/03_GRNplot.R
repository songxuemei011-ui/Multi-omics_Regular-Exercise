# ============================================
# CD8_CTL TF Analysis
#   1. RUNX3 Regulatory Network
#   2. AUCell Heatmap for selected TFs
# ============================================

library(readr)
library(dplyr)
library(ggraph)
library(igraph)
library(ggplot2)
library(reshape2)
library(AUCell)
library(loomR)

# ============================================
# USER SETTINGS - MODIFY THESE PATHS
# ============================================

# Path to GRN results
GRN_FILE <- "./CD8_CTL_scenic_grn.csv"

# Path to AUCell loom file
AUC_LOOM_FILE <- "/media/AnalysisTempDisk2/Songxuemei/运动/pyscenic/CD8CTL_scenic_auc.loom"

# Group information (modify according to your data)
GROUP_INFO <- c("SED", "EX")  # or read from metadata

# ============================================
# Part 1: RUNX3 Regulatory Network (CD8_CTL)
# ============================================

cat("=", rep("=", 60), "\n")
cat("Part 1: RUNX3 Regulatory Network (CD8_CTL)\n")
cat("=", rep("=", 60), "\n")

df <- read_csv(GRN_FILE)

selected_immune_genes <- c(
  "IRF1","PARP1","LCK","MICB","SMARCE1","TPI1","CRIP1",
  "VPS26A","CYTH2","KDM2B","SETBP1","CKAP2","BCL11B",
  "PIK3R2","PKM","JMJD1C","CERK","TUBB","UCP2","SELPLG",
  "GZMB","TBX21","TNFRSF1B","SMAD7","CD8A","ID2",
  "KLF6","ELF4","MAPK9","CD226","EFHD2"
)

filtered_df_final <- df %>%
  filter(target %in% selected_immune_genes) %>%
  filter(importance > 0.1) %>%
  filter(TF == "RUNX3")

if(nrow(filtered_df_final) > 0) {
  
  edges <- filtered_df_final %>%
    rename(from = TF, to = target, weight = importance)
  
  target_nodes <- unique(edges$to)
  nodes <- data.frame(
    name = c("RUNX3", target_nodes),
    type = c("TF", rep("Target", length(target_nodes)))
  )
  
  gr <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
  
  set.seed(42)
  
  p1 <- ggraph(gr, layout = "fr") +
    geom_edge_link(
      arrow = arrow(length = unit(2, "mm"), type = "closed"),
      end_cap = circle(7, "mm"),
      start_cap = circle(5, "mm"),
      edge_colour = "grey55",
      edge_width = 0.45,
      edge_alpha = 0.75
    ) +
    geom_node_point(
      data = function(x) subset(x, type == "Target"),
      aes(x = x, y = y),
      shape = 21, size = 14,
      fill = "#7ECCE8", colour = "grey70", stroke = 0.4
    ) +
    geom_node_point(
      data = function(x) subset(x, type == "TF"),
      aes(x = x, y = y),
      shape = 22, size = 20,
      fill = "#F5E17A", colour = "grey55", stroke = 1
    ) +
    geom_node_text(
      data = function(x) subset(x, type == "Target"),
      aes(x = x, y = y, label = name),
      size = 2.5, colour = "black", fontface = "plain"
    ) +
    geom_node_text(
      data = function(x) subset(x, type == "TF"),
      aes(x = x, y = y, label = name),
      size = 4.5, colour = "black", fontface = "bold"
    ) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    coord_fixed()
  
  print(p1)
  ggsave("CD8_CTL_RUNX3_network.pdf", p1, width = 9, height = 9)
  ggsave("CD8_CTL_RUNX3_network.png", p1, width = 9, height = 9, dpi = 300)
  cat("Saved: CD8_CTL_RUNX3_network.pdf\n")
  
} else {
  cat("No data found for RUNX3 in CD8_CTL\n")
}

# ============================================
# Part 2: AUCell Heatmap for CD8_CTL
# ============================================

cat("\n", "=", rep("=", 60), "\n")
cat("Part 2: AUCell Heatmap (CD8_CTL)\n")
cat("=", rep("=", 60), "\n")

# Read AUCell results from loom file
cat("Reading AUCell loom file:", AUC_LOOM_FILE, "\n")
loom <- open_loom(AUC_LOOM_FILE)

# Extract AUC matrix (regulons x cells)
auc_matrix <- as.matrix(loom[["/AUC"]][, ])
rownames(auc_matrix) <- loom[["/row_attrs/Regulons"]][]
colnames(auc_matrix) <- loom[["/col_attrs/CellID"]][]

close_loom(loom)

# Transpose to cells x regulons
auc_df <- as.data.frame(t(auc_matrix))
auc_df$group <- GROUP_INFO  # Replace with actual group info from your data

# Calculate mean AUC per group
auc_avg <- auc_df %>%
  group_by(group) %>%
  summarise(across(where(is.numeric), mean)) %>%
  arrange(match(group, c("SED", "EX")))

auc_avg_mat <- as.data.frame(auc_avg)
rownames(auc_avg_mat) <- auc_avg_mat$group
auc_avg_mat <- auc_avg_mat[, -1]

# Select TFs of interest
top_TFs <- c("RUNX3", "FOSB", "KLF6", "EOMES", "TBX21", "ETS1",
             "IRF7", "STAT1", "STAT2", "FOS", "TBX3", "SOX4", "ZEB2")

# Match with available regulons (remove (+) suffix if needed)
colnames_clean <- gsub("\\(\\+\\)", "", colnames(auc_avg_mat))
top_TFs_found <- top_TFs[top_TFs %in% colnames_clean]
cat("Found TFs in CD8_CTL:", paste(top_TFs_found, collapse = ", "), "\n")

# Subset and rename columns
auc_avg_sub <- auc_avg_mat[, gsub("\\(\\+\\)", "", colnames(auc_avg_mat)) %in% top_TFs_found, drop = FALSE]
colnames(auc_avg_sub) <- gsub("\\(\\+\\)", "", colnames(auc_avg_sub))

# Keep order
auc_avg_sub <- auc_avg_sub[, top_TFs_found[top_TFs_found %in% colnames(auc_avg_sub)], drop = FALSE]

# Z-score
auc_scaled <- as.data.frame(scale(auc_avg_sub))
auc_scaled$Group <- rownames(auc_scaled)

# Melt for ggplot
auc_long <- melt(as.matrix(auc_scaled[, -ncol(auc_scaled)]),
                 varnames = c("Group", "TF"),
                 value.name = "Activity")
auc_long$Group <- factor(auc_long$Group, levels = c("SED", "EX"))
auc_long$TF <- factor(auc_long$TF, levels = top_TFs_found)

# Heatmap
p2 <- ggplot(auc_long, aes(x = TF, y = Group, fill = Activity)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "#99CC00",      # green
    mid = "white",
    high = "#FF9999",     # pink
    midpoint = 0,
    name = "Z-score",
    guide = guide_colorbar(barwidth = 0.8, barheight = 4)
  ) +
  labs(x = "", y = "", title = "CD8_CTL AUCell Heatmap") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.background = element_rect(fill = "white", colour = NA)
  )

print(p2)
ggsave("CD8_CTL_AUC_heatmap.pdf", p2, width = 10, height = 2.5)
ggsave("CD8_CTL_AUC_heatmap.png", p2, width = 10, height = 2.5, dpi = 300)
cat("Saved: CD8_CTL_AUC_heatmap.pdf\n")

cat("\n", "=", rep("=", 60), "\n")
cat("CD8_CTL TF Analysis Complete!\n")
cat("=", rep("=", 60), "\n")
