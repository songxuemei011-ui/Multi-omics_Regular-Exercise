library(readr)
library(dplyr)
library(ggraph)
library(igraph)
library(ggplot2)

# ============================================
# USER SETTINGS - MODIFY THIS PATH
# ============================================

# Change this to your file path
INPUT_FILE <- "./CD8_CTL_scenic_grn.csv"

# ============================================
# DO NOT MODIFY BELOW
# ============================================

# Read data
df <- read_csv(INPUT_FILE)

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

# Build graph
edges <- filtered_df_final %>%
  rename(from = TF, to = target, weight = importance)

target_nodes <- unique(edges$to)
nodes <- data.frame(
  name = c("RUNX3", target_nodes),
  type = c("TF", rep("Target", length(target_nodes)))
)

gr <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

# Layout
set.seed(42)

# Plot
p <- ggraph(gr, layout = "fr") +
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

print(p)

# Save outputs
ggsave("RUNX3_network.pdf", p, width = 9, height = 9)
ggsave("RUNX3_network.png", p, width = 9, height = 9, dpi = 300)
