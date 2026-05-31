# ============================================
# Metabolomics and Lipidomics Analysis
# ============================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ============================================
# 1. Volcano Plot
# ============================================

lipid_data <- read.csv("path/to/lipid_data.csv", check.names = FALSE)
metabo_data <- read.csv("path/to/metabo_data.csv", check.names = FALSE)

metabo_data <- metabo_data %>%
  mutate(Type = "Polar Metabolites")

lipid_data <- lipid_data %>%
  mutate(Type = "Complex Lipids")

combined_data <- bind_rows(metabo_data, lipid_data)

combined_data_volcano <- combined_data %>%
  mutate(
    `-log10(pvalue)` = -log10(p.value),
    Significance = case_when(
      abs(`log2(FC)`) >= 0.15 & p.value < 0.05 & `log2(FC)` > 0 ~ "Up-regulated",
      abs(`log2(FC)`) >= 0.15 & p.value < 0.05 & `log2(FC)` < 0 ~ "Down-regulated",
      TRUE ~ "Not Significant"
    )
  )

integrated_volcano_plot <- ggplot(
  data = combined_data_volcano,
  aes(
    x = `log2(FC)`,
    y = `-log10(pvalue)`,
    color = Significance,
    shape = Type,
    size = Type
  )
) +
  geom_point(alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Down-regulated" = "#3372A6",
      "Up-regulated" = "#367B34",
      "Not Significant" = "grey80"
    ),
    name = "Regulation"
  ) +
  scale_shape_manual(
    values = c(
      "Complex Lipids" = 17,
      "Polar Metabolites" = 16
    ),
    name = ""
  ) +
  scale_size_manual(
    values = c(
      "Complex Lipids" = 4.5,
      "Polar Metabolites" = 4.0
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(-1.5, 1.5),
    oob = scales::oob_squish
  ) +
  scale_y_continuous(
    limits = c(0, NA),
    oob = scales::oob_squish
  ) +
  labs(
    title = "",
    x = "log2(Fold Change)",
    y = "-log10(p-value)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(2, 2, 2, 2),
    legend.box.margin = margin(2, 2, 2, 2),
    legend.key = element_rect(fill = "white"),
    legend.background = element_rect(colour = "black", linewidth = 0.2)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3),
      order = 1
    ),
    shape = guide_legend(
      override.aes = list(size = 3),
      order = 2
    )
  ) +
  geom_vline(xintercept = c(-0.15, 0, 0.15), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5)

print(integrated_volcano_plot)

# ============================================
# 2. Donut Chart
# ============================================

significant_molecules <- combined_data_volcano %>%
  filter(Significance != "Not Significant") %>%
  mutate(Molecule_Type = ifelse(Type == "Polar Metabolites", "Metabolite", "Lipid"))

donut_data <- significant_molecules %>%
  group_by(Molecule_Type) %>%
  summarise(Count = n()) %>%
  mutate(
    Percentage = Count / sum(Count) * 100,
    Label = paste0(round(Percentage, 1), "%"),
    ymax = cumsum(Percentage),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

donut_plot <- ggplot(donut_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = Molecule_Type)) +
  geom_rect(alpha = 0.8) +
  scale_fill_manual(
    values = c("Metabolite" = "#D2B48C", "Lipid" = "#D8BFD8"),
    labels = c("Metabolite" = "Metabolites", "Lipid" = "Lipids")
  ) +
  coord_polar(theta = "y") +
  xlim(c(0.5, 4)) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 14, face = "plain")
  ) +
  geom_text(
    aes(x = 3, y = label_pos, label = paste0(Count, "\n", Label)),
    color = "black",
    size = 3.5,
    fontface = "plain"
  ) +
  annotate(
    "text", x = 0.5, y = 0,
    label = paste0(sum(donut_data$Count), "\nmolecules"),
    size = 5,
    fontface = "plain",
    color = "black"
  )

print(donut_plot)
ggsave("donut_chart.pdf", width = 5, height = 4)

# ============================================
# 3. Pathway Enrichment Plot
# ============================================

pathway_df <- data.frame(
  Pathway = c(
    "Glycine, serine and threonine metabolism",
    "Phenylalanine, tyrosine and tryptophan biosynthesis",
    "Phenylalanine metabolism",
    "Cysteine and methionine metabolism",
    "Tyrosine metabolism",
    "Pantothenate and CoA biosynthesis",
    "Valine, leucine and isoleucine biosynthesis",
    "Glyoxylate and dicarboxylate metabolism",
    "Citrate cycle (TCA cycle)",
    "Propanoate metabolism",
    "Pyruvate metabolism",
    "Alanine, aspartate and glutamate metabolism",
    "Glutathione metabolism",
    "Lipoic acid metabolism"
  ),
  Raw_p = c(
    3.86e-10, 4.46e-06, 6.07e-05, 0.000343, 0.000883, 0.00114, 0.00309, 0.00416,
    0.0194, 0.0213, 0.0254, 0.0367, 0.0367, 0.0367
  )
)

pathway_df <- pathway_df %>%
  mutate(
    Log10P = -log10(Raw_p),
    Pathway = factor(Pathway, levels = rev(Pathway))
  ) %>%
  filter(Raw_p < 0.05)

enrichment_plot <- ggplot(pathway_df, aes(x = Log10P, y = Pathway)) +
  geom_col(fill = "#87CEEB", width = 0.7) +
  geom_text(aes(x = 0.1, label = Pathway), hjust = 0, size = 3) +
  labs(
    x = "-Log₁₀(P-value)",
    y = NULL,
    title = "Metabolic Pathway Enrichment Analysis"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

print(enrichment_plot)
ggsave("pathway_enrichment.pdf", width = 8, height = 4)
# ============================================
# 4. Dot Plot with Error Bars for Significant TAG Lipids
# ============================================

library(tidyverse)
library(ggpubr)

# Extract significant TAG lipids
TAG_lipids <- significant_lipids %>%
  filter(str_detect(Molecule_Name, "^TAG")) %>%
  pull(Molecule_Name)

cat("Found", length(TAG_lipids), "significant TAG lipids\n")
print(TAG_lipids)

# Extract these TAG lipids from lipid_long data
selected_data <- lipid_long %>%
  select(Sample_ID, Group, all_of(TAG_lipids))

# Convert to long format for plotting
plot_data <- selected_data %>%
  pivot_longer(
    cols = -c(Sample_ID, Group),
    names_to = "Lipid",
    values_to = "Intensity"
  ) %>%
  group_by(Lipid, Group) %>%
  summarise(
    mean_intensity = mean(Intensity, na.rm = TRUE),
    se_intensity = sd(Intensity, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Plot
dot_error_plot <- ggplot(plot_data, aes(x = Lipid, y = mean_intensity, color = Group)) +
  geom_point(
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = mean_intensity - se_intensity, ymax = mean_intensity + se_intensity),
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  geom_line(
    aes(group = Lipid),
    position = position_dodge(width = 0.5),
    color = "gray50",
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c("SED" = "#3372A6", "EX" = "#367B34"),
    labels = c("SED", "EX")
  ) +
  theme_bw() +
  labs(
    x = "",
    y = "Mean Intensity ± SE",
    color = "Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5)
  )

print(dot_error_plot)

# Save
ggsave("dot_error_plot_TAG.pdf", width = 10, height = 6)



