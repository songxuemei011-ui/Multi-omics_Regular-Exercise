#!/usr/bin/env Rscript

# ============================================
# Metabolomics and Lipidomics Analysis
# ANCOVA + Volcano Plot + Donut Chart + Pathway Enrichment + TAG Dot Plot
# ============================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(tidyverse)
library(ggpubr)
library(broom)
library(readr)

# ============================================
# Part 1: ANCOVA Differential Analysis
# ============================================

run_ancova <- function(data, data_name) {
  
  non_measure_cols <- c("Sample_name", "exercise_group_new", "sex", "Age", "BMI")
  biomarkers <- setdiff(names(data), non_measure_cols)
  
  cat("\n", data_name, "- testing", length(biomarkers), "features\n")
  
  results <- data.frame(
    biomarker = character(),
    p_value = numeric(),
    log2FC = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (biomarker in biomarkers) {
    formula <- as.formula(paste("`", biomarker, "` ~ exercise_group_new + sex + Age + BMI"))
    model <- aov(formula, data = data)
    model_summary <- broom::tidy(model)
    p_value_row <- model_summary %>% filter(term == "exercise_group_new")
    
    if (nrow(p_value_row) > 0) {
      p_value <- p_value_row$p.value
      mean_sed <- mean(data[[biomarker]][data$exercise_group_new == "SED"], na.rm = TRUE)
      mean_ex <- mean(data[[biomarker]][data$exercise_group_new == "EX"], na.rm = TRUE)
      log2fc <- log2(mean_ex / mean_sed)
      
      results <- rbind(results, data.frame(
        biomarker = biomarker,
        p_value = p_value,
        log2FC = log2fc,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  results$q_value <- p.adjust(results$p_value, method = "fdr")
  results <- results %>%
    mutate(
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      significant = ifelse(p_value < 0.05 & abs(log2FC) >= 0.15, TRUE, FALSE)
    ) %>%
    arrange(p_value)
  
  return(results)
}

# Load data (modify paths)
metabolite_data <- read_csv("path/to/preprocessed_metabolites.csv")
lipid_data <- read_csv("path/to/preprocessed_lipids.csv")

cat("Metabolites data:", dim(metabolite_data), "\n")
cat("Lipids data:", dim(lipid_data), "\n")

# Run ANCOVA
metabolite_results <- run_ancova(metabolite_data, "Metabolites")
lipid_results <- run_ancova(lipid_data, "Lipids")

# Save results
write.csv(metabolite_results, "differential_metabolites_ancova.csv", row.names = FALSE)
write.csv(lipid_results, "differential_lipids_ancova.csv", row.names = FALSE)

# Save significant only
metabolite_sig <- metabolite_results %>% filter(significant == TRUE)
lipid_sig <- lipid_results %>% filter(significant == TRUE)

write.csv(metabolite_sig, "differential_metabolites_significant.csv", row.names = FALSE)
write.csv(lipid_sig, "differential_lipids_significant.csv", row.names = FALSE)

cat("Metabolites - total:", nrow(metabolite_results), "| significant:", nrow(metabolite_sig), "\n")
cat("Lipids - total:", nrow(lipid_results), "| significant:", nrow(lipid_sig), "\n")

# ============================================
# Part 2: Volcano Plot
# ============================================

# Add type column
metabolite_results$Type <- "Polar Metabolites"
lipid_results$Type <- "Complex Lipids"

combined_results <- bind_rows(metabolite_results, lipid_results)

combined_volcano <- combined_results %>%
  mutate(
    neg_log10_q = -log10(q_value),
    regulation = case_when(
      significant & log2FC > 0 ~ "Up-regulated",
      significant & log2FC < 0 ~ "Down-regulated",
      TRUE ~ "Not Significant"
    )
  )

volcano_plot <- ggplot(combined_volcano, aes(x = log2FC, y = neg_log10_q, color = regulation, shape = Type)) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(
    values = c(
      "Up-regulated" = "#367B34",
      "Down-regulated" = "#3372A6",
      "Not Significant" = "grey70"
    )
  ) +
  scale_shape_manual(values = c("Polar Metabolites" = 16, "Complex Lipids" = 17)) +
  geom_vline(xintercept = c(-0.15, 0.15), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    x = "log2 Fold Change (EX vs SED)",
    y = "-log10(q-value)",
    title = "Differential Metabolites & Lipids"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("volcano_plot_combined.pdf", volcano_plot, width = 8, height = 6)

# ============================================
# Part 3: Donut Chart
# ============================================

significant_combined <- combined_results %>% filter(significant == TRUE)

donut_data <- significant_combined %>%
  group_by(Type) %>%
  summarise(Count = n()) %>%
  mutate(
    Percentage = Count / sum(Count) * 100,
    Label = paste0(round(Percentage, 1), "%"),
    ymax = cumsum(Percentage),
    ymin = c(0, head(ymax, n = -1)),
    label_pos = (ymax + ymin) / 2
  )

donut_plot <- ggplot(donut_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2, fill = Type)) +
  geom_rect(alpha = 0.8) +
  scale_fill_manual(
    values = c("Polar Metabolites" = "#D2B48C", "Complex Lipids" = "#D8BFD8"),
    labels = c("Polar Metabolites" = "Metabolites", "Complex Lipids" = "Lipids")
  ) +
  coord_polar(theta = "y") +
  xlim(c(0.5, 4)) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank()
  ) +
  geom_text(aes(x = 3, y = label_pos, label = paste0(Count, "\n", Label)), size = 3.5) +
  annotate("text", x = 0.5, y = 0,
           label = paste0(sum(donut_data$Count), "\nmolecules"),
           size = 5)

ggsave("donut_chart.pdf", donut_plot, width = 5, height = 4)

# ============================================
# Part 4: Pathway Enrichment Plot
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
  filter(Raw_p < 0.05) %>%
  mutate(row_num = row_number())

# Define pathway categories
pathway_categories <- list(
  "Amino acid" = c(
    "Glycine, serine and threonine metabolism",
    "Phenylalanine, tyrosine and tryptophan biosynthesis",
    "Phenylalanine metabolism",
    "Cysteine and methionine metabolism",
    "Tyrosine metabolism",
    "Valine, leucine and isoleucine biosynthesis",
    "Alanine, aspartate and glutamate metabolism"
  ),
  "Energy" = c(
    "Glyoxylate and dicarboxylate metabolism",
    "Citrate cycle (TCA cycle)",
    "Propanoate metabolism",
    "Pyruvate metabolism"
  ),
  "Cofactors and vitamins" = c(
    "Pantothenate and CoA biosynthesis",
    "Lipoic acid metabolism"
  ),
  "Peptide" = c(
    "Glutathione metabolism"
  )
)

# Assign categories
pathway_df$Category <- NA
for (category in names(pathway_categories)) {
  pathway_df$Category[pathway_df$Pathway %in% pathway_categories[[category]]] <- category
}

pathway_df <- pathway_df %>% filter(!is.na(Category))

# Order categories
category_order <- pathway_df %>%
  group_by(Category) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count)) %>%
  pull(Category)

pathway_df$Category <- factor(pathway_df$Category, levels = category_order)

# Gradient color function
get_category_color <- function(category, value, max_value) {
  intensity <- value / max_value
  switch(as.character(category),
    "Amino acid" = rgb(0.8 - intensity * 0.3, 0.6 - intensity * 0.2, 0.7 - intensity * 0.2),
    "Energy" = rgb(0.5 - intensity * 0.2, 0.6 - intensity * 0.2, 0.8 - intensity * 0.2),
    "Cofactors and vitamins" = rgb(0.6 - intensity * 0.2, 0.5 - intensity * 0.2, 0.7 - intensity * 0.2),
    "Peptide" = rgb(0.5 - intensity * 0.2, 0.7 - intensity * 0.2, 0.5 - intensity * 0.2),
    rgb(0.5, 0.5, 0.5)
  )
}

# Apply colors
max_log10p <- max(pathway_df$Log10P)
pathway_df$bar_color <- mapply(get_category_color, pathway_df$Category, pathway_df$Log10P, max_log10p)

# Color bar data
colorbar_data <- pathway_df %>%
  group_by(Category) %>%
  summarise(
    start = min(row_num),
    end = max(row_num),
    mid = (start + end) / 2,
    avg_log10p = mean(Log10P)
  ) %>%
  mutate(
    bar_color = mapply(get_category_color, Category, avg_log10p, max_log10p),
    start_adj = start + 0.3,
    end_adj = end - 0.3
  )

# Main plot
p_main <- ggplot(pathway_df, aes(x = Log10P, y = reorder(Pathway, Log10P))) +
  geom_col(aes(fill = I(bar_color)), width = 0.8) +
  geom_text(aes(x = 0.1, label = Pathway), hjust = 0, color = "black", size = 3.2) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, hjust = 0.5),
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    plot.margin = margin(l = 60),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = "-Log₁₀(P-value)", y = NULL, title = "Metabolic Pathway Enrichment Analysis") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

# Color bar
p_colorbar <- ggplot(colorbar_data) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = start_adj, ymax = end_adj, fill = I(bar_color)), alpha = 0.9) +
  geom_text(aes(x = 0.5, y = mid, label = Category), color = "white", size = 3.5,
            angle = 90, hjust = 0.5, vjust = 0.5) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(r = 5)) +
  scale_y_continuous(limits = c(0.5, nrow(pathway_df) + 0.5)) +
  coord_cartesian(expand = FALSE, xlim = c(0, 1))

# Combine
final_plot <- p_colorbar + p_main + plot_layout(widths = c(1, 15))

print(final_plot)
ggsave("pathway_enrichment.pdf", final_plot, width = 10, height = 6)


# ============================================
# 5. Dot Plot with Error Bars for Significant TAG Lipids
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
