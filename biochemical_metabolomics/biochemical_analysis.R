library(ggplot2)
library(ggbeeswarm)
library(ggsignif)
library(dplyr)

# ============================================================
# User settings - change these to your own paths
# ============================================================

output_dir <- "./results/blood_analysis"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(file.path(output_dir, "individual_plots"))) {
  dir.create(file.path(output_dir, "individual_plots"), recursive = TRUE)
}

# ============================================================
# Biomarker units
# ============================================================

biomarker_units <- list(
  "ALT" = "U/L",
  "AST" = "U/L",
  "AST/ALT" = "ratio",
  "DBIL" = "\u00b5mol/L",
  "IBIL" = "\u00b5mol/L",
  "Tbil" = "\u00b5mol/L",
  "ALB" = "g/L",
  "GLOB" = "g/L",
  "A/G" = "ratio",
  "GGT" = "U/L",
  "Cr" = "\u00b5mol/L",
  "UR" = "mmol/L",
  "UA" = "\u00b5mol/L",
  "CHOL" = "mmol/L",
  "TG" = "mmol/L",
  "LDL-C" = "mmol/L",
  "HDL-C" = "mmol/L",
  "TP" = "g/L",
  "GLU" = "mmol/L"
)

biomarkers <- c("ALT", "AST", "AST/ALT", "DBIL", "IBIL", "Tbil", "ALB", "GLOB", 
                "A/G", "GGT", "Cr", "UR", "UA", "CHOL", "TG", "LDL-C", "HDL-C", "TP", "GLU")

safe_filename <- function(x) gsub("[/\\\\]", "_ratio_", x)

common_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      plot.title = element_text(hjust = 0.5, size = 11),
      plot.margin = unit(c(2, 2, 2, 2), "mm")
    )
}

calc_annotation_pos <- function(values) {
  y_max <- max(values, na.rm = TRUE)
  y_range <- diff(range(values, na.rm = TRUE))
  text_offset <- y_range * 0.08
  line_offset <- y_range * 0.06
  list(
    text = y_max + text_offset,
    line = y_max + line_offset,
    tip_length = min(0.02, y_range * 0.03)
  )
}

stats_table <- data.frame()

pdf(file.path(output_dir, "all_biomarkers.pdf"), width = 2.5, height = 3)

for (biomarker in biomarkers) {
  
  plot_data <- merged_data[!is.na(merged_data[[biomarker]]), ]
  plot_data$group <- factor(plot_data$group, levels = c("SED", "EX"))
  
  pos <- calc_annotation_pos(plot_data[[biomarker]])
  
  formula <- as.formula(paste0("`", biomarker, "` ~ group + sex + Age + BMI"))
  model <- aov(formula, data = plot_data)
  p_value <- summary(model)[[1]]$`Pr(>F)`[1]
  q_value <- p.adjust(p_value, method = "BH")
  
  p <- ggplot(plot_data, aes(x = group, y = .data[[biomarker]])) +
    
    geom_boxplot(
      aes(fill = group),
      width = 0.3,
      outlier.shape = NA,
      linewidth = 0.25,
      alpha = 0.3,
      show.legend = FALSE
    ) +
    
    stat_summary(
      geom = "crossbar",
      width = 0.3,
      fatten = 2,
      fun = median,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    
    geom_beeswarm(
      aes(color = group),
      size = 1.6,
      cex = 3,
      alpha = 0.7,
      dodge.width = 0.3,
      show.legend = FALSE
    ) +
    
    geom_segment(
      aes(x = 1, xend = 2, y = pos$line, yend = pos$line),
      color = "black",
      linewidth = 0.15,
      inherit.aes = FALSE
    ) +
    
    geom_segment(
      aes(x = 1, xend = 1, y = pos$line - pos$tip_length, yend = pos$line),
      color = "black",
      linewidth = 0.15,
      inherit.aes = FALSE
    ) +
    
    geom_segment(
      aes(x = 2, xend = 2, y = pos$line - pos$tip_length, yend = pos$line),
      color = "black",
      linewidth = 0.15,
      inherit.aes = FALSE
    ) +
    
    annotate(
      "text",
      x = 1.5,
      y = pos$text,
      label = if(q_value < 0.0001) "q<0.0001" else sprintf("q=%.4f", q_value),
      size = 3.5,
      vjust = 0
    ) +
    
    scale_color_manual(values = c("SED" = "#3372A6", "EX" = "#367B34")) +
    scale_fill_manual(values = c("SED" = "#3372A6", "EX" = "#367B34")) +
    
    common_theme() +
    
    labs(
      x = "",
      y = ifelse(biomarker %in% names(biomarker_units),
                 paste0("Concentration (", biomarker_units[[biomarker]], ")"),
                 biomarker),
      title = biomarker
    )
  
  ggsave(file.path(output_dir, "individual_plots", paste0(safe_filename(biomarker), ".pdf")),
         plot = p, width = 2.5, height = 3, device = "pdf")
  
  print(p)
  
  stats_table <- rbind(stats_table, data.frame(
    Biomarker = biomarker,
    P_value = p_value,
    Q_value = q_value
  ))
}

dev.off()

write.csv(stats_table, file.path(output_dir, "stat_results.csv"), row.names = FALSE)

