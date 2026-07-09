# ============================================================
# Fig 2 - ACE model source-data plotting script
# Usage:
#   source("source_data/fig2.R")  # from project root
#   or setwd("source_data/"); source("fig2.R")
# Data: fig2.RData
# ============================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(ggsignif)
  library(Cairo)
  library(grid)
  library(tibble)
})

script_dir <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
}, error = function(e) getwd())

load(file.path(script_dir, "fig2.RData"), verbose = TRUE)
# Objects used here:
#   fig2_ACE_heritability, fig2_ACE_composition, fig2_ACE_shared_genus, fig2f_data

out_dir <- file.path(script_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w, h, dpi = 600) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
}

# ============================================================
cat("\n=== Fig 2A: ACE A distribution ===\n")
df <- fig2_ACE_heritability %>%
  filter(!is.na(group), !is.na(h2)) %>%
  mutate(
    h2 = as.numeric(h2),
    group = factor(group, levels = c("infant", "child", "adult"))
  )

stats <- df %>%
  group_by(group) %>%
  summarise(n = n(), median = median(h2), mean = mean(h2), .groups = "drop")
cat("  Data:", nrow(df), "rows |", paste(stats$group, "n=", stats$n, "median=", round(stats$median, 3), collapse = " | "), "\n")

pA <- ggplot(df, aes(group, h2, fill = group, colour = group)) +
  geom_violin(alpha = 0.45, linewidth = 0.6, trim = TRUE) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5,
               fill = "black", colour = "black", stroke = 0.6) +
  geom_signif(
    comparisons = list(c("infant", "child"), c("infant", "adult"), c("child", "adult")),
    test = "wilcox.test",
    p.adjust.method = "bonferroni",
    map_signif_level = TRUE,
    textsize = 2.8,
    linewidth = 0.35,
    tip_length = 0.01,
    y_position = c(1.04, 1.10, 1.16)
  ) +
  scale_fill_manual(values = c("#cccccc", "#6fafae", "#005352")) +
  scale_colour_manual(values = c("#cccccc", "#6fafae", "#005352")) +
  scale_x_discrete(labels = c("Infant", "Child", "Adult")) +
  labs(y = expression("ACE additive genetic effect"), x = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(margin = margin(t = 4)),
    axis.title.y = element_text(size = 9.5, margin = margin(r = 8)),
    axis.ticks = element_line(linewidth = 0.4),
    axis.line = element_line(linewidth = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(8, 8, 12, 4, unit = "pt"),
    legend.position = "none"
  )
save_plot(pA, "Fig2A_nature_boxviolin", 3.4, 3.0)
cat("  2A done.\n")

# ============================================================
cat("\n=== Fig 2B: pooled ACE composition ===\n")
df_comp <- fig2_ACE_composition
df_long <- df_comp %>%
  mutate(stage_raw = factor(tolower(stage), levels = c("infant", "child", "adult"))) %>%
  arrange(stage_raw) %>%
  pivot_longer(cols = c(A, C, E), names_to = "component", values_to = "value") %>%
  mutate(
    component = factor(component, levels = c("E", "C", "A")),
    x_label = case_when(
      stage_raw == "infant" ~ paste0("Infant\n", stage_n_MZ, " MZ / ", stage_n_DZ, " DZ"),
      stage_raw == "child"  ~ paste0("Child\u2020\n", stage_n_MZ, " MZ / ", stage_n_DZ, " DZ"),
      stage_raw == "adult"  ~ paste0("Adult\n", stage_n_MZ, " MZ / ", stage_n_DZ, " DZ")
    )
  )
df_long$x_label <- factor(df_long$x_label, levels = unique(df_long$x_label))
df_label <- df_long %>%
  filter(value >= 0.075) %>%
  mutate(label_fontface = ifelse(component == "A", "bold", "plain"))

pB <- ggplot(df_long, aes(x = x_label, y = value, fill = component)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.8) +
  geom_text(
    data = df_label,
    aes(label = sprintf("%.2f", value), fontface = label_fontface),
    position = position_stack(vjust = 0.5),
    size = 2.6,
    color = "white",
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(A = "#D6604D", C = "#66809d", E = "#418383"),
    labels = c(expression(italic("A")), expression(italic("C")), expression(italic("E")))
  ) +
  scale_y_continuous(limits = c(0, 1.04), expand = c(0, 0)) +
  labs(y = "Pooled ACE fraction", x = NULL) +
  guides(fill = guide_legend(title = NULL, keyheight = unit(0.35, "cm"), keywidth = unit(0.35, "cm"), order = 3)) +
  annotate(
    "text", x = 0.03, y = -0.12,
    label = "\u2020 Childhood stratum has fewer twin pairs; interpret this estimate cautiously.",
    size = 2.4, hjust = 0, color = "grey40"
  ) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    axis.text.x = element_text(size = 8.5, color = "black", margin = margin(t = 4)),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 9.5, margin = margin(r = 8)),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(8, 8, 20, 4, unit = "pt"),
    legend.position = "top",
    legend.text = element_text(size = 8.5, margin = margin(r = 6)),
    legend.box.spacing = unit(2, "pt")
  )
save_plot(pB, "Fig2B_nature_ACE_composition", 3.6, 3.0)
cat("  2B done.\n")

# ============================================================
cat("\n=== Fig 2E: ACE heatmap ===\n")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
  cat("  Skipped: ComplexHeatmap and/or circlize is not installed.\n")
} else {
Heatmap <- ComplexHeatmap::Heatmap
rowAnnotation <- ComplexHeatmap::rowAnnotation
anno_points <- ComplexHeatmap::anno_points
draw <- ComplexHeatmap::draw
colorRamp2 <- circlize::colorRamp2
shared_group <- fig2_ACE_shared_genus
important_genera <- shared_group %>%
  group_by(genus) %>%
  filter(n() >= 2) %>%
  mutate(is_important = (shared_ratio > 0.5) | (spearman_cor > 0.7 & prevalence > 0.5)) %>%
  summarise(
    n_groups = n(),
    n_important = sum(is_important),
    mean_shared = mean(shared_ratio),
    mean_h2 = mean(h2),
    max_h2 = max(h2),
    mean_cor = mean(spearman_cor),
    mean_prev = mean(prevalence),
    .groups = "drop"
  ) %>%
  filter(n_important >= 1) %>%
  arrange(desc(n_important), desc(mean_h2), desc(mean_shared)) %>%
  slice_head(n = 40)

selected_data <- shared_group %>% filter(genus %in% important_genera$genus)
category_table <- selected_data %>%
  group_by(genus) %>%
  summarise(
    high_h2 = any(h2 > 0.2, na.rm = TRUE),
    high_shared = sum(shared_ratio > 0.7, na.rm = TRUE) >= 2,
    max_A = max(h2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(category = case_when(
    high_h2 ~ "High h2",
    high_shared ~ "High Shared",
    TRUE ~ "Other"
  )) %>%
  arrange(factor(category, levels = c("High h2", "High Shared", "Other")), desc(max_A), genus)
genus_order <- category_table$genus

create_matrix <- function(data, value_col, row_order) {
  mat <- data %>%
    select(genus, group, !!sym(value_col)) %>%
    pivot_wider(names_from = group, values_from = !!sym(value_col)) %>%
    column_to_rownames("genus") %>%
    select(infant, child, adult) %>%
    as.matrix()
  mat[row_order, , drop = FALSE]
}

mat_shared <- create_matrix(selected_data, "shared_ratio", genus_order)
mat_h2 <- create_matrix(selected_data, "h2", genus_order)
mat_cor <- create_matrix(selected_data, "spearman_cor", genus_order)
mat_prev <- create_matrix(selected_data, "prevalence", genus_order)
mat_abund <- create_matrix(selected_data, "mean_abundance", genus_order)
stopifnot(identical(rownames(mat_shared), rownames(mat_h2)),
          identical(rownames(mat_h2), rownames(mat_cor)),
          identical(rownames(mat_cor), rownames(mat_prev)),
          identical(rownames(mat_prev), rownames(mat_abund)))

genus_labels <- rownames(mat_shared)
genus_category <- setNames(category_table$category, category_table$genus)[genus_labels]
label_colors <- unname(c("High h2" = "#e5722e", "High Shared" = "#6984a2", "Other" = "black")[genus_category])

left_anno <- rowAnnotation(
  Category = genus_category,
  col = list(Category = c("High h2" = "#e5722e", "High Shared" = "#6984a2", "Other" = "grey80")),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  simple_anno_size = unit(0.4, "cm"),
  width = unit(0.8, "cm")
)
right_anno <- rowAnnotation(
  "Max A" = anno_points(apply(mat_h2, 1, max, na.rm = TRUE), pch = 16, size = unit(2, "mm"),
                        gp = gpar(col = "#E31A1C"), ylim = c(0, 1),
                        axis_param = list(gp = gpar(fontsize = 8))),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  annotation_name_rot = 0
)

ht_list <- Heatmap(
  mat_shared,
  name = "Shared\nRatio",
  col = colorRamp2(c(0, 0.5, 1), c("#ffffff", "#d9d9d9", "#969696")),
  cluster_rows = TRUE,
  row_split = factor(genus_category, levels = c("High h2", "High Shared", "Other")),
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9, fontface = "italic", col = label_colors),
  left_annotation = left_anno,
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_title = "Shared Ratio",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  width = unit(4, "cm"),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(mat_shared[i, j])) {
      if (mat_shared[i, j] > 0.7) {
        grid.points(x, y, pch = 16, size = unit(3, "mm"), gp = gpar(col = "black"))
      } else if (mat_shared[i, j] > 0.5) {
        grid.text("*", x, y, gp = gpar(fontsize = 12, col = "black"))
      }
    }
  }
) +
  Heatmap(mat_h2, name = "ACE A", col = colorRamp2(c(0, 0.2, 0.5), c("#ffffff", "#b5d6e0", "#6cb3c5")),
          cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
          column_title = "ACE Heritability", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
          column_names_gp = gpar(fontsize = 11, fontface = "bold"), column_names_rot = 0,
          column_names_centered = TRUE, width = unit(4, "cm"), border = TRUE, rect_gp = gpar(col = "white", lwd = 1)) +
  Heatmap(mat_cor, name = "Twin\nCorrelation", col = colorRamp2(c(0, 0.5, 1), c("#ffffff", "#c3c6d5", "#929ab5")),
          cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
          column_title = "Twin Correlation", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
          column_names_gp = gpar(fontsize = 11, fontface = "bold"), column_names_rot = 0,
          column_names_centered = TRUE, width = unit(4, "cm"), border = TRUE, rect_gp = gpar(col = "white", lwd = 1)) +
  Heatmap(mat_prev, name = "Prevalence", col = colorRamp2(c(0, 0.5, 1), c("#FFFFCC", "#A1DAB4", "#2C7FB8")),
          cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
          column_title = "Prevalence", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
          column_names_gp = gpar(fontsize = 11, fontface = "bold"), column_names_rot = 0,
          column_names_centered = TRUE, width = unit(4, "cm"), border = TRUE, rect_gp = gpar(col = "white", lwd = 1)) +
  Heatmap(mat_abund, name = "Mean\nAbundance",
          col = colorRamp2(c(0, max(mat_abund, na.rm = TRUE) / 10, max(mat_abund, na.rm = TRUE) / 1.8),
                           c("#ffffff", "#cbe3de", "#66a797")),
          cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
          column_title = "Mean Abundance", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
          column_names_gp = gpar(fontsize = 11, fontface = "bold"), column_names_rot = 0,
          column_names_centered = TRUE, width = unit(4, "cm"), border = TRUE, rect_gp = gpar(col = "white", lwd = 1))

CairoPDF(file.path(out_dir, "Fig2E_revised_ACE_heatmap.pdf"), width = 18, height = 12, family = "Arial")
draw(ht_list + right_anno, padding = unit(c(2, 2, 2, 10), "mm"),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
grid.text("* Shared ratio > 0.5  \u25cf Shared ratio > 0.7", x = 0.5, y = 0.02,
          gp = gpar(fontsize = 10, fontface = "italic"))
dev.off()
png(file.path(out_dir, "Fig2E_revised_ACE_heatmap.png"), width = 18, height = 12, units = "in", res = 150, bg = "white")
draw(ht_list + right_anno, padding = unit(c(2, 2, 2, 10), "mm"),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
grid.text("* Shared ratio > 0.5  \u25cf Shared ratio > 0.7", x = 0.5, y = 0.02,
          gp = gpar(fontsize = 10, fontface = "italic"))
dev.off()
cat("  2E done.\n")
}

# ============================================================
cat("\n=== Fig 2F: ACE A and non-genetic effect axes ===\n")
d <- fig2f_data %>%
  filter(!is.na(ACE_A), !is.na(twin_sim), !is.na(cohens_d), !is.na(age_rho), !is.na(prevalence))
axes <- list(
  list(col = "twin_sim", label = "Twin similarity", color = "#418383"),
  list(col = "cohens_d", label = "Geography (|d|)", color = "#66809d"),
  list(col = "age_rho", label = "Age effect (|\u03c1|)", color = "#4393C3")
)

panel_list <- lapply(axes, function(ax) {
  sp <- cor.test(d$ACE_A, d[[ax$col]], method = "spearman")
  p_str <- if (sp$p.value < 0.001) sprintf("p = %.1e", sp$p.value) else sprintf("p = %.2f", sp$p.value)
  plot_df <- data.frame(ACE_A = d$ACE_A, y_val = d[[ax$col]], prevalence = d$prevalence)
  ggplot(plot_df, aes(x = ACE_A, y = y_val, color = prevalence)) +
    geom_point(size = 1.8, alpha = 0.78, stroke = 0.2) +
    scale_color_viridis_c(option = "D", direction = -1, name = "Prevalence") +
    geom_smooth(method = "lm", se = FALSE, color = ax$color, linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("\u03c1 = %.2f\n%s", sp$estimate, p_str),
             hjust = -0.05, vjust = 1.5, size = 2.8, color = "grey30") +
    labs(x = expression("ACE " * italic("A")), y = ax$label) +
    theme_classic(base_size = 9, base_family = "Arial") +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      panel.grid = element_blank(),
      legend.position = "none"
    )
})
panel_list[[3]] <- panel_list[[3]] +
  theme(legend.position = "right") +
  guides(color = guide_colorbar(
    barwidth = unit(0.25, "cm"),
    barheight = unit(2.5, "cm"),
    label.theme = element_text(size = 7),
    title.theme = element_text(size = 8),
    title.position = "top"
  ))
pF <- ggarrange(plotlist = panel_list, ncol = 3, nrow = 1, widths = c(1, 1, 1.25), align = "h")
save_plot(pF, "Fig2F_nature_effect_axis_correlation", 7.2, 2.4)
cat("  2F done.\n")

cat("\nAll revised Fig 2 panels generated.\n")
