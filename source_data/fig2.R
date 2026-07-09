# ============================================================
# Fig 2 — ACE version (cleaned for RData loading)
# Usage:
#   source("source_data/fig2.R")  # from project root
#   or setwd("source_data/"); source("fig2.R")
# Data: fig2.RData
# ============================================================
library(ggplot2); library(ggridges); library(dplyr); library(scales)
library(patchwork); library(tibble); library(ggbeeswarm)
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(Cairo); library(tidyr) })

# Resolve paths relative to this script when sourced, with getwd() fallback.
script_dir <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
}, error = function(e) getwd())

# ---- Load data ----
load(file.path(script_dir, "fig2.RData"), verbose = TRUE)
# Objects: fig2_ACE_heritability, fig2_ACE_plot_data, fig2_ACE_shared_genus,
#          fig2_trajectory, fig2f_data

out_dir <- script_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sp <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white")}

age_colors   <- c("infant" = "#E8B4A8", "child" = "#8BB8C8", "adult" = "#A8A4C4")
env_colors   <- c("Twin similarity" = "#418383", "Geography (|d|)" = "#66809d",
                   "Age effect (|rho|)" = "#4393C3", "Heritability (ACE A)" = "#D6604D")

# ============================================================
cat("\n=== Fig 2A: ACE A density plot ===\n")
h.merged <- fig2_ACE_heritability
h.merged$group <- factor(h.merged$group, levels = c("infant","child","adult"))
h2_stats <- h.merged %>% group_by(group) %>% summarise(median = median(h2, na.rm = TRUE), mean = mean(h2, na.rm = TRUE))
pA <- ggplot(h.merged, aes(x = h2, y = group, fill = group)) +
  ggridges::geom_density_ridges(scale = 0.9, alpha = 0.8, color = "white", linewidth = 0.8, rel_min_height = 0.01) +
  geom_segment(data = h2_stats, aes(x = median, xend = median, y = as.numeric(group), yend = as.numeric(group) + 0.85),
    color = "grey20", linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE) +
  labs(x = "ACE Heritability (A)", y = NULL) +
  theme_classic(base_size = 11) + theme(
    axis.line = element_line(color = "black", linewidth = 0.5), axis.line.y = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5), axis.ticks.y = element_blank(),
    axis.text = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, color = "black", face = "bold", margin = margin(t = 8)),
    legend.position = "none", panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank()) +
  scale_fill_manual(values = age_colors) +
  scale_y_discrete(labels = c("Infant","Child","Adult")) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)), breaks = seq(0, 1, 0.2), limits = c(0, 1))
sp(pA, "Fig2A_ACE_density", 7, 5); cat("  2A done.\n")

# ============================================================
cat("\n=== Fig 2B: C/E bar chart ===\n")
plot_data <- fig2_ACE_plot_data
plot_data$group <- factor(plot_data$group, levels = c("infant","child","adult"))
pB <- ggplot(plot_data, aes(x = group, y = mean_value, fill = env_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = "white", linewidth = 0.3) +
  labs(x = NULL, y = "Average Variation", fill = NULL) +
  theme_classic(base_size = 11) + theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(size = 10, color = "black", face = "bold"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 11, color = "black", face = "bold", margin = margin(r = 8)),
    legend.position = "top", legend.text = element_text(size = 9), legend.key.size = unit(4, "mm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3), panel.grid.major.x = element_blank()) +
  scale_fill_manual(values = c("Shared Environment" = "#66809d", "Unique Environment" = "#418383")) +
  scale_x_discrete(labels = c("Infant","Child","Adult"))
sp(pB, "Fig2B_ACE_bar", 7, 5); cat("  2B done.\n")

# ============================================================
cat("\n=== Fig 2E: ACE composite heatmap ===\n")
shared_group <- fig2_ACE_shared_genus
important_genera <- shared_group %>% group_by(genus) %>% filter(n() >= 2) %>%
  mutate(is_important = (shared_ratio > 0.5) | (spearman_cor > 0.7 & prevalence > 0.5)) %>%
  summarise(n_groups = n(), n_important = sum(is_important), mean_shared = mean(shared_ratio),
    mean_h2 = mean(h2), max_h2 = max(h2), mean_cor = mean(spearman_cor), mean_prev = mean(prevalence)) %>%
  filter(n_important >= 1) %>% arrange(desc(n_important), desc(mean_h2), desc(mean_shared)) %>% slice_head(n = 40)
selected_data <- shared_group %>% filter(genus %in% important_genera$genus)

create_matrix <- function(data, value_col) {
  data %>% select(genus, group, !!sym(value_col)) %>%
    pivot_wider(names_from = group, values_from = !!sym(value_col)) %>%
    column_to_rownames("genus") %>% select(infant, child, adult) %>% as.matrix()}
mat_shared <- create_matrix(selected_data, "shared_ratio")
mat_h2 <- create_matrix(selected_data, "h2")
mat_cor <- create_matrix(selected_data, "spearman_cor")
mat_prev <- create_matrix(selected_data, "prevalence")
mat_abund <- create_matrix(selected_data, "mean_abundance")

high_h2_genera <- selected_data %>% filter(h2 > 0.2) %>% pull(genus) %>% unique()
high_shared_genera <- selected_data %>% filter(shared_ratio > 0.7) %>% group_by(genus) %>% filter(n() >= 2) %>% pull(genus) %>% unique()
core_genera <- selected_data %>% filter(prevalence > 0.5) %>% group_by(genus) %>% filter(n() == 3) %>% pull(genus) %>% unique()
genus_labels <- rownames(mat_shared)
label_colors <- rep("black", length(genus_labels))
label_colors[genus_labels %in% high_h2_genera] <- "#e5722e"
label_colors[genus_labels %in% high_shared_genera] <- "#6984a2"
label_colors[genus_labels %in% core_genera] <- "#0e646c"
genus_category <- ifelse(genus_labels %in% high_h2_genera, "High h2",
  ifelse(genus_labels %in% high_shared_genera, "High Shared",
    ifelse(genus_labels %in% core_genera, "Core", "Other")))

left_anno <- rowAnnotation(Category = genus_category,
  col = list(Category = c("High h2" = "#e5722e", "High Shared" = "#6984a2", "Core" = "#0e646c", "Other" = "grey80")),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"), simple_anno_size = unit(0.4, "cm"))
right_anno <- rowAnnotation(
  "Max A" = anno_points(apply(mat_h2, 1, max, na.rm = TRUE), pch = 16, size = unit(2, "mm"),
    gp = gpar(col = "#E31A1C"), ylim = c(0, 1), axis_param = list(gp = gpar(fontsize = 8))),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_rot = 0)

col_shared <- colorRamp2(c(0, 0.5, 1), c("#ffffff", "#d9d9d9", "#969696"))
col_h2 <- colorRamp2(c(0, 0.2, 0.5), c("#ffffff", "#b5d6e0", "#6cb3c5"))
col_cor <- colorRamp2(c(0, 0.5, 1), c("#ffffff", "#c3c6d5", "#929ab5"))
col_prev <- colorRamp2(c(0, 0.5, 1), c("#FFFFCC", "#A1DAB4", "#2C7FB8"))
col_abund <- colorRamp2(c(0, max(mat_abund)/10, max(mat_abund)/1.8), c("#ffffff", "#cbe3de", "#66a797"))

ht_list <- Heatmap(mat_shared, name = "Shared\nRatio", col = col_shared, cluster_rows = TRUE, cluster_columns = FALSE,
  show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 9, fontface = "italic", col = label_colors),
  column_names_gp = gpar(fontsize = 11, fontface = "bold"), column_names_rot = 0, column_names_centered = TRUE,
  column_title = "Shared Ratio", column_title_gp = gpar(fontsize = 12, fontface = "bold"), width = unit(4, "cm"), border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(!is.na(mat_shared[i, j])) {
      if(mat_shared[i, j] > 0.7) grid.points(x, y, pch = 16, size = unit(3, "mm"), gp = gpar(col = "black"))
      else if(mat_shared[i, j] > 0.5) grid.text("*", x, y, gp = gpar(fontsize = 12, col = "black"))
    }
  }) +
  Heatmap(mat_h2, name = "ACE A", col = col_h2, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, column_title = "ACE Heritability", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_names_rot = 0, column_names_centered = TRUE, width = unit(4, "cm"), border = TRUE, rect_gp = gpar(col = "white", lwd = 1)) +
  Heatmap(mat_cor, name = "Twin\nCorrelation", col = col_cor, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, column_title = "Twin Correlation", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_names_rot = 0, column_names_centered = TRUE, width = unit(4, "cm"), border = TRUE, rect_gp = gpar(col = "white", lwd = 1)) +
  Heatmap(mat_prev, name = "Prevalence", col = col_prev, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, column_title = "Prevalence", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_names_rot = 0, column_names_centered = TRUE, width = unit(4, "cm"), border = TRUE, rect_gp = gpar(col = "white", lwd = 1)) +
  Heatmap(mat_abund, name = "Mean\nAbundance", col = col_abund, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = FALSE, column_title = "Mean Abundance", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    column_names_rot = 0, column_names_centered = TRUE, width = unit(4, "cm"), border = TRUE, rect_gp = gpar(col = "white", lwd = 1))

cairo_pdf(file.path(out_dir, "Fig2E_ACE_heatmap.pdf"), width = 18, height = 12, family = "Arial")
draw(left_anno + ht_list + right_anno, padding = unit(c(2, 2, 2, 10), "mm"),
  heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
grid.text("* Shared ratio > 0.5  ● Shared ratio > 0.7", x = 0.5, y = 0.02, gp = gpar(fontsize = 10, fontface = "italic"))
dev.off()
png(file.path(out_dir, "Fig2E_ACE_heatmap.png"), width = 18, height = 12, units = "in", res = 150, bg = "white")
draw(left_anno + ht_list + right_anno, padding = unit(c(2, 2, 2, 10), "mm"),
  heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
grid.text("* Shared ratio > 0.5  ● Shared ratio > 0.7", x = 0.5, y = 0.02, gp = gpar(fontsize = 10, fontface = "italic"))
dev.off()
cat("  2E done.\n")

# ============================================================
cat("\n=== Fig 2F: Final square ===\n")
d <- fig2f_data
d <- d[!is.na(d$twin_sim) & !is.na(d$cohens_d) & !is.na(d$age_rho), ]
# Rename merged prevalence column
if(!"prevalence" %in% names(d) && "prevalence.y" %in% names(d)) d$prevalence <- d$prevalence.y

make_scat <- function(dt, xvar, yvar, xlab, ylab) {
  ct <- cor.test(dt[[xvar]], dt[[yvar]], method = "spearman")
  ggplot(dt, aes(x = .data[[xvar]], y = .data[[yvar]], fill = prevalence)) +
    geom_point(shape = 21, size = 1.6, alpha = 0.75, color = "white", stroke = 0.15) +
    geom_smooth(method = "lm", se = TRUE, color = "grey50", linewidth = 0.4, alpha = 0.15) +
    scale_fill_viridis_c(option = "D", direction = -1, name = "Prevalence", labels = scales::percent) +
    labs(x = xlab, y = ylab, subtitle = sprintf("rho = %.2f  p = %.2g", ct$estimate, ct$p.value)) +
    theme_classic(base_size = 8) + theme(text = element_text(family = "Arial"),
      panel.grid = element_blank(), plot.subtitle = element_text(size = 7, color = "grey40"),
      legend.position = "right", legend.key.height = unit(3, "mm"),
      legend.title = element_text(size = 7), legend.text = element_text(size = 6))
}
p1 <- make_scat(d, "ACE_A", "twin_sim", "Heritability (ACE A)", "Twin similarity")
p2 <- make_scat(d, "ACE_A", "cohens_d", "Heritability (ACE A)", "Geography (|d|)")
p3 <- make_scat(d, "ACE_A", "age_rho",  "Heritability (ACE A)", "Age effect (|rho|)")
top_row <- p1 + p2 + p3 + plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "bottom")

wt_t <- wilcox.test(d$twin_sim, d$ACE_A, paired = TRUE)
wt_g <- wilcox.test(d$cohens_d, d$ACE_A, paired = TRUE)
wt_a <- wilcox.test(d$age_rho, d$ACE_A, paired = TRUE)
diff_data <- rbind(
  data.frame(genus = d$genus, diff = d$twin_sim - d$ACE_A, force = "Twin similarity", pval = wt_t$p.value),
  data.frame(genus = d$genus, diff = d$cohens_d - d$ACE_A, force = "Geography (|d|)", pval = wt_g$p.value),
  data.frame(genus = d$genus, diff = d$age_rho - d$ACE_A, force = "Age effect (|rho|)", pval = wt_a$p.value))
diff_data$force <- factor(diff_data$force, levels = c("Twin similarity","Geography (|d|)","Age effect (|rho|)"))
diff_data$p_label <- sprintf("p = %.1e\nmed diff = %.2f", diff_data$pval, ave(diff_data$diff, diff_data$force, FUN = median))
diff_data$sig <- ifelse(diff_data$pval < 0.001, "***", ifelse(diff_data$pval < 0.01, "**", ifelse(diff_data$pval < 0.05, "*", "ns")))

p_bot <- ggplot(diff_data, aes(x = force, y = diff, fill = force)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_violin(scale = "width", alpha = 0.4, linewidth = 0.2, color = NA, trim = TRUE) +
  geom_beeswarm(size = 0.35, alpha = 0.25, cex = 0.7) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5, fill = "white", color = "black", stroke = 0.6) +
  geom_text(data = unique(diff_data[, c("force","p_label","sig")]),
    aes(label = sprintf("%s\n%s", p_label, sig)), y = 1.0, size = 2.6, vjust = 0, lineheight = 0.9, color = "grey20", fontface = "bold") +
  scale_fill_manual(values = c("Twin similarity" = "#418383", "Geography (|d|)" = "#66809d", "Age effect (|rho|)" = "#4393C3"), guide = "none") +
  scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0)) +
  labs(x = NULL, y = expression(Delta * " (env - heritability)")) +
  theme_classic(base_size = 9) + theme(text = element_text(family = "Arial"),
    panel.grid = element_blank(), panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    axis.text.x = element_text(size = 7.5, face = "bold"))

pF <- wrap_elements(top_row) / p_bot + plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    title = "F. ACE A is weakly coupled to non-genetic effect axes",
    subtitle = sprintf("Top: Spearman rank correlations across genera (n=%d, colored by prevalence).\nBottom: exploratory native-scale differences; interpret as relative effect-axis summary, not variance decomposition.", nrow(d)),
    theme = theme(plot.title = element_text(size = 10, face = "bold", family = "Arial"),
      plot.subtitle = element_text(size = 7.5, color = "grey40", family = "Arial")))
sp(pF, "Fig2F_final_square", 7, 7)
cat("  2F done.\n")

cat("\nAll Fig 2 panels generated.\n")
