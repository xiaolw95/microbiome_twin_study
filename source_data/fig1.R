# ============================================================
# Fig 1 — Twin microbiome overview
# Usage: setwd("source_data/"); source("fig1.R")
# Data: fig1.RData
# ============================================================
library(ggplot2); library(dplyr); library(patchwork)
library(reshape2); library(data.table); library(viridis)
library(tidyr); library(mgcv); library(RColorBrewer); library(tibble)
library(ggpubr)
library(scales)

# ---- Load data ----
script_dir <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (is.null(ofile)) {
    file_arg <- commandArgs(trailingOnly = FALSE)
    file_arg <- file_arg[grepl("^--file=", file_arg)]
    if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
  } else dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
}, error = function(e) getwd())
out_dir <- file.path(script_dir, "figures")
if (!file.exists(file.path(script_dir, "fig1.RData")) &&
    file.exists(file.path(getwd(), "source_data", "fig1.RData"))) {
  script_dir <- file.path(getwd(), "source_data")
  out_dir <- file.path(script_dir, "figures")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load(file.path(script_dir, "fig1.RData"), verbose = TRUE)
# Objects: fig1_pcoa, fig1_connections, fig1_age_stats, fig1_ACE_per_stage,
#          fig1_multicontrol, fig1e_twin_distance, fig1c_random, fig1c_twin, fig1d_similarity,
#          fig1_R3_1c_C_hist, fig1_R3_1c_C_summary, fig1_R3_1c_E_violin, fig1_R3_1c_E_summary

# ---- Shared settings ----
custom_colors  <- c("#E0E0E0","#C8C8C8","#B0B0B0","#98B5B4","#86c2c1","#418e88","#1b7b82","#095f67")
contrast_colors <- c("#1b7b82", "#cccccc")
age_colors      <- c("infant"="#E8B4A8","child"="#8BB8C8","adult"="#A8A4C4")
get_twin_colors <- function(n) colorRampPalette(custom_colors)(n)

sp <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white")
}

# ============================================================
# Fig 1B — Twin PCoA with connections
# ============================================================
cat("\n=== Fig 1B: Twin PCoA ===\n")
p_twin_pcoa <- ggplot() +
  geom_segment(data = fig1_connections,
    aes(x = x1, y = y1, xend = x2, yend = y2, color = zygosity, alpha = 1 - distance),
    linewidth = 0.8) +
  geom_point(data = fig1_pcoa,
    aes(x = PC1, y = PC2, fill = age, shape = zygosity, size = shannon),
    color = "black", stroke = 0.0, alpha = 0.8) +
  geom_smooth(data = fig1_pcoa, aes(x = PC1, y = PC2, linetype = population),
    method = "gam", formula = y ~ s(x, k = 5), se = TRUE, alpha = 0.2, color = "gray20") +
  scale_color_manual(name = "Twin Type",
    values = c("MZ" = contrast_colors[2], "DZ" = contrast_colors[1]),
    labels = c("Monozygotic", "Dizygotic")) +
  scale_fill_gradientn(name = "Age\n(years)", colors = get_twin_colors(7),
    breaks = c(0,3,10,18,30,50,70), labels = c("0","3","10","18","30","50","70+")) +
  scale_shape_manual(name = "Twin Type", values = c("MZ" = 21, "DZ" = 22),
    labels = c("Monozygotic", "Dizygotic")) +
  scale_linetype_manual(name = "Population",
    values = c("western" = "solid", "non-western" = "dashed"),
    labels = c("Western", "Non-Western")) +
  scale_size_continuous(name = "Shannon\nDiversity", range = c(3, 6)) +
  scale_alpha_identity() +
  labs(x = "PC1", y = "PC2") +
  theme_classic() +
  theme(legend.position = "right")
sp(p_twin_pcoa, "Fig1B_twin_PCoA", 8, 6)
cat("  1B done.\n")

# ============================================================
# Fig 1C — Distance histogram (3-group: Twin / Same-geo / Random)
# ============================================================
cat("\n=== Fig 1C: Distance histogram (R3_1c) ===\n")
hist_colors <- c("Twin pairs" = "#1b7b82",
                 "Random same country and age" = "#e8a87c",
                 "Random different country and age" = "#cccccc")
plot_dt <- as.data.table(fig1_R3_1c_C_hist)
plot_dt <- plot_dt[pair_type %in% names(hist_colors)]
plot_dt[, pair_type := factor(pair_type, levels = names(hist_colors))]
medians <- plot_dt[, .(median_BC = median(distance)), by = pair_type]
counts  <- plot_dt[, .N, by = pair_type]

p_distance_histogram <- ggplot(plot_dt, aes(x = distance)) +
  # Histograms: bottom layer to top
  geom_histogram(data = subset(plot_dt, pair_type == "Random different country and age"),
    aes(y = after_stat(density), fill = pair_type), alpha = 0.4, bins = 30, color = NA) +
  geom_histogram(data = subset(plot_dt, pair_type == "Random same country and age"),
    aes(y = after_stat(density), fill = pair_type), alpha = 0.55, bins = 30, color = NA) +
  geom_histogram(data = subset(plot_dt, pair_type == "Twin pairs"),
    aes(y = after_stat(density), fill = pair_type), alpha = 0.8, bins = 30, color = NA) +
  # Density lines
  geom_density(data = subset(plot_dt, pair_type == "Random different country and age"),
    aes(color = pair_type), linewidth = 1.0) +
  geom_density(data = subset(plot_dt, pair_type == "Random same country and age"),
    aes(color = pair_type), linewidth = 1.0) +
  geom_density(data = subset(plot_dt, pair_type == "Twin pairs"),
    aes(color = pair_type), linewidth = 1.2) +
  # Median dashed lines
  geom_vline(data = medians, aes(xintercept = median_BC, color = pair_type),
    linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
  scale_fill_manual(values = hist_colors, name = NULL) +
  scale_color_manual(values = hist_colors, name = NULL) +
  scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  labs(x = "Bray-Curtis distance", y = "Density",
       subtitle = paste0(
         "Twin (n=", counts[pair_type == "Twin pairs", N],
         ", med=", sprintf("%.3f", medians[pair_type == "Twin pairs", median_BC]), ") | ",
         "Same country/age (n=", counts[pair_type == "Random same country and age", N],
         ", med=", sprintf("%.3f", medians[pair_type == "Random same country and age", median_BC]), ") | ",
         "Different country/age (n=", counts[pair_type == "Random different country and age", N],
         ", med=", sprintf("%.3f", medians[pair_type == "Random different country and age", median_BC]), ")")) +
  theme_classic(base_size = 11) +
  theme(text = element_text(family = "Arial"),
        plot.subtitle = element_text(size = 8.5, color = "grey40"),
        legend.position = c(0.75, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.key.size = unit(0.8, "lines"),
        legend.text = element_text(size = 8.5))
sp(p_distance_histogram, "Fig1C_distance_histogram", 8, 5)
cat("  1C done.\n")

# ============================================================
# Fig 1D — Co-twin excess similarity over matched random pairs
# ============================================================
cat("\n=== Fig 1D: Excess co-twin similarity ===\n")
age_order <- c("0-10y", "10-30y", "30-60y", "60+y")
excess_dt <- as.data.table(fig1D_excess_similarity)
excess_dt <- excess_dt[age_group %in% age_order]
excess_dt[, age_group := factor(age_group, levels = age_order)]
all_dt <- excess_dt[contrast == "All twins"]
context_dt <- excess_dt[contrast %in% c("MZ twins", "DZ twins")]

p_twin_excess <- ggplot() +
  geom_hline(yintercept = 0, color = "#BDBDBD", linewidth = 0.35) +
  geom_line(data = context_dt,
    aes(x = age_group, y = excess_similarity, group = contrast, color = contrast),
    linewidth = 0.55, alpha = 0.55) +
  geom_point(data = context_dt,
    aes(x = age_group, y = excess_similarity, color = contrast),
    size = 1.8, alpha = 0.55) +
  geom_errorbar(data = all_dt,
    aes(x = age_group, ymin = ci_low, ymax = ci_high),
    width = 0.12, color = "#2B7A78", linewidth = 0.45) +
  geom_line(data = all_dt, aes(x = age_group, y = excess_similarity, group = 1),
    color = "#2B7A78", linewidth = 0.95) +
  geom_point(data = all_dt, aes(x = age_group, y = excess_similarity),
    color = "#2B7A78", size = 2.4) +
  geom_text(data = all_dt,
    aes(x = age_group, y = ci_high + 0.018, label = paste0("n=", n_twin)),
    size = 2.2, color = "grey30") +
  scale_color_manual(values = c("MZ twins" = "#8FA6C1", "DZ twins" = "#A9D5D1"), name = NULL,
                     labels = c("MZ", "DZ")) +
  coord_cartesian(ylim = c(0, max(excess_dt$ci_high, na.rm = TRUE) + 0.05)) +
  labs(x = "Age group", y = "Excess co-twin similarity") +
  theme_classic(base_size = 10) +
  theme(text = element_text(family = "Arial"),
        legend.position = c(0.78, 0.80),
        legend.background = element_rect(fill = alpha("white", 0.75), color = NA),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.grid = element_blank())
sp(p_twin_excess, "Fig1D_excess_similarity", 3.4, 2.4)
cat("  1D done.\n")

# ============================================================
# Fig 1E — Age-stratified MZ/DZ/matched-random violin
# ============================================================
cat("\n=== Fig 1E: Age-stratified violin (R3_1c) ===\n")
cols_e <- c("MZ twins" = "#095f67", "DZ twins" = "#418e88",
            "Random same country and age" = "#e8a87c")

comparisons <- list(
  c("MZ twins", "DZ twins"),
  c("MZ twins", "Random same country and age"),
  c("DZ twins", "Random same country and age")
)

plot_e <- as.data.table(fig1_R3_1c_E_violin)
plot_e <- plot_e[pair_type %in% names(cols_e)]
plot_e[, pair_type := factor(pair_type, levels = names(cols_e))]
plot_e[, age_group := factor(age_group, levels = age_order)]

p1_violin <- ggplot(plot_e, aes(x = pair_type, y = distance, fill = pair_type)) +
  geom_violin(alpha = 0.75, trim = FALSE, scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.3) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",
    label = "p.signif", size = 2.8, bracket.size = 0.3,
    step.increase = 0.08, tip.length = 0.01) +
  scale_fill_manual(values = cols_e, name = NULL) +
  facet_wrap(~ age_group, nrow = 1) +
  coord_cartesian(ylim = c(0, 1.15)) +
  labs(x = NULL, y = "Bray-Curtis dissimilarity") +
  theme_classic(base_size = 10) +
  theme(text = element_text(family = "Arial"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(size = 9, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.position = "bottom")
sp(p1_violin, "Fig1E_distance_violin", 10, 5)
cat("  1E done.\n")

cat("\n")

cat("\nAll Fig 1 panels generated.\n")
