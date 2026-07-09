# ============================================================
# Fig 1 — Twin microbiome overview
# Usage: setwd("source_data/"); source("fig1.R")
# Data: fig1.RData
# ============================================================
library(ggplot2); library(dplyr); library(patchwork)
library(reshape2); library(data.table); library(viridis)
library(tidyr); library(mgcv); library(RColorBrewer); library(tibble)
library(ggraph); library(igraph); library(tidygraph); library(ggpubr)

# ---- Load data ----
load("fig1.RData", verbose = TRUE)
# Objects: fig1_pcoa, fig1_connections, fig1_age_stats, fig1_ACE_per_stage,
#          fig1_multicontrol, fig1e_twin_distance, fig1c_random, fig1c_twin, fig1d_similarity,
#          fig1_R3_1c_C_hist, fig1_R3_1c_C_summary, fig1_R3_1c_E_violin, fig1_R3_1c_E_summary

# ---- Shared settings ----
custom_colors  <- c("#E0E0E0","#C8C8C8","#B0B0B0","#98B5B4","#86c2c1","#418e88","#1b7b82","#095f67")
contrast_colors <- c("#1b7b82", "#cccccc")
age_colors      <- c("infant"="#E8B4A8","child"="#8BB8C8","adult"="#A8A4C4")
get_twin_colors <- function(n) colorRampPalette(custom_colors)(n)

sp <- function(p, name, w = 7, h = 5) {
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0(name, ".pdf"), p, width = w, height = h, device = cairo_pdf, bg = "white")
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
                 "Random (same geography)" = "#e8a87c",
                 "Random pairs" = "#cccccc")
plot_dt <- fig1_R3_1c_C_hist
medians <- plot_dt[, .(median_BC = median(distance)), by = type]
counts  <- plot_dt[, .N, by = type]

p_distance_histogram <- ggplot(plot_dt, aes(x = distance)) +
  # Histograms: bottom layer to top
  geom_histogram(data = subset(plot_dt, type == "Random pairs"),
    aes(y = after_stat(density), fill = type), alpha = 0.4, bins = 30, color = NA) +
  geom_histogram(data = subset(plot_dt, type == "Random (same geography)"),
    aes(y = after_stat(density), fill = type), alpha = 0.55, bins = 30, color = NA) +
  geom_histogram(data = subset(plot_dt, type == "Twin pairs"),
    aes(y = after_stat(density), fill = type), alpha = 0.8, bins = 30, color = NA) +
  # Density lines
  geom_density(data = subset(plot_dt, type == "Random pairs"),
    aes(color = type), linewidth = 1.0) +
  geom_density(data = subset(plot_dt, type == "Random (same geography)"),
    aes(color = type), linewidth = 1.0) +
  geom_density(data = subset(plot_dt, type == "Twin pairs"),
    aes(color = type), linewidth = 1.2) +
  # Median dashed lines
  geom_vline(data = medians, aes(xintercept = median_BC, color = type),
    linetype = "dashed", linewidth = 0.5, alpha = 0.7) +
  scale_fill_manual(values = hist_colors, name = NULL) +
  scale_color_manual(values = hist_colors, name = NULL) +
  scale_x_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0.02)) +
  labs(x = "Bray-Curtis distance", y = "Density",
       subtitle = paste0(
         "Twin (n=", counts[type == "Twin pairs", N],
         ", med=", sprintf("%.3f", medians[type == "Twin pairs", median_BC]), ") | ",
         "Same-geo (n=", counts[type == "Random (same geography)", N],
         ", med=", sprintf("%.3f", medians[type == "Random (same geography)", median_BC]), ") | ",
         "Random (n=", counts[type == "Random pairs", N],
         ", med=", sprintf("%.3f", medians[type == "Random pairs", median_BC]), ")")) +
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
# Fig 1D — Twin similarity heatmap by age group
# ============================================================
cat("\n=== Fig 1D: Age similarity heatmap ===\n")
sim_long <- fig1d_similarity
sim_long$Age1 <- factor(sim_long$Age1, levels = unique(sim_long$Age1))
sim_long$Age2 <- factor(sim_long$Age2, levels = unique(sim_long$Age2))

p_twin_age_heatmap <- ggplot(sim_long, aes(x = Age1, y = Age2, fill = Twin_Similarity)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(name = "Twin Microbiome\nSimilarity",
    colors = get_twin_colors(9), breaks = scales::pretty_breaks(n = 4)) +
  labs(x = "Age Group (years)", y = "Age Group (years)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9), panel.grid = element_blank())
sp(p_twin_age_heatmap, "Fig1D_age_heatmap", 7, 6)
cat("  1D done.\n")

# ============================================================
# Fig 1D (network) — Age-group network graph
# ============================================================
cat("\n=== Fig 1D: Age network ===\n")
node_data <- fig1_age_stats %>%
  mutate(node_size = n_pairs, within_similarity = mean_similarity,
         mz_ratio = mz_n / (mz_n + dz_n), age_group = age_group,
         age_numeric = case_when(
           age_group == "0-0.5" ~ 0.25, age_group == "0.5-1" ~ 0.75,
           age_group == "1-2" ~ 1.5,   age_group == "2-3" ~ 2.5,
           age_group == "3-6" ~ 4.5,   age_group == "6-12" ~ 9,
           age_group == "12-18" ~ 15,  age_group == "18-30" ~ 24,
           age_group == "30-50" ~ 40,  age_group == "50-60" ~ 55,
           age_group == "60-70" ~ 65,  age_group == "70+" ~ 75),
         life_stage = case_when(
           age_numeric <= 3 ~ "Infancy", age_numeric <= 18 ~ "Childhood", TRUE ~ "Adulthood"))

similarity_threshold <- 0.4
edge_data <- sim_long %>%
  filter(Twin_Similarity >= similarity_threshold, Age1 != Age2) %>%
  mutate(edge_weight = Twin_Similarity,
         similarity_level = cut(Twin_Similarity, breaks = c(0, 0.4, 0.6, 0.8, 1.0),
           labels = c("Moderate","High","Very High","Exceptional"), include.lowest = TRUE),
         age_diff = abs(node_data$age_numeric[match(Age1, node_data$age_group)] -
                        node_data$age_numeric[match(Age2, node_data$age_group)]))

network_graph <- graph_from_data_frame(
  d = edge_data[, c("Age1","Age2","edge_weight","similarity_level","age_diff")],
  vertices = node_data[, c("age_group","node_size","within_similarity","mz_ratio","life_stage","age_numeric")],
  directed = FALSE)
tidy_network <- as_tbl_graph(network_graph)

node_colors <- c("Infancy" = "#E8601C", "Childhood" = "#FD9A44", "Adulthood" = "#1b7b82")
edge_colors <- c("Moderate" = "#B5CEDB", "High" = "#7FADC3", "Very High" = "#4891AD", "Exceptional" = "#2180A1")

p_network_linear <- tidy_network %>%
  ggraph(layout = "linear", circular = FALSE) +
  geom_edge_arc(aes(width = edge_weight, color = similarity_level, alpha = 1/(1 + age_diff/5)),
    strength = 0.3, show.legend = c(color = TRUE, width = TRUE, alpha = FALSE)) +
  geom_node_point(aes(size = node_size, fill = within_similarity, shape = life_stage),
    color = "white", stroke = 1.5, alpha = 0.9) +
  geom_node_text(aes(label = name), size = 3, fontface = "bold", angle = 45, hjust = 1, vjust = 0) +
  scale_edge_color_manual(name = "Similarity Level", values = edge_colors) +
  scale_edge_width_continuous(name = "Similarity Strength", range = c(0.5, 3)) +
  scale_fill_gradientn(name = "Within-Age\nSimilarity", colors = get_twin_colors(7)) +
  scale_shape_manual(name = "Life Stage", values = c("Infancy" = 21, "Childhood" = 22, "Adulthood" = 24)) +
  scale_size_continuous(name = "Twin Pairs", range = c(8, 20)) +
  theme_graph() + theme(legend.position = "right")
sp(p_network_linear, "Fig1D_network", 12, 5)
cat("  1D network done.\n")

# ============================================================
# Fig 1E — Age-stratified 4-group violin (MZ/DZ/SameGeo/Complete)
# ============================================================
cat("\n=== Fig 1E: Age-stratified violin (R3_1c) ===\n")
cols_e <- c("MZ twins" = "#095f67", "DZ twins" = "#418e88",
            "Random (same geography)" = "#e8a87c", "Random (complete)" = "#c9c6c2")

comparisons <- list(
  c("MZ twins", "DZ twins"),
  c("MZ twins", "Random (same geography)"),
  c("MZ twins", "Random (complete)"),
  c("DZ twins", "Random (same geography)"),
  c("DZ twins", "Random (complete)"),
  c("Random (same geography)", "Random (complete)")
)

p1_violin <- ggplot(fig1_R3_1c_E_violin, aes(x = pair_type, y = value, fill = pair_type)) +
  geom_violin(alpha = 0.75, trim = FALSE, scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.3) +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",
    label = "p.signif", size = 2.8, bracket.size = 0.3,
    step.increase = 0.07, tip.length = 0.01) +
  scale_fill_manual(values = cols_e, name = NULL) +
  facet_wrap(~ age_group, nrow = 1) +
  coord_cartesian(ylim = c(0, 1.35)) +
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
