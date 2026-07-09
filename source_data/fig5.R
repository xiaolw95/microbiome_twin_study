# ============================================================
# Fig 5 — Network & functional analysis
# Usage: setwd("source_data/"); source("fig5.R")
# Data: fig5.RData
# ============================================================
library(tidyverse); library(ggplot2); library(patchwork); library(dplyr)
library(ggtext); library(scales)

# ---- Load data ----
load("fig5.RData", verbose = TRUE)
# Objects: fig5_hub_tvr, fig5_hub_3g, fig5_net_tvr, fig5_net_3g,
#          fig5_path_tvr, fig5_path_3g, fig5_shared_tvr, fig5_shared_3g,
#          fig5_net_param, fig5_pathways

# ---- Shared settings ----
colors_main <- c("Random" = "#B0B0B0", "Twin" = "#E64B35", "DZ" = "#4DBBD5", "MZ" = "#00A087")
age_colors  <- c("infant" = "#E8B4A8", "child" = "#8BB8C8", "adult" = "#A8A4C4")

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) + theme(
    text = element_text(family = "Arial", color = "black"),
    plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0, vjust = 1),
    axis.title = element_text(size = base_size, face = "bold"),
    axis.text = element_text(size = base_size - 0.5, color = "black"),
    legend.title = element_text(size = base_size, face = "bold"),
    legend.text = element_text(size = base_size - 0.5),
    legend.key.size = unit(0.4, "cm"), legend.background = element_blank(),
    legend.position = "top", strip.background = element_blank(),
    strip.text = element_text(size = base_size, face = "bold", hjust = 0),
    plot.margin = margin(5, 5, 5, 5), panel.spacing = unit(0.3, "cm"),
    panel.grid = element_blank(), axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5), axis.ticks.length = unit(0.1, "cm"))
}

theme_publication <- function(base_size = 11) {
  theme_classic(base_size = base_size) + theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40", margin = margin(b = 15)),
    axis.title = element_text(size = 11, face = "bold"), axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.5), axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.position = "none", panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(10))
}

sp <- function(p, name, w = 7, h = 5) {
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0(name, ".pdf"), p, width = w, height = h, device = cairo_pdf, bg = "white")
}

# Factor ordering helper
set_age_factors <- function(...) {
  for(d in list(...)) {
    if("Age_Group" %in% names(d)) d$Age_Group <- factor(d$Age_Group, levels = c("Infant","Child","Adult"))
    if("Group" %in% names(d)) d$Group <- factor(d$Group, levels = c("Random","DZ","MZ"))
    if("Pair_Type" %in% names(d)) d$Pair_Type <- factor(d$Pair_Type, levels = c("Random","Twin"))
  }
}
set_age_factors(fig5_hub_tvr, fig5_hub_3g, fig5_net_tvr, fig5_net_3g,
                fig5_path_tvr, fig5_path_3g, fig5_shared_tvr, fig5_shared_3g)

# ============================================================
# Fig 5A — Effect size heatmap (Twin vs Random)
# ============================================================
cat("\n=== Fig 5A: Effect size heatmap ===\n")
calc_effect_size <- function(data_tvr) {
  data_tvr %>% group_by(Age_Group) %>%
    summarise(Mean_Random = Mean[Pair_Type == "Random"], Mean_Twin = Mean[Pair_Type == "Twin"],
              SD_Random = SD[Pair_Type == "Random"], .groups = "drop") %>%
    mutate(Effect_Size = (Mean_Twin - Mean_Random) / SD_Random)
}

effect_all <- bind_rows(
  calc_effect_size(fig5_hub_tvr) %>% mutate(Metric = "Hub Concordance"),
  calc_effect_size(fig5_net_tvr) %>% mutate(Metric = "Network Similarity"),
  calc_effect_size(fig5_path_tvr) %>% mutate(Metric = "Pathway Similarity")) %>%
  mutate(Metric = factor(Metric, levels = c("Network Similarity","Hub Concordance","Pathway Similarity")))

p5A <- ggplot(effect_all, aes(x = Age_Group, y = Metric, fill = Effect_Size)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Effect_Size)), size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0,
    limits = c(-0.5, 3), name = "Effect size\n(Cohen's d)") +
  labs(x = NULL, y = NULL, title = "a") +
  theme_nature() + theme(legend.position = "right", axis.text.x = element_text(angle = 0, hjust = 0.5))
sp(p5A, "Fig5A_effect_heatmap", 5, 3)
cat("  5A done.\n")

# ============================================================
# Fig 5B — Twin advantage trend line
# ============================================================
cat("\n=== Fig 5B: Twin advantage trend ===\n")
calc_improvement <- function(data_tvr) {
  data_tvr %>% group_by(Age_Group) %>%
    summarise(Mean_Random = Mean[Pair_Type == "Random"], Mean_Twin = Mean[Pair_Type == "Twin"],
              SE_Random = SD[Pair_Type == "Random"] / sqrt(N[Pair_Type == "Random"]),
              SE_Twin   = SD[Pair_Type == "Twin"] / sqrt(N[Pair_Type == "Twin"]), .groups = "drop") %>%
    mutate(Improvement = ((Mean_Twin - Mean_Random) / Mean_Random) * 100,
           SE_Improvement = abs(Improvement) * sqrt((SE_Twin/Mean_Twin)^2 + (SE_Random/Mean_Random)^2))
}

improve_all <- bind_rows(
  calc_improvement(fig5_hub_tvr) %>% mutate(Metric = "Hub Concordance"),
  calc_improvement(fig5_net_tvr) %>% mutate(Metric = "Network Similarity"),
  calc_improvement(fig5_path_tvr) %>% mutate(Metric = "Pathway Similarity")) %>%
  mutate(Metric = factor(Metric, levels = c("Network Similarity","Hub Concordance","Pathway Similarity")))

p5B <- ggplot(improve_all, aes(x = Age_Group, y = Improvement, color = Metric, group = Metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, shape = 21, fill = "white", stroke = 1.2) +
  geom_errorbar(aes(ymin = Improvement - SE_Improvement, ymax = Improvement + SE_Improvement), width = 0.15, linewidth = 0.6) +
  scale_color_manual(values = c("Network Similarity" = "#4DBBD5", "Hub Concordance" = "#E64B35", "Pathway Similarity" = "#00A087"), name = NULL) +
  labs(x = "Age group", y = "Twin advantage over random (%)", title = "b") +
  theme_nature() + theme(legend.position = c(0.7, 0.2), legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3))
sp(p5B, "Fig5B_trend", 5, 4)
cat("  5B done.\n")

# ============================================================
# Fig 5C — Three-group gradient bar chart
# ============================================================
cat("\n=== Fig 5C: Three-group gradient ===\n")
prepare_gradient <- function(data_3g, metric_name) {
  data_3g %>% mutate(SE = SD / sqrt(N), Metric = metric_name) %>% select(Age_Group, Group, Mean, SE, Metric)
}

gradient_all <- bind_rows(
  prepare_gradient(fig5_hub_3g, "Hub\nConcordance"),
  prepare_gradient(fig5_net_3g, "Network\nSimilarity"),
  prepare_gradient(fig5_path_3g, "Pathway\nSimilarity")) %>%
  mutate(Metric = factor(Metric, levels = c("Network\nSimilarity","Hub\nConcordance","Pathway\nSimilarity")))

p5C <- ggplot(gradient_all, aes(x = Metric, y = Mean, fill = Group)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(0.8), width = 0.25, linewidth = 0.5) +
  facet_wrap(~Age_Group, nrow = 1) +
  scale_fill_manual(values = colors_main[c("Random","DZ","MZ")], labels = c("Random","DZ twins","MZ twins"), name = NULL) +
  labs(x = NULL, y = "Concordance (mean ± SE)", title = "d") +
  theme_nature() + theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 0.5, size = 6.5), strip.text = element_text(hjust = 0.5))
sp(p5C, "Fig5C_gradient", 7, 4)
cat("  5C done.\n")

# ============================================================
# Fig 5D — Network architecture parameters
# ============================================================
cat("\n=== Fig 5D: Network parameters ===\n")
net_colors <- c("Clustering" = "#1f77b4", "Mean_Degree" = "#ff7f0e", "Positive_Ratio" = "#2ca02c", "Modularity" = "#d62728")
net_data <- fig5_net_param
net_data$Dataset <- factor(net_data$Dataset, levels = c("infant","child","adult"))

main_plot <- net_data %>% filter(variable != "Modularity") %>%
  ggplot(aes(x = Dataset, y = value, group = variable)) +
  geom_smooth(aes(color = variable), method = "loess", se = TRUE, alpha = 0.15, linewidth = 0.8, span = 1.5) +
  geom_line(aes(color = variable), linewidth = 2, alpha = 0.9) +
  geom_point(aes(color = variable, shape = variable), size = 5, stroke = 1.5, alpha = 0.95) +
  scale_color_manual(values = net_colors[1:3], labels = c("Clustering Coefficient","Mean Degree","Positive Ratio")) +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Clustering Coefficient","Mean Degree","Positive Ratio")) +
  labs(y = "Network Index", color = "", shape = "") +
  theme_classic(base_size = 12) + theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 11), axis.title.y = element_text(size = 13, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.8),
    legend.position = "top", legend.justification = "left", legend.text = element_text(size = 11),
    legend.key.size = unit(1.2, "cm"), legend.margin = margin(b = 15),
    panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
    plot.margin = margin(t = 10, r = 15, b = 5, l = 15))

bar_plot <- net_data %>% filter(variable == "Modularity") %>%
  ggplot(aes(x = Dataset, y = value)) +
  geom_col(fill = net_colors["Modularity"], alpha = 0.8, width = 0.6, color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.5, size = 4, color = "black", fontface = "bold") +
  labs(x = "**Developmental Stage**", y = "**Modularity**") +
  scale_y_continuous(limits = c(0, max(net_data$value[net_data$variable == "Modularity"]) * 1.15), expand = c(0, 0)) +
  theme_classic(base_size = 12) + theme(
    axis.text = element_text(size = 11), axis.title.x = element_markdown(size = 13),
    axis.title.y = element_markdown(size = 13), axis.line = element_line(color = "black", linewidth = 0.8),
    panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
    plot.margin = margin(t = 5, r = 15, b = 10, l = 15))

p5D <- main_plot / bar_plot + plot_layout(heights = c(2.2, 1.3)) +
  plot_annotation(title = "**Network Architecture Across Human Development**",
    subtitle = "*Connectivity patterns from infancy to adulthood*",
    theme = theme(plot.title = element_markdown(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_markdown(size = 12, hjust = 0.5, color = "grey40", margin = margin(b = 15))))
sp(p5D, "Fig5D_network", 8, 8)
cat("  5D done.\n")

# ============================================================
# Fig 5F — Top pathways by age
# ============================================================
cat("\n=== Fig 5F: Top pathways ===\n")
fig5_pathways$Age_Group <- factor(fig5_pathways$Age_Group, levels = c("infant","child","adult"))

p5F <- ggplot(fig5_pathways, aes(x = Mean_Abundance, y = reorder(Enhanced_Category, Mean_Abundance))) +
  geom_segment(aes(x = 0, xend = Mean_Abundance, y = Enhanced_Category, yend = Enhanced_Category, color = Age_Group), linewidth = 1.2, alpha = 0.7) +
  geom_point(aes(fill = Age_Group, size = Count, alpha = Mean_Prevalence), shape = 21, color = "black", stroke = 0.8) +
  facet_wrap(~ Age_Group, ncol = 3, scales = "free") +
  scale_color_manual(values = age_colors) + scale_fill_manual(values = age_colors) +
  scale_size_continuous(range = c(3, 8), name = "Count") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Mean Pathway Abundance", y = NULL) +
  theme_publication() + theme(legend.position = "right",
    strip.background = element_rect(fill = "grey95"), panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    axis.text.y = element_text(size = 9)) +
  guides(color = "none", fill = "none", size = guide_legend(ncol = 1))
sp(p5F, "Fig5F_pathways", 10, 6)
cat("  5F done.\n")

cat("\nAll Fig 5 panels generated.\n")
