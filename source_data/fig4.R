# ============================================================
# Fig 4 — Structural variation analysis
# Usage: setwd("source_data/"); source("fig4.R")
# Data: fig4.RData
# ============================================================
library(tidyr); library(tidyverse); library(ggpubr); library(rstatix)
library(RColorBrewer); library(scales); library(ggrepel); library(patchwork)
library(circlize); library(ComplexHeatmap); library(Cairo)

# ---- Load data ----
load("fig4.RData", verbose = TRUE)
# Objects: fig4_dist_raw, fig4_dist_trim, fig4b_dsv, fig4b_vsv,
#          fig4b_diff_dsgv, fig4b_diff_vsgv, fig4b_dsgv_map, fig4b_vsgv_map,
#          fig4c_dsv, fig4c_vsv, fig4c_heatmap, fig4d_sig_func,
#          fig4e_species, fig4e_func_info, fig4e_link_mat, fig4e_heatmap_mat, fig4e_row_anno

# ---- Shared settings ----
sp <- function(p, name, w = 7, h = 5) {
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0(name, ".pdf"), p, width = w, height = h, device = cairo_pdf, bg = "white")
}

zygosity_colors <- c("MZ" = "#364d78", "DZ" = "#1b8d79")
stage_colors    <- c("adult" = "#248e79", "child" = "#354c7a", "infant" = "#c44b36")

theme_violin <- function(base_size = 14) {
  theme_classic(base_size = base_size) + theme(
    axis.title = element_text(size = 16, face = "bold", color = "black", margin = margin(t = 12, r = 12)),
    axis.text  = element_text(size = c(14, 13), color = "black", face = c("bold", "plain")),
    axis.line  = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(), legend.position = "none",
    plot.margin = margin(25))
}

# ---- Reusable violin+boxplot builder ----
make_violin_box <- function(data, x, y, fill_colors, x_labels = NULL, ylab,
                            stat_data = NULL, stat_label = "p", y_expand = c(0.05, 0.18)) {
  p <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = "black", linewidth = 0.8) +
    geom_boxplot(width = 0.12, alpha = 0.9, outlier.shape = NA, color = "black", linewidth = 0.6, fill = "white") +
    scale_fill_manual(values = fill_colors) +
    scale_y_continuous(expand = expansion(mult = y_expand), breaks = scales::pretty_breaks(n = 6)) +
    labs(x = NULL, y = ylab) + theme_violin()
  if(!is.null(x_labels)) p <- p + scale_x_discrete(labels = x_labels)
  if(!is.null(stat_data)) p <- p + stat_pvalue_manual(stat_data, label = stat_label, tip.length = 0.02,
    bracket.nudge.y = 0, step.increase = 0.10, size = 5, fontface = "bold", bracket.size = 0.8)
  p
}

# ============================================================
# Fig 4A — Top species with differential SVs
# ============================================================
cat("\n=== Fig 4A: Species distribution ===\n")
top20 <- fig4_dist_trim %>%
  filter(taxonomic_level == "species") %>% arrange(desc(total_differential_svs)) %>% head(21) %>%
  mutate(taxon_name = factor(taxon_name, levels = rev(taxon_name)),
         dsgv_prop = dsgv_count / total_differential_svs, vsgv_prop = vsgv_count / total_differential_svs)

top20_long <- top20 %>%
  select(taxon_name, total_differential_svs, dsgv_count, vsgv_count) %>%
  pivot_longer(cols = c(dsgv_count, vsgv_count), names_to = "sv_type", values_to = "count") %>%
  mutate(sv_type = factor(sv_type, levels = c("vsgv_count","dsgv_count"), labels = c("vSGV","dSGV")))

p4A <- ggplot(top20_long, aes(x = taxon_name, y = count, fill = sv_type)) +
  geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.3) +
  scale_fill_manual(name = "SV Type", values = c("dSGV" = "#2E86AB", "vSGV" = "#A23B72")) +
  geom_text(data = top20, aes(x = taxon_name, y = total_differential_svs + max(total_differential_svs) * 0.02,
    label = total_differential_svs), inherit.aes = FALSE, hjust = 0, size = 3.5, color = "black", fontface = "bold") +
  coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = scales::pretty_breaks(n = 6)) +
  labs(title = "Top 20 Species with Differential SVs", x = "Species", y = "Number of Differential SVs") +
  theme_classic(base_size = 12) + theme(
    axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = c(11, 10), color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8), axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.background = element_rect(fill = "white"), panel.grid.major.x = element_line(color = "gray90", linewidth = 0.5),
    legend.position = "top", legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), plot.margin = margin(20))
sp(p4A, "Fig4A_species_SV", 8, 7)
cat("  4A done.\n")

# ============================================================
# Fig 4B — dSV & vSV by zygosity
# ============================================================
cat("\n=== Fig 4B: SV by zygosity ===\n")
stat_dsv <- fig4b_dsv %>% wilcox_test(discordance_rate ~ factor_level, p.adjust.method = "bonferroni") %>%
  add_significance("p") %>% add_xy_position(x = "factor_level")
p4B_dsv <- make_violin_box(fig4b_dsv, "factor_level", "discordance_rate", zygosity_colors,
  c("MZ" = "MZ", "DZ" = "DZ"), "Discordance Rate in dSVs", stat_dsv, "p")

stat_vsv <- fig4b_vsv %>% wilcox_test(mean_rel_difference ~ factor_level, p.adjust.method = "bonferroni") %>%
  add_significance("p") %>% add_xy_position(x = "factor_level")
p4B_vsv <- make_violin_box(fig4b_vsv, "factor_level", "mean_rel_difference", zygosity_colors,
  c("MZ" = "MZ", "DZ" = "DZ"), "Mean Relative Difference in vSVs", stat_vsv, "p")

sp(p4B_dsv, "Fig4B_dSV_zygosity", 5, 6)
sp(p4B_vsv, "Fig4B_vSV_zygosity", 5, 6)
cat("  4B done.\n")

# ============================================================
# Fig 4C — dSV & vSV by life stage
# ============================================================
cat("\n=== Fig 4C: SV by life stage ===\n")
fig4c_dsv$factor_level <- factor(fig4c_dsv$factor_level, levels = c("infant","child","adult"))
fig4c_vsv$factor_level <- factor(fig4c_vsv$factor_level, levels = c("infant","child","adult"))

stat_c_dsv <- fig4c_dsv %>% wilcox_test(discordance_rate ~ factor_level, p.adjust.method = "bonferroni") %>%
  add_significance("p") %>% add_xy_position(x = "factor_level")
p4C_dsv <- make_violin_box(fig4c_dsv, "factor_level", "discordance_rate", stage_colors,
  function(x) tools::toTitleCase(x), "Discordance Rate in dSVs", stat_c_dsv, "p.adj.signif")

stat_c_vsv <- fig4c_vsv %>% wilcox_test(mean_rel_difference ~ factor_level, p.adjust.method = "bonferroni") %>%
  add_significance("p") %>% add_xy_position(x = "factor_level")
p4C_vsv <- make_violin_box(fig4c_vsv, "factor_level", "mean_rel_difference", stage_colors,
  function(x) tools::toTitleCase(x), "Mean Relative Difference in vSVs", stat_c_vsv, "p.adj.signif")

sp(p4C_dsv, "Fig4C_dSV_stage", 5, 6)
sp(p4C_vsv, "Fig4C_vSV_stage", 5, 6)
cat("  4C done.\n")

# ============================================================
# Fig 4D — Functional enrichment forest plot
# ============================================================
cat("\n=== Fig 4D: Functional enrichment ===\n")
sig_func <- fig4d_sig_func
# Ensure derived columns exist (may be pre-computed in CSV or need creation)
if(!"function_label" %in% names(sig_func))
  sig_func <- sig_func %>% mutate(function_label = paste0("[", annotation_type, "] ", str_trunc(function_category, 40)))
if(!"enrichment_score_mirror" %in% names(sig_func))
  sig_func$enrichment_score_mirror <- -sig_func$enrichment_score
if(!"neg_log10_p" %in% names(sig_func))
  sig_func$neg_log10_p <- pmin(-log10(sig_func$p_adj), 10)
# sv_count always derived from column 'a' (SV count per function)
sig_func$sv_count <- sig_func$a

custom_colors <- rev(c("#E0E0E0","#C8C8C8","#B0B0B0","#98B5B4","#86c2c1","#418e88","#1b7b82","#095f67"))
p4D <- ggplot(sig_func, aes(y = function_label)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray70", linewidth = 0.3) +
  geom_point(aes(x = enrichment_score_mirror, size = sv_count, fill = neg_log10_p, color = sv_type), shape = 23, stroke = 1) +
  scale_fill_gradientn(colors = custom_colors, name = expression(-log[10](italic(p)[adj]))) +
  scale_color_manual(values = c("dSGV" = "#095f67", "vSGV" = "#7F7F7F"), name = "SV type") +
  scale_size_continuous(range = c(4, 12), name = "SV count") +
  scale_x_continuous(breaks = seq(-6, 6, 2), labels = function(x) abs(x)) +
  labs(title = "Functional enrichment forest plot", subtitle = "dSGV ← | → vSGV",
    x = "Enrichment score (|log2(OR)|)", y = NULL) +
  theme_classic(base_size = 10) + theme(axis.text.y = element_text(size = 9), panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3))
sp(p4D, "Fig4D_forest", 8, 6)
cat("  4D done.\n")

# ============================================================
# Fig 4E — Circos chord diagram (species ↔ function)
# ============================================================
cat("\n=== Fig 4E: Circos diagram ===\n")
species_info <- fig4e_species
function_info <- fig4e_func_info
link_matrix  <- fig4e_link_mat

circos.clear()
circos.par(start.degree = 90,
  gap.degree = c(rep(2, nrow(species_info) - 1), 10, rep(2, nrow(function_info) - 1), 10),
  track.margin = c(0.01, 0.01))

sector_order <- c(species_info$species_short, function_info$function_label)
all_nodes <- rbind(
  data.frame(node_name = species_info$species_short, node_color = species_info$node_color),
  data.frame(node_name = function_info$function_label, node_color = function_info$node_color))

chordDiagram(link_matrix[, c("from","to","value")],
  order = sector_order, grid.col = setNames(all_nodes$node_color, all_nodes$node_name),
  transparency = 0.5, directional = 1, direction.type = "arrows", link.arr.type = "big.arrow",
  link.sort = TRUE, link.decreasing = TRUE,
  annotationTrack = c("grid","name"), annotationTrackHeight = c(0.03, 0.03),
  preAllocateTracks = list(track.height = 0.1))

circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  sn <- get.cell.meta.data("sector.index"); xlim <- get.cell.meta.data("xlim"); ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1], sn, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
    cex = ifelse(sn %in% species_info$species_short, 0.6, 0.5),
    font = ifelse(sn %in% species_info$species_short, 3, 1))
}, bg.border = NA)

cairo_pdf("Fig4E_circos.pdf", width = 10, height = 10, family = "Arial")
# Redraw circos for PDF
circos.clear()
circos.par(start.degree = 90,
  gap.degree = c(rep(2, nrow(species_info) - 1), 10, rep(2, nrow(function_info) - 1), 10),
  track.margin = c(0.01, 0.01))
chordDiagram(link_matrix[, c("from","to","value")],
  order = sector_order, grid.col = setNames(all_nodes$node_color, all_nodes$node_name),
  transparency = 0.5, directional = 1, direction.type = "arrows", link.arr.type = "big.arrow",
  link.sort = TRUE, link.decreasing = TRUE,
  annotationTrack = c("grid","name"), annotationTrackHeight = c(0.03, 0.03),
  preAllocateTracks = list(track.height = 0.1))
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  sn <- get.cell.meta.data("sector.index"); xlim <- get.cell.meta.data("xlim"); ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1], sn, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
    cex = ifelse(sn %in% species_info$species_short, 0.6, 0.5),
    font = ifelse(sn %in% species_info$species_short, 3, 1))
}, bg.border = NA)
dev.off()

png("Fig4E_circos.png", width = 10, height = 10, units = "in", res = 300, bg = "white")
circos.clear()
circos.par(start.degree = 90,
  gap.degree = c(rep(2, nrow(species_info) - 1), 10, rep(2, nrow(function_info) - 1), 10),
  track.margin = c(0.01, 0.01))
chordDiagram(link_matrix[, c("from","to","value")],
  order = sector_order, grid.col = setNames(all_nodes$node_color, all_nodes$node_name),
  transparency = 0.5, directional = 1, direction.type = "arrows", link.arr.type = "big.arrow",
  link.sort = TRUE, link.decreasing = TRUE,
  annotationTrack = c("grid","name"), annotationTrackHeight = c(0.03, 0.03),
  preAllocateTracks = list(track.height = 0.1))
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  sn <- get.cell.meta.data("sector.index"); xlim <- get.cell.meta.data("xlim"); ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1], sn, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
    cex = ifelse(sn %in% species_info$species_short, 0.6, 0.5),
    font = ifelse(sn %in% species_info$species_short, 3, 1))
}, bg.border = NA)
dev.off()
circos.clear()
cat("  4E done.\n")

cat("\nAll Fig 4 panels generated.\n")
