# ============================================================
# Fig 3 — Strain-level effects & heritability dynamics
# Usage: setwd("source_data/"); source("fig3.R")
# Data: fig3.RData
# ============================================================
library(ggplot2); library(dplyr); library(patchwork); library(scales)
library(ggbeeswarm); library(ggpubr); library(rstatix); library(readr)
library(RColorBrewer); library(viridis); library(Cairo); library(tidyr)
library(ggrepel)

# Resolve paths relative to this script when sourced, with getwd() fallback.
script_dir <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
}, error = function(e) getwd())
out_dir <- file.path(script_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ----
load(file.path(script_dir, "fig3.RData"), verbose = TRUE)
# Objects: fig3_strain, fig3_radar_data, fig3_radar_summary, fig3_mz_dz,
#          fig3c_twin_effect, fig3_heatmap, fig3_pie

# ---- Shared settings ----
journal_colors <- list(effects = c(
  "Heritability" = "#D6604D", "Twin Effect" = "#6cb3c5",
  "Geographic Effect" = "#929ab5", "Age Effect" = "#66a797"))
phylum_colors <- c(
  "Bacillota" = "#dd4226", "Bacteroidota" = "#289cb9", "Actinomycetota" = "#2b4578",
  "Proteobacteria" = "#2ECC71", "Verrucomicrobia" = "#F39C12", "Other" = "#95A5A6")

theme_journal <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace% theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.4),
    axis.title = element_text(size = base_size, colour = "black", face = "bold"),
    axis.text = element_text(size = base_size - 1, colour = "black"),
    plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = base_size - 1, colour = "grey30", hjust = 0.5),
    legend.position = "bottom", plot.margin = margin(8, 8, 8, 8))
}
theme_set(theme_journal())

sp <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf, bg = "white")
}

# ============================================================
# Fig 3A — Four-dimensional effect size (strain-level)
# ============================================================
cat("\n=== Fig 3A: Strain effect plot ===\n")
results_df <- fig3_strain
if(!"species_name" %in% colnames(results_df)) results_df$species_name <- results_df$strain_id

strain_order <- c(
  "Hungatella hathewayi","Paraprevotella clara","Prevotella copri","Streptococcus salivarius",
  "Limosilactobacillus fermentum","Romboutsia timonensis","Intestinimonas butyriciproducens",
  "Anaerobutyricum hallii","Ruminococcus bicirculans","Lachnospira eligens",
  "Eubacterium ventriosum","Roseburia intestinalis","Blautia wexlerae",
  "Gordonibacter pamelaeae","Bifidobacterium scardovii","Bifidobacterium breve",
  "Bifidobacterium bifidum","Bifidobacterium pseudocatenulatum","Bifidobacterium adolescentis",
  "Odoribacter splanchnicus","Alistipes putredinis","Alistipes onderdonkii",
  "Parabacteroides merdae","Parabacteroides distasonis","Phocaeicola plebeius",
  "Phocaeicola coprophilus","Phocaeicola dorei","Phocaeicola vulgatus",
  "Bacteroides caccae","Bacteroides caecimuris","Bacteroides fragilis",
  "Bacteroides thetaiotaomicron","Bacteroides ovatus","Bacteroides uniformis")
available <- intersect(strain_order, results_df$species_name)
strain_order <- rev(strain_order[strain_order %in% available])
plot_data <- results_df %>% filter(species_name %in% strain_order) %>% mutate(species_name = factor(species_name, levels = strain_order))

make_effect_plot <- function(data, y_var, sig_var, title, y_label, color, plot_type = "point", ref_lines = NULL, y_limits = NULL, show_y = TRUE) {
  valid <- data %>% filter(!is.na(.data[[y_var]]))
  if(nrow(valid) == 0) return(ggplot() + theme_void() + labs(title = title, subtitle = "No data"))
  p <- ggplot(valid, aes(x = species_name, y = .data[[y_var]]))
  if(!is.null(ref_lines)) for(l in ref_lines) p <- p + geom_hline(yintercept = l, linetype = if(l == 0) "solid" else "dotted", color = if(l == 0) "grey40" else "grey60", linewidth = 0.3, alpha = 0.7)
  if(plot_type == "point") {
    p <- p + geom_point(aes(color = .data[[sig_var]]), size = 2.5, alpha = 0.8) +
      scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = color[[1]]), name = "Significant", labels = c("No","Yes"))
  } else {
    p <- p + geom_col(aes(fill = .data[[sig_var]]), alpha = 0.8, width = 0.7) +
      scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = color[[1]]), name = "Significant", labels = c("No","Yes"))
  }
  p <- p + coord_flip()
  if(!is.null(y_limits)) p <- p + scale_y_continuous(limits = y_limits)
  p + labs(title = title, subtitle = sprintf("Sig: %d/%d", sum(valid[[sig_var]], na.rm = TRUE), sum(!is.na(valid[[y_var]]))),
           x = "", y = y_label) +
    theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8),
      axis.text.y = if(show_y) element_text(size = 7.5, face = "italic") else element_blank(),
      axis.ticks.y = if(show_y) element_line() else element_blank())
}

p1 <- make_effect_plot(plot_data, "twin_cohens_d", "twin_significant", "Twin Effect", "Cohen's d",
  journal_colors$effects["Twin Effect"], "point", c(0, 0.2, 0.5, 0.8), show_y = TRUE)
p2 <- make_effect_plot(plot_data, "age_correlation", "age_significant", "Age Effect", "Spearman ρ",
  journal_colors$effects["Age Effect"], "bar", c(0), show_y = FALSE)
p3 <- make_effect_plot(plot_data, "geographic_effect_size", "geographic_significant", "Geographic Effect", "Cohen's d",
  journal_colors$effects["Geographic Effect"], "bar", c(0), show_y = FALSE)
p4 <- make_effect_plot(plot_data, "heritability", "heritability_significant", "Heritability", "h²",
  journal_colors$effects["Heritability"], "point", c(0, 0.3, 0.6), c(0, 1), show_y = FALSE)

p3A <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4, widths = c(1.3, 1, 1, 1)) +
  plot_annotation(title = "Four-Dimensional Effect Size Analysis",
    subtitle = sprintf("Analysis of %d strains", nrow(plot_data)),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30")))
sp(p3A, "Fig3A_strain_effects", 14, 8)
cat("  3A done.\n")

# ============================================================
# Fig 3B — Radar plot (genus-level effect profiles)
# ============================================================
cat("\n=== Fig 3B: Radar plot ===\n")
radar_coords <- fig3_radar_data
top_genera   <- fig3_radar_summary

extract_phylum <- function(genus_names) {
  phylum_db <- data.frame(genus = c("Bacteroides","Phocaeicola","Parabacteroides","Prevotella","Paraprevotella",
    "Alistipes","Odoribacter","Bifidobacterium","Gordonibacter","Blautia","Roseburia","Eubacterium",
    "Lachnospira","Ruminococcus","Anaerobutyricum","Romboutsia","Intestinimonas","Limosilactobacillus",
    "Streptococcus","Hungatella","Escherichia","Klebsiella","Akkermansia"),
    phylum = c(rep("Bacteroidota",7), rep("Actinomycetota",2), rep("Bacillota",11), rep("Proteobacteria",2), "Verrucomicrobia"),
    stringsAsFactors = FALSE)
  sapply(genus_names, function(g) { p <- phylum_db$phylum[phylum_db$genus == g]; if(length(p) == 0) "Other" else p[1] })
}

if(!"phylum" %in% colnames(top_genera)) top_genera <- top_genera %>% mutate(phylum = extract_phylum(genus))
if(!"phylum" %in% colnames(radar_coords)) radar_coords <- radar_coords %>% left_join(top_genera %>% select(genus, phylum), by = "genus")

genera_with_phylum <- top_genera %>% arrange(phylum, genus) %>%
  mutate(phylum = factor(phylum, levels = names(phylum_colors)[names(phylum_colors) %in% unique(phylum)]))
colors_genus <- c()
for(p in unique(genera_with_phylum$phylum)) {
  g_in_p <- genera_with_phylum$genus[genera_with_phylum$phylum == p]
  n_g <- length(g_in_p)
  base <- phylum_colors[as.character(p)]; if(is.na(base)) base <- phylum_colors["Other"]
  pal <- colorRampPalette(c(adjustcolor(base, alpha.f = 0.5), base, adjustcolor(base, red.f = 0.7, green.f = 0.7, blue.f = 0.7)))(n_g)
  names(pal) <- g_in_p; colors_genus <- c(colors_genus, pal)
}

grid_lines <- do.call(rbind, lapply(c(0.2, 0.4, 0.6, 0.8, 1.0), function(r)
  data.frame(x = r * cos(seq(0, 2*pi, length.out = 100)), y = r * sin(seq(0, 2*pi, length.out = 100)), radius = r)))
effect_types <- c("Twin","Age","Geographic","Heritability")
angles <- seq(0, 2*pi, length.out = 5)[1:4]
axis_lines <- data.frame(x = 0, y = 0, xend = 1.1 * cos(angles), yend = 1.1 * sin(angles))
label_data <- data.frame(x = 1.25 * cos(angles), y = 1.25 * sin(angles), label = effect_types,
  hjust = ifelse(cos(angles) > 0.1, 0, ifelse(cos(angles) < -0.1, 1, 0.5)),
  vjust = ifelse(sin(angles) > 0.1, 0, ifelse(sin(angles) < -0.1, 1, 0.5)))
legend_info <- genera_with_phylum %>% arrange(phylum, genus) %>% mutate(legend_label = paste0(genus, " (", phylum, ", n=", n_strains, ")"))

p_radar <- ggplot() +
  geom_path(data = grid_lines, aes(x = x, y = y, group = radius), color = "grey90", linewidth = 0.3) +
  geom_segment(data = axis_lines, aes(x = x, y = y, xend = xend, yend = yend), color = "grey70", linewidth = 0.4) +
  geom_polygon(data = radar_coords, aes(x = x, y = y, fill = genus, group = genus), alpha = 0.2) +
  geom_path(data = radar_coords, aes(x = x, y = y, color = genus, group = genus), linewidth = 1.2, alpha = 0.9) +
  geom_point(data = radar_coords[!duplicated(paste(radar_coords$genus, radar_coords$effect_type)), ],
    aes(x = x, y = y, color = genus), size = 2.5, alpha = 0.9) +
  geom_text(data = label_data, aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust), size = 4.5, fontface = "bold") +
  annotate("text", x = c(0.2, 0.4, 0.6, 0.8, 1.0), y = 0, label = c("0.2","0.4","0.6","0.8","1.0"), size = 2.5, color = "grey60") +
  scale_color_manual(values = colors_genus, name = "Genus (Phylum, Strains)", labels = setNames(legend_info$legend_label, legend_info$genus)) +
  scale_fill_manual(values = colors_genus, name = "Genus (Phylum, Strains)", labels = setNames(legend_info$legend_label, legend_info$genus)) +
  coord_fixed(ratio = 1, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +
  theme_void() + theme(legend.position = "right", legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9), legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"), plot.margin = margin(20)) +
  labs(title = "Genus-Level Effect Size Profiles (Colored by Phylum)",
    subtitle = sprintf("Mean effect sizes for %d genera", nrow(top_genera)))
sp(p_radar, "Fig3B_radar", 10, 8)
cat("  3B done.\n")

# ============================================================
# Fig 3C — Twin effect & MZ vs DZ per species
# ============================================================
cat("\n=== Fig 3C: Twin effect & MZ vs DZ ===\n")
species_list <- c("Gordonibacter_pamelaeae", "Bifidobacterium_breve")
twin_data  <- fig3c_twin_effect
mz_dz_data <- fig3_mz_dz

for (sp_name in species_list) {
  cat("  Processing:", sp_name, "\n")
  sp_dir <- gsub("_", " ", sp_name)
  # Twin effect
  sp_twin <- twin_data[twin_data$Species == sp_name, ]
  if(nrow(sp_twin) > 0) {
    stat_t <- sp_twin %>% wilcox_test(Distance ~ Effect_Type) %>% add_significance() %>% add_xy_position(x = "Effect_Type")
    p_twin <- ggplot(sp_twin, aes(x = Effect_Type, y = Distance, fill = Effect_Type)) +
      geom_boxplot(alpha = 0.8, width = 0.5, color = "black", linewidth = 0.5) +
      scale_fill_manual(values = c("Within-Twin" = "#c9c4c0", "Between-Twin" = "#00575c")) +
      stat_pvalue_manual(stat_t, label = "p = {p}", size = 3.5, bracket.nudge.y = 0.02, tip.length = 0.01) +
      theme_journal() + labs(title = gsub("_", " ", sp_name), subtitle = "Twin Effect Analysis", x = "", y = "Phylogenetic Distance")
    sp(p_twin, paste0("Fig3C_", sp_dir, "_twin"), 4.5, 5)
  }
  # MZ vs DZ
  sp_mzdz <- mz_dz_data[mz_dz_data$Species == sp_name, ]
  if(nrow(sp_mzdz) > 0) {
    stat_m <- sp_mzdz %>% wilcox_test(Distance ~ Zygosity) %>% add_significance() %>% add_xy_position(x = "Zygosity")
    p_mzdz <- ggplot(sp_mzdz, aes(x = Zygosity, y = Distance, fill = Zygosity)) +
      geom_boxplot(alpha = 0.8, width = 0.5, color = "black", linewidth = 0.5) +
      scale_fill_manual(values = c("MZ" = "#00575c", "DZ" = "#6daaa5")) +
      stat_pvalue_manual(stat_m, label = "p = {p}", size = 3.5, bracket.nudge.y = 0.02, tip.length = 0.01) +
      theme_journal() + labs(title = gsub("_", " ", sp_name), x = "Zygosity", y = "Within-Family Distance")
    sp(p_mzdz, paste0("Fig3C_", sp_dir, "_MZvsDZ"), 4.5, 5)
  }
}
cat("  3C done.\n")

# ============================================================
# Fig 3D - distance-metric robustness of MZ-DZ strain distances
# ============================================================
cat("\n=== Fig 3D: distance-metric robustness ===\n")
long_df <- fig3d_pairwise_distances_long %>%
  filter(!is.na(distance), !is.na(zygosity))
summary_df <- fig3d_metric_summary

metric_labels <- c("p-distance", "Jukes-Cantor", "Kimura-2P")
long_df$metric <- factor(long_df$metric, levels = metric_labels)
long_df$zygosity <- factor(long_df$zygosity, levels = c("MZ", "DZ"))
summary_df$metric <- factor(summary_df$metric, levels = metric_labels)

p3D <- ggplot(long_df, aes(x = metric, y = distance, fill = zygosity)) +
  geom_violin(
    aes(color = zygosity),
    position = position_dodge(width = 0.75),
    width = 0.6,
    alpha = 0.35,
    linewidth = 0.4,
    show.legend = FALSE
  ) +
  geom_boxplot(
    aes(group = interaction(metric, zygosity)),
    position = position_dodge(width = 0.75),
    width = 0.04,
    outlier.shape = NA,
    alpha = 0.9,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  stat_summary(
    aes(group = interaction(metric, zygosity), color = zygosity),
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.75),
    size = 1.8,
    stroke = 0.3,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c(MZ = "#00575c", DZ = "#6daaa5")) +
  scale_color_manual(values = c(MZ = "#00575c", DZ = "#6daaa5")) +
  scale_y_continuous(limits = c(-0.03, 1.10), expand = c(0, 0)) +
  labs(x = NULL, y = "Within-pair strain distance") +
  geom_text(
    data = summary_df,
    aes(x = metric, y = pmax(MZ_median, DZ_median) + 0.10, label = sprintf("p = %.2f", p_wilcox)),
    inherit.aes = FALSE,
    size = 2.6,
    color = "black"
  ) +
  geom_segment(
    data = summary_df,
    aes(x = as.numeric(metric) - 0.32, xend = as.numeric(metric) + 0.32,
        y = pmax(MZ_median, DZ_median) + 0.085, yend = pmax(MZ_median, DZ_median) + 0.085),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.35
  ) +
  theme_classic(base_size = 9, base_family = "Arial") +
  theme(
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 9, margin = margin(r = 6)),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    panel.grid.major.y = element_line(color = "#EDEDED", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    legend.box.spacing = unit(2, "pt"),
    plot.margin = margin(8, 8, 4, 4, unit = "pt")
  )
sp(p3D, "Fig3D_revised_distance_metric_robustness", 5.6, 2.75)
cat("  3D done.\n")

cat("\nAll Fig 3 panels generated.\n")
