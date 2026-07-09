# ============================================================
# _build_RData.R — Rebuild all fig*.RData from source CSV/TXT
# Usage: setwd("source_data/"); source("_build_RData.R")
# ============================================================
library(data.table); library(readr)

script_dir <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE))
}, error = function(e) getwd())

here <- function(...) file.path(script_dir, ...)

read_first <- function(..., label = NULL) {
  paths <- c(...)
  hit <- paths[file.exists(paths)][1]
  if (is.na(hit)) {
    stop(sprintf("Missing source file%s. Tried:\n  %s",
                 ifelse(is.null(label), "", paste0(" for ", label)),
                 paste(paths, collapse = "\n  ")))
  }
  fread(hit)
}

get_existing_object <- function(rdata_path, object_name) {
  if (!file.exists(rdata_path)) return(NULL)
  e <- new.env()
  load(rdata_path, envir = e)
  if (exists(object_name, envir = e)) get(object_name, envir = e) else NULL
}

cat("Building RData files from source CSV/TXT...\n")

# ---- fig1 ----
cat("  fig1.RData ...")
fig1_pcoa         <- fread(here("fig1b_pcoa_df.csv"))
fig1_connections  <- fread(here("fig1b_twin_connections.csv"))
fig1_age_stats    <- fread(here("fig1d_twin_within_age_stats.csv"))
fig1_ACE_per_stage <- fread(here("fig1_ACE_per_stage.csv"))
fig1_multicontrol <- fread(here("fig1e_twin_distance_5group.csv"))
fig1e_twin_distance <- fread(here("fig1e_twin_distance.csv"))
fig1c_random      <- fread(here("fig1c_random_dist_df.csv"))
fig1c_twin        <- fread(here("fig1c_within_twin_dist_df.csv"))
fig1d_similarity  <- fread(here("fig1_twin_similarity_long.txt"))

# R3_1c preprocessed data (requires twin_metadata.csv for family matching)
meta <- fread(here("..", "data", "twin_metadata.csv"), encoding = "UTF-8")
setnames(meta, 1, "id")
map <- meta[, .(id, family, subject, zygosity)]

dt1e <- copy(fig1e_twin_distance)
dt1e <- merge(dt1e, map[, .(id, fam1 = family, sub1 = subject, zyg1 = zygosity)],
              by.x = "Var1", by.y = "id", all.x = TRUE, sort = FALSE)
dt1e <- merge(dt1e, map[, .(id, fam2 = family, sub2 = subject, zyg2 = zygosity)],
              by.x = "Var2", by.y = "id", all.x = TRUE, sort = FALSE)
dt1e[, is_cotwin := (!is.na(fam1) & fam1 == fam2 & sub1 != sub2)]
dt1e[, is_random := (!is.na(fam1) & !is.na(fam2) & fam1 != fam2)]
dt1e[, same_month := (month1 == month2)]
dt1e[, same_nation := (nation1 == nation2)]
dt1e[, same_geo := (population1 == population2)]

dedup_key <- function(d) {
  d[, pkey := ifelse(Var1 < Var2, paste(Var1, Var2), paste(Var2, Var1))]
  d[!duplicated(pkey)]
}

# Fig1C: 3-group histogram (Twin / Same-geo random / Random)
twin_c <- dt1e[is_cotwin == TRUE & same_month == TRUE]; twin_c <- dedup_key(twin_c)
twin_c[, type := "Twin pairs"]
rand_sg <- dt1e[is_random == TRUE & same_geo == TRUE]; rand_sg <- dedup_key(rand_sg)
rand_sg[, type := "Random (same geography)"]
rand_orig <- copy(fig1c_random); rand_orig[, type := "Random pairs"]
fig1_R3_1c_C_hist <- rbindlist(list(
  twin_c[, .(distance = value, type)],
  rand_sg[, .(distance = value, type)],
  rand_orig[, .(distance, type)]), use.names = TRUE)
fig1_R3_1c_C_hist[, type := factor(type, levels = c("Twin pairs", "Random (same geography)", "Random pairs"))]

# Fig1E: 4-group violin (MZ/DZ/SameGeo/Complete) by age
twin_e <- dt1e[is_cotwin == TRUE & same_month == TRUE & zyg1 %in% c("MZ", "DZ")]
twin_e <- dedup_key(twin_e)
twin_e[, pair_type := ifelse(zyg1 == "MZ", "MZ twins", "DZ twins")]
rand_e <- dt1e[is_random == TRUE]; rand_e <- dedup_key(rand_e)
rand_sg_e <- rand_e[same_geo == TRUE]
age_breaks <- c(0, 120, 360, 720, 1059); age_labels <- c("0-10y", "10-30y", "30-60y", "60+y")
assign_age <- function(m) cut(m, breaks = age_breaks, labels = age_labels, include.lowest = TRUE)
twin_e[, age_group := assign_age(month1)]
rand_sg_e[, age_group := assign_age(month1)]
rand_e[, age_group := assign_age(month1)]
fig1_R3_1c_E_violin <- rbindlist(list(
  twin_e[, .(value, pair_type, age_group)],
  rand_sg_e[, .(value, pair_type = "Random (same geography)", age_group)],
  rand_e[, .(value, pair_type = "Random (complete)", age_group)]), use.names = TRUE)
fig1_R3_1c_E_violin <- fig1_R3_1c_E_violin[!is.na(age_group)]
set.seed(1)
fig1_R3_1c_E_violin <- fig1_R3_1c_E_violin[, .SD[sample(.N, min(.N, 2500))], by = .(pair_type, age_group)]
fig1_R3_1c_E_violin[, pair_type := factor(pair_type, levels = c("MZ twins", "DZ twins", "Random (same geography)", "Random (complete)"))]

# Summary tables (from pre-computed CSVs if available, otherwise compute)
fig1_R3_1c_C_summary <- tryCatch(
  fread(here("..", "data", "result_R3_1c_Fig1C_summary.csv")),
  error = function(e) fig1_R3_1c_C_hist[, .(median_BC = median(distance), mean_BC = mean(distance), n = .N), by = type])
fig1_R3_1c_E_summary <- tryCatch(
  fread(here("..", "data", "result_R3_1c_Fig1E_summary.csv")),
  error = function(e) fig1_R3_1c_E_violin[, .(median_BC = median(value), mean_BC = mean(value), n = .N), by = .(pair_type, age_group)])

save(fig1_pcoa, fig1_connections, fig1_age_stats, fig1_ACE_per_stage, fig1_multicontrol,
     fig1e_twin_distance, fig1c_random, fig1c_twin, fig1d_similarity,
     fig1_R3_1c_C_hist, fig1_R3_1c_C_summary, fig1_R3_1c_E_violin, fig1_R3_1c_E_summary,
     file = here("fig1.RData"), compress = "xz")
cat(" OK\n")

# ---- fig2 ----
cat("  fig2.RData ...")
existing_fig2 <- here("fig2.RData")
fig2_ACE_heritability <- read_first(
  here("..", "revised_figure", "figure2", "panel_A_ACE_A_distribution", "plotdata_fig2A_ACE_A_distribution.csv"),
  here("fig2a_ACE_heritability.csv"),
  here("fig2a_ACE_heritability_updated.csv"),
  here("..", "data", "fig2a_ACE_heritability_updated.csv"),
  label = "Fig2A ACE heritability")
fig2_ACE_plot_data <- read_first(
  here("fig2a_ACE_plot_data.csv"),
  here("fig2a_ACE_plot_data_updated.csv"),
  here("..", "data", "fig2a_ACE_plot_data_updated.csv"),
  label = "Fig2B ACE C/E plot data")
fig2_ACE_shared_genus <- tryCatch(
  read_first(
    here("..", "revised_figure", "figure2", "panel_E_ACE_heatmap", "plotdata_fig2E_full_ACE_shared_genus.csv"),
    here("fig2e_ACE_shared_genus.csv"),
    here("..", "data", "fig2e_ACE_shared_genus.csv"),
    label = "Fig2E ACE shared genus"),
  error = function(e) {
    old <- get_existing_object(existing_fig2, "fig2_ACE_shared_genus")
    if (is.null(old)) stop(e)
    message("    Fig2E source CSV not found; preserving fig2_ACE_shared_genus from existing fig2.RData")
    old
  })
fig2_ACE_composition <- read_first(
  here("..", "revised_figure", "figure2", "panel_B_ACE_composition", "plotdata_fig2B_pooled_ACE_composition.csv"),
  here("Pooled_ACE_by_stage.csv"),
  here("..", "data", "Pooled_ACE_by_stage.csv"),
  label = "Fig2B pooled ACE composition")
fig2_trajectory <- read_first(
  here("Pooled_ACE_by_stage.csv"),
  here("..", "data", "Pooled_ACE_by_stage.csv"),
  label = "Fig2 ACE trajectory")
fig2f_data <- read_first(
  here("..", "revised_figure", "figure2", "panel_F_effect_axis_correlation", "plotdata_fig2F_effect_axis_correlation.csv"),
  here("fig2f_three_panel.csv"),
  here("..", "data", "fig2f_three_panel.csv"),
  label = "Fig2F effect-axis data")
save(fig2_ACE_heritability, fig2_ACE_plot_data, fig2_ACE_shared_genus, fig2_ACE_composition,
     fig2_trajectory, fig2f_data,
     file = here("fig2.RData"), compress = "xz")
cat(" OK\n")

# ---- fig3 ----
cat("  fig3.RData ...")
fig3_strain        <- fread(here("fig3a_corrected_strain_effects_results.csv"))
fig3_radar_data    <- fread(here("fig3b_radar_plot_data.csv"))
fig3_radar_summary <- fread(here("fig3b_radar_plot_summary.csv"))
fig3_mz_dz         <- fread(here("fig3c_mz_dz_comparison_data.csv"))
fig3c_twin_effect  <- fread(here("fig3c_twin_effect_data.csv"))
fig3_heatmap       <- fread(here("fig3d_plot_data_heatmap.csv"))
fig3_pie           <- fread(here("fig3d_plot_data_pie.csv"))
fig3d_pairwise_distances_long <- read_first(
  here("..", "revised_figure", "figure3", "panel_D_distance_metric_robustness", "plotdata_fig3D_pairwise_distances_long.csv"),
  here("..", "data", "reproducibility", "R1_2c_strain_distance_metrics", "plotdata_pairwise_distances.csv"),
  label = "Fig3D distance-metric pairwise distances")
if ("p_distance" %in% names(fig3d_pairwise_distances_long)) {
  metric_labels <- c(p_distance = "p-distance", jukes_cantor = "Jukes-Cantor", kimura_2p = "Kimura-2P")
  fig3d_pairwise_distances_long <- rbindlist(lapply(names(metric_labels), function(m) {
    tmp <- fig3d_pairwise_distances_long[, .(sgb, zygosity, distance = get(m))]
    tmp[, metric := metric_labels[[m]]]
    tmp
  }), use.names = TRUE)
}
fig3d_metric_summary <- read_first(
  here("..", "revised_figure", "figure3", "panel_D_distance_metric_robustness", "summary_fig3D_MZ_vs_DZ_by_metric.csv"),
  here("..", "data", "reproducibility", "R1_2c_strain_distance_metrics", "plotdata_MZ_vs_DZ_by_metric.csv"),
  label = "Fig3D MZ vs DZ metric summary")
save(fig3_strain, fig3_radar_data, fig3_radar_summary, fig3_mz_dz, fig3c_twin_effect,
     fig3_heatmap, fig3_pie, fig3d_pairwise_distances_long, fig3d_metric_summary,
     file = here("fig3.RData"), compress = "xz")
cat(" OK\n")

# ---- fig4 ----
cat("  fig4.RData ...")
fig4_dist_raw      <- fread(here("fig4a_species_distribution_raw.csv"))
fig4_dist_trim     <- fread(here("fig4a_species_distribution_trim.csv"))
fig4b_dsv          <- fread(here("fig4b_dSV.csv"))
fig4b_vsv          <- fread(here("fig4b_vSV.csv"))
fig4b_diff_dsgv    <- fread(here("fig4b_differential_dSGV.csv"))
fig4b_diff_vsgv    <- fread(here("fig4b_differential_vSGV.csv"))
fig4b_dsgv_map     <- fread(here("fig4b_dSGV_mapping.csv"))
fig4b_vsgv_map     <- fread(here("fig4b_vSGV_mapping.csv"))
fig4c_dsv          <- fread(here("fig4c_dSV.csv"))
fig4c_vsv          <- fread(here("fig4c_vSV.csv"))
fig4c_heatmap      <- fread(here("fig4c_heatmap.csv"))
fig4d_sig_func     <- fread(here("fig4d_significant_func.csv"))
fig4e_species      <- fread(here("fig4e_species_info.csv"))
fig4e_func_info    <- fread(here("fig4e_function_info.csv"))
fig4e_link_mat     <- fread(here("fig4e_link_matrix.csv"))
fig4e_heatmap_mat  <- fread(here("fig4e_heatmap_matrix.csv"))
fig4e_row_anno     <- fread(here("fig4e_row_annotation.csv"))
save(fig4_dist_raw, fig4_dist_trim,
     fig4b_dsv, fig4b_vsv, fig4b_diff_dsgv, fig4b_diff_vsgv, fig4b_dsgv_map, fig4b_vsgv_map,
     fig4c_dsv, fig4c_vsv, fig4c_heatmap,
     fig4d_sig_func,
     fig4e_species, fig4e_func_info, fig4e_link_mat, fig4e_heatmap_mat, fig4e_row_anno,
     file = here("fig4.RData"), compress = "xz")
cat(" OK\n")

# ---- fig5 ----
cat("  fig5.RData ...")
fig5_hub_tvr    <- fread(here("fig5a_Summary_Hub_TwinVsRandom.txt"))
fig5_hub_3g     <- fread(here("fig5a_Summary_Hub_ThreeGroups.txt"))
fig5_net_tvr    <- fread(here("fig5a_Summary_Network_TwinVsRandom.txt"))
fig5_net_3g     <- fread(here("fig5a_Summary_Network_ThreeGroups.txt"))
fig5_path_tvr   <- fread(here("fig5a_Summary_Pathway_TwinVsRandom.txt"))
fig5_path_3g    <- fread(here("fig5a_Summary_Pathway_ThreeGroups.txt"))
fig5_shared_tvr <- fread(here("fig5a_Summary_Shared_TwinVsRandom.txt"))
fig5_shared_3g  <- fread(here("fig5a_Summary_Shared_ThreeGroups.txt"))
fig5_net_param  <- fread(here("fig5d_net_parameter.csv"))
fig5_pathways   <- fread(here("fig5f_top_pathways_plot.csv"))
save(fig5_hub_tvr, fig5_hub_3g, fig5_net_tvr, fig5_net_3g,
     fig5_path_tvr, fig5_path_3g, fig5_shared_tvr, fig5_shared_3g,
     fig5_net_param, fig5_pathways,
     file = here("fig5.RData"), compress = "xz")
cat(" OK\n")

cat("\nAll RData files rebuilt.\n")
