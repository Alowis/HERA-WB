###############################################################################
# COMBINED VALIDATION SUMMARY PLOT
#
# Reads Spearman rho results from all four validated variables:
#   - AET (LISFLOOD vs GLEAM)
#   - SWE (LISFLOOD vs GlobSnow)
#   - Soil Moisture (LISFLOOD vs ESA CCI)
#   - River Discharge (HERA vs observed)
#
# Produces:
#   Fig A: Faceted violin + box (variable × aggregation)
#   Fig B: Climate-stratified heatmap (median rho per variable × climate × window)
#   Fig C: Radar/spider chart showing median rho per climate per variable (Daily)
#   Fig D: Slope chart showing how rho changes from Daily <U+2192> Monthly per climate
#   Tables: summary CSVs
###############################################################################

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(exactextractr)

# --- Paths ---
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
sm_stats_dir <- "Z:/ClimateRun4/nahaUsers/tilloal/SoilMositure_2026_06_03/2.Diego_Analysis/0.Stats_time_windows/"

# ===========================================================================
# 1. LOAD STATS FROM EACH VARIABLE
# ===========================================================================

# --- AET ---
aet_stats <- fread(file.path(base_dir, "output/aet_diego/2.stats/stats_all_windows.csv"))
aet_stats$variable <- "AET"
aet_stats[scenario == "15-Day", scenario := "15-Day"]

# --- SWE ---
swe_stats <- fread(file.path(base_dir, "output/swe_diego/2.stats/stats_all_windows_snowcat.csv"))
swe_stats$variable <- "SWE"
swe_stats[scenario == "15-Day", scenario := "15-Day"]

# --- Soil Moisture (from RDS on network drive) ---
if (dir.exists(sm_stats_dir)) {
  sm_daily <- readRDS(file.path(sm_stats_dir, "stats_daily_snowmasked.rds"))
  sm_7d <- readRDS(file.path(sm_stats_dir, "stats_7d_snowmasked.rds"))
  sm_15d <- readRDS(file.path(sm_stats_dir, "stats_15d_snowmasked.rds"))
  sm_month <- readRDS(file.path(sm_stats_dir, "stats_30d_snowmasked.rds"))
  sm_daily$scenario <- "Daily"
  sm_7d$scenario <- "7-Day"
  sm_15d$scenario <- "15-Day"
  sm_month$scenario <- "Monthly"
  sm_stats <- rbindlist(list(sm_daily, sm_7d, sm_15d, sm_month), fill = TRUE)
  sm_stats$variable <- "Soil Moisture"
  # Remove leading "X" from catch_id (R prepends X to numeric column names)
  sm_stats[, catch_id := sub("^X", "", catch_id)]
} else {
  cat("  WARNING: SM stats directory not found. Skipping Soil Moisture.\n")
  sm_stats <- data.table(
    catch_id = character(), rho = numeric(),
    scenario = character(), variable = character()
  )
}

# --- River Discharge ---
dis_path <- file.path(base_dir, "output/discharge_validation_all_windows.csv")
if (file.exists(dis_path)) {
  dis_stats <- fread(dis_path)
  dis_stats[, rho := spearman_rho]
  dis_stats$variable <- "Discharge"
  dis_stats <- dis_stats[, .(catch_id = catch_id, rho = rho, scenario = scenario, variable = variable)]
} else {
  cat("  WARNING: Discharge stats file not found. Skipping.\n")
  dis_stats <- data.table(
    catch_id = character(), rho = numeric(),
    scenario = character(), variable = character()
  )
}

# ===========================================================================
# 2. COMBINE ALL VARIABLES
# ===========================================================================
all_stats <- rbindlist(list(
  aet_stats[, .(catch_id, rho, scenario, variable)],
  swe_stats[, .(catch_id, rho, scenario, variable)],
  sm_stats[, .(catch_id, rho, scenario, variable)],
  dis_stats[, .(catch_id, rho, scenario, variable)]
), fill = TRUE)

all_stats <- all_stats[!is.na(rho)]
all_stats$scenario <- factor(all_stats$scenario,
  levels = c("Daily", "7-Day", "15-Day", "Monthly")
)
all_stats$variable <- factor(all_stats$variable,
  levels = c("Discharge", "AET", "Soil Moisture", "SWE")
)

# ===========================================================================
# 3. ATTACH CLIMATE ZONES TO CATCHMENTS
# ===========================================================================
cat("  Enriching with climate zones...\n")

file_shp <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
path_clim <- file.path(base_dir, "data", "koppen_geiger_0p1.tif")

shp <- st_read(file_shp, quiet = TRUE)

norm_id <- function(x) {
  x <- as.character(x)
  suppressWarnings(ifelse(grepl("^\\d+$", x), as.character(as.integer(x)), x))
}
shp$join_id <- norm_id(shp$catch_id)

clim_lookup <- c(
  setNames(rep("Tropical", 3), as.character(1:3)),
  setNames(rep("Arid", 4), as.character(4:7)),
  setNames(rep("Temperate", 9), as.character(8:16)),
  setNames(rep("Cold", 12), as.character(17:28)),
  setNames(rep("Polar", 2), as.character(29:30))
)

r_clim <- terra::rast(path_clim)
shp$clim_class <- exactextractr::exact_extract(
  r_clim, sf::st_transform(shp, terra::crs(r_clim)),
  fun = function(values, coverage_fractions) {
    if (all(is.na(values))) {
      return(NA_character_)
    }
    maj <- as.character(names(which.max(table(values[!is.na(values)]))))
    unname(clim_lookup[maj])
  }
)

clim_lut <- st_drop_geometry(shp)[, c("join_id", "clim_class")]
clim_lut <- as.data.table(clim_lut)

# Merge climate into all_stats
all_stats[, join_id := norm_id(catch_id)]
all_stats <- merge(all_stats, clim_lut, by = "join_id", all.x = TRUE)
all_stats$clim_class <- factor(all_stats$clim_class,
  levels = c("Polar", "Cold", "Temperate", "Arid", "Tropical")
)

cat(sprintf(
  "  Combined stats: %d rows, %d with climate info\n",
  nrow(all_stats), sum(!is.na(all_stats$clim_class))
))

# ===========================================================================
# 4. PALETTES AND THEME
# ===========================================================================
clim_pal <- c(
  "Polar" = "#7fbfff", "Cold" = "#2166ac", "Temperate" = "#4dac26",
  "Arid" = "#d6604d", "Tropical" = "#8e0152"
)
var_pal <- c(
  "Discharge" = "#4575b4", "AET" = "#91cf60",
  "Soil Moisture" = "#d73027", "SWE" = "#fc8d59"
)
scen_pal <- c(
  "Daily" = "#bdbdbd", "7-Day" = "#74c4e4",
  "14-Day" = "#2c9e4b", "Monthly" = "#1a3f7a"
)

theme_nature <- function(base_size = 9) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "grey96", colour = "grey60"),
      strip.text = element_text(face = "bold", size = base_size),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      legend.key = element_blank()
    )
}

# ===========================================================================
# 5. FIGURE A: Faceted violin + box (variable × aggregation)
# ===========================================================================
message("Fig A: faceted violin...")
fig_a <- ggplot(all_stats, aes(x = scenario, y = rho, fill = scenario)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.2) +
  geom_violin(alpha = 0.4, trim = FALSE, scale = "width", colour = "grey50", linewidth = 0.3) +
  geom_boxplot(
    width = 0.15, outlier.shape = NA, fill = "white", fatten = 1.5,
    colour = "grey30", linewidth = 0.4
  ) +
  scale_fill_manual(values = scen_pal, guide = "none") +
  scale_y_continuous(limits = c(-0.4, 1.05), breaks = seq(-0.25, 1, 0.25)) +
  facet_wrap(~variable, nrow = 1) +
  labs(
    title = "Cross-variable validation: Spearman \u03c1 by temporal aggregation",
    subtitle = "Dashed line = \u03c1 = 0.5",
    x = "Temporal aggregation", y = "Spearman \u03c1"
  ) +
  theme_nature(10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
    strip.text = element_text(size = 11)
  )




# ===========================================================================
# 6. FIGURE B: Climate-stratified heatmap (median rho per variable x climate x window)
# ===========================================================================
message("Fig B: climate heatmap...")
heat_dt <- all_stats[!is.na(clim_class), .(
  median_rho = round(median(rho, na.rm = TRUE), 2)
), by = .(variable, clim_class, scenario)]

palet <- hcl.colors(11, palette = "RdYlBu", rev = F)

fig_b <- ggplot(heat_dt, aes(x = scenario, y = clim_class, fill = median_rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", median_rho)), size = 2.5, fontface = "bold") +
  scale_fill_gradientn(
    colours = palet,
    limits = c(0, 1), oob = squish, name = "Median \u03c1"
  ) +
  facet_wrap(~variable, nrow = 1) +
  labs(
    title = "Median Spearman \u03c1 by climate zone, variable, and aggregation",
    x = "Temporal aggregation", y = "Climate zone"
  ) +
  theme_bw(15) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 14, color = "white"), # Optional: Makes text white for contrast
    strip.background = element_rect(fill = "#2c3e50"), # Changes the background behind "discharge", "AET", etc.
    panel.grid = element_blank()
  )

# ===========================================================================
# 7. FIGURE C: Slope chart - Daily to Monthly gain per climate & variable
# ===========================================================================
message("Fig C: slope chart...")
slope_dt <- all_stats[!is.na(clim_class) & scenario %in% c("Daily", "Monthly"), .(
  median_rho = median(rho, na.rm = TRUE)
), by = .(variable, clim_class, scenario)]

fig_c <- ggplot(slope_dt, aes(
  x = scenario, y = median_rho,
  group = interaction(variable, clim_class),
  colour = clim_class
)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = clim_pal, name = "Climate") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  facet_wrap(~variable, nrow = 1) +
  labs(
    title = "Aggregation effect: Daily \u2192 Monthly (median \u03c1 per climate zone)",
    x = NULL, y = "Median Spearman \u03c1"
  ) +
  theme_nature(10) +
  theme(legend.position = "bottom", panel.spacing = unit(1, "lines"))

# ===========================================================================
# 8. FIGURE D: Boxplots by climate zone, faceted by variable (Daily only)
# ===========================================================================
message("Fig D: climate boxplots...")
daily_clim <- all_stats[scenario == "Daily" & !is.na(clim_class)]

fig_d <- ggplot(daily_clim, aes(x = clim_class, y = rho, fill = variable)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_boxplot(
    position = position_dodge(width = 0.8), width = 0.7,
    outlier.size = 0.3, outlier.alpha = 0.3, linewidth = 0.35
  ) +
  scale_fill_manual(values = var_pal, name = "Variable") +
  scale_y_continuous(limits = c(-0.3, 1.05), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "Daily Spearman \u03c1 by climate zone and variable",
    x = "Climate zone", y = "Spearman \u03c1"
  ) +
  theme_nature(10) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

# ===========================================================================
# 9. FIGURE E: Ridge density by climate (all variables overlaid)
# ===========================================================================
if (requireNamespace("ggridges", quietly = TRUE)) {
  message("Fig E: ridge density by climate...")
  daily_clim_valid <- all_stats[scenario == "Daily" & !is.na(clim_class)]

  fig_e <- ggplot(
    daily_clim_valid,
    aes(
      x = rho, y = forcats::fct_rev(clim_class),
      fill = variable, colour = variable
    )
  ) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55") +
    ggridges::geom_density_ridges(
      alpha = 0.35, scale = 1.1, bandwidth = 0.04,
      linewidth = 0.4
    ) +
    scale_fill_manual(values = var_pal, name = "Variable") +
    scale_colour_manual(values = var_pal, guide = "none") +
    scale_x_continuous(limits = c(-0.3, 1), breaks = seq(-0.25, 1, 0.25)) +
    labs(
      title = "Daily \u03c1 distribution per climate zone (all variables)",
      x = "Spearman \u03c1", y = NULL
    ) +
    theme_nature(10) +
    theme(legend.position = "bottom")
}

# ===========================================================================
# 10. SUMMARY TABLE
# ===========================================================================
summary_tbl <- all_stats[!is.na(clim_class), .(
  n = .N,
  median_rho = round(median(rho, na.rm = TRUE), 3),
  pct_ge_0.5 = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 1)
), by = .(variable, scenario, clim_class)]

# ===========================================================================
# 11. SAVE FIGURES
# ===========================================================================
path_out <- file.path(base_dir, "output", "combined_validation")
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(path_out, "FigA_combined_rho_violin.png"), fig_a,
  width = 28, height = 12, units = "cm", dpi = 300, bg = "white"
)
ggsave(file.path(path_out, "FigB_climate_heatmap.png"), fig_b,
  width = 30, height = 12, units = "cm", dpi = 300, bg = "white"
)
ggsave(file.path(path_out, "FigC_slope_daily_monthly.png"), fig_c,
  width = 28, height = 12, units = "cm", dpi = 300, bg = "white"
)
ggsave(file.path(path_out, "FigD_climate_boxplots_daily.png"), fig_d,
  width = 22, height = 14, units = "cm", dpi = 300, bg = "white"
)
if (exists("fig_e")) {
  ggsave(file.path(path_out, "FigE_ridge_climate_all_vars.png"), fig_e,
    width = 20, height = 14, units = "cm", dpi = 300, bg = "white"
  )
}

# Combined multi-panel figure for paper
fig_paper <- (fig_a / fig_b / fig_c) +
  plot_annotation(
    tag_levels = "A",
    title = "HERA cross-variable validation summary",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )
ggsave(file.path(path_out, "Fig_combined_paper_v2.png"), fig_paper,
  width = 30, height = 36, units = "cm", dpi = 300, bg = "white"
)

fwrite(summary_tbl, file.path(path_out, "summary_rho_by_variable_climate_window.csv"))

cat(sprintf("\nDone. Figures saved in: %s\n", path_out))
