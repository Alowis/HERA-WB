###############################################################################
# LISFLOOD vs ESA CCI SOIL MOISTURE — EXTENDED CROSS-COMPARISON ANALYSIS
#
# DATA INPUTS (confirmed from diagnostic output):
#   sm_timeseries  keys : "daily" / "7d" / "15d" / "monthly", each $lisf + $esa
#   monthly date col    : "year_month"
#   stats_* columns     : catch_id, rho, p_val, status  (+  area_km2, nesting cols)
#   stats_month extra   : catch_id.y  (duplicate — dropped on load)
#   scenario labels     : "Daily" / "7-Day" / "15-Day" / "Monthly"
#
# OUTPUTS  (written to path_out):
#   Fig0  Spatial rho maps (4-panel, your existing function)
#   Fig1  rho violin + box across aggregations
#   Fig2  Scatter — 4 aggregations, coloured by climate class
#   Fig3  Scatter matrix — climate rows x aggregation columns
#   Fig4  Stratified diagnostics — climate / area / elevation boxplots
#   Fig5  Ridge density of rho per climate class
#   Fig6  Aggregation-gain heatmap (<U+0394>rho vs daily)
#   Table_rho_by_climate/area/elevation.csv
###############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(exactextractr)
library(ggplot2)
library(patchwork)
library(zoo)
library(scales)
library(ggridges)
library(forcats)
library(grid)


# =============================================================================
# 0.  PATHS
# =============================================================================
path_stats <- "Z:\\ClimateRun4\\nahaUsers\\tilloal\\SoilMositure_2026_06_03\\2.Diego_Analysis\\0.Stats_time_windows\\"
path_out <- file.path(base_dir, "output", "figures")
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)
file_shp <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
path_clim <- file.path(base_dir, "data", "koppen_geiger_0p1.tif")
path_dem <- file.path(base_dir, "data", "dem.nc")


# =============================================================================
# 1.  LOAD DATA
# =============================================================================
message("Loading data...")

sm_timeseries <- readRDS(file.path(path_stats, "SM_time_series_masked_all_scenarios.rds"))

drop_geom <- function(x) if (inherits(x, "sf")) sf::st_drop_geometry(x) else x

# Drop duplicate catch_id.y that stats_month carries from its internal join
clean_stats <- function(path) {
  drop_geom(readRDS(path)) %>% dplyr::select(-any_of("catch_id.y"))
}

stats_daily <- clean_stats(file.path(path_stats, "stats_daily_snowmasked.rds"))
stats_7d <- clean_stats(file.path(path_stats, "stats_7d_snowmasked.rds"))
stats_15d <- clean_stats(file.path(path_stats, "stats_15d_snowmasked.rds"))
stats_month <- clean_stats(file.path(path_stats, "stats_30d_snowmasked.rds"))

shp <- st_read(file_shp, quiet = TRUE)

strip_id <- function(df) {
  df$join_id <- as.character(gsub("^X", "", df$catch_id))
  df
}

map_d <- left_join(shp, strip_id(stats_daily), by = c("catch_id" = "join_id"))
map_7 <- left_join(shp, strip_id(stats_7d), by = c("catch_id" = "join_id"))
map_15 <- left_join(shp, strip_id(stats_15d), by = c("catch_id" = "join_id"))
map_m <- left_join(shp, strip_id(stats_month), by = c("catch_id" = "join_id"))

eu_mask <- st_read("D:\\gomezdi\\My Data\\0. Projects\\3.WildFireVulnerab\\Datasets\\TRACE\\Luc_delivery\\NUTS3_wgs84LatLon.shp",
  quiet = TRUE
)
EU30 <- c(
  "HR", "DE", "BE", "AT", "ES", "BG", "FR", "CZ", "HU", "EL", "FI", "CY", "DK", "EE",
  "NL", "LV", "MT", "LT", "IT", "PT", "PL", "SE", "IE", "LU", "RO", "UK", "SK", "SI", "NO", "CH"
)

message("Data loaded.")


# =============================================================================
# 2.  CLIMATE + ELEVATION + AREA ENRICHMENT
# =============================================================================
message("=== Step 2: Climate / elevation / area enrichment ===")

clim_lookup <- c(
  setNames(rep("Tropical", 3), as.character(1:3)),
  setNames(rep("Arid", 4), as.character(4:7)),
  setNames(rep("Temperate", 9), as.character(8:16)),
  setNames(rep("Cold", 12), as.character(17:28)),
  setNames(rep("Polar", 2), as.character(29:30))
)

r_clim <- rast(path_clim)
shp$clim_class <- exact_extract(
  r_clim, shp,
  fun = function(values, coverage_fractions) {
    if (all(is.na(values))) {
      return(NA_character_)
    }
    maj_val <- as.character(names(which.max(table(values[!is.na(values)]))))
    unname(clim_lookup[maj_val])
  }
)

# GMTED raster is WGS84 — reproject shapefile to match before extraction
r_dem <- rast(path_dem)
shp_wgs84 <- st_transform(shp, crs = 4326)
shp$elev_m <- exact_extract(r_dem, shp_wgs84, fun = "mean")
message(
  "  elev_m range: ", round(min(shp$elev_m, na.rm = TRUE), 1),
  " to ", round(max(shp$elev_m, na.rm = TRUE), 1), " m  |  NAs: ", sum(is.na(shp$elev_m))
)

shp <- shp %>%
  mutate(
    clim_class = factor(clim_class, levels = c("Polar", "Cold", "Temperate", "Arid", "Tropical")),
    elev_class = cut(
      elev_m,
      breaks = c(-Inf, 200, 500, 1000, 2000, Inf),
      labels = c("< 200 m", "200-500 m", "500-1000 m", "1000-2000 m", "> 2000 m"),
      right  = TRUE
    ),
    area_class = cut(
      area_km2,
      breaks         = quantile(area_km2, probs = c(0, .25, .50, .75, 1), na.rm = TRUE),
      labels         = c("Q1 (smallest)", "Q2", "Q3", "Q4 (largest)"),
      include.lowest = TRUE
    )
  )

message("Enrichment done.")


# =============================================================================
# 3.  ALL-SCENARIO META TABLE  (rho + catchment attributes)
# =============================================================================
message("=== Step 3: Assembling stratified rho table ===")

meta_lut <- st_drop_geometry(shp) %>%
  dplyr::select(catch_id, clim_class, area_km2, area_class, elev_m, elev_class)

attach_meta <- function(stats_df, scenario_label) {
  stats_df %>%
    mutate(catch_id = as.character(gsub("^X", "", catch_id))) %>%
    left_join(meta_lut, by = "catch_id") %>%
    filter(!is.na(rho)) %>%
    mutate(scenario = scenario_label)
}

all_meta <- bind_rows(
  attach_meta(stats_daily, "Daily"),
  attach_meta(stats_7d, "7-Day"),
  attach_meta(stats_15d, "15-Day"),
  attach_meta(stats_month, "Monthly")
) %>%
  mutate(scenario = factor(scenario, levels = c("Daily", "7-Day", "15-Day", "Monthly")))

message(
  "all_meta rows: ", nrow(all_meta),
  "  | scenarios: ", paste(levels(all_meta$scenario), collapse = ", ")
)


# =============================================================================
# 4.  LONG SCATTER TABLE  (sm_mod x sm_obs per catchment x time-step)
# =============================================================================
message("=== Step 4: Building scatter table ===")

clim_lut_dt <- as.data.table(meta_lut)[, .(catch_id, clim_class)]
clim_lut_dt[, catch_id_x := paste0("X", catch_id)]

build_scatter_dt <- function(lisf_wide, esa_wide, scenario_label, date_col = "date") {
  lisf_wide <- as.data.table(lisf_wide)
  esa_wide <- as.data.table(esa_wide)
  l <- melt(lisf_wide, id.vars = date_col, variable.name = "catch_id_x", value.name = "sm_mod")
  e <- melt(esa_wide, id.vars = date_col, variable.name = "catch_id_x", value.name = "sm_obs")
  dt <- merge(l, e, by = c(date_col, "catch_id_x"))
  dt <- dt[!is.na(sm_mod) & !is.na(sm_obs)]
  dt <- merge(dt, clim_lut_dt[, .(catch_id_x, clim_class)], by = "catch_id_x", all.x = TRUE)
  dt[, scenario := scenario_label]
  dt
}

# Monthly: rename year_month -> date
# lisf_m_dt <- as.data.table(sm_timeseries$monthly$obs)
# esa_m_dt  <- as.data.table(sm_timeseries$monthly$mod)
# setnames(lisf_m_dt, "month", "date")
# setnames(esa_m_dt,  "month", "date")

plot_dt <- rbind(
  build_scatter_dt(lisf_wide = sm_timeseries$daily$mod, esa_wide = sm_timeseries$daily$obs, "a) Daily"),
  build_scatter_dt(sm_timeseries$`7d`$mod, sm_timeseries$`7d`$obs, "b) 7-Day"),
  build_scatter_dt(sm_timeseries$`15d`$mod, sm_timeseries$`15d`$obs, "c) 15-Day"),
  build_scatter_dt(lisf_wide = sm_timeseries$monthly$mod, esa_wide = sm_timeseries$monthly$obs, "d) Monthly"),
  fill = TRUE
)

plot_dt[, scenario := factor(scenario, levels = c("a) Daily", "b) 7-Day", "c) 15-Day", "d) Monthly"))]
plot_dt[, clim_class := factor(clim_class, levels = c("Polar", "Cold", "Temperate", "Arid", "Tropical"))]

ax_lim <- range(
  quantile(plot_dt$sm_mod, probs = c(0.02, 0.99), na.rm = TRUE),
  quantile(plot_dt$sm_obs, probs = c(0.02, 0.99), na.rm = TRUE)
)

# Per-panel annotation: rho + n
ann <- plot_dt[!is.na(clim_class), .(
  rho = round(cor(sm_mod, sm_obs, method = "spearman", use = "complete.obs"), 2),
  n   = .N
), by = .(scenario, clim_class)]
ann[, label := paste0("\u03c1 = ", rho, "\nn = ", format(n, big.mark = ","))]

# Stratified subsample for rendering (max 3000 pts per climate x scenario)
set.seed(42)

plot_dt_sub <- plot_dt[!is.na(clim_class),
  {
    tmp <- copy(.SD)
    tmp[, xbin := cut(sm_mod, breaks = 100)]
    tmp[, .SD[sample(.N, min(.N, 100L))], by = xbin]
  },
  by = .(scenario, clim_class)
]
message("Scatter table ready. Total rows: ", nrow(plot_dt_sub))


# =============================================================================
# 5.  SHARED THEME + PALETTES
# =============================================================================

theme_nature <- function(base_size = 9) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
      strip.background = element_rect(fill = "grey96", colour = "grey60", linewidth = 0.4),
      strip.text = element_text(face = "bold", size = base_size),
      axis.ticks = element_line(colour = "grey40", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      plot.title = element_text(
        face = "bold", size = base_size + 2, hjust = 0,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = base_size - 0.5, colour = "grey35",
        hjust = 0, margin = margin(b = 6)
      ),
      plot.caption = element_text(
        size = base_size - 1.5, colour = "grey50",
        hjust = 0, face = "italic"
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      plot.margin = margin(6, 8, 4, 6)
    )
}

clim_pal <- c(
  "Polar"     = "#7fbfff",
  "Cold"      = "#2166ac",
  "Temperate" = "#4dac26",
  "Arid"      = "#d6604d",
  "Tropical"  = "#8e0152"
)

scen_pal <- c(
  "Daily"   = "#bdbdbd",
  "7-Day"   = "#74c4e4",
  "15-Day"  = "#2c9e4b",
  "Monthly" = "#1a3f7a"
)


# =============================================================================
# FIG 0  |  SPATIAL RHO MAPS  (your existing function, unmodified)
# =============================================================================
message("=== Fig 0: Spatial rho maps ===")

message("Fig 0: spatial rho maps...")
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
basemap <- sf::st_transform(world, crs = 3035)

build_map <- function(spatial_df, plot_title) {
  sp <- sf::st_transform(spatial_df, 3035)
  nco <- sf::st_coordinates(sf::st_centroid(sp))
  ggplot() +
    geom_sf(data = basemap, fill = "beige", color = NA) +
    geom_sf(data = sp, fill = "grey70", color = "transparent", linewidth = 0.05) +
    geom_sf(
      data = sp, aes(fill = rho),
      color = "transparent", linewidth = 0
    ) +
    # geom_sf(
    #     data = dplyr::filter(sp, status == "Not Significant"),
    #     fill = "grey75", color = "transparent", linewidth = 0
    # ) +
    geom_sf(data = basemap, fill = "transparent", linewidth = 0.2) +
    scale_fill_gradientn(
      colors = c("#ffffff", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"),
      values = c(0.0, 0.25, 0.50, 0.75, 0.90, 1.0),
      limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1.0),
      oob = scales::squish, name = "Spearman \u03c1", na.value = "transparent"
    ) +
    coord_sf(
      crs = 3035,
      xlim = c(min(nco[, 1]), max(nco[, 1])),
      ylim = c(min(nco[, 2]), max(nco[, 2])), expand = FALSE
    ) +
    labs(subtitle = plot_title) +
    theme_void(base_size = 10) +
    theme(
      plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.6, "cm"),
      legend.key.height = grid::unit(0.4, "cm")
    )
}

p0_1 <- build_map(map_d, "a) Daily Resolution")
p0_2 <- build_map(map_7, "b) 7-Day Moving Mean")
p0_3 <- build_map(map_15, "c) 15-Day Moving Mean")
p0_4 <- build_map(map_m, "d) Natural Calendar Month")

fig0 <- (p0_1 + plot_spacer() + p0_2 + plot_spacer() + p0_3 + plot_spacer() + p0_4) +
  plot_layout(ncol = 7, widths = c(1, 0.15, 1, 0.15, 1, 0.15, 1), guides = "collect") &
  theme(legend.position = "bottom")

fig0 <- fig0 +
  plot_annotation(
    title = "Cross-Comparison Performance Validation: LISFLOOD vs ESA CCI Surface SM",
    subtitle = "Grey = insignificant (p \u2265 0.05)  |  White = no data",
    caption = "*Negative \u03c1 values squished to 0.0",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9, colour = "grey30", face = "italic", hjust = 0.5, margin = margin(b = 6)),
      plot.caption  = element_text(size = 8, colour = "grey40", face = "italic", hjust = 0.01, margin = margin(t = 4))
    )
  )

ggsave(fig0,
  filename = file.path(path_out, "Fig0_SM_Aggregation_Spearman_Composite_1Row.png"),
  width = 20, height = 8, units = "cm", dpi = 300, bg = "white"
)
message("Fig 0 saved.")


# =============================================================================
# FIG 1  |  rho VIOLIN + BOX across 4 temporal aggregations
# =============================================================================
message("=== Fig 1: rho violin + boxplot ===")

fig1_ann <- all_meta %>%
  group_by(scenario) %>%
  summarise(
    med = round(median(rho, na.rm = TRUE), 2),
    pct = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 0),
    .groups = "drop"
  )

fig1 <- ggplot(all_meta, aes(x = scenario, y = rho, fill = scenario)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey70", linewidth = 0.3) +
  geom_violin(alpha = 0.35, linewidth = 0.3, trim = FALSE, scale = "width", colour = "grey40") +
  geom_boxplot(
    width = 0.12, outlier.shape = NA, linewidth = 0.45,
    colour = "grey25", fill = "white", fatten = 2
  ) +
  geom_text(
    data = fig1_ann,
    aes(x = scenario, y = med, label = paste0("md=", med)),
    vjust = -0.5, size = 2.2, colour = "grey15", fontface = "bold", inherit.aes = FALSE
  ) +
  geom_text(
    data = fig1_ann,
    aes(x = scenario, y = 0.97, label = paste0(pct, "% \u2265 0.5")),
    size = 3.2, colour = "grey30", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = scen_pal, guide = "none") +
  scale_y_continuous(
    limits = c(-0.35, 1.02),
    breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1.0), expand = c(0, 0)
  ) +
  labs(
    title = "Fig. 1 | Spearman \u03c1 distributions across four temporal aggregation windows",
    subtitle = paste0(
      "Spearman \u03c1 across all European catchments (n = ",
      nrow(filter(all_meta, scenario == "Daily")), ")"
    ),
    x = NULL, y = "Spearman \u03c1",
    caption = "Dashed: \u03c1 = 0.5. Boxes: IQR + median. Whiskers: 1.5 \u00d7 IQR."
  ) +
  theme_nature(base_size = 9)

ggsave(fig1,
  filename = file.path(path_out, "Fig1_rho_violin_aggregations.png"),
  width = 16, height = 12, units = "cm", dpi = 300, bg = "white"
)
message("Fig 1 saved.")


# =============================================================================
# FIG 2  |  SCATTER PANELS — 4 aggregations, coloured by climate class
# =============================================================================
message("=== Fig 2: Scatter panels by climate class ===")

ann_overall <- plot_dt[!is.na(clim_class), .(
  rho_all = round(cor(sm_mod, sm_obs, method = "spearman", use = "complete.obs"), 2)
), by = scenario]

fig2 <- ggplot(plot_dt_sub, aes(x = sm_mod, y = sm_obs, colour = clim_class)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "longdash", linewidth = 0.25) +
  geom_point(alpha = 0.45, size = 0.6, stroke = 0, shape = 16) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.75) +
  geom_text(
    data = ann_overall,
    aes(
      x = ax_lim[1] + diff(ax_lim) * 0.04,
      y = ax_lim[2] - diff(ax_lim) * 0.06,
      label = paste0("\u03c1 = ", rho_all)
    ),
    hjust = 0, vjust = 1, size = 2.5, colour = "grey15",
    fontface = "bold", inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = clim_pal, name = "Climate class",
    guide = guide_legend(override.aes = list(alpha = 1, size = 2), nrow = 1)
  ) +
  scale_x_continuous(limits = ax_lim, labels = label_number(accuracy = 0.01)) +
  scale_y_continuous(limits = ax_lim, labels = label_number(accuracy = 0.01)) +
  facet_wrap(~scenario, ncol = 4) +
  coord_fixed(ratio = 1) +
  labs(
    title    = "Fig. 2 | SM agreement by temporal aggregation and climate class",
    subtitle = "Each point = one catchment \u00d7 time-step. Lines = per-climate OLS. Dashed = 1:1.",
    x        = "LISFLOOD SM (m\u00b3/m\u00b3)",
    y        = "ESA CCI SM (m\u00b3/m\u00b3)",
    caption  = "Subsample: max 3 000 pts per climate \u00d7 scenario. \u03c1 = overall Spearman."
  ) +
  theme_nature(base_size = 9) +
  theme(legend.position = "bottom", panel.spacing = unit(1.0, "lines"))

ggsave(fig2,
  filename = file.path(path_out, "Fig2_scatter_climate_4panels.png"),
  width = 20, height = 10, units = "cm", dpi = 300, bg = "white"
)
message("Fig 2 saved.")


# =============================================================================
# FIG 3  |  SCATTER MATRIX — climate rows x aggregation columns
# =============================================================================
message("=== Fig 3: Scatter matrix ===")

fig3 <- ggplot(plot_dt_sub, aes(x = sm_mod, y = sm_obs, colour = clim_class)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "longdash", linewidth = 0.3) +
  geom_point(alpha = 0.50, size = 0.8, stroke = 0, shape = 16) +
  geom_smooth(aes(group = clim_class),
    method = "lm", formula = y ~ x,
    se = FALSE, colour = "grey15", linewidth = 0.7
  ) +
  geom_text(
    data = ann, aes(label = label),
    x = ax_lim[1] + diff(ax_lim) * 0.04,
    y = ax_lim[2] - diff(ax_lim) * 0.04,
    hjust = 0, vjust = 1, size = 2.8, colour = "grey15",
    lineheight = 1.3, inherit.aes = FALSE
  ) +
  scale_colour_manual(values = clim_pal, guide = "none") +
  scale_x_continuous(
    limits = ax_lim, labels = label_number(accuracy = 0.01),
    breaks = pretty(ax_lim, n = 3)
  ) +
  scale_y_continuous(
    limits = ax_lim, labels = label_number(accuracy = 0.01),
    breaks = pretty(ax_lim, n = 3)
  ) +
  coord_fixed(ratio = 1, clip = "on") +
  facet_grid(clim_class ~ scenario) +
  labs(
    title    = "Fig. 3 | Agreement matrix: climate class \u00d7 temporal aggregation",
    subtitle = "Black line = OLS per panel. Dashed = 1:1. \u03c1 and n annotated per panel.",
    x        = "LISFLOOD SM (m\u00b3/m\u00b3)",
    y        = "ESA CCI SM (m\u00b3/m\u00b3)",
    caption  = "\u03c1 = Spearman. n = catchment \u00d7 time-step pairs in subsample."
  ) +
  theme_nature(base_size = 9) +
  theme(
    panel.spacing  = unit(0.5, "lines"),
    axis.text      = element_text(size = 7.5),
    strip.text.x   = element_text(size = 9, face = "bold"),
    strip.text.y   = element_text(size = 9, face = "bold", angle = 0),
    plot.margin    = margin(4, 4, 4, 4) # tight outer margins
  )

# Height computed so each of 5 rows gets ~equal square space:
# 4 cols x square panels -> height proportional to 5/4 of width, plus margins for strips/labels
ggsave(fig3,
  filename = file.path(path_out, "Fig3_scatter_matrix_5x4.png"),
  width = 26, height = 24, units = "cm", dpi = 600, bg = "white"
)
message("Fig 3 saved.")


# =============================================================================
# FIG 4  |  STRATIFIED DIAGNOSTICS  (3-panel boxplot composite)
# =============================================================================
message("=== Fig 4: Stratified diagnostics ===")

make_strat_box <- function(data, x_var, fill_var, fill_pal, fill_name,
                           x_lab, panel_tag, rotate_x = FALSE) {
  p <- ggplot(data, aes(x = .data[[x_var]], y = rho, fill = .data[[fill_var]])) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", colour = "grey75", linewidth = 0.25) +
    geom_boxplot(
      position = position_dodge(width = 0.8), width = 0.68,
      outlier.size = 0.2, outlier.alpha = 0.25, outlier.shape = 1,
      linewidth = 0.35, fatten = 2
    ) +
    scale_fill_manual(values = fill_pal, name = fill_name) +
    scale_y_continuous(
      limits = c(-0.3, 1.05),
      breaks = c(0, 0.25, 0.5, 0.75, 1.0), expand = c(0, 0)
    ) +
    labs(tag = panel_tag, x = x_lab, y = "Spearman \u03c1") +
    theme_nature(base_size = 8.5) +
    theme(
      legend.position = "bottom", legend.key.size = unit(0.35, "cm"),
      plot.tag = element_text(face = "bold", size = 9)
    )
  if (rotate_x) p <- p + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  p
}

p4a <- make_strat_box(
  filter(all_meta, !is.na(clim_class)),
  "scenario", "clim_class", clim_pal, "Climate class", NULL, "a"
)
p4b <- make_strat_box(filter(all_meta, !is.na(area_class)),
  "area_class", "scenario", scen_pal, "Temporal aggregation",
  "Catchment area quartile", "b",
  rotate_x = TRUE
)
p4c <- make_strat_box(filter(all_meta, !is.na(elev_class)),
  "elev_class", "scenario", scen_pal, "Temporal aggregation",
  "Elevation class", "c",
  rotate_x = TRUE
)

fig4 <- (p4a / (p4b + p4c)) +
  plot_layout(heights = c(1, 1), guides = "keep") +
  plot_annotation(
    title = "Fig. 4 | LISFLOOD-ESA CCI agreement stratified by catchment characteristics",
    caption = "Boxes: IQR; whiskers: 1.5 \u00d7 IQR; dashed: \u03c1 = 0.5.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      plot.caption = element_text(size = 7, colour = "grey50", hjust = 0)
    )
  )

ggsave(fig4,
  filename = file.path(path_out, "Fig4_stratified_diagnostics.png"),
  width = 20, height = 20, units = "cm", dpi = 600, bg = "white"
)
message("Fig 4 saved.")


# =============================================================================
# FIG 5  |  RIDGE DENSITY — rho per climate class, faceted by aggregation
# =============================================================================
message("=== Fig 5: Ridge density ===")

fig5 <- ggplot(
  filter(all_meta, !is.na(clim_class)),
  aes(x = rho, y = fct_rev(clim_class), fill = clim_class, colour = clim_class)
) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "solid", colour = "grey70", linewidth = 0.25) +
  ggridges::geom_density_ridges(
    alpha = 0.55, scale = 1.15, bandwidth = 0.04, linewidth = 0.35,
    jittered_points = TRUE, point_alpha = 0.12, point_size = 0.4,
    position = position_points_jitter(width = 0, height = 0)
  ) +
  stat_summary(fun = median, geom = "point", shape = 124, size = 4, colour = "grey15") +
  scale_fill_manual(values = clim_pal, guide = "none") +
  scale_colour_manual(values = clim_pal, guide = "none") +
  scale_x_continuous(
    limits = c(-0.3, 1.0),
    breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1.0)
  ) +
  facet_wrap(~scenario, ncol = 4) +
  labs(
    title = "Fig. 5 | Distribution of Spearman \u03c1 per climate class",
    subtitle = "Vertical bar = median. Dashed: \u03c1 = 0.5.",
    x = "Spearman \u03c1", y = NULL,
    caption = "Jittered points = individual catchment values."
  ) +
  theme_nature(base_size = 9) +
  theme(axis.text.y = element_text(size = 8.5), panel.spacing = unit(1.0, "lines"))

ggsave(fig5,
  filename = file.path(path_out, "Fig5_ridge_density_climate.png"),
  width = 22, height = 10, units = "cm", dpi = 600, bg = "white"
)
message("Fig 5 saved.")


# =============================================================================
# FIG 6  |  AGGREGATION-GAIN HEATMAP  (Delta rho vs daily baseline)
# =============================================================================
message("=== Fig 6: Aggregation-gain heatmap ===")

# Confirmed column names from diagnostic:
#   pivot_wider produces: "stratum" "15-Day" "7-Day" "Daily" "Monthly"  (alphabetical)
# Strategy: compute deltas row-wise BEFORE pivot, avoiding column-name issues entirely.

compute_gain <- function(data, group_var, group_label) {
  # Long table of median rho per stratum x scenario
  long <- data %>%
    filter(!is.na(.data[[group_var]])) %>%
    group_by(stratum = as.character(.data[[group_var]]), scenario) %>%
    summarise(med_rho = median(rho, na.rm = TRUE), .groups = "drop")

  # Extract daily baseline as a named vector
  daily_baseline <- long %>%
    filter(scenario == "Daily") %>%
    dplyr::select(stratum, baseline = med_rho)

  # Join baseline and compute delta, then keep only non-Daily rows
  long %>%
    filter(scenario != "Daily") %>%
    left_join(daily_baseline, by = "stratum") %>%
    mutate(
      delta_rho = med_rho - baseline,
      stratum_type = group_label
    ) %>%
    dplyr::select(stratum, stratum_type, aggregation = scenario, delta_rho)
}

gain_dt <- bind_rows(
  compute_gain(all_meta, "clim_class", "Climate class"),
  compute_gain(all_meta, "area_class", "Catchment area"),
  compute_gain(all_meta, "elev_class", "Elevation")
) %>%
  mutate(
    aggregation = factor(aggregation, levels = c("7-Day", "15-Day", "Monthly")),
    stratum_type = factor(stratum_type, levels = c("Climate class", "Catchment area", "Elevation")),
    stratum = factor(stratum, levels = c(
      "Polar", "Cold", "Temperate", "Arid", "Tropical",
      "Q1 (smallest)", "Q2", "Q3", "Q4 (largest)",
      "< 200 m", "200-500 m", "500-1000 m", "1000-2000 m", "> 2000 m"
    ))
  )

delta_max <- max(abs(gain_dt$delta_rho), na.rm = TRUE)

fig6 <- ggplot(gain_dt, aes(x = aggregation, y = stratum, fill = delta_rho)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(
    aes(
      label = sprintf("%+.2f", delta_rho),
      colour = abs(delta_rho) > delta_max * 0.55
    ),
    size = 2.8, fontface = "bold"
  ) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "grey20"), guide = "none") +
  scale_fill_gradientn(
    colours = c("#c2523c", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#4393c3", "#1a3f7a"),
    values = scales::rescale(
      c(-delta_max, -delta_max / 2, -0.02, 0, 0.02, delta_max / 2, delta_max)
    ),
    limits = c(-delta_max, delta_max),
    oob = squish,
    name = "\u0394\u03c1 vs daily",
    guide = guide_colourbar(barwidth = unit(0.4, "cm"), barheight = unit(5, "cm"), ticks = TRUE)
  ) +
  facet_grid(stratum_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    title = "Fig. 6 | Aggregation gain (\u0394\u03c1 vs daily) across catchment strata",
    subtitle = "Blue = improvement; red = degradation vs daily baseline.",
    x = "Temporal aggregation",
    y = NULL,
    caption = paste0(
      "Each cell = median \u03c1 at that aggregation minus median \u03c1 at daily resolution for the same stratum.\n",
      "Blue cells indicate temporal smoothing improves model-satellite agreement; red cells indicate degradation.\n",
      "Rows with near-zero \u0394\u03c1 across all columns (e.g. Arid) suggest structural disagreement ",
      "not resolved by smoothing.\n",
      "Large positive values in Temperate/Cold rows indicate daily noise is the dominant source of disagreement."
    )
  ) +
  theme_nature(base_size = 9) +
  theme(
    strip.placement   = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 8.5),
    panel.spacing     = unit(0.4, "lines"),
    axis.text.y       = element_text(size = 8),
    legend.position   = "right"
  )

ggsave(fig6,
  filename = file.path(path_out, "Fig6_aggregation_gain_heatmap.png"),
  width = 18, height = 16, units = "cm", dpi = 600, bg = "white"
)
message("Fig 6 saved.")


# =============================================================================
# SUMMARY TABLES
# =============================================================================
message("=== Exporting summary tables ===")

summarise_rho <- function(data, group_var) {
  data %>%
    filter(!is.na(.data[[group_var]])) %>%
    group_by(stratum = .data[[group_var]], scenario) %>%
    summarise(
      n_catchments = n(),
      median_rho = round(median(rho, na.rm = TRUE), 3),
      mean_rho = round(mean(rho, na.rm = TRUE), 3),
      sd_rho = round(sd(rho, na.rm = TRUE), 3),
      pct_above_0.5 = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 1),
      pct_above_0.7 = round(mean(rho >= 0.7, na.rm = TRUE) * 100, 1),
      .groups = "drop"
    )
}

write.csv(summarise_rho(all_meta, "clim_class"),
  file.path(path_out, "Table_rho_by_climate.csv"),
  row.names = FALSE
)
write.csv(summarise_rho(all_meta, "area_class"),
  file.path(path_out, "Table_rho_by_area.csv"),
  row.names = FALSE
)
write.csv(summarise_rho(all_meta, "elev_class"),
  file.path(path_out, "Table_rho_by_elevation.csv"),
  row.names = FALSE
)

message("Tables written.")
message("\n=== ALL DONE. Outputs in: ", path_out, " ===")
