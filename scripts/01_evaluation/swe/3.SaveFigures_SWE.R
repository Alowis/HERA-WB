###############################################################################
# SWE CROSS-COMPARISON - STEP 3: FIGURES  (GlobSnow vs LISFLOOD)
#
# Mirrors Diego_script/3.SaveFigures.R, adapted for SWE. Produces:
#   Fig0  Spatial rho maps (Daily / 7-Day / 15-Day / Monthly)
#   Fig1  rho violin + box across the four aggregations
#   Fig2  Scatter panels (4 aggregations) coloured by climate class
#   Fig3  Scatter matrix (climate rows x aggregation columns)
#   Fig4  Stratified boxplots (climate / area / elevation)
#   Fig5  Ridge density of rho per climate class
#   Fig6  Aggregation-gain heatmap (delta-rho vs daily)
#   Tables: rho by climate / area / elevation
#
# External rasters (now available in data/):
#   Climate : data/koppen_geiger_0p1.tif
#   DEM     : data/dem.nc
# The NUTS basemap from Diego's script is replaced by rnaturalearth.
###############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(exactextractr)
library(ggplot2)
library(patchwork)
library(scales)
library(rnaturalearth)

# 0. PATHS ---------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
path_stats <- file.path(base_dir, "output", "swe_diego", "2.stats")
path_out <- file.path(base_dir, "output", "swe_diego", "3.figures")
file_shp <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
path_clim <- file.path(base_dir, "data", "koppen_geiger_0p1.tif")
path_dem <- file.path(base_dir, "data", "dem.nc")
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

norm_id <- function(x) {
    x <- as.character(x)
    suppressWarnings(ifelse(grepl("^\\d+$", x), as.character(as.integer(x)), x))
}

# 1. LOAD STATS + TIME SERIES + CATCHMENTS -------------------------
message("Loading stats, time series and catchments...")
swe_ts <- readRDS(file.path(path_stats, "SWE_time_series_all_scenarios.rds"))
stats_daily <- readRDS(file.path(path_stats, "stats_daily.rds"))
stats_7d <- readRDS(file.path(path_stats, "stats_7d.rds"))
stats_15d <- readRDS(file.path(path_stats, "stats_15d.rds"))
stats_month <- readRDS(file.path(path_stats, "stats_month.rds"))

shp <- st_read(file_shp, quiet = TRUE)
shp$join_id <- norm_id(shp$catch_id)

strip_id <- function(df) {
    df <- as.data.frame(df)
    df$join_id <- norm_id(df$catch_id)
    df
}

map_d <- left_join(shp, strip_id(stats_daily), by = "join_id")
map_7 <- left_join(shp, strip_id(stats_7d), by = "join_id")
map_15 <- left_join(shp, strip_id(stats_15d), by = "join_id")
map_m <- left_join(shp, strip_id(stats_month), by = "join_id")

hist(map_d$rho)
length(which(!is.na(map_d$rho)))
length(which(map_d$status=="Not Significant"))


# 2. CLIMATE + ELEVATION + AREA ENRICHMENT -------------------------
message("Enriching catchments with climate / elevation / area...")

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

r_dem <- terra::rast(path_dem)
shp$elev_m <- exactextractr::exact_extract(
    r_dem, sf::st_transform(shp, terra::crs(r_dem)),
    fun = "mean"
)

shp <- shp %>%
    mutate(
        clim_class = factor(
            clim_class,
            levels = c("Polar", "Cold", "Temperate", "Arid", "Tropical")
        ),
        elev_class = cut(
            elev_m,
            breaks = c(-Inf, 200, 500, 1000, 2000, Inf),
            labels = c("< 200 m", "200-500 m", "500-1000 m", "1000-2000 m", "> 2000 m")
        ),
        area_class = cut(
            area_km2,
            breaks = quantile(area_km2, probs = c(0, .25, .50, .75, 1), na.rm = TRUE),
            labels = c("Q1 (smallest)", "Q2", "Q3", "Q4 (largest)"),
            include.lowest = TRUE
        )
    )

meta_lut <- st_drop_geometry(shp) %>%
    dplyr::select(join_id, clim_class, area_km2, area_class, elev_m, elev_class)

# 3. LONG rho TABLE ACROSS SCENARIOS -------------------------------
attach_meta <- function(stats_df, label) {
    strip_id(stats_df) %>%
        left_join(meta_lut, by = "join_id") %>%
        filter(!is.na(rho)) %>%
        mutate(scenario = label)
}
all_meta <- bind_rows(
    attach_meta(stats_daily, "Daily"),
    attach_meta(stats_7d, "7-Day"),
    attach_meta(stats_15d, "15-Day"),
    attach_meta(stats_month, "Monthly")
) %>%
    mutate(scenario = factor(scenario, levels = c("Daily", "7-Day", "15-Day", "Monthly")))

# 4. LONG SCATTER TABLE  (obs x mod per catchment x time step) -----
message("Building scatter table...")
clim_lut_dt <- as.data.table(meta_lut)[, .(id_key = join_id, clim_class)]

build_scatter_dt <- function(obs_wide, mod_wide, label) {
    obs_wide <- as.data.table(obs_wide)
    mod_wide <- as.data.table(mod_wide)
    o <- melt(obs_wide, id.vars = "date", variable.name = "id_key", value.name = "obs")
    m <- melt(mod_wide, id.vars = "date", variable.name = "id_key", value.name = "mod")
    o[, id_key := as.character(id_key)]
    m[, id_key := as.character(id_key)]
    dt <- merge(m, o, by = c("date", "id_key"))
    dt <- dt[!is.na(obs) & !is.na(mod)]
    dt <- merge(dt, clim_lut_dt, by = "id_key", all.x = TRUE)
    dt[, date := as.character(date)] # unify type across windows
    dt[, scenario := label]
    dt
}

plot_dt <- rbind(
    build_scatter_dt(swe_ts$daily$obs, swe_ts$daily$mod, "a) Daily"),
    build_scatter_dt(swe_ts$`7d`$obs, swe_ts$`7d`$mod, "b) 7-Day"),
    build_scatter_dt(swe_ts$`15d`$obs, swe_ts$`15d`$mod, "c) 15-Day"),
    build_scatter_dt(swe_ts$monthly$obs, swe_ts$monthly$mod, "d) Monthly")
)
plot_dt[, scenario := factor(scenario, levels = c("a) Daily", "b) 7-Day", "c) 15-Day", "d) Monthly"))]
plot_dt[, clim_class := factor(clim_class, levels = c("Polar", "Cold", "Temperate", "Arid", "Tropical"))]

ax_lim <- range(
    quantile(plot_dt$mod, probs = c(0.02, 0.98), na.rm = TRUE),
    quantile(plot_dt$obs, probs = c(0.02, 0.98), na.rm = TRUE)
)

ann <- plot_dt[!is.na(clim_class), .(
    rho = round(cor(mod, obs, method = "spearman", use = "complete.obs"), 2),
    n = .N
), by = .(scenario, clim_class)]
ann[, label := paste0("\u03c1 = ", rho, "\nn = ", format(n, big.mark = ","))]

set.seed(42)
plot_dt_sub <- plot_dt[!is.na(clim_class),
    .SD[sample(.N, min(.N, 3000L))],
    by = .(scenario, clim_class)
]

# 5. THEME + PALETTES ----------------------------------------------
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
clim_pal <- c(
    "Polar" = "#7fbfff", "Cold" = "#2166ac", "Temperate" = "#4dac26",
    "Arid" = "#d6604d", "Tropical" = "#8e0152"
)
scen_pal <- c(
    "Daily" = "#bdbdbd", "7-Day" = "#74c4e4",
    "15-Day" = "#2c9e4b", "Monthly" = "#1a3f7a"
)

# 6. FIG 0 | SPATIAL rho MAPS --------------------------------------
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

fig0 <- (build_map(map_d, "a) Daily") + build_map(map_7, "b) 7-Day Moving Mean") +
    build_map(map_15, "c) 15-Day Moving Mean") + build_map(map_m, "d) Monthly")) +
    plot_layout(ncol = 4, guides = "collect") &
    theme(legend.position = "bottom")
fig0 <- fig0 + plot_annotation(
    title = "SWE cross-comparison: LISFLOOD vs GlobSnow (Spearman \u03c1)",
    subtitle = "Grey = no Globsnow",
    theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0.5)
    )
)
ggsave(file.path(path_out, "Fig0_SWE_spatial_rho.png"), fig0,
    width = 40, height = 13, units = "cm", dpi = 300, bg = "white"
)

  # 7. FIG 1 | rho VIOLIN + BOX --------------------------------------
message("Fig 1: rho violin...")
fig1_ann <- all_meta %>%
    group_by(scenario) %>%
    summarise(
        med = round(median(rho, na.rm = TRUE), 2),
        pct = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 0), .groups = "drop"
    )
fig1 <- ggplot(all_meta, aes(scenario, rho, fill = scenario)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", colour = "grey40") +
    geom_boxplot(width = 0.12, outlier.shape = NA, colour = "grey25", fill = "white", fatten = 2) +
    geom_text(
        data = fig1_ann, aes(scenario, med, label = paste0("md=", med)),
        vjust = -0.6, size = 2.6, fontface = "bold", inherit.aes = FALSE
    ) +
    geom_text(
        data = fig1_ann, aes(scenario, 0.97, label = paste0(pct, "% \u2265 0.5")),
        size = 3, colour = "grey30", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = scen_pal, guide = "none") +
    scale_y_continuous(limits = c(-0.35, 1.02), breaks = seq(-0.25, 1, 0.25)) +
    labs(title = "Fig. 1 | Spearman \u03c1 across temporal aggregations (SWE)", x = NULL, y = "Spearman \u03c1") +
    theme_nature()
ggsave(file.path(path_out, "Fig1_SWE_rho_violin.png"), fig1,
    width = 16, height = 12, units = "cm", dpi = 300, bg = "white"
)

# 8. FIG 2 | SCATTER PANELS BY CLIMATE -----------------------------
message("Fig 2: scatter panels by climate...")
ann_overall <- plot_dt[!is.na(clim_class), .(
    rho_all = round(cor(mod, obs, method = "spearman", use = "complete.obs"), 2)
), by = scenario]

fig2 <- ggplot(plot_dt_sub, aes(mod, obs, colour = clim_class)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "longdash", linewidth = 0.25) +
    geom_point(alpha = 0.45, size = 0.6, stroke = 0, shape = 16) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.75) +
    geom_text(
        data = ann_overall,
        aes(
            x = ax_lim[1] + diff(ax_lim) * 0.04, y = ax_lim[2] - diff(ax_lim) * 0.06,
            label = paste0("\u03c1 = ", rho_all)
        ),
        hjust = 0, vjust = 1, size = 2.5, colour = "grey15", fontface = "bold", inherit.aes = FALSE
    ) +
    scale_colour_manual(
        values = clim_pal, name = "Climate class",
        guide = guide_legend(override.aes = list(alpha = 1, size = 2), nrow = 1)
    ) +
    scale_x_continuous(limits = ax_lim) +
    scale_y_continuous(limits = ax_lim) +
    facet_wrap(~scenario, ncol = 4) +
    coord_fixed(ratio = 1) +
    labs(
        title = "Fig. 2 | SWE agreement by aggregation and climate class",
        subtitle = "Each point = one catchment \u00d7 time-step. Lines = per-climate OLS. Dashed = 1:1.",
        x = "LISFLOOD SWE (mm)", y = "GlobSnow SWE (mm)"
    ) +
    theme_nature() +
    theme(legend.position = "bottom", panel.spacing = unit(1.0, "lines"))
ggsave(file.path(path_out, "Fig2_SWE_scatter_climate.png"), fig2,
    width = 22, height = 10, units = "cm", dpi = 300, bg = "white"
)

# 9. FIG 3 | SCATTER MATRIX climate x aggregation ------------------
message("Fig 3: scatter matrix...")
fig3 <- ggplot(plot_dt_sub, aes(mod, obs, colour = clim_class)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "longdash", linewidth = 0.3) +
    geom_point(alpha = 0.50, size = 0.8, stroke = 0, shape = 16) +
    geom_smooth(aes(group = clim_class), method = "lm", formula = y ~ x, se = FALSE, colour = "grey15", linewidth = 0.7) +
    geom_text(
        data = ann, aes(label = label),
        x = ax_lim[1] + diff(ax_lim) * 0.04, y = ax_lim[2] - diff(ax_lim) * 0.04,
        hjust = 0, vjust = 1, size = 2.8, colour = "grey15", lineheight = 1.3, inherit.aes = FALSE
    ) +
    scale_colour_manual(values = clim_pal, guide = "none") +
    scale_x_continuous(limits = ax_lim, breaks = pretty(ax_lim, n = 3)) +
    scale_y_continuous(limits = ax_lim, breaks = pretty(ax_lim, n = 3)) +
    coord_fixed(ratio = 1, clip = "on") +
    facet_grid(clim_class ~ scenario) +
    labs(
        title = "Fig. 3 | Agreement matrix: climate class \u00d7 aggregation",
        x = "LISFLOOD SWE (mm)", y = "GlobSnow SWE (mm)"
    ) +
    theme_nature() +
    theme(
        panel.spacing = unit(0.5, "lines"),
        strip.text.y = element_text(angle = 0, face = "bold")
    )
ggsave(file.path(path_out, "Fig3_SWE_scatter_matrix.png"), fig3,
    width = 26, height = 24, units = "cm", dpi = 300, bg = "white"
)

# 10. FIG 4 | STRATIFIED BOXPLOTS (climate / area / elevation) -----
message("Fig 4: stratified diagnostics...")
make_strat_box <- function(data, x_var, fill_var, fill_pal, fill_name, x_lab, tag, rotate = FALSE) {
    p <- ggplot(data, aes(x = .data[[x_var]], y = rho, fill = .data[[fill_var]])) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
        geom_boxplot(
            position = position_dodge(width = 0.8), width = 0.68,
            outlier.size = 0.2, outlier.alpha = 0.25, outlier.shape = 1, linewidth = 0.35
        ) +
        scale_fill_manual(values = fill_pal, name = fill_name) +
        scale_y_continuous(limits = c(-0.3, 1.05), breaks = seq(0, 1, 0.25)) +
        labs(tag = tag, x = x_lab, y = "Spearman \u03c1") +
        theme_nature(8.5) +
        theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))
    if (rotate) p <- p + theme(axis.text.x = element_text(angle = 25, hjust = 1))
    p
}
p4a <- make_strat_box(filter(all_meta, !is.na(clim_class)), "scenario", "clim_class", clim_pal, "Climate class", NULL, "a")
p4b <- make_strat_box(filter(all_meta, !is.na(area_class)), "area_class", "scenario", scen_pal, "Aggregation", "Catchment area quartile", "b", TRUE)
p4c <- make_strat_box(filter(all_meta, !is.na(elev_class)), "elev_class", "scenario", scen_pal, "Aggregation", "Elevation class", "c", TRUE)
fig4 <- (p4a / (p4b + p4c)) +
    plot_annotation(title = "Fig. 4 | SWE agreement stratified by catchment characteristics")
ggsave(file.path(path_out, "Fig4_SWE_stratified.png"), fig4,
    width = 20, height = 20, units = "cm", dpi = 300, bg = "white"
)

# 11. FIG 5 | RIDGE DENSITY of rho per climate ---------------------
if (requireNamespace("ggridges", quietly = TRUE)) {
    message("Fig 5: ridge density...")
    fig5 <- ggplot(
        filter(all_meta, !is.na(clim_class)),
        aes(x = rho, y = forcats::fct_rev(clim_class), fill = clim_class, colour = clim_class)
    ) +
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55") +
        geom_vline(xintercept = 0, colour = "grey70") +
        ggridges::geom_density_ridges(alpha = 0.55, scale = 1.15, bandwidth = 0.04, linewidth = 0.35) +
        stat_summary(fun = median, geom = "point", shape = 124, size = 4, colour = "grey15") +
        scale_fill_manual(values = clim_pal, guide = "none") +
        scale_colour_manual(values = clim_pal, guide = "none") +
        scale_x_continuous(limits = c(-0.3, 1.0), breaks = seq(-0.25, 1, 0.25)) +
        facet_wrap(~scenario, ncol = 4) +
        labs(title = "Fig. 5 | Distribution of Spearman \u03c1 per climate class", x = "Spearman \u03c1", y = NULL) +
        theme_nature()
    ggsave(file.path(path_out, "Fig5_SWE_ridge_climate.png"), fig5,
        width = 22, height = 10, units = "cm", dpi = 300, bg = "white"
    )
} else {
    message("Fig 5 skipped (install 'ggridges' to enable).")
}

# 12. FIG 6 | AGGREGATION-GAIN HEATMAP (delta-rho vs daily) ---------
message("Fig 6: aggregation-gain heatmap...")
compute_gain <- function(data, group_var, group_label) {
    long <- data %>%
        filter(!is.na(.data[[group_var]])) %>%
        group_by(stratum = as.character(.data[[group_var]]), scenario) %>%
        summarise(med_rho = median(rho, na.rm = TRUE), .groups = "drop")
    baseline <- long %>%
        filter(scenario == "Daily") %>%
        dplyr::select(stratum, baseline = med_rho)
    long %>%
        filter(scenario != "Daily") %>%
        left_join(baseline, by = "stratum") %>%
        mutate(delta_rho = med_rho - baseline, stratum_type = group_label) %>%
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
fig6 <- ggplot(gain_dt, aes(aggregation, stratum, fill = delta_rho)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%+.2f", delta_rho)), size = 2.8, fontface = "bold") +
    scale_fill_gradientn(
        colours = c("#c2523c", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#4393c3", "#1a3f7a"),
        values = scales::rescale(c(-delta_max, -delta_max / 2, -0.02, 0, 0.02, delta_max / 2, delta_max)),
        limits = c(-delta_max, delta_max), oob = squish, name = "\u0394\u03c1 vs daily"
    ) +
    facet_grid(stratum_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(
        title = "Fig. 6 | Aggregation gain (\u0394\u03c1 vs daily) across strata",
        x = "Temporal aggregation", y = NULL
    ) +
    theme_nature() +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"))
ggsave(file.path(path_out, "Fig6_SWE_aggregation_gain.png"), fig6,
    width = 18, height = 16, units = "cm", dpi = 300, bg = "white"
)

# --- Scatterplot: yearly mean SWE vs rho (correlation) ------------------------

# Extract LISFLOOD monthly SWE time series (row 1 = GlobSnow, row 2 = LISFLOOD)
# Adjust row index depending on which scenario is LISFLOOD
lf_swe_monthly <- swe_ts$monthly[[1]]  # row 2 = LISFLOOD (adjust if needed)

date_col <- names(lf_swe_monthly)[1]
# Parse dates and compute annual sums, then average across years
lf_swe_monthly[, date := as.Date(paste0(get(date_col),"-01"))]
lf_swe_monthly[, year := year(date)]

catch_cols_swe <- setdiff(names(lf_swe_monthly), date_col)

# Annual sum per catchment per year, then mean across years
annual_sums <- lf_swe_monthly[, lapply(.SD, sum, na.rm = TRUE),
                              by = year, .SDcols = catch_cols_swe]
mean_swe <- colMeans(annual_sums[, ..catch_cols_swe], na.rm = TRUE)

mean_swe_df <- data.frame(
  join_id = norm_id(sub("^X", "", names(mean_swe))),
  mean_swe = as.numeric(mean_swe),
  stringsAsFactors = FALSE
)

# Get rho from stats (check column name — likely "rho" or "spearman_r")
stats_df <- strip_id(stats_month)  # or stats_daily, depending on which rho you want
# Identify the rho column
rho_col <- intersect(c("rho", "spearman_r", "r", "cor"), names(stats_df))
cat("Rho column found:", rho_col, "\n")

# Merge
scatter_df <- merge(mean_swe_df, stats_df[, c("join_id", rho_col[1])],
                    by = "join_id", all.x = TRUE)
names(scatter_df)[3] <- "rho"

scatter_df=scatter_df[which(!is.na(scatter_df$rho)),]

nonan=scatter_df$join_id
# Plot
library(ggplot2)

p_scatter <- ggplot(scatter_df, aes(x = mean_swe, y = rho)) +
  geom_point(alpha = 0.5, size = 1.5, color = "steelblue") +
  geom_smooth(method = "loess", se = TRUE, color = "darkblue") +
  labs(
    title = "Correlation (rho) vs HERA-WC mean yearly SWE per catchment",
    x = "Mean yearly SWE (mm)",
    y = expression("Monthly" ~ rho)
  ) +
 scale_x_log10(breaks=c(0.001,0.1,1,10,100,1000,10000))+
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_scatter

ggsave(file.path(path_out, "FigS_rho_vs_SWE.png"), p_scatter,
       width = 18, height = 16, units = "cm", dpi = 300, bg = "white"
)


# 2. Extract mean elevation std per catchment from NetCDF raster

elev_std_path <- file.path(base_dir, "data", "elvstd_European_01min.nc")  # adjust path if needed
r_elev_std <- rast(elev_std_path)

# Ensure catchments are in WGS84 for extraction
shp_wgs <- if (st_crs(shp) == st_crs(4326)) shp else st_transform(shp, 4326)

# Area-weighted mean of elevation std within each catchment
shp_wgs$elev_std <- exact_extract(r_elev_std, shp_wgs, "mean")


# Filter: elevation std < 50
flat_catches <- norm_id(shp_wgs$catch_id[shp_wgs$elev_std < 50])


# --- Filter catchments: snow every year + low elevation variability -----------

# 1. Identify catchments with snow every year (annual max SWE > 0 for all years)
gsw_swe_monthly <- swe_ts$monthly[[1]]  # row 2 = LISFLOOD (adjust if needed)
gsw_swe_monthly[, date := as.Date(paste0(get(date_col),"-01"))]
gsw_swe_monthly[, year := year(date)]
catch_cols_swe <- setdiff(names(gsw_swe_monthly), date_col)

annual_max_swe <- gsw_swe_monthly[, lapply(.SD, max, na.rm = TRUE),
                                 by = year, .SDcols = catch_cols_swe]

# Count years with zero snow per catchment
years_no_snow <- colSums(annual_max_swe[, ..catch_cols_swe] <= 0, na.rm = TRUE)
# Keep only catchments with snow every year
catches_with_snow <- names(years_no_snow[years_no_snow == 0])


# 3. Intersect both criteria
snow_flat <- intersect(norm_id(sub("^X", "", catches_with_snow)), flat_catches)
cat("Catchments with snow every year AND elev_std < 50:", length(snow_flat), "\n")

subsample=which(!is.na(match(nonan,snow_flat)))

scat_sample<-scatter_df[subsample,]
median(scat_sample$rho)

# 13. TABLES | median rho by stratum -------------------------------
message("Writing stratified tables...")
tbl <- function(group_var) {
    all_meta %>%
        filter(!is.na(.data[[group_var]])) %>%
        group_by(scenario, stratum = .data[[group_var]]) %>%
        summarise(
            n = dplyr::n(),
            median_rho = round(median(rho, na.rm = TRUE), 3),
            pct_ge_0.5 = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 1),
            .groups = "drop"
        )
}
data.table::fwrite(tbl("clim_class"), file.path(path_out, "Table_rho_by_climate.csv"))
data.table::fwrite(tbl("area_class"), file.path(path_out, "Table_rho_by_area.csv"))
data.table::fwrite(tbl("elev_class"), file.path(path_out, "Table_rho_by_elevation.csv"))

message("Step 3 done. Figures + tables in: ", path_out)
