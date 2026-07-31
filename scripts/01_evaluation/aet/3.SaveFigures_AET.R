###############################################################################
# AET DIEGO-STYLE - STEP 3: FIGURES  (GLEAM vs LISFLOOD)
#
# Produces Fig0-Fig6 + stratified tables (climate, area, elevation),
# mirroring the SWE / SM Diego workflow but for actual evapotranspiration.
#
# External rasters:
#   Climate : data/koppen_geiger_0p1.tif
#   DEM     : data/dem.nc
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
path_stats <- file.path(base_dir, "output", "aet_diego", "2.stats")
path_out <- file.path(base_dir, "output", "aet_diego", "3.figures")
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
aet_ts <- readRDS(file.path(path_stats, "AET_time_series_all_scenarios.rds"))
stats_daily <- readRDS(file.path(path_stats, "stats_daily.rds"))
stats_7d <- readRDS(file.path(path_stats, "pearson_extremes_daily.rds"))
stats_7d$rho <- stats_7d$pearson_top20
stats_15d <- readRDS(file.path(path_stats, "pearson_extremes_daily.rds"))
stats_15d$rho <- stats_15d$pearson_bot20
stats_month <- readRDS(file.path(path_stats, "stats_regime.rds"))
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

# 2. CLIMATE + ELEVATION + AREA ENRICHMENT -------------------------
message("Enriching catchments...")
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
shp <- shp %>% mutate(
    clim_class = factor(clim_class, levels = c("Polar", "Cold", "Temperate", "Arid", "Tropical")),
    elev_class = cut(elev_m,
        breaks = c(-Inf, 200, 500, 1000, 2000, Inf),
        labels = c("< 200 m", "200-500 m", "500-1000 m", "1000-2000 m", "> 2000 m")
    ),
    area_class = cut(area_km2,
        breaks = quantile(area_km2, probs = c(0, .25, .5, .75, 1), na.rm = TRUE),
        labels = c("Q1 (smallest)", "Q2", "Q3", "Q4 (largest)"), include.lowest = TRUE
    )
)
meta_lut <- st_drop_geometry(shp) %>%
    dplyr::select(join_id, clim_class, area_km2, area_class, elev_m, elev_class)

# 3. LONG rho TABLE ------------------------------------------------
attach_meta <- function(stats_df, label) {
    strip_id(stats_df) %>%
        left_join(meta_lut, by = "join_id") %>%
        filter(!is.na(rho)) %>%
        mutate(scenario = label)
}
all_meta <- bind_rows(
    attach_meta(stats_daily, "Daily"), attach_meta(stats_7d, "Top 20%"),
    attach_meta(stats_15d, "Bottom 20%"), attach_meta(stats_month, "Regime")
) %>% mutate(scenario = factor(scenario, levels = c("Daily", "Top 20%", "Bottom 20%", "Regime")))

# 4. SCATTER TABLE -------------------------------------------------
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
    dt[, date := as.character(date)]
    dt[, scenario := label]
    dt
}
plot_dt <- rbind(
    build_scatter_dt(aet_ts$daily$obs, aet_ts$daily$mod, "a) Daily"),
    build_scatter_dt(aet_ts$`7d`$obs, aet_ts$`7d`$mod, "b) 7-Day"),
    build_scatter_dt(aet_ts$`15d`$obs, aet_ts$`15d`$mod, "c) 15-Day"),
    build_scatter_dt(aet_ts$monthly$obs, aet_ts$monthly$mod, "d) Monthly")
)
plot_dt[, scenario := factor(scenario, levels = c("a) Daily", "b) Top 20%", "c) Bottom 20%", "d) Regime"))]
plot_dt[, clim_class := factor(clim_class, levels = c("Polar", "Cold", "Temperate", "Arid", "Tropical"))]
ax_lim <- range(
    quantile(plot_dt$mod, c(0.02, 0.98), na.rm = TRUE),
    quantile(plot_dt$obs, c(0.02, 0.98), na.rm = TRUE)
)
ann <- plot_dt[!is.na(clim_class), .(
    rho = round(cor(mod, obs, method = "spearman", use = "complete.obs"), 2), n = .N
), by = .(scenario, clim_class)]
ann[, label := paste0("\u03c1 = ", rho, "\nn = ", format(n, big.mark = ","))]
set.seed(42)
plot_dt_sub <- plot_dt[!is.na(clim_class), .SD[sample(.N, min(.N, 3000L))], by = .(scenario, clim_class)]

rm(aet_ts)
gc()
# 5. THEME + PALETTES ----------------------------------------------
theme_nature <- function(base_size = 9) {
    theme_bw(base_size = base_size) %+replace% theme(
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey96", colour = "grey60"),
        strip.text = element_text(face = "bold"), legend.key = element_blank()
    )
}
clim_pal <- c(
    "Polar" = "#7fbfff", "Cold" = "#2166ac", "Temperate" = "#4dac26",
    "Arid" = "#d6604d", "Tropical" = "#8e0152"
)
scen_pal <- c("Daily" = "#bdbdbd", "7-Day" = "#74c4e4", "15-Day" = "#2c9e4b", "Monthly" = "#1a3f7a")

scen_pal <- c("Daily" = "#bdbdbd", "Top 20%" = "#74c4e4", "Bottom 20%" = "#2c9e4b", "Regime" = "#1a3f7a")
# 6. FIG 0 | SPATIAL rho MAPS --------------------------------------
message("Fig 0...")
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
basemap <- sf::st_transform(world, crs = 3035)
build_map <- function(spatial_df, ttl) {
    sp <- sf::st_transform(spatial_df, 3035)
    nco <- sf::st_coordinates(sf::st_centroid(sp))
    ggplot() +
        geom_sf(data = basemap, fill = "#f7f7f7", color = NA) +
        geom_sf(data = sp, fill = "white", color = "grey92", linewidth = 0.05) +
        geom_sf(
            data = sp, aes(fill = rho),
            color = "transparent", linewidth = 0
        ) +
        # geom_sf(
        #     data = dplyr::filter(sp, status == "Not Significant"),
        #     fill = "grey75", color = "transparent"
        # ) +
        geom_sf(data = basemap, fill = NA, color = "grey20", linewidth = 0.2) +
        scale_fill_gradientn(
            colors = c("#ffffff", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"),
            values = c(0, 0.25, 0.5, 0.75, 0.9, 1), limits = c(0, 1),
            oob = squish, name = "Spearman \u03c1", na.value = "transparent"
        ) +
        coord_sf(
            crs = 3035, xlim = c(min(nco[, 1]), max(nco[, 1])),
            ylim = c(min(nco[, 2]), max(nco[, 2])), expand = FALSE
        ) +
        labs(subtitle = ttl) +
        theme_void(10) +
        theme(
            plot.subtitle = element_text(face = "bold", hjust = 0.5),
            legend.position = "bottom", legend.key.width = grid::unit(1.6, "cm")
        )
}
fig0 <- (build_map(map_d, "a) Daily") + build_map(map_7, "b) Top 20%") +
    build_map(map_15, "c) Bottom 20%") + build_map(map_m, "d) Regime")) +
    plot_layout(ncol = 4, guides = "collect") & theme(legend.position = "bottom")
fig0 <- fig0 + plot_annotation(
    title = "AET cross-comparison: HERA-WB vs GLEAM (Spearman \u03c1)",
    subtitle = "Grey = insignificant (p \u2265 0.05) | White = no data"
)
ggsave(file.path(path_out, "Fig0_AET_spatial_performances.png"), fig0,
    width = 40, height = 13, units = "cm", dpi = 300, bg = "white"
)

# 7. FIG 1 | VIOLIN ------------------------------------------------
message("Fig 1...")
fig1_ann <- all_meta %>%
    group_by(scenario) %>%
    summarise(
        med = round(median(rho, na.rm = TRUE), 2),
        pct = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 0), .groups = "drop"
    )
fig1 <- ggplot(all_meta, aes(scenario, rho, fill = scenario)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", colour = "grey40") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", fatten = 2) +
    geom_text(
        data = fig1_ann, aes(scenario, med, label = paste0("md=", med)),
        vjust = -0.6, size = 2.6, fontface = "bold", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = scen_pal, guide = "none") +
    scale_y_continuous(limits = c(-0.35, 1.02), breaks = seq(-0.25, 1, 0.25)) +
    labs(title = "Fig. 1 | Spearman \u03c1 across aggregations (AET)", x = NULL, y = "Spearman \u03c1") +
    theme_nature()
ggsave(file.path(path_out, "Fig1_AET_perf_violin.png"), fig1,
    width = 16, height = 12, units = "cm", dpi = 300, bg = "white"
)

# 8. FIG 2 | SCATTER BY CLIMATE ------------------------------------
message("Fig 2...")
ann_ov <- plot_dt[!is.na(clim_class), .(rho_all = round(cor(mod, obs, method = "spearman", use = "complete.obs"), 2)), by = scenario]
fig2 <- ggplot(plot_dt_sub, aes(mod, obs, colour = clim_class)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "longdash", linewidth = 0.25) +
    geom_point(alpha = 0.45, size = 0.6, shape = 16) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.75) +
    scale_colour_manual(
        values = clim_pal, name = "Climate",
        guide = guide_legend(override.aes = list(alpha = 1, size = 2), nrow = 1)
    ) +
    scale_x_continuous(limits = ax_lim) +
    scale_y_continuous(limits = ax_lim) +
    facet_wrap(~scenario, ncol = 4) +
    coord_fixed() +
    labs(
        title = "Fig. 2 | AET agreement by climate class",
        x = "LISFLOOD AET (mm/day)", y = "GLEAM AET (mm/day)"
    ) +
    theme_nature() +
    theme(legend.position = "bottom")
ggsave(file.path(path_out, "Fig2_AET_scatter_climate.png"), fig2,
    width = 22, height = 10, units = "cm", dpi = 300, bg = "white"
)

# 9. FIG 3 | SCATTER MATRIX ----------------------------------------
message("Fig 3...")
fig3 <- ggplot(plot_dt_sub, aes(mod, obs, colour = clim_class)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = "longdash") +
    geom_point(alpha = 0.5, size = 0.8, shape = 16) +
    geom_smooth(aes(group = clim_class), method = "lm", formula = y ~ x, se = FALSE, colour = "grey15") +
    geom_text(
        data = ann, aes(label = label), x = ax_lim[1] + diff(ax_lim) * 0.04,
        y = ax_lim[2] - diff(ax_lim) * 0.04, hjust = 0, vjust = 1, size = 2.8,
        colour = "grey15", lineheight = 1.3, inherit.aes = FALSE
    ) +
    scale_colour_manual(values = clim_pal, guide = "none") +
    scale_x_continuous(limits = ax_lim) +
    scale_y_continuous(limits = ax_lim) +
    coord_fixed(clip = "on") +
    facet_grid(clim_class ~ scenario) +
    labs(title = "Fig. 3 | AET scatter matrix", x = "LISFLOOD AET", y = "GLEAM AET") +
    theme_nature() +
    theme(strip.text.y = element_text(angle = 0, face = "bold"))
ggsave(file.path(path_out, "Fig3_AET_scatter_matrix.png"), fig3,
    width = 26, height = 24, units = "cm", dpi = 300, bg = "white"
)

# 10. FIG 4 | STRATIFIED BOXPLOTS ----------------------------------
message("Fig 4...")
make_box <- function(data, x_var, fill_var, fill_pal, fill_name, x_lab, tag, rot = FALSE) {
    p <- ggplot(data, aes(.data[[x_var]], rho, fill = .data[[fill_var]])) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
        geom_boxplot(position = position_dodge(0.8), width = 0.68, outlier.size = 0.2, outlier.shape = 1) +
        scale_fill_manual(values = fill_pal, name = fill_name) +
        scale_y_continuous(limits = c(-0.3, 1.05), breaks = seq(0, 1, 0.25)) +
        labs(tag = tag, x = x_lab, y = "Spearman \u03c1") +
        theme_nature(8.5) +
        theme(legend.position = "bottom")
    if (rot) p <- p + theme(axis.text.x = element_text(angle = 25, hjust = 1))
    p
}
p4a <- make_box(filter(all_meta, !is.na(clim_class)), "scenario", "clim_class", clim_pal, "Climate", NULL, "a")
p4b <- make_box(filter(all_meta, !is.na(area_class)), "area_class", "scenario", scen_pal, "Aggregation", "Area quartile", "b", TRUE)
p4c <- make_box(filter(all_meta, !is.na(elev_class)), "elev_class", "scenario", scen_pal, "Aggregation", "Elevation", "c", TRUE)
fig4 <- (p4a / (p4b + p4c)) + plot_annotation(title = "Fig. 4 | AET stratified diagnostics")
ggsave(file.path(path_out, "Fig4_AET_stratified.png"), fig4,
    width = 20, height = 20, units = "cm", dpi = 300, bg = "white"
)

# 11. FIG 5 | RIDGE DENSITY ----------------------------------------
if (requireNamespace("ggridges", quietly = TRUE)) {
    message("Fig 5...")
    fig5 <- ggplot(
        filter(all_meta, !is.na(clim_class)),
        aes(rho, forcats::fct_rev(clim_class), fill = clim_class, colour = clim_class)
    ) +
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55") +
        ggridges::geom_density_ridges(alpha = 0.55, scale = 1.15, bandwidth = 0.04) +
        stat_summary(fun = median, geom = "point", shape = 124, size = 4, colour = "grey15") +
        scale_fill_manual(values = clim_pal, guide = "none") +
        scale_colour_manual(values = clim_pal, guide = "none") +
        scale_x_continuous(limits = c(-0.3, 1), breaks = seq(-0.25, 1, 0.25)) +
        facet_wrap(~scenario, ncol = 4) +
        labs(title = "Fig. 5 | \u03c1 density per climate (AET)", x = "Spearman \u03c1", y = NULL) +
        theme_nature()
    ggsave(file.path(path_out, "Fig5_AET_ridge_climate.png"), fig5,
        width = 22, height = 10, units = "cm", dpi = 300, bg = "white"
    )
}

# 12. FIG 6 | AGGREGATION-GAIN HEATMAP -----------------------------
message("Fig 6...")
compute_gain <- function(data, gvar, glabel) {
    long <- data %>%
        filter(!is.na(.data[[gvar]])) %>%
        group_by(stratum = as.character(.data[[gvar]]), scenario) %>%
        summarise(med_rho = median(rho, na.rm = TRUE), .groups = "drop")
    baseline <- long %>%
        filter(scenario == "Daily") %>%
        dplyr::select(stratum, baseline = med_rho)
    long %>%
        filter(scenario != "Daily") %>%
        left_join(baseline, by = "stratum") %>%
        mutate(delta_rho = med_rho - baseline, stratum_type = glabel) %>%
        dplyr::select(stratum, stratum_type, aggregation = scenario, delta_rho)
}
gain_dt <- bind_rows(
    compute_gain(all_meta, "clim_class", "Climate class"),
    compute_gain(all_meta, "area_class", "Catchment area"),
    compute_gain(all_meta, "elev_class", "Elevation")
) %>% mutate(
    aggregation = factor(aggregation, levels = c("7-Day", "15-Day", "Monthly")),
    stratum_type = factor(stratum_type, levels = c("Climate class", "Catchment area", "Elevation")),
    stratum = factor(stratum, levels = c(
        "Polar", "Cold", "Temperate", "Arid", "Tropical",
        "Q1 (smallest)", "Q2", "Q3", "Q4 (largest)",
        "< 200 m", "200-500 m", "500-1000 m", "1000-2000 m", "> 2000 m"
    ))
)
dmax <- max(abs(gain_dt$delta_rho), na.rm = TRUE)
fig6 <- ggplot(gain_dt, aes(aggregation, stratum, fill = delta_rho)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%+.2f", delta_rho)), size = 2.8, fontface = "bold") +
    scale_fill_gradientn(
        colours = c("#c2523c", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#4393c3", "#1a3f7a"),
        values = rescale(c(-dmax, -dmax / 2, -0.02, 0, 0.02, dmax / 2, dmax)),
        limits = c(-dmax, dmax), oob = squish, name = "\u0394\u03c1 vs daily"
    ) +
    facet_grid(stratum_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = "Fig. 6 | AET aggregation gain", x = "Aggregation", y = NULL) +
    theme_nature() +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"))
ggsave(file.path(path_out, "Fig6_AET_aggregation_gain.png"), fig6,
    width = 18, height = 16, units = "cm", dpi = 300, bg = "white"
)

# 13. TABLES -------------------------------------------------------
message("Tables...")
tbl <- function(gvar) {
    all_meta %>%
        filter(!is.na(.data[[gvar]])) %>%
        group_by(scenario, stratum = .data[[gvar]]) %>%
        summarise(
            n = n(), median_rho = round(median(rho, na.rm = TRUE), 3),
            pct_ge_0.5 = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 1), .groups = "drop"
        )
}
fwrite(tbl("clim_class"), file.path(path_out, "Table_rho_by_climate.csv"))
fwrite(tbl("area_class"), file.path(path_out, "Table_rho_by_area.csv"))
fwrite(tbl("elev_class"), file.path(path_out, "Table_rho_by_elevation.csv"))

message("Step 3 done. Figures + tables in: ", path_out)
