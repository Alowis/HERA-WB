# =============================================================================
# Trend maps: spatial patterns of linear trends across all catchments
#   Panel A: Yearly AET (from monthly aggregates)
#   Panel B: Winter surface soil moisture (DJFM, from monthly aggregates)
#   Panel C: Winter rainfall (DJFM, from TSS postprocess)
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(ggplot2)
library(sf)
library(dplyr)
library(lubridate)
library(scales)
library(cowplot)
library(rnaturalearth)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
agg_dir <- file.path(base_dir, "data", "aggregates")
tss_dir <- file.path(base_dir, "data", "tss_postprocess")
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
out_dir <- file.path(base_dir, "output", "temporal_evolution")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load catchments ----------------------------------------------------------
cat("Loading catchments...\n")
catchments <- st_read(gpkg_path, quiet = TRUE)
catchments_3035 <- st_transform(catchments, 3035)

# Basemap
basemap <- ne_countries(scale = "medium", returnclass = "sf") |>
    st_transform(3035)
bbox <- st_bbox(catchments_3035)

# --- Helper: compute linear trend per catchment -------------------------------
compute_trends <- function(annual_dt, catch_cols) {
    # annual_dt: data.table with 'year' column + catchment columns
    # Returns data.table with catch_id, slope, p_value
    results <- data.table(
        catch_id = catch_cols,
        slope = NA_real_,
        p_value = NA_real_
    )
    for (i in seq_along(catch_cols)) {
        cid <- catch_cols[i]
        vals <- annual_dt[[cid]]
        years <- annual_dt$year
        valid <- !is.na(vals)
        if (sum(valid) < 10) next
        fit <- lm(vals[valid] ~ years[valid])
        results[i, slope := coef(fit)[2]]
        results[i, p_value := summary(fit)$coefficients[2, 4]]
    }
    results[!is.na(slope)]
}

# --- Helper: make trend map ---------------------------------------------------
make_trend_map <- function(trend_dt, title_text, subtitle_text = NULL,
                           palette = "BrBG", unit_label = "Trend") {
    # Merge with geometries
    cats_map <- catchments_3035[as.character(as.numeric(catchments_3035$catch_id)) %in%
        trend_dt$catch_id, ]
    cats_map <- merge(cats_map, trend_dt,
        by.x = "catch_id", by.y = "catch_id",
        all.x = FALSE
    )
    cats_map$significant <- cats_map$p_value < 0.05

    # Symmetric color limits (95th percentile to avoid outlier domination)
    slope_lim <- quantile(abs(cats_map$slope), 0.95, na.rm = TRUE)
    palet <- hcl.colors(11, palette = palette, rev = FALSE)

    p <- ggplot() +
        geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
        geom_sf(
            data = cats_map[!cats_map$significant, ],
            aes(fill = slope), color = NA, alpha = 0.3
        ) +
        geom_sf(
            data = cats_map[cats_map$significant, ],
            aes(fill = slope), color = NA, alpha = 0.9
        ) +
        scale_fill_gradientn(
            colors = palet,
            limits = c(-slope_lim, slope_lim),
            oob = squish,
            name = unit_label
        ) +
        coord_sf(
            xlim = c(bbox["xmin"], bbox["xmax"]),
            ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
        ) +
        labs(title = title_text, subtitle = subtitle_text) +
        theme_minimal(base_size = 12) +
        theme(
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 10),
            legend.position = "right"
        )
    return(p)
}

# =============================================================================
# PANEL A: Yearly AET trends
# =============================================================================
cat("[A] Computing yearly AET trends...\n")

aet_path <- file.path(agg_dir, "ActEvapo", "ActEvapo_monthly_all_years.csv")
aet_dt <- fread(aet_path)
aet_dt[, date := as.Date(period_start)]
aet_dt[, n_days := as.numeric(as.Date(period_end) - as.Date(period_start))]

meta_cols <- c("month_idx", "period_start", "period_end", "date", "n_days")
catch_cols_aet <- setdiff(names(aet_dt), meta_cols)

# Filter to winter months (Dec, Jan, Feb, Mar)
aet_dt <- aet_dt[month(date) %in% c(12, 1, 2, 3)]
# Assign Dec to next winter year
aet_dt[, winter_year := fifelse(month(date) == 12, year(date) + 1L, year(date))]
# Remove incomplete edge winters
aet_dt <- aet_dt[winter_year > min(winter_year) & winter_year < max(winter_year)]

# Annual total AET per catchment (sum of monthly values, already in mm)
aet_annual <- aet_dt[, lapply(.SD, sum, na.rm = TRUE),
    by = .(year = year(date)), .SDcols = catch_cols_aet
]

trend_aet <- compute_trends(aet_annual, catch_cols_aet)
cat(
    "  AET trends: ", nrow(trend_aet), " catchments,",
    sum(trend_aet$p_value < 0.05), "significant (p<0.05)\n"
)

pA <- make_trend_map(trend_aet,
    title_text = "A) Trend in winter AET (1951–2020)",
    subtitle_text = "Linear slope per catchment. Faded = non-significant (p=0.05)",
    palette = "BrBG",
    unit_label = "Trend\n(mm/yr²)"
)

# =============================================================================
# PANEL B: Winter surface soil moisture trends (DJFM)
# =============================================================================
cat("[B] Computing winter surface SM trends...\n")

sm_path <- file.path(
    agg_dir, "surface_soil_moisture",
    "surface_soil_moisture_monthly_all_years.csv"
)
sm_dt <- fread(sm_path)
sm_dt[, date := as.Date(period_start)]

meta_cols_sm <- c("month_idx", "period_start", "period_end", "date")
catch_cols_sm <- setdiff(names(sm_dt), meta_cols_sm)

# Filter to winter months (Dec, Jan, Feb, Mar)
sm_dt <- sm_dt[month(date) %in% c(12, 1, 2, 3)]
# Assign Dec to next winter year
sm_dt[, winter_year := fifelse(month(date) == 12, year(date) + 1L, year(date))]
# Remove incomplete edge winters
sm_dt <- sm_dt[winter_year > min(winter_year) & winter_year < max(winter_year)]

# Winter mean per catchment per year
sm_winter <- sm_dt[, lapply(.SD, mean, na.rm = TRUE),
    by = .(year = winter_year), .SDcols = catch_cols_sm
]

trend_sm <- compute_trends(sm_winter, catch_cols_sm)
cat(
    "  Winter SM trends: ", nrow(trend_sm), " catchments,",
    sum(trend_sm$p_value < 0.05), "significant (p<0.05)\n"
)

pB <- make_trend_map(trend_sm,
    title_text = "B) Trend in winter surface soil moisture (DJFM, 1951–2020)",
    subtitle_text = "Linear slope per catchment. Faded = non-significant (p=0.05)",
    palette = "BrBG",
    unit_label = "Trend\n(mm/yr²)"
)

# =============================================================================
# PANEL C: Winter rainfall trends (DJFM)
# =============================================================================
cat("[C] Loading rainfall data (this may take a few minutes)...\n")

rain_path <- file.path(tss_dir, "rainUpsX_nested_1951_2020.csv")
rain_dt <- fread(rain_path, header = TRUE)

# Assign dates
n_rows <- nrow(rain_dt)
n_days <- as.integer(as.Date("2020-12-31") - as.Date("1951-01-01")) + 1

if (n_rows == n_days) {
    rain_dt[, date := seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")]
} else if (n_rows == n_days * 4) {
    # 6-hourly: aggregate to daily first
    rain_dt[, day_idx := rep(seq_len(n_days), each = 4)]
    date_seq <- seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")
    rain_dt[, date := date_seq[day_idx]]
    catch_cols_rain <- setdiff(names(rain_dt), c("day_idx", "date"))
    rain_dt <- rain_dt[, lapply(.SD, mean, na.rm = TRUE),
        by = date, .SDcols = catch_cols_rain
    ]
} else {
    warning("Unexpected row count for rainfall (", n_rows, "). Assuming daily.")
    rain_dt[, date := seq.Date(as.Date("1951-01-01"), length.out = n_rows, by = "day")]
}

catch_cols_rain <- setdiff(names(rain_dt), c("date", "day_idx"))

# Filter to winter months (Dec, Jan, Feb, Mar)
rain_dt <- rain_dt[month(date) %in% c(12, 1, 2, 3)]
rain_dt[, winter_year := fifelse(month(date) == 12, year(date) + 1L, year(date))]
rain_dt <- rain_dt[winter_year > min(winter_year) & winter_year < max(winter_year)]

# Winter total rainfall per catchment per year (sum of daily values)
cat("  Aggregating to winter totals...\n")
rain_winter <- rain_dt[, lapply(.SD, sum, na.rm = TRUE),
    by = .(year = winter_year), .SDcols = catch_cols_rain
]

trend_rain <- compute_trends(rain_winter, catch_cols_rain)
cat(
    "  Winter rainfall trends: ", nrow(trend_rain), " catchments,",
    sum(trend_rain$p_value < 0.05), "significant (p<0.05)\n"
)

pC <- make_trend_map(trend_rain,
    title_text = "C) Trend in winter rainfall (DJFM, 1951–2020)",
    subtitle_text = "Linear slope per catchment. Faded = non-significant (p=0.05)",
    palette = "BrBG",
    unit_label = "Trend\n(mm/yr²)"
)

# =============================================================================
# COMPOSE AND SAVE
# =============================================================================
cat("Composing trend maps figure...\n")

fig_trends <- plot_grid(pA, pB, pC,
    ncol = 1, align = "v",
    rel_heights = c(1, 1, 1)
)

ggsave(file.path(out_dir, "Figure_trend_maps.png"),
    fig_trends,
    width = 10, height = 20, dpi = 200
)
cat("  -> Saved Figure_trend_maps.png\n")

# Also save individual panels
ggsave(file.path(out_dir, "trend_map_AET.png"), pA,
    width = 10, height = 8, dpi = 200
)
ggsave(file.path(out_dir, "trend_map_winter_SM.png"), pB,
    width = 10, height = 8, dpi = 200
)
ggsave(file.path(out_dir, "trend_map_winter_rainfall.png"), pC,
    width = 10, height = 8, dpi = 200
)

cat("\nDone! Trend maps saved to:", out_dir, "\n")
