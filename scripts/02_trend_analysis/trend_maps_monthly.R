# =============================================================================
# Monthly trend maps: one map per month (12 panels) for each variable
#   Figure A: AET trends by month
#   Figure B: Surface soil moisture trends by month
#   Figure C: Rainfall trends by month
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

# --- Start year (align with satellite products: GLEAM=1980, ESA CCI=1978) -----
start_year <- 1951

# --- Load catchments ----------------------------------------------------------
cat("Loading catchments...\n")
catchments <- st_read(gpkg_path, quiet = TRUE)
catchments_3035 <- st_transform(catchments, 3035)

basemap <- ne_countries(scale = "medium", returnclass = "sf") |>
    st_transform(3035)
bbox <- st_bbox(catchments_3035)

month_names <- c(
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)

# --- Helper: compute Sen's slope + Mann-Kendall p-value per catchment ---------
compute_trends <- function(annual_dt, catch_cols) {
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
        y <- vals[valid]
        x <- years[valid]
        n <- length(y)
        # Sen's slope: median of all pairwise slopes
        slopes <- outer(seq_len(n), seq_len(n), function(i, j) {
            (y[i] - y[j]) / (x[i] - x[j])
        })
        slopes <- slopes[upper.tri(slopes)]
        results[i, slope := median(slopes, na.rm = TRUE)]
        # Mann-Kendall p-value
        mk <- cor.test(x, y, method = "kendall")
        results[i, p_value := mk$p.value]
    }
    results[!is.na(slope)]
}

# --- Helper: make trend map for one month -------------------------------------
make_month_map <- function(trend_dt, month_label, slope_lim, palette = "BrBG") {
    cats_map <- catchments_3035[as.character(as.numeric(catchments_3035$catch_id)) %in%
        trend_dt$catch_id, ]
    cats_map <- merge(cats_map, trend_dt,
        by.x = "catch_id", by.y = "catch_id",
        all.x = FALSE
    )
    cats_map$significant <- cats_map$p_value < 0.05

    palet <- hcl.colors(11, palette = palette, rev = FALSE)

    ggplot() +
        geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.15) +
        # geom_sf(
        #     data = cats_map[!cats_map$significant, ],
        #     aes(fill = slope), color = NA, alpha = 0.3
        # ) +
        geom_sf(
            data = cats_map,
            aes(fill = slope), color = NA, alpha = 0.6
        ) +
        scale_fill_gradientn(
            colors = palet,
            limits = c(-slope_lim, slope_lim),
            oob = squish,
            name = NULL
        ) +
        coord_sf(
            xlim = c(bbox["xmin"], bbox["xmax"]),
            ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
        ) +
        labs(title = month_label) +
        theme_void(base_size = 9) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
            legend.position = "none"
        )
}

# --- Helper: make 12-panel figure with shared legend --------------------------
make_monthly_figure <- function(monthly_data, catch_cols, fig_title,
                                unit_label = "Trend\n(mm/yr²)",
                                palette = "BrBG") {
    # Compute trend for each month
    cat("  Computing trends for 12 months...\n")
    trend_list <- list()
    for (m in 1:12) {
        dt_m <- monthly_data[month == m]
        if (nrow(dt_m) < 10) {
            trend_list[[m]] <- data.table(
                catch_id = character(0),
                slope = numeric(0),
                p_value = numeric(0)
            )
            next
        }
        trend_list[[m]] <- compute_trends(dt_m, catch_cols)
        cat("    Month", m, ":", nrow(trend_list[[m]]), "catchments\n")
    }

    # Global color limit (95th percentile across all months)
    all_slopes <- unlist(lapply(trend_list, function(x) x$slope))
    slope_lim <- quantile(abs(all_slopes), 0.98, na.rm = TRUE)

    # Generate 12 maps (all without legend)
    plots <- lapply(1:12, function(m) {
        make_month_map(trend_dt=trend_list[[m]], month_label=month_names[m], slope_lim, palette)
    })

    # Add legend to the last panel (Dec) by rebuilding it with legend visible
    palet <- hcl.colors(11, palette = palette, rev = FALSE)
    plots[[12]] <- plots[[12]] +
        scale_fill_gradientn(
            colors = palet,
            limits = c(-slope_lim, slope_lim),
            oob = squish,
            name = unit_label,
            guide = guide_colorbar(barwidth = 1, barheight = 8)
        ) +
        theme(
            legend.position = "right",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8)
        )

    # Compose: 4 cols x 3 rows + title
    grid <- plot_grid(plotlist = plots, ncol = 4, align = "hv")
    title_gg <- ggdraw() +
        draw_label(fig_title, fontface = "bold", size = 14, x = 0.5)

    plot_grid(title_gg, grid,
        ncol = 1,
        rel_heights = c(0.05, 1)
    )
}

# =============================================================================
# FIGURE A: AET trends by month
# =============================================================================
cat("[A] Loading AET data...\n")

aet_path <- file.path(agg_dir, "ActEvapo", "ActEvapo_monthly_all_years.csv")
aet_dt <- fread(aet_path)
aet_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
aet_dt[, month := month(date)]
aet_dt[, year := year(date)]

meta_cols <- c("month_idx", "period_start", "period_end", "date", "month", "year")
catch_cols_aet <- setdiff(names(aet_dt), meta_cols)

# Keep only needed columns, filter to start_year
aet_monthly <- aet_dt[year >= start_year, c("year", "month", catch_cols_aet), with = FALSE]

fig_aet <- make_monthly_figure(aet_monthly, catch_cols_aet,
    fig_title = paste0("Trends in monthly AET (", start_year, "–2020)"),
    unit_label = "Sen's Slope \n(mm/yr)",
    palette = "BrBG"
)

ggsave(file.path(out_dir, paste0("trend_maps_AET_by_month_",start_year,".png")),
    fig_aet,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_AET_by_month.png\n")

# =============================================================================
# FIGURE B: Surface soil moisture trends by month
# =============================================================================
cat("[B] Loading surface soil moisture data...\n")

sm_path <- file.path(
    agg_dir, "surface_soil_moisture",
    "surface_soil_moisture_monthly_all_years.csv"
)


sm_dt <- fread(sm_path)
# Rows are chronological: row 1 = Jan 1951
sm_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]


sm_dt[, month := month(date)]
sm_dt[, year := year(date)]

plot(diff(sm_dt[["date"]]))
meta_cols_sm <- c("month_idx", "period_start", "period_end", "date", "month", "year")
catch_cols_sm <- setdiff(names(sm_dt), meta_cols_sm)

sm_monthly <- sm_dt[year >= start_year, c("year", "month", catch_cols_sm), with = FALSE]

october=sm_monthly[month == 10,]
fig_sm <- make_monthly_figure(sm_monthly, catch_cols_sm,
    fig_title = paste0("Trends in monthly surface soil moisture (", start_year, "–2020)"),
    unit_label = "Sen's Slope\n(mm3/mm3/yr)",
    palette = "BrBG"
)

ggsave(file.path(out_dir, paste0("trend_maps_SM_by_month",start_year,".png")),
    fig_sm,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_SM_by_month.png\n")

# =============================================================================
# FIGURE C: Rainfall trends by month
# =============================================================================
cat("[C] Loading rainfall data (large file, may take a few minutes)...\n")

rain_path <- file.path(tss_dir, "rainUpsX_nested_1951_2020.csv")
rain_dt <- fread(rain_path, header = TRUE)

# Assign dates
n_rows <- nrow(rain_dt)
n_days <- as.integer(as.Date("2020-12-31") - as.Date("1951-01-01")) + 1

if (n_rows == n_days) {
    rain_dt[, date := seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")]
} else if (n_rows == n_days * 4) {
    rain_dt[, day_idx := rep(seq_len(n_days), each = 4)]
    date_seq <- seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")
    rain_dt[, date := date_seq[day_idx]]
    catch_cols_rain <- setdiff(names(rain_dt), c("day_idx", "date"))
    cat("  Aggregating 6-hourly to daily...\n")
    rain_dt <- rain_dt[, lapply(.SD, sum, na.rm = TRUE),
        by = date, .SDcols = catch_cols_rain
    ]
} else {
    rain_dt[, date := seq.Date(as.Date("1951-01-01"), length.out = n_rows, by = "day")]
}

catch_cols_rain <- setdiff(names(rain_dt), c("date", "day_idx"))

# Aggregate daily to monthly totals
cat("  Aggregating to monthly totals...\n")
rain_dt[, month := month(date)]
rain_dt[, year := year(date)]

rain_monthly <- rain_dt[year >= start_year, lapply(.SD, sum, na.rm = TRUE),
    by = .(year, month), .SDcols = catch_cols_rain
]

fig_rain <- make_monthly_figure(rain_monthly, catch_cols_rain,
    fig_title = paste0("Trends in monthly rainfall (", start_year, "–2020)"),
    unit_label = "Sen's Slope\n(mm/yr)",
    palette = "BrBG"
)

ggsave(file.path(out_dir, paste0("trend_maps_rainfall_by_month_",start_year,".png")),
    fig_rain,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_rainfall_by_month.png\n")

# =============================================================================
# FIGURE D: Root zone soil moisture trends by month
# =============================================================================
cat("[D] Loading root zone soil moisture data...\n")

rsm_path <- file.path(
    agg_dir, "root_soil_moisture",
    "root_soil_moisture_monthly_all_years.csv"
)
rsm_dt <- fread(rsm_path)
rsm_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
rsm_dt[, month := month(date)]
rsm_dt[, year := year(date)]

meta_cols_rsm <- c("month_idx", "period_start", "period_end", "date", "month", "year")
catch_cols_rsm <- setdiff(names(rsm_dt), meta_cols_rsm)

rsm_monthly <- rsm_dt[year >= start_year, c("year", "month", catch_cols_rsm), with = FALSE]

fig_rsm <- make_monthly_figure(rsm_monthly, catch_cols_rsm,
    fig_title = paste0("Trends in monthly root zone soil moisture (", start_year, "–2020)"),
    unit_label = "Sen's Slope\n(mm/yr)",
    palette = "BrBG"
)

ggsave(file.path(out_dir, paste0("trend_maps_root_SM_by_month_",start_year,".png")),
    fig_rsm,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_root_SM_by_month.png\n")

# =============================================================================
# Helper: load monthly aggregate (from pre-computed CSV or raw TSS fallback)
# =============================================================================
load_tss_monthly <- function(var_name, tss_filename, agg_method = "sum") {
    # Try pre-computed aggregate first
    agg_path <- file.path(agg_dir, var_name, paste0(var_name, "_monthly_all_years.csv"))
    if (file.exists(agg_path) & agg_method=="sum") {
        cat("  Loading pre-computed aggregate:", agg_path, "\n")
        dt <- fread(agg_path)
        dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
        dt[, month := month(date)]
        dt[, year := year(date)]
        meta <- c("month_idx", "period_start", "period_end", "date", "month", "year","V1")
        catch_cols <- setdiff(names(dt), meta)
        monthly <- dt[year >= start_year, c("year", "month", catch_cols), with = FALSE]
        return(list(data = monthly, cols = catch_cols))
    }

    # Fallback: read raw TSS
    cat("  Aggregate not found, reading raw TSS:", tss_filename, "...\n")
    tss_path <- file.path(tss_dir, tss_filename)
    dt <- fread(tss_path, header = TRUE)
    n_rows_v <- nrow(dt)
    n_days_v <- as.integer(as.Date("2020-12-31") - as.Date("1951-01-01")) + 1

    if (n_rows_v == n_days_v) {
        dt[, date := seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")]
    } else if (n_rows_v == n_days_v * 4) {
        dt[, day_idx := rep(seq_len(n_days_v), each = 4)]
        dt[, date := seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")[day_idx]]
        catch_cols <- setdiff(names(dt), c("day_idx", "date"))
        if (agg_method == "sum") {
            dt <- dt[, lapply(.SD, sum, na.rm = TRUE), by = date, .SDcols = catch_cols]
        }
        if (agg_method == "max") {
          dt <- dt[, lapply(.SD, max, na.rm = TRUE), by = date, .SDcols = catch_cols]
        } else {
            dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = date, .SDcols = catch_cols]
        }
    } else {
        dt[, date := seq.Date(as.Date("1951-01-01"), length.out = n_rows_v, by = "day")]
    }
    catch_cols <- setdiff(names(dt), c("date", "day_idx"))
    dt[, month := month(date)]
    dt[, year := year(date)]

    if (agg_method == "sum") {
        monthly <- dt[year >= start_year, lapply(.SD, sum, na.rm = TRUE),
            by = .(year, month), .SDcols = catch_cols
        ]
    } else {
        monthly <- dt[year >= start_year, lapply(.SD, mean, na.rm = TRUE),
            by = .(year, month), .SDcols = catch_cols
        ]
    }
    return(list(data = monthly, cols = catch_cols))
}

# =============================================================================
# FIGURE E: Snowfall trends by month
# =============================================================================
cat("[E] Loading snowfall data...\n")
snow_res <- load_tss_monthly(var_name="snowfall", tss_filename = "snowUpsX_nested_1951_2020.csv", agg_method =  "sum")

fig_snow <- make_monthly_figure(snow_res$data, snow_res$cols,
    fig_title = paste0("Trends in monthly snowfall (", start_year, "\u20132020)"),
    unit_label = "Sen's Slope\n(mm/yr)", palette = "BrBG"
)
ggsave(file.path(out_dir, "trend_maps_snowfall_by_month.png"),
    fig_snow,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_snowfall_by_month.png\n")

# =============================================================================
# FIGURE F: Snowmelt trends by month
# =============================================================================
cat("[F] Loading snowmelt data...\n")
melt_res <- load_tss_monthly("snowmelt", "snowMeltUpsX_nested_1951_2020.csv", "sum")

fig_melt <- make_monthly_figure(melt_res$data, melt_res$cols,
    fig_title = paste0("Trends in monthly snowmelt (", start_year, "\u20132020)"),
    unit_label = "Sen's Slope\n(mm/yr)", palette = "BrBG"
)
ggsave(file.path(out_dir, "trend_maps_snowmelt_by_month.png"),
    fig_melt,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_snowmelt_by_month.png\n")

# =============================================================================
# FIGURE G: Deep soil layer (theta3) trends by month
# =============================================================================
cat("[G] Loading deep soil layer (theta3) data...\n")
theta3_res <- load_tss_monthly("theta3", "theta3totalX_nested_1951_2020.csv", "mean")

op=data.frame(theta3_res$data)
plot(op$X57326)
fig_theta3 <- make_monthly_figure(theta3_res$data, theta3_res$cols,
    fig_title = paste0("Trends in monthly deep soil moisture (theta3) (", start_year, "\u20132020)"),
    unit_label = "Sen's Slope\n(mm/yr)", palette = "BrBG"
)
ggsave(file.path(out_dir, "trend_maps_theta3_by_month.png"),
    fig_theta3,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_theta3_by_month.png\n")

# =============================================================================
# FIGURE H: Percolation UZ to LZ trends by month
# =============================================================================
cat("[H] Loading percolation (UZ->LZ) data...\n")
perc_res <- load_tss_monthly("percolation", "percUZLZUpsX_nested_1951_2020.csv", "sum")

fig_perc <- make_monthly_figure(perc_res$data, perc_res$cols,
    fig_title = paste0("Trends in monthly percolation UZ\u2192LZ (", start_year, "\u20132020)"),
    unit_label = "Sen's Slope\n(mm/yr)", palette = "BrBG"
)
ggsave(file.path(out_dir, "trend_maps_percolation_by_month.png"),
    fig_perc,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_percolation_by_month.png\n")

cat("\nDone! All monthly trend maps saved to:", out_dir, "\n")

# =============================================================================
# FIGURE I: Direct runnoff by month
# =============================================================================
cat("[H] Loading direct runoff data...\n")
perc_res <- load_tss_monthly("runoff", "surfaceRunoffUpsX_nested_1951_2020.csv", "max")

fig_run <- make_monthly_figure(perc_res$data, perc_res$cols,
                                fig_title = paste0("Trends in monthly max surface runoff UZ\u2192LZ (", start_year, "\u20132020)"),
                                unit_label = "Sen's Slope\n(mm/yr)", palette = "BrBG"
)
ggsave(file.path(out_dir, "trend_maps_,axrunoff_by_month.png"),
       fig_run,
       width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_runoff_by_month.png\n")

cat("\nDone! All monthly trend maps saved to:", out_dir, "\n")

# =============================================================================
# FIGURE J: Direct runnoff by month
# =============================================================================

cat("[J] Loading preferential flow data...\n")
perc_res <- load_tss_monthly("preferential flow", "prefFlowUpsX_nested_1951_2020.csv", "sum")

fig_pf <- make_monthly_figure(perc_res$data, perc_res$cols,
                               fig_title = paste0("Trends in monthly preferential flow UZ\u2192LZ (", start_year, "\u20132020)"),
                               unit_label = "Sen's Slope\n(mm/yr)", palette = "BrBG"
)
ggsave(file.path(out_dir, "trend_maps_prefflow_by_month.png"),
       fig_run,
       width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_runoff_by_month.png\n")

cat("\nDone! All monthly trend maps saved to:", out_dir, "\n")

# =============================================================================
# FIGURE J: Infiltration by month
# =============================================================================

cat("[J] Loading infiltration data...\n")
perc_res <- load_tss_monthly("infiltration", "InfUpsX_nested_1951_2020.csv", "sum")

fig_inf <- make_monthly_figure(perc_res$data, perc_res$cols[-1],
                               fig_title = paste0("Trends in monthly infiltration UZ\u2192LZ (", start_year, "\u20132020)"),
                               unit_label = "Sen's Slope\n(mm/yr)", palette = "BrBG"
)
ggsave(file.path(out_dir, "trend_maps_inf_by_month.png"),
       fig_inf,
       width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_runoff_by_month.png\n")

cat("\nDone! All monthly trend maps saved to:", out_dir, "\n")