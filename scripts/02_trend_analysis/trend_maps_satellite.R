# =============================================================================
# Monthly trend maps from satellite products
#   Figure A: GLEAM AET trends by month (12 panels)
#   Figure B: ESA CCI soil moisture trends by month (12 panels)
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
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
out_dir <- file.path(base_dir, "output", "temporal_evolution")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Start year (align with GLEAM=1980)
start_year <- 1980

# Input data (homogenized monthly CSVs)
gleam_path <- file.path(
    base_dir, "output", "aet_diego",
    "1.homogenized", "gleam_monthly_homog.csv"
)
cci_path <- file.path(
    base_dir, "output", "soil_moisture_diego",
    "1.Diego_Merged", "esacci_homogenized.csv"
)

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

start_year <- 1991
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
                                unit_label = "Trend", palette = "BrBG") {
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
    slope_lim <- quantile(abs(all_slopes), 0.95, na.rm = TRUE)

    # Generate 12 maps (all without legend)
    plots <- lapply(1:12, function(m) {
        make_month_map(trend_list[[m]], month_names[m], slope_lim, palette)
    })

    # Add legend to the last panel (Dec)
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

    # Compose: title + 4x3 grid
    grid <- plot_grid(plotlist = plots, ncol = 4, align = "hv")
    title_gg <- ggdraw() +
        draw_label(fig_title, fontface = "bold", size = 14, x = 0.5)

    plot_grid(title_gg, grid,
        ncol = 1,
        rel_heights = c(0.05, 1)
    )
}

# =============================================================================
# FIGURE A: GLEAM AET trends by month
# =============================================================================
cat("[A] Loading GLEAM AET data...\n")

gleam_dt <- fread(gleam_path, header = TRUE)
gleam_dt[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
gleam_dt[, month := month(date)]
gleam_dt[, year := year(date)]

# Column names have catchment IDs (some numeric, some with no prefix)
gleam_catch_cols <- setdiff(names(gleam_dt), c("date", "month", "year"))

# Convert mm/day to mm/month (multiply by days in each month)
gleam_dt[, n_days := days_in_month(date)]
gleam_dt[, (gleam_catch_cols) := lapply(.SD, function(x) x * n_days),
    .SDcols = gleam_catch_cols
]
gleam_dt[, n_days := NULL]

# Keep only needed columns
gleam_monthly <- gleam_dt[, c("year", "month", gleam_catch_cols), with = FALSE]

# Normalize column names (remove X prefix if present for matching)
clean_ids <- sub("^X", "", gleam_catch_cols)
setnames(gleam_monthly, gleam_catch_cols, clean_ids)

fig_gleam <- make_monthly_figure(gleam_monthly, clean_ids,
    fig_title = "Trends in monthly GLEAM AET (1980–2020)",
    unit_label = "Sen's Slope\n(mm/mo/yr)",
    palette = "BrBG"
)

ggsave(file.path(out_dir, "trend_maps_GLEAM_AET_by_month.png"),
    fig_gleam,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_GLEAM_AET_by_month.png\n")

# =============================================================================
# FIGURE B: ESA CCI soil moisture trends by month

# =============================================================================
# FIGURE B: ESA CCI soil moisture trends by month
# =============================================================================
cat("[B] Loading ESA CCI soil moisture data...\n")

# Load the SM time series (nested data.table: 2 rows x 4 temporal windows)
# Row 1 = ESA CCI, Row 2 = LISFLOOD
# Columns: daily, 7d, 15d, monthly (each a nested data.table)
cci_ts_path <- file.path(
    base_dir, "output", "soil_moisture_diego",
    "2.Diego_Analysis", "0.Stats_time_windows", "SM_time_series_all_scenarios.rds"
)

cci_ts_all <- readRDS(cci_ts_path)

# Extract ESA CCI monthly data (row 1)
cci_monthly_dt <- cci_ts_all$monthly[[2]]
cat(
    "  ESA CCI monthly:", nrow(cci_monthly_dt), "rows x",
    ncol(cci_monthly_dt), "cols\n"
)

# First column is the date
date_col <- names(cci_monthly_dt)[1]
cci_monthly_dt[, date := as.Date(paste0(get(date_col),"-01"))] 
cci_monthly_dt[, month := month(date)]
cci_monthly_dt[, year := year(date)]

# Catchment columns
meta_cci <- c(date_col, "date", "month", "year")
catch_cols_cci <- setdiff(names(cci_monthly_dt), meta_cci)

# Remove X prefix if present
clean_ids_cci <- sub("^X", "", catch_cols_cci)
if (any(catch_cols_cci != clean_ids_cci)) {
    setnames(cci_monthly_dt, catch_cols_cci, clean_ids_cci)
    catch_cols_cci <- clean_ids_cci
}

cci_monthly <- cci_monthly_dt[year >= start_year,
    c("year", "month", catch_cols_cci),
    with = FALSE
]

fig_cci <- make_monthly_figure(cci_monthly, catch_cols_cci,
    fig_title = paste0("Trends in monthly ESA CCI soil moisture (", start_year, "\u20132020)"),
    unit_label = "Sen's Slope\n(m3/m3/yr)",
    palette = "BrBG"
)

ggsave(file.path(out_dir, paste0("trend_maps_ESACCI_SM_by_month",start_year,".png")),
    fig_cci,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved trend_maps_ESACCI_SM_by_month.png\n")

cat("\nDone! Satellite trend maps saved to:", out_dir, "\n")
