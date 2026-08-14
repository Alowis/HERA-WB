# =============================================================================
# EXAMPLE: Monthly Soil Moisture Trend Maps (Sen's slope + Mann-Kendall)
#
# Produces a 12-panel figure (4 cols x 3 rows) showing per-catchment trends
# in surface soil moisture for each calendar month, following the same logic
# as the original trend_maps_monthly.R script.
#
# Before running:
#   - Edit config/paths.R to match your environment
#   - Ensure the monthly aggregate CSV exists (run preprocessing first)
#
# Output:
#   A single PNG with 12 monthly maps + shared legend
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
source("config/paths.R")

out_dir <- file.path(base_dir, "output", "examples")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Configuration ------------------------------------------------------------
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

# =============================================================================
# Helper: compute Sen's slope + Mann-Kendall p-value per catchment
# =============================================================================
compute_trends <- function(annual_dt, catch_cols) {
    results <- data.table(
        catch_id = catch_cols,
        slope    = NA_real_,
        p_value  = NA_real_
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

        # Mann-Kendall p-value (via Kendall correlation test)
        mk <- cor.test(x, y, method = "kendall")
        results[i, p_value := mk$p.value]
    }

    results[!is.na(slope)]
}

# =============================================================================
# Helper: make trend map for one month
# =============================================================================
make_month_map <- function(trend_dt, month_label, slope_lim, palette = "BrBG") {
    cats_map <- catchments_3035[
        as.character(as.numeric(catchments_3035$catch_id)) %in% trend_dt$catch_id,
    ]
    cats_map <- merge(
        cats_map, trend_dt,
        by.x = "catch_id", by.y = "catch_id",
        all.x = FALSE
    )
    cats_map$significant <- cats_map$p_value < 0.05

    palet <- hcl.colors(11, palette = palette, rev = FALSE)

    ggplot() +
        geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.15) +
        geom_sf(
            data = cats_map,
            aes(fill = slope), color = NA, alpha = 0.6
        ) +
        scale_fill_gradientn(
            colors = palet,
            limits = c(-slope_lim, slope_lim),
            oob    = squish,
            name   = NULL
        ) +
        coord_sf(
            xlim = c(bbox["xmin"], bbox["xmax"]),
            ylim = c(bbox["ymin"], bbox["ymax"]),
            expand = FALSE
        ) +
        labs(title = month_label) +
        theme_void(base_size = 9) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
            legend.position = "none"
        )
}

# =============================================================================
# Helper: assemble 12-panel figure with shared legend
# =============================================================================
make_monthly_figure <- function(monthly_data, catch_cols, fig_title,
                                unit_label = "Sen's Slope\n(mm/yr)",
                                palette = "BrBG") {
    # Compute trend for each month
    cat("  Computing trends for 12 months...\n")
    trend_list <- list()

    for (m in 1:12) {
        dt_m <- monthly_data[month == m]
        if (nrow(dt_m) < 10) {
            trend_list[[m]] <- data.table(
                catch_id = character(0),
                slope    = numeric(0),
                p_value  = numeric(0)
            )
            next
        }
        trend_list[[m]] <- compute_trends(dt_m, catch_cols)
        cat("    Month", m, ":", nrow(trend_list[[m]]), "catchments\n")
    }

    # Global color limit (98th percentile across all months)
    all_slopes <- unlist(lapply(trend_list, function(x) x$slope))
    slope_lim <- quantile(abs(all_slopes), 0.98, na.rm = TRUE)

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
            oob    = squish,
            name   = unit_label,
            guide  = guide_colorbar(barwidth = 1, barheight = 8)
        ) +
        theme(
            legend.position = "right",
            legend.title    = element_text(size = 9),
            legend.text     = element_text(size = 8)
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
# LOAD SURFACE SOIL MOISTURE AND PRODUCE THE FIGURE
# =============================================================================
cat("Loading surface soil moisture data...\n")

sm_path <- file.path(
    agg_dir, "surface_soil_moisture",
    "surface_soil_moisture_monthly_all_years.csv"
)

if (!file.exists(sm_path)) {
    stop(
        "Monthly aggregate not found at:\n  ", sm_path,
        "\nRun the preprocessing pipeline first (scripts/00_preprocessing/)."
    )
}

sm_dt <- fread(sm_path)

# Rows are chronological: row 1 = Jan 1951
sm_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
sm_dt[, month := month(date)]
sm_dt[, year := year(date)]

meta_cols <- c("month_idx", "period_start", "period_end", "date", "month", "year")
catch_cols_sm <- setdiff(names(sm_dt), meta_cols)

sm_monthly <- sm_dt[year >= start_year, c("year", "month", catch_cols_sm), with = FALSE]

cat("  Catchments:", length(catch_cols_sm), "\n")
cat("  Years:", start_year, "- 2020\n")
cat("  Months per year: 12 x", length(unique(sm_monthly$year)), "years\n\n")

# --- Generate figure ----------------------------------------------------------
fig_sm <- make_monthly_figure(
    sm_monthly, catch_cols_sm,
    fig_title = paste0("Trends in monthly surface soil moisture (", start_year, "\u20132020)"),
    unit_label = "Sen's Slope\n(mm\u00b3/mm\u00b3/yr)",
    palette = "BrBG"
)

ggsave(
    file.path(out_dir, paste0("trend_maps_SM_by_month_", start_year, ".png")),
    fig_sm,
    width = 14, height = 12, dpi = 200
)

cat("\n=== Done ===\n")
cat("Figure saved to:", file.path(out_dir, paste0("trend_maps_SM_by_month_", start_year, ".png")), "\n")
