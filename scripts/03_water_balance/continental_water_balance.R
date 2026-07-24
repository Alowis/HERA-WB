# =============================================================================
# Continental water balance evolution (area-weighted)
#
# Panel A: Inputs — rainfall and snowmelt (positive, stacked area)
# Panel B: Fluxes — preferential flow, direct runoff (positive),
#                   infiltration, actual ET (negative/downward)
#                   percolation UZ->LZ (positive)
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(ggplot2)
library(sf)
library(lubridate)
library(cowplot)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
agg_dir <- file.path(base_dir, "data", "aggregates")
tss_dir <- file.path(base_dir, "data", "tss_postprocess")
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
out_dir <- file.path(base_dir, "output", "temporal_evolution")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load catchments for area weighting ---------------------------------------
cat("Loading catchments...\n")
catchments <- st_read(gpkg_path, quiet = TRUE)

# --- Helper: load monthly aggregate -------------------------------------------
load_monthly_agg <- function(var_name, tss_filename, agg_method = "sum") {
    agg_path <- file.path(agg_dir, var_name, paste0(var_name, "_monthly_all_years.csv"))
    if (file.exists(agg_path)) {
        dt <- fread(agg_path)
    } else {
        cat("  Aggregate not found for", var_name, "— reading TSS (slow)...\n")
        dt <- fread(file.path(tss_dir, tss_filename), header = TRUE)
        n_rows <- nrow(dt)
        n_days <- as.integer(as.Date("2020-12-31") - as.Date("1951-01-01")) + 1
        if (n_rows == n_days) {
            dt[, date := seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")]
        } else if (n_rows == n_days * 4) {
            dt[, day_idx := rep(seq_len(n_days), each = 4)]
            dt[, date := seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")[day_idx]]
            catch_cols <- setdiff(names(dt), c("day_idx", "date"))
            if (agg_method == "sum") {
                dt <- dt[, lapply(.SD, sum, na.rm = TRUE), by = date, .SDcols = catch_cols]
            } else {
                dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = date, .SDcols = catch_cols]
            }
        } else {
            dt[, date := seq.Date(as.Date("1951-01-01"), length.out = n_rows, by = "day")]
        }
        catch_cols <- setdiff(names(dt), c("date", "day_idx"))
        dt[, year := year(date)]
        dt[, month := month(date)]
        if (agg_method == "sum") {
            dt <- dt[, lapply(.SD, sum, na.rm = TRUE), by = .(year, month), .SDcols = catch_cols]
        } else {
            dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = .(year, month), .SDcols = catch_cols]
        }
        dt[, date := as.Date(paste(year, month, "15", sep = "-"))]
    }
    # Assign reliable dates from row position
    dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
    dt[, year := year(date)]
    return(dt)
}

# --- Helper: compute area-weighted annual total -------------------------------
compute_continental_annual <- function(dt, catchments_sf) {
    meta <- c("month_idx", "period_start", "period_end", "date", "year", "month", "window")
    catch_cols <- setdiff(names(dt), meta)

    # Match catchments with areas
    common_ids <- intersect(catch_cols, as.character(as.numeric(catchments_sf$catch_id)))
    area_vec <- catchments_sf$residual_area_km2[match(
        common_ids,
        as.character(as.numeric(catchments_sf$catch_id))
    )]
    weights <- area_vec / sum(area_vec, na.rm = TRUE)

    # Area-weighted mean per month (NA-robust)
    mat <- as.matrix(dt[, ..common_ids])
    monthly_wmean <- apply(mat, 1, function(row) {
        valid <- !is.na(row)
        if (sum(valid) == 0) {
            return(NA_real_)
        }
        sum(row[valid] * weights[valid]) / sum(weights[valid])
    })

    # Annual sum of monthly weighted means
    annual <- data.table(
        year = dt$year,
        monthly_val = monthly_wmean
    )[, .(annual_total = sum(monthly_val, na.rm = TRUE)), by = year]

    return(annual)
}

# --- Load all variables -------------------------------------------------------
cat("Loading variables...\n")

rain_dt <- load_monthly_agg("rainfall", "rainUpsX_nested_1951_2020.csv", "sum")
snowmelt_dt <- load_monthly_agg("snowmelt", "snowMeltUpsX_nested_1951_2020.csv", "sum")
prefflow_dt <- load_monthly_agg("prefflow", "prefFlowUpsX_nested_1951_2020.csv", "sum")
aet_dt <- load_monthly_agg("ActEvapo", "ActEvapo_nested_1951_2020.csv", "sum")
runoff_dt <- load_monthly_agg("runoff", "surfaceRunoffUpsX_nested_1951_2020.csv", "sum")
infil_dt <- load_monthly_agg("infiltration", "infUpsX_nested_1951_2020.csv", "sum")
perc_dt <- load_monthly_agg("dSubToUz", "dSubToUzUpsX_nested_1951_2020.csv", "sum")

# --- Compute continental annual totals ----------------------------------------
cat("Computing area-weighted annual totals...\n")

ann_rain <- compute_continental_annual(rain_dt, catchments)
ann_snowmelt <- compute_continental_annual(snowmelt_dt, catchments)
ann_prefflow <- compute_continental_annual(prefflow_dt, catchments)
ann_aet <- compute_continental_annual(aet_dt, catchments)
ann_runoff <- compute_continental_annual(runoff_dt, catchments)
ann_infil <- compute_continental_annual(infil_dt, catchments)
ann_perc <- compute_continental_annual(perc_dt, catchments)

# --- Panel A: Inputs (rainfall + snowmelt) ------------------------------------
inputs <- merge(ann_rain, ann_snowmelt, by = "year", suffixes = c("_rain", "_melt"))
setnames(inputs, c("year", "Rainfall", "Snowmelt"))

inputs_long <- melt(inputs,
    id.vars = "year",
    variable.name = "component", value.name = "value"
)
inputs_long[, component := factor(component, levels = c("Snowmelt", "Rainfall"))]

pA <- ggplot(inputs_long, aes(x = year, y = value, fill = component)) +
    geom_area(alpha = 0.8) +
    scale_fill_manual(
        values = c("Rainfall" = "#4393c3", "Snowmelt" = "#92c5de"),
        name = NULL
    ) +
    labs(
        title = "A) Continental inputs: Precipitation & Snowmelt",
        x = NULL, y = "Annual total (mm/year)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

# --- Panel B: Fluxes (positive & negative) ------------------------------------
fluxes <- Reduce(
    function(a, b) merge(a, b, by = "year"),
    list(
        setnames(ann_runoff, "annual_total", "Direct runoff"),
        setnames(ann_prefflow, "annual_total", "Preferential flow"),
        setnames(ann_perc, "annual_total", "Soil to UZ recharge"),
        setnames(ann_infil, "annual_total", "Infiltration"),
        setnames(ann_aet, "annual_total", "Actual ET")
    )
)

# Make infiltration and AET negative (losses going downward)
fluxes[, Infiltration := -Infiltration]
fluxes[, `Actual ET` := -`Actual ET`]

fluxes_long <- melt(fluxes,
    id.vars = "year",
    variable.name = "component", value.name = "value"
)

# Separate positive and negative for plotting
fluxes_long[, direction := fifelse(value >= 0, "positive", "negative")]

# Colors
flux_colors <- c(
    "Direct runoff" = "#ff7f0e",
    "Preferential flow" = "#9467bd",
    "Soil to UZ recharge" = "#1f77b4",
    "Infiltration" = "#2ca02c",
    "Actual ET" = "#d62728"
)

pB <- ggplot(fluxes_long, aes(x = year, y = value, fill = component)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_col(position = "stack", width = 0.8) +
    scale_fill_manual(values = flux_colors, name = NULL) +
    labs(
        title = "B) Continental fluxes (negative = losses to soil/atmosphere)",
        x = NULL, y = "Annual total (mm/year)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

# --- Compose and save ---------------------------------------------------------
fig <- plot_grid(pA, pB, ncol = 1, align = "v", rel_heights = c(1, 1.2))

ggsave(file.path(out_dir, "Figure_continental_water_balance.png"),
    fig,
    width = 12, height = 10, dpi = 200
)
cat("  -> Saved Figure_continental_water_balance.png\n")

cat("\nDone!\n")
