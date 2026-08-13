# =============================================================================
# Water balance check per catchment (REFACTORED — uses shared workspace)
#
# Water balance equation:
#   P = AET + Q + GW_loss + dS/dt
#   Residual = P - AET - Q - GW_loss  (should ˜ 0 over long term)
#
# Outputs:
#   - CSV with annual water balance per catchment
#   - Map of long-term mean residual (closure error)
#   - Map of residual as % of precipitation
#   - Temporal water balance plot (continental + example catchments)
#   - Storage trend maps (soil + GW)
# =============================================================================

# --- Load shared workspace (replaces ~50 lines of loading code) ---------------
source("R/load_workspace.R")

library(ggplot2)
library(scales)
library(cowplot)

out_dir <- file.path(base_dir, "output", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# LOAD WATER BALANCE COMPONENTS (from workspace)
# =============================================================================
cat("Extracting water balance components from workspace...\n")

# These are already loaded with correct dates, year, month columns
rain_dt <- monthly$rainfall # NULL if not yet aggregated
snow_dt <- monthly$snowfall
aet_dt <- monthly$ActEvapo
quz_dt <- monthly$quz
qlz_dt <- monthly$qlz
runoff_dt <- monthly$runoff

# GW loss may not be in workspace — load separately if needed
gwloss_path <- file.path(agg_dir, "gwloss", "gwloss_monthly_all_years.csv")
if (file.exists(gwloss_path)) {
    gwloss_dt <- fread(gwloss_path)
    gwloss_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
    gwloss_dt[, year := year(date)]
} else {
    # Try TSS fallback
    gwloss_tss <- file.path(tss_dir, "lossUpsX_nested_1951_2020.csv")
    if (file.exists(gwloss_tss)) {
        cat("  Loading GW loss from TSS (slow)...\n")
        gwloss_dt <- fread(gwloss_tss, header = TRUE)
        gwloss_dt[, date := seq.Date(as.Date("1951-01-01"), by = "day", length.out = .N)]
        gwloss_dt[, year := year(date)]
        gwloss_dt[, month := month(date)]
        catch_cols_gw <- setdiff(names(gwloss_dt), c("date", "year", "month"))
        gwloss_dt <- gwloss_dt[, lapply(.SD, sum, na.rm = TRUE),
            by = .(year, month), .SDcols = catch_cols_gw
        ]
        gwloss_dt[, date := as.Date(paste(year, month, "15", sep = "-"))]
    } else {
        gwloss_dt <- NULL
        cat("  WARNING: No GW loss data found. Proceeding without.\n")
    }
}

# --- Identify common catchment columns ----------------------------------------
meta <- c("month_idx", "period_start", "period_end", "date", "year", "month", "window")
get_catch_cols <- function(dt) setdiff(names(dt), meta)

# Check which components are available
available <- list(
    rain = rain_dt, snow = snow_dt, aet = aet_dt,
    quz = quz_dt, qlz = qlz_dt, runoff = runoff_dt
)
available <- available[!sapply(available, is.null)]

if (length(available) < 5) {
    cat("  Available components:", paste(names(available), collapse = ", "), "\n")
    cat("  Missing some. Rainfall may need preprocess_tss_aggregates.R first.\n")
}

all_catch_sets <- lapply(available, get_catch_cols)
if (!is.null(gwloss_dt)) all_catch_sets$gwloss <- get_catch_cols(gwloss_dt)

common_catches <- Reduce(intersect, all_catch_sets)
cat("  Common catchments across all variables:", length(common_catches), "\n")

# =============================================================================
# COMPUTE ANNUAL WATER BALANCE
# =============================================================================
cat("Computing annual water balance...\n")

annual_sum <- function(dt, cols) {
    dt[, lapply(.SD, sum, na.rm = TRUE), by = year, .SDcols = cols]
}

ann_rain <- annual_sum(rain_dt, common_catches)
ann_snow <- annual_sum(snow_dt, common_catches)
ann_aet <- annual_sum(aet_dt, common_catches)
ann_quz <- annual_sum(quz_dt, common_catches)
ann_qlz <- annual_sum(qlz_dt, common_catches)
ann_runoff <- annual_sum(runoff_dt, common_catches)
ann_gwloss <- if (!is.null(gwloss_dt)) annual_sum(gwloss_dt, common_catches) else NULL

# Matrix approach
years <- ann_rain$year
mat_P <- as.matrix(ann_rain[, ..common_catches]) + as.matrix(ann_snow[, ..common_catches])
mat_AET <- as.matrix(ann_aet[, ..common_catches])
mat_Q <- as.matrix(ann_quz[, ..common_catches]) +
    as.matrix(ann_qlz[, ..common_catches]) +
    as.matrix(ann_runoff[, ..common_catches])

if (!is.null(ann_gwloss)) {
    mat_GWloss <- as.matrix(ann_gwloss[, ..common_catches])
    mat_residual <- mat_P - mat_AET - mat_Q - mat_GWloss
    mean_GWloss <- colMeans(mat_GWloss, na.rm = TRUE)
} else {
    mat_GWloss <- matrix(0, nrow = nrow(mat_P), ncol = ncol(mat_P))
    mat_residual <- mat_P - mat_AET - mat_Q
    mean_GWloss <- rep(0, ncol(mat_P))
}

# Long-term means
mean_P <- colMeans(mat_P, na.rm = TRUE)
mean_AET <- colMeans(mat_AET, na.rm = TRUE)
mean_Q <- colMeans(mat_Q, na.rm = TRUE)
mean_residual <- colMeans(mat_residual, na.rm = TRUE)
mean_resid_pct <- 100 * mean_residual / mean_P

# Summary table
wb_summary <- data.table(
    catch_id = common_catches,
    mean_P_mm = mean_P,
    mean_AET_mm = mean_AET,
    mean_Q_mm = mean_Q,
    mean_GWloss_mm = mean_GWloss,
    mean_residual_mm = mean_residual,
    residual_pct_P = mean_resid_pct
)

fwrite(wb_summary, file.path(out_dir, "water_balance_summary.csv"))
cat("  Saved water_balance_summary.csv (", nrow(wb_summary), "catchments)\n")

# Print continental stats
cat("\n--- Continental water balance ---\n")
cat("  Mean P:        ", round(mean(mean_P, na.rm = TRUE), 1), "mm/year\n")
cat("  Mean AET:      ", round(mean(mean_AET, na.rm = TRUE), 1), "mm/year\n")
cat("  Mean Q:        ", round(mean(mean_Q, na.rm = TRUE), 1), "mm/year\n")
cat("  Mean GW loss:  ", round(mean(mean_GWloss, na.rm = TRUE), 1), "mm/year\n")
cat("  Mean residual: ", round(mean(mean_residual, na.rm = TRUE), 1), "mm/year\n")
cat("  Mean |resid|/P:", round(mean(abs(mean_resid_pct), na.rm = TRUE), 1), "%\n")

# =============================================================================
# TEMPORAL WATER BALANCE PLOT (continental stacked bars)
# =============================================================================
cat("Generating temporal water balance plot...\n")

# Area-weighted continental means per year (use workspace weights)
w_wb <- weights[common_catches]
w_wb <- w_wb / sum(w_wb, na.rm = TRUE)

continental_wb <- data.table(
    year = years,
    Rainfall = as.numeric(as.matrix(ann_rain[, ..common_catches]) %*% w_wb),
    Snowfall = as.numeric(as.matrix(ann_snow[, ..common_catches]) %*% w_wb),
    P = as.numeric(mat_P %*% w_wb),
    AET = as.numeric(mat_AET %*% w_wb),
    Q = as.numeric(mat_Q %*% w_wb),
    GW_loss = as.numeric(mat_GWloss %*% w_wb),
    Residual = as.numeric(mat_residual %*% w_wb)
)

# Stacked bar: Rain+Snow positive, AET/Q/GWloss negative
continental_wb[, AET_neg := -AET]
continental_wb[, Q_neg := -Q]
continental_wb[, GW_loss_neg := -GW_loss]

bar_data <- melt(continental_wb,
    id.vars = c("year", "Residual"),
    measure.vars = c("Rainfall", "Snowfall", "AET_neg", "Q_neg", "GW_loss_neg"),
    variable.name = "component", value.name = "value"
)
bar_data[, component := factor(component,
    levels = c("Snowfall", "Rainfall", "AET_neg", "Q_neg", "GW_loss_neg"),
    labels = c("Snowfall", "Rainfall", "AET", "Discharge (Q)", "GW loss")
)]

bar_colors <- c(
    "Rainfall" = "steelblue", "Snowfall" = "#92c5de",
    "AET" = "#d62728", "Discharge (Q)" = "#2ca02c", "GW loss" = "#9467bd"
)

p_continental <- ggplot() +
    geom_col(
        data = bar_data, aes(x = year, y = value, fill = component),
        width = 0.8, position = "stack"
    ) +
    geom_line(
        data = continental_wb, aes(x = year, y = Residual),
        color = "black", linewidth = 0.9
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = bar_colors, name = NULL) +
    labs(
        title = "Continental water balance (area-weighted mean)",
        subtitle = "Bars: Rain+Snow (up), AET+Q+GW_loss (down). Black line: residual",
        x = NULL, y = "Annual total (mm/year)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

ggsave(file.path(out_dir, "water_balance_temporal.png"),
    p_continental,
    width = 12, height = 6, dpi = 200
)
cat("  -> Saved water_balance_temporal.png\n")

# =============================================================================
# MAPS (using pre-loaded basemap and catchments_3035 from workspace)
# =============================================================================
cat("Generating maps...\n")

# Merge summary with geometry
cats_wb <- catchments_3035[as.character(catchments_3035$catch_id) %in% common_catches, ]
cats_wb <- merge(cats_wb, wb_summary, by.x = "catch_id", by.y = "catch_id", all.x = FALSE)

palet_div <- hcl.colors(11, palette = "RdBu", rev = TRUE)

# Residual in mm/year
resid_lim <- quantile(abs(wb_summary$mean_residual_mm), 0.95, na.rm = TRUE)

p_resid <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_wb, aes(fill = mean_residual_mm), color = NA) +
    scale_fill_gradientn(
        colors = palet_div, limits = c(-resid_lim, resid_lim),
        oob = squish, name = "Residual\n(mm/yr)"
    ) +
    coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(
        title = "Water balance residual: P - AET - Q - GW_loss",
        subtitle = "Positive = excess input. Negative = deficit."
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "right"
    )

ggsave(file.path(out_dir, "map_residual_mm.png"), p_resid,
    width = 10, height = 8, dpi = 200
)
cat("  -> Saved map_residual_mm.png\n")

cat("\nDone! Results saved to:", out_dir, "\n")
