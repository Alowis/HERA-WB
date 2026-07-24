# =============================================================================
# Water balance check per catchment
#
# Water balance equation:
#   P = AET + Q + dS/dt
#   where P = precipitation (rain + snowfall)
#         Q = total discharge (quz + qlz + direct runoff)
#         AET = actual evapotranspiration
#         dS/dt = change in storage (soil moisture + groundwater)
#
# Residual = P - AET - Q  (should ˜ dS/dt ˜ 0 over long term)
#
# Outputs:
#   - CSV with annual water balance per catchment
#   - Map of long-term mean residual (closure error)
#   - Map of residual as % of precipitation
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(ggplot2)
library(sf)
library(lubridate)
library(scales)
library(cowplot)
library(rnaturalearth)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
agg_dir <- file.path(base_dir, "data", "aggregates")
tss_dir <- file.path(base_dir, "data", "tss_postprocess")
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
out_dir <- file.path(base_dir, "output", "water_balance")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load catchments ----------------------------------------------------------
cat("Loading catchments...\n")
catchments <- st_read(gpkg_path, quiet = TRUE)

# --- Helper: load monthly aggregate ------------------------------------------
load_monthly <- function(var_name, tss_filename = NULL, agg_method = "sum") {
    agg_path <- file.path(agg_dir, var_name, paste0(var_name, "_monthly_all_years.csv"))
    if (file.exists(agg_path)) {
        cat("  Loading:", var_name, "(from aggregate)\n")
        dt <- fread(agg_path)
    } else if (!is.null(tss_filename)) {
        cat("  Loading:", var_name, "(from TSS — slow)...\n")
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
    } else {
        stop("No aggregate and no TSS file specified for: ", var_name)
    }
    # Assign reliable dates
    dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
    dt[, year := year(date)]
    return(dt)
}

# --- Load all water balance components ----------------------------------------
cat("Loading water balance components...\n")

# Precipitation = rainfall + snowfall
rain_dt <- load_monthly("rainfall", "rainUpsX_nested_1951_2020.csv", "sum")
snow_dt <- load_monthly("snowfall", "snowUpsX_nested_1951_2020.csv", "sum")

# AET
aet_dt <- load_monthly("ActEvapo", NULL, "sum")

# Q components: quz + qlz + direct runoff (= total channel discharge)
quz_dt <- load_monthly("quz", NULL, "sum")
qlz_dt <- load_monthly("qlz", NULL, "sum")
runoff_dt <- load_monthly("runoff", "surfaceRunoffUpsX_nested_1951_2020.csv", "sum")

# Groundwater losses (water lost to deep aquifers, not reaching the outlet)
gwloss_dt <- load_monthly("gwloss", "lossUpsX_nested_1951_2020.csv", "sum")

# --- Identify common catchment columns ---------------------------------------
meta <- c("month_idx", "period_start", "period_end", "date", "year", "month", "window")

get_catch_cols <- function(dt) setdiff(names(dt), meta)

# Find catchments present in all datasets
all_catch_sets <- list(
    get_catch_cols(rain_dt),
    get_catch_cols(snow_dt),
    get_catch_cols(aet_dt),
    get_catch_cols(quz_dt),
    get_catch_cols(qlz_dt),
    get_catch_cols(runoff_dt),
    get_catch_cols(gwloss_dt)
)
common_catches <- Reduce(intersect, all_catch_sets)
cat("  Common catchments across all variables:", length(common_catches), "\n")

# --- Compute annual water balance per catchment --------------------------------
cat("Computing annual water balance...\n")

# Annual sums per catchment
annual_sum <- function(dt, cols) {
    dt[, lapply(.SD, sum, na.rm = TRUE), by = year, .SDcols = cols]
}

ann_rain <- annual_sum(rain_dt, common_catches)
ann_snow <- annual_sum(snow_dt, common_catches)
ann_aet <- annual_sum(aet_dt, common_catches)
ann_quz <- annual_sum(quz_dt, common_catches)
ann_qlz <- annual_sum(qlz_dt, common_catches)
ann_runoff <- annual_sum(runoff_dt, common_catches)
ann_gwloss <- annual_sum(gwloss_dt, common_catches)

test=ann_quz[["164208"]]
plot(test)
# Compute P, Q, GW loss, and residual for each year and catchment
# P = AET + Q + GW_loss + dS/dt  =>  Residual = P - AET - Q - GW_loss
years <- ann_rain$year
n_years <- length(years)
n_catch <- length(common_catches)

# Matrix approach for speed
mat_P <- as.matrix(ann_rain[, ..common_catches]) + as.matrix(ann_snow[, ..common_catches])
mat_R <- as.matrix(ann_rain[, ..common_catches])
mat_S <- as.matrix(ann_snow[, ..common_catches])
mat_AET <- as.matrix(ann_aet[, ..common_catches])
mat_Q <- as.matrix(ann_quz[, ..common_catches]) +
    as.matrix(ann_qlz[, ..common_catches]) +
    as.matrix(ann_runoff[, ..common_catches])
mat_GWloss <- as.matrix(ann_gwloss[, ..common_catches])
mat_residual <- mat_P - mat_AET - mat_Q - mat_GWloss

test=mat_residual[,which(colnames(mat_residual)=="164208")]
plot(test)
# --- Long-term mean balance per catchment -------------------------------------
mean_P <- colMeans(mat_P, na.rm = TRUE)
mean_R <- colMeans(mat_R, na.rm = TRUE)
mean_S <- colMeans(mat_S, na.rm = TRUE)
mean_AET <- colMeans(mat_AET, na.rm = TRUE)
mean_Q <- colMeans(mat_Q, na.rm = TRUE)
mean_GWloss <- colMeans(mat_GWloss, na.rm = TRUE)
mean_residual <- colMeans(mat_residual, na.rm = TRUE)
mean_resid_pct <- 100 * mean_residual / mean_P # as % of precipitation

wb_summary <- data.table(
    catch_id = common_catches,
    mean_P_mm = mean_P,
    mean_R_mm = mean_R,
    mean_S_mm = mean_S,
    mean_AET_mm = mean_AET,
    mean_Q_mm = mean_Q,
    mean_GWloss_mm = mean_GWloss,
    mean_residual_mm = mean_residual,
    residual_pct_P = mean_resid_pct
)

# Save summary CSV
fwrite(wb_summary, file.path(out_dir, "water_balance_summary.csv"))
cat("  Saved water_balance_summary.csv (", nrow(wb_summary), "catchments)\n")

# Print continental stats
cat("\n--- Continental water balance (all catchments) ---\n")
cat("  Mean P:        ", round(mean(mean_P, na.rm = TRUE), 1), "mm/year\n")
cat("  Mean AET:      ", round(mean(mean_AET, na.rm = TRUE), 1), "mm/year\n")
cat("  Mean Q:        ", round(mean(mean_Q, na.rm = TRUE), 1), "mm/year\n")
cat("  Mean GW loss:  ", round(mean(mean_GWloss, na.rm = TRUE), 1), "mm/year\n")
cat("  Mean residual: ", round(mean(mean_residual, na.rm = TRUE), 1), "mm/year\n")
cat("  Mean |residual|/P:", round(mean(abs(mean_resid_pct), na.rm = TRUE), 1), "%\n")

# --- Temporal water balance plot (continental + example catchments) ------------
cat("Generating temporal water balance plots...\n")

# Continental: area-weighted mean of each component per year
area_vec_wb <- catchments$residual_area_km2[match(
    common_catches,
    as.character(as.numeric(catchments$catch_id))
)]
weights_wb <- area_vec_wb / sum(area_vec_wb, na.rm = TRUE)

continental_wb <- data.table(
    year = years,
    Rainfall = as.numeric(as.matrix(ann_rain[, ..common_catches]) %*% weights_wb),
    Snowfall = as.numeric(as.matrix(ann_snow[, ..common_catches]) %*% weights_wb),
    P = as.numeric(mat_P %*% weights_wb),
    AET = as.numeric(mat_AET %*% weights_wb),
    Q = as.numeric(mat_Q %*% weights_wb),
    GW_loss = as.numeric(mat_GWloss %*% weights_wb),
    Residual = as.numeric(mat_residual %*% weights_wb)
)

# Prepare data for stacked bar plot: Rain+Snow positive, AET and Q negative
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
    "Rainfall" = "steelblue",
    "Snowfall" = "#92c5de",
    "AET" = "#d62728",
    "Discharge (Q)" = "#2ca02c",
    "GW loss" = "#9467bd"
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
        subtitle = "Bars: Rain+Snow (positive), AET+Q+GW_loss (negative). Black line: residual",
        x = NULL, y = "Annual total (mm/year)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

# --- Partitioning plot: fluxes as % of P -------------------------------------
continental_wb[, pct_Rainfall := 100 * Rainfall / P]
continental_wb[, pct_Snowfall := 100 * Snowfall / P]
continental_wb[, pct_AET := 100 * AET / P]
continental_wb[, pct_Q := 100 * Q / P]
continental_wb[, pct_GWloss := 100 * GW_loss / P]

pct_data <- melt(continental_wb,
    id.vars = "year",
    measure.vars = c("pct_AET", "pct_Q", "pct_GWloss"),
    variable.name = "component", value.name = "pct"
)
pct_data[, component := factor(component,
    levels = c("pct_AET", "pct_Q", "pct_GWloss"),
    labels = c("AET / P", "Q / P", "GW loss / P")
)]

pct_colors <- c(
    "AET / P" = "#d62728",
    "Q / P" = "#2ca02c",
    "GW loss / P" = "#9467bd"
)

p_partitioning <- ggplot(pct_data, aes(x = year, y = pct, fill = component)) +
    geom_col(width = 0.8, position = "stack") +
    scale_fill_manual(values = pct_colors, name = NULL) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey30") +
    labs(
        title = "Partitioning of precipitation into fluxes (% of P)",
        subtitle = "Stacked: AET + Q + GW_loss as fraction of total precipitation",
        x = NULL, y = "Fraction of P (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

# Also show snow fraction evolution
pct_snow_data <- melt(continental_wb,
    id.vars = "year",
    measure.vars = c("pct_Snowfall", "pct_Rainfall"),
    variable.name = "component", value.name = "pct"
)
pct_snow_data[, component := factor(component,
    levels = c("pct_Snowfall", "pct_Rainfall"),
    labels = c("Snowfall / P", "Rainfall / P")
)]

p_precip_partition <- ggplot(pct_snow_data, aes(x = year, y = pct, fill = component)) +
    geom_col(width = 0.8, position = "stack") +
    scale_fill_manual(
        values = c("Snowfall / P" = "#92c5de", "Rainfall / P" = "steelblue"),
        name = NULL
    ) +
    labs(
        title = "Precipitation partitioning: rainfall vs snowfall (% of P)",
        x = NULL, y = "Fraction of P (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

# Save partitioning figures
fig_partition <- plot_grid(p_precip_partition, p_partitioning, ncol = 1, align = "v")
ggsave(file.path(out_dir, "water_balance_partitioning.png"),
    fig_partition,
    width = 12, height = 10, dpi = 200
)
cat("  -> Saved water_balance_partitioning.png\n")

# Example catchments (pick 4 diverse ones)
example_catches <- c("224577", "204769", "335731", "290666")
example_names <- c("Paris", "Upper Vistula", "Cega", "Ticino")

catch_panels <- list()
for (k in seq_along(example_catches)) {
    cid <- example_catches[k]
    if (!cid %in% common_catches) next
    cidx <- which(common_catches == cid)

    catch_wb <- data.table(
        year = years,
        Rainfall = mat_R[, cidx],
        Snowfall = mat_S[, cidx],
        AET = mat_AET[, cidx],
        Q = mat_Q[, cidx],
        GW_loss = mat_GWloss[, cidx],
        Residual = mat_residual[, cidx]
    )
    catch_wb[, AET_neg := -AET]
    catch_wb[, Q_neg := -Q]
    catch_wb[, GW_loss_neg := -GW_loss]

    catch_bar <- melt(catch_wb,
        id.vars = c("year", "Residual"),
        measure.vars = c("Rainfall", "Snowfall", "AET_neg", "Q_neg", "GW_loss_neg"),
        variable.name = "component", value.name = "value"
    )
    catch_bar[, component := factor(component,
        levels = c("Snowfall", "Rainfall", "AET_neg", "Q_neg", "GW_loss_neg"),
        labels = c("Snowfall", "Rainfall", "AET", "Discharge (Q)", "GW loss")
    )]

    catch_panels[[k]] <- ggplot() +
        geom_col(
            data = catch_bar, aes(x = year, y = value, fill = component),
            width = 0.8, position = "stack"
        ) +
        geom_line(
            data = catch_wb, aes(x = year, y = Residual),
            color = "black", linewidth = 0.7
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        scale_fill_manual(values = bar_colors, name = NULL) +
        labs(
            title = paste0(example_names[k], " (", cid, ")"),
            x = NULL, y = "mm/year"
        ) +
        theme_minimal(base_size = 10) +
        theme(
            plot.title = element_text(face = "bold", size = 10),
            legend.position = "none"
        )
}

# Compose: continental (with legend) on top, 4 catchments below
grid_catches <- plot_grid(plotlist = catch_panels, ncol = 2, align = "hv")
fig_wb_temporal <- plot_grid(
    p_continental, grid_catches,
    ncol = 1, rel_heights = c(0.45, 0.55)
)

ggsave(file.path(out_dir, "water_balance_temporal.png"),
    fig_wb_temporal,
    width = 12, height = 14, dpi = 200
)
cat("  -> Saved water_balance_temporal.png\n")

# --- Maps ---------------------------------------------------------------------
cat("Generating maps...\n")

# Merge with geometry
cats_wb <- catchments[as.character(as.numeric(catchments$catch_id)) %in% common_catches, ]
cats_wb <- merge(cats_wb, wb_summary,
    by.x = "catch_id", by.y = "catch_id", all.x = FALSE
)
cats_wb <- st_transform(cats_wb, 3035)

basemap <- ne_countries(scale = "medium", returnclass = "sf") |> st_transform(3035)
bbox <- st_bbox(cats_wb)
palet_div <- hcl.colors(11, palette = "RdBu", rev = TRUE)

# --- Map 1: Residual in mm/year ----------------------------------------------
resid_lim <- quantile(abs(wb_summary$mean_residual_mm), 0.95, na.rm = TRUE)

p_resid_mm <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_wb, aes(fill = mean_residual_mm), color = NA) +
    scale_fill_gradientn(
        colors = palet_div,
        limits = c(-resid_lim, resid_lim),
        oob = squish,
        name = "Residual\n(mm/yr)"
    ) +
    coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(
        title = "Water balance residual: P - AET - Q (mm/year)",
        subtitle = "Positive = excess input (storage gain). Negative = deficit."
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "right"
    )

# --- Map 2: Residual as % of P -----------------------------------------------
pct_lim <- quantile(abs(wb_summary$residual_pct_P), 0.95, na.rm = TRUE)

p_resid_pct <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_wb, aes(fill = residual_pct_P), color = NA) +
    scale_fill_gradientn(
        colors = palet_div,
        limits = c(-pct_lim, pct_lim),
        oob = squish,
        name = "Residual\n(% of P)"
    ) +
    coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(
        title = "Water balance residual as % of precipitation",
        subtitle = "Long-term mean (P - AET - Q) / P × 100"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "right"
    )

# --- Map 3: Mean annual P ----------------------------------------------------
p_lim <- quantile(wb_summary$mean_P_mm, c(0.05, 0.95), na.rm = TRUE)

p_precip <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_wb, aes(fill = mean_P_mm), color = NA) +
    scale_fill_viridis_c(
        option = "viridis", name = "P\n(mm/yr)",
        limits = p_lim, oob = squish
    ) +
    coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(title = "Mean annual precipitation (rain + snowfall)") +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "right"
    )

# --- Compose and save ---------------------------------------------------------
fig_wb <- plot_grid(p_precip, p_resid_mm, p_resid_pct, ncol = 1)

ggsave(file.path(out_dir, "water_balance_maps.png"),
    fig_wb,
    width = 10, height = 20, dpi = 200
)
cat("  -> Saved water_balance_maps.png\n")

# Also save individual maps
ggsave(file.path(out_dir, "map_residual_mm.png"), p_resid_mm, width = 10, height = 8, dpi = 200)
ggsave(file.path(out_dir, "map_residual_pct.png"), p_resid_pct, width = 10, height = 8, dpi = 200)
ggsave(file.path(out_dir, "map_precipitation.png"), p_precip, width = 10, height = 8, dpi = 200)

# --- Histogram of closure errors ----------------------------------------------
p_hist <- ggplot(wb_summary, aes(x = residual_pct_P)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
        title = "Distribution of water balance closure errors",
        x = "Residual (% of P)", y = "Number of catchments"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "histogram_closure_error.png"),
    p_hist,
    width = 8, height = 5, dpi = 200
)
cat("  -> Saved histogram_closure_error.png\n")

cat("\nDone! Results saved to:", out_dir, "\n")

# =============================================================================
# Map of trend in groundwater storage change (dS_GW/dt)
# GW storage change per year = (dSubToUz + prefFlow) - (quz + qlz) - GW_loss
# i.e., water entering GW minus water leaving GW
# =============================================================================
cat("Computing groundwater storage trends...\n")

# Load the additional components needed
prefflow_dt <- load_monthly("prefflow", "prefFlowUpsX_nested_1951_2020.csv", "sum")
dsubtuz_dt <- load_monthly("dSubToUz", "dSubToUzUpsX_nested_1951_2020.csv", "sum")

# Annual sums
ann_prefflow <- annual_sum(prefflow_dt, common_catches)
ann_dsubtuz <- annual_sum(dsubtuz_dt, common_catches)

# GW storage change = inflows to GW - outflows from GW
# Inflows: dSubToUz + prefFlow (water entering upper zone from soil + bypass)
# Outflows: quz + qlz + GW_loss (water leaving GW as discharge or deep loss)
mat_GW_in <- as.matrix(ann_dsubtuz[, ..common_catches]) +
    as.matrix(ann_prefflow[, ..common_catches])
mat_GW_out <- as.matrix(ann_quz[, ..common_catches]) +
    as.matrix(ann_qlz[, ..common_catches]) +
    as.matrix(ann_gwloss[, ..common_catches])
mat_dS_GW <- mat_GW_in - mat_GW_out # positive = GW gaining storage

# Compute Sen's slope trend per catchment (over 1951-2020)
cat("  Computing Sen's slope for GW storage change...\n")
gw_trend <- data.table(
    catch_id = common_catches,
    slope = NA_real_,
    p_value = NA_real_
)

for (i in seq_along(common_catches)) {
    y <- mat_dS_GW[, i]
    x <- years
    valid <- !is.na(y)
    if (sum(valid) < 10) next
    yv <- y[valid]
    xv <- x[valid]
    n <- length(yv)
    slopes <- outer(seq_len(n), seq_len(n), function(a, b) {
        (yv[a] - yv[b]) / (xv[a] - xv[b])
    })
    gw_trend[i, slope := median(slopes[upper.tri(slopes)], na.rm = TRUE)]
    mk <- cor.test(xv, yv, method = "kendall")
    gw_trend[i, p_value := mk$p.value]
}
gw_trend <- gw_trend[!is.na(slope)]

cat("  GW storage trends computed for", nrow(gw_trend), "catchments\n")
cat("  Significant (p<0.05) increasing:", sum(gw_trend$slope > 0 & gw_trend$p_value < 0.05), "\n")
cat("  Significant (p<0.05) decreasing:", sum(gw_trend$slope < 0 & gw_trend$p_value < 0.05), "\n")

# Map
cats_gw <- catchments[as.character(as.numeric(catchments$catch_id)) %in% gw_trend$catch_id, ]
cats_gw <- merge(cats_gw, gw_trend, by.x = "catch_id", by.y = "catch_id", all.x = FALSE)
cats_gw <- st_transform(cats_gw, 3035)

slope_lim_gw <- quantile(abs(gw_trend$slope), 0.95, na.rm = TRUE)
palet_gw <- hcl.colors(11, palette = "BrBG", rev = FALSE)

cats_gw$significant <- cats_gw$p_value < 0.05

p_gw_trend <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(
        data = cats_gw[!cats_gw$significant, ],
        aes(fill = slope), color = NA, alpha = 0.3
    ) +
    geom_sf(
        data = cats_gw[cats_gw$significant, ],
        aes(fill = slope), color = NA, alpha = 0.9
    ) +
    scale_fill_gradientn(
        colors = palet_gw,
        limits = c(-slope_lim_gw, slope_lim_gw),
        oob = squish,
        name = "Trend\n(mm/yr\u00b2)"
    ) +
    coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(
        title = "Trend in annual groundwater storage change (1951\u20132020)",
        subtitle = "Sen's slope of (GW inflows \u2212 GW outflows). Faded = non-significant (p\u22650.05)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
    )

ggsave(file.path(out_dir, "map_GW_storage_trend.png"),
    p_gw_trend,
    width = 10, height = 8, dpi = 200
)
cat("  -> Saved map_GW_storage_trend.png\n")

# =============================================================================
# Map of trend in SOIL storage change (dS_soil/dt)
# Soil storage change = Infiltration - AET - dSubToUz
# =============================================================================
cat("Computing soil storage trends...\n")

# Load infiltration
infil_dt <- load_monthly("infiltration", "infUpsX_nested_1951_2020.csv", "sum")
ann_infil <- annual_sum(infil_dt, common_catches)

# Soil storage change = infiltration - AET - dSubToUz
mat_infil <- as.matrix(ann_infil[, ..common_catches])
mat_dS_soil <- mat_infil - mat_AET - as.matrix(ann_dsubtuz[, ..common_catches])

# Compute Sen's slope trend per catchment
cat("  Computing Sen's slope for soil storage change...\n")
soil_trend <- data.table(
    catch_id = common_catches,
    slope = NA_real_,
    p_value = NA_real_
)

for (i in seq_along(common_catches)) {
    y <- mat_dS_soil[, i]
    x <- years
    valid <- !is.na(y)
    if (sum(valid) < 10) next
    yv <- y[valid]
    xv <- x[valid]
    n <- length(yv)
    slopes_s <- outer(seq_len(n), seq_len(n), function(a, b) {
        (yv[a] - yv[b]) / (xv[a] - xv[b])
    })
    soil_trend[i, slope := median(slopes_s[upper.tri(slopes_s)], na.rm = TRUE)]
    mk_s <- cor.test(xv, yv, method = "kendall")
    soil_trend[i, p_value := mk_s$p.value]
}
soil_trend <- soil_trend[!is.na(slope)]

cat("  Soil storage trends:", nrow(soil_trend), "catchments\n")
cat("  Significant (p<0.05) wetting:", sum(soil_trend$slope > 0 & soil_trend$p_value < 0.05), "\n")
cat("  Significant (p<0.05) drying: ", sum(soil_trend$slope < 0 & soil_trend$p_value < 0.05), "\n")

# Map: soil storage trend
cats_soil <- catchments[as.character(as.numeric(catchments$catch_id)) %in% soil_trend$catch_id, ]
cats_soil <- merge(cats_soil, soil_trend, by.x = "catch_id", by.y = "catch_id", all.x = FALSE)
cats_soil <- st_transform(cats_soil, 3035)

slope_lim_soil <- quantile(abs(soil_trend$slope), 0.95, na.rm = TRUE)
cats_soil$significant <- cats_soil$p_value < 0.05

p_soil_trend <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(
        data = cats_soil[!cats_soil$significant, ],
        aes(fill = slope), color = NA, alpha = 0.3
    ) +
    geom_sf(
        data = cats_soil[cats_soil$significant, ],
        aes(fill = slope), color = NA, alpha = 0.9
    ) +
    scale_fill_gradientn(
        colors = palet_gw,
        limits = c(-slope_lim_soil, slope_lim_soil),
        oob = squish,
        name = "Trend\n(mm/yr\u00b2)"
    ) +
    coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(
        title = "Trend in annual soil storage change (1951\u20132020)",
        subtitle = "Sen's slope of (Infiltration \u2212 AET \u2212 dSubToUz). Faded = non-significant"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
    )

ggsave(file.path(out_dir, "map_soil_storage_trend.png"),
    p_soil_trend,
    width = 10, height = 8, dpi = 200
)
cat("  -> Saved map_soil_storage_trend.png\n")

# Combined figure: soil + GW trends side by side
fig_storage <- plot_grid(p_soil_trend, p_gw_trend, ncol = 1, align = "v")
ggsave(file.path(out_dir, "map_storage_trends_combined.png"),
    fig_storage,
    width = 10, height = 14, dpi = 200
)
cat("  -> Saved map_storage_trends_combined.png\n")

# =============================================================================
# Maps of MEAN storage residuals (averaged over entire period)
# =============================================================================
cat("Computing mean storage residuals...\n")

# Mean annual soil storage change per catchment
mean_dS_soil <- colMeans(mat_dS_soil, na.rm = TRUE)
# Mean annual GW storage change per catchment
mean_dS_GW <- colMeans(mat_dS_GW, na.rm = TRUE)

storage_summary <- data.table(
    catch_id = common_catches,
    mean_dS_soil_mm = mean_dS_soil,
    mean_dS_GW_mm = mean_dS_GW
)

# Merge with geometry
cats_storage <- catchments[as.character(as.numeric(catchments$catch_id)) %in% common_catches, ]
cats_storage <- merge(cats_storage, storage_summary,
    by.x = "catch_id", by.y = "catch_id", all.x = FALSE
)
cats_storage <- st_transform(cats_storage, 3035)

soil_lim <- quantile(abs(storage_summary$mean_dS_soil_mm), 0.95, na.rm = TRUE)
gw_lim <- quantile(abs(storage_summary$mean_dS_GW_mm), 0.95, na.rm = TRUE)

# --- Map: Mean soil storage residual ------------------------------------------
p_soil_mean <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_storage, aes(fill = mean_dS_soil_mm), color = NA) +
    scale_fill_gradientn(
        colors = palet_gw,
        limits = c(-soil_lim, soil_lim),
        oob = squish,
        name = "\u0394S soil\n(mm/yr)"
    ) +
    coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(
        title = "Mean annual soil storage residual (1951\u20132020)",
        subtitle = "Infiltration \u2212 AET \u2212 dSubToUz. Positive = soil gaining water."
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
    )

# --- Map: Mean GW storage residual --------------------------------------------
p_gw_mean <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_storage, aes(fill = mean_dS_GW_mm), color = NA) +
    scale_fill_gradientn(
        colors = palet_gw,
        limits = c(-gw_lim, gw_lim),
        oob = squish,
        name = "\u0394S GW\n(mm/yr)"
    ) +
    coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(
        title = "Mean annual groundwater storage residual (1951\u20132020)",
        subtitle = "(dSubToUz + prefFlow) \u2212 (quz + qlz + GW_loss). Positive = GW filling."
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
    )

# Save
fig_mean_storage <- plot_grid(p_soil_mean, p_gw_mean, ncol = 1, align = "v")
ggsave(file.path(out_dir, "map_mean_storage_residuals.png"),
    fig_mean_storage,
    width = 10, height = 14, dpi = 200
)
ggsave(file.path(out_dir, "map_mean_soil_storage.png"),
    p_soil_mean,
    width = 10, height = 8, dpi = 200
)
ggsave(file.path(out_dir, "map_mean_GW_storage.png"),
    p_gw_mean,
    width = 10, height = 8, dpi = 200
)
cat("  -> Saved mean storage residual maps\n")

cat("\n--- Mean storage residuals (continental) ---\n")
cat("  Mean dS_soil:", round(mean(mean_dS_soil, na.rm = TRUE), 2), "mm/yr\n")
cat("  Mean dS_GW:  ", round(mean(mean_dS_GW, na.rm = TRUE), 2), "mm/yr\n")

cat("\nAll done!\n")
