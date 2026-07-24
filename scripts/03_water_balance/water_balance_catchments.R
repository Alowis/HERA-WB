# =============================================================================
# Water balance comparison per catchment
# Inputs: precipitation (rain + snowmelt)
# Outputs: preferential flow, actual evapotranspiration, direct runoff,
#          infiltration, recharge from soil to groundwater (percolation UZ->LZ)
#
# Plot style: stacked area for outputs, line for total input,
#             one panel per catchment
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(ggplot2)
library(lubridate)
library(cowplot)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
agg_dir  <- file.path(base_dir, "data", "aggregates")
tss_dir  <- file.path(base_dir, "data", "tss_postprocess")
out_dir  <- file.path(base_dir, "output", "temporal_evolution")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Catchment choices (modify these) -----------------------------------------
catchments_to_plot <- list(
    list(id = "307920", name = "Avignon"),
    list(id = "325150", name = "Aquila"),
    list(id = "191110", name = "Brussels"),
    list(id = "289805", name = "Aosta")
)

# --- Helper: find column in header --------------------------------------------
find_col <- function(catch_id, col_names) {
    if (catch_id %in% col_names) return(catch_id)
    x_catch <- paste0("X", catch_id)
    if (x_catch %in% col_names) return(x_catch)
    stop("Catchment ", catch_id, " not found.")
}

# --- Helper: load monthly aggregate or TSS fallback ---------------------------
load_monthly <- function(var_name, tss_filename, agg_method = "sum") {
    agg_path <- file.path(agg_dir, var_name, paste0(var_name, "_monthly_all_years.csv"))
    if (file.exists(agg_path)) {
        dt <- fread(agg_path)
        dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
        dt[, year := year(date)]
        return(dt)
    }
    # Fallback: read TSS
    cat("  Aggregate not found for", var_name, "— reading TSS...\n")
    tss_path <- file.path(tss_dir, tss_filename)
    dt <- fread(tss_path, header = TRUE)
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
        monthly <- dt[, lapply(.SD, sum, na.rm = TRUE), by = .(year, month), .SDcols = catch_cols]
    } else {
        monthly <- dt[, lapply(.SD, mean, na.rm = TRUE), by = .(year, month), .SDcols = catch_cols]
    }
    monthly[, date := as.Date(paste(year, month, "15", sep = "-"))]
    return(monthly)
}

# --- Load all datasets --------------------------------------------------------
cat("Loading datasets...\n")

# Inputs
rain_dt     <- load_monthly("rainfall", "rainUpsX_nested_1951_2020.csv", "sum")
snowmelt_dt <- load_monthly("snowmelt", "snowMeltUpsX_nested_1951_2020.csv", "sum")

# Outputs
prefflow_dt <- load_monthly("prefflow", "prefFlowUpsX_nested_1951_2020.csv", "sum")
aet_dt      <- load_monthly("ActEvapo", "ActEvapo_nested_1951_2020.csv", "sum")
runoff_dt   <- load_monthly("runoff", "surfaceRunoffUpsX_nested_1951_2020.csv", "sum")
infil_dt    <- load_monthly("infiltration", "infUpsX_nested_1951_2020.csv", "sum")
perc_dt     <- load_monthly("percolation", "percUZLZUpsX_nested_1951_2020.csv", "sum")

# --- Compute annual totals per catchment and plot -----------------------------
make_water_balance_panel <- function(catch_id, catch_name, panel_label) {
    cat("  Processing:", catch_name, "(", catch_id, ")\n")

    meta <- c("month_idx", "period_start", "period_end", "date", "year", "month", "window")

    # Get column name for this catchment in each dataset
    rain_col <- find_col(catch_id, setdiff(names(rain_dt), meta))
    melt_col <- find_col(catch_id, setdiff(names(snowmelt_dt), meta))
    pref_col <- find_col(catch_id, setdiff(names(prefflow_dt), meta))
    aet_col  <- find_col(catch_id, setdiff(names(aet_dt), meta))
    ro_col   <- find_col(catch_id, setdiff(names(runoff_dt), meta))
    inf_col  <- find_col(catch_id, setdiff(names(infil_dt), meta))
    perc_col <- find_col(catch_id, setdiff(names(perc_dt), meta))

    # Annual sums
    ann_rain <- rain_dt[, .(rain = sum(get(rain_col), na.rm = TRUE)), by = year]
    ann_melt <- snowmelt_dt[, .(snowmelt = sum(get(melt_col), na.rm = TRUE)), by = year]
    ann_pref <- prefflow_dt[, .(pref_flow = sum(get(pref_col), na.rm = TRUE)), by = year]
    ann_aet  <- aet_dt[, .(aet = sum(get(aet_col), na.rm = TRUE)), by = year]
    ann_ro   <- runoff_dt[, .(runoff = sum(get(ro_col), na.rm = TRUE)), by = year]
    ann_inf  <- infil_dt[, .(infiltration = sum(get(inf_col), na.rm = TRUE)), by = year]
    ann_perc <- perc_dt[, .(percolation = sum(get(perc_col), na.rm = TRUE)), by = year]

    # Merge
    wb <- Reduce(function(a, b) merge(a, b, by = "year", all = TRUE),
        list(ann_rain, ann_melt, ann_pref, ann_aet, ann_ro, ann_inf, ann_perc))

    # Total input
    wb[, total_input := rain + snowmelt]

    # Melt outputs for stacked area
    out_long <- melt(wb,
        id.vars = "year",
        measure.vars = c("aet", "infiltration", "runoff", "pref_flow", "percolation"),
        variable.name = "component", value.name = "value"
    )

    # Factor order (bottom to top in stacked area)
    out_long[, component := factor(component,
        levels = c("aet", "infiltration", "percolation", "pref_flow", "runoff"),
        labels = c("Actual ET", "Infiltration", "Percolation UZ\u2192LZ",
                   "Preferential flow", "Direct runoff")
    )]

    # Colors
    comp_colors <- c(
        "Actual ET" = "#d62728",
        "Infiltration" = "#2ca02c",
        "Percolation UZ\u2192LZ" = "#1f77b4",
        "Preferential flow" = "#9467bd",
        "Direct runoff" = "#ff7f0e"
    )

    # Plot: stacked area for outputs + line for inputs
    p <- ggplot() +
        geom_area(data = out_long, aes(x = year, y = value, fill = component),
                  alpha = 0.7, position = "stack") +
        geom_line(data = wb, aes(x = year, y = total_input),
                  color = "black", linewidth = 0.9, linetype = "solid") +
        geom_line(data = wb, aes(x = year, y = rain),
                  color = "grey40", linewidth = 0.6, linetype = "dashed") +
        geom_line(data = wb, aes(x = year, y = snowmelt),
                  color = "steelblue", linewidth = 0.6, linetype = "dotted") +
        scale_fill_manual(values = comp_colors, name = "Outputs") +
        labs(
            title = paste0(panel_label, ") ", catch_name, " (", catch_id, ")"),
            x = NULL, y = "Annual total (mm/year)"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", size = 12),
            legend.position = "none",
            panel.grid.minor = element_blank()
        )

    return(p)
}

# --- Generate panels ----------------------------------------------------------
cat("Building panels...\n")

panels <- list()
labels <- c("a", "b", "c", "d")

for (i in seq_along(catchments_to_plot)) {
    catch <- catchments_to_plot[[i]]
    panels[[i]] <- make_water_balance_panel(catch$id, catch$name, labels[i])
}

# Shared legend from first panel
p_for_legend <- panels[[1]] +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 10)) +
    # Add input lines to legend
    geom_line(aes(x = 0, y = 0, linetype = "Total input (P + snowmelt)"), color = "black") +
    geom_line(aes(x = 0, y = 0, linetype = "Rainfall"), color = "grey40") +
    geom_line(aes(x = 0, y = 0, linetype = "Snowmelt"), color = "steelblue") +
    scale_linetype_manual(values = c(
        "Total input (P + snowmelt)" = "solid",
        "Rainfall" = "dashed",
        "Snowmelt" = "dotted"
    ), name = "Inputs")

legend <- get_legend(p_for_legend)

# Compose
grid <- plot_grid(plotlist = panels, ncol = 2, align = "hv")

# Add a manual legend annotation
# Simpler approach: just add legend text below
legend_text <- ggdraw() +
    draw_label("Lines: black solid = total input | grey dashed = rainfall | blue dotted = snowmelt",
               size = 10, x = 0.5)

fig <- plot_grid(grid, legend, legend_text,
    ncol = 1, rel_heights = c(1, 0.08, 0.03))

ggsave(file.path(out_dir, "Figure_water_balance_catchments.png"),
    fig,
    width = 14, height = 10, dpi = 200
)
cat("  -> Saved Figure_water_balance_catchments.png\n")

cat("\nDone!\n")
