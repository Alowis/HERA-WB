# =============================================================================
# Multi-variable evolution plot per catchment
# Style: multiple smoothed curves with confidence bands on the same panel,
#        one panel per catchment, each variable has its own y-axis color.
#
# Variables:
#   1. Total annual snowmelt (green) — from snowmelt aggregates
#   2. Maximum monthly soil moisture (blue) — from surface_soil_moisture aggregates
#   3. Maximum 7-day direct runoff (purple) — from runoff 7-day aggregates
#   4. Mean spring AET (orange) — from ActEvapo aggregates (Mar-May)
#
# Each variable is standardized (z-score) so they can share the same y-axis.
# The LOESS smoother + CI band gives the visual style from the reference figure.
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(ggplot2)
library(sf)
library(dplyr)
library(lubridate)
library(scales)
library(cowplot)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
agg_dir <- file.path(base_dir, "data", "aggregates")
tss_dir <- file.path(base_dir, "data", "tss_postprocess")
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
out_dir <- file.path(base_dir, "output", "temporal_evolution")
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
    if (catch_id %in% col_names) {
        return(catch_id)
    }
    x_catch <- paste0("X", catch_id)
    if (x_catch %in% col_names) {
        return(x_catch)
    }
    stop("Catchment ", catch_id, " not found.")
}

# --- Load all datasets once ---------------------------------------------------
cat("Loading datasets...\n")

# 1. Snowmelt (monthly aggregates or from TSS preprocessing)
melt_path <- file.path(agg_dir, "snowmelt", "snowmelt_monthly_all_years.csv")
if (!file.exists(melt_path)) {
    # Fallback: use snow_water_equivalent monthly as proxy
    melt_path <- file.path(
        agg_dir, "snow_water_equivalent",
        "snow_water_equivalent_monthly_all_years.csv"
    )
    cat("  Snowmelt aggregate not found, using SWE as fallback.\n")
}
melt_dt <- fread(melt_path)
melt_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
melt_dt[, year := year(date)]
melt_dt[, month := month(date)]

# 2. Surface soil moisture (monthly)
sm_path <- file.path(
    agg_dir, "surface_soil_moisture",
    "surface_soil_moisture_monthly_all_years.csv"
)
sm_dt <- fread(sm_path)
sm_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
sm_dt[, year := year(date)]
sm_dt[, month := month(date)]

# 3. Direct runoff (7-day)
runoff7_path <- file.path(agg_dir, "runoff", "runoff_7day_all_years.csv")
runoff7_dt <- fread(runoff7_path)
runoff7_dt[, date := seq.Date(as.Date("1951-01-01"), by = "7 days", length.out = .N)]
runoff7_dt[, year := year(date)]

# 4. AET (monthly)
aet_path <- file.path(agg_dir, "ActEvapo", "ActEvapo_monthly_all_years.csv")
aet_dt <- fread(aet_path)
aet_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
aet_dt[, year := year(date)]
aet_dt[, month := month(date)]

# --- Compute annual indicators per catchment ----------------------------------
compute_annual_indicators <- function(catch_id) {
    cat("  Processing catchment:", catch_id, "\n")

    # Meta columns to exclude
    meta_melt <- c("month_idx", "period_start", "period_end", "date", "year", "month", "window")
    meta_sm <- meta_melt
    meta_r7 <- c("window", "period_start", "period_end", "date", "year")
    meta_aet <- meta_melt

    # 1. Total annual snowmelt (sum all months)
    melt_col <- find_col(catch_id, setdiff(names(melt_dt), meta_melt))
    annual_melt <- melt_dt[, .(snowmelt = sum(get(melt_col), na.rm = TRUE)),
        by = .(year)
    ]

    # 2. Maximum monthly soil moisture per year
    sm_col <- find_col(catch_id, setdiff(names(sm_dt), meta_sm))
    annual_sm <- sm_dt[, .(max_sm = max(get(sm_col), na.rm = TRUE)),
        by = .(year)
    ]

    # 3. Maximum 7-day direct runoff per year
    r7_col <- find_col(catch_id, setdiff(names(runoff7_dt), meta_r7))
    annual_runoff <- runoff7_dt[, .(max_runoff7 = max(get(r7_col), na.rm = TRUE)),
        by = .(year)
    ]

    # 4. Mean spring AET (Mar–May)
    aet_col <- find_col(catch_id, setdiff(names(aet_dt), meta_aet))
    annual_aet <- aet_dt[month %in% 3:5, .(spring_aet = mean(get(aet_col), na.rm = TRUE)),
        by = .(year)
    ]

    # Merge all
    out <- merge(annual_melt, annual_sm, by = "year", all = TRUE)
    out <- merge(out, annual_runoff, by = "year", all = TRUE)
    out <- merge(out, annual_aet, by = "year", all = TRUE)

    return(out)
}

# --- Standardize (z-score) for overlay plotting -------------------------------
standardize <- function(x) {
    mu <- mean(x, na.rm = TRUE)
    sigma <- sd(x, na.rm = TRUE)
    if (is.na(sigma) || sigma == 0) {
        return(rep(0, length(x)))
    }
    (x - mu) / sigma
}

# --- Build panel plot for one catchment ---------------------------------------
make_panel <- function(annual_dt, panel_label, catch_name) {
    # Standardize each variable
    dt <- copy(annual_dt)
    dt[, snowmelt_z := standardize(snowmelt)]
    dt[, max_sm_z := standardize(max_sm)]
    dt[, max_runoff7_z := standardize(max_runoff7)]
    dt[, spring_aet_z := standardize(spring_aet)]

    # Melt to long format
    dt_long <- melt(annual_dt,
        id.vars = "year",
        measure.vars = c("snowmelt", "max_sm", "max_runoff7", "spring_aet"),
        variable.name = "variable", value.name = "value"
    )

    # Labels and colors matching reference style
    var_labels <- c(
        "snowmelt" = "Total snowmelt",
        "max_sm" = "Max monthly SM",
        "max_runoff7" = "Max 7-d runoff",
        "spring_aet" = "Mean spring AET"
    )
    var_colors <- c(
        "snowmelt" = "purple", # green
        "max_sm" = "darkgreen", # blue
        "max_runoff7" = "blue", # purple
        "spring_aet" = "#d67e00" # orange
    )

    dt_long[, variable := factor(variable, levels = names(var_labels))]

    ggplot(dt_long, aes(x = year, y = value, color = variable, fill = variable)) +
        geom_smooth(
            method = "loess", span = 0.35, se = TRUE,
            alpha = 0.15, linewidth = 0.9
        ) +
        geom_line(linewidth = 0.3, alpha = 0.5) +
        scale_color_manual(values = var_colors, labels = var_labels, name = NULL) +
        scale_fill_manual(values = var_colors, labels = var_labels, name = NULL) +
        labs(
            title = paste0(panel_label, "  ", catch_name),
            x = NULL, y = "Standardized anomaly"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", size = 12),
            legend.position = "none",
            panel.grid.minor = element_blank()
        )
}

# --- Compute and plot for all catchments --------------------------------------
cat("Computing annual indicators...\n")

panels <- list()
labels <- c("a", "b", "c", "d")

for (i in seq_along(catchments_to_plot)) {
    catch <- catchments_to_plot[[i]]
    annual <- compute_annual_indicators(catch$id)
    panels[[i]] <- make_panel(annual, labels[i], catch$name)
}

# Create a shared legend from one panel
p_legend_src <- panels[[1]] +
    theme(
        legend.position = "bottom",
        legend.text = element_text(size = 10)
    )
legend <- get_legend(p_legend_src)

# Compose figure
grid <- plot_grid(plotlist = panels, ncol = 2, align = "hv")
fig <- plot_grid(grid, legend, ncol = 1, rel_heights = c(1, 0.08))

ggsave(file.path(out_dir, "Figure_multivar_catchments.png"),
    fig,
    width = 14, height = 10, dpi = 200
)
cat("  -> Saved Figure_multivar_catchments.png\n")

cat("\nDone!\n")
