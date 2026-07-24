###############################################################################
# CROSS-COMPARISON ANALYSIS: LISFLOOD vs ESA CCI SOIL MOISTURE
#
# This script:
# 1. Creates one scatter panel for each temporal aggregation
# 2. Uses HEX bins instead of millions of points
# 3. Colours catchments by Area
# 4. Computes catchment-level performance statistics
# 5. Produces diagnostic plots (rho vs Area)
###############################################################################

library(data.table)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(patchwork)
library(zoo)
library(data.table)
library(scales)
library(viridis)


# 1. DEFINE PATHS
path_out  <- "Z:\\ClimateRun4\\nahaUsers\\tilloal\\SoilMositure_2026_06_03\\1.Diego_Merged\\"
file_shp  <- "Z:\\ClimateRun4\\nahaUsers\\tilloal\\SoilMositure_2026_06_03\\catchments_analysis_final_v3.gpkg"

# 2. LOAD HOMOGENIZED ARRAYS
lisf_mat <- fread(paste0(path_out, "lisflood_homogenized.csv"))
esa_mat  <- fread(paste0(path_out, "esacci_homogenized.csv"))

# Parse exact Date type for time aggregation steps
lisf_mat[, date := as.IDate(date)]
esa_mat[, date := as.IDate(date)]

catchment_cols <- setdiff(names(lisf_mat), "date")

# 3. STATS CORE EVALUATOR FUNCTION
compute_spearman <- function(mat_lisf, mat_esa, ids) {
  rbindlist(lapply(ids, function(id) {
    x <- mat_lisf[[id]]
    y <- mat_esa[[id]]
    
    # Complete cases filter to drop paired NAs dynamically per catchment
    valid_idx <- which(!is.na(x) & !is.na(y))
    
    # Condition 1: Entirely missing data -> White
    if (length(valid_idx) < 3 || all(is.na(y))) {
      return(list(catch_id = id, rho = NA_real_, p_val = NA_real_, status = "No Data"))
    }
    
    res <- cor.test(x[valid_idx], y[valid_idx], method = "spearman", exact = FALSE)
    
    # Condition 2: Evaluate significance status -> Insignificant becomes Grey
    status <- if (res$p.value >= 0.05) "Not Significant" else "Significant"
    
    list(catch_id = id, rho = round(res$estimate, 2), p_val = res$p.value, status = status)
  }))
}
compute_spearman <- function(mat_lisf, mat_esa, ids) {
  rbindlist(lapply(ids, function(id) {
    x <- mat_lisf[[id]]
    y <- mat_esa[[id]]
    valid_idx <- which(!is.na(x) & !is.na(y))
    
    if (length(valid_idx) < 10 || all(is.na(y)))
      return(list(catch_id = id, rho = NA_real_, p_val = NA_real_,
                  n_eff = NA_real_, status = "No Data"))
    
    xv <- x[valid_idx]
    yv <- y[valid_idx]
    n  <- length(xv)
    
    # Spearman rho (unchanged)
    rho <- cor(xv, yv, method = "spearman")
    
    # Effective sample size via lag-1 autocorrelation of both series
    r1x <- cor(xv[-n], xv[-1], use = "complete.obs")
    r1y <- cor(yv[-n], yv[-1], use = "complete.obs")
    # Use average lag-1 ACF of the two series
    r1  <- mean(c(r1x, r1y), na.rm = TRUE)
    r1  <- max(0, r1)   # floor at 0; negative ACF would inflate n_eff (conservative)
    
    n_eff <- max(3, floor(n * (1 - r1) / (1 + r1)))
    
    # t-test with corrected df
    t_stat <- rho * sqrt((n_eff - 2) / (1 - rho^2))
    p_corr <- 2 * pt(abs(t_stat), df = n_eff - 2, lower.tail = FALSE)
    
    status <- if (is.na(p_corr) || p_corr >= 0.05) "Not Significant" else "Significant"
    
    list(catch_id = id, rho = round(rho, 2), p_val = p_corr,
         n_eff = n_eff, status = status)
  }))
} # Effective sample size correction due to possible autocorrelation — it keeps the same ρ but replaces the naive p-value with a corrected one

# -------------------------------------------------------------------------
# SCENARIO A: DAILY DATA CRUNCH
# -------------------------------------------------------------------------
message("Evaluating Daily Data...")
stats_daily <- compute_spearman(lisf_mat, esa_mat, catchment_cols)

# -------------------------------------------------------------------------
# SCENARIO B: 7-DAY ROLLING AVERAGE
# -------------------------------------------------------------------------
message("Evaluating 7-Day Running Window...")
lisf_7d <- copy(lisf_mat)[, (catchment_cols) := lapply(.SD, function(v) rollmean(v, 7, fill = NA, align = "right", na.rm = TRUE)), .SDcols = catchment_cols]
esa_7d  <- copy(esa_mat)[, (catchment_cols) := lapply(.SD, function(v) rollmean(v, 7, fill = NA, align = "right", na.rm = TRUE)), .SDcols = catchment_cols]
stats_7d <- compute_spearman(lisf_7d, esa_7d, catchment_cols)

# -------------------------------------------------------------------------
# SCENARIO C: 15-DAY ROLLING AVERAGE
# -------------------------------------------------------------------------
message("Evaluating 15-Day Running Window...")
lisf_15d <- copy(lisf_mat)[, (catchment_cols) := lapply(.SD, function(v) rollmean(v, 15, fill = NA, align = "right", na.rm = TRUE)), .SDcols = catchment_cols]
esa_15d  <- copy(esa_mat)[, (catchment_cols) := lapply(.SD, function(v) rollmean(v, 15, fill = NA, align = "right", na.rm = TRUE)), .SDcols = catchment_cols]
stats_15d <- compute_spearman(lisf_15d, esa_15d, catchment_cols)

# -------------------------------------------------------------------------
# SCENARIO D: NATURAL CALENDAR MONTH AVERAGE
# -------------------------------------------------------------------------
message("Evaluating Calendar Month Step...")
lisf_month <- copy(lisf_mat)[, year_month := format(date, "%Y-%m")][, lapply(.SD, mean, na.rm = TRUE), by = year_month, .SDcols = catchment_cols]
esa_month  <- copy(esa_mat)[, year_month := format(date, "%Y-%m")][, lapply(.SD, mean, na.rm = TRUE), by = year_month, .SDcols = catchment_cols]
stats_month <- compute_spearman(lisf_month, esa_month, catchment_cols)


# -------------------------------------------------------------------------
# Save objects
# -------------------------------------------------------------------------
path_stats <- "Z:\\ClimateRun4\\nahaUsers\\tilloal\\SoilMositure_2026_06_03\\2.Diego_Analysis\\0.Stats_time_windows\\"
sm_timeseries <- list(
  daily   = list(lisf = lisf_mat,   esa = esa_mat),
  `7d`    = list(lisf = lisf_7d,    esa = esa_7d),
  `15d`   = list(lisf = lisf_15d,   esa = esa_15d),
  monthly = list(lisf = lisf_month, esa = esa_month)
)
saveRDS(sm_timeseries, file.path(path_stats, "SM_time_series_all_scenarios.rds"))
saveRDS(stats_daily, file.path(path_stats, "stats_daily.rds"))
saveRDS(stats_7d, file.path(path_stats, "stats_7d.rds"))
saveRDS(stats_15d, file.path(path_stats, "stats_15d.rds"))
saveRDS(stats_month, file.path(path_stats, "stats_month.rds"))

