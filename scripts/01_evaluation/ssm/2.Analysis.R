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
path_out <- "Z:\\ClimateRun4\\nahaUsers\\tilloal\\SoilMositure_2026_06_03\\1.Diego_Merged\\"
file_shp <- "Z:\\ClimateRun4\\nahaUsers\\tilloal\\SoilMositure_2026_06_03\\catchments_analysis_final_v3.gpkg"

# 2. LOAD HOMOGENIZED ARRAYS
mod_d <- fread(paste0(path_out, "lisflood_homogenized.csv"))
obs_d <- fread(paste0(path_out, "esacci_homogenized.csv"))

# Parse exact Date type for time aggregation steps
mod_d[, date := as.IDate(date)]
obs_d[, date := as.IDate(date)]

catchment_cols <- setdiff(names(mod_d), "date")

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

    if (length(valid_idx) < 10 || all(is.na(y))) {
      return(list(
        catch_id = id, rho = NA_real_, p_val = NA_real_,
        n_eff = NA_real_, status = "No Data"
      ))
    }

    xv <- x[valid_idx]
    yv <- y[valid_idx]
    n <- length(xv)

    # Spearman rho (unchanged)
    rho <- cor(xv, yv, method = "spearman")

    # Effective sample size via lag-1 autocorrelation of both series
    r1x <- cor(xv[-n], xv[-1], use = "complete.obs")
    r1y <- cor(yv[-n], yv[-1], use = "complete.obs")
    # Use average lag-1 ACF of the two series
    r1 <- mean(c(r1x, r1y), na.rm = TRUE)
    r1 <- max(0, r1) # floor at 0; negative ACF would inflate n_eff (conservative)

    n_eff <- max(3, floor(n * (1 - r1) / (1 + r1)))

    # t-test with corrected df
    t_stat <- rho * sqrt((n_eff - 2) / (1 - rho^2))
    p_corr <- 2 * pt(abs(t_stat), df = n_eff - 2, lower.tail = FALSE)

    status <- if (is.na(p_corr) || p_corr >= 0.05) "Not Significant" else "Significant"

    list(
      catch_id = id, rho = round(rho, 2), p_val = p_corr,
      n_eff = n_eff, status = status
    )
  }))
} # Effective sample size correction due to possible autocorrelation — it keeps the same <U+03C1> but replaces the naive p-value with a corrected one


compute_spearman <- function(mat_obs, mat_mod, ids) {
  rbindlist(lapply(ids, function(id) {
    x <- mat_mod[[id]]
    y <- mat_obs[[id]]
    valid_idx <- which(!is.na(x) & !is.na(y))
    
    if (length(valid_idx) < 10 || all(is.na(y)) || all(is.na(x))) {
      return(list(
        catch_id = id, rho = NA_real_, p_val = NA_real_,
        n_eff = NA_real_, status = "No Data"
      ))
    }
    
    xv <- x[valid_idx]
    yv <- y[valid_idx]
    n <- length(xv)
    
    if (sd(xv) == 0 || sd(yv) == 0) {
      return(list(
        catch_id = id, rho = NA_real_, p_val = NA_real_,
        n_eff = NA_real_, status = "No Variance"
      ))
    }
    
    rho <- cor(xv, yv, method = "spearman")
    
    r1x <- cor(xv[-n], xv[-1], use = "complete.obs")
    r1y <- cor(yv[-n], yv[-1], use = "complete.obs")
    r1 <- mean(c(r1x, r1y), na.rm = TRUE)
    r1 <- max(0, r1)
    
    n_eff <- max(3, floor(n * (1 - r1) / (1 + r1)))
    t_stat <- rho * sqrt((n_eff - 2) / (1 - rho^2))
    p_corr <- 2 * pt(abs(t_stat), df = n_eff - 2, lower.tail = FALSE)
    
    status <- if (is.na(p_corr) || p_corr >= 0.05) "Not Significant" else "Significant"
    
    list(
      catch_id = id, rho = round(rho, 2), p_val = p_corr,
      n_eff = n_eff, status = status
    )
  }))
}

# -------------------------------------------------------------------------
# SCENARIO A: DAILY DATA CRUNCH
# -------------------------------------------------------------------------
message("Evaluating Daily Data...")
stats_daily <- compute_spearman(mod_d, obs_d, catchment_cols)

# -------------------------------------------------------------------------
# SCENARIO B: 7-DAY ROLLING AVERAGE
# -------------------------------------------------------------------------
message("Evaluating 7-Day Running Window...")

# --- 7-day block averages (non-overlapping weeks) ---
obs_d[, week := (seq_len(.N) - 1) %/% 7]
mod_d[, week := (seq_len(.N) - 1) %/% 7]

obs_7d <- obs_d[, lapply(.SD, mean, na.rm = TRUE), by = week]
mod_7d <- mod_d[, lapply(.SD, mean, na.rm = TRUE), by = week]


stats_7d <- compute_spearman(mod_7d, obs_7d, catchment_cols)

# -------------------------------------------------------------------------
# SCENARIO C: 15-DAY ROLLING AVERAGE
# -------------------------------------------------------------------------
message("Evaluating 15-Day Running Window...")

obs_d[, biweek := (seq_len(.N) - 1) %/% 15]
mod_d[, biweek := (seq_len(.N) - 1) %/% 15]

obs_15d <- obs_d[, lapply(.SD, mean, na.rm = TRUE), by = biweek]
mod_15d <- mod_d[, lapply(.SD, mean, na.rm = TRUE), by = biweek]

message("15-day block means: ", nrow(obs_15d), " biweeks")
stats_15d <- compute_spearman(obs_15d, mod_15d, catchment_cols)

# -------------------------------------------------------------------------
# SCENARIO D: NATURAL CALENDAR MONTH AVERAGE
# -------------------------------------------------------------------------
message("Evaluating Calendar Month Step...")

obs_d[, month := (seq_len(.N) - 1) %/% 30]
mod_d[, month := (seq_len(.N) - 1) %/% 30]

obs_30d <- obs_d[, lapply(.SD, mean, na.rm = TRUE), by = month]
mod_30d <- mod_d[, lapply(.SD, mean, na.rm = TRUE), by = month]

message("15-day block means: ", nrow(obs_15d), " biweeks")
stats_30d <- compute_spearman(obs_30d, mod_30d, catchment_cols)


# -------------------------------------------------------------------------
# Save objects
# -------------------------------------------------------------------------
path_stats <- "Z:\\ClimateRun4\\nahaUsers\\tilloal\\SoilMositure_2026_06_03\\2.Diego_Analysis\\0.Stats_time_windows\\"
sm_timeseries <- list(
  daily   = list(mod = mod_d, obs = obs_d),
  `7d`    = list(mod = mod_7d, obs = obs_7d),
  `15d`   = list(mod = mod_15d, obs = obs_15d),
  monthly = list(mod = mod_30d, obs = obs_30d)
)
saveRDS(sm_timeseries, file.path(path_stats, "SM_time_series_all_scenarios.rds"))
saveRDS(stats_daily, file.path(path_stats, "stats_daily.rds"))
saveRDS(stats_7d, file.path(path_stats, "stats_7d.rds"))
saveRDS(stats_15d, file.path(path_stats, "stats_15d.rds"))
saveRDS(stats_month, file.path(path_stats, "stats_month.rds"))



# =============================================================================
# SNOW-MASKED ANALYSIS: Remove snow-covered days from SM comparison
# =============================================================================
message("Applying snow mask to soil moisture analysis...")

base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
start_year <- 1991

# Load daily SWE from LISFLOOD to identify snow-covered days
swe_dir <- file.path(base_dir, "output", "swe_diego", "1.homogenized")
swe_daily_path <- file.path(swe_dir, "lisflood_daily_homog.csv")

if (file.exists(swe_daily_path)) {
  swe_daily <- fread(swe_daily_path, header = TRUE)
  swe_daily[, date := as.IDate(date)]
  swe_daily[, year := year(date)]

  # Find common catchments and dates between SWE and SM
  swe_cols <- setdiff(names(swe_daily), "date")
  sm_cols_clean <- sub("^X", "", catchment_cols)
  common_mask_cols <- intersect(swe_cols, sm_cols_clean)

  common_dates <- intersect(swe_daily$date, obs_d$date)
  cat("  Common catchments for masking:", length(common_mask_cols), "\n")
  cat("  Common dates:", length(common_dates), "\n")

  # Create masked copies
  mod_masked <- copy(mod_d[date %in% common_dates])
  obs_masked <- copy(obs_d[date %in% common_dates])
  obs_masked[, year := year(date)]
  mod_masked[, year := year(date)]
  
  swe_sub <- swe_daily[date %in% common_dates]
  
  obs_masked <- obs_masked[year >= start_year, ]
  mod_masked <- mod_masked[year >= start_year, ]
  swe_sub <- swe_sub[year >= start_year, ]


  # Apply mask: set SM to NA where SWE > 0 (per catchment)
  for (col in common_mask_cols) {
    # Find matching column names (may have X prefix in SM data)
    sm_col <- if (col %in% names(mod_masked)) col else paste0("X", col)
    if (sm_col %in% names(mod_masked) && col %in% names(swe_sub)) {
      snow_mask <- swe_sub[[col]] > 5
      mod_masked[snow_mask, (sm_col) := NA_real_]
      obs_masked[snow_mask, (sm_col) := NA_real_]
    }
  }

  # Recompute stats on snow-masked data
  masked_cols <- intersect(catchment_cols, names(mod_masked))

  message("  Recomputing daily stats (snow-masked)...")
  max(obs_masked$date)
  max(mod_masked$date)
  plot(mod_masked[["X1029"]])
  points(obs_masked[["X1029"]],col=2)
  stats_daily_masked <- compute_spearman(mod_masked, obs_masked, masked_cols)

  median(stats_daily$rho,na.rm=T)
  
  message("  Recomputing 7-day stats (snow-masked)...")
  obs_masked[, week := (seq_len(.N) - 1) %/% 7]
  mod_masked[, week := (seq_len(.N) - 1) %/% 7]
  
  obs_7dm <- obs_masked[, lapply(.SD, mean, na.rm = TRUE), by = week]
  mod_7dm <- mod_masked[, lapply(.SD, mean, na.rm = TRUE), by = week]
  
  stats_7d_masked <- compute_spearman(mod_7dm, obs_7dm, masked_cols)

  median (stats_7d$rho,na.rm=T)
  median (stats_7d_masked$rho,na.rm=T)
  plot(stats_7d_masked$rho,stats_7d$rho)
  abline(a=0,b=1)
  
  message("  Recomputing 15-day stats (snow-masked)...")
  obs_masked[, biweek := (seq_len(.N) - 1) %/% 15]
  mod_masked[, biweek := (seq_len(.N) - 1) %/% 15]
  
  obs_15dm <- obs_masked[, lapply(.SD, mean, na.rm = TRUE), by = biweek]
  mod_15dm <- mod_masked[, lapply(.SD, mean, na.rm = TRUE), by = biweek]
  
  stats_15d_masked <- compute_spearman(mod_15dm, obs_15dm, masked_cols)
  
  
  message("  Recomputing monthly stats (snow-masked)...")
  obs_masked[, month := (seq_len(.N) - 1) %/% 30]
  mod_masked[, month := (seq_len(.N) - 1) %/% 30]
  
  obs_30dm <- obs_masked[, lapply(.SD, mean, na.rm = TRUE), by = month]
  mod_30dm <- mod_masked[, lapply(.SD, mean, na.rm = TRUE), by = month]
  
  stats_30d_masked <- compute_spearman(mod_30dm, obs_30dm, masked_cols)
  
  
  # Save snow-masked stats
  sm_timeseries_masked <- list(
    daily   = list(mod = mod_masked, obs = obs_masked),
    `7d`    = list(mod = mod_7dm, obs = obs_7dm),
    `15d`   = list(mod = mod_15dm, obs = obs_15dm),
    monthly = list(mod = mod_30dm, obs = obs_30dm)
  )
  saveRDS(sm_timeseries_masked, file.path(path_stats, "SM_time_series_masked_all_scenarios.rds"))
  
  saveRDS(stats_daily_masked, file.path(path_stats, "stats_daily_snowmasked.rds"))
  saveRDS(stats_7d_masked, file.path(path_stats, "stats_7d_snowmasked.rds"))
  saveRDS(stats_15d_masked, file.path(path_stats, "stats_15d_snowmasked.rds"))
  saveRDS(stats_30d_masked, file.path(path_stats, "stats_30d_snowmasked.rds"))

  cat("  Median rho (daily, all):", median(stats_daily$rho, na.rm = TRUE), "\n")
  cat("  Median rho (daily, snow-masked):", median(stats_daily_masked$rho, na.rm = TRUE), "\n")
  cat("  Median rho (monthly, all):", median(stats_month$rho, na.rm = TRUE), "\n")
  cat("  Median rho (monthly, snow-masked):", median(stats_30d_masked$rho, na.rm = TRUE), "\n")

  message("Snow-masked SM analysis complete.")
} else {
  message("  SWE daily file not found at: ", swe_daily_path)
  message("  Skipping snow-masked analysis. Run SWE homogenization first.")
}
