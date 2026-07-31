###############################################################################
# SWE CROSS-COMPARISON - STEP 2: ANALYSIS
#
# Mirrors Diego_script/2.Analysis.R but for SWE (GlobSnow vs LISFLOOD).
# Computes per-catchment Spearman rho across four temporal windows:
#   Daily / 7-Day / 15-Day (rolling means from daily) / Monthly (reused output)
# Uses an effective-sample-size correction for the p-value (autocorrelation).
#
# Inputs : output/swe_diego/1.homogenized/*.csv
# Outputs: output/swe_diego/2.stats/*.rds  and  stats_all_windows.csv
###############################################################################

# --- Load shared workspace ----------------------------------------------------
source("R/load_workspace.R")
library(zoo)
library(data.table)
library(terra)
library(sf)
library(exactextractr)
# 1. PATHS ---------------------------------------------------------
in_dir <- file.path(base_dir, "output", "swe_diego", "1.homogenized")
out_dir <- file.path(base_dir, "output", "swe_diego", "2.stats")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. LOAD HOMOGENIZED MATRICES -------------------------------------
obs_d <- fread(file.path(in_dir, "globsnow_daily_homog.csv"), header = TRUE)
mod_d <- fread(file.path(in_dir, "lisflood_daily_homog.csv"), header = TRUE)
obs_m <- fread(file.path(in_dir, "globsnow_monthly_homog.csv"), header = TRUE)
mod_m <- fread(file.path(in_dir, "lisflood_monthly_homog.csv"), header = TRUE)

obs_d[, date := as.IDate(date)]
mod_d[, date := as.IDate(date)]

plot(mod_d$date)
daily_cols <- setdiff(names(mod_d), "date")
month_cols <- setdiff(names(mod_m), "date")

# 3. SPEARMAN WITH EFFECTIVE-SAMPLE-SIZE CORRECTION ----------------
# (obs = GlobSnow, mod = LISFLOOD). Keeps rho; corrects the p-value for
# temporal autocorrelation via lag-1 ACF (Diego's approach).
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

rollmean_dt <- function(dt, cols, k) {
    out <- copy(dt)
    out[, (cols) := lapply(.SD, function(v) {
        zoo::rollmean(v, k, fill = NA, align = "center", na.rm = TRUE)
    }), .SDcols = cols]
    out
}
# 
# rollmean_dt2 <- function(dt, cols, k) {
#   out <- copy(dt)
#   out[, (cols) := lapply(.SD, function(v) {
#     data.table::frollmean(v, k, align = "center", na.rm = TRUE)
#   }), .SDcols = cols]
#   out
# }




norm_id <- function(x) {
  x <- as.character(x)
  suppressWarnings(ifelse(grepl("^\\d+$", x), as.character(as.integer(x)), x))
}
# 4. RUN EACH WINDOW ----------------------------------------------
message("Daily...")
stats_daily <- compute_spearman(obs_d, mod_d, daily_cols)


message("7-Day rolling mean...")
# obs_7 <- rollmean_dt(obs_d, daily_cols, 7)
# mod_7 <- rollmean_dt(mod_d, daily_cols, 7)
# stats_7d <- compute_spearman(obs_7, mod_7, daily_cols)
# 
# mean(stats_7d$rho,na.rm=T)

# obs_72 <- rollmean_dt2(obs_d, daily_cols, 7)
# mod_72 <- rollmean_dt2(mod_d, daily_cols, 7)
# stats_7d2 <- compute_spearman(obs_72, mod_72, daily_cols)
# 
# median(stats_7d2$rho,na.rm=T)
# 
# obs_72 <- rollmean_fast(obs_d, daily_cols, 7)
# mod_72 <- rollmean_fast(mod_d, daily_cols, 7)
# stats_7d2 <- compute_spearman(obs_72, mod_72, daily_cols)
# 
# median(stats_7d2$rho,na.rm=T)
# 
# message("15-Day rolling mean...")
# obs_15 <- rollmean_dt2(obs_d, daily_cols, 15)
# mod_15 <- rollmean_dt2(mod_d, daily_cols, 15)
# stats_15d <- compute_spearman(obs_15, mod_15, daily_cols)
# 
# median(stats_15d$rho,na.rm=T)

# --- 7-day block averages (non-overlapping weeks) ---
obs_d[, week := (seq_len(.N) - 1) %/% 7]
mod_d[, week := (seq_len(.N) - 1) %/% 7]

obs_7d <- obs_d[, lapply(.SD, mean, na.rm = TRUE), by = week]
mod_7d <- mod_d[, lapply(.SD, mean, na.rm = TRUE), by = week]

message("7-day block means: ", nrow(obs_7d), " weeks")
stats_7d <- compute_spearman(obs_7d, mod_7d, daily_cols)

# --- 14-day block averages (non-overlapping biweeks) ---
obs_d[, biweek := (seq_len(.N) - 1) %/% 15]
mod_d[, biweek := (seq_len(.N) - 1) %/% 15]

obs_15d <- obs_d[, lapply(.SD, mean, na.rm = TRUE), by = biweek]
mod_15d <- mod_d[, lapply(.SD, mean, na.rm = TRUE), by = biweek]

message("15-day block means: ", nrow(obs_15d), " biweeks")
stats_15d <- compute_spearman(obs_15d, mod_15d, daily_cols)


d <- as.Date(obs_15d$date)

plot(d[100:700], obs_15d[["303662"]][100:700],
    type = "l", col = "#1b9e77",
    xlab = "Date", ylab = "SWE (mm)", main = "Catchment 303662 - 15-day mean"
)
lines(d, mod_15d[["303662"]], col = "#d95f02")
legend("topright", c("GlobSnow", "LISFLOOD"), col = c("#1b9e77", "#d95f02"), lty = 1)

cor(obs_m[["303662"]], mod_m[["303662"]], use = "complete.obs")


obs_d[["date"]]
mod_d[["date"]]
plot(obs_d[["86849"]], mod_d[["86849"]])
abline(a=0,b=1)
day <- as.Date(obs_d$date)
plot(day[100:700], obs_d[["86849"]][100:700],
    type = "l", col = "#1b9e77",
    xlab = "Date", ylab = "SWE (mm)", main = "Catchment 303662 - 15-day mean"
)
lines(day, mod_d[["86849"]], col = "#d95f02")
lines(d, obs_15d[["86849"]], col = "blue")
lines(d, mod_15d[["86849"]], col = "red")

legend("topright", c("GlobSnow", "LISFLOOD"), col = c("#1b9e77", "#d95f02"), lty = 1)






message("Monthly (reused monthly outputs)...")
# Determine overlapping period as set intersection
lf_months <- mod_m$date
glowsnow_months <- obs_m$date
overlap_dates <- intersect(lf_months, glowsnow_months)
# Subset both data frames to overlapping dates
put <- mod_m$date %in% overlap_dates
unique(mod_m$date)
obs_m <- obs_m[obs_m$date %in% overlap_dates, ]
mod_m <- mod_m[mod_m$date %in% overlap_dates, ]

stats_month <- compute_spearman(obs_m, mod_m, month_cols)

# 5. SAVE ----------------------------------------------------------
swe_timeseries <- list(
    daily   = list(obs = obs_d, mod = mod_d),
    `7d`    = list(obs = obs_7d, mod = mod_7d),
    `15d`   = list(obs = obs_15d, mod = mod_15d),
    monthly = list(obs = obs_m, mod = mod_m)
)
saveRDS(swe_timeseries, file.path(out_dir, "SWE_time_series_all_scenarios.rds"))
saveRDS(stats_daily, file.path(out_dir, "stats_daily.rds"))
saveRDS(stats_7d, file.path(out_dir, "stats_7d.rds"))
saveRDS(stats_15d, file.path(out_dir, "stats_15d.rds"))
saveRDS(stats_month, file.path(out_dir, "stats_month.rds"))

# Combined tidy CSV for quick inspection
all_stats <- rbindlist(list(
    cbind(stats_daily, scenario = "Daily"),
    cbind(stats_7d, scenario = "7-Day"),
    cbind(stats_15d, scenario = "15-Day"),
    cbind(stats_month, scenario = "Monthly")
), use.names = TRUE, fill = TRUE)
data.table::fwrite(all_stats, file.path(out_dir, "stats_all_windows.csv"))

med <- all_stats[, .(
    median_rho = round(median(rho, na.rm = TRUE), 3),
    pct_ge_0.5 = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 1)
),
by = scenario
]
print(med)
message("Step 2 done. Stats in: ", out_dir)


# =============================================================================
# SNOW-FILTERED ANALYSIS: Keep only catchments with >= 5% snow days (SWE > 5mm)
# =============================================================================
message("Applying snow catchment filter...")

mod_d[, year := year(date)]

# Identify catchments with meaningful snow cover
n_days_total <- nrow(mod_d)
snow_pct <- mod_d[, lapply(.SD, function(x) sum(x > 5, na.rm = TRUE) / n_days_total),
    .SDcols = daily_cols
]

snow_pct2 <- obs_d[, lapply(.SD, function(x) sum(x > 5, na.rm = TRUE) / n_days_total),
                  .SDcols = daily_cols
]

catch_cols_swe <- setdiff(names(mod_d), c("date", "year"))

annual_max_swe <- mod_d[, lapply(.SD, max, na.rm = TRUE),
                          by = year, .SDcols = catch_cols_swe
]

years_no_snow <- colSums(annual_max_swe[, ..catch_cols_swe] <= 5, na.rm = TRUE)
# Keep only catchments with snow every year
snow_catches <- names(years_no_snow[years_no_snow <= 20])
snow_catches <- names(snow_pct)[as.numeric(snow_pct[1, ]) > 0.01]

cat(
    "  Catchments with >= 1% snow days (SWE > 5mm):", length(snow_catches),
    "out of", length(daily_cols), "\n"
)


#filter catchments based on standard deviation of elevation

elev_std_path <- file.path(base_dir, "data", "elvstd_European_01min.nc") # adjust path if needed
r_elev_std <- rast(elev_std_path)

dem_path <- file.path(base_dir, "data", "dem.nc") # adjust path if needed
r_dem <- rast(dem_path)

# Ensure catchments are in WGS84 for extraction
shp <- catchments
shp_wgs <- if (st_crs(shp) == st_crs(4326)) shp else st_transform(shp, 4326)

# Area-weighted mean of elevation std within each catchment
shp_wgs$elev_std <- exact_extract(r_elev_std, shp_wgs, "mean")
shp_wgs$elev_std2 <- exact_extract(r_dem, shp_wgs, "stdev")

# Filter: elevation std < 50
flat_catches <- norm_id(shp_wgs$catch_id[shp_wgs$elev_std2 < 200])


catches_iceland <- ws$catch_ids[-match(ws$iceland_ids,ws$catch_ids)]



# 3. Intersect both criteria
snow_flat <- intersect(norm_id(sub("^X", "", snow_catches)), flat_catches)
snow_flat <- intersect(catches_iceland, snow_flat)

# Subset to snow-relevant catchments
obs_d_snow <- obs_d[, c("date", snow_flat), with = FALSE]
mod_d_snow <- mod_d[, c("date", snow_flat), with = FALSE]


# Recompute stats on filtered catchments
message("  Recomputing daily stats (snow catchments only)...")
stats_daily_snow <- compute_spearman(obs_d_snow, mod_d_snow, snow_flat)


# --- 7-day block averages (non-overlapping weeks) ---

message("7-day block means: ", nrow(obs_7d), " weeks")
stats_7d_snow <- compute_spearman(obs_7d, mod_7d, snow_flat)

# --- 14-day block averages (non-overlapping biweeks) ---)

message("  Recomputing monthly stats...")
snow_month_cols <- intersect(snow_flat, month_cols)
obs_m_snow <- obs_m[, c("date", snow_month_cols), with = FALSE]
mod_m_snow <- mod_m[, c("date", snow_month_cols), with = FALSE]
stats_month_snow <- compute_spearman(obs_m_snow, mod_m_snow, snow_month_cols)

# Save snow-filtered stats
saveRDS(stats_daily_snow, file.path(out_dir, "stats_daily_snow_filtered.rds"))
saveRDS(stats_7d_snow, file.path(out_dir, "stats_7d_snow_filtered.rds"))
saveRDS(stats_15d_snow, file.path(out_dir, "stats_15d_snow_filtered.rds"))
saveRDS(stats_month_snow, file.path(out_dir, "stats_month_snow_filtered.rds"))

saveRDS(snow_flat, file.path(out_dir, "snow_filter.rds"))

cat("  Median rho (monthly, all):", median(stats_month$rho, na.rm = TRUE), "\n")
cat("  Median rho (monthly, snow-filtered):", median(stats_month_snow$rho, na.rm = TRUE), "\n")


# Combined tidy CSV for quick inspection
all_stats_snow <- rbindlist(list(
  cbind(stats_daily_snow, scenario = "Daily"),
  cbind(stats_7d_snow, scenario = "7-Day"),
  cbind(stats_15d_snow, scenario = "15-Day"),
  cbind(stats_month_snow, scenario = "Monthly")
), use.names = TRUE, fill = TRUE)
data.table::fwrite(all_stats_snow, file.path(out_dir, "stats_all_windows_snowcat.csv"))


med <- all_stats_snow[, .(
  median_rho = round(median(rho, na.rm = TRUE), 3),
  pct_ge_0.5 = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 1)
),
by = scenario
]
med
message("Snow-filtered SWE analysis complete.")
