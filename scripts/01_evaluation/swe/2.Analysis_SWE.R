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

library(data.table)
library(zoo)

# 1. PATHS ---------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
in_dir <- file.path(base_dir, "output", "swe_diego", "1.homogenized")
out_dir <- file.path(base_dir, "output", "swe_diego", "2.stats")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. LOAD HOMOGENIZED MATRICES -------------------------------------
obs_d <- data.table::fread(file.path(in_dir, "globsnow_daily_homog.csv"))
mod_d <- data.table::fread(file.path(in_dir, "lisflood_daily_homog.csv"),header = T)
obs_m <- data.table::fread(file.path(in_dir, "globsnow_monthly_homog.csv"),header=T)
mod_m <- data.table::fread(file.path(in_dir, "lisflood_monthly_homog.csv"),header=T)

obs_d[, date := as.IDate(date)]
mod_d[, date := as.IDate(date)]

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

# 4. RUN EACH WINDOW ----------------------------------------------
message("Daily...")
stats_daily <- compute_spearman(obs_d, mod_d, daily_cols)


message("7-Day rolling mean...")
obs_7 <- rollmean_dt(obs_d, daily_cols, 7)
mod_7 <- rollmean_dt(mod_d, daily_cols, 7)
stats_7d <- compute_spearman(obs_7, mod_7, daily_cols)

message("15-Day rolling mean...")
obs_15 <- rollmean_dt(obs_d, daily_cols, 15)
mod_15 <- rollmean_dt(mod_d, daily_cols, 15)
stats_15d <- compute_spearman(obs_15, mod_15, daily_cols)



d <- as.Date(obs_15$date)

plot(d[100:700], obs_15[["303662"]][100:700], type = "l", col = "#1b9e77",
     xlab = "Date", ylab = "SWE (mm)", main = "Catchment 303662 - 15-day mean")
lines(d, mod_15[["303662"]], col = "#d95f02")
legend("topright", c("GlobSnow", "LISFLOOD"), col = c("#1b9e77", "#d95f02"), lty = 1)

cor(obs_m[["303662"]],mod_m[["303662"]],use="complete.obs")


obs_d[["date"]]
mod_d[["date"]]
plot(obs_d[["86849"]],mod_d[["86849"]])

day<-as.Date(obs_d$date)
plot(day[100:700], obs_d[["86849"]][100:700], type = "l", col = "#1b9e77",
     xlab = "Date", ylab = "SWE (mm)", main = "Catchment 303662 - 15-day mean")
lines(day, mod_d[["86849"]], col = "#d95f02")
lines(d, obs_15[["86849"]], col = "blue")
lines(d, mod_15[["86849"]], col = "red")

legend("topright", c("GlobSnow", "LISFLOOD"), col = c("#1b9e77", "#d95f02"), lty = 1)






message("Monthly (reused monthly outputs)...")
# Determine overlapping period as set intersection
lf_months=mod_m$date
glowsnow_months=obs_m$date
overlap_dates <- intersect(lf_months, glowsnow_months)
# Subset both data frames to overlapping dates
put=mod_m$date %in% overlap_dates
unique(mod_m$date)
obs_m <- obs_m[obs_m$date %in% overlap_dates, ]
mod_m <- mod_m[mod_m$date %in% overlap_dates, ]

stats_month <- compute_spearman(obs_m, mod_m, month_cols)

# 5. SAVE ----------------------------------------------------------
swe_timeseries <- list(
    daily   = list(obs = obs_d, mod = mod_d),
    `7d`    = list(obs = obs_7, mod = mod_7),
    `15d`   = list(obs = obs_15, mod = mod_15),
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

med <- all_stats[, .(median_rho = round(median(rho, na.rm = TRUE), 3),
    pct_ge_0.5 = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 1)),
by = scenario
]
print(med)
message("Step 2 done. Stats in: ", out_dir)
