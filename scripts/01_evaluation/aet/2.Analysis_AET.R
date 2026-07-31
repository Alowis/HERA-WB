###############################################################################
# AET DIEGO-STYLE - STEP 2: ANALYSIS
#
# Computes per-catchment Spearman rho (with autocorrelation-corrected p-value)
# across four temporal aggregation windows: Daily / 7-Day / 15-Day / Monthly.
#
# Inputs : output/aet_diego/1.homogenized/*.csv
# Outputs: output/aet_diego/2.stats/*.rds + stats_all_windows.csv
###############################################################################

library(data.table)
library(zoo)

# 1. PATHS ---------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
in_dir <- file.path(base_dir, "output", "aet_diego", "1.homogenized")
out_dir <- file.path(base_dir, "output", "aet_diego", "2.stats")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. LOAD HOMOGENIZED MATRICES -------------------------------------
obs_d <- data.table::fread(file.path(in_dir, "gleam_daily_homog.csv"))
mod_d <- data.table::fread(file.path(in_dir, "lisflood_daily_homog.csv"))
obs_m <- data.table::fread(file.path(in_dir, "gleam_monthly_homog.csv"), header = T)
mod_m <- data.table::fread(file.path(in_dir, "lisflood_monthly_homog.csv"), header = T)

obs_d[, date := as.IDate(date)]
mod_d[, date := as.IDate(date)]

daily_cols <- setdiff(names(mod_d), "date")
month_cols <- setdiff(names(mod_m), "date")

# 3. SPEARMAN WITH EFFECTIVE-SAMPLE-SIZE CORRECTION ----------------
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
        zoo::rollmean(v, k, fill = NA, align = "right", na.rm = TRUE)
    }), .SDcols = cols]
    out
}

# 4. RUN EACH WINDOW -----------------------------------------------
message("Daily...")
stats_daily <- compute_spearman(obs_d, mod_d, daily_cols)


# --- 7-day block averages (non-overlapping weeks) ---
obs_d[, week := (seq_len(.N) - 1) %/% 7]
mod_d[, week := (seq_len(.N) - 1) %/% 7]

obs_7d <- obs_d[, lapply(.SD, mean, na.rm = TRUE), by = week]
mod_7d <- mod_d[, lapply(.SD, mean, na.rm = TRUE), by = week]

message("7-day block means: ", nrow(obs_7d), " weeks")
stats_7d <- compute_spearman(obs_7d, mod_7d, daily_cols)

# --- 15-day block averages (non-overlapping biweeks) ---
obs_d[, biweek := (seq_len(.N) - 1) %/% 15]
mod_d[, biweek := (seq_len(.N) - 1) %/% 15]

obs_15d <- obs_d[, lapply(.SD, mean, na.rm = TRUE), by = biweek]
mod_15d <- mod_d[, lapply(.SD, mean, na.rm = TRUE), by = biweek]

message("15-day block means: ", nrow(obs_15d), " biweeks")
stats_15d <- compute_spearman(obs_15d, mod_15d, daily_cols)

message("Evaluating Calendar Month Step...")

obs_d[, month := (seq_len(.N) - 1) %/% 30]
mod_d[, month := (seq_len(.N) - 1) %/% 30]

obs_m <- obs_d[, lapply(.SD, mean, na.rm = TRUE), by = month]
mod_m <- mod_d[, lapply(.SD, mean, na.rm = TRUE), by = month]

message("30-day block means: ", nrow(obs_m), " months")
stats_month <- compute_spearman(obs_m, mod_m, daily_cols)


# message("7-Day rolling mean...")
# obs_7 <- rollmean_dt(obs_d, daily_cols, 7)
# mod_7 <- rollmean_dt(mod_d, daily_cols, 7)
# stats_7d <- compute_spearman(obs_7, mod_7, daily_cols)
#
# message("15-Day rolling mean...")
# obs_15 <- rollmean_dt(obs_d, daily_cols, 15)
# mod_15 <- rollmean_dt(mod_d, daily_cols, 15)
# stats_15d <- compute_spearman(obs_15, mod_15, daily_cols)
#
# message("Monthly...")
# obs_m <- rollmean_dt(obs_d, daily_cols, 30)
# mod_m <- rollmean_dt(mod_d, daily_cols, 30)
# stats_month <- compute_spearman(obs_m, mod_m, month_cols)

# 5. SAVE ----------------------------------------------------------
aet_timeseries <- list(
    daily   = list(obs = obs_d, mod = mod_d),
    `7d`    = list(obs = obs_7d, mod = mod_7d),
    `15d`   = list(obs = obs_15d, mod = mod_15d),
    monthly = list(obs = obs_m, mod = mod_m)
)
saveRDS(aet_timeseries, file.path(out_dir, "AET_time_series_all_scenarios_v2.rds"))
saveRDS(stats_daily, file.path(out_dir, "stats_daily.rds"))
saveRDS(stats_7d, file.path(out_dir, "stats_7d.rds"))
saveRDS(stats_15d, file.path(out_dir, "stats_15d.rds"))
saveRDS(stats_month, file.path(out_dir, "stats_month.rds"))

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

# 5b. REGIME SPEARMAN ----------------------------------------------
# Regime = mean day-of-year cycle (365 values per catchment).
# Compare obs vs mod regime with Spearman rho.
message("Computing regime (mean DOY cycle)...")

RegimeFast <- function(data) {
    names(data) <- c("date", "Q")
    jours <- as.numeric(format(data$date, "%j"))
    ind.j <- tapply(seq(length(jours)), jours, c)
    ind.j <- ind.j[-366] # remove day 366 (leap years)
    Qc <- data.frame(
        date = as.numeric(names(ind.j)),
        mean = sapply(ind.j, function(x) mean(data$Q[x], na.rm = TRUE))
    )
    return(Qc)
}

# Compute regime for each catchment and evaluate Spearman
stats_regime <- rbindlist(lapply(daily_cols, function(id) {
    obs_ts <- data.frame(date = obs_d[["date"]], Q = obs_d[[id]])
    mod_ts <- data.frame(date = mod_d[["date"]], Q = mod_d[[id]])

    obs_regime <- tryCatch(RegimeFast(obs_ts), error = function(e) NULL)
    mod_regime <- tryCatch(RegimeFast(mod_ts), error = function(e) NULL)

    if (is.null(obs_regime) || is.null(mod_regime)) {
        return(list(
            catch_id = id, rho = NA_real_, p_val = NA_real_,
            n_eff = NA_real_, status = "No Data"
        ))
    }

    # Merge on DOY to align
    regime <- merge(obs_regime, mod_regime, by = "date", suffixes = c("_obs", "_mod"))
    valid_idx <- which(!is.na(regime$mean_obs) & !is.na(regime$mean_mod))

    if (length(valid_idx) < 10) {
        return(list(
            catch_id = id, rho = NA_real_, p_val = NA_real_,
            n_eff = NA_real_, status = "No Data"
        ))
    }

    xv <- regime$mean_mod[valid_idx]
    yv <- regime$mean_obs[valid_idx]

    if (sd(xv) == 0 || sd(yv) == 0) {
        return(list(
            catch_id = id, rho = NA_real_, p_val = NA_real_,
            n_eff = NA_real_, status = "No Variance"
        ))
    }

    rho <- cor(xv, yv, method = "spearman")
    n <- length(xv)

    # Autocorrelation correction
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

saveRDS(stats_regime, file.path(out_dir, "stats_regime.rds"))

# Add regime to all_stats
all_stats <- rbindlist(list(
    all_stats,
    cbind(stats_regime, scenario = "Regime")
), use.names = TRUE, fill = TRUE)
data.table::fwrite(all_stats, file.path(out_dir, "stats_all_windows.csv"))

med_regime <- stats_regime[, .(
    median_rho = round(median(rho, na.rm = TRUE), 3),
    pct_ge_0.5 = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 1)
)]
message(
    "Regime: median rho = ", med_regime$median_rho,
    ", % >= 0.5 = ", med_regime$pct_ge_0.5
)

# 6. ANOMALY COMPUTATION -------------------------------------------
# Normalize data to approximate normality using inverse normal transformation
# (quantile mapping to standard normal), then compute anomalies as deviations.
message("Normalizing data and computing anomalies...")

compute_anomalies <- function(dt, cols) {
    out <- copy(dt)
    out[, (cols) := lapply(.SD, function(v) {
        valid <- !is.na(v)
        n_valid <- sum(valid)
        if (n_valid < 3) {
            return(rep(NA_real_, length(v)))
        }
        r <- rep(NA_real_, length(v))
        # Map to uniform via empirical CDF, then to normal via qnorm
        p <- rank(v[valid], ties.method = "average") / (n_valid + 1)
        r[valid] <- qnorm(p)
        r
    }), .SDcols = cols]
    out
}

# Compute anomalies for each temporal window
anom_obs_d <- compute_anomalies(obs_d, daily_cols)
anom_mod_d <- compute_anomalies(mod_d, daily_cols)

anom_obs_7 <- compute_anomalies(obs_7d, daily_cols)
anom_mod_7 <- compute_anomalies(mod_7d, daily_cols)

anom_obs_15 <- compute_anomalies(obs_15d, daily_cols)
anom_mod_15 <- compute_anomalies(mod_15d, daily_cols)

anom_obs_m <- compute_anomalies(obs_m, month_cols)
anom_mod_m <- compute_anomalies(mod_m, month_cols)

plot(anom_mod_d[["233248"]])
# 7. PEARSON ON TOP/BOTTOM 20% OF ANOMALIES ------------------------
# For each catchment, select the time steps where observed anomalies fall
# in the top 20% (high extremes) or bottom 20% (low extremes), then compute
# the Pearson correlation on those subsets.
message("Computing Pearson on top/bottom 20% anomalies...")

compute_pearson_extremes <- function(mat_obs_anom, mat_mod_anom, mat_obs, mat_mod, ids, quantile_lo = 0.2, quantile_hi = 0.8) {
    rbindlist(lapply(ids, function(id) {
        obs_a <- mat_obs_anom[[id]]
        mod_a <- mat_mod_anom[[id]]

        obs_b <- mat_obs[[id]]
        mod_b <- mat_mod[[id]]

        valid_idx <- which(!is.na(obs_a) & !is.na(mod_a))

        if (length(valid_idx) < 10) {
            return(list(
                catch_id = id,
                pearson_top20 = NA_real_,
                pearson_bot20 = NA_real_,
                n_top20 = NA_integer_,
                n_bot20 = NA_integer_,
                status = "No Data"
            ))
        }

        obs_v <- obs_a[valid_idx]
        mod_v <- mod_a[valid_idx]



        # Quantile thresholds based on observed anomalies
        q_lo <- quantile(obs_v, quantile_lo, na.rm = TRUE)
        q_hi <- quantile(obs_v, quantile_hi, na.rm = TRUE)

        # Bottom 20%
        idx_bot <- which(obs_v <= q_lo)
        # Top 20%
        idx_top <- which(obs_v >= q_hi)

        pearson_top <- if (length(idx_top) >= 5 && sd(obs_v[idx_top]) > 0 && sd(mod_v[idx_top]) > 0) {
            cor(obs_v[idx_top], mod_v[idx_top], method = "spearman")
        } else {
            NA_real_
        }

        pearson_bot <- if (length(idx_bot) >= 5 && sd(obs_v[idx_bot]) > 0 && sd(mod_v[idx_bot]) > 0) {
            cor(obs_v[idx_bot], mod_v[idx_bot], method = "spearman")
        } else {
            NA_real_
        }

        list(
            catch_id = id,
            pearson_top20 = round(pearson_top, 4),
            pearson_bot20 = round(pearson_bot, 4),
            n_top20 = length(idx_top),
            n_bot20 = length(idx_bot),
            status = "OK"
        )
    }))
}

pearson_ext_daily <- compute_pearson_extremes(anom_obs_d, anom_mod_d, obs_d, mod_d, daily_cols)
pearson_ext_7d <- compute_pearson_extremes(anom_obs_7, anom_mod_7, obs_7d, mod_7d, daily_cols)
pearson_ext_15d <- compute_pearson_extremes(anom_obs_15, anom_mod_15, obs_15d, mod_15d, daily_cols)
pearson_ext_month <- compute_pearson_extremes(anom_obs_m, anom_mod_m, obs_m, mod_m, month_cols)

# 8. NORMALIZED RMSE ------------------------------------------------
# NRMSE = RMSE / sd(obs), computed per catchment
message("Computing Normalized RMSE...")

compute_nrmse <- function(mat_obs, mat_mod, ids) {
    rbindlist(lapply(ids, function(id) {
        x <- mat_mod[[id]]
        y <- mat_obs[[id]]
        valid_idx <- which(!is.na(x) & !is.na(y))

        if (length(valid_idx) < 10) {
            return(list(catch_id = id, rmse = NA_real_, nrmse = NA_real_, status = "No Data"))
        }

        xv <- x[valid_idx]
        yv <- y[valid_idx]

        rmse_val <- sqrt(mean((xv - yv)^2))
        sd_obs <- sd(yv)

        nrmse_val <- if (sd_obs > 0) rmse_val / sd_obs else NA_real_

        list(
            catch_id = id,
            rmse = round(rmse_val, 4),
            nrmse = round(nrmse_val, 4),
            status = "OK"
        )
    }))
}

nrmse_daily <- compute_nrmse(obs_d, mod_d, daily_cols)
nrmse_7d <- compute_nrmse(obs_7d, mod_7d, daily_cols)
nrmse_15d <- compute_nrmse(obs_15d, mod_15d, daily_cols)
nrmse_month <- compute_nrmse(obs_m, mod_m, month_cols)

# 9. SAVE ANOMALY RESULTS ------------------------------------------
message("Saving anomaly and NRMSE results...")

# Save Pearson extremes
saveRDS(pearson_ext_daily, file.path(out_dir, "pearson_extremes_daily.rds"))
saveRDS(pearson_ext_7d, file.path(out_dir, "pearson_extremes_7d.rds"))
saveRDS(pearson_ext_15d, file.path(out_dir, "pearson_extremes_15d.rds"))
saveRDS(pearson_ext_month, file.path(out_dir, "pearson_extremes_month.rds"))

all_pearson_ext <- rbindlist(list(
    cbind(pearson_ext_daily, scenario = "Daily"),
    cbind(pearson_ext_7d, scenario = "7-Day"),
    cbind(pearson_ext_15d, scenario = "15-Day"),
    cbind(pearson_ext_month, scenario = "Monthly")
), use.names = TRUE, fill = TRUE)
data.table::fwrite(all_pearson_ext, file.path(out_dir, "pearson_extremes_all_windows.csv"))

# Save NRMSE
saveRDS(nrmse_daily, file.path(out_dir, "nrmse_daily.rds"))
saveRDS(nrmse_7d, file.path(out_dir, "nrmse_7d.rds"))
saveRDS(nrmse_15d, file.path(out_dir, "nrmse_15d.rds"))
saveRDS(nrmse_month, file.path(out_dir, "nrmse_month.rds"))

all_nrmse <- rbindlist(list(
    cbind(nrmse_daily, scenario = "Daily"),
    cbind(nrmse_7d, scenario = "7-Day"),
    cbind(nrmse_15d, scenario = "15-Day"),
    cbind(nrmse_month, scenario = "Monthly")
), use.names = TRUE, fill = TRUE)
data.table::fwrite(all_nrmse, file.path(out_dir, "nrmse_all_windows.csv"))

# Summary
med_pearson <- all_pearson_ext[, .(
    median_pearson_top20 = round(median(pearson_top20, na.rm = TRUE), 3),
    median_pearson_bot20 = round(median(pearson_bot20, na.rm = TRUE), 3)
), by = scenario]
print(med_pearson)

med_nrmse <- all_nrmse[, .(
    median_nrmse = round(median(nrmse, na.rm = TRUE), 3)
), by = scenario]
print(med_nrmse)

message("All done. Anomaly stats and NRMSE saved in: ", out_dir)
