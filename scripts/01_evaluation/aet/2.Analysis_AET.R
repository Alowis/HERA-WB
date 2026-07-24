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
obs_m <- data.table::fread(file.path(in_dir, "gleam_monthly_homog.csv"),header=T)
mod_m <- data.table::fread(file.path(in_dir, "lisflood_monthly_homog.csv"),header=T)

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

message("7-Day rolling mean...")
obs_7 <- rollmean_dt(obs_d, daily_cols, 7)
mod_7 <- rollmean_dt(mod_d, daily_cols, 7)
stats_7d <- compute_spearman(obs_7, mod_7, daily_cols)

message("15-Day rolling mean...")
obs_15 <- rollmean_dt(obs_d, daily_cols, 15)
mod_15 <- rollmean_dt(mod_d, daily_cols, 15)
stats_15d <- compute_spearman(obs_15, mod_15, daily_cols)

message("Monthly...")
obs_m <- rollmean_dt(obs_d, daily_cols, 30)
mod_m <- rollmean_dt(mod_d, daily_cols, 30)
stats_month <- compute_spearman(obs_m, mod_m, month_cols)

# 5. SAVE ----------------------------------------------------------
aet_timeseries <- list(
    daily   = list(obs = obs_d, mod = mod_d),
    `7d`    = list(obs = obs_7, mod = mod_7),
    `15d`   = list(obs = obs_15, mod = mod_15),
    monthly = list(obs = obs_m, mod = mod_m)
)
saveRDS(aet_timeseries, file.path(out_dir, "AET_time_series_all_scenarios.rds"))
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
