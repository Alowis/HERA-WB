###############################################################################
# AET DIEGO-STYLE - STEP 1: HOMOGENIZE
#
# Aligns GLEAM (observation) and LISFLOOD (model) AET on shared dates and
# shared catchment IDs for daily and monthly resolution.
#
# Inputs (pre-existing):
#   GLEAM daily wide: data/gleam_aet_daily_wide.csv  (from step 0)
#   LISFLOOD daily:   data/aggregates/ActEvapo/ActEvapo_daily_all_years.csv
#   LISFLOOD monthly: data/aggregates/ActEvapo/ActEvapo_monthly_all_years.csv
#
# Outputs: output/aet_diego/1.homogenized/{gleam,lisflood}_{daily,monthly}_homog.csv
###############################################################################

library(data.table)

# 1. PATHS ---------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"

gleam_daily_path <- file.path(base_dir,  "output", "aet_diego", "0.extracted", "gleam_aet_daily_wide.csv")
lf_daily_path <- file.path(
    base_dir, "data", "aggregates", "ActEvapo", "ActEvapo_daily_all_years.csv"
)
lf_month_path <- file.path(
    base_dir, "data", "aggregates", "ActEvapo", "ActEvapo_monthly_all_years.csv"
)

out_dir <- file.path(base_dir, "output", "aet_diego", "1.homogenized")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. HELPERS -------------------------------------------------------
norm_id <- function(x) {
    x <- as.character(x)
    suppressWarnings(ifelse(grepl("^\\d+$", x), as.character(as.integer(x)), x))
}

align_wide <- function(obs, mod) {
    id_obs <- setdiff(names(obs), "date")
    id_mod <- setdiff(names(mod), "date")
    key_obs <- norm_id(id_obs)
    key_mod <- norm_id(id_mod)
    common_keys <- intersect(key_obs, key_mod)
    if (length(common_keys) == 0) stop("No shared catchment ids.")
    obs_sel <- id_obs[match(common_keys, key_obs)]
    mod_sel <- id_mod[match(common_keys, key_mod)]
    common_dates <- intersect(obs$date, mod$date)
    if (length(common_dates) == 0) stop("No shared dates.")
    obs2 <- obs[date %in% common_dates, c("date", obs_sel), with = FALSE]
    mod2 <- mod[date %in% common_dates, c("date", mod_sel), with = FALSE]
    data.table::setnames(obs2, c("date", common_keys))
    data.table::setnames(mod2, c("date", common_keys))
    data.table::setorder(obs2, date)
    data.table::setorder(mod2, date)
    list(obs = obs2, mod = mod2)
}

round_cols <- function(dt, digits = 3) {
    num_cols <- setdiff(names(dt), "date")
    for (col in num_cols) set(dt, j = col, value = round(as.numeric(dt[[col]]), digits))
    dt
}

# 3. DAILY ---------------------------------------------------------
message("Homogenizing DAILY AET...")
gleam_d <- data.table::fread(gleam_daily_path)
lf_d <- data.table::fread(lf_daily_path)
if ("year" %in% names(lf_d)) lf_d[, year := NULL]

daily <- align_wide(gleam_d, lf_d)
daily$obs <- round_cols(daily$obs)
daily$mod <- round_cols(daily$mod)

data.table::fwrite(daily$obs, file.path(out_dir, "gleam_daily_homog.csv"))
data.table::fwrite(daily$mod, file.path(out_dir, "lisflood_daily_homog.csv"))
message("  daily: ", nrow(daily$obs), " dates x ", ncol(daily$obs) - 1, " catchments")

# 4. MONTHLY -------------------------------------------------------
message("Homogenizing MONTHLY AET...")
lf_m <- data.table::fread(lf_month_path)
# LISFLOOD monthly has month_idx, period_start, period_end + id cols
lf_m[, date := substr(as.character(period_end), 1, 7)]
lf_m[, c("month_idx", "period_start", "period_end") := NULL]

# GLEAM monthly: aggregate GLEAM daily wide to monthly means
gleam_m <- copy(gleam_d)
gleam_m[, date := substr(date, 1, 7)]
num_cols <- setdiff(names(gleam_m), "date")
gleam_m <- gleam_m[, lapply(.SD, mean, na.rm = TRUE), by = date, .SDcols = num_cols]

monthly <- align_wide(gleam_m, lf_m)
monthly$obs <- round_cols(monthly$obs)
monthly$mod <- round_cols(monthly$mod)

ptn<-monthly$obs
data.table::fwrite(monthly$obs, file.path(out_dir, "gleam_monthly_homog.csv"))
data.table::fwrite(monthly$mod, file.path(out_dir, "lisflood_monthly_homog.csv"))
message("  monthly: ", nrow(monthly$obs), " dates x ", ncol(monthly$obs) - 1, " catchments")

message("Step 1 done. Homogenized files in: ", out_dir)
