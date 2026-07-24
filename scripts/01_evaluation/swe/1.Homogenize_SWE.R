###############################################################################
# SWE CROSS-COMPARISON - STEP 1: HOMOGENIZE
#
# Mirrors Diego_script/1.Homogenize.R but for SWE (GlobSnow vs LISFLOOD).
# Reuses already-saved daily and monthly outputs (no re-extraction):
#   Observation (GlobSnow):
#     data/globsnow_swe_catchment_daily.csv   (date=YYYY-MM-DD + id cols)
#     data/globsnow_swe_catchment.csv         (date=YYYY-MM    + id cols)
#   Model (LISFLOOD):
#     data/aggregates/snow_water_equivalent/snow_water_equivalent_daily_all_years.csv
#         (year, date=YYYY-MM-DD + id cols)
#     data/aggregates/snow_water_equivalent/snow_water_equivalent_monthly_all_years.csv
#         (month_idx, period_start, period_end + id cols)
#
# Produces aligned wide matrices (same dates, same catchments) for daily and
# monthly resolution, ready for the analysis step.
###############################################################################

library(data.table)

# 1. PATHS ---------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"

gs_daily_path <- file.path(base_dir, "data", "globsnow_swe_catchment_daily.csv")
gs_month_path <- file.path(base_dir, "data", "globsnow_swe_catchment.csv")
lf_daily_path <- file.path(
    base_dir, "data", "aggregates", "snow_water_equivalent",
    "snow_water_equivalent_daily_all_years.csv"
)
lf_month_path <- file.path(
    base_dir, "data", "aggregates", "snow_water_equivalent",
    "snow_water_equivalent_monthly_all_years.csv"
)

out_dir <- file.path(base_dir, "output", "swe_diego", "1.homogenized")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 2. HELPERS -------------------------------------------------------
# Normalise catchment ids for matching (strip leading zeros)
norm_id <- function(x) {
    x <- as.character(x)
    suppressWarnings(ifelse(grepl("^\\d+$", x), as.character(as.integer(x)), x))
}

# Align two wide tables (obs, mod) on shared `date` values and shared
# catchment-id columns; returns the two aligned data.tables.
align_wide <- function(obs, mod) {
    id_obs <- setdiff(names(obs), "date")
    id_mod <- setdiff(names(mod), "date")

    # Map original -> normalised id, then find the common set
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

    # Rename both to the normalised keys so columns line up exactly
    data.table::setnames(obs2, c("date", common_keys))
    data.table::setnames(mod2, c("date", common_keys))
    data.table::setorder(obs2, date)
    data.table::setorder(mod2, date)
    list(obs = obs2, mod = mod2)
}

round_cols <- function(dt, digits = 3) {
    num_cols <- setdiff(names(dt), "date")
    for (col in num_cols) {
        set(dt, j = col, value = round(as.numeric(dt[[col]]), digits))
    }
    dt
}

# 3. DAILY ---------------------------------------------------------
message("Homogenizing DAILY SWE...")
gs_d <- data.table::fread(gs_daily_path)
lf_d <- data.table::fread(lf_daily_path)
if ("year" %in% names(lf_d)) lf_d[, year := NULL] # drop redundant year col

daily <- align_wide(gs_d, lf_d)
daily$obs <- round_cols(daily$obs)
daily$mod <- round_cols(daily$mod)

data.table::fwrite(daily$obs, file.path(out_dir, "globsnow_daily_homog.csv"))
data.table::fwrite(daily$mod, file.path(out_dir, "lisflood_daily_homog.csv"))
message(
    "  daily: ", nrow(daily$obs), " dates x ",
    ncol(daily$obs) - 1, " catchments"
)

# 4. MONTHLY -------------------------------------------------------
message("Homogenizing MONTHLY SWE...")
gs_m <- data.table::fread(gs_month_path,header = T) # date = YYYY-MM
lf_m <- data.table::fread(lf_month_path) # month_idx, period_start, period_end

# LISFLOOD monthly date -> YYYY-MM from period_start; drop index columns
lf_m[, date := substr(as.character(period_end), 1, 7)]

lisflood_dates <- format(
  seq.Date(
    as.Date("1951-01-01"),
    as.Date("2020-12-01"),
    by = "month"
  ),
  "%Y-%m"
)

lf_m[, c("month_idx", "period_start", "period_end") := NULL]

monthly <- align_wide(gs_m, lf_m)
monthly$obs <- round_cols(monthly$obs)
monthly$mod <- round_cols(monthly$mod)

data.table::fwrite(monthly$obs, file.path(out_dir, "globsnow_monthly_homog.csv"))
data.table::fwrite(monthly$mod, file.path(out_dir, "lisflood_monthly_homog.csv"))
message(
    "  monthly: ", nrow(monthly$obs), " dates x ",
    ncol(monthly$obs) - 1, " catchments"
)

message("Step 1 done. Homogenized files in: ", out_dir)

