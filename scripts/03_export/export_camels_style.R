# =============================================================================
# Export HERA data in CAMELS-style CSV format
#
# Input:
#   data/SHARE/timeseries_var/  — one CSV per variable (columns = catchments)
#
# Output:
#   data/SHARE/timeseries_cat/  — one CSV per catchment (columns = variables)
#
# Strategy:
#   Process one variable at a time. For each variable file, read it once and
#   append that column to every catchment's output file. This uses minimal RAM
#   (only one variable in memory at a time) while still reading each large file
#   only once.
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
tss_dir <- file.path(base_dir, "data", "SHARE", "timeseries_var")
out_dir <- file.path(base_dir, "data", "SHARE", "timeseries_cat")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Variable definitions -----------------------------------------------------
variables <- list(
    list(name = "RF", file = "RF_nested_1951_2020.csv"),
    list(name = "SF", file = "SF_nested_1951_2020.csv"),
    list(name = "SNM", file = "SNM_nested_1951_2020.csv"),
    list(name = "INF", file = "INF_nested_1951_2020.csv"),
    list(name = "AET", file = "AET_nested_1951_2020.csv"),
    list(name = "SRF", file = "SRF_nested_1951_2020.csv"),
    list(name = "ULF", file = "ULF_nested_1951_2020.csv"),
    list(name = "SGW", file = "SGW_nested_1951_2020.csv"),
    list(name = "QUZ", file = "QUZ_nested_1951_2020.csv"),
    list(name = "QLZ", file = "QLZ_nested_1951_2020.csv"),
    list(name = "GWL", file = "GWL_nested_1951_2020.csv"),
    list(name = "GWR", file = "GWR_nested_1951_2020.csv"),
    list(name = "SSM", file = "SSM_nested_1951_2020.csv"),
    list(name = "RSM", file = "RSM_nested_1951_2020.csv"),
    list(name = "LSM", file = "LSM_nested_1951_2020.csv"),
    list(name = "SWE", file = "SWE_nested_1951_2020.csv"),
    list(name = "Q", file = "Q_nested_1951_2020.csv")
)

# --- Check which files exist --------------------------------------------------
cat("Checking available files in:", tss_dir, "\n")
available_vars <- list()
for (v in variables) {
    fpath <- file.path(tss_dir, v$file)
    if (file.exists(fpath)) {
        available_vars[[v$name]] <- v
    } else {
        cat("  MISSING:", v$file, "\n")
    }
}
cat("Available:", length(available_vars), "/", length(variables), "variables\n\n")

if (length(available_vars) == 0) stop("No input files found in ", tss_dir)

# --- Determine catchment IDs and temporal resolution --------------------------
first_file <- file.path(tss_dir, available_vars[[1]]$file)
catch_ids <- names(fread(first_file, nrows = 0, header = TRUE))
n_catch <- length(catch_ids)

n_rows <- nrow(fread(first_file, select = 1L, header = TRUE))
n_days <- as.integer(as.Date("2020-12-31") - as.Date("1951-01-01")) + 1

if (n_rows == n_days) {
    cat("Data is DAILY.\n")
    dates <- format(seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"),
        by = "day"
    ), "%Y-%m-%d")
    hours <- rep(0L, n_days)
} else if (n_rows == n_days * 4) {
    cat("Data is 6-HOURLY.\n")
    dates <- format(rep(seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"),
        by = "day"
    ), each = 4), "%Y-%m-%d")
    hours <- rep(c(0L, 6L, 12L, 18L), times = n_days)
} else {
    cat("Unexpected row count:", n_rows, ". Assuming daily.\n")
    dates <- format(seq.Date(as.Date("1951-01-01"),
        length.out = n_rows,
        by = "day"
    ), "%Y-%m-%d")
    hours <- rep(0L, n_rows)
}

cat("Catchments:", n_catch, " | Timesteps:", n_rows, "\n\n")

# =============================================================================
# PASS 1: Write date/hour columns to all catchment files
# =============================================================================
cat("Pass 0: Writing date/hour to all catchment files...\n")
time_dt <- data.table(date = dates, hour = hours)

for (ic in seq_along(catch_ids)) {
    cid <- catch_ids[ic]
    out_path <- file.path(out_dir, paste0(cid, "_timeseries.csv"))
    fwrite(time_dt, out_path)
}
cat("  Done.\n\n")

# =============================================================================
# PASS 2: For each variable, read file once, append column to each catchment
# =============================================================================
var_names <- names(available_vars)

for (iv in seq_along(var_names)) {
    vname <- var_names[iv]
    vinfo <- available_vars[[vname]]
    fpath <- file.path(tss_dir, vinfo$file)

    cat(sprintf(
        "[%d/%d] Reading %s (%s)...\n",
        iv, length(var_names), vname, vinfo$file
    ))

    # Read full file — one variable, all catchments, all timesteps
    dt <- fread(fpath, header = TRUE)
    gc()

    cat(sprintf("  Writing column to %d catchment files...\n", n_catch))

    for (ic in seq_along(catch_ids)) {
        cid <- catch_ids[ic]
        out_path <- file.path(out_dir, paste0(cid, "_timeseries.csv"))

        # Read existing file, add column, overwrite
        ts_dt <- fread(out_path)
        set(ts_dt, j = vname, value = dt[[ic]])
        fwrite(ts_dt, out_path)
    }

    rm(dt)
    gc()
    cat("  Done.\n\n")
}

cat("Export complete!", n_catch, "files in:", out_dir, "\n")
