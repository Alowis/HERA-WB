# =============================================================================
# Convert variable-oriented CSVs to Parquet format
#
# Input:
#   data/SHARE/timeseries_var/*.csv  — one CSV per variable (cols = catchments)
#
# Output:
#   data/SHARE/timeseries_var_parquet/*.parquet — same structure, compressed
#
# Parquet advantages:
#   - ~5-10x smaller than CSV (columnar + zstd compression)
#   - Faster to read (no text parsing)
#   - Column-selective reads (load only the catchments you need)
#   - Portable across R, Python, Julia, Spark, etc.
#
# Usage to read back:
#   library(arrow)
#   dt <- read_parquet("file.parquet")
#   # Or specific columns:
#   dt <- read_parquet("file.parquet", col_select = c("time", "12345"))
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(arrow)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
csv_dir <- file.path(base_dir, "data", "SHARE", "timeseries_var")
pq_dir <- file.path(base_dir, "data", "SHARE", "timeseries_var_parquet")
dir.create(pq_dir, showWarnings = FALSE, recursive = TRUE)

# --- Determine time vector ----------------------------------------------------
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) stop("No CSV files found in: ", csv_dir)

n_rows <- nrow(fread(csv_files[1], select = 1L, header = TRUE))
n_days <- as.integer(as.Date("2020-12-31") - as.Date("1951-01-01")) + 1

if (n_rows == n_days) {
    cat("Data is DAILY.\n")
    time_vec <- seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")
} else if (n_rows == n_days * 4) {
    cat("Data is 6-HOURLY.\n")
    time_vec <- seq(
        as.POSIXct("1951-01-01 00:00:00", tz = "UTC"),
        as.POSIXct("2020-12-31 18:00:00", tz = "UTC"),
        by = "6 hours"
    )
} else {
    cat("Unexpected row count:", n_rows, ". Assuming daily.\n")
    time_vec <- seq.Date(as.Date("1951-01-01"), length.out = n_rows, by = "day")
}

# --- Convert each CSV to Parquet ----------------------------------------------
cat(sprintf("\nConverting %d CSV files to Parquet...\n\n", length(csv_files)))

total_csv_mb <- 0
total_pq_mb <- 0

for (csv_path in csv_files) {
    fname <- tools::file_path_sans_ext(basename(csv_path))
    pq_path <- file.path(pq_dir, paste0(fname, ".parquet"))

    # Skip if already converted
    if (file.exists(pq_path)) {
        cat(sprintf("  SKIP (exists): %s\n", fname))
        next
    }

    cat(sprintf("  Reading: %s ...", fname))
    dt <- fread(csv_path, header = TRUE)

    # Prepend time column
    dt[, time := time_vec]
    setcolorder(dt, c("time", setdiff(names(dt), "time")))

    # Write as Parquet with zstd compression
    cat(" writing parquet...")
    write_parquet(dt, pq_path, compression = "zstd", compression_level = 3)

    # Report size reduction
    csv_size <- file.size(csv_path) / 1e6
    pq_size <- file.size(pq_path) / 1e6
    ratio <- csv_size / pq_size
    total_csv_mb <- total_csv_mb + csv_size
    total_pq_mb <- total_pq_mb + pq_size

    cat(sprintf(
        " done. (%.0f MB -> %.0f MB, %.1fx smaller)\n",
        csv_size, pq_size, ratio
    ))

    rm(dt)
    gc()
}

cat(sprintf(
    "\nTotal: %.1f GB -> %.1f GB (%.1fx compression)\n",
    total_csv_mb / 1000, total_pq_mb / 1000,
    total_csv_mb / total_pq_mb
))
cat("Output:", pq_dir, "\n")
