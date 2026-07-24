###############################################################################
# AET CROSS-COMPARISON - STEP 0: EXTRACT GLEAM DAILY AET BY CATCHMENT
#
# Builds on extract_gleam_aet_daily.R but produces a WIDE-format output
# (date + one column per catchment) compatible with the Diego-style pipeline.
# Processes year by year to keep memory low and saves one CSV per year into
# output/swe_diego_aet/0.extracted/. A final concatenation step merges all
# yearly files into a single daily matrix.
#
# Requires: GLEAM NetCDF files in data/GLEAM/AET_daily/E_YYYY_GLEAM_v4.3a.nc
###############################################################################

library(sf)
library(terra)
library(exactextractr)

# Path configuration -----------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"

catchments_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
aet_daily_dir <- file.path(base_dir, "data", "GLEAM", "AET_daily")
out_dir <- file.path(base_dir, "output", "aet_diego", "0.extracted")
out_merged <- file.path(base_dir, "output", "aet_diego", "0.extracted", "gleam_aet_daily_wide.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

id_col <- "catch_id"

# Load and align catchments ----------------------------------------
cat("[1/3] Loading catchments...\n")
catchments <- st_read(catchments_path, quiet = TRUE)

if (!isTRUE(st_crs(catchments) == st_crs(4326))) {
    catchments <- st_transform(catchments, 4326)
}
catchments$Id <- catchments[[id_col]]

# Crop window from catchment centroids (+1 deg buffer)
nco <- sf::st_coordinates(sf::st_centroid(catchments))
buf <- 1
crop_win <- terra::ext(
    min(nco[, 1]) - buf, max(nco[, 1]) + buf,
    min(nco[, 2]) - buf, max(nco[, 2]) + buf
)

catch_ids <- as.character(catchments$Id)

# Extract one daily layer: returns named numeric vector of catchment means
extract_day <- function(sm, v_sf) {
    exactextractr::exact_extract(sm, v_sf, "mean", progress = FALSE)
}

# Process ONE yearly file -> wide data.frame (date + catchment columns) ----
process_year <- function(file, v_sf, win, ids) {
    meta <- terra::rast(file)
    n_days <- terra::nlyr(meta)

    dates <- as.Date(terra::time(meta))
    if (all(is.na(dates))) {
        yr <- as.integer(regmatches(basename(file), regexpr("\\d{4}", basename(file))))
        dates <- seq(as.Date(paste0(yr, "-01-01")), by = "day", length.out = n_days)
    }

    # Pre-allocate matrix: days x catchments
    mat <- matrix(NA_real_, nrow = n_days, ncol = length(ids))
    colnames(mat) <- ids

    for (d in seq_len(n_days)) {
        if (d %% 30 == 1) cat("    day", d, "/", n_days, "\n")
        sm <- terra::crop(terra::rast(file, lyrs = d), win)
        mat[d, ] <- extract_day(sm, v_sf)
    }

    out <- as.data.frame(mat)
    out$date <- format(dates, "%Y-%m-%d")
    out <- out[, c("date", ids)]
    out
}

# Run over all yearly files, saving one CSV per year ----------------
cat("[2/3] Extracting GLEAM daily AET per catchment (year by year)...\n")

nc_files <- list.files(
    aet_daily_dir,
    pattern = "^E_\\d{4}_GLEAM_v4\\.3a\\.nc$",
    full.names = TRUE
)
nc_files <- sort(nc_files)
if (length(nc_files) == 0) {
    stop("No GLEAM daily NetCDF files found in: ", aet_daily_dir)
}

n_files <- length(nc_files)
yearly_csvs <- character(n_files)

for (i in seq_along(nc_files)) {
    f <- nc_files[i]
    yr <- regmatches(basename(f), regexpr("\\d{4}", basename(f)))
    yr_csv <- file.path(out_dir, paste0("gleam_aet_", yr, ".csv"))
    yearly_csvs[i] <- yr_csv

    if (file.exists(yr_csv)) {
        cat("  [", i, "/", n_files, "] ", basename(f), " -> already exists, skipping\n")
        next
    }

    cat("  [", i, "/", n_files, "] ", basename(f), "\n")
    df <- process_year(f, catchments, crop_win, catch_ids)
    write.csv(df, yr_csv, row.names = FALSE)
}

# Merge all yearly CSVs into a single daily matrix ------------------
# only if all 40 years were computed

lyr=length(yearly_csvs)
print(lyr)
cat("[3/3] Merging yearly CSVs into", basename(out_merged), "...\n")



all_dfs <- lapply(yearly_csvs, function(f) {
    read.csv(f, check.names = FALSE)
})
merged <- do.call(rbind, all_dfs)
merged <- merged[order(merged$date), ]
rownames(merged) <- NULL
write.csv(merged, out_merged, row.names = FALSE)

cat(
    "\n=== Done ===\n",
    "Years processed: ", n_files, "\n",
    "Total days:      ", nrow(merged), "\n",
    "Catchments:      ", ncol(merged) - 1, "\n",
    "Output:          ", out_merged, "\n",
    sep = ""
)
