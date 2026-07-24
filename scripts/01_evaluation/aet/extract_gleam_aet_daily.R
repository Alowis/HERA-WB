# Extract GLEAM daily AET and aggregate to catchment level ----------
# Each GLEAM daily file (data/GLEAM/AET_daily/E_YYYY_GLEAM_v4.3a.nc) is a
# yearly stack of ~365 global daily layers. To keep memory low we NEVER
# load the whole stack: we read ONE day at a time and crop it to the
# catchment region before extracting. For every day and every catchment
# we compute area-weighted stats and write one tidy row per catchment/day.

# Library calling --------------------------------------------------
library(sf)
library(terra)
library(exactextractr)

# Path configuration -----------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"

catchments_path <- file.path(
    base_dir, "data", "catchments_analysis_final_v3.gpkg"
)
aet_daily_dir <- file.path(base_dir, "data", "GLEAM", "AET_daily")
out_csv <- file.path(base_dir, "data", "gleam_aet_catchment_daily.csv")

# Catchment identifier column in the gpkg
id_col <- "catch_id"

# Load and align catchments ----------------------------------------
cat("[1/3] Loading catchments...\n")
catchments <- st_read(catchments_path, quiet = TRUE)

# GLEAM grid is WGS84 (EPSG:4326); match it
if (!isTRUE(st_crs(catchments) == st_crs(4326))) {
    catchments <- st_transform(catchments, 4326)
}

# The summary function below refers to v_sf$Id, so expose that name
catchments$Id <- catchments[[id_col]]

# Region of interest: bounding box around the catchment centroids --
# (+/- buffer so edge catchments are fully covered by the crop window)
nco <- sf::st_coordinates(sf::st_centroid(catchments))
buf <- 1 # degrees
crop_win <- terra::ext(
    min(nco[, 1]) - buf, max(nco[, 1]) + buf,
    min(nco[, 2]) - buf, max(nco[, 2]) + buf
)

# Area-weighted stats for ONE daily raster (the requested logic) ---
extract_layer <- function(sm, date, v_sf) {
    stats <- exactextractr::exact_extract(
        x = sm,
        y = v_sf,
        fun = function(df) {
            ok <- !is.na(df$value)
            if (sum(ok) == 0) {
                return(data.frame(
                    n_pixels = nrow(df), n_valid = 0L, n_na = nrow(df),
                    mean_w = NA_real_, sd = NA_real_,
                    p2.5 = NA_real_, p50 = NA_real_, p97.5 = NA_real_
                ))
            }
            v <- df$value[ok]
            w <- df$coverage_fraction[ok]
            w <- w / sum(w)
            mw <- weighted.mean(v, w)
            data.frame(
                n_pixels = nrow(df), n_valid = sum(ok), n_na = sum(!ok),
                mean_w = round(mw, 4),
                sd = round(sqrt(weighted.mean((v - mw)^2, w)), 4),
                p2.5 = round(quantile(v, 0.025, na.rm = TRUE), 4),
                p50 = round(quantile(v, 0.500, na.rm = TRUE), 4),
                p97.5 = round(quantile(v, 0.975, na.rm = TRUE), 4)
            )
        },
        summarize_df = TRUE,
        progress = FALSE
    )
    stats$date <- date
    stats$Id <- v_sf$Id
    # date and Id first, then the statistics
    stats[, c(
        "date", "Id", "n_pixels", "n_valid", "n_na",
        "mean_w", "sd", "p2.5", "p50", "p97.5"
    )]
}

# Process ONE yearly file, one day at a time -> data.frame ---------
process_one_file <- function(file, v_sf, win) {
    # Metadata only (lazy); no pixel data loaded here
    meta <- terra::rast(file)
    n_days <- terra::nlyr(meta)

    # Daily dates from the NetCDF time dimension; fall back to a day
    # sequence built from the year in the filename if it is missing.
    dates <- as.Date(terra::time(meta))
    if (all(is.na(dates))) {
        yr <- as.integer(regmatches(basename(file), regexpr("\\d{4}", basename(file))))
        dates <- seq(as.Date(paste0(yr, "-01-01")), by = "day", length.out = n_days)
    }

    day_list <- vector("list", n_days)
    for (d in seq_len(n_days)) {
        print(d)
        # Read ONLY this day's layer, then crop to the region of interest.
        # crop() reads just the windowed pixels -> minimal memory use.
        sm <- terra::crop(terra::rast(file, lyrs = d), win)
        day_list[[d]] <- extract_layer(sm, dates[d], v_sf)
    }
    do.call(rbind, day_list)
}

# Run over all daily files, writing incrementally ------------------
cat("[2/3] Extracting GLEAM daily AET per catchment...\n")

nc_files <- list.files(
    aet_daily_dir,
    pattern = "^E_\\d{4}_GLEAM_v4\\.3a\\.nc$",
    full.names = TRUE
)
if (length(nc_files) == 0) {
    stop("No GLEAM daily NetCDF files found in: ", aet_daily_dir)
}

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
if (file.exists(out_csv)) file.remove(out_csv)

n_files <- length(nc_files)
for (i in seq_along(nc_files)) {
    f <- nc_files[i]
    cat("  file", i, "of", n_files, "-", basename(f), "\n")

    df <- process_one_file(f, catchments, crop_win)

    # Append to the single CSV (header only for the first file)
    write.table(
        df, out_csv,
        sep = ",", row.names = FALSE,
        col.names = (i == 1), append = (i > 1)
    )
}

# Completion summary -----------------------------------------------
cat("[3/3] Done.\n")
cat(
    "\n=== Completion Summary ===\n",
    "Catchments:        ", nrow(catchments), "\n",
    "GLEAM files:       ", n_files, "\n",
    "Output CSV:        ", out_csv, "\n",
    sep = ""
)
