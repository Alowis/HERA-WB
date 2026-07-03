# Optimized daily GlobSnow SWE extraction per catchment ----------------
# Reads daily GlobSnow NetCDF files, masks flags (water/-1, mountain/-2),
# and extracts area-weighted mean SWE per catchment.
#
# Optimizations vs. the monthly per-file loop:
#   1. Catchments reprojected to the raster CRS ONCE (EASE-Grid, EPSG:3408)
#   2. Rasters cropped to the European catchment extent before extraction
#   3. Files stacked per year -> a SINGLE exact_extract call per year
#      (polygon/grid intersection weights are computed once, reused for
#       every layer, instead of recomputed for every file)
#   4. Optional parallelism across years via foreach %dopar%

# Library calling --------------------------------------------------
library(sf)
library(terra)
library(exactextractr)

# Set to TRUE to parallelize across years (requires foreach + doParallel)
use_parallel <- TRUE
n_workers <- 4

# Path configuration -----------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"

catchments_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
globsnow_dir <- file.path(base_dir, "data", "globsnow_swe_daily")
globsnow_csv_out <- file.path(base_dir, "data", "globsnow_swe_catchment_daily.csv")

# Daily GlobSnow filename pattern. ADJUST if your files differ.
# Assumes the date is the leading 8 digits: YYYYMMDD_...nc
file_pattern <- "^\\d{8}.*\\.nc$"

# Load catchments --------------------------------------------------
cat("[1/4] Loading catchments...\n")
catchments <- st_read(catchments_path, quiet = TRUE)

# List daily files -------------------------------------------------
nc_files <- list.files(globsnow_dir, pattern = file_pattern, full.names = TRUE)
if (length(nc_files) == 0) {
  stop("No daily GlobSnow NetCDF files found in: ", globsnow_dir)
}
n_files <- length(nc_files)
cat("[2/4] Found", n_files, "daily files.\n")

# Parse dates (leading YYYYMMDD) and group by year
date_tokens <- substr(basename(nc_files), 1, 8)
file_dates <- as.Date(date_tokens, format = "%Y%m%d")
if (any(is.na(file_dates))) {
  stop(
    "Could not parse YYYYMMDD dates from ", sum(is.na(file_dates)),
    " filename(s). Adjust 'file_pattern' / date parsing for your naming."
  )
}
file_years <- format(file_dates, "%Y")
files_by_year <- split(
  data.frame(path = nc_files, date = file_dates, stringsAsFactors = FALSE),
  file_years
)

# One-time spatial setup: reproject catchments to the raster CRS ---
crs_ease <- terra::crs(terra::rast(nc_files[1]))
catch_ease <- sf::st_transform(catchments, crs_ease)
crop_ext <- terra::ext(terra::vect(catch_ease)) * 1.05 # small buffer
# Plain numeric window (SpatExtent can't be sent to parallel workers)
crop_vec <- as.vector(crop_ext)

# Per-year extraction function -------------------------------------
extract_year <- function(year_df, catch_ease, crop_vec) {
  year_df <- year_df[order(year_df$date), ]
  crop_ext <- terra::ext(crop_vec) # rebuild SpatExtent inside worker

  # Read each file once, select the SWE layer, crop early to Europe
  layers <- lapply(year_df$path, function(f) {
    r <- tryCatch(terra::rast(f), error = function(e) NULL)
    if (is.null(r)) {
      return(NULL)
    }
    if (terra::nlyr(r) > 1) r <- r[[1]] # SWE is the first field
    terra::crop(r, crop_ext)
  })

  keep <- !vapply(layers, is.null, logical(1))
  layers <- layers[keep]
  dates_k <- year_df$date[keep]
  if (length(layers) == 0) {
    return(NULL)
  }

  # Stack and mask flags (negative = water/mountain/no-data)
  r_stack <- terra::rast(layers)
  r_stack[r_stack < 0] <- NA

  # SINGLE exact_extract call for the whole year's stack
  means <- exactextractr::exact_extract(
    r_stack, catch_ease, "mean",
    progress = FALSE
  )
  # means: rows = catchments, cols = days -> transpose to days x catchments
  out <- as.data.frame(t(as.matrix(means)))
  colnames(out) <- catch_ease$catch_id
  out <- cbind(date = format(dates_k, "%Y-%m-%d"), out)
  rownames(out) <- NULL
  out
}

# Run extraction ---------------------------------------------------
cat(
  "[3/4] Extracting SWE by year",
  if (use_parallel) "(parallel)" else "(sequential)", "...\n"
)

if (use_parallel) {
  library(foreach)
  library(doParallel)
  cl <- parallel::makeCluster(n_workers)
  doParallel::registerDoParallel(cl)
  results <- foreach::foreach(
    yr = files_by_year,
    .packages = c("terra", "sf", "exactextractr")
  ) %dopar% {
    extract_year(yr, catch_ease, crop_vec)
  }
  parallel::stopCluster(cl)
} else {
  results <- lapply(
    files_by_year, extract_year,
    catch_ease = catch_ease, crop_vec = crop_vec
  )
}

results <- results[!vapply(results, is.null, logical(1))]
globsnow_df <- do.call(rbind, results)
globsnow_df <- globsnow_df[order(globsnow_df$date), ]
rownames(globsnow_df) <- NULL

# Sanity check: western catchments (lon < 0) must have data --------
west_ids <- catchments$catch_id[
  sf::st_coordinates(sf::st_centroid(sf::st_transform(catchments, 4326)))[, 1] < 0
]
west_cols <- intersect(as.character(west_ids), names(globsnow_df))
if (length(west_cols) > 0) {
  west_vals <- as.matrix(globsnow_df[, west_cols, drop = FALSE])
  n_west_with_data <- sum(colSums(!is.na(west_vals)) > 0)
  cat(
    "[check] Western catchments (lon < 0): ", length(west_cols),
    " | with >=1 non-NA SWE value: ", n_west_with_data, "\n",
    sep = ""
  )
  if (n_west_with_data == 0) {
    stop(
      "No SWE data west of 0 deg longitude - likely a CRS mismatch between ",
      "the GlobSnow raster (EASE-Grid) and catchments."
    )
  }
}

# Save daily CSV ---------------------------------------------------
cat("[4/4] Saving daily GlobSnow catchment CSV...\n")
dir.create(dirname(globsnow_csv_out), recursive = TRUE, showWarnings = FALSE)
write.csv(globsnow_df, globsnow_csv_out, row.names = FALSE)

cat(
  "\n=== Done ===\n",
  "Days extracted:  ", nrow(globsnow_df), "\n",
  "Catchments:      ", ncol(globsnow_df) - 1, "\n",
  "Output:          ", globsnow_csv_out, "\n",
  sep = ""
)
