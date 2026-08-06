# =============================================================================
# Aggregate land use fractions at the catchment level
#
# For each land use type (forest, irrigated, other, rice, sealed, water),
# extracts the area-weighted mean fraction per catchment per year using
# exact_extract. Saves results as one CSV per land use type.
#
# Output: data/aggregates/landuse/{type}_catchment_yearly.csv
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(terra)
library(sf)
library(exactextractr)
library(data.table)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
landuse_dir <- "D:/tilloal/Documents/06_Floodrivers/landuse/"
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
out_dir <- file.path(base_dir, "data", "aggregates", "landuse")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Land use types -----------------------------------------------------------
lu_types <- c(
    "fracforest", "fracirrigated", "fracother",
    "fracrice", "fracsealed", "fracwater"
)

lu_long_names <- c(
    fracforest = "Forest fraction",
    fracirrigated = "Irrigated cropland fraction",
    fracother = "Other land use fraction",
    fracrice = "Rice paddy fraction",
    fracsealed = "Sealed (impervious) fraction",
    fracwater = "Water body fraction"
)

# --- Load catchments ----------------------------------------------------------
cat("Loading catchments...\n")
catchments <- st_read(gpkg_path, quiet = TRUE)

# Ensure WGS84 for extraction
if (!isTRUE(st_crs(catchments) == st_crs(4326))) {
    catchments <- st_transform(catchments, 4326)
}
catch_ids <- as.character(catchments$catch_id)
n_catch <- length(catch_ids)
cat("  Catchments:", n_catch, "\n")

# --- Process each land use type -----------------------------------------------
for (lu_type in lu_types) {
    cat("\n[", lu_type, "] Processing...\n")

    out_path <- file.path(out_dir, paste0(lu_type, "_catchment_yearly.csv"))

    # Skip if already computed
    if (file.exists(out_path)) {
        cat("  Already exists, skipping.\n")
        next
    }

    # Find all yearly NetCDF files for this type
    pattern <- paste0("^", lu_type, "_European_01min_\\d{4}\\.nc$")
    nc_files <- list.files(landuse_dir, pattern = pattern, full.names = TRUE)

    if (length(nc_files) == 0) {
        cat("  No files found! Skipping.\n")
        next
    }

    # Extract years from filenames
    years <- as.integer(sub(
        paste0(".*", lu_type, "_European_01min_(\\d{4})\\.nc$"),
        "\\1", nc_files
    ))
    ord <- order(years)
    nc_files <- nc_files[ord]
    years <- years[ord]
    cat("  Years:", min(years), "-", max(years), "(", length(years), "files)\n")

    # Result matrix: years x catchments
    result <- matrix(NA_real_, nrow = length(years), ncol = n_catch)

    for (i in seq_along(nc_files)) {
        if (i %% 10 == 1) {
            cat("  Extracting year", years[i], "(", i, "/", length(nc_files), ")\n")
        }

        r <- tryCatch(terra::rast(nc_files[i]), error = function(e) {
            warning("  Could not read: ", basename(nc_files[i]))
            NULL
        })
        if (is.null(r)) next

        # Use first layer if multi-band
        if (nlyr(r) > 1) r <- r[[1]]

        # Extract area-weighted mean per catchment
        means <- exact_extract(r, catchments, "mean")
        result[i, ] <- means
    }

    # Build output data.table
    out_dt <- data.table(year = years)
    for (j in seq_along(catch_ids)) {
        out_dt[, (catch_ids[j]) := result[, j]]
    }

    # Save
    fwrite(out_dt, out_path)
    cat("  Saved:", out_path, "\n")
}

cat("\n--- Summary ---\n")
for (lu_type in lu_types) {
    out_path <- file.path(out_dir, paste0(lu_type, "_catchment_yearly.csv"))
    if (file.exists(out_path)) {
        dt <- fread(out_path, nrows = 1)
        cat(
            "  ", lu_type, ":", ncol(dt) - 1, "catchments,",
            nrow(fread(out_path)), "years\n"
        )
    }
}

cat("\nDone! Land use aggregates saved to:", out_dir, "\n")
