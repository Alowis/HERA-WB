# =============================================================================
# Export land use fractions in CAMELS-style CSV format
#
# One CSV per catchment with yearly land use fractions (1950–2020).
# Uses pre-computed aggregates from aggregate_landuse_catchments.R if available,
# otherwise extracts directly from NetCDF rasters.
#
# Output:
#   data/HERA_CAMELS/landuse/
#     {catch_id}_landuse.csv          — one per catchment
#     landuse_variable_description.csv — metadata
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(sf)
library(terra)
library(exactextractr)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
landuse_dir <- "D:/tilloal/Documents/06_Floodrivers/landuse/"
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
agg_lu_dir <- file.path(base_dir, "data", "aggregates", "landuse")
out_dir <- file.path(base_dir, "data", "HERA_CAMELS", "landuse")
out_attr <- file.path(base_dir, "data", "HERA_CAMELS", "attributes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_attr, showWarnings = FALSE, recursive = TRUE)

# --- Land use type definitions ------------------------------------------------
lu_types <- list(
    list(
        name = "fracforest", col_name = "forest",
        long_name = "Forest fraction", units = "-",
        description = "Fraction of catchment covered by forest"
    ),
    list(
        name = "fracirrigated", col_name = "irrigated",
        long_name = "Irrigated cropland fraction", units = "-",
        description = "Fraction of catchment with irrigated agriculture"
    ),
    list(
        name = "fracother", col_name = "other",
        long_name = "Other land use fraction", units = "-",
        description = "Fraction of catchment with other land use (rainfed crops, grassland, etc.)"
    ),
    list(
        name = "fracrice", col_name = "rice",
        long_name = "Rice paddy fraction", units = "-",
        description = "Fraction of catchment with rice paddies"
    ),
    list(
        name = "fracsealed", col_name = "sealed",
        long_name = "Sealed (impervious) fraction", units = "-",
        description = "Fraction of catchment with impervious surfaces (urban, roads)"
    ),
    list(
        name = "fracwater", col_name = "water",
        long_name = "Water body fraction", units = "-",
        description = "Fraction of catchment covered by water bodies (lakes, reservoirs)"
    )
)

# --- Save variable description ------------------------------------------------
var_desc <- data.table(
    variable_name = sapply(lu_types, `[[`, "col_name"),
    long_name = sapply(lu_types, `[[`, "long_name"),
    units = sapply(lu_types, `[[`, "units"),
    description = sapply(lu_types, `[[`, "description"),
    source_pattern = paste0(sapply(lu_types, `[[`, "name"), "_European_01min_{YEAR}.nc"),
    temporal_resolution = "yearly (1950-2020)",
    spatial_aggregation = "area-weighted mean over catchment polygon"
)
fwrite(var_desc, file.path(out_attr, "landuse_variable_description.csv"))
cat("Saved landuse_variable_description.csv\n")

# --- Try loading from pre-computed aggregates ---------------------------------
cat("Checking for pre-computed land use aggregates...\n")

lu_data <- list()
use_precomputed <- TRUE

for (lu in lu_types) {
    agg_path <- file.path(agg_lu_dir, paste0(lu$name, "_catchment_yearly.csv"))
    if (file.exists(agg_path)) {
        lu_data[[lu$col_name]] <- fread(agg_path)
        cat("  Loaded:", lu$name, "\n")
    } else {
        use_precomputed <- FALSE
        cat("  NOT FOUND:", agg_path, "\n")
        break
    }
}

if (use_precomputed) {
    # --- Export from pre-computed aggregates -----------------------------------
    cat("\nUsing pre-computed aggregates.\n")

    # Get catchment IDs (columns except 'year')
    catch_ids <- setdiff(names(lu_data[[1]]), "year")
    years <- lu_data[[1]]$year
    n_catch <- length(catch_ids)
    cat("Catchments:", n_catch, "| Years:", min(years), "-", max(years), "\n")

    cat("Exporting CAMELS-style land use CSVs...\n")
    for (ic in seq_along(catch_ids)) {
        cid <- catch_ids[ic]

        if (ic %% 200 == 0 || ic == 1) {
            cat("  [", ic, "/", n_catch, "]", cid, "\n")
        }

        out_path <- file.path(out_dir, paste0(cid, "_landuse.csv"))
        if (file.exists(out_path)) next

        # Build one data.table with year + all land use columns
        ts_dt <- data.table(year = years)
        for (lu in lu_types) {
            ts_dt[, (lu$col_name) := lu_data[[lu$col_name]][[cid]]]
        }

        fwrite(ts_dt, out_path)
    }
} else {
    # --- Extract directly from NetCDF rasters ---------------------------------
    cat("\nPre-computed aggregates not available. Extracting from rasters...\n")

    # Load catchments
    catchments <- st_read(gpkg_path, quiet = TRUE)
    if (!isTRUE(st_crs(catchments) == st_crs(4326))) {
        catchments <- st_transform(catchments, 4326)
    }
    catch_ids <- as.character(catchments$catch_id)
    n_catch <- length(catch_ids)

    # Get years from first land use type
    pattern <- paste0("^", lu_types[[1]]$name, "_European_01min_\\d{4}\\.nc$")
    nc_files <- list.files(landuse_dir, pattern = pattern, full.names = TRUE)
    years <- sort(as.integer(sub(".*_(\\d{4})\\.nc$", "\\1", nc_files)))
    n_years <- length(years)
    cat("Years:", min(years), "-", max(years), "(", n_years, ")\n")
    cat("Catchments:", n_catch, "\n")
    cat("This will take a while...\n\n")

    # Initialize storage: list of matrices (years x catchments) per LU type
    lu_mats <- list()
    for (lu in lu_types) {
        lu_mats[[lu$col_name]] <- matrix(NA_real_, nrow = n_years, ncol = n_catch)
    }

    # Extract year by year
    for (iy in seq_along(years)) {
        yr <- years[iy]
        if (iy %% 5 == 1) cat("  Year", yr, "(", iy, "/", n_years, ")\n")

        for (lu in lu_types) {
            nc_path <- file.path(
                landuse_dir,
                paste0(lu$name, "_European_01min_", yr, ".nc")
            )
            if (!file.exists(nc_path)) next

            r <- tryCatch(terra::rast(nc_path), error = function(e) NULL)
            if (is.null(r)) next
            if (nlyr(r) > 1) r <- r[[1]]

            means <- exact_extract(r, catchments, "mean")
            lu_mats[[lu$col_name]][iy, ] <- means
        }
    }

    # Export per catchment
    cat("\nExporting CSVs...\n")
    for (ic in seq_along(catch_ids)) {
        cid <- catch_ids[ic]

        if (ic %% 200 == 0 || ic == 1) {
            cat("  [", ic, "/", n_catch, "]", cid, "\n")
        }

        out_path <- file.path(out_dir, paste0(cid, "_landuse.csv"))
        if (file.exists(out_path)) next

        ts_dt <- data.table(year = years)
        for (lu in lu_types) {
            ts_dt[, (lu$col_name) := lu_mats[[lu$col_name]][, ic]]
        }

        fwrite(ts_dt, out_path)
    }
}

cat("\nDone! Land use CSVs exported to:", out_dir, "\n")
cat("Structure:\n")
cat("  {catch_id}_landuse.csv with columns: year, forest, irrigated, other, rice, sealed, water\n")
