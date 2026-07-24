# =============================================================================
# Export HERA data in CAMELS-style CSV format
#
# Structure (following CAMELS-DE conventions):
#   timeseries/
#     {catch_id}_6hourly.csv     — 6-hourly time series (all variables)
#   attributes/
#     catchment_attributes.csv   — static attributes (area, coordinates, etc.)
#     variable_description.csv   — metadata for each time series variable
#
# Each timeseries CSV has columns:
#   date, hour, rainfall, snowfall, snowmelt, infiltration, ActEvapo,
#   surfaceRunoff, prefFlow, dSubToUz, percUZLZ, quz, qlz, GWloss,
#   theta1, theta2, theta3, snowSWE
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(sf)
library(lubridate)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
tss_dir <- file.path(base_dir, "data", "tss_postprocess")
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
out_base <- file.path(base_dir, "data", "HERA_CAMELS")
out_ts <- file.path(out_base, "timeseries")
out_attr <- file.path(out_base, "attributes")
dir.create(out_ts, showWarnings = FALSE, recursive = TRUE)
dir.create(out_attr, showWarnings = FALSE, recursive = TRUE)

# --- Variable definitions -----------------------------------------------------
variables <- list(
    list(
        name = "rainfall", file = "rainUpsX_nested_1951_2020.csv",
        long_name = "Rainfall", units = "mm/6h",
        description = "Liquid precipitation reaching the surface"
    ),
    list(
        name = "snowfall", file = "snowUpsX_nested_1951_2020.csv",
        long_name = "Snowfall", units = "mm/6h",
        description = "Solid precipitation (snow)"
    ),
    list(
        name = "snowmelt", file = "snowMeltUpsX_nested_1951_2020.csv",
        long_name = "Snowmelt", units = "mm/6h",
        description = "Melt from snowpack"
    ),
    list(
        name = "infiltration", file = "infUpsX_nested_1951_2020.csv",
        long_name = "Infiltration", units = "mm/6h",
        description = "Water infiltrating from surface into soil layers"
    ),
    list(
        name = "ActEvapo", file = "ActEvapo_nested_1951_2020.csv",
        long_name = "Actual evapotranspiration", units = "mm/6h",
        description = "Actual evapotranspiration (soil evap + transpiration + interception)"
    ),
    list(
        name = "surfaceRunoff", file = "surfaceRunoffUpsX_nested_1951_2020.csv",
        long_name = "Direct surface runoff", units = "mm/6h",
        description = "Water running off the surface directly to the channel"
    ),
    list(
        name = "prefFlow", file = "prefFlowUpsX_nested_1951_2020.csv",
        long_name = "Preferential flow", units = "mm/6h",
        description = "Bypass flow from surface directly to upper groundwater zone"
    ),
    list(
        name = "dSubToUz", file = "dSubToUzUpsX_nested_1951_2020.csv",
        long_name = "Soil to upper zone drainage", units = "mm/6h",
        description = "Drainage from lower soil layer to upper groundwater zone"
    ),
    list(
        name = "percUZLZ", file = "percUZLZUpsX_nested_1951_2020.csv",
        long_name = "Percolation UZ to LZ", units = "mm/6h",
        description = "Percolation from upper to lower groundwater zone"
    ),
    list(
        name = "quz", file = "qUzUpsX_nested_1951_2020.csv",
        long_name = "Upper zone outflow", units = "mm/6h",
        description = "Subsurface flow from upper groundwater zone to channel"
    ),
    list(
        name = "qlz", file = "qLZUpsX_nested_1951_2020.csv",
        long_name = "Lower zone outflow (baseflow)", units = "mm/6h",
        description = "Baseflow from lower groundwater zone to channel"
    ),
    list(
        name = "GWloss", file = "lossUpsX_nested_1951_2020.csv",
        long_name = "Groundwater loss", units = "mm/6h",
        description = "Water lost from lower zone to deep aquifers (not reaching outlet)"
    ),
    list(
        name = "theta1", file = "theta1totalX_nested_1951_2020.csv",
        long_name = "Upper soil moisture", units = "mm",
        description = "Water stored in upper (superficial) soil layer"
    ),
    list(
        name = "theta2", file = "theta2totalX_nested_1951_2020.csv",
        long_name = "Middle soil moisture", units = "mm",
        description = "Water stored in middle (upper) soil layer"
    ),
    list(
        name = "theta3", file = "theta3totalX_nested_1951_2020.csv",
        long_name = "Deep soil moisture", units = "mm",
        description = "Water stored in deep (lower) soil layer"
    ),
    list(
        name = "snowSWE", file = "scovUps_nested_1951_2020.csv",
        long_name = "Snow water equivalent", units = "mm",
        description = "Water stored in the snowpack"
    )
)

# --- Save variable description CSV --------------------------------------------
var_desc <- data.table(
    variable_name = sapply(variables, `[[`, "name"),
    long_name = sapply(variables, `[[`, "long_name"),
    units = sapply(variables, `[[`, "units"),
    description = sapply(variables, `[[`, "description"),
    source_file = sapply(variables, `[[`, "file"),
    temporal_resolution = "6-hourly",
    spatial_aggregation = "upstream catchment area-weighted mean"
)
fwrite(var_desc, file.path(out_attr, "variable_description.csv"))
cat("Saved variable_description.csv\n")

# --- Load catchment attributes and save ---------------------------------------
cat("Loading catchment attributes...\n")
catchments <- st_read(gpkg_path, quiet = TRUE)

# Compute centroids for lat/lon
centroids <- st_centroid(st_transform(catchments, 4326))
coords <- st_coordinates(centroids)

catch_attr <- data.table(
    catch_id = catchments$catch_id,
    longitude = coords[, 1],
    latitude = coords[, 2],
    area_km2 = catchments$area_km2,
    residual_area_km2 = catchments$residual_area_km2
)

# Add any other available fields from the gpkg
extra_fields <- intersect(
    c("immediate_nested_ids", "n_nested", "is_headwater"),
    names(catchments)
)
for (f in extra_fields) {
    catch_attr[, (f) := catchments[[f]]]
}

fwrite(catch_attr, file.path(out_attr, "catchment_attributes.csv"))
cat("Saved catchment_attributes.csv (", nrow(catch_attr), "catchments)\n")

# --- Check which TSS files exist ----------------------------------------------
cat("Checking available TSS files...\n")
available_vars <- list()
for (v in variables) {
    fpath <- file.path(tss_dir, v$file)
    if (file.exists(fpath)) {
        available_vars[[v$name]] <- v
    } else {
        cat("  MISSING:", v$file, "\n")
    }
}
cat("Available:", length(available_vars), "/", length(variables), "variables\n")

# --- Get catchment IDs from first file header ---------------------------------
first_file <- file.path(tss_dir, available_vars[[1]]$file)
header <- names(fread(first_file, nrows = 0, header = TRUE))
catch_ids <- header
n_catch <- length(catch_ids)
cat("Catchments in data:", n_catch, "\n")

# --- Time vector (6-hourly) ---------------------------------------------------
n_days <- as.integer(as.Date("2020-12-31") - as.Date("1951-01-01")) + 1
n_steps_expected <- n_days * 4

# Check actual row count
test_nrows <- nrow(fread(first_file, select = 1L, header = TRUE))
if (test_nrows == n_days) {
    cat("Data is DAILY (not 6-hourly).\n")
    dates <- seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")
    hours <- rep(0L, length(dates))
    n_steps <- n_days
} else {
    cat("Data is 6-HOURLY.\n")
    dates <- rep(seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day"), each = 4)
    hours <- rep(c(0L, 6L, 12L, 18L), times = n_days)
    n_steps <- n_steps_expected
}

# --- Export one CSV per catchment ---------------------------------------------
cat("\nExporting CAMELS-style CSVs (", n_catch, "catchments)...\n")
cat("Reading", length(available_vars), "TSS files per catchment.\n\n")

for (ic in seq_along(catch_ids)) {
    cid <- catch_ids[ic]

    if (ic %% 50 == 0 || ic == 1) {
        cat("  [", ic, "/", n_catch, "]", cid, "\n")
    }

    out_path <- file.path(out_ts, paste0(cid, "_timeseries.csv"))

    # Skip if already exists
    if (file.exists(out_path)) next

    # Start with date/hour columns
    ts_dt <- data.table(date = dates, hour = hours)

    # Read each variable's column
    for (vname in names(available_vars)) {
        vinfo <- available_vars[[vname]]
        fpath <- file.path(tss_dir, vinfo$file)

        col_data <- tryCatch(
            fread(fpath, select = cid, header = TRUE)[[1]],
            error = function(e) {
                tryCatch(
                    fread(fpath, select = paste0("X", cid), header = TRUE)[[1]],
                    error = function(e2) rep(NA_real_, n_steps)
                )
            }
        )

        if (length(col_data) != n_steps) {
            col_data <- rep(NA_real_, n_steps)
        }

        ts_dt[, (vname) := col_data]
    }

    # Write CSV
    fwrite(ts_dt, out_path)
}

cat("\nDone! CAMELS-style dataset exported to:", out_base, "\n")
cat("Structure:\n")
cat("  ", out_base, "/\n")
cat("    attributes/\n")
cat("      catchment_attributes.csv\n")
cat("      variable_description.csv\n")
cat("    timeseries/\n")
cat("      {catch_id}_timeseries.csv  (one per catchment)\n")
