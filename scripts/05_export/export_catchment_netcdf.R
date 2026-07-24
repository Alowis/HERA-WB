# =============================================================================
# Export all hydrological variables into one NetCDF file per catchment
#
# Each NetCDF contains monthly time series of all available variables with
# proper metadata (long name, units, description).
#
# Output: data/netcdf_by_catchment/{catch_id}.nc
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(ncdf4)
library(lubridate)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
agg_dir <- file.path(base_dir, "data", "aggregates")
tss_dir <- file.path(base_dir, "data", "tss_postprocess")
out_dir <- file.path(base_dir, "data", "netcdf_by_catchment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Variable definitions -----------------------------------------------------
# name: folder name in aggregates (or TSS fallback)
# file_monthly: filename of the monthly aggregate CSV
# tss_file: raw TSS file (fallback if no monthly aggregate)
# long_name: CF-style long name
# units: CF-style units
# agg_method: how sub-daily is aggregated ("sum" for fluxes, "mean" for states)
variables <- list(
    list(
        name = "ActEvapo", file_monthly = "ActEvapo_monthly_all_years.csv",
        long_name = "Actual evapotranspiration", units = "mm/month", agg_method = "sum"
    ),
    list(
        name = "surface_soil_moisture", file_monthly = "surface_soil_moisture_monthly_all_years.csv",
        long_name = "Surface soil moisture", units = "mm", agg_method = "mean"
    ),
    list(
        name = "root_soil_moisture", file_monthly = "root_soil_moisture_monthly_all_years.csv",
        long_name = "Root zone soil moisture", units = "mm", agg_method = "mean"
    ),
    list(
        name = "snow_water_equivalent", file_monthly = "snow_water_equivalent_monthly_all_years.csv",
        long_name = "Snow water equivalent", units = "mm", agg_method = "mean"
    ),
    list(
        name = "runoff", file_monthly = "runoff_monthly_all_years.csv",
        long_name = "Direct surface runoff", units = "mm/month", agg_method = "sum"
    ),
    list(
        name = "qlz", file_monthly = "qlz_monthly_all_years.csv",
        long_name = "Baseflow from lower groundwater zone", units = "mm/month", agg_method = "sum"
    ),
    list(
        name = "quz", file_monthly = "quz_monthly_all_years.csv",
        long_name = "Flow from upper groundwater zone", units = "mm/month", agg_method = "sum"
    ),
    list(
        name = "qutl", file_monthly = "qutl_monthly_all_years.csv",
        long_name = "Percolation from upper to lower zone", units = "mm/month", agg_method = "sum"
    ),
    list(
        name = "Q", file_monthly = "Q_monthly_all_years.csv",
        long_name = "Total channel discharge", units = "mm/month", agg_method = "sum"
    ),
    list(
        name = "snowfall", file_monthly = "snowfall_monthly_all_years.csv",
        tss_file = "snowUpsX_nested_1951_2020.csv",
        long_name = "Snowfall", units = "mm/month", agg_method = "sum"
    ),
    list(
        name = "snowmelt", file_monthly = "snowmelt_monthly_all_years.csv",
        tss_file = "snowMeltUpsX_nested_1951_2020.csv",
        long_name = "Snowmelt", units = "mm/month", agg_method = "sum"
    ),
    list(
        name = "theta3", file_monthly = "theta3_monthly_all_years.csv",
        tss_file = "theta3totalX_nested_1951_2020.csv",
        long_name = "Deep soil layer moisture (theta3)", units = "mm", agg_method = "mean"
    ),
    list(
        name = "percolation", file_monthly = "percolation_monthly_all_years.csv",
        tss_file = "percUZLZUpsX_nested_1951_2020.csv",
        long_name = "Percolation from upper to lower groundwater zone", units = "mm/month", agg_method = "sum"
    ),
    list(
        name = "tha1", file_monthly = "tha1_monthly_all_years.csv",
        long_name = "Upper soil layer moisture (theta1)", units = "mm", agg_method = "mean"
    )
)

# --- Load all monthly datasets ------------------------------------------------
cat("Loading monthly aggregate datasets...\n")

loaded_vars <- list()
common_catches <- NULL

for (v in variables) {
    agg_path <- file.path(agg_dir, v$name, v$file_monthly)
    if (!file.exists(agg_path)) {
        cat("  SKIP:", v$name, "(file not found:", agg_path, ")\n")
        next
    }
    cat("  Loading:", v$name, "\n")
    dt <- fread(agg_path)
    # Assign reliable dates from row position
    dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]

    meta <- c("month_idx", "period_start", "period_end", "date", "window")
    catch_cols <- setdiff(names(dt), meta)

    # Track common catchments
    if (is.null(common_catches)) {
        common_catches <- catch_cols
    } else {
        common_catches <- intersect(common_catches, catch_cols)
    }

    loaded_vars[[v$name]] <- list(
        data = dt,
        catch_cols = catch_cols,
        long_name = v$long_name,
        units = v$units
    )
}

cat("\nVariables loaded:", length(loaded_vars), "\n")
cat("Common catchments:", length(common_catches), "\n")

# --- Time dimension -----------------------------------------------------------
# All datasets should have same number of months (840 = 70 years * 12)
dates <- loaded_vars[[1]]$data$date
n_time <- length(dates)
cat("Time steps:", n_time, "months (", year(dates[1]), "-", year(dates[n_time]), ")\n")

# Time as days since 1951-01-01 (CF convention)
time_vals <- as.numeric(dates - as.Date("1951-01-01"))

# --- Export one NetCDF per catchment ------------------------------------------
cat("\nExporting NetCDF files...\n")
n_catch <- length(common_catches)

for (ic in seq_along(common_catches)) {
    cid <- common_catches[ic]

    if (ic %% 100 == 0 || ic == 1) {
        cat("  Progress:", ic, "/", n_catch, "(", cid, ")\n")
    }

    nc_path <- file.path(out_dir, paste0(cid, ".nc"))

    # Define time dimension
    time_dim <- ncdim_def("time", "days since 1951-01-01",
        vals = time_vals, unlim = TRUE,
        calendar = "standard"
    )

    # Define variables
    nc_vars <- list()
    for (vname in names(loaded_vars)) {
        vinfo <- loaded_vars[[vname]]
        nc_vars[[vname]] <- ncvar_def(
            name = vname,
            units = vinfo$units,
            dim = list(time_dim),
            missval = -9999,
            longname = vinfo$long_name
        )
    }

    # Create NetCDF file
    nc <- nc_create(nc_path, nc_vars, force_v4 = TRUE)

    # Write data for each variable
    for (vname in names(loaded_vars)) {
        vinfo <- loaded_vars[[vname]]
        vals <- vinfo$data[[cid]]
        if (is.null(vals)) {
            vals <- rep(NA_real_, n_time)
        }
        # Replace NA with missval for cleanliness
        vals[is.na(vals)] <- -9999
        ncvar_put(nc, nc_vars[[vname]], vals)
    }

    # Global attributes
    ncatt_put(
        nc, 0, "title",
        paste("LISFLOOD monthly hydrological variables — Catchment", cid)
    )
    ncatt_put(nc, 0, "institution", "JRC — Joint Research Centre")
    ncatt_put(nc, 0, "source", "HERA pan-European hydrological reanalysis (1951-2020)")
    ncatt_put(nc, 0, "catchment_id", cid)
    ncatt_put(nc, 0, "temporal_resolution", "monthly")
    ncatt_put(nc, 0, "spatial_coverage", "upstream catchment area-weighted mean")
    ncatt_put(
        nc, 0, "history",
        paste("Created on", Sys.Date(), "by export_catchment_netcdf.R")
    )
    ncatt_put(nc, 0, "conventions", "CF-1.8")

    # Time attributes
    ncatt_put(nc, "time", "calendar", "standard")
    ncatt_put(nc, "time", "axis", "T")

    nc_close(nc)
}

cat("\nDone! Exported", n_catch, "NetCDF files to:", out_dir, "\n")
cat("Each file contains", length(loaded_vars), "variables ×", n_time, "months.\n")
