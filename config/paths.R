# =============================================================================
# Path configuration — edit these to match your local/HPC environment
# Source this file at the top of any script: source("config/paths.R")
# =============================================================================

# --- Base project directory ---
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"

# --- LISFLOOD utilities and static data ---
hydro_dir <- "D:/tilloal/Documents/LFRuns_utils/data/"

# --- Land use rasters ---
landuse_dir <- "D:/tilloal/Documents/06_Floodrivers/landuse/"

# --- Derived paths (do not edit) ---
agg_dir <- file.path(base_dir, "data", "aggregates")
tss_dir <- file.path(base_dir, "data", "tss_postprocess")
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
out_dir <- file.path(base_dir, "output")
plot_dir <- file.path(base_dir, "plots")

# --- Satellite data paths ---
gleam_daily_dir <- file.path(base_dir, "data", "GLEAM", "AET_daily")
gleam_monthly_dir <- file.path(base_dir, "data", "GLEAM", "AET_monthly")
globsnow_dir <- file.path(base_dir, "data", "globsnow_swe_daily")

# --- Climate / elevation rasters ---
path_clim <- file.path(base_dir, "data", "koppen_geiger_0p1.tif")
path_dem <- file.path(base_dir, "data", "dem.nc")

# --- Validation reference data ---
discharge_dir <- file.path(base_dir, "data", "river_discharge")

# --- Analysis parameters ---
start_year <- 1951
end_year <- 2020
year_range <- start_year:end_year
