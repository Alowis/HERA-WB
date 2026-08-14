# =============================================================================
# Path configuration
# Source this file at the top of any script: source("config/paths.R")
#
# Edit the three user-defined paths below. All other paths are derived from
# base_dir automatically.
#
# Expected layout of base_dir:
#   HERA-WB_data/
#     +-- aggregates/
#     +-- tss_postprocess/
#     +-- attributes/
#     ¦   +-- catchment_attributes.gpkg
#     +-- GLEAM/
#     +-- globsnow_swe_daily/
#     +-- river_discharge/
#     +-- koppen_geiger_0p1.tif
#     +-- dem.nc
#     +-- output/
#     +-- plots/
# =============================================================================

# =============================================================================
# USER-DEFINED PATHS — edit these three lines to match your environment
# =============================================================================

# Root of the HERA-WB data folder (contains aggregates/, attributes/, etc.)
# data_dir <- "/path/to/HERA-WB_data"
data_dir <- "D:/tilloal/Documents/01_Projects/HERA-WB_data"

# base_dir <- "/path/to/folder"
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/data"


# LISFLOOD utilities and static data
# hydro_dir <- "/path/to/LFRuns_utils/data"
hydro_dir <- "D:/tilloal/Documents/LFRuns_utils/data"

# Land use fraction rasters
# landuse_dir <- "/path/to/landuse"
landuse_dir <- file.path(data_dir, "landuse")

# =============================================================================
# DERIVED PATHS (do not edit below this line)
# =============================================================================
agg_dir <- file.path(base_dir, "aggregates")
tss_dir <- file.path(base_dir, "tss_postprocess")
tss_raw <- file.path(base_dir, "tss")
gpkg_path <- file.path(data_dir, "attributes", "catchment_attributes.gpkg")
out_dir <- file.path(base_dir, "output")
plot_dir <- file.path(base_dir, "plots")

# --- Satellite data paths ---
gleam_daily_dir <- file.path(base_dir, "GLEAM", "AET_daily")
gleam_monthly_dir <- file.path(base_dir, "GLEAM", "AET_monthly")
globsnow_dir <- file.path(base_dir, "globsnow_swe_daily")

# --- Validation reference data ---
discharge_dir <- file.path(base_dir, "river_discharge")

# --- Climate / elevation rasters ---
path_clim <- file.path(data_dir, "koppen_geiger_0p1.tif")
path_dem <- file.path(data_dir, "dem.nc")


# --- Analysis parameters ---
start_year <- 1951
end_year <- 2020
year_range <- start_year:end_year
