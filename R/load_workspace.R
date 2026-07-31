# =============================================================================
# Shared workspace loader
# Source this at the top of any analysis script:
#   source("R/load_workspace.R")
#
# Provides:
#   ws            — full workspace list
#   base_dir      — project root path
#   catchments    — sf object (WGS84)
#   catchments_3035 — sf object (EPSG:3035)
#   basemap       — European basemap (EPSG:3035)
#   bbox          — bounding box for maps
#   catch_ids     — character vector of catchment IDs
#   weights       — named numeric vector of normalized area weights
#   monthly       — named list of monthly aggregate data.tables
#   homog         — named list of homogenized satellite data.tables
#   snow_mask     — monthly snow presence (logical data.table)
# =============================================================================

# Load configuration
# Requires RStudio project or setwd() to RegimeShift_codes/
# (the .Rproj file ensures this when opening in RStudio)

source("config/paths.R")

# Load workspace
ws_path <- file.path(base_dir, "output", "workspace.rds")
if (!file.exists(ws_path)) {
    stop(
        "Workspace not found at: ", ws_path,
        "\nRun scripts/00_preprocessing/build_workspace.R first."
    )
}

cat("Loading workspace...\n")
ws <- readRDS(ws_path)

# Unpack commonly used objects into global environment
catchments <- ws$catchments
catchments_3035 <- ws$catchments_3035
basemap <- ws$basemap
bbox <- ws$bbox
catch_ids <- ws$catch_ids
area_vec <- ws$area_vec
weights <- ws$weights
monthly <- ws$monthly
homog <- ws$homog
snow_mask <- ws$snow_mask
swe_daily <- ws$swe_daily

cat(
    "  Workspace loaded:", ws$meta$n_catchments, "catchments,",
    ws$meta$n_variables, "variables\n"
)
