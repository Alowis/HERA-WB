# =============================================================================
# Build shared analysis workspace
#
# Run this ONCE after preprocessing (deaggregation + aggregation).
# All downstream scripts load the result with:
#   ws <- readRDS(file.path(base_dir, "output", "workspace.rds"))
#
# Contains:
#   - Catchments (sf + projected) with Iceland removed
#   - Area weights (residual_area_km2, normalized)
#   - All monthly aggregates (with reliable dates)
#   - Homogenized satellite products (GLEAM, ESA CCI, GlobSnow)
#   - Snow mask (monthly: which catchments/months have snow)
#   - Continental area-weighted annual means (all variables)
#   - Basemap for plotting
#
# Output: output/workspace.rds (~200-500 MB depending on data)
# =============================================================================

library(data.table)
library(sf)
library(lubridate)
library(rnaturalearth)

# --- Configuration -----------------------------------------------------------
source("config/paths.R")

cat("=== Building shared workspace ===\n")
cat("Base dir:", base_dir, "\n\n")

# =============================================================================
# 1. CATCHMENTS + AREA WEIGHTS
# =============================================================================
cat("[1/7] Loading catchments...\n")

catchments <- st_read(gpkg_path, quiet = TRUE)

# Remove Iceland (for snow analysis only — stored as separate mask)
centroids_wgs <- st_coordinates(st_centroid(st_transform(catchments, 4326)))
iceland_mask <- centroids_wgs[, 1] < -13 & centroids_wgs[, 2] > 63
iceland_ids <- catch_ids[iceland_mask]

# Keep ALL catchments in the workspace (Iceland included)
catchments_3035 <- st_transform(catchments, 3035)

# Catch IDs and area weights
# Catch IDs and area weights (all catchments, Iceland included)
catch_ids <- as.character(as.numeric(catchments$catch_id))
area_vec <- setNames(catchments$residual_area_km2, catch_ids)
weights <- area_vec / sum(area_vec, na.rm = TRUE)

cat("  Catchments:", length(catch_ids), "(all, including Iceland)\n")
cat("  Iceland catchments:", length(iceland_ids), "(excluded in snow analysis)\n")
cat("  Total area:", round(sum(area_vec, na.rm = TRUE)), "km²\n")

# =============================================================================
# 2. BASEMAP
# =============================================================================
cat("[2/7] Preparing basemap...\n")

basemap <- ne_countries(scale = "medium", returnclass = "sf") |>
    st_transform(3035)
bbox <- st_bbox(catchments_3035)

# =============================================================================
# 3. MONTHLY AGGREGATES (all variables, with reliable dates)
# =============================================================================
cat("[3/7] Loading monthly aggregates...\n")

load_monthly_agg <- function(var_name) {
    path <- file.path(agg_dir, var_name, paste0(var_name, "_monthly_all_years.csv"))
    if (!file.exists(path)) {
        cat("    SKIP:", var_name, "(not found)\n")
        return(NULL)
    }
    dt <- fread(path)
    # Assign reliable dates from row position (month_idx unreliable)
    dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
    dt[, year := year(date)]
    dt[, month := month(date)]
    cat("    OK:", var_name, "(", nrow(dt), "months)\n")
    dt
}

monthly_vars <- c(
    "ActEvapo", "surface_soil_moisture", "root_soil_moisture",
    "snow_water_equivalent", "runoff", "qlz", "quz", "qutl",
    "snowfall", "snowmelt", "theta3", "percolation", "Q","rainfall","prefflow","GWloss"
)

monthly <- setNames(lapply(monthly_vars, load_monthly_agg), monthly_vars)
# Remove NULLs
monthly <- monthly[!sapply(monthly, is.null)]
cat("  Loaded:", length(monthly), "variables\n")

# =============================================================================
# 4. HOMOGENIZED SATELLITE DATA
# =============================================================================
cat("[4/7] Loading homogenized satellite products...\n")

load_homog <- function(path, date_fmt = "YYYY-MM") {
    if (!file.exists(path)) {
        cat("    SKIP:", basename(path), "\n")
        return(NULL)
    }
    dt <- fread(path, header = TRUE)
    # Parse date column
    if ("date" %in% names(dt)) {
        if (date_fmt == "YYYY-MM") {
            dt[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
        } else {
            dt[, date := as.Date(date)]
        }
        dt[, year := year(date)]
        dt[, month := month(date)]
    }
    cat("    OK:", basename(path), "(", nrow(dt), "rows)\n")
    dt
}

homog <- list(
    gleam_monthly = load_homog(file.path(base_dir, "output/aet_diego/1.homogenized/gleam_monthly_homog.csv")),
    lisflood_aet_monthly = load_homog(file.path(base_dir, "output/aet_diego/1.homogenized/lisflood_monthly_homog.csv")),
    cci_monthly = load_homog(file.path(base_dir, "output/soil_moisture_diego/1.Diego_Merged/esacci_homogenized.csv"), "date"),
    lisflood_sm_monthly = load_homog(file.path(base_dir, "output/soil_moisture_diego/1.Diego_Merged/lisflood_homogenized.csv"), "date"),
    globsnow_monthly = load_homog(file.path(base_dir, "output/swe_diego/1.homogenized/globsnow_monthly_homog.csv")),
    lisflood_swe_monthly = load_homog(file.path(base_dir, "output/swe_diego/1.homogenized/lisflood_monthly_homog.csv"))
)
homog <- homog[!sapply(homog, is.null)]
cat("  Loaded:", length(homog), "products\n")

# =============================================================================
# 5. SNOW MASK (monthly + daily SWE for detailed masking)
# =============================================================================
cat("[5/7] Computing snow mask + loading daily SWE...\n")

snow_mask <- NULL
swe_daily <- NULL

if ("snow_water_equivalent" %in% names(monthly)) {
    swe_dt <- monthly$snow_water_equivalent
    meta_cols <- c("month_idx", "period_start", "period_end", "date", "year", "month")
    swe_catch_cols <- intersect(catch_ids, setdiff(names(swe_dt), meta_cols))

    # Monthly snow mask: TRUE where SWE > 0
    snow_mask <- swe_dt[, lapply(.SD, function(x) x > 0), .SDcols = swe_catch_cols]
    snow_mask[, date := swe_dt$date]
    snow_mask[, year := swe_dt$year]
    snow_mask[, month := swe_dt$month]

    snow_pct <- colMeans(snow_mask[, ..swe_catch_cols], na.rm = TRUE)
    cat("  Monthly snow mask:", length(swe_catch_cols), "catchments\n")
    cat("  Catchments with >5% snow months:", sum(snow_pct > 0.05), "\n")
}

# Load daily SWE (for per-day snow masking in SM/discharge scripts)
swe_daily_path <- file.path(base_dir, "output/swe_diego/1.homogenized/lisflood_daily_homog.csv")
if (!file.exists(swe_daily_path)) {
    swe_daily_path <- file.path(
        agg_dir, "snow_water_equivalent",
        "snow_water_equivalent_daily_all_years.csv"
    )
}
if (file.exists(swe_daily_path)) {
    cat("  Loading daily SWE:", basename(swe_daily_path), "\n")
    swe_daily <- fread(swe_daily_path, header = TRUE)
    if ("date" %in% names(swe_daily)) {
        swe_daily[, date := as.IDate(date)]
    } else if ("period_start" %in% names(swe_daily)) {
        swe_daily[, date := as.IDate(period_start)]
    } else {
        swe_daily[, date := seq.Date(as.Date("1951-01-01"), by = "day", length.out = .N)]
    }
    cat("  Daily SWE:", nrow(swe_daily), "days x", ncol(swe_daily) - 1, "catchments\n")
} else {
    cat("  WARNING: No daily SWE found. Daily snow masking unavailable.\n")
}

# =============================================================================
# 6. CONTINENTAL AREA-WEIGHTED ANNUAL MEANS
# =============================================================================
cat("[6/7] Computing continental annual means...\n")

compute_continental_annual <- function(dt, var_name) {
    meta_cols <- c("month_idx", "period_start", "period_end", "date", "year", "month")
    cols <- intersect(catch_ids, setdiff(names(dt), meta_cols))
    if (length(cols) < 10) {
        return(NULL)
    }

    w <- weights[cols]
    w <- w / sum(w, na.rm = TRUE)

    # Compute weighted mean per timestep (NA-robust)
    mat <- as.matrix(dt[, ..cols])
    w_mat <- matrix(w, nrow = nrow(mat), ncol = length(w), byrow = TRUE)
    w_mat[is.na(mat)] <- 0
    mat[is.na(mat)] <- 0
    weighted_mean <- rowSums(mat * w_mat) / rowSums(w_mat)

    # Annual mean
    annual <- data.table(
        year = dt$year,
        monthly_val = weighted_mean
    )[, .(annual_mean = mean(monthly_val, na.rm = TRUE)), by = year]
    annual[, variable := var_name]
    annual
}

continental_annual <- rbindlist(
    lapply(names(monthly), function(nm) compute_continental_annual(monthly[[nm]], nm)),
    fill = TRUE
)
continental_annual <- continental_annual[!is.na(annual_mean)]
cat("  Variables with continental means:", length(unique(continental_annual$variable)), "\n")
cat("  Years:", range(continental_annual$year), "\n")

# =============================================================================
# 7. SAVE WORKSPACE
# =============================================================================
cat("[7/7] Saving workspace...\n")

workspace <- list(
    # Spatial
    catchments = catchments,
    catchments_3035 = catchments_3035,
    basemap = basemap,
    bbox = bbox,

    # IDs and weights
    catch_ids = catch_ids,
    iceland_ids = iceland_ids,
    area_vec = area_vec,
    weights = weights,

    # Monthly time series (all variables)
    monthly = monthly,

    # Homogenized satellite products
    homog = homog,

    # Snow mask
    snow_mask = snow_mask,
    swe_daily = swe_daily,

    # Precomputed continental means
    continental_annual = continental_annual,

    # Metadata
    meta = list(
        created = Sys.time(),
        n_catchments = length(catch_ids),
        n_variables = length(monthly),
        year_range = c(1951, 2020)
    )
)

out_path <- file.path(base_dir, "output", "workspace.rds")
saveRDS(workspace, out_path)

file_size_mb <- round(file.size(out_path) / 1e6, 1)
cat("\n=== Workspace saved ===\n")
cat("  Path:", out_path, "\n")
cat("  Size:", file_size_mb, "MB\n")
cat("  Contents:\n")
cat("    - catchments:", length(catch_ids), "\n")
cat("    - monthly variables:", length(monthly), "\n")
cat("    - homogenized products:", length(homog), "\n")
cat("    - snow mask:", if (!is.null(snow_mask)) "yes" else "no", "\n")
cat("    - continental annual means:", nrow(continental_annual), "rows\n")
cat("\n  Usage in downstream scripts:\n")
cat("    ws <- readRDS('output/workspace.rds')\n")
cat("    catchments <- ws$catchments\n")
cat("    weights <- ws$weights\n")
cat("    monthly_aet <- ws$monthly$ActEvapo\n")
cat("    snow_mask <- ws$snow_mask\n")

