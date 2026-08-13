# =============================================================================
# deaggregate_and_metrics.R
# =============================================================================
# Clean re-implementation of Statevar_ana_nested.R
#
# Purpose:
#   Deaggregate upstream-averaged LISFLOOD outputs to residual (inter-catchment)
#   values, compute derived climate/hydrologic metrics, produce maps, and
#   batch-export deaggregated CSVs.
#
# Concept:
#   TSS files store spatially-aggregated means over the full upstream area.
#   For a nested outlet i with immediate children j1, j2, ...:
#
#     residual_i = (Value_i * Area_i - Sum(Value_j * Area_j)) / ResArea_i
#
#   For discharge (m3/s), volume mass balance is used instead:
#     residual_i = (Q_i * 86400 - Sum(Q_j * 86400)) / 86400
#
#   Headwater catchments (no nested children) keep their original values.
# =============================================================================

# =============================================================================
# SECTION 0: Libraries
# =============================================================================
library(sf)
library(data.table)
library(lubridate)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)
library(ncdf4)
library(terra)
library(exactextractr)
library(scales)
library(cowplot)

# =============================================================================
# SECTION 1: Configuration
# =============================================================================
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
hydro_dir <- "D:/tilloal/Documents/LFRuns_utils/data/"

# Source helper functions (preprocess_in, process_data, process_precip, etc.)
source(paste0(base_dir, "RegimeShift_codes/R/functions_regime.R"))

# --- Input paths ---
gpkg_path <- paste0(base_dir, "data/catchments_analysis_final_v3.gpkg")
outlet_name <- "outletsv8_hybas07_01min"
uparea_file <- "upArea_European_01min.nc"
efas_file <- "efas_rnet_100km_01min"
path_clim <- paste0(base_dir, "data/koppen_geiger_0p1.tif")
path_dem <- paste0(base_dir, "data/dem.nc")

# --- TSS file paths ---
rain_file <- paste0(base_dir, "data/tss/HERA_SocCF/RainUpsX_1951_2020.csv")
snow_file <- paste0(base_dir, "data/tss/HERA_Histo/SnowUpsX_1951_2020.csv")
aet_file <- paste0(hydro_dir, "tss/HERA_Histo/ActEvapo_1951_2020.csv")
pet_file <- paste0(hydro_dir, "tss/HERA_Histo/etUpsX_1951_2020.csv")

# --- Output paths ---
plot_dir <- file.path(base_dir, "output", "figures")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# --- Analysis parameters ---
year_range <- 1951:2020
dt_step <- 4 # 6-hourly time steps per day

# --- Variables for batch deaggregation ---
batch_vars <- c(
  "rainUpsX", "snowMeltUpsX", "snowUpsX", "scovUps",
  "qUzUpsX", "qLZUpsX", "surfaceRunoffUpsX", "infUpsX", "ActEvapo", "disWin",
  "percUZLZUpsX", "dSubToUzUpsX", "prefFlowUpsX", "lossUpsX", "theta3UpsX",
  "theta2UpsX", "theta1UpsX", "etUpsX"
)

# --- Plot settings ---
tsize <- 12
osize <- 12

# =============================================================================
# SECTION 2: Core Deaggregation Function
# =============================================================================

#' Deaggregate upstream-averaged values to residual (inter-catchment) values
#'
#' @param dt_raw       data.table [nTime x nOutlets], columns named by catch_id
#' @param imm_children Named list: catch_id -> character vector of child ids
#' @param area_lkp     Named numeric vector of upstream area (km2), key = catch_id
#' @param resarea_lkp  Named numeric vector of residual area (km2), key = catch_id
#' @param Q            Logical. If TRUE, treats values as discharge (m3/s) and
#'                     uses volume mass balance instead of area-weighted.
#' @return data.table same shape as dt_raw with residual values
deaggregate_to_residual <- function(dt_raw, imm_children, area_lkp,
                                    resarea_lkp, Q = FALSE) {
  res_dt <- copy(dt_raw)
  col_names <- names(dt_raw)

  # Only process outlets that appear in both the nesting list and data
  outlets_in_data <- intersect(names(imm_children), col_names)

  for (p_id in outlets_in_data) {
    child_ids <- imm_children[[p_id]]

    # Headwater: nothing to subtract
    if (length(child_ids) == 0) next

    p_area <- area_lkp[p_id]
    res_area <- resarea_lkp[p_id]

    if (is.na(res_area) || res_area <= 0) next

    # Keep only children present in this data.table
    valid_children <- intersect(child_ids, col_names)
    if (length(valid_children) == 0) next

    # --- Mass balance ---
    if (Q) {
      # Discharge: convert m3/s to daily volume (m3) for balance
      factor <- 86400
      child_vols <- dt_raw[, .SD * factor, .SDcols = valid_children]
      tot_vol <- dt_raw[[p_id]] * factor
    } else {
      # Area-weighted: value * area
      child_areas <- area_lkp[valid_children]
      child_vols <- dt_raw[, mapply(`*`, .SD, child_areas),
        .SDcols = valid_children
      ]
      tot_vol <- dt_raw[[p_id]] * p_area
    }

    total_child_vol <- if (length(valid_children) > 1) {
      rowSums(child_vols, na.rm = TRUE)
    } else {
      as.vector(child_vols)
    }

    # Compute residual
    if (Q) {
      res_series <- (tot_vol - total_child_vol) / factor
    } else {
      res_series <- (tot_vol - total_child_vol) / res_area
    }

    set(res_dt, j = p_id, value = res_series)
  }

  return(res_dt)
}

# =============================================================================
# SECTION 2b: Lag-Corrected Deaggregation for Discharge
# =============================================================================

#' Estimate optimal lag (in time steps) between parent and child discharge series
#' using cross-correlation. The lag represents the travel time from the child
#' outlet to the parent outlet.
#'
#' @param parent_q   Numeric vector of parent discharge (m3/s)
#' @param child_q    Numeric vector of child discharge (m3/s)
#' @param max_lag    Maximum lag to search (in time steps, default 10 days)
#' @return Integer lag (>= 0) that maximises cross-correlation
estimate_lag <- function(parent_q, child_q, max_lag = 10L) {
  # Remove NAs for correlation computation
  valid <- !is.na(parent_q) & !is.na(child_q)
  if (sum(valid) < 100) {
    return(0L)
  } # not enough data, assume no lag

  p <- parent_q[valid]
  c <- child_q[valid]

  # Normalise to zero-mean for cross-correlation
  p <- p - mean(p)
  c <- c - mean(c)

  # Compute cross-correlation at positive lags only (child leads parent)
  # A positive lag k means: child_q at time (t - k) correlates with parent_q at time t
  best_lag <- 0L
  best_cor <- cor(p, c)

  n <- length(p)
  for (k in seq(0, max_lag)) {
    if (k >= n) break
    # Shift child back by k steps: compare parent[k+1:n] with child[1:(n-k)]
    r <- cor(p[(k + 1):n], c[1:(n - k)])
    if (!is.na(r) && r > best_cor) {
      best_cor <- r
      best_lag <- k
    }
  }

  return(best_lag)
}

#' Shift a numeric vector forward by `lag` positions (pad with NA at the start)
#' This effectively delays the child series so it aligns with when its water
#' arrives at the parent outlet.
#'
#' @param x   Numeric vector
#' @param lag Integer >= 0
#' @return Shifted numeric vector (same length, leading NAs)
shift_series <- function(x, lag) {
  if (lag <= 0L) {
    return(x)
  }
  n <- length(x)
  c(rep(NA_real_, lag), x[1:(n - lag)])
}

#' Deaggregate discharge using lag-corrected mass balance.
#'
#' For each parent-child pair, the optimal routing lag is estimated via
#' cross-correlation. The child series is then shifted forward by that lag
#' before subtraction, so that the volumes are temporally aligned.
#'
#' @param dt_raw       data.table [nTime x nOutlets], columns named by catch_id
#' @param imm_children Named list: catch_id -> character vector of child ids
#' @param max_lag      Maximum lag to search in time steps (default 10)
#' @param verbose      Logical. Print lag estimates per pair.
#' @return A list with:
#'   - dt_residual: data.table same shape as dt_raw with lag-corrected residuals
#'   - lag_table:   data.table with columns parent_id, child_id, lag_steps
deaggregate_discharge_lagcorr <- function(dt_raw, imm_children,
                                          max_lag = 14L, verbose = TRUE) {
  res_dt <- copy(dt_raw)
  col_names <- names(dt_raw)
  factor <- 86400

  # Store estimated lags for diagnostics
  lag_records <- list()

  outlets_in_data <- intersect(names(imm_children), col_names)

  for (p_id in outlets_in_data) {
    child_ids <- imm_children[[p_id]]
    if (length(child_ids) == 0) next

    valid_children <- intersect(child_ids, col_names)
    if (length(valid_children) == 0) next

    parent_q <- dt_raw[[p_id]]
    tot_vol <- parent_q * factor

    # For each child, estimate lag and shift before summing
    shifted_child_vols <- matrix(0,
      nrow = nrow(dt_raw),
      ncol = length(valid_children)
    )

    for (i in seq_along(valid_children)) {
      ch_id <- valid_children[i]
      child_q <- dt_raw[[ch_id]]

      # Estimate lag (child -> parent travel time)
      lag_k <- estimate_lag(parent_q, child_q, max_lag = max_lag)

      if (verbose) {
        message(sprintf("  Lag %s -> %s: %d days", ch_id, p_id, lag_k))
      }

      lag_records[[length(lag_records) + 1]] <- data.table(
        parent_id = p_id, child_id = ch_id, lag_steps = lag_k
      )

      # Shift child series forward by lag, then convert to volume
      shifted_q <- shift_series(child_q, lag_k)
      shifted_child_vols[, i] <- shifted_q * factor
    }

    # Sum shifted children volumes (handle NAs from shifting)
    total_child_vol <- rowSums(shifted_child_vols, na.rm = TRUE)

    # Where all children are NA (first few rows due to shift), mark as NA
    all_na_rows <- rowSums(!is.na(shifted_child_vols)) == 0
    total_child_vol[all_na_rows] <- NA_real_

    # Compute lag-corrected residual discharge
    res_series <- (tot_vol - total_child_vol) / factor

    set(res_dt, j = p_id, value = res_series)
  }

  lag_table <- rbindlist(lag_records)

  if (verbose) {
    message(sprintf(
      "Lag correction applied to %d parent-child pairs. Mean lag: %.1f days.",
      nrow(lag_table), mean(lag_table$lag_steps)
    ))
  }

  return(list(dt_residual = res_dt, lag_table = lag_table))
}

# =============================================================================
# SECTION 3: Load Spatial Data and Nesting Structure
# =============================================================================
message("Reading catchments GeoPackage ...")
catchments_gpkg <- st_read(gpkg_path, quiet = TRUE)

stopifnot(all(c(
  "catch_id", "area_km2", "residual_area_km2",
  "immediate_nested_ids"
) %in% names(catchments_gpkg)))

# --- Parse nesting relationships ---
parse_ids <- function(x) {
  if (is.na(x) || trimws(x) == "" || trimws(x) == "NA") {
    return(character(0))
  }
  trimws(strsplit(as.character(x), ",")[[1]])
}

imm_children_list <- setNames(
  lapply(catchments_gpkg$immediate_nested_ids, parse_ids),
  as.character(catchments_gpkg$catch_id)
)

message(sprintf(
  "Loaded %d catchments (%d headwaters, %d nested).",
  nrow(catchments_gpkg),
  sum(lengths(imm_children_list) == 0),
  sum(lengths(imm_children_list) > 0)
))

# =============================================================================
# SECTION 4: Load Upstream Area (grid-based) and Build Lookups
# =============================================================================
message("Loading outlet and upstream area data ...")

# Load outlet coordinates
outhybas_raw <- outletopen(hydro_dir, outlet_name)
outhybas_raw$latlong <- paste(round(outhybas_raw$Var1, 4),
  round(outhybas_raw$Var2, 4),
  sep = " "
)
outhybas_raw$idlalo <- paste(outhybas_raw$idlo, outhybas_raw$idla, sep = " ")

# Load upstream area (km2 from NetCDF)
UpArea_full <- UpAopen(hydro_dir, uparea_file, outhybas_raw)

# Filter to European domain (EFAS network)
out_efas <- outletopen(hydro_dir, efas_file)
out_efas$latlong <- paste(round(out_efas$Var1, 4),
  round(out_efas$Var2, 4),
  sep = " "
)
outhybas_eu <- merge(out_efas, outhybas_raw, by = "latlong")

# Match to upstream area
matcat <- match(outhybas_eu$latlong, UpArea_full$latlong)
UpArea <- UpArea_full[matcat, ]

# --- Column selector for TSS files ---
header_ref <- fread(snow_file, nrows = 0, header = TRUE)
matcol <- match(UpArea_full$outlets, as.numeric(colnames(header_ref)))
matcol <- matcol[!is.na(matcol)]
cnames <- as.character(UpArea$outlets)

# --- Build area lookups from grid-based upstream area ---
# Using grid areas ensures mass-balance consistency with LISFLOOD
area_lookup <- setNames(UpArea$upa, as.character(UpArea$outlets))

# Residual area = parent upstream area - sum of children upstream areas
resarea_lookup <- setNames(
  sapply(as.character(catchments_gpkg$catch_id), function(p_id) {
    p_upa <- UpArea$upa[match(as.numeric(p_id), UpArea$outlets)]
    children <- imm_children_list[[p_id]]
    if (length(children) == 0) {
      return(p_upa)
    }
    children_upa <- sum(UpArea$upa[match(children, UpArea$outlets)], na.rm = TRUE)
    p_upa - children_upa
  }),
  as.character(catchments_gpkg$catch_id)
)

# --- Basemap setup for plotting ---
world <- ne_countries(scale = "medium", returnclass = "sf")
basemap <- st_transform(world, crs = 3035)

# Plot bounding box from outlet coordinates
cord.dec <- SpatialPoints(outhybas_raw[, c(2, 3)],
  proj4string = CRS("+proj=longlat")
)
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco <- cord.UTM@coords

# =============================================================================
# SECTION 2: Batch Deaggregation of All Variables
# =============================================================================
# Deaggregates each variable and saves as CSV with "_nested" suffix.
# Discharge (disWin) uses volume-based mass balance (Q = TRUE).
message("Starting batch deaggregation ...")

batch_vars <- c("disWin")
for (var in batch_vars) {
  is_discharge <- (var == "disWin")
  message("  Deaggregating: ", var)

  in_file <- paste0(base_dir, "data/tss/HERA_Histo/", var, "_1951_2020.csv")
  if (!file.exists(in_file)) {
    message("    File not found, skipping: ", in_file)
    next
  }

  Vari <- fread(in_file, header = TRUE)
  time_v <- Vari$V1
  Vari <- Vari[order(time_v), ]
  matcol_v <- match(UpArea$outlets, as.numeric(colnames(Vari)))
  matcol_v <- matcol_v[!is.na(matcol_v)]
  Vari <- Vari[, .SD, .SDcols = matcol_v]

  if (is_discharge) {
    # Use lag-corrected deaggregation for discharge
    message("    Using lag-corrected mass balance for discharge ...")
    lag_result <- deaggregate_discharge_lagcorr(
      Vari, imm_children_list,
      max_lag = 10L, verbose = TRUE
    )
    Vari_res <- lag_result$dt_residual

    # Save lag diagnostics table
    lag_out <- paste0(base_dir, "data/discharge_lag_estimates.csv")
    fwrite(lag_result$lag_table, lag_out)
    message("    Lag table saved to: ", lag_out)

    # Clamp remaining negatives (small residuals from imperfect lag estimate)
    Vari_res <- Vari_res[, lapply(.SD, function(x) pmax(0, x))]
  } else {
    Vari_res <- deaggregate_to_residual(Vari, imm_children_list,
      area_lookup, resarea_lookup,
      Q = FALSE
    )
    # Clip negative values for non-discharge variables
    Vari_res <- Vari_res[, lapply(.SD, function(x) pmax(0, x))]
  }

  # Prepend sorted timestamp column
  Vari_res <- data.table(time = time_v[order(time_v)], Vari_res)

  out_file <- paste0(base_dir, "data/", var, "_nested_1951_2020_v2.csv")
  fwrite(Vari_res, out_file,
    sep = ",", na = "-9999", row.names = FALSE,
    quote = FALSE, showProgress = TRUE, nThread = 2
  )

  rm(Vari, Vari_res)
  gc()
}

message("deaggregate_and_metrics.R completed successfully.")
