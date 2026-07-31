###############################################################################
# DISCHARGE VALIDATION - DIEGO-STYLE (Station-Specific)
#
# Computes per-station Spearman rho (with autocorrelation-corrected p-value)
# between HERA simulated discharge at the station pixel and observed daily
# discharge, then aggregates to catchment level using upstream-area-weighted
# averaging when multiple stations are matched to the same catchment.
#
# Temporal windows: Daily / 7-Day / 14-Day / Monthly (30-day rolling mean)
#
# Inputs:
#   output/station_catchment_matches.csv (from match_stations_catchments.R)
#   data/river_discharge/HERA_Val4_19502020.csv (simulated Q at station pixel)
#   data/river_discharge/HERA_CorStatv3_19512020.csv (corrected sim, overrides)
#   data/river_discharge/Q_19502020.csv (observed Q at station)
#
# Outputs:
#   output/discharge_station_metrics.csv  (per-station results, daily)
#   output/discharge_validation_all_windows.csv (all scenarios combined)
#   output/discharge_validation_summary.csv (per-catchment aggregated results)
#   output/discharge_diego/3.figures/Fig0_discharge_spatial_rho.png
#   output/discharge_diego/3.figures/Fig1_discharge_rho_violin.png
#   output/discharge_diego/3.figures/Fig2_discharge_scatter_climate.png
#   output/discharge_diego/3.figures/Fig4_discharge_stratified.png
#
# External rasters:
#   Climate : data/koppen_geiger_0p1.tif
#   DEM     : data/dem.nc
###############################################################################

# ===========================================================================
# SECTION 1: Library loading, paths, matching table
# ===========================================================================
library(sf)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(rnaturalearth)
library(terra)
library(exactextractr)


base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"

# Read matching table
matching_path <- file.path(base_dir, "output", "station_catchment_matches.csv")
if (!file.exists(matching_path)) {
    stop(
        "Matching table not found: ", matching_path,
        "\nRun match_stations_catchments.R first."
    )
}
matches <- read.csv(matching_path, stringsAsFactors = FALSE)
length(unique(matches$catch_id[which(matches$match_type=="relaxed")]))

cat("[1/5] Loading matching table and discharge data...\n")
cat(sprintf("  Matched pairs loaded: %d\n", nrow(matches)))

# Load station metadata for coordinates (needed for point-based map)
stations_meta_path <- file.path(base_dir, "Stations_ValidationF.csv")
stations_raw <- read.csv(stations_meta_path, stringsAsFactors = FALSE)

# ===========================================================================
# SECTION 2: Load simulated and observed discharge (station-specific)
# ===========================================================================

# --- Load observed discharge at station ---
Q_obs_path <- file.path(base_dir, "data", "river_discharge", "Q_19502020.csv")
if (!file.exists(Q_obs_path)) stop("Observed discharge file not found: ", Q_obs_path)

Q_data <- read.csv(Q_obs_path, header = FALSE)
time <- as.Date(as.numeric(Q_data$V2) - 1, origin = "0000-01-01")[-c(1, 2)]
Q_data <- Q_data[-1, -1] # remove first row (original header) and first column (row index)
Station_obs_IDs <- as.numeric(as.vector(t(Q_data[1, ])))[-1] # first remaining row = station IDs
Q_data <- Q_data[-1, ] # remove station ID row, keep data only
Q_data <- as.data.frame(lapply(Q_data, as.numeric))
rownames(Q_data) <- NULL


rm_obs <- (which(year(time) == 1950))

# --- Load simulated discharge (HERA_Val4 at station pixel) ---
Q_sim_path <- file.path(base_dir, "data", "river_discharge", "HERA_Val4_19502020.csv")
if (!file.exists(Q_sim_path)) stop("Simulated discharge file not found: ", Q_sim_path)

Q_sim <- read.csv(Q_sim_path, header = FALSE)
Q_sim <- Q_sim[-1, -1] # remove first row and first column
Station_sim_IDs <- as.numeric(as.vector(t(Q_sim[1, ]))) # first remaining row = station IDs
Q_sim <- Q_sim[-1, ] # remove station ID row
Q_sim <- as.data.frame(lapply(Q_sim, as.numeric))
rownames(Q_sim) <- NULL

# --- Load corrected stations (override HERA_Val4 for matching IDs) ---
Q_cor_path <- file.path(base_dir, "data", "river_discharge", "HERA_CorStatv3_19512020.csv")
has_corrected <- file.exists(Q_cor_path)
if (has_corrected) {
    Q_cor <- read.csv(Q_cor_path, header = FALSE)
    Q_cor <- Q_cor[-1, -1]
    Station_cor_IDs <- as.numeric(as.vector(t(Q_cor[1, ])))
    Q_cor <- Q_cor[-1, ]
    Q_cor <- as.data.frame(lapply(Q_cor, as.numeric))
    rownames(Q_cor) <- NULL
    cat(sprintf("  Corrected stations loaded: %d IDs\n", length(Station_cor_IDs)))
} else {
    Station_cor_IDs <- numeric(0)
    cat("  No corrected station file found, using HERA_Val4 only.\n")
}

# --- Date vectors ---
date_obs <- seq(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")
date_sim <- date_obs # HERA_Val4 starts 1950-01-01
date_cor <- seq(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")

cat(sprintf(
    "  Observed stations: %d | Simulated stations: %d\n",
    length(Station_obs_IDs), length(Station_sim_IDs)
))

# ===========================================================================
# SECTION 3: Per-station validation loop (with temporal aggregation)
# ===========================================================================
cat("[2/5] Computing per-station Spearman correlations...\n")

# Spearman with autocorrelation-corrected p-value
compute_spearman <- function(x, y) {
    valid_idx <- which(!is.na(x) & !is.na(y))
    n <- length(valid_idx)
    if (n < 365) {
        return(list(
            rho = NA, p_val = NA, n_obs = n, n_eff = NA,
            mean_obs = NA, mean_sim = NA, status = "Insufficient Data"
        ))
    }
    xv <- x[valid_idx]
    yv <- y[valid_idx]
    if (sd(xv) == 0 || sd(yv) == 0) {
        return(list(
            rho = NA, p_val = NA, n_obs = n, n_eff = NA,
            mean_obs = mean(yv), mean_sim = mean(xv), status = "No Variance"
        ))
    }
    rho <- cor(xv, yv, method = "spearman", use = "complete.obs")
    r1x <- cor(xv[-n], xv[-1], use = "complete.obs")
    r1y <- cor(yv[-n], yv[-1], use = "complete.obs")
    r1 <- max(0, mean(c(r1x, r1y), na.rm = TRUE))
    n_eff <- max(3, floor(n * (1 - r1) / (1 + r1)))
    t_stat <- rho * sqrt((n_eff - 2) / (1 - rho^2))
    p_val <- 2 * pt(abs(t_stat), df = n_eff - 2, lower.tail = FALSE)
    status <- if (is.na(p_val) || p_val >= 0.05) "Not Significant" else "Significant"
    list(
        rho = rho, p_val = p_val, n_obs = n, n_eff = n_eff,
        mean_obs = mean(yv), mean_sim = mean(xv), status = status
    )
}

compute_spearman <- function(mat_obs, mat_mod, ids) {
  rbindlist(lapply(ids, function(id) {
    x <- mat_mod[[id]]
    y <- mat_obs[[id]]
    valid_idx <- which(!is.na(x) & !is.na(y))
    n=length(valid_idx)
    
    if (length(valid_idx) < 10 || all(is.na(y)) || all(is.na(x))) {
      return(list(
        catch_id = id, rho = NA_real_, p_val = NA_real_,
        n_eff = NA_real_, status = "No Data"
      ))
    }
    
    xv <- x[valid_idx]
    yv <- y[valid_idx]
    n <- length(xv)
    
    if (sd(xv) == 0 || sd(yv) == 0) {
      return(list(
        catch_id = id, rho = NA_real_, p_val = NA_real_,
        n_eff = NA_real_, status = "No Variance"
      ))
    }
    
    rho <- cor(xv, yv, method = "spearman")
    
    r1x <- cor(xv[-n], xv[-1], use = "complete.obs")
    r1y <- cor(yv[-n], yv[-1], use = "complete.obs")
    r1 <- mean(c(r1x, r1y), na.rm = TRUE)
    r1 <- max(0, r1)
    
    n_eff <- max(3, floor(n * (1 - r1) / (1 + r1)))
    t_stat <- rho * sqrt((n_eff - 2) / (1 - rho^2))
    p_corr <- 2 * pt(abs(t_stat), df = n_eff - 2, lower.tail = FALSE)
    
    status <- if (is.na(p_corr) || p_corr >= 0.05) "Not Significant" else "Significant"
    
    list(
      catch_id = id, rho = round(rho, 2), p_val = p_corr, n_obs = n,
      mean_obs = mean(yv), mean_sim = mean(xv), n_eff = n_eff, status = status
    )
  }))
}
# ---------------------------------------------------------------------------
# Function to compute per-station correlations with optional smoothing
# k = smoothing window (1 = daily, 7 = 7-day, 14 = 14-day, 30 = monthly)
# ---------------------------------------------------------------------------
compute_all_stations <- function(matches, Q_sim, Station_sim_IDs, Q_cor,
                                 Station_cor_IDs, Q_data, Station_obs_IDs,
                                 date_sim, date_obs, date_cor,
                                 rm_obs, has_corrected, k = 1) {
    station_results <- vector("list", nrow(matches))

    for (i in seq_len(nrow(matches))) {
        cid <- as.character(matches$catch_id[i])
        sid <- as.numeric(matches$station_id[i])
        mtype <- matches$match_type[i]
        station_upa <- matches$station_upa[i]
        catchment_upa <- matches$catchment_upa[i]

        # Determine simulated source: corrected overrides HERA_Val4
        use_corrected <- has_corrected && (sid %in% Station_cor_IDs)

        # --- Find station in simulated data ---
        if (use_corrected) {
            sim_col_idx <- which(Station_cor_IDs == sid)
        } else {
            sim_col_idx <- which(Station_sim_IDs == sid)
        }

        if (length(sim_col_idx) == 0) {
            station_results[[i]] <- data.frame(
                catch_id = cid, station_id = sid, match_type = mtype,
                station_upa = station_upa, catchment_upa = catchment_upa,
                spearman_rho = NA, p_value = NA, n_obs = NA, n_eff = NA,
                mean_obs = NA, mean_sim = NA, status = "Station Not In Sim",
                stringsAsFactors = FALSE
            )
            next
        }

        # --- Find station in observed data ---
        obs_col_idx <- which(Station_obs_IDs == sid)

        if (length(obs_col_idx) == 0) {
            station_results[[i]] <- data.frame(
                catch_id = cid, station_id = sid, match_type = mtype,
                station_upa = station_upa, catchment_upa = catchment_upa,
                spearman_rho = NA, p_value = NA, n_obs = NA, n_eff = NA,
                mean_obs = NA, mean_sim = NA, status = "Station Not In Obs",
                stringsAsFactors = FALSE
            )
            next
        }

        # --- Extract time series ---
        if (use_corrected) {
            sim_vec <- Q_cor[, sim_col_idx[1]]
            sim_dates <- date_cor
        } else {
            sim_vec <- Q_sim[, sim_col_idx[1]]
            sim_dates <- date_sim
        }

        obs_vec <- Q_data[-rm_obs, obs_col_idx[1] + 1]
        obs_dates <- date_obs

        # Build data.tables and merge on date
        dt_sim <- data.table(date = sim_dates, sim = as.numeric(sim_vec))
        dt_obs <- data.table(date = obs_dates, obs = as.numeric(obs_vec))
        paired <- merge(dt_sim, dt_obs, by = "date", all = FALSE)
        pairedX <- paired[!is.na(sim) & !is.na(obs)]
        n_paired <- nrow(pairedX)

        if (n_paired < 365) {
            station_results[[i]] <- data.frame(
                catch_id = cid, station_id = sid, match_type = mtype,
                station_upa = station_upa, catchment_upa = catchment_upa,
                spearman_rho = NA, p_value = NA, n_obs = n_paired, n_eff = NA,
                mean_obs = NA, mean_sim = NA, status = "Insufficient Data",
                stringsAsFactors = FALSE
            )
            next
        }

        # Apply rolling mean smoothing if k > 1
        sim_final <- data.table(paired$sim)
        obs_final <- data.table(paired$obs) 
        if (k > 1) {
            obs_final[, agg := (seq_len(.N) - 1) %/% k]
            sim_final[, agg := (seq_len(.N) - 1) %/% k]
            obs_final <- obs_final[, lapply(.SD, mean, na.rm = TRUE), by = agg]
            sim_final <- sim_final[, lapply(.SD, mean, na.rm = TRUE), by = agg]
            # sim_final <- data.table::frollmean(sim_final, n = k, align = "right", na.rm = TRUE)
            # obs_final <- data.table::frollmean(obs_final, n = k, align = "right", na.rm = TRUE)
        }

        # Compute Spearman
        daily_cols <- setdiff(names(sim_final), "agg")
        sp <- compute_spearman(sim_final, obs_final,daily_cols)
        sp$catch_id=cid

        station_results[[i]] <- data.frame(
            catch_id = cid, station_id = sid, match_type = mtype,
            station_upa = station_upa, catchment_upa = catchment_upa,
            spearman_rho = sp$rho, p_value = sp$p_val,
            n_obs = sp$n_obs, n_eff = sp$n_eff,
            mean_obs = sp$mean_obs, mean_sim = sp$mean_sim,
            status = sp$status,
            stringsAsFactors = FALSE
        )
    }

    # Combine per-station results
    result_df <- do.call(rbind, station_results)
    rownames(result_df) <- NULL
    result_df
}

# --- Run all temporal windows ---
cat("  Computing Daily...\n")
station_df <- compute_all_stations(
    matches, Q_sim, Station_sim_IDs, Q_cor, Station_cor_IDs,
    Q_data, Station_obs_IDs, date_sim, date_obs, date_cor,
    rm_obs, has_corrected,
    k = 1
)

cat("  Computing 7-day rolling mean...\n")
station_df_7d <- compute_all_stations(
    matches, Q_sim, Station_sim_IDs, Q_cor, Station_cor_IDs,
    Q_data, Station_obs_IDs, date_sim, date_obs, date_cor,
    rm_obs, has_corrected,
    k = 7
)

cat("  Computing 15-day rolling mean...\n")
station_df_15d <- compute_all_stations(
    matches, Q_sim, Station_sim_IDs, Q_cor, Station_cor_IDs,
    Q_data, Station_obs_IDs, date_sim, date_obs, date_cor,
    rm_obs, has_corrected,
    k = 15
)

cat("  Computing 30-day rolling mean...\n")
station_df_30d <- compute_all_stations(
    matches, Q_sim, Station_sim_IDs, Q_cor, Station_cor_IDs,
    Q_data, Station_obs_IDs, date_sim, date_obs, date_cor,
    rm_obs, has_corrected,
    k = 30
)

# Add scenario labels
station_df$scenario <- "Daily"
station_df_7d$scenario <- "7-Day"
station_df_15d$scenario <- "15-Day"
station_df_30d$scenario <- "Monthly"

# Combine all windows
station_df_all <- rbind(station_df, station_df_7d, station_df_15d, station_df_30d)
station_df_all$scenario <- factor(station_df_all$scenario,
    levels = c("Daily", "7-Day", "15-Day", "Monthly")
)

median(station_df_15d$spearman_rho, na.rm = TRUE)
cat(sprintf(
    "  Stations with valid rho (daily): %d / %d\n",
    sum(!is.na(station_df$spearman_rho)), nrow(station_df)
))

# ===========================================================================
# SECTION 4: Aggregate to catchment level (weighted average) - all windows
# ===========================================================================
cat("[3/5] Aggregating to catchment level...\n")

# Function to aggregate station results to catchment level
aggregate_to_catchment <- function(sdf, match) {
    sdf$upa_ratio <- sdf$station_upa / sdf$catchment_upa

    # if (match=="strict"){
    #   sdf=sdf[which(sdf$match_type="strict"),]
    # }
    # if (match=="relaxed"){
    #   sdf=sdf[which(sdf$match_type="relaxed"),]
    # }
    catchment_summary <- sdf %>%
        group_by(catch_id) %>%
        summarise(
            n_stations = n(),
            match_type = paste(unique(match_type), collapse = "/"),
            spearman_rho = {
                valid <- !is.na(spearman_rho)
                if (sum(valid) == 0) {
                    NA_real_
                } else if (sum(valid) == 1) {
                    spearman_rho[valid]
                } else {
                    weighted.mean(spearman_rho[valid], w = upa_ratio[valid], na.rm = TRUE)
                }
            },
            p_value = {
                valid <- !is.na(p_value)
                if (sum(valid) == 0) {
                    NA_real_
                } else {
                    min(p_value[valid], na.rm = TRUE)
                }
            },
            n_obs = sum(n_obs, na.rm = TRUE),
            n_eff = {
                valid <- !is.na(n_eff)
                if (sum(valid) == 0) {
                    NA_real_
                } else {
                    sum(n_eff[valid])
                }
            },
            mean_obs = {
                valid <- !is.na(mean_obs)
                if (sum(valid) == 0) {
                    NA_real_
                } else if (sum(valid) == 1) {
                    mean_obs[valid]
                } else {
                    weighted.mean(mean_obs[valid], w = upa_ratio[valid], na.rm = TRUE)
                }
            },
            mean_sim = {
                valid <- !is.na(mean_sim)
                if (sum(valid) == 0) {
                    NA_real_
                } else if (sum(valid) == 1) {
                    mean_sim[valid]
                } else {
                    weighted.mean(mean_sim[valid], w = upa_ratio[valid], na.rm = TRUE)
                }
            },
            status = {
                valid <- !is.na(spearman_rho)
                if (sum(valid) == 0) {
                    "No Valid Stations"
                } else if (!is.na(min(p_value[valid], na.rm = TRUE)) &&
                    min(p_value[valid], na.rm = TRUE) < 0.05) {
                    "Significant"
                } else {
                    "Not Significant"
                }
            },
            .groups = "drop"
        ) %>%
        as.data.frame()

    catchment_summary
}


# Aggregate each window
catchment_summary <- aggregate_to_catchment(station_df)
catchment_summary_7d <- aggregate_to_catchment(station_df_7d)
catchment_summary_15d <- aggregate_to_catchment(station_df_15d)
catchment_summary_30d <- aggregate_to_catchment(station_df_30d)

# Label and combine
catchment_summary$scenario <- "Daily"
catchment_summary_7d$scenario <- "7-Day"
catchment_summary_15d$scenario <- "15-Day"
catchment_summary_30d$scenario <- "Monthly"

catchment_summary_all <- rbind(
    catchment_summary, catchment_summary_7d,
    catchment_summary_14d, catchment_summary_30d
)
catchment_summary_all$scenario <- factor(catchment_summary_all$scenario,
    levels = c("Daily", "7-Day", "15-Day", "Monthly")
)

cat(sprintf(
    "  Catchments with valid rho (daily): %d / %d\n",
    sum(!is.na(catchment_summary$spearman_rho)), nrow(catchment_summary)
))
cat(sprintf(
    "  Multi-station catchments: %d\n",
    sum(catchment_summary$n_stations > 1)
))

median(catchment_summary$spearman_rho[which(catchment_summary$match_type=="strict")])
length(which(catchment_summary$spearman_rho[which(catchment_summary$match_type=="strict")]>0.5))/
  length(catchment_summary$spearman_rho[which(catchment_summary$match_type=="strict")])


median(catchment_summary$spearman_rho[which(catchment_summary$match_type=="relaxed")])
length(which(catchment_summary$spearman_rho[which(catchment_summary$match_type=="relaxed")]>0.5))/
  length(catchment_summary$spearman_rho[which(catchment_summary$match_type=="relaxed")])
# ===========================================================================
# SECTION 5: Export CSVs
# ===========================================================================
cat("[4/5] Exporting results...\n")

output_dir <- file.path(base_dir, "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Per-station results (daily only, backward compatible)
write.csv(station_df,
    file = file.path(output_dir, "discharge_station_metrics.csv"),
    row.names = FALSE, fileEncoding = "UTF-8"
)

# All windows combined (station-level)
write.csv(station_df_all,
    file = file.path(output_dir, "discharge_validation_all_windows.csv"),
    row.names = FALSE, fileEncoding = "UTF-8"
)

# Per-catchment aggregated results (daily)
write.csv(catchment_summary,
    file = file.path(output_dir, "discharge_validation_summary.csv"),
    row.names = FALSE, fileEncoding = "UTF-8"
)

# Print summary statistics
valid_rows <- catchment_summary[!is.na(catchment_summary$spearman_rho), ]
n_validated <- nrow(valid_rows)
median_rho <- median(valid_rows$spearman_rho, na.rm = TRUE)
pct_ge05 <- round(mean(valid_rows$spearman_rho >= 0.5, na.rm = TRUE) * 100, 1)

cat(sprintf("\n  === VALIDATION SUMMARY ===\n"))
cat(sprintf("  Catchments validated: %d / %d\n", n_validated, nrow(catchment_summary)))
cat(sprintf("  Median Spearman rho: %.3f\n", median_rho))
cat(sprintf("  Percent with rho >= 0.5: %.1f%%\n", pct_ge05))
cat(sprintf(
    "  Per-station output: %s\n",
    file.path(output_dir, "discharge_station_metrics.csv")
))
cat(sprintf(
    "  Per-catchment output: %s\n",
    file.path(output_dir, "discharge_validation_summary.csv")
))

# Summary by scenario
cat("\n  === PER-SCENARIO SUMMARY ===\n")
for (sc in levels(catchment_summary_all$scenario)) {
    sub <- catchment_summary_all[catchment_summary_all$scenario == sc, ]
    valid_sub <- sub[!is.na(sub$spearman_rho), ]
    cat(sprintf(
        "  %s: median rho = %.3f, pct >= 0.5: %.1f%%\n",
        sc, median(valid_sub$spearman_rho, na.rm = TRUE),
        round(mean(valid_sub$spearman_rho >= 0.5, na.rm = TRUE) * 100, 1)
    ))
}

# ===========================================================================
# SECTION 6: Diego-style diagnostic figures (matches SWE script exactly)
# ===========================================================================
cat("[5/5] Generating diagnostic figures...\n")

# --- Create output dir, load classification data ---
path_out <- file.path(base_dir, "output", "discharge_diego", "3.figures")
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

# Load catchment polygons
file_shp <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
path_clim <- file.path(base_dir, "data", "koppen_geiger_0p1.tif")
path_dem <- file.path(base_dir, "data", "dem.nc")

shp <- st_read(file_shp, quiet = TRUE)

norm_id <- function(x) {
    x <- as.character(x)
    suppressWarnings(ifelse(grepl("^\\d+$", x), as.character(as.integer(x)), x))
}
shp$join_id <- norm_id(shp$catch_id)

clim_lookup <- c(
    setNames(rep("Tropical", 3), as.character(1:3)),
    setNames(rep("Arid", 4), as.character(4:7)),
    setNames(rep("Temperate", 9), as.character(8:16)),
    setNames(rep("Cold", 12), as.character(17:28)),
    setNames(rep("Polar", 2), as.character(29:30))
)

r_clim <- terra::rast(path_clim)
shp$clim_class <- exactextractr::exact_extract(
    r_clim, sf::st_transform(shp, terra::crs(r_clim)),
    fun = function(values, coverage_fractions) {
        if (all(is.na(values))) {
            return(NA_character_)
        }
        maj <- as.character(names(which.max(table(values[!is.na(values)]))))
        unname(clim_lookup[maj])
    }
)

# Extract mean elevation from DEM
r_dem <- terra::rast(path_dem)
shp$elev_m <- exactextractr::exact_extract(
    r_dem, sf::st_transform(shp, terra::crs(r_dem)),
    fun = "mean"
)

shp <- shp %>% mutate(
    clim_class = factor(clim_class, levels = c("Polar", "Cold", "Temperate", "Arid", "Tropical")),
    elev_class = cut(elev_m,
        breaks = c(-Inf, 200, 500, 1000, 2000, Inf),
        labels = c("< 200 m", "200-500 m", "500-1000 m", "1000-2000 m", "> 2000 m")
    ),
    area_class = cut(area_km2,
        breaks = quantile(area_km2, probs = c(0, .25, .5, .75, 1), na.rm = TRUE),
        labels = c("Q1 (smallest)", "Q2", "Q3", "Q4 (largest)"), include.lowest = TRUE
    )
)
meta_lut <- st_drop_geometry(shp) %>%
    dplyr::select(join_id, clim_class, area_km2, area_class, elev_m, elev_class)

# Build long rho table (all_meta) matching SWE format
strip_id <- function(df) {
    df <- as.data.frame(df)
    df$join_id <- norm_id(df$catch_id)
    df$rho <- df$spearman_rho
    df
}
attach_meta <- function(stats_df, label) {
    strip_id(stats_df) %>%
        left_join(meta_lut, by = "join_id") %>%
        filter(!is.na(rho)) %>%
        mutate(scenario = label)
}
all_meta <- bind_rows(
    attach_meta(catchment_summary, "Daily"),
    attach_meta(catchment_summary_7d, "7-Day"),
    attach_meta(catchment_summary_14d, "14-Day"),
    attach_meta(catchment_summary_30d, "Monthly")
) %>% mutate(scenario = factor(scenario, levels = c("Daily", "7-Day", "14-Day", "Monthly")))

# Theme + palettes (identical to SWE)
theme_nature <- function(base_size = 9) {
    theme_bw(base_size = base_size) %+replace%
        theme(
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
            strip.background = element_rect(fill = "grey96", colour = "grey60"),
            strip.text = element_text(face = "bold", size = base_size),
            plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
            legend.key = element_blank()
        )
}
clim_pal <- c(
    "Polar" = "#7fbfff", "Cold" = "#2166ac", "Temperate" = "#4dac26",
    "Arid" = "#d6604d", "Tropical" = "#8e0152"
)
scen_pal <- c(
    "Daily" = "#bdbdbd", "7-Day" = "#74c4e4",
    "14-Day" = "#2c9e4b", "Monthly" = "#1a3f7a"
)

# --- Fig0: Point-based spatial rho maps (4-panel, improved) ---
message("Fig 0: spatial rho maps...")
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
basemap <- sf::st_transform(world, crs = 3035)
palet <- c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = TRUE, fixup = TRUE))
tsize <- 12
osize <- 12

station_coords <- stations_raw[, c("V1", "Var1", "Var2", "upa")]
names(station_coords)[1] <- "station_id"
station_coords$station_id <- as.numeric(station_coords$station_id)

build_point_map <- function(sdf, ttl) {
    plot_df <- merge(sdf[!is.na(sdf$spearman_rho), ],
        station_coords[, c("station_id", "Var1", "Var2")],
        by = "station_id", all.x = TRUE
    )
    plot_df <- plot_df[!is.na(plot_df$Var1) & !is.na(plot_df$Var2), ]
    plot_df$upa <- plot_df$station_upa
    parpl <- st_as_sf(plot_df, coords = c("Var1", "Var2"), crs = 4326)
    parpl <- st_transform(parpl, crs = 3035)
    nco <- st_coordinates(parpl)
    ggplot(basemap) +
        geom_sf(fill = "gray95", color = NA) +
        geom_sf(
            data = parpl, aes(fill = spearman_rho, size = upa),
            color = "transparent", alpha = 0.7, shape = 21, stroke = 0
        ) +
        geom_sf(
            data = parpl, aes(size = upa),
            col = "grey20", alpha = 1, stroke = 0.05, shape = 1
        ) +
        geom_sf(data = basemap, fill = NA, color = "grey20") +
        scale_size(range = c(1, 3), trans = "sqrt") +
        scale_fill_gradientn(
            colors = palet, n.breaks = 5, oob = scales::squish,
            limits = c(0, 1), name = "Spearman \u03c1"
        ) +
        coord_sf(
            xlim = c(min(nco[, 1]), max(nco[, 1])),
            ylim = c(min(nco[, 2]), max(nco[, 2]))
        ) +
        labs(x = "Longitude", y = "Latitude", subtitle = ttl) +
        guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12), size = "none") +
        theme(
            axis.title = element_text(size = tsize),
            panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
            panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
            legend.title = element_text(size = tsize),
            legend.text = element_text(size = osize),
            legend.position = "right",
            panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
            panel.grid.minor = element_line(colour = "grey90"),
            legend.key = element_rect(fill = "transparent", colour = "transparent"),
            legend.key.size = unit(0.8, "cm")
        )
}
fig0 <- (build_point_map(station_df, "a) Daily") +
    build_point_map(station_df_7d, "b) 7-Day Moving Mean")) /
    (build_point_map(station_df_14d, "c) 14-Day Moving Mean") +
        build_point_map(station_df_30d, "d) Monthly")) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
fig0 <- fig0 + plot_annotation(
    title = "Discharge cross-comparison: HERA vs Observed (Spearman \u03c1)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
)
ggsave(file.path(path_out, "Fig0_discharge_spatial_rho.png"), fig0,
    width = 34, height = 26, units = "cm", dpi = 300, bg = "white"
)

# 7. FIG 1 | rho VIOLIN + BOX (identical to SWE) -------------------
message("Fig 1: rho violin...")
fig1_ann <- all_meta %>%
    group_by(scenario) %>%
    summarise(
        med = round(median(rho, na.rm = TRUE), 2),
        pct = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 0), .groups = "drop"
    )
fig1 <- ggplot(all_meta, aes(scenario, rho, fill = scenario)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", colour = "grey40") +
    geom_boxplot(width = 0.12, outlier.shape = NA, colour = "grey25", fill = "white", fatten = 2) +
    geom_text(
        data = fig1_ann, aes(scenario, med, label = paste0("md=", med)),
        vjust = -0.6, size = 2.6, fontface = "bold", inherit.aes = FALSE
    ) +
    geom_text(
        data = fig1_ann, aes(scenario, 0.97, label = paste0(pct, "% \u2265 0.5")),
        size = 3, colour = "grey30", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = scen_pal, guide = "none") +
    scale_y_continuous(limits = c(-0.35, 1.02), breaks = seq(-0.25, 1, 0.25)) +
    labs(
        title = "Fig. 1 | Spearman \u03c1 across temporal aggregations (Discharge)",
        x = NULL, y = "Spearman \u03c1"
    ) +
    theme_nature()
ggsave(file.path(path_out, "Fig1_discharge_rho_violin.png"), fig1,
    width = 16, height = 12, units = "cm", dpi = 300, bg = "white"
)

# 10. FIG 4 | STRATIFIED BOXPLOTS (identical to SWE) ---------------
message("Fig 4: stratified diagnostics...")
make_strat_box <- function(data, x_var, fill_var, fill_pal, fill_name, x_lab, tag, rotate = FALSE) {
    p <- ggplot(data, aes(x = .data[[x_var]], y = rho, fill = .data[[fill_var]])) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
        geom_boxplot(
            position = position_dodge(width = 0.8), width = 0.68,
            outlier.size = 0.2, outlier.alpha = 0.25, outlier.shape = 1, linewidth = 0.35
        ) +
        scale_fill_manual(values = fill_pal, name = fill_name) +
        scale_y_continuous(limits = c(-0.3, 1.05), breaks = seq(0, 1, 0.25)) +
        labs(tag = tag, x = x_lab, y = "Spearman \u03c1") +
        theme_nature(8.5) +
        theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))
    if (rotate) p <- p + theme(axis.text.x = element_text(angle = 25, hjust = 1))
    p
}
p4a <- make_strat_box(filter(all_meta, !is.na(clim_class)), "scenario", "clim_class", clim_pal, "Climate class", NULL, "a")
p4b <- make_strat_box(filter(all_meta, !is.na(area_class)), "area_class", "scenario", scen_pal, "Aggregation", "Catchment area quartile", "b", TRUE)
p4c <- make_strat_box(filter(all_meta, !is.na(elev_class)), "elev_class", "scenario", scen_pal, "Aggregation", "Elevation class", "c", TRUE)
fig4 <- (p4a / (p4b + p4c)) +
    plot_annotation(title = "Fig. 4 | Discharge agreement stratified by catchment characteristics")
ggsave(file.path(path_out, "Fig4_discharge_stratified.png"), fig4,
    width = 20, height = 20, units = "cm", dpi = 300, bg = "white"
)

# 11. FIG 5 | RIDGE DENSITY (identical to SWE) ---------------------
if (requireNamespace("ggridges", quietly = TRUE)) {
    message("Fig 5: ridge density...")
    fig5 <- ggplot(
        filter(all_meta, !is.na(clim_class)),
        aes(x = rho, y = forcats::fct_rev(clim_class), fill = clim_class, colour = clim_class)
    ) +
        geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55") +
        geom_vline(xintercept = 0, colour = "grey70") +
        ggridges::geom_density_ridges(alpha = 0.55, scale = 1.15, bandwidth = 0.04, linewidth = 0.35) +
        stat_summary(fun = median, geom = "point", shape = 124, size = 4, colour = "grey15") +
        scale_fill_manual(values = clim_pal, guide = "none") +
        scale_colour_manual(values = clim_pal, guide = "none") +
        scale_x_continuous(limits = c(-0.3, 1.0), breaks = seq(-0.25, 1, 0.25)) +
        facet_wrap(~scenario, ncol = 4) +
        labs(
            title = "Fig. 5 | Distribution of Spearman \u03c1 per climate class",
            x = "Spearman \u03c1", y = NULL
        ) +
        theme_nature()
    ggsave(file.path(path_out, "Fig5_discharge_ridge_climate.png"), fig5,
        width = 22, height = 10, units = "cm", dpi = 300, bg = "white"
    )
} else {
    message("Fig 5 skipped (install 'ggridges' to enable).")
}

# 12. FIG 6 | AGGREGATION-GAIN HEATMAP (identical to SWE) ----------
message("Fig 6: aggregation-gain heatmap...")
compute_gain <- function(data, group_var, group_label) {
    long <- data %>%
        filter(!is.na(.data[[group_var]])) %>%
        group_by(stratum = as.character(.data[[group_var]]), scenario) %>%
        summarise(med_rho = median(rho, na.rm = TRUE), .groups = "drop")
    baseline <- long %>%
        filter(scenario == "Daily") %>%
        dplyr::select(stratum, baseline = med_rho)
    long %>%
        filter(scenario != "Daily") %>%
        left_join(baseline, by = "stratum") %>%
        mutate(delta_rho = med_rho - baseline, stratum_type = group_label) %>%
        dplyr::select(stratum, stratum_type, aggregation = scenario, delta_rho)
}
gain_dt <- bind_rows(
    compute_gain(all_meta, "clim_class", "Climate class"),
    compute_gain(all_meta, "area_class", "Catchment area"),
    compute_gain(all_meta, "elev_class", "Elevation")
) %>% mutate(
    aggregation = factor(aggregation, levels = c("7-Day", "14-Day", "Monthly")),
    stratum_type = factor(stratum_type, levels = c("Climate class", "Catchment area", "Elevation")),
    stratum = factor(stratum, levels = c(
        "Polar", "Cold", "Temperate", "Arid", "Tropical",
        "Q1 (smallest)", "Q2", "Q3", "Q4 (largest)",
        "< 200 m", "200-500 m", "500-1000 m", "1000-2000 m", "> 2000 m"
    ))
)
delta_max <- max(abs(gain_dt$delta_rho), na.rm = TRUE)
fig6 <- ggplot(gain_dt, aes(aggregation, stratum, fill = delta_rho)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%+.2f", delta_rho)), size = 2.8, fontface = "bold") +
    scale_fill_gradientn(
        colours = c("#c2523c", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#4393c3", "#1a3f7a"),
        values = scales::rescale(c(-delta_max, -delta_max / 2, -0.02, 0, 0.02, delta_max / 2, delta_max)),
        limits = c(-delta_max, delta_max), oob = squish, name = "\u0394\u03c1 vs daily"
    ) +
    facet_grid(stratum_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(
        title = "Fig. 6 | Aggregation gain (\u0394\u03c1 vs daily) across strata",
        x = "Temporal aggregation", y = NULL
    ) +
    theme_nature() +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"))
ggsave(file.path(path_out, "Fig6_discharge_aggregation_gain.png"), fig6,
    width = 18, height = 16, units = "cm", dpi = 300, bg = "white"
)

# 13. TABLES (identical to SWE) ------------------------------------
message("Writing stratified tables...")
tbl <- function(group_var) {
    all_meta %>%
        filter(!is.na(.data[[group_var]])) %>%
        group_by(scenario, stratum = .data[[group_var]]) %>%
        summarise(
            n = dplyr::n(), median_rho = round(median(rho, na.rm = TRUE), 3),
            pct_ge_0.5 = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 1), .groups = "drop"
        )
}
data.table::fwrite(tbl("clim_class"), file.path(path_out, "Table_rho_by_climate.csv"))
data.table::fwrite(tbl("area_class"), file.path(path_out, "Table_rho_by_area.csv"))
data.table::fwrite(tbl("elev_class"), file.path(path_out, "Table_rho_by_elevation.csv"))

message("Done. Figures + tables in: ", path_out)
