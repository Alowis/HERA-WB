# =============================================================================
# Homogenize discharge data: create paired obs/sim monthly CSVs
# similar to the GLEAM/LISFLOOD AET homogenized files.
#
# Output (in output/discharge_diego/1.homogenized/):
#   - obs_monthly_strict.csv   (observed Q, strict-match stations only)
#   - obs_monthly_relaxed.csv  (observed Q, all matched stations)
#   - sim_monthly.csv          (simulated Q at matched catchments)
#
# Format: date column (YYYY-MM) + one column per catch_id (matched)
# =============================================================================

library(data.table)
library(lubridate)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
dis_dir <- file.path(base_dir, "data", "river_discharge")
match_path <- file.path(base_dir, "output", "station_catchment_matches.csv")
out_dir <- file.path(base_dir, "output", "discharge_diego", "1.homogenized")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load station-catchment matches -------------------------------------------
cat("Loading station-catchment matches...\n")
matches <- fread(match_path)
cat("  Total matches:", nrow(matches), "\n")
cat("  Strict:", sum(matches$match_type == "strict"), "\n")
cat("  Relaxed:", sum(matches$match_type == "relaxed"), "\n")

# --- Load observed discharge --------------------------------------------------
cat("Loading observed discharge...\n")
Q_obs_raw <- read.csv(file.path(dis_dir, "Q_19502020.csv"), header = FALSE)

# First row after header = station IDs, data starts after
time_obs <- as.Date(as.numeric(Q_obs_raw[-1, 2]) - 1, origin = "0000-01-01")
Q_obs_raw <- Q_obs_raw[-1, -1] # remove header row + row index column
Station_obs_IDs <- as.numeric(as.vector(t(Q_obs_raw[1, ])))
Q_obs_raw <- Q_obs_raw[-1, ]
Q_obs_raw <- as.data.frame(lapply(Q_obs_raw, as.numeric))
rownames(Q_obs_raw) <- NULL

# Remove 1950
keep_idx <- which(year(time_obs) >= 1951)
time_obs <- time_obs[keep_idx]
Q_obs_raw <- Q_obs_raw[keep_idx, ]

cat("  Observed: ", ncol(Q_obs_raw), "stations,", nrow(Q_obs_raw), "days\n")

# --- Load simulated discharge (HERA_Val4 + corrected) -------------------------
cat("Loading simulated discharge...\n")
Q_sim_raw <- read.csv(file.path(dis_dir, "HERA_Val4_19502020.csv"), header = FALSE)
Q_sim_raw <- Q_sim_raw[-1, -1]
Station_sim_IDs <- as.numeric(as.vector(t(Q_sim_raw[1, ])))
Q_sim_raw <- Q_sim_raw[-1, ]
Q_sim_raw <- as.data.frame(lapply(Q_sim_raw, as.numeric))
rownames(Q_sim_raw) <- NULL

date_sim <- seq(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")
# Trim to match if needed
if (nrow(Q_sim_raw) > length(date_sim)) {
    Q_sim_raw <- Q_sim_raw[1:length(date_sim), ]
}

# Load corrected stations (override)
Q_cor_path <- file.path(dis_dir, "HERA_CorStatv3_19512020.csv")
has_corrected <- file.exists(Q_cor_path)
if (has_corrected) {
    Q_cor_raw <- read.csv(Q_cor_path, header = FALSE)
    Q_cor_raw <- Q_cor_raw[-1, -1]
    Station_cor_IDs <- as.numeric(as.vector(t(Q_cor_raw[1, ])))
    Q_cor_raw <- Q_cor_raw[-1, ]
    Q_cor_raw <- as.data.frame(lapply(Q_cor_raw, as.numeric))
    date_cor <- seq(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")
    cat("  Corrected stations available:", length(Station_cor_IDs), "\n")
}

# --- Extract paired time series per match -------------------------------------
cat("Extracting paired daily time series...\n")

# Common date range
date_common <- seq(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")

extract_paired <- function(match_subset) {
    n_matches <- nrow(match_subset)

    # Initialize daily data.tables
    obs_daily <- data.table(date = date_common)
    sim_daily <- data.table(date = date_common)

    for (i in seq_len(n_matches)) {
        sid <- match_subset$station_id[i]
        cid <- as.character(match_subset$catch_id[i])

        # --- Observed ---
        obs_col_idx <- which(Station_obs_IDs == sid)
        if (length(obs_col_idx) == 0) next

        obs_vec <- Q_obs_raw[, obs_col_idx[1]]
        # Align to common dates (obs may have different length)
        if (length(obs_vec) == length(time_obs)) {
            obs_aligned <- rep(NA_real_, length(date_common))
            idx_match <- match(date_common, time_obs)
            obs_aligned[!is.na(idx_match)] <- obs_vec[idx_match[!is.na(idx_match)]]
        } else {
            obs_aligned <- rep(NA_real_, length(date_common))
        }

        # --- Simulated (use corrected if available) ---
        if (has_corrected && sid %in% Station_cor_IDs) {
            sim_col_idx <- which(Station_cor_IDs == sid)
            sim_vec <- Q_cor_raw[, sim_col_idx[1]]
            sim_dates <- date_cor
        } else {
            sim_col_idx <- which(Station_sim_IDs == sid)
            if (length(sim_col_idx) == 0) next
            sim_vec <- Q_sim_raw[, sim_col_idx[1]]
            sim_dates <- date_sim
        }

        sim_aligned <- rep(NA_real_, length(date_common))
        idx_sim <- match(date_common, sim_dates)
        sim_aligned[!is.na(idx_sim)] <- sim_vec[idx_sim[!is.na(idx_sim)]]

        # Use catch_id as column name
        obs_daily[, (cid) := obs_aligned]
        sim_daily[, (cid) := sim_aligned]
    }

    return(list(obs = obs_daily, sim = sim_daily))
}

# --- Extract for strict matches -----------------------------------------------
cat("  Extracting strict matches...\n")
strict_matches <- matches[match_type == "strict"]
paired_strict <- extract_paired(strict_matches)

# --- Extract for all matches (strict + relaxed) -------------------------------
cat("  Extracting all matches (strict + relaxed)...\n")
paired_all <- extract_paired(matches)

# --- Aggregate to monthly means -----------------------------------------------
cat("Aggregating to monthly...\n")

aggregate_monthly <- function(daily_dt) {
    dt <- copy(daily_dt)
    dt[, year_month := format(date, "%Y-%m")]
    catch_cols <- setdiff(names(dt), c("date", "year_month"))
    monthly <- dt[, lapply(.SD, mean, na.rm = TRUE),
        by = year_month, .SDcols = catch_cols
    ]
    setnames(monthly, "year_month", "date")
    # Replace NaN with NA
    for (col in catch_cols) {
        monthly[is.nan(get(col)), (col) := NA_real_]
    }
    return(monthly)
}

obs_monthly_strict <- aggregate_monthly(paired_strict$obs)
obs_monthly_relaxed <- aggregate_monthly(paired_all$obs)
sim_monthly <- aggregate_monthly(paired_all$sim)

# Also create sim for strict subset (same columns as obs_strict)
strict_cols <- setdiff(names(obs_monthly_strict), "date")
sim_monthly_strict <- sim_monthly[, c("date", intersect(strict_cols, names(sim_monthly))),
    with = FALSE
]

# --- Save homogenized files ---------------------------------------------------
cat("Saving homogenized files...\n")

fwrite(obs_monthly_strict, file.path(out_dir, "obs_monthly_strict.csv"))
fwrite(obs_monthly_relaxed, file.path(out_dir, "obs_monthly_relaxed.csv"))
fwrite(sim_monthly, file.path(out_dir, "sim_monthly.csv"))
fwrite(sim_monthly_strict, file.path(out_dir, "sim_monthly_strict.csv"))

# Also save daily for detailed analysis
fwrite(paired_strict$obs, file.path(out_dir, "obs_daily_strict.csv"))
fwrite(paired_strict$sim, file.path(out_dir, "sim_daily_strict.csv"))
fwrite(paired_all$obs, file.path(out_dir, "obs_daily_relaxed.csv"))
fwrite(paired_all$sim, file.path(out_dir, "sim_daily.csv"))

cat("\n--- Summary ---\n")
cat("  Strict matches:", ncol(obs_monthly_strict) - 1, "catchments\n")
cat("  All matches:", ncol(obs_monthly_relaxed) - 1, "catchments\n")
cat("  Monthly rows:", nrow(sim_monthly), "\n")
cat("  Daily rows:", nrow(paired_all$sim), "\n")
cat("\nSaved to:", out_dir, "\n")
cat("Files:\n")
cat("  obs_monthly_strict.csv   — observed Q (strict matches only)\n")
cat("  obs_monthly_relaxed.csv  — observed Q (all matches)\n")
cat("  sim_monthly.csv          — simulated Q (all matches)\n")
cat("  sim_monthly_strict.csv   — simulated Q (strict matches only)\n")
cat("  obs_daily_strict.csv     — daily observed (strict)\n")
cat("  sim_daily_strict.csv     — daily simulated (strict)\n")
cat("  obs_daily_relaxed.csv    — daily observed (all)\n")
cat("  sim_daily.csv            — daily simulated (all)\n")


# =============================================================================
# Map of matched catchments (strict vs relaxed)
# =============================================================================
cat("Generating map of matched catchments...\n")

library(sf)
library(ggplot2)
library(rnaturalearth)

# Load catchments
catchments <- st_read(file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg"),
    quiet = TRUE
)
catchments_3035 <- st_transform(catchments, 3035)

# Basemap
basemap <- ne_countries(scale = "medium", returnclass = "sf") |>
    st_transform(3035)
bbox <- st_bbox(catchments_3035)

# Tag matched catchments
catchments_3035$match_type <- "Not matched"
strict_ids <- as.character(matches[match_type == "strict"]$catch_id)
relaxed_ids <- as.character(matches[match_type == "relaxed"]$catch_id)
catchments_3035$match_type[catchments_3035$catch_id %in% strict_ids] <- "Strict"
catchments_3035$match_type[catchments_3035$catch_id %in% relaxed_ids] <- "Relaxed"

catchments_3035$match_type <- factor(catchments_3035$match_type,
    levels = c("Not matched", "Relaxed", "Strict")
)

# Separate for layered plotting
matched <- catchments_3035[catchments_3035$match_type != "Not matched", ]

p_map <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(
        data = catchments_3035[catchments_3035$match_type == "Not matched", ],
        fill = "grey85", color = NA, alpha = 0.4
    ) +
    geom_sf(data = matched, aes(fill = match_type), color = NA, alpha = 0.8) +
    scale_fill_manual(
        values = c("Strict" = "#1b7837", "Relaxed" = "#fdae61"),
        name = "Match type"
    ) +
    coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(
        title = "Discharge validation: matched catchments",
        subtitle = sprintf(
            "Strict: %d | Relaxed: %d | Total: %d",
            length(strict_ids), length(relaxed_ids),
            length(strict_ids) + length(relaxed_ids)
        )
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom"
    )

ggsave(file.path(out_dir, "map_matched_catchments.png"),
    p_map,
    width = 10, height = 8, dpi = 200
)
cat("  -> Saved map_matched_catchments.png\n")
