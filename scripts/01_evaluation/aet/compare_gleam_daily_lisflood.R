# Compare catchment daily AET: GLEAM vs LISFLOOD --------------------
# GLEAM: one CSV per catchment (data/gleam_aet_by_catchment/gleam_aet_<Id>.csv)
#        with a daily `mean_w` (area-weighted mean AET, mm/day).
# LISFLOOD: data/ActEvapo_nested_1951_2020.csv, 2214 catchment-id columns,
#        one row per day starting 1951-01-01 (no date column).
# Per catchment we align on date over the overlapping period and compute
# Spearman correlation plus bias/RMSE, then map + boxplot the result.

# Library calling --------------------------------------------------
library(sf)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rnaturalearth)
library(cowplot)
library(scales)

# Path configuration -----------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"

catchments_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
gleam_dir <- file.path(base_dir, "data", "gleam_aet_by_catchment")
# LISFLOOD AET is 6-hourly. The nested CSV holds the values (one column per
# catchment, no time column); the raw file supplies the matching timestamps
# (same row order). Values are aggregated to daily totals (sum of 4 steps).
lisflood_nested <- file.path(base_dir, "data", "ActEvapo_nested_1951_2020.csv")
lisflood_raw <- file.path(
    base_dir, "data", "tss", "HERA_Histo", "ActEvapo_1951_2020.csv"
)

summary_out <- file.path(base_dir, "output", "aet_daily_correlation_summary.csv")
plots_dir <- file.path(base_dir, "output", "plots")

# Value column used from the GLEAM per-catchment files
gleam_value_col <- "mean_w"

# Helpers ----------------------------------------------------------
# Normalise ids for matching (strip leading zeros: "01029" -> "1029")
norm_id <- function(x) {
    x <- as.character(x)
    suppressWarnings(ifelse(grepl("^\\d+$", x), as.character(as.integer(x)), x))
}

# Load catchments --------------------------------------------------
cat("[1/5] Loading catchments...\n")
catchments <- st_read(catchments_path, quiet = TRUE)
catchments$id_key <- norm_id(catchments$catch_id)

# Load GLEAM per-catchment daily series ----------------------------
cat("[2/5] Loading GLEAM per-catchment daily AET...\n")
gleam_files <- list.files(
    gleam_dir,
    pattern = "^gleam_aet_.*\\.csv$", full.names = TRUE
)
if (length(gleam_files) == 0) {
    stop(
        "No GLEAM per-catchment files in: ", gleam_dir,
        " (run extract_gleam_aet_daily.R first)."
    )
}

# Each file -> data.frame(id_key, date, gleam)
gleam_list <- lapply(gleam_files, function(f) {
    d <- read.csv(f, colClasses = "character")
    data.frame(
        id_key = norm_id(d$Id),
        date = as.Date(d$date),
        gleam = as.numeric(d[[gleam_value_col]]),
        stringsAsFactors = FALSE
    )
})
gleam_long <- do.call(rbind, gleam_list)
gleam_ids <- unique(gleam_long$id_key)
cat(
    "  GLEAM catchments:", length(gleam_ids), "| date range:",
    format(min(gleam_long$date)), "to", format(max(gleam_long$date)), "\n"
)

# Load LISFLOOD 6-hourly AET and aggregate to daily ----------------
cat("[3/5] Loading LISFLOOD 6-hourly AET (nested) and aggregating to daily...\n")

# Timestamps come from the raw file's leading time column (6-hourly).
# The nested CSV rows are in the SAME order, so time aligns positionally.
time_vec <- data.table::fread(lisflood_raw, select = 1L)[[1]]
time_vec <- time_vec[order(time_vec)]
time_vec<-time_vec[-which(is.na(time_vec))]
# Read only the nested columns matching a GLEAM catchment (memory saver)
# nested_header <- (data.table::fread(lisflood_nested))
# lf_keys <- norm_id(nested_header)
# keep_idx <- which(lf_keys %in% gleam_ids)
# if (length(keep_idx) == 0) {
#     stop("No catchment ids shared between GLEAM files and LISFLOOD nested header.")
# }
nested <- data.table::fread(lisflood_nested, header = TRUE)

if (nrow(nested) != length(time_vec)) {
    stop(
        "Row mismatch: raw time has ", length(time_vec),
        " steps but nested CSV has ", nrow(nested), " rows."
    )
}

# 6-hourly (00/06/12/18) -> daily totals: sum the four steps per day.
nested[, day := as.Date(substr(time_vec, 1, 10))]
val_cols <- setdiff(names(nested), "day")
daily <- nested[, lapply(.SD, sum, na.rm = TRUE), by = day, .SDcols = val_cols]
data.table::setorder(daily, day)

# Long form for the join; normalise ids to match GLEAM
lf_long <- data.table::melt(
    daily,
    id.vars = "day", variable.name = "orig_id", value.name = "lisflood"
)
lf_long <- as.data.frame(lf_long)
lf_long$id_key <- norm_id(lf_long$orig_id)
lf_long$date <- lf_long$day
lf_long <- lf_long[, c("id_key", "date", "lisflood")]

# Align and compute per-catchment metrics --------------------------
cat("[4/5] Computing daily Spearman correlation and error metrics...\n")

paired <- dplyr::inner_join(gleam_long, lf_long, by = c("id_key", "date"))
paired <- paired[!is.na(paired$gleam) & !is.na(paired$lisflood), ]

metrics <- paired |>
    dplyr::group_by(id_key) |>
    dplyr::summarise(
        n_obs = dplyr::n(),
        mean_gleam = mean(gleam),
        mean_lisflood = mean(lisflood),
        bias = mean(lisflood)/mean(gleam), # GLEAM minus LISFLOOD
        rmse = sqrt(mean((gleam - lisflood)^2)),
        spearman_r = if (dplyr::n() >= 10 && sd(gleam) > 0 && sd(lisflood) > 0) {
            suppressWarnings(cor(gleam, lisflood, method = "spearman"))
        } else {
            NA_real_
        },
        .groups = "drop"
    )

# Quality category mirrors the monthly comparison script
metrics$category <- dplyr::case_when(
    is.na(metrics$spearman_r) & metrics$n_obs < 10 ~ "few_obs",
    is.na(metrics$spearman_r) ~ "no_variance",
    TRUE ~ "valid"
)

# Save summary -----------------------------------------------------
cat("[5/5] Saving summary and plots...\n")
dir.create(dirname(summary_out), recursive = TRUE, showWarnings = FALSE)
write.csv(metrics, summary_out, row.names = FALSE)

median_r <- median(metrics$spearman_r, na.rm = TRUE)
median_b <- median(metrics$bias, na.rm = TRUE)
cat(
    "\n=== Completion Summary ===\n",
    "Catchments compared: ", nrow(metrics), "\n",
    "Paired daily obs:    ", sum(metrics$n_obs), "\n",
    "Median Spearman r:   ", round(median_r, 4), "\n",
    "Median bias:         ", round(median_b, 4), "\n",
    "Summary CSV:         ", summary_out, "\n",
    sep = ""
)

# --- Map + boxplot of Spearman r (template style) -----------------
catch_cor <- dplyr::left_join(
    catchments, metrics,
    by = "id_key"
)
catch_cor_3035 <- sf::st_transform(catch_cor, crs = 3035)
nco <- sf::st_coordinates(sf::st_centroid(catch_cor_3035))

palet <- hcl.colors(9, palette = "viridis", rev = TRUE, fixup = TRUE)
limi <- c(0, 1)

parpl_valid <- catch_cor_3035[
    !is.na(catch_cor_3035$spearman_r),
]
parpl_missing <- catch_cor_3035[
    is.na(catch_cor_3035$spearman_r),
]

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
basemap <- sf::st_transform(world, crs = 3035)

cor_map <- ggplot2::ggplot(basemap) +
    ggplot2::geom_sf(fill = "gray95", color = NA) +
    ggplot2::geom_sf(
        data = parpl_missing, ggplot2::aes(geometry = geom),
        fill = "grey70", color = "transparent", alpha = 0.7
    ) +
    ggplot2::geom_sf(
        data = parpl_valid,
        ggplot2::aes(geometry = geom, fill = spearman_r),
        color = "transparent", alpha = 0.85
    ) +
    ggplot2::geom_sf(fill = NA, color = "grey20") +
    ggplot2::scale_fill_gradientn(
        colors = palet, n.breaks = 5, oob = scales::squish,
        limits = limi, name = "Spearman r"
    ) +
    ggplot2::coord_sf(
        xlim = c(min(nco[, 1]), max(nco[, 1])),
        ylim = c(min(nco[, 2]), max(nco[, 2]))
    ) +
    ggplot2::labs(
        x = "Longitude", y = "Latitude",
        title = "Daily AET correlation: GLEAM vs LISFLOOD"
    ) +
    ggplot2::theme_minimal()

cor_map

limi=c(0.8,1.2)
bias_map <- ggplot2::ggplot(basemap) +
  ggplot2::geom_sf(fill = "gray95", color = NA) +
  ggplot2::geom_sf(
    data = parpl_missing, ggplot2::aes(geometry = geom),
    fill = "grey70", color = "transparent", alpha = 0.7
  ) +
  ggplot2::geom_sf(
    data = parpl_valid,
    ggplot2::aes(geometry = geom, fill = bias),
    color = "transparent", alpha = 0.85
  ) +
  ggplot2::geom_sf(fill = NA, color = "grey20") +
  ggplot2::scale_fill_gradientn(
    colors = palet, n.breaks = 5, oob = scales::squish,
    limits = limi, name = "bias ratio"
  ) +
  ggplot2::coord_sf(
    xlim = c(min(nco[, 1]), max(nco[, 1])),
    ylim = c(min(nco[, 2]), max(nco[, 2]))
  ) +
  ggplot2::labs(
    x = "Longitude", y = "Latitude",
    title = "Daily AET correlation: GLEAM vs LISFLOOD"
  ) +
  ggplot2::theme_minimal()

bias_map

r_vals <- metrics$spearman_r[!is.na(metrics$spearman_r)]
box_df <- data.frame(spearman_r = r_vals)
box_plot <- ggplot2::ggplot(box_df, ggplot2::aes(x = "", y = spearman_r)) +
    ggplot2::geom_boxplot(fill = "grey85", color = "grey30", width = 0.4) +
    ggplot2::stat_summary(
        fun = median, geom = "point", shape = 18, size = 4, color = "#d95f02"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(-1, 1, 0.25)) +
    ggplot2::labs(x = NULL, y = "Spearman r", title = paste0("n = ", length(r_vals))) +
    ggplot2::theme_minimal()

combined_plot <- cowplot::plot_grid(
    cor_map, box_plot,
    ncol = 2, rel_widths = c(3, 1), align = "h", axis = "tb"
)

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
    file.path(plots_dir, "aet_daily_correlation_map.png"),
    combined_plot,
    width = 14, height = 8, dpi = 150
)
cat("Saved: output/plots/aet_daily_correlation_map.png\n")

# --- Per-catchment daily time-series plot (optional helper) -------
plot_aet_daily <- function(id, save = TRUE) {
    key <- norm_id(id)
    g <- gleam_long[gleam_long$id_key == key, c("date", "gleam")]
    l <- lf_long[lf_long$id_key == key, c("date", "lisflood")]
    if (nrow(g) == 0) stop("Catchment '", id, "' not found in GLEAM data.")
    df <- dplyr::inner_join(g, l, by = "date")
    long <- tidyr::pivot_longer(
        df,
        cols = c("gleam", "lisflood"),
        names_to = "source", values_to = "aet"
    )
    p <- ggplot2::ggplot(long, ggplot2::aes(date, aet, color = source)) +
        ggplot2::geom_line(na.rm = TRUE) +
        ggplot2::scale_color_manual(
            values = c("gleam" = "#1b9e77", "lisflood" = "#d95f02")
        ) +
        ggplot2::labs(
            title = paste("Daily AET - Catchment", id),
            x = "Date", y = "AET (mm/day)", color = "Source"
        ) +
        ggplot2::theme_minimal()
    if (save) {
        dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave(
            file.path(plots_dir, paste0("aet_daily_ts_", key, ".png")),
            p,
            width = 10, height = 5, dpi = 150
        )
    }
    p
}

plot1=plot_aet_daily(303662,save=F)
plot1
