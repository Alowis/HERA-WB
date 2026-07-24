# Compare catchment daily SWE: GlobSnow vs LISFLOOD ----------------
# GlobSnow: data/globsnow_swe_catchment_daily.csv (wide: a `date` column +
#           one column per catchment with daily mean SWE).
# LISFLOOD: data/scovUps_nested_1951_2020.csv, one column per catchment,
#           6-hourly rows (no time column). Timestamps come from the raw
#           file (same row order). SWE is a STATE variable, so 6-hourly ->
#           daily aggregation uses the MEAN of the four steps (not the sum).
# Per catchment we align on date over the overlapping period and compute
# Spearman correlation plus a bias ratio and RMSE, then map + boxplot.

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
globsnow_csv <- file.path(base_dir, "data", "globsnow_swe_catchment_daily.csv")
# LISFLOOD SWE (scovUps): nested values + raw file for the 6-hourly timestamps
lisflood_nested <- file.path(base_dir, "data", "scovUps_nested_1951_2020.csv")
lisflood_raw <- file.path(
    base_dir, "data", "tss", "HERA_Histo", "scovUps_1951_2020.csv"
)

summary_out <- file.path(base_dir, "output", "swe_daily_correlation_summary.csv")
plots_dir <- file.path(base_dir, "output", "plots")

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

# Load GlobSnow daily SWE (wide) -----------------------------------
cat("[2/5] Loading GlobSnow daily SWE...\n")
if (!file.exists(globsnow_csv)) {
    stop(
        "GlobSnow CSV not found: ", globsnow_csv,
        " (run extract_globsnow_daily.R first)."
    )
}
gs <- data.table::fread(globsnow_csv, header = TRUE)
gs_long <- data.table::melt(
    gs,
    id.vars = "date", variable.name = "orig_id", value.name = "globsnow"
)
gs_long <- as.data.frame(gs_long)
gs_long$id_key <- norm_id(gs_long$orig_id)
gs_long$date <- as.Date(gs_long$date)
gs_long <- gs_long[, c("id_key", "date", "globsnow")]
globsnow_ids <- unique(gs_long$id_key)
cat(
    "  GlobSnow catchments:", length(globsnow_ids), "| date range:",
    format(min(gs_long$date)), "to", format(max(gs_long$date)), "\n"
)

# Load LISFLOOD 6-hourly SWE and aggregate to daily ----------------
cat("[3/5] Loading LISFLOOD 6-hourly SWE (nested) and aggregating to daily...\n")

# Timestamps from the raw file's leading time column (6-hourly);
# nested CSV rows are in the SAME order, so time aligns positionally.
time_vec <- data.table::fread(lisflood_raw, select = 1L)[[1]][-1]
time_vec <- time_vec[order(time_vec)]
nested <- data.table::fread(lisflood_nested, header = TRUE)

if (nrow(nested) != length(time_vec)) {
    stop(
        "Row mismatch: raw time has ", length(time_vec),
        " steps but nested CSV has ", nrow(nested), " rows."
    )
}

# 6-hourly -> daily MEAN (SWE is a state variable, so we average the four
# sub-daily steps rather than summing them). Grouping by day is order-free.
nested[, day := as.Date(substr(time_vec, 1, 10))]
val_cols <- setdiff(names(nested), "day")
daily <- nested[!is.na(day),
    lapply(.SD, mean, na.rm = TRUE),
    by = day, .SDcols = val_cols
]
data.table::setorder(daily, day)
daily <- daily[, -c("V1")]

# Long form for the join; normalise ids to match GlobSnow
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

paired <- dplyr::inner_join(gs_long, lf_long, by = c("id_key", "date"))
paired <- paired[!is.na(paired$globsnow) & !is.na(paired$lisflood), ]

metrics <- paired |>
    dplyr::group_by(id_key) |>
    dplyr::summarise(
        n_obs = dplyr::n(),
        mean_globsnow = mean(globsnow),
        mean_lisflood = mean(lisflood),
        bias = mean(lisflood) / mean(globsnow), # model / obs ratio
        rmse = sqrt(mean((globsnow - lisflood)^2)),
        spearman_r = if (dplyr::n() >= 10 && sd(globsnow) > 0 && sd(lisflood) > 0) {
            suppressWarnings(cor(globsnow, lisflood, method = "spearman"))
        } else {
            NA_real_
        },
        .groups = "drop"
    )

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
    "Median bias ratio:   ", round(median_b, 4), "\n",
    "Summary CSV:         ", summary_out, "\n",
    sep = ""
)

# --- Map + boxplot of Spearman r ----------------------------------
catch_cor <- dplyr::left_join(catchments, metrics, by = "id_key")
catch_cor_3035 <- sf::st_transform(catch_cor, crs = 3035)
nco <- sf::st_coordinates(sf::st_centroid(catch_cor_3035))

palet <- hcl.colors(9, palette = "viridis", rev = TRUE, fixup = TRUE)
limi <- c(0, 1)

parpl_valid <- catch_cor_3035[!is.na(catch_cor_3035$spearman_r), ]
parpl_missing <- catch_cor_3035[is.na(catch_cor_3035$spearman_r), ]

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
        title = "Daily SWE correlation: GlobSnow vs LISFLOOD (scovUps)"
    ) +
    ggplot2::theme_minimal()

cor_map

# --- Bias ratio map -----------------------------------------------
palet_b <- hcl.colors(9, palette = "RdYlBu", rev = TRUE, fixup = TRUE)
limi_b <- c(0, 2)
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
        colors = palet_b, n.breaks = 5, oob = scales::squish,
        limits = limi_b, name = "bias ratio"
    ) +
    ggplot2::coord_sf(
        xlim = c(min(nco[, 1]), max(nco[, 1])),
        ylim = c(min(nco[, 2]), max(nco[, 2]))
    ) +
    ggplot2::labs(
        x = "Longitude", y = "Latitude",
        title = "Daily SWE bias ratio (LISFLOOD / GlobSnow)"
    ) +
    ggplot2::theme_minimal()
bias_map




plot(metrics$mean_lisflood / metrics$mean_globsnow, ylim = c(0, 10))


# --- snow in lisflood -----------------


lf_stats <- lf_long |>
    dplyr::group_by(id_key) |>
    dplyr::summarise(
        n_obs = dplyr::n(),
        mean_lisflood = mean(lisflood),
        sd = sd(lisflood),
        q1 = quantile(lisflood, 0.05),
        q2 = quantile(lisflood, 0.95),
        med = median(lisflood),
        .groups = "drop"
    )
lf_stat_catch <- dplyr::left_join(catchments, lf_stats, by = "id_key")
lf_stat_catch_3035 <- sf::st_transform(lf_stat_catch, crs = 3035)

limi_lf <- c(0.1, 200)
palet_l <- hcl.colors(9, palette = "BuPu", rev = TRUE, fixup = TRUE)
lf_map <- ggplot2::ggplot(basemap) +
    ggplot2::geom_sf(fill = "gray95", color = NA) +
    ggplot2::geom_sf(
        data = lf_stat_catch_3035,
        ggplot2::aes(geometry = geom, fill = mean_lisflood),
        color = "transparent", alpha = 0.85
    ) +
    ggplot2::geom_sf(fill = NA, color = "grey20") +
    ggplot2::scale_fill_gradientn(
        colors = palet_l, breaks = c(0.1, 1, 10, 100), oob = scales::squish, trans = "log",
        limits = limi_lf, name = "mean SWE - lisflood"
    ) +
    ggplot2::coord_sf(
        xlim = c(min(nco[, 1]), max(nco[, 1])),
        ylim = c(min(nco[, 2]), max(nco[, 2]))
    ) +
    ggplot2::labs(
        x = "Longitude", y = "Latitude",
        title = "Mean SWE (LISFLOOD) - 1951-2020"
    ) +
    ggplot2::theme_minimal()

lf_map

r_vals <- metrics$spearman_r[!is.na(metrics$spearman_r)]
box_df <- data.frame(spearman_r = r_vals)
box_plot <- ggplot2::ggplot(box_df, ggplot2::aes(x = "", y = spearman_r)) +
    ggplot2::geom_boxplot(fill = "grey85", color = "grey30", width = 0.4) +
    ggplot2::stat_summary(
        fun = median, geom = "point", shape = 18, size = 4, color = "#d95f02"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(x = NULL, y = "Spearman r", title = paste0("n = ", length(r_vals))) +
    ggplot2::theme_minimal()

combined_plot <- cowplot::plot_grid(
    cor_map, box_plot,
    ncol = 2, rel_widths = c(3, 1), align = "h", axis = "tb"
)

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
    file.path(plots_dir, "swe_daily_correlation_map.png"),
    combined_plot,
    width = 14, height = 8, dpi = 150
)
ggplot2::ggsave(
    file.path(plots_dir, "swe_daily_bias_map.png"),
    bias_map,
    width = 10, height = 8, dpi = 150
)
cat("Saved: output/plots/swe_daily_correlation_map.png\n")

# --- Per-catchment daily time-series plot (optional helper) -------
# `year` (optional): if given, only that calendar year is plotted.
plot_swe_daily <- function(id, year = NULL, save = TRUE) {
    key <- norm_id(id)
    g <- gs_long[gs_long$id_key == key, c("date", "globsnow")]
    l <- lf_long[lf_long$id_key == key, c("date", "lisflood")]
    if (nrow(g) == 0) stop("Catchment '", id, "' not found in GlobSnow data.")
    df <- dplyr::inner_join(g, l, by = "date")

    # Optional filter to a single year
    if (!is.null(year)) {
        df <- df[as.integer(format(df$date, "%Y")) == year, ]
        if (nrow(df) == 0) {
            stop("No data for catchment '", id, "' in year ", year, ".")
        }
    }

    long <- tidyr::pivot_longer(
        df,
        cols = c("globsnow", "lisflood"),
        names_to = "source", values_to = "swe"
    )
    ttl <- if (is.null(year)) {
        paste("Daily SWE - Catchment", id)
    } else {
        paste0("Daily SWE - Catchment ", id, " (", year, ")")
    }
    p <- ggplot2::ggplot(long, ggplot2::aes(date, swe, color = source)) +
        ggplot2::geom_line(na.rm = TRUE) +
        ggplot2::scale_color_manual(
            values = c("globsnow" = "#1b9e77", "lisflood" = "#d95f02")
        ) +
        ggplot2::labs(
            title = ttl,
            x = "Date", y = "SWE", color = "Source"
        ) +
        ggplot2::theme_minimal()
    if (save) {
        dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
        suffix <- if (is.null(year)) key else paste0(key, "_", year)
        ggplot2::ggsave(
            file.path(plots_dir, paste0("swe_daily_ts_", suffix, ".png")),
            p,
            width = 10, height = 5, dpi = 150
        )
    }
    p
}

p1=plot_swe_daily(303662,1983,save=F)
p1
