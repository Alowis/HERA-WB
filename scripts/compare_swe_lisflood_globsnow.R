# Comparison of SWE between LISFLOOD model output and GlobSnow v3.0 -----
# Extracts area-weighted mean SWE from GlobSnow NetCDF rasters -----
# Computes spearman correlation per catchment over overlapping period -----

# Library calling --------------------------------------------------
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(rnaturalearth)
library(cowplot)

# Path configuration -----------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"


# Input paths
catchments_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
globsnow_dir <- file.path(base_dir, "data", "globsnow_swe_monthly")
lisflood_path <- file.path(base_dir, "data", "aggregates", "snow_water_equivalent", "snow_water_equivalent_monthly_all_years.csv")

# Output paths
globsnow_csv_path <- file.path(base_dir, "data", "globsnow_swe_catchment.csv")
correlation_out_path <- file.path(base_dir, "output", "swe_correlation_summary.csv")

# Load and reproject catchments ------------------------------------
cat("[1/6] Loading and reprojecting catchment polygons...\n")

catchments <- st_read(catchments_path, quiet = TRUE)
catchments$catch_id=as.character(as.numeric(catchments$catch_id))
ids_before <- catchments$catch_id

if (!isTRUE(st_crs(catchments) == st_crs(4326))) {
  catchments <- st_transform(catchments, 4326)
}

ids_after <- catchments$catch_id
missing_ids <- setdiff(ids_before, ids_after)
if (length(missing_ids) > 0) {
  stop(
    "Catchment IDs lost during reprojection: ",
    paste(missing_ids, collapse = ", ")
  )
}

# Extract GlobSnow SWE per catchment --------------------------------
cat("[2/6] Extracting GlobSnow SWE...\n")

nc_files <- list.files(
  globsnow_dir,
  pattern = "^\\d{6}_northern_hemisphere_monthly_swe_0\\.25grid\\.nc$",
  full.names = TRUE
)

if (length(nc_files) == 0) {
  stop("No GlobSnow NetCDF files found in: ", globsnow_dir)
}

n_files <- length(nc_files)
results_list <- vector("list", n_files)

for (i in seq_along(nc_files)) {
  print(i)
  f <- nc_files[i]
  fname <- basename(f)
  cat(
    "[2/6] Extracting GlobSnow SWE (file", i, "of",
    n_files, ")...\n"
  )

  r <- tryCatch(
    terra::rast(f),
    error = function(e) {
      warning("Could not read file: ", fname, " - skipping")
      NULL
    }
  )
  if (is.null(r)) next

  r[r < 0] <- NA
  yr <- substr(fname, 1, 4)
  mo <- substr(fname, 5, 6)
  date_str <- paste0(yr, "-", mo)

  # Reproject catchments to the raster CRS (EASE-Grid) before extraction.
  # exact_extract requires matching CRS; otherwise catchment coordinates
  # (degrees) are misread in the raster's metre space, collapsing all
  # polygons near the pole and producing a false split at 0 deg longitude.
  catchments_r <- sf::st_transform(catchments, terra::crs(r))

  means <- exact_extract(r, catchments_r, "mean")
  row_df <- data.frame(
    date = date_str,
    t(means),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(row_df) <- c("date", catchments$catch_id)

  results_list[[i]] <- row_df
}

# Remove NULLs from failed reads
results_list <- results_list[!vapply(
  results_list, is.null, logical(1)
)]

globsnow_df <- do.call(rbind, results_list)
globsnow_df <- globsnow_df[order(globsnow_df$date), ]
rownames(globsnow_df) <- NULL

# Sanity check: verify catchments west of 0 deg longitude have SWE data ----
# Guards against the CRS-mismatch bug that produced an artificial data
# boundary at the prime meridian (all-NA west of 0 deg).
west_ids <- catchments$catch_id[
  sf::st_coordinates(sf::st_centroid(sf::st_transform(catchments, 4326)))[, 1] < 0
]
west_ids <- as.character(west_ids)
west_cols <- intersect(west_ids, names(globsnow_df))

if (length(west_cols) > 0) {
  west_vals <- as.matrix(globsnow_df[, west_cols, drop = FALSE])
  n_west_with_data <- sum(colSums(!is.na(west_vals)) > 0)
  cat(
    "[check] Western catchments (lon < 0): ", length(west_cols),
    " | with at least one non-NA SWE value: ", n_west_with_data, "\n",
    sep = ""
  )
  if (n_west_with_data == 0) {
    stop(
      "No SWE data extracted for any catchment west of 0 deg longitude. ",
      "This usually indicates a CRS mismatch between the GlobSnow raster ",
      "(EASE-Grid, EPSG:3408) and the catchment polygons. Ensure catchments ",
      "are reprojected to terra::crs(r) before exact_extract()."
    )
  }
} else {
  cat("[check] No catchments west of 0 deg longitude found.\n")
}

# Save GlobSnow catchment CSV --------------------------------------
cat("[3/6] Saving GlobSnow catchment CSV...\n")

dir.create(dirname(globsnow_csv_path), recursive = TRUE, showWarnings = FALSE)
write.csv(globsnow_df, globsnow_csv_path, row.names = FALSE)


# start here if snow data already generated

globsnow_df <- read.csv(globsnow_csv_path)
# Load LISFLOOD SWE and align temporally ----------------------------
cat("[4/6] Loading LISFLOOD SWE and aligning temporally...\n")

lisflood_df <- read.csv(lisflood_path, check.names = FALSE)

lisflood_dates <- format(
  seq.Date(
    as.Date("1951-01-01"),
    as.Date("2020-12-01"),
    by = "month"
  ),
  "%Y-%m"
)

lisflood_df$date <- lisflood_dates

# Determine overlapping period as set intersection
overlap_dates <- intersect(lisflood_dates, globsnow_df$date)

# Compute spearman correlation per catchment -------------------------
cat("[5/6] Computing spearman correlations...\n")

catch_ids <- catchments$catch_id
n_catch <- length(catch_ids)


cor_results <- data.frame(
  catchment_id = catch_ids,
  spearman_r = NA_real_,
  n_obs = NA_integer_,
  p_value = NA_real_,
  # GlobSnow data quality category per catchment:
  #   "valid"      - has non-zero SWE observations (correlation computable)
  #   "all_zero"   - all GlobSnow values are 0 (maritime/low-snow catchment)
  #   "no_data"    - all GlobSnow values are NA (fully masked: water/mountain)
  #   "few_obs"    - some non-NA/non-zero but < 10 paired obs for correlation
  gs_category = NA_character_,
  stringsAsFactors = FALSE
)

# Subset both data frames to overlapping dates
gs_sub <- globsnow_df[globsnow_df$date %in% overlap_dates, ]
names(gs_sub)[-1] <- sub("^X", "", names(gs_sub)[-1])
lf_sub <- lisflood_df[lisflood_df$date %in% overlap_dates, ]

# matching col names

mcol <- match(colnames(gs_sub), colnames(lf_sub))

colnames(lf_sub)[mcol[2]]
colnames(gs_sub)[2]
for (i in seq_along(catch_ids)) {
  cid <- as.character(catch_ids[i])

  gs_vals <- gs_sub[[cid]]
  lf_vals <- lf_sub[[cid]]

  # Classify GlobSnow data quality for this catchment
  gs_all_na <- all(is.na(gs_vals))
  gs_all_zero <- !gs_all_na && all(gs_vals == 0, na.rm = TRUE)

  # Filter to pairwise complete (non-NA) observations
  complete <- !is.na(gs_vals) & !is.na(lf_vals)
  gs_complete <- gs_vals[complete]
  lf_complete <- lf_vals[complete]

  n <- length(gs_complete)
  cor_results$n_obs[i] <- n

  if (gs_all_na) {
    cor_results$gs_category[i] <- "no_data"
    cor_results$spearman_r[i] <- NA_real_
    cor_results$p_value[i] <- NA_real_
  } else if (gs_all_zero) {
    cor_results$gs_category[i] <- "all_zero"
    cor_results$spearman_r[i] <- NA_real_
    cor_results$p_value[i] <- NA_real_
  } else if (n < 10) {
    cor_results$gs_category[i] <- "few_obs"
    cor_results$spearman_r[i] <- NA_real_
    cor_results$p_value[i] <- NA_real_
  } else {
    cor_results$gs_category[i] <- "valid"
    ct <- cor.test(lf_complete, gs_complete, method = "spearman")
    cor_results$spearman_r[i] <- ct$estimate
    cor_results$p_value[i] <- ct$p.value
  }
}

# Save correlation summary and print completion ---------------------
cat("[6/6] Saving correlation summary...\n")

dir.create(
  dirname(correlation_out_path),
  recursive = TRUE,
  showWarnings = FALSE
)
write.csv(cor_results, correlation_out_path, row.names = FALSE)

# Validate at least one catchment and one GlobSnow file processed
if (n_catch < 1 || n_files < 1) {
  stop(
    "Invalid completion state: ",
    n_catch, " catchments and ",
    n_files, " GlobSnow files processed. ",
    "At least one of each is required."
  )
}

# Print completion summary
median_r <- median(cor_results$spearman_r, na.rm = TRUE)
cat(
  "\n=== Completion Summary ===\n",
  "Catchments processed: ", n_catch, "\n",
  "GlobSnow files processed: ", n_files, "\n",
  "Median spearman correlation: ", round(median_r, 4), "\n",
  sep = ""
)



# Section 7: Diagnostic Plots -----------------------------------------
cat("[7/7] Generating diagnostic plots...\n")

plot_swe_timeseries <- function(catchment_id,
                                gs_data = gs_sub,
                                lf_data = lf_sub,
                                output_dir = file.path(base_dir, "output", "plots"),
                                save = TRUE) {
  # --- Input validation ---
  # Validate catchment_id is a non-NULL, non-NA character scalar

  if (is.null(catchment_id) || length(catchment_id) != 1 ||
    !is.character(catchment_id) || is.na(catchment_id)) {
    stop("catchment_id must be a single valid character string (non-NULL, non-NA).")
  }

  # Validate catchment_id exists in GlobSnow data (check first per Req 2.3)
  if (!catchment_id %in% names(gs_data)) {
    stop("Catchment '", catchment_id, "' not found in GlobSnow data.")
  }

  # Validate catchment_id exists in LISFLOOD data
  if (!catchment_id %in% names(lf_data)) {
    stop("Catchment '", catchment_id, "' not found in LISFLOOD data.")
  }

  # --- Data preparation ---
  # Compute overlapping dates (set intersection)
  common_dates <- intersect(gs_data$date, lf_data$date)

  # Subset both data frames to overlapping dates
  gs_overlap <- gs_data[gs_data$date %in% common_dates, ]
  lf_overlap <- lf_data[lf_data$date %in% common_dates, ]

  # Build plot data frame with only overlapping dates
  plot_df <- data.frame(
    date = as.Date(paste0(gs_overlap$date, "-01")),
    GlobSnow = gs_overlap[[catchment_id]],
    LISFLOOD = lf_overlap[[catchment_id]]
  )

  # Pivot to long format
  plot_long <- tidyr::pivot_longer(
    plot_df,
    cols = c("GlobSnow", "LISFLOOD"),
    names_to = "source",
    values_to = "swe"
  )

  # --- Plot construction ---
  p <- ggplot2::ggplot(
    plot_long,
    ggplot2::aes(x = date, y = swe, color = source)
  ) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = c("GlobSnow" = "#1b9e77", "LISFLOOD" = "#d95f02")
    ) +
    ggplot2::labs(
      title = paste("SWE Comparison \u2013 Catchment", catchment_id),
      x = "Date",
      y = "SWE (mm)",
      color = "Source"
    ) +
    ggplot2::theme_minimal()

  # --- Save ---
  if (save) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    out_path <- file.path(
      output_dir,
      paste0("swe_timeseries_", catchment_id, ".png")
    )
    ggplot2::ggsave(out_path, p, width = 10, height = 5, dpi = 150)
    cat("Saved:", out_path, "\n")
  }

  invisible(p)
}


# --- Correlation Map and Boxplot (template style) -----------------

# Build basemap (world countries in EPSG:3035)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
basemap <- sf::st_transform(world, crs = 3035)

# Prepare catchment data in EPSG:3035
catch_cor <- dplyr::left_join(
  catchments,
  cor_results,
  by = c("catch_id" = "catchment_id")
)
catch_cor_3035 <- sf::st_transform(catch_cor, crs = 3035)

# Derive map extent from catchment centroids (mirrors nco pattern)
nco <- sf::st_coordinates(sf::st_centroid(catch_cor_3035))

# Plot parameters (consistent with your other scripts)
tsize <- 12
osize <- 10
palet <- hcl.colors(9, palette = "viridis", rev = TRUE, fixup = TRUE)
limi <- c(0, 1)
metric <- "Spearman r"

# Split catchments by GlobSnow category
parpl_valid <- catch_cor_3035[
  !is.na(catch_cor_3035$gs_category) & catch_cor_3035$gs_category == "valid",
]
parpl_allzero <- catch_cor_3035[
  !is.na(catch_cor_3035$gs_category) & catch_cor_3035$gs_category == "all_zero",
]
parpl_nodata <- catch_cor_3035[
  is.na(catch_cor_3035$gs_category) |
    catch_cor_3035$gs_category %in% c("no_data", "few_obs"),
]

# --- Map ---
cor_map <- ggplot2::ggplot(basemap) +
  # Background countries
  ggplot2::geom_sf(fill = "gray95", color = NA) +
  # No data / too few obs — grey
  ggplot2::geom_sf(
    data = parpl_nodata,
    ggplot2::aes(geometry = geom),
    fill = "gray22", color = "transparent", alpha = 0.7, shape = 21, stroke = 0
  ) +
  # All-zero SWE — yellow
  ggplot2::geom_sf(
    data = parpl_allzero,
    ggplot2::aes(geometry = geom),
    fill = "gray56", color = "transparent", alpha = 0.7, shape = 21, stroke = 0
  ) +
  # Valid: coloured by Spearman r
  ggplot2::geom_sf(
    data = parpl_valid,
    ggplot2::aes(geometry = geom, fill = spearman_r),
    color = "transparent", alpha = 0.7, shape = 21, stroke = 0
  ) +
  # Outline rings for valid catchments
  # ggplot2::geom_sf(
  #   data = parpl_valid,
  #   ggplot2::aes(geometry = geom),
  #   col = "grey20", alpha = 1, stroke = 0.05, shape = 1
  # ) +
  # Country borders on top
  ggplot2::geom_sf(fill = NA, color = "grey20") +
  ggplot2::scale_x_continuous(breaks = seq(-30, 40, by = 10)) +
  ggplot2::scale_fill_gradientn(
    colors   = palet,
    n.breaks = 5,
    oob      = scales::squish,
    limits   = limi,
    name     = metric
  ) +
  ggplot2::coord_sf(
    xlim = c(min(nco[, 1]), max(nco[, 1])),
    ylim = c(min(nco[, 2]), max(nco[, 2]))
  ) +
  ggplot2::labs(x = "Longitude", y = "Latitude") +
  ggplot2::guides(
    fill = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 12, reverse = FALSE)
  ) +
  ggplot2::theme(
    axis.title       = ggplot2::element_text(size = tsize),
    panel.background = ggplot2::element_rect(fill = "aliceblue", colour = "grey1"),
    panel.border     = ggplot2::element_rect(linetype = "solid", fill = NA, colour = "black"),
    legend.title     = ggplot2::element_text(size = tsize),
    legend.text      = ggplot2::element_text(size = osize),
    legend.position  = "right",
    panel.grid.major = ggplot2::element_line(colour = "grey85", linetype = "dashed"),
    panel.grid.minor = ggplot2::element_line(colour = "grey90"),
    legend.key       = ggplot2::element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size  = ggplot2::unit(0.8, "cm")
  )

cor_map
# --- Boxplot of Spearman r (valid catchments only) ---
r_vals <- cor_results$spearman_r[!is.na(cor_results$spearman_r)]
box_df <- data.frame(spearman_r = r_vals)

box_plot <- ggplot2::ggplot(box_df, ggplot2::aes(x = "", y = spearman_r)) +
  ggplot2::geom_boxplot(
    fill = "grey85", color = "grey30",
    width = 0.4, outlier.shape = 16, outlier.size = 1.5
  ) +
  ggplot2::stat_summary(
    fun = mean, geom = "point",
    shape = 18, size = 4, color = "#d95f02"
  ) +
  ggplot2::stat_summary(
    fun = mean, geom = "text",
    ggplot2::aes(label = paste0("mean=", round(ggplot2::after_stat(y), 2))),
    hjust = -0.2, size = 3.5, color = "#d95f02"
  ) +
  ggplot2::stat_summary(
    fun = median, geom = "text",
    ggplot2::aes(label = paste0("med=", round(ggplot2::after_stat(y), 2))),
    hjust = -0.2, vjust = 2, size = 3.5, color = "grey30"
  ) +
  ggplot2::scale_y_continuous(limits = c(-0.2, 1), breaks = seq(-1, 1, 0.1)) +
  ggplot2::labs(
    x = NULL,
    y = "Spearman r",
    title = paste0("n = ", length(r_vals))
  ) +
  ggplot2::theme_minimal(base_size = tsize) +
  ggplot2::theme(
    axis.text.x        = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    plot.title         = ggplot2::element_text(size = osize, hjust = 0.5)
  )

# --- Combine map + boxplot side by side ---
combined_plot <- cowplot::plot_grid(
  cor_map, box_plot,
  ncol = 2,
  rel_widths = c(3, 1),
  align = "h",
  axis = "tb"
)

# --- Save ---
plots_dir <- file.path(base_dir, "output", "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(
  file.path(plots_dir, "swe_correlation_map.png"),
  combined_plot,
  width  = 14,
  height = 8,
  dpi    = 150
)
cat("Saved: output/plots/swe_correlation_map.png\n")

# --- Example usage ---
# Uncomment below to generate time series plots for specific catchments:
# plot_swe_timeseries("01029")
# plot_swe_timeseries("01245")
#
# To get the plot object without saving:
# p <- plot_swe_timeseries("01029", save = FALSE)
# p + ggplot2::theme_bw()
