# Comparison of AET between LISFLOOD model output and GLEAM v4.3a -----
# Extracts area-weighted mean AET from GLEAM yearly NetCDF rasters -----
# Computes Spearman correlation per catchment over overlapping period -----

# Library calling --------------------------------------------------
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rnaturalearth)
library(cowplot)
library(scales)

# Path configuration -----------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"

# Input paths
catchments_path <- file.path(
  base_dir, "data", "catchments_analysis_final_v3.gpkg"
)
gleam_dir <- file.path(base_dir, "data", "GLEAM", "AET_monthly")
lisflood_path <- file.path(
  base_dir, "data", "ActEvapo_nested_1951_2020.csv"
)

# Output paths
gleam_csv_path <- file.path(base_dir, "data", "gleam_aet_catchment.csv")
correlation_out_path <- file.path(
  base_dir, "output", "aet_correlation_summary.csv"
)

# Load and reproject catchments ------------------------------------
cat("[1/6] Loading and reprojecting catchment polygons...\n")

catchments <- st_read(catchments_path, quiet = TRUE)
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



# =============================================================================
# SEASONAL CYCLE COMPARISON: Mean monthly AET (GLEAM vs LISFLOOD)
# Uses homogenized monthly data from Diego pipeline
# =============================================================================
cat("[7/7] Computing seasonal cycle comparison...\n")

library(data.table)

# --- Load homogenized monthly data --------------------------------------------
in_dir <- file.path(base_dir, "output", "aet_diego", "1.homogenized")
# 
# obs_d <- data.table::fread(file.path(in_dir, "gleam_daily_homog.csv"))
# mod_d <- data.table::fread(file.path(in_dir, "lisflood_daily_homog.csv"))
obs_m <- data.table::fread(file.path(in_dir, "gleam_monthly_homog.csv"), header = TRUE)
mod_m <- data.table::fread(file.path(in_dir, "lisflood_monthly_homog.csv"), header = TRUE)
# 
# obs_d[, date := as.IDate(date)]
# mod_d[, date := as.IDate(date)]

daily_cols <- setdiff(names(mod_d), "date")
month_cols <- setdiff(names(mod_m), "date")

# Parse dates
obs_m[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
obs_m[, month := month(date)]

mod_m[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
mod_m[, month := month(date)]

# Catchment columns (common between obs and mod)
meta_cols <- c("date", "month")
obs_catches <- setdiff(names(obs_m), meta_cols)
mod_catches <- setdiff(names(mod_m), meta_cols)
common_cols <- intersect(obs_catches, mod_catches)

cat("  Common catchments:", length(common_cols), "\n")
cat("  Months:", nrow(obs_m), "(GLEAM),", nrow(mod_m), "(LISFLOOD)\n")

# --- Area weights -------------------------------------------------------------
if ("residual_area_km2" %in% names(catchments)) {
  area_vec <- catchments$residual_area_km2[match(
    common_cols, as.character(catchments$catch_id)
  )]
} else {
  area_vec <- rep(1, length(common_cols))
}
area_vec[is.na(area_vec)] <- 1
weights <- area_vec / sum(area_vec, na.rm = TRUE)

# --- Continental mean seasonal cycle (area-weighted) --------------------------
seasonal_gleam <- sapply(1:12, function(m) {
  mat <- as.matrix(obs_m[month == m, ..common_cols])
  monthly_means <- colMeans(mat, na.rm = TRUE)
  sum(monthly_means * weights, na.rm = TRUE)
})

seasonal_lf <- sapply(1:12, function(m) {
  mat <- as.matrix(mod_m[month == m, ..common_cols])
  monthly_means <- colMeans(mat, na.rm = TRUE)
  sum(monthly_means * weights, na.rm = TRUE)
})

seasonal_df <- data.frame(
  month = rep(1:12, 2),
  AET = c(seasonal_gleam, seasonal_lf),
  source = rep(c("GLEAM v4.3a", "HERA-WB"), each = 12)
)

# --- Plot: Continental mean seasonal cycle ------------------------------------
month_labels <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

p_seasonal <- ggplot(seasonal_df, aes(x = month, y = AET, color = source)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_color_manual(values = c("HERA-WB" = "steelblue", "GLEAM v4.3a" = "firebrick")) +
  labs(
    title = "Mean monthly AET: HERA-WB vs GLEAM (continental average)",
    subtitle = paste0("Area-weighted mean over ", length(common_cols), " catchments"),
    x = "Month", y = "Mean AET (mm/month)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_seasonal
ggsave(file.path(plots_dir, "aet_seasonal_cycle_continental.png"),
  p_seasonal,
  width = 8, height = 5, dpi = 200
)
cat("  -> Saved aet_seasonal_cycle_continental.png\n")

# --- Per-catchment seasonal bias (month by month) -----------------------------
seasonal_bias <- data.table(catch_id = common_cols)

for (m in 1:12) {
  gleam_m <- colMeans(as.matrix(obs_m[month == m, ..common_cols]), na.rm = TRUE)
  lf_m <- colMeans(as.matrix(mod_m[month == m, ..common_cols]), na.rm = TRUE)
  seasonal_bias[, paste0("bias_m", sprintf("%02d", m)) := lf_m - gleam_m]
}

# --- Plot: Seasonal bias boxplots ---------------------------------------------
bias_long <- melt(seasonal_bias,
  id.vars = "catch_id",
  variable.name = "month_col", value.name = "bias"
)
bias_long[, month := as.integer(sub("bias_m", "", month_col))]

p_bias_seasonal <- ggplot(bias_long, aes(x = factor(month), y = bias)) +
  geom_violin(alpha = 0.35, trim = FALSE, scale = "width", fill = "royalblue") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", fatten = 2) +
  # geom_boxplot(fill = "lightyellow", outlier.size = 0.5, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_x_discrete(labels = month_labels) +
  labs(
    title = "Monthly AET bias (HERA-WB - GLEAM) across all catchments",
    x = "Month", y = "Bias (mm/month)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(plots_dir, "aet_seasonal_bias_boxplot.png"),
  p_bias_seasonal,
  width = 8, height = 5, dpi = 200
)
cat("  -> Saved aet_seasonal_bias_boxplot.png\n")

# --- Combined figure ----------------------------------------------------------
fig_seasonal <- plot_grid(p_seasonal, p_bias_seasonal, ncol = 1, align = "v")
ggsave(file.path(plots_dir, "aet_seasonal_comparison.png"),
  fig_seasonal,
  width = 9, height = 18, dpi = 200
)
cat("  -> Saved aet_seasonal_comparison.png\n")

cat("Seasonal cycle comparison complete.\n")
