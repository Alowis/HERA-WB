# =============================================================================
# Seasonal cycle comparison for all validation variables
# Uses homogenized monthly data from Diego pipelines
#
# Variables:
#   - AET: GLEAM v4.3a vs LISFLOOD
#   - SWE: GlobSnow v3.0 vs LISFLOOD
#   - Soil Moisture: ESA CCI vs LISFLOOD
#   - Discharge: Observed vs LISFLOOD (from matched stations)
#
# Produces for each variable:
#   1. Continental mean seasonal cycle (2 curves)
#   2. Per-catchment monthly bias boxplots
# =============================================================================

library(data.table)
library(ggplot2)
library(sf)
library(lubridate)
library(cowplot)
library(terra)
library(exactextractr)
library(scales)
library(rnaturalearth)
source("R/load_workspace.R")
# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
plots_dir <- file.path(base_dir, "output", "figures")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Homogenized data directories
aet_dir <- file.path(base_dir, "output", "aet_diego", "1.homogenized")
swe_dir <- file.path(base_dir, "output", "swe_diego", "1.homogenized")
sm_dir <- file.path(base_dir, "output", "soil_moisture_diego", "1.Diego_Merged")
dis_dir <- file.path(base_dir, "data", "aggregate", "Q")

# --- Load catchments for area weighting ---------------------------------------
catchments <- st_read(gpkg_path, quiet = TRUE)
catchments_3035 <- st_transform(catchments, 3035)

# Basemap for maps
basemap <- ne_countries(scale = "medium", returnclass = "sf") |> st_transform(3035)
bbox <- st_bbox(catchments_3035)

# Remove Iceland
centroids_wgs <- st_coordinates(st_centroid(st_transform(catchments_3035, 4326)))
iceland_mask <- centroids_wgs[, 1] < -13 & centroids_wgs[, 2] > 63
catchments_3035 <- catchments_3035[!iceland_mask, ]
catchments <- catchments[!iceland_mask, ]

# --- Month labels -------------------------------------------------------------
month_labels <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

# =============================================================================
# Helper function: compute and plot seasonal comparison
# =============================================================================
seasonal_comparison <- function(obs_dt, mod_dt, var_name, obs_label, mod_label,
                                y_label, color_obs = "firebrick",
                                color_mod = "steelblue") {
  # Find common catchments
  meta <- c("date", "month", "year")
  obs_cols <- setdiff(names(obs_dt), meta)
  mod_cols <- setdiff(names(mod_dt), meta)
  common_cols <- intersect(obs_cols, mod_cols)

  if (length(common_cols) < 10) {
    cat("  WARNING: Only", length(common_cols), "common catchments for", var_name, "\n")
  }
  cat("  Common catchments:", length(common_cols), "\n")

  # Area weights
  catch_ids_clean <- sub("^X", "", common_cols)
  area_vec <- catchments$residual_area_km2[match(
    catch_ids_clean, as.character(catchments$catch_id)
  )]
  area_vec[is.na(area_vec)] <- 1
  weights <- area_vec / sum(area_vec, na.rm = TRUE)

  # Continental mean seasonal cycle
  seasonal_obs <- sapply(1:12, function(m) {
    mat <- as.matrix(obs_dt[month == m, ..common_cols])
    means <- colMeans(mat, na.rm = TRUE)
    sum(means * weights, na.rm = TRUE)
  })

  seasonal_mod <- sapply(1:12, function(m) {
    mat <- as.matrix(mod_dt[month == m, ..common_cols])
    means <- colMeans(mat, na.rm = TRUE)
    sum(means * weights, na.rm = TRUE)
  })

  seasonal_df <- data.frame(
    month = rep(1:12, 2),
    value = c(seasonal_obs, seasonal_mod),
    source = rep(c(obs_label, mod_label), each = 12)
  )

  # Plot 1: Seasonal cycle
  p_cycle <- ggplot(seasonal_df, aes(x = month, y = value, color = source)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_x_continuous(breaks = 1:12, labels = month_labels) +
    scale_color_manual(values = setNames(
      c(color_mod, color_obs), c(mod_label, obs_label)
    )) +
    labs(
      title = paste0("Mean monthly ", var_name, ": ", mod_label, " vs ", obs_label),
      subtitle = paste0("Area-weighted mean over ", length(common_cols), " catchments"),
      x = "Month", y = y_label, color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  # Compute per-catchment bias per month
  bias_dt <- data.table(catch_id = common_cols)
  for (m in 1:12) {
    obs_means <- colMeans(as.matrix(obs_dt[month == m, ..common_cols]), na.rm = TRUE)
    mod_means <- colMeans(as.matrix(mod_dt[month == m, ..common_cols]), na.rm = TRUE)
    bias_dt[, paste0("m", sprintf("%02d", m)) := mod_means - obs_means]
  }

  # Plot 2: Bias violin + boxplot
  bias_long <- melt(bias_dt,
    id.vars = "catch_id",
    variable.name = "month_col", value.name = "bias"
  )
  bias_long[, month := as.integer(sub("m", "", month_col))]

  p_bias <- ggplot(bias_long, aes(x = factor(month), y = bias)) +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", fill = "royalblue") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", fatten = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_x_discrete(labels = month_labels) +
    labs(
      title = paste0("Monthly ", var_name, " bias (", mod_label, " - ", obs_label, ")"),
      x = "Month", y = paste0("Bias (", y_label, ")")
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))


  # --- Plot 3: Bias map (mean annual bias per catchment) --------------------
  bias_annual <- bias_dt[, .(mean_bias = rowMeans(.SD, na.rm = TRUE)),
    .SDcols = paste0("m", sprintf("%02d", 1:12))
  ]
  bias_annual$catch_id <- sub("^X", "", bias_dt$catch_id)

  cats_bias <- catchments_3035[as.character(catchments_3035$catch_id) %in% bias_annual$catch_id, ]
  cats_bias <- merge(cats_bias, bias_annual,
    by.x = "catch_id", by.y = "catch_id", all.x = FALSE
  )

  bias_lim <- quantile(abs(cats_bias$mean_bias), 0.95, na.rm = TRUE)
  palet_bias <- hcl.colors(11, palette = "RdBu", rev = TRUE)

  p_map <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_bias, aes(fill = mean_bias), color = NA) +
    scale_fill_gradientn(
      colors = palet_bias,
      limits = c(-bias_lim, bias_lim), oob = squish,
      name = paste0("Mean bias\n(", y_label, ")")
    ) +
    coord_sf(
      xlim = c(bbox["xmin"], bbox["xmax"]),
      ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(
      title = paste0("Mean annual bias: ", mod_label, " - ", obs_label),
      subtitle = var_name
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "right"
    )

  # --- Plot 4: Continental temporal evolution (annual means) -----------------
  obs_dt[, year := year(date)]
  mod_dt[, year := year(date)]

  # Remove incomplete years (keep only years with 12 months of data)


  complete_years_obs <- obs_dt[, .N, by = year][N >= 360, year]
  complete_years_mod <- mod_dt[, .N, by = year][N >= 360, year]
  complete_years <- intersect(complete_years_obs, complete_years_mod)


  annual_obs <- obs_dt[year %in% complete_years, .(annual_mean = sum(
    sapply(common_cols, function(col) mean(.SD[[col]], na.rm = TRUE)) * weights,
    na.rm = TRUE
  )), by = year, .SDcols = common_cols]

  annual_mod <- mod_dt[year %in% complete_years, .(annual_mean = sum(
    sapply(common_cols, function(col) mean(.SD[[col]], na.rm = TRUE)) * weights,
    na.rm = TRUE
  )), by = year, .SDcols = common_cols]


  # Simpler approach: compute weighted mean per year
  annual_obs2 <- obs_dt[,
    {
      mat <- as.matrix(.SD)
      means <- colMeans(mat, na.rm = TRUE)
      .(val = sum(means * weights, na.rm = TRUE))
    },
    by = year,
    .SDcols = common_cols
  ]
  annual_obs[, source := obs_label]

  annual_mod2 <- mod_dt[,
    {
      mat <- as.matrix(.SD)
      means <- colMeans(mat, na.rm = TRUE)
      .(val = sum(means * weights, na.rm = TRUE))
    },
    by = year,
    .SDcols = common_cols
  ]
  annual_mod[, source := mod_label]

  annual_combined <- rbind(annual_obs, annual_mod)

  p_temporal <- ggplot(annual_combined, aes(x = year, y = annual_mean, color = source)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 1) +
    scale_color_manual(values = setNames(
      c(color_mod, color_obs), c(mod_label, obs_label)
    )) +
    labs(
      title = paste0("Continental ", var_name, ": annual evolution"),
      x = NULL, y = y_label, color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  # Combined
  fig <- plot_grid(p_cycle, p_bias, p_map, p_temporal, ncol = 2, align = "h")

  return(list(
    cycle = p_cycle, bias = p_bias, map = p_map,
    temporal = p_temporal, combined = fig
  ))
}


bias_trends_comparison <- function(obs_dt, mod_dt, var_name, obs_label, mod_label,
                                   y_label, color_obs = "firebrick",
                                   color_mod = "steelblue", dt = "monthly") {
  # Find common catchments
  meta <- c("date", "month", "year")
  obs_cols <- setdiff(names(obs_dt), meta)
  mod_cols <- setdiff(names(mod_dt), meta)
  common_cols <- intersect(obs_cols, mod_cols)

  if (length(common_cols) < 10) {
    cat("  WARNING: Only", length(common_cols), "common catchments for", var_name, "\n")
  }
  cat("  Common catchments:", length(common_cols), "\n")

  # Area weights
  catch_ids_clean <- sub("^X", "", common_cols)
  area_vec <- catchments$residual_area_km2[match(
    catch_ids_clean, as.character(catchments$catch_id)
  )]
  area_vec[is.na(area_vec)] <- 1
  weights <- area_vec / sum(area_vec, na.rm = TRUE)


  # Compute per-catchment bias per month
  bias_dt <- data.table(catch_id = common_cols)
  for (m in 1:12) {
    obs_means <- colMeans(as.matrix(obs_dt[month == m, ..common_cols]), na.rm = TRUE)
    mod_means <- colMeans(as.matrix(mod_dt[month == m, ..common_cols]), na.rm = TRUE)
    bias_dt[, paste0("m", sprintf("%02d", m)) := mod_means - obs_means]
  }

  # Plot 2: Bias violin + boxplot
  bias_long <- melt(bias_dt,
    id.vars = "catch_id",
    variable.name = "month_col", value.name = "bias"
  )
  bias_long[, month := as.integer(sub("m", "", month_col))]

  p_bias <- ggplot(bias_long, aes(x = factor(month), y = bias)) +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", fill = "royalblue") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", fatten = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_x_discrete(labels = month_labels) +
    labs(
      title = paste0("(b)"),
      x = "Month", y = paste0("Bias (", y_label, ")")
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))


  # --- Plot 3: Bias map (mean annual bias per catchment) --------------------
  bias_annual <- bias_dt[, .(mean_bias = rowMeans(.SD, na.rm = TRUE)),
    .SDcols = paste0("m", sprintf("%02d", 1:12))
  ]
  bias_annual$catch_id <- as.character(as.numeric(sub("^X", "", bias_dt$catch_id)))
  cats_bias <- catchments_3035[as.character(as.numeric(catchments_3035$catch_id)) %in% bias_annual$catch_id, ]
  cats_bias$catch_id <- as.character(as.numeric(cats_bias$catch_id))
  cats_bias <- merge(cats_bias, bias_annual,
    by.x = "catch_id", by.y = "catch_id", all.x = FALSE
  )

  bias_lim <- quantile(abs(cats_bias$mean_bias), 0.95, na.rm = TRUE)
  palet_bias <- hcl.colors(11, palette = "RdBu", rev = F)

  p_map <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_bias, aes(fill = mean_bias), color = "gray44", linewidth = 0.02) +
    scale_fill_gradientn(
      colors = palet_bias,
      limits = c(-bias_lim, bias_lim), oob = squish,
      guide = guide_colorbar(
        direction = "vertical", title.position = "top",
        barwidth = .5, barheight = 10
      ),
      name = paste0("Mean bias (", y_label, ")")
    ) +
    coord_sf(
      xlim = c(bbox["xmin"], bbox["xmax"]),
      ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
    ) +
    labs(
      title = paste0("(c)"),
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "right"
    )

  # --- Plot 4: Continental temporal evolution (annual means) -----------------
  obs_dt[, year := year(date)]
  mod_dt[, year := year(date)]

  obs_dt[["year"]]
  # Remove incomplete years (keep only years with 12 months of data)

  if (dt == "monthly") nm <- 12
  if (dt == "daily") nm <- 360

  if (var_name == "SWE") nm <- nm * 0.6

  complete_years_obs <- obs_dt[, .N, by = year][N >= nm, year]
  complete_years_mod <- mod_dt[, .N, by = year][N >= nm, year]
  complete_years <- intersect(complete_years_obs, complete_years_mod)


  annual_obs <- obs_dt[year %in% complete_years, .(annual_mean = sum(
    sapply(common_cols, function(col) mean(.SD[[col]], na.rm = TRUE)) * weights,
    na.rm = TRUE
  )), by = year, .SDcols = common_cols]

  annual_mod <- mod_dt[year %in% complete_years, .(annual_mean = sum(
    sapply(common_cols, function(col) mean(.SD[[col]], na.rm = TRUE)) * weights,
    na.rm = TRUE
  )), by = year, .SDcols = common_cols]


  # Simpler approach: compute weighted mean per year
  annual_obs2 <- obs_dt[,
    {
      mat <- as.matrix(.SD)
      means <- colMeans(mat, na.rm = TRUE)
      .(val = sum(means * weights, na.rm = TRUE))
    },
    by = year,
    .SDcols = common_cols
  ]
  annual_obs[, source := obs_label]

  annual_mod2 <- mod_dt[,
    {
      mat <- as.matrix(.SD)
      means <- colMeans(mat, na.rm = TRUE)
      .(val = sum(means * weights, na.rm = TRUE))
    },
    by = year,
    .SDcols = common_cols
  ]
  annual_mod[, source := mod_label]

  annual_combined <- rbind(annual_obs, annual_mod)
  maxo <- max(annual_obs$annual_mean)
  maxm <- max(annual_mod$annual_mean)
  bias <- mean(annual_mod$annual_mean) / mean(annual_obs$annual_mean)
  bias <- (bias - 1) * 100

  p_temporal <- ggplot(annual_combined, aes(x = year, y = annual_mean, color = source)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 1) +
    scale_color_manual(values = setNames(
      c(color_mod, color_obs), c(mod_label, obs_label)
    )) +
    # scale_y_continuous(limits = c(0,max(maxo,maxm)))+
    labs(
      title = paste0("(a)"),
      x = NULL, y = y_label, color = NULL
    ) +
    annotate("text",
      x = min(annual_obs$year) + 3, y = max(annual_obs$annual_mean),
      label = paste0("mean bias = ", round(bias, 1), "%"), hjust = 0, size = 4
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  # Combined
  bottom_row <- plot_grid(p_bias, p_map, ncol = 2, align = "h")
  fig <- plot_grid(p_temporal, bottom_row, nrow = 2, rel_heights = c(1, 1))


  return(list(
    bias = p_bias, map = p_map,
    temporal = p_temporal, combined = fig
  ))
}

norm_id <- function(x) {
  x <- as.character(x)
  suppressWarnings(ifelse(grepl("^\\d+$", x), as.character(as.integer(x)), x))
}
# =============================================================================
# 1. AET: GLEAM vs LISFLOOD
# =============================================================================
cat("[1/4] AET seasonal comparison...\n")

obs_aet <- fread(file.path(aet_dir, "gleam_monthly_homog.csv"), header = TRUE)
mod_aet <- fread(file.path(aet_dir, "lisflood_monthly_homog.csv"), header = TRUE)

obs_aet[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
obs_aet[, month := month(date)]
mod_aet[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
mod_aet[, month := month(date)]

res_aet <- bias_trends_comparison(
  obs_dt = obs_aet, mod_dt = mod_aet, var_name = "AET",
  obs_label = "GLEAM v4.3a", mod_label = "HERA-WB",
  y_label = "mm/month", dt = "monthly"
)

res_aet$combined
ggsave(file.path(plots_dir, "Supplement_Figure_AET.png"), res_aet$combined,
  width = 14, height = 11, dpi = 400
)

# =============================================================================
# 2. SWE: GlobSnow vs LISFLOOD
# =============================================================================
cat("[2/4] SWE seasonal comparison...\n")
start_year <- 1978
obs_swe <- fread(file.path(swe_dir, "globsnow_monthly_homog.csv"), header = TRUE)
mod_swe <- fread(file.path(swe_dir, "lisflood_monthly_homog.csv"), header = TRUE)

obs_swe[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
obs_swe[, month := month(date)]
mod_swe[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
mod_swe[, month := month(date)]

res_swe <- bias_trends_comparison(
  obs_dt = obs_swe, mod_dt = mod_swe, "SWE",
  obs_label = "GlobSnow v3.0", mod_label = "HERA-WB",
  y_label = "mm/month", dt = "monthly"
)

res_swe$combined
ggsave(file.path(plots_dir, "seasonal_SWE.png"), res_swe$combined,
  width = 9, height = 9, dpi = 200
)
ggsave(file.path(plots_dir, "seasonal_SWE_map.png"), res_swe$map,
  width = 10, height = 8, dpi = 200
)
ggsave(file.path(plots_dir, "seasonal_SWE_temporal.png"), res_swe$temporal,
  width = 10, height = 5, dpi = 200
)

mod_swe_days <- fread(file.path(swe_dir, "lisflood_daily_homog.csv"), header = TRUE)
mod_swe_days[, year := year(date)]
mod_swe_days[, month := month(date)]
mod_swe_days <- mod_swe_days[year >= start_year, ]

# 1. Identify catchments with snow every year (annual max SWE > 0 for all years)

date_col <- names(obs_swe)[1]
obs_swe[, year := year(date)]
mod_swe[, year := year(date)]

catch_cols_swe <- setdiff(names(mod_swe), c(date_col, "month", "year"))

annual_max_swe <- mod_swe[, lapply(.SD, max, na.rm = TRUE),
  by = year, .SDcols = catch_cols_swe
]


message("Applying snow catchment filter...")

mod_swe_days[, year := year(date)]

# Identify catchments with meaningful snow cover
catch_cols_swe <- setdiff(names(mod_swe_days), c("date", "year"))

n_days_total <- nrow(mod_swe_days)
snow_pct <- mod_swe_days[, lapply(.SD, function(x) sum(x > 5, na.rm = TRUE) / n_days_total),
  .SDcols = catch_cols_swe
]



annual_max_swe <- mod_swe_days[, lapply(.SD, max, na.rm = TRUE),
  by = year, .SDcols = catch_cols_swe
]

years_no_snow <- colSums(annual_max_swe[, ..catch_cols_swe] <= 5, na.rm = TRUE)
# Keep only catchments with snow every year
snow_catches <- names(years_no_snow[years_no_snow <= 20])
snow_catches <- names(snow_pct)[as.numeric(snow_pct[1, ]) > 0.01]

cat(
  "  Catchments with >= 1% snow days (SWE > 5mm):", length(snow_catches),
  "out of", length(daily_cols), "\n"
)


# filter catchments based on standard deviation of elevation

elev_std_path <- file.path(base_dir, "data", "elvstd_European_01min.nc") # adjust path if needed
r_elev_std <- rast(elev_std_path)

dem_path <- file.path(base_dir, "data", "dem.nc") # adjust path if needed
r_dem <- rast(dem_path)

# Ensure catchments are in WGS84 for extraction
shp <- catchments
shp_wgs <- if (st_crs(shp) == st_crs(4326)) shp else st_transform(shp, 4326)

# Area-weighted mean of elevation std within each catchment
shp_wgs$elev_std <- exact_extract(r_elev_std, shp_wgs, "mean")
shp_wgs$elev_std2 <- exact_extract(r_dem, shp_wgs, "stdev")

# Filter: elevation std < 50
flat_catches <- norm_id(shp_wgs$catch_id[shp_wgs$elev_std2 < 200])


catches_iceland <- ws$catch_ids[-match(ws$iceland_ids, ws$catch_ids)]



# 3. Intersect both criteria
snow_flat <- intersect(norm_id(sub("^X", "", snow_catches)), flat_catches)
snow_flat <- intersect(catches_iceland, snow_flat)

cat("Catchments with frequent snow AND elev_std < 200:", length(snow_flat), "\n")

obs_sweX <- obs_swe[, .SD, .SDcols = c(date_col, snow_flat)]
mod_sweX <- mod_swe[, .SD, .SDcols = c(date_col, snow_flat)]
obs_sweX[, month := month(date)]
mod_sweX[, month := month(date)]

res_swe_locs <- bias_trends_comparison(obs_sweX, mod_sweX, "SWE",
  obs_label = "GlobSnow v3.0", mod_label = "HERA-WB",
  y_label = "mm/month", dt = "monthly"
)

res_swe_locs$combined
ggsave(file.path(plots_dir, "Supplement_Figure_SWE.png"), res_swe_locs$combined,
  width = 14, height = 11, dpi = 400
)




cat("Loading catchments...\n")
catchments <- st_read(gpkg_path, quiet = TRUE)
catchments_3035 <- st_transform(catchments, 3035)
library(rnaturalearth)
basemap <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(3035)
bbox <- st_bbox(catchments_3035)




catchments_3035$xsnow <- snow_flat[
  match(as.numeric(catchments_3035$catch_id), as.numeric((snow_flat)))
]


# subsample
catchments_snow <- catchments_3035[which(!is.na(catchments_3035$xsnow)), ]


cats_map <- catchments_snow

mean_swe_mod <- mod_sweX[, lapply(.SD, mean, na.rm = TRUE),
  .SDcols = snow_flat
]

mean_swe_obs <- obs_sweX[, lapply(.SD, mean, na.rm = TRUE),
  .SDcols = snow_flat
]

mean(as.numeric(mean_swe_obs), na.rm = T)
bias <- as.numeric(mean_swe_mod - mean_swe_obs)
median(bias, na.rm = T)
cats_map$bias <- bias[
  match(as.numeric(cats_map$catch_id), as.numeric((snow_flat)))
]
cats_map$mod <- as.numeric(mean_swe_mod)[
  match(as.numeric(cats_map$catch_id), as.numeric((snow_flat)))
]
cats_map$obs <- as.numeric(mean_swe_obs)[
  match(as.numeric(cats_map$catch_id), as.numeric((snow_flat)))
]
cats_map$catch_id <- as.numeric(cats_map$catch_id)
# map of no_snow days
palet <- hcl.colors(11, palette = "RdBu", rev = F)
map <- ggplot() +
  geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.15) +
  # geom_sf(
  #     data = cats_map[!cats_map$significant, ],
  #     aes(fill = slope), color = NA, alpha = 0.3
  # ) +
  geom_sf(
    data = cats_map, aes(fill = bias), color = "gray44", alpha = 0.6
  ) +
  scale_fill_gradientn(
    colors = palet,
    limits = c(-10, 10),
    oob = scales::squish,
    name = "bias (mm)"
  ) +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
  ) +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title.position = "top"
  )

# =============================================================================
# 3. Soil Moisture: ESA CCI vs LISFLOOD
# =============================================================================
cat("[3/4] Soil Moisture seasonal comparison...\n")
start_year <- 1991
obs_sm <- fread(file.path(sm_dir, "esacci_homogenized.csv"), header = TRUE)
mod_sm <- fread(file.path(sm_dir, "lisflood_homogenized.csv"), header = TRUE)



obs_sm[, date := as.Date(date)]
obs_sm[, month := month(date)]
obs_sm[, year := year(date)]

mod_sm[, date := as.Date(date)]
mod_sm[, month := month(date)]
mod_sm[, year := year(date)]

obs_sm <- obs_sm[year >= start_year, ]
mod_sm <- mod_sm[year >= start_year, ]


# Remove X prefix from column names for matching
obs_sm_cols <- setdiff(names(obs_sm), c("date", "month"))
clean_obs <- sub("^X", "", obs_sm_cols)
setnames(obs_sm, obs_sm_cols, clean_obs)

mod_sm_cols <- setdiff(names(mod_sm), c("date", "month"))
clean_mod <- sub("^X", "", mod_sm_cols)
setnames(mod_sm, mod_sm_cols, clean_mod)

res_sm <- bias_trends_comparison(obs_sm, mod_sm, "Soil Moisture",
  obs_label = "ESA CCI", mod_label = "HERA-WB",
  y_label = "m\u00b3/m\u00b3", dt = "daily"
)
res_sm$combined


ggsave(file.path(plots_dir, "seasonal_SM.png"), res_sm$combined,
  width = 9, height = 9, dpi = 200
)
ggsave(file.path(plots_dir, "seasonal_SM_map.png"), res_sm$map,
  width = 10, height = 8, dpi = 200
)
ggsave(file.path(plots_dir, "seasonal_SM_temporal.png"), res_sm$temporal,
  width = 10, height = 5, dpi = 200
)


# mask snow days

mod_swe_days <- fread(file.path(swe_dir, "lisflood_daily_homog.csv"), header = TRUE)
mod_swe_days[, year := year(date)]
mod_swe_days[, month := month(date)]
mod_swe_days <- mod_swe_days[year >= start_year, ]

common_dates <- intersect(mod_swe_days$date, mod_sm$date)


# Identify days with snow (SWE > 0) per catchment — returns logical matrix
catch_cols_swe <- setdiff(names(mod_swe_days), c("date", "month"))

days_with_snow <- mod_swe_days[, lapply(.SD, function(x) x > 0), .SDcols = catch_cols_swe]


# Count snow days per year
mod_swe_days[, year := year(as.Date(date))]
no_snow_days_per_year <- mod_swe_days[, lapply(.SD, function(x) sum(x == 0, na.rm = TRUE)),
  by = year, .SDcols = catch_cols_swe
]

no_snow_days <- mod_swe_days[, lapply(.SD, function(x) sum(x == 0, na.rm = TRUE)),
  .SDcols = catch_cols_swe
]

# mask snow days
# For each catchment, set obs/mod to NA where SWE > 5
catch_cols <- intersect(names(mod_swe_days), names(obs_sm))

obs_sm <- obs_sm[date %in% common_dates, c(catch_cols), with = FALSE]
mod_sm <- mod_sm[date %in% common_dates, c(catch_cols), with = FALSE]
mod_swe_days <- mod_swe_days[date %in% common_dates, c(catch_cols), with = FALSE]

catch_cols_na <- setdiff(catch_cols, c("date", "month", "year"))

for (col in catch_cols_na) {
  snow_mask <- mod_swe_days[[col]] > 5
  obs_sm[snow_mask, (col) := NA_real_]
  mod_sm[snow_mask, (col) := NA_real_]
}

max(obs_sm$date)

res_sm <- bias_trends_comparison(obs_sm, mod_sm, "Soil Moisture - snow mask",
  obs_label = "ESA CCI", mod_label = "HERA-WB",
  y_label = "m\u00b3/m\u00b3", dt = "daily"
)

res_sm$combined

# ggsave(file.path(plots_dir, "seasonal_SM_snowmasked.png"), res_sm$combined,
#   width = 9, height = 9, dpi = 200
# )
# ggsave(file.path(plots_dir, "seasonal_SM_snowmasked_map.png"), res_sm$map,
#   width = 10, height = 8, dpi = 200
# )
ggsave(file.path(plots_dir, "Supplement_Figure_SM_snowmasked.png"), res_sm$combined,
  width = 14, height = 11, dpi = 400
)






cat("Loading catchments...\n")
catchments <- st_read(gpkg_path, quiet = TRUE)
catchments_3035 <- st_transform(catchments, 3035)
library(rnaturalearth)
basemap <- ne_countries(scale = "medium", returnclass = "sf") |>
  st_transform(3035)
bbox <- st_bbox(catchments_3035)


# No-snow days per catchment (scalar per catchment)
no_snow_count <- colSums(mod_swe_days[, ..catch_cols_swe] <= 1e-2 |
  is.na(mod_swe_days[, ..catch_cols_swe]))

catchments_3035$no_snow_days <- no_snow_count[
  match(as.numeric(catchments_3035$catch_id), as.numeric(names(no_snow_count)))
]
# As percentage of total days
catchments_3035$pct_snow_free <- 100 * catchments_3035$no_snow_days / nrow(mod_swe_days)


cats_map <- catchments_3035
# map of no_snow days
palet <- hcl.colors(11, palette = "Blues", rev = T)
ggplot() +
  geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.15) +
  # geom_sf(
  #     data = cats_map[!cats_map$significant, ],
  #     aes(fill = slope), color = NA, alpha = 0.3
  # ) +
  geom_sf(
    data = cats_map,
    aes(fill = pct_snow_free), color = "gray44", alpha = 0.9
  ) +
  scale_fill_gradientn(
    colors = palet,
    limits = c(0, 100),
    oob = scales::squish,
    name = "Percentage of no snow days"
  ) +
  coord_sf(
    xlim = c(bbox["xmin"], bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"]), expand = FALSE
  ) +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title.position = "top"
  )



# Root SM

agg_dir <- file.path(base_dir, "data", "aggregates")
rsm_path <- file.path(
  agg_dir, "surface_soil_moisture",
  "surface_soil_moisture_monthly_all_years.csv"
)
rsm_dt <- fread(rsm_path)
rsm_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
rsm_dt[, month := month(date)]
rsm_dt[, year := year(date)]

meta_cols_rsm <- c("month_idx", "period_start", "period_end", "date", "month", "year")
catch_cols_rsm <- setdiff(names(rsm_dt), meta_cols_rsm)

rsm_monthly <- rsm_dt[year >= start_year, c("year", "month", catch_cols_rsm), with = FALSE]

# Area weights
catch_ids_clean <- sub("^X", "", catch_cols_rsm)
area_vec <- catchments$residual_area_km2[match(
  catch_ids_clean, as.character(catchments$catch_id)
)]
area_vec[is.na(area_vec)] <- 1
weights <- area_vec / sum(area_vec, na.rm = TRUE)

# Continental mean seasonal cycle
seasonal_rsm <- sapply(1:12, function(m) {
  mat <- as.matrix(rsm_monthly[month == m, ..catch_cols_rsm])
  means <- colMeans(mat, na.rm = TRUE)
  sum(means * weights, na.rm = TRUE)
})



seasonal_df <- data.frame(
  month = c(1:12),
  value = c(seasonal_rsm),
  source = rep("Root SM", each = 12)
)

# Plot 1: Seasonal cycle
p_cycle <- ggplot(seasonal_df, aes(x = month, y = value, color = source)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  scale_y_continuous(limits = c(0, 0.6)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_cycle

# =============================================================================
# 4. Discharge: Observed vs LISFLOOD
# =============================================================================
cat("[4/4] Discharge seasonal comparison...\n")


# Read matching table
matching_path <- file.path(base_dir, "output", "station_catchment_matches.csv")
if (!file.exists(matching_path)) {
  stop(
    "Matching table not found: ", matching_path,
    "\nRun match_stations_catchments.R first."
  )
}
matches <- read.csv(matching_path, stringsAsFactors = FALSE)
length(unique(matches$catch_id[which(matches$match_type == "relaxed")]))


# Load homogenized discharge files (from homogenize_discharge.R)
dis_homog_dir <- file.path(base_dir, "output", "discharge_diego", "1.homogenized")

dis_obs_strict_path <- file.path(dis_homog_dir, "obs_monthly_strict.csv")
dis_obs_relaxed_path <- file.path(dis_homog_dir, "obs_monthly_relaxed.csv")
dis_sim_path <- file.path(dis_homog_dir, "sim_monthly.csv")

if (file.exists(dis_obs_strict_path) && file.exists(dis_sim_path)) {
  obs_dis <- fread(dis_obs_strict_path, header = TRUE)
  mod_dis <- fread(file.path(dis_homog_dir, "sim_monthly_strict.csv"), header = TRUE)

  # Parse dates (format: YYYY-MM)
  obs_dis[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
  obs_dis[, month := month(date)]
  mod_dis[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
  mod_dis[, month := month(date)]

  # transform into specific discharge
  # Transform into specific discharge (mm/day)
  # Q_specific = Q (m³/s) * 86400 (s/day) / (area_km² * 1e6 m²/km²) * 1000 (mm/m)
  #            = Q * 86.4 / area_km²

  catch_cols_dis <- setdiff(names(obs_dis), c("date", "month"))

  # Get upstream area for each matched catchment
  area_vec_dis <- catchments$area_km2[match(
    as.character(catch_cols_dis), as.character(catchments$catch_id)
  )]

  # Divide each column by its catchment area, convert to mm/day
  obs_dis[, (catch_cols_dis) := Map(`*`, .SD, 86.4 / area_vec_dis), .SDcols = catch_cols_dis]
  mod_dis[, (catch_cols_dis) := Map(`*`, .SD, 86.4 / area_vec_dis), .SDcols = catch_cols_dis]


  res_dis <- seasonal_comparison(obs_dis, mod_dis, "Discharge (strict)",
    obs_label = "Observed", mod_label = "HERA-WB",
    y_label = "m\u00b3/s/km\u00b2", color_obs = "#2ca02c"
  )

  ggsave(file.path(plots_dir, "seasonal_Discharge_strict.png"), res_dis$combined,
    width = 9, height = 9, dpi = 200
  )

  # Also do relaxed matches
  obs_dis_r <- fread(dis_obs_relaxed_path, header = TRUE)
  mod_dis_r <- fread(dis_sim_path, header = TRUE)
  obs_dis_r[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
  obs_dis_r[, month := month(date)]
  mod_dis_r[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]
  mod_dis_r[, month := month(date)]


  # transform into specific discharge
  # Transform into specific discharge (mm/day)
  # Q_specific = Q (m³/s) * 86400 (s/day) / (area_km² * 1e6 m²/km²) * 1000 (mm/m)
  #            = Q * 86.4 / area_km²

  catch_cols_dis <- setdiff(names(obs_dis_r), c("date", "month"))

  # Get upstream area for each matched catchment
  area_vec_dis <- catchments$area_km2[match(
    as.character(catch_cols_dis), as.character(catchments$catch_id)
  )]

  # Divide each column by its catchment area, convert to mm/day
  obs_dis_r[, (catch_cols_dis) := Map(`*`, .SD, 86.4 / area_vec_dis), .SDcols = catch_cols_dis]
  mod_dis_r[, (catch_cols_dis) := Map(`*`, .SD, 86.4 / area_vec_dis), .SDcols = catch_cols_dis]
  res_dis_relaxed <- bias_trends_comparison(obs_dis_r, mod_dis_r, "Discharge (all)",
    obs_label = "Observed", mod_label = "HERA-WB",
    y_label = "m\u00b3/s/km\u00b2"
  )

  ggsave(file.path(plots_dir, "Supplement_Figure_Discharge.png"), res_dis_relaxed$combined,
    width = 14, height = 11, dpi = 400
  )
} else {
  cat("  Homogenized discharge files not found. Run homogenize_discharge.R first.\n")
  cat("  Expected:", dis_obs_strict_path, "\n")
  res_dis <- NULL
}

# =============================================================================
# Combined 4-variable figure (seasonal cycles only)
# =============================================================================
cat("Composing combined seasonal cycle figure...\n")

panels_cycle <- list(res_aet$cycle, res_swe$cycle, res_sm$cycle)
if (!is.null(res_dis)) panels_cycle[[4]] <- res_dis$cycle

fig_all_cycles <- plot_grid(plotlist = panels_cycle, ncol = 2, align = "hv")
ggsave(file.path(plots_dir, "seasonal_cycles_all_variables.png"),
  fig_all_cycles,
  width = 14, height = 10, dpi = 200
)

# Combined bias figure
panels_bias <- list(res_aet$bias, res_swe$bias, res_sm$bias)
if (!is.null(res_dis)) panels_bias[[4]] <- res_dis$bias

fig_all_bias <- plot_grid(plotlist = panels_bias, ncol = 2, align = "hv")
ggsave(file.path(plots_dir, "seasonal_bias_all_variables.png"),
  fig_all_bias,
  width = 14, height = 10, dpi = 400
)

cat("\nDone! Seasonal comparison figures saved to:", plots_dir, "\n")
