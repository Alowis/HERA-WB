# =============================================================================
# Temporal evolution plots
# 1. Monthly infiltration for a single catchment
# 2. Monthly SWE for a single catchment
# 3. Continental area-weighted variable (climate stripes style)
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(ggplot2)
library(sf)
library(dplyr)
library(lubridate)
library(scales)
library(terra)
library(exactextractr)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
agg_dir <- file.path(base_dir, "data", "aggregates")
tss_dir <- file.path(base_dir, "data", "tss_postprocess")
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
out_dir <- file.path(base_dir, "output", "temporal_evolution")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- User choices -------------------------------------------------------------
# Catchment IDs to plot (change as needed)
catch_infiltration <- "224577" # Paris
catch_infiltration_name <- "Paris"
catch_infiltration <- "184168" # London
catch_infiltration <- "307920" # Avignon
catch_infiltration_name <- "Avignon"
catch_infiltration <- "191110" # bruxellss
catch_infiltration <- "162248" # birmimgham
catch_infiltration <- "292302" # Milano
catch_swe <- "290666" # catchment for SWE plot
# Variable for the continental stripes (must match a folder in data/aggregates/)
stripes_variable <- "snow_water_equivalent"
# Land use sealed fraction
landuse_dir <- "D:/tilloal/Documents/06_Floodrivers/landuse/"
catch_sealed <- catch_infiltration # catchments for sealed fraction plot

# =============================================================================
# 1. MONTHLY INFILTRATION RATIO — SINGLE CATCHMENT
#    Ratio = infiltration / (rainfall + snowfall)
# =============================================================================
cat("[1/4] Loading infiltration ratio data for catchment", catch_infiltration, "...\n")

# The data is daily in tss_postprocess (very large files).
# Use fread with select to read only the needed column.
# Note: TSS postprocess files have catchment IDs as the header row (first row),
# no timestamp column — rows are sequential timesteps.
inf_path <- file.path(tss_dir, "infUpsX_nested_1951_2020.csv")
rain_path <- file.path(tss_dir, "rainUpsX_nested_1951_2020.csv")
snow_path <- file.path(tss_dir, "snowUpsX_nested_1951_2020.csv")

# Read header to find column names (force header=TRUE since first row
# contains numeric catchment IDs that fread may not auto-detect as header)
inf_header <- fread(inf_path, nrows = 0, header = TRUE)
inf_cols <- names(inf_header)

# Find target catchment column (no timestamp column in these files)
if (catch_infiltration %in% inf_cols) {
  col_to_read <- catch_infiltration
} else {
  # Try with X prefix (R's check.names default)
  x_catch <- paste0("X", catch_infiltration)
  if (x_catch %in% inf_cols) {
    col_to_read <- x_catch
  } else {
    stop("Catchment ", catch_infiltration, " not found in infiltration CSV.")
  }
}

cat("  Reading column:", col_to_read, "\n")
inf_dt <- fread(inf_path, select = col_to_read, header = TRUE)
rain_dt <- fread(rain_path, select = col_to_read, header = TRUE)
snow_dt <- fread(snow_path, select = col_to_read, header = TRUE)

# Combine into single data.table
inf_dt[, infiltration := inf_dt[[1]]]
inf_dt[, precipitation := rain_dt[[1]] + snow_dt[[1]]]
inf_dt[, value := fifelse(precipitation > 0, infiltration / precipitation, NA_real_)]

# Assign date: daily from 1951-01-01 to 2020-12-31
# LISFLOOD TSS files are 6-hourly (4 steps/day) but the nested/CC versions
# are already daily aggregates. Check row count to determine timestep.
n_rows <- nrow(inf_dt)
n_days <- as.integer(as.Date("2020-12-31") - as.Date("1951-01-01")) + 1

if (n_rows == n_days) {
  inf_dt[, date := seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")]
} else if (n_rows == n_days * 4) {
  # 6-hourly: aggregate to daily first
  inf_dt[, day_idx := rep(seq_len(n_days), each = 4)]
  date_seq <- seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")
  inf_dt[, date := date_seq[day_idx]]
  precip_dt <- inf_dt[, .(precip = mean(precipitation, na.rm = TRUE)), by = date]
  inf_dt <- inf_dt[, .(infil = mean(infiltration, na.rm = TRUE)), by = date]
} else {
  # Monthly data from tss_postprocess (840 months)
  n_months <- length(seq.Date(as.Date("1951-01-01"), as.Date("2020-12-01"), by = "month"))
  if (n_rows == n_months) {
    inf_dt[, date := seq.Date(as.Date("1951-01-01"), as.Date("2020-12-01"), by = "month")]
  } else {
    warning("Unexpected row count (", n_rows, "). Assuming daily with possible offset.")
    inf_dt[, date := seq.Date(as.Date("1951-01-01"), length.out = n_rows, by = "day")]
  }
}

# Aggregate to monthly means
inf_monthly <- inf_dt[, .(
  infil = mean(infil, na.rm = TRUE)
), by = .(year = year(date), month = month(date))]
inf_monthly[, date := as.Date(paste(year, month, "15", sep = "-"))]

precip_monthly <- precip_dt[, .(
  infil = mean(precip, na.rm = TRUE)
), by = .(year = year(date), month = month(date))]
precip_monthly[, date := as.Date(paste(year, month, "15", sep = "-"))]

inf_ratio_monthly <- inf_monthly
inf_ratio_monthly[, inf_ratio := inf_monthly[[3]] / precip_monthly[[3]]]
inf_ratio_monthly[, inf_ratio := frollmean(inf_ratio, 3, align = "center")]
# Plot
p1 <- ggplot(inf_ratio_monthly, aes(x = date, y = inf_ratio)) +
  geom_line(color = "steelblue", linewidth = 0.4) +
  geom_smooth(
    method = "loess", span = 0.5, se = T,
    color = "darkblue", linewidth = 0.8
  ) +
  labs(
    title = paste("Monthly infiltration ratio — Catchment", catch_infiltration),
    x = NULL,
    y = "Infiltration / Precipitation (-)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p1

ggsave(file.path(out_dir, paste0("infil_ratio_monthly_", catch_infiltration_name, ".png")),
  p1,
  width = 10, height = 4, dpi = 150
)
cat("  -> Saved infiltration ratio plot.\n")



# =============================================================================
# 4. SEALED LAND USE FRACTION — EVOLUTION PER CATCHMENT
# =============================================================================
cat("[4/4] Extracting sealed fraction evolution for selected catchments...\n")

# List available fracsealed NetCDF files (one per year)
sealed_files <- list.files(
  landuse_dir,
  pattern = "^fracsealed_European_01min_\\d{4}\\.nc$",
  full.names = TRUE
)

if (length(sealed_files) == 0) {
  stop("No fracsealed NetCDF files found in: ", landuse_dir)
}

# Extract year from filename
sealed_years <- as.integer(sub(
  ".*fracsealed_European_01min_(\\d{4})\\.nc$", "\\1", sealed_files
))
sealed_order <- order(sealed_years)
sealed_files <- sealed_files[sealed_order]
sealed_years <- sealed_years[sealed_order]

# Load catchments if not already loaded
if (!exists("catchments")) {
  catchments <- st_read(gpkg_path, quiet = TRUE)
}
# Ensure WGS84 for extraction
if (!isTRUE(st_crs(catchments) == st_crs(4326))) {
  catchments <- st_transform(catchments, 4326)
}

# Subset to selected catchments
catch_idx <- which(as.character(catchments$catch_id) %in% catch_sealed)
if (length(catch_idx) == 0) {
  stop(
    "None of the selected catchments found in the GeoPackage: ",
    paste(catch_sealed, collapse = ", ")
  )
}
catch_sub <- catchments[catch_idx, ]

# Extract area-weighted mean sealed fraction per catchment per year
sealed_results <- data.frame(
  year = rep(sealed_years, each = nrow(catch_sub)),
  catch_id = rep(as.character(catch_sub$catch_id), times = length(sealed_years)),
  frac_sealed = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(sealed_files)) {
  cat(
    "  Extracting year", sealed_years[i],
    "(", i, "/", length(sealed_files), ")\n"
  )

  r <- tryCatch(
    terra::rast(sealed_files[i]),
    error = function(e) {
      warning("Could not read: ", basename(sealed_files[i]), " - skipping")
      NULL
    }
  )
  if (is.null(r)) next

  # Use first layer if multi-band
  if (nlyr(r) > 1) r <- r[[1]]

  # Extract area-weighted mean
  means <- exact_extract(r, catch_sub, "mean")

  row_idx <- ((i - 1) * nrow(catch_sub) + 1):(i * nrow(catch_sub))
  sealed_results$frac_sealed[row_idx] <- means
}

# Remove failed extractions
sealed_results <- sealed_results[!is.na(sealed_results$frac_sealed), ]

# Plot: line per catchment
p4 <- ggplot(sealed_results, aes(
  x = year, y = frac_sealed,
  color = catch_id
)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.2, alpha = 0.6) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  labs(
    title = "Sealed land use fraction — evolution per catchment",
    x = NULL,
    y = "Sealed fraction",
    color = "Catchment ID"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )


p4


ggsave(file.path(out_dir, paste0("sealed_fraction_", catch_infiltration_name, ".png")),
  p4,
  width = 10, height = 5, dpi = 150
)
# =============================================================================
# 2. MONTHLY SWE — SINGLE CATCHMENT
# =============================================================================
cat("[2/4] Loading SWE data for catchment", catch_swe, "...\n")

catch_swe <- "290666" # ticino
catch_swe <- "303662" # ardeche
catch_swe <- "289805" # aosta
catch_swe <- "325150" #Aquila

swe_path <- file.path(
  agg_dir, "snow_water_equivalent",
  "snow_water_equivalent_monthly_all_years.csv"
)
swe_dt <- fread(swe_path)

# Parse date
swe_dt[, date := as.Date(period_start)]

# Extract catchment column
swe_col <- catch_swe
if (!swe_col %in% names(swe_dt)) {
  swe_col <- paste0("X", catch_swe)
}
if (!swe_col %in% names(swe_dt)) {
  stop("Catchment ", catch_swe, " not found in SWE aggregates.")
}

swe_catch <- swe_dt[, .(date, swe = get(swe_col))]

# Keep only first 6 months (Jan–Jun): snow accumulation & melt season
swe_catch <- swe_catch[month(date) <= 6]

# Compute annual means for stripes
swe_annual <- swe_catch[, .(
  annual_sum = sum(swe, na.rm = TRUE)
), by = .(year = year(date))]

# Plot
p2 <- ggplot(swe_annual, aes(x = year, y = annual_sum)) +
  geom_line(color = "dodgerblue3", linewidth = 0.3) +
  geom_smooth(
    method = "lm", se = T,
    color = "navy", linewidth = 0.8
  ) +
  labs(
    title = paste("Monthly mean SWE (Jan–Jun) — Catchment", catch_swe),
    x = NULL,
    y = "Snow Water Equivalent (mm)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p2

ggsave(file.path(out_dir, paste0("swe_monthly_", catch_swe, ".png")),
  p2,
  width = 10, height = 4, dpi = 150
)
cat("  -> Saved SWE plot.\n")

# =============================================================================
# 3. CONTINENTAL CLIMATE STRIPES — AREA-WEIGHTED
# =============================================================================
cat("[3/4] Computing continental area-weighted stripes for:", stripes_variable, "\n")

stripes_variable<-"ActEvapo"
limi=c(-10, 10)
# Load the monthly aggregates for the chosen variable
stripes_path <- file.path(
  agg_dir, stripes_variable,
  paste0(stripes_variable, "_monthly_all_years.csv")
)
if (!file.exists(stripes_path)) {
  stop("Aggregates file not found: ", stripes_path)
}
stripes_dt <- fread(stripes_path)
stripes_dt[, date := as.Date(period_end)]
if(stripes_variable=="snow_water_equivalent") stripes_dt <- stripes_dt[month(date) <= 6]
# Identify catchment columns (numeric IDs)
meta_cols <- c("month_idx", "period_start", "period_end", "date")
catch_cols <- setdiff(names(stripes_dt), meta_cols)

# Load catchments to get upstream area for weighting
cat("  Loading catchment areas...\n")
catchments <- st_read(gpkg_path, quiet = TRUE)

common_ids <- intersect(catch_cols, as.character(as.numeric(catchments$catch_id)))

if (length(common_ids) < 10) {
  warning(
    "Few matching catchments (", length(common_ids),
    "). Check catch_id format."
  )
}

# Build weight vector (normalized)
area_vec <- catchments$residual_area_km2[match(common_ids, as.character(as.numeric(catchments$catch_id)))]
weights <- area_vec / sum(area_vec, na.rm = TRUE)

which(is.na(weights))
# Compute area-weighted mean per timestep
cat("  Computing area-weighted mean (", length(common_ids), " catchments)...\n")
mat <- (as.matrix(stripes_dt[, ..common_ids]))
is.na(mat)
stripes_dt[, weighted_mean := as.numeric(mat %*% weights)]



# Compute annual means for stripes
stripes_annual <- stripes_dt[, .(
  annual = max(weighted_mean, na.rm = TRUE)
), by = .(year = year(date))]

# Climate stripes plot
# Color each year bar by its anomaly relative to the long-term mean
long_term_mean <- mean(stripes_annual$annual, na.rm = TRUE)
stripes_annual[, anomaly := annual - long_term_mean]

palet=hcl.colors(9, palette = "RdBu", alpha = NULL, rev = T, fixup = TRUE)

p3 <- ggplot(stripes_annual, aes(x = year, y = 1, fill = anomaly)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_gradientn(
    colors = palet, limits = limi, oob = squish,
    breaks = seq(limi[1],limi[2],limi[2]/5) ,
    name = paste0(stripes_variable, "\nanomaly (mm)")
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(1955, 2020, by = 10)
  ) +
  labs(
    title = paste(
      "Continental", gsub("_", " ", stripes_variable),
      "— area-weighted max anomaly (1951–2020)"
    ),
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.text.x = element_text(size = 10, margin = margin(t = 4)),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

p3
ggsave(file.path(out_dir, paste0("stripes_AMAX_", stripes_variable, ".png")),
  p3,
  width = 12, height = 3, dpi = 150
)
cat("  -> Saved climate stripes plot.\n")

cat("  -> Saved sealed fraction plot.\n")

cat("\nAll plots saved to:", out_dir, "\n")
