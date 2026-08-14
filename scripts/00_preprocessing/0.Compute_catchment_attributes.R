# =============================================================================
# Compute static catchment attributes and save as catchment_attributes.csv
#
# Attributes:
#   - climate class (Köppen-Geiger majority)
#   - mean elevation (m)
#   - elevation standard deviation (m)
#   - mean gradient (m/m)
#   - upstream area (km²)
#   - residual upstream area (km²)
#   - outlet elevation (m) [elevation at the outlet pixel]
#   - total population (sum within catchment)
#   - ksat1_forest (mm/day, top soil layer, forest)
#   - ksat1_other (mm/day, top soil layer, non-forest)
#   - soildepth1a_forest (layer 1a depth, forest)
#   - soildepth1a_other (layer 1a depth, non-forest)
#   - soildepth1b_forest (layer 1b depth, forest)
#   - soildepth1b_other (layer 1b depth, non-forest)
#   - soildepth2_forest (layer 2 depth, forest)
#   - soildepth2_other (layer 2 depth, non-forest)
#   - longitude, latitude (centroid)
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(data.table)
library(sf)
library(terra)
library(exactextractr)
library(dplyr)
library(ncdf4)


outletopen <- function(dir, outletname, nrspace = rep(NA, 5)) {
  ncbassin <- paste0(dir, "/", outletname, ".nc")
  ncb <- nc_open(ncbassin)
  name.vb <- names(ncb[["var"]])
  namev <- name.vb[1]
  name.lon <- "lon"
  name.lat <- "lat"
  if (!is.na(nrspace[1])) {
    start <- as.numeric(nrspace[c(2, 4)])
    count <- as.numeric(nrspace[c(3, 5)]) - start + 1
  } else {
    londat <- ncvar_get(ncb, name.lon)
    llo <- length(londat)
    latdat <- ncvar_get(ncb, name.lat)
    lla <- length(latdat)
    start <- c(1, 1)
    count <- c(llo, lla)
  }

  londat <- ncvar_get(ncb, name.lon, start = start[1], count = count[1])
  llo <- length(londat)
  latdat <- ncvar_get(ncb, name.lat, start = start[2], count = count[2])
  lla <- length(latdat)
  outlets <- ncvar_get(ncb, namev, start = start, count = count)
  outlets <- as.vector(outlets)
  outll <- expand.grid(londat, latdat)
  lonlatloop <- expand.grid(c(1:llo), c(1:lla))
  outll$idlo <- lonlatloop$Var1
  outll$idla <- lonlatloop$Var2

  outll <- outll[which(!is.na(outlets)), ]
  outlets <- outlets[which(!is.na(outlets))]
  outll <- data.frame(outlets, outll)
  return(outll)
}

DemOpen <- function(dir, outletname, Sloc_final) {
  ncbassin <- paste0(dir, outletname)
  ncb <- nc_open(ncbassin)
  name.vb <- names(ncb[["var"]])
  namev <- name.vb[2]
  # time <- ncvar_get(ncb,"time")

  # timestamp corretion
  name.lon <- "lon"
  name.lat <- "lat"
  londat <- ncvar_get(ncb, name.lon)
  llo <- length(londat)
  latdat <- ncvar_get(ncb, name.lat)
  lla <- length(latdat)
  start <- c(1, 1)
  count <- c(llo, lla)


  londat <- ncvar_get(ncb, name.lon, start = start[1], count = count[1])
  llo <- length(londat)
  latdat <- ncvar_get(ncb, name.lat, start = start[2], count = count[2])
  lla <- length(latdat)
  outlets <- ncvar_get(ncb, namev, start = start, count = count)
  outlets <- as.vector(outlets)
  outll <- expand.grid(londat, latdat)
  lonlatloop <- expand.grid(c(1:llo), c(1:lla))
  outll$elev <- outlets
  outll$idlo <- lonlatloop$Var1
  outll$idla <- lonlatloop$Var2

  # outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outll$latlong <- paste(round(outll$Var1, 4), round(outll$Var2, 4), sep = " ")
  outfinal <- inner_join(outll, Sloc_final, by = "latlong")
  return(outfinal)
}

# --- Paths --------------------------------------------------------------------
source("config/paths.R")
attr_out_dir <- file.path(data_dir, "attributes")
dir.create(attr_out_dir, showWarnings = FALSE, recursive = TRUE)

# Raster paths
path_elvstd <- file.path(hydro_dir, "mapscal", "elvstd_European_01min.nc")
path_gradient <- file.path(hydro_dir, "mapscal", "gradient_European_01min.nc")
path_pop <- file.path(hydro_dir, "mapscal", "population_European_01min.nc")
path_ksat1_f <- file.path(hydro_dir, "mapscal", "soilhyd", "ksat1_f_European_01min.nc")
path_ksat1_o <- file.path(hydro_dir, "mapscal", "soilhyd", "ksat1_o_European_01min.nc")
path_outlets <- hydro_dir

# Soil depth paths (layer 1a, 1b, 2 × forest/other)
path_soildepth1_f <- file.path(flood_dir, "mapscal", "table2maps", "soildepth1_f_European_01min.nc")
path_soildepth1_o <- file.path(flood_dir, "mapscal", "table2maps", "soildepth1_o_European_01min.nc")
path_soildepth2_f <- file.path(flood_dir, "mapscal", "table2maps", "soildepth2_f_European_01min.nc")
path_soildepth2_o <- file.path(flood_dir, "mapscal", "table2maps", "soildepth2_o_European_01min.nc")
path_soildepth3_f <- file.path(flood_dir, "mapscal", "table2maps", "soildepth3_f_European_01min.nc")
path_soildepth3_o <- file.path(flood_dir, "mapscal", "table2maps", "soildepth3_o_European_01min.nc")
# --- Load catchments ----------------------------------------------------------
cat("Loading catchments...\n")
catchments <- st_read(gpkg_path, quiet = TRUE)
n_catch <- nrow(catchments)
cat("  Catchments:", n_catch, "\n")

# Ensure WGS84 for extraction
catch_wgs <- st_transform(catchments, 4326)

# --- Compute centroids --------------------------------------------------------
centroids <- st_centroid(catch_wgs)
coords <- st_coordinates(centroids)

# --- Extract raster attributes ------------------------------------------------
cat("Extracting raster attributes...\n")

# Mean elevation
cat("  Mean elevation...\n")
r_dem <- rast(path_dem)
catch_dem <- st_transform(catch_wgs, crs(r_dem))
elev_mean <- exact_extract(r_dem, catch_dem, "mean")

elev_std <- exact_extract(r_dem, catch_dem, "stdev")

# Mean gradient
cat("  Mean gradient...\n")
r_grad <- rast(path_gradient)
catch_grad <- st_transform(catch_wgs, crs(r_grad))
gradient_mean <- exact_extract(r_grad, catch_grad, "mean")


# Ksat top layer - forest
cat("  Ksat1 (forest)...\n")
r_ksat_f <- rast(path_ksat1_f)
catch_ksat <- st_transform(catch_wgs, crs(r_ksat_f))
ksat1_forest <- exact_extract(r_ksat_f, catch_ksat, "mean")

# Ksat top layer - other
cat("  Ksat1 (other)...\n")
r_ksat_o <- rast(path_ksat1_o)
ksat1_other <- exact_extract(r_ksat_o, catch_ksat, "mean")

# Soil depth layer 1a (top) - forest
cat("  Soil depth layer 1a (forest)...\n")
r_sd1_f <- rast(path_soildepth1_f)
catch_sd <- st_transform(catch_wgs, crs(r_sd1_f))
soildepth1a_forest <- exact_extract(r_sd1_f, catch_sd, "mean")

# Soil depth layer 1a (top) - other
cat("  Soil depth layer 1a (other)...\n")
r_sd1_o <- rast(path_soildepth1_o)
soildepth1a_other <- exact_extract(r_sd1_o, catch_sd, "mean")

# Soil depth layer 1b - forest
cat("  Soil depth layer 1b (forest)...\n")
r_sd2_f <- rast(path_soildepth2_f)
soildepth1b_forest <- exact_extract(r_sd2_f, catch_sd, "mean")

# Soil depth layer 1b - other
cat("  Soil depth layer 1b (other)...\n")
r_sd2_o <- rast(path_soildepth2_o)
soildepth1b_other <- exact_extract(r_sd2_o, catch_sd, "mean")

# Soil depth layer 2 - forest
cat("  Soil depth layer 2 (forest)...\n")
r_sd3_f <- rast(path_soildepth3_f)
soildepth2_forest <- exact_extract(r_sd3_f, catch_sd, "mean")

# Soil depth layer 2 - other
cat("  Soil depth layer 2 (other)...\n")
r_sd3_o <- rast(path_soildepth3_o)
soildepth2_other <- exact_extract(r_sd3_o, catch_sd, "mean")

# Climate class (majority)
cat("  Climate class...\n")
clim_lookup <- c(
  setNames(rep("Tropical", 3), as.character(1:3)),
  setNames(rep("Arid", 4), as.character(4:7)),
  setNames(rep("Temperate", 9), as.character(8:16)),
  setNames(rep("Cold", 12), as.character(17:28)),
  setNames(rep("Polar", 2), as.character(29:30))
)

r_clim <- rast(path_clim)
catch_clim <- st_transform(catch_wgs, crs(r_clim))
clim_class <- exact_extract(r_clim, catch_clim,
  fun = function(values, coverage_fractions) {
    if (all(is.na(values))) {
      return(NA_character_)
    }
    maj <- as.character(names(which.max(table(values[!is.na(values)]))))
    unname(clim_lookup[maj])
  }
)

# Outlet elevation (elevation at centroid/outlet point)
cat("  Outlet elevation...\n")
outlets <- outletopen(path_outlets, "outletsv8_hybas07_01min")

outlets <- as.data.frame(outlets)
outlets$latlong <- paste(round(outlets$Var1, 4), round(outlets$Var2, 4), sep = " ")
outlets_elev <- DemOpen(dir = file.path(base_dir, "data"), outletname = "/dem.nc", outlets)

outlets_elev2 <- outlets_elev[which(!is.na(match(outlets_elev$outlets, as.numeric(catchments$catch_id)))), ]

# --- Build attributes table ---------------------------------------------------
cat("Building attributes table...\n")

attr_dt <- data.table(
  catch_id           = catchments$catch_id,
  longitude          = round(coords[, 1], 5),
  latitude           = round(coords[, 2], 5),
  area_km2           = catchments$area_km2,
  residual_area_km2  = catchments$residual_area_km2,
  climate_class      = clim_class,
  elev_mean_m        = round(elev_mean, 1),
  elev_std_m         = round(elev_std, 1),
  outlet_elev_m      = round(outlets_elev2$elev, 1),
  gradient_mean      = round(gradient_mean, 5),
  ksat1_forest_mmday = round(ksat1_forest, 2),
  ksat1_other_mmday  = round(ksat1_other, 2),
  soildepth1a_forest = round(soildepth1a_forest, 2),
  soildepth1a_other  = round(soildepth1a_other, 2),
  soildepth1b_forest = round(soildepth1b_forest, 2),
  soildepth1b_other  = round(soildepth1b_other, 2),
  soildepth2_forest  = round(soildepth2_forest, 2),
  soildepth2_other   = round(soildepth2_other, 2)
)

# --- Save ---------------------------------------------------------------------
out_path <- file.path(out_dir, "catchment_attributes.csv")
fwrite(attr_dt, out_path)
cat("\nSaved:", out_path, "\n")
cat("  Rows:", nrow(attr_dt), "| Columns:", ncol(attr_dt), "\n")

catchments2 <- catchments[, c(1, 3, 4, 7, 8)]
attr_sf <- inner_join(catchments2, attr_dt, by = c("catch_id"))

st_write(attr_sf, file.path(out_dir, "catchment_attributes.gpkg"), delete_dsn = TRUE)

# --- Summary ------------------------------------------------------------------
cat("\n--- Summary ---\n")
cat("  Elevation range:", range(attr_dt$elev_mean_m, na.rm = TRUE), "m\n")
cat("  Area range:", range(attr_dt$residual_area_km2, na.rm = TRUE), "km²\n")
cat("  Climate classes:", paste(unique(na.omit(attr_dt$climate_class)), collapse = ", "), "\n")
cat("  Ksat1 (forest) range:", range(attr_dt$ksat1_forest_mmday, na.rm = TRUE), "mm/day\n")
cat("  Ksat1 (other) range:", range(attr_dt$ksat1_other_mmday, na.rm = TRUE), "mm/day\n")
cat("  Soil depth 1a (forest) range:", range(attr_dt$soildepth1a_forest, na.rm = TRUE), "\n")
cat("  Soil depth 1a (other) range:", range(attr_dt$soildepth1a_other, na.rm = TRUE), "\n")
cat("  Soil depth 1b (forest) range:", range(attr_dt$soildepth1b_forest, na.rm = TRUE), "\n")
cat("  Soil depth 1b (other) range:", range(attr_dt$soildepth1b_other, na.rm = TRUE), "\n")
cat("  Soil depth 2 (forest) range:", range(attr_dt$soildepth2_forest, na.rm = TRUE), "\n")
cat("  Soil depth 2 (other) range:", range(attr_dt$soildepth2_other, na.rm = TRUE), "\n")
cat("  Outlet elevation:", range(attr_dt$outlet_elev_m, na.rm = TRUE), "\n")

cat("\nDone!\n")
