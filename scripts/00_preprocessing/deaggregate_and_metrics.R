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
source(paste0(base_dir, "RegimeShift_codes/functions_regime.R"))

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
plot_dir <- paste0(base_dir, "plots/")

# --- Analysis parameters ---
year_range <- 1951:2020
dt_step <- 4 # 6-hourly time steps per day

# --- Variables for batch deaggregation ---
batch_vars <- c(
    "rainUpsX", "snowMeltUpsX", "qUzUpsX", "qLZUpsX",
    "surfaceRunoffUpsX", "infUpsX", "ActEvapo", "disWin",
    "percUZLZUpsX", "dSubToUzUpsX", "prefFlowUpsX", "lossUpsX"
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
# SECTION 5: Deaggregate Rain and Snow, Compute Precipitation & Snow Fraction
# =============================================================================
message("Loading and deaggregating Rain ...")
Rain_raw <- fread(rain_file, header = TRUE)
time <- Rain_raw$V1
timeStampX <- time[order(time)]
Rain_raw <- Rain_raw[order(time), .SD, .SDcols = matcol]

message("Loading and deaggregating Snow ...")
Snow_raw <- fread(snow_file, header = TRUE)
Snow_raw <- Snow_raw[order(time), .SD, .SDcols = matcol]

Rain_res <- deaggregate_to_residual(
    Rain_raw, imm_children_list,
    area_lookup, resarea_lookup
)
Snow_res <- deaggregate_to_residual(
    Snow_raw, imm_children_list,
    area_lookup, resarea_lookup
)
Snow_res <- Snow_res[, lapply(.SD, function(x) pmax(0, x))]

# Derived fields
Precipitation_res <- Rain_res + Snow_res
Snowfraction_res <- Snow_res / Precipitation_res # NaN where P = 0 is acceptable

rm(Rain_raw, Snow_raw)
gc()

# =============================================================================
# SECTION 6: Compute Yearly Snow Fraction
# =============================================================================
message("Computing yearly snow fraction ...")
dsel <- hour(timeStampX)

SaveAgg_SF <- data.table(year = year_range)
for (col_name in cnames) {
    Trun <- preprocess_frac(
        time = timeStampX, input_var = Snowfraction_res,
        col_name, dsel, dt_step
    )
    Yagg <- process_frac(Trun)
    SaveAgg_SF[, (col_name) := Yagg$val]
}

save(SaveAgg_SF, file = paste0(base_dir, "SnowFraction_nested2.Rdata"))

# =============================================================================
# SECTION 7: Load AET and PET, Deaggregate, Compute Yearly Aggregates
# =============================================================================
message("Loading and deaggregating AET ...")
ActEvapo_raw <- fread(aet_file, header = TRUE)
colnames(ActEvapo_raw)[1] <- "time"
time <- ActEvapo_raw$time
timeStampX <- time[order(time)]
ActEvapo_raw <- ActEvapo_raw[order(time), .SD, .SDcols = matcol]

ActEvapo_res <- deaggregate_to_residual(
    ActEvapo_raw, imm_children_list,
    area_lookup, resarea_lookup
)
ActEvapo_res <- ActEvapo_res[, lapply(.SD, function(x) pmax(0, x))]
rm(ActEvapo_raw)
gc()

message("Loading and deaggregating PET ...")
PEvapo_raw <- fread(pet_file, header = TRUE)
PEvapo_raw <- PEvapo_raw[order(time), .SD, .SDcols = matcol]

PEvapo_res <- deaggregate_to_residual(
    PEvapo_raw, imm_children_list,
    area_lookup, resarea_lookup
)
PEvapo_res <- PEvapo_res[, lapply(.SD, function(x) pmax(0, x))]
rm(PEvapo_raw)
gc()

# --- Yearly aggregates ---
message("Computing yearly aggregates for AET, PET, Precipitation ...")
SaveE0_res <- data.table(year = year_range)
SaveAET_res <- data.table(year = year_range)
SavePrecip_res <- NULL
SeasonPrecip_res <- NULL

for (col_name in cnames) {
    AET <- preprocess_in(
        time = timeStampX, input_var = ActEvapo_res,
        col_name, dsel, dt_step
    )
    Precip <- preprocess_in(
        time = timeStampX, input_var = Precipitation_res,
        col_name, dsel, dt_step
    )
    TE0 <- preprocess_in(
        time = timeStampX, input_var = PEvapo_res,
        col_name, dsel, dt_step
    )

    E0Acc <- process_data(TE0)
    AETAcc <- process_data(AET)
    PrecipAcc <- process_precip(Precip)

    PrecipAcc$seaonal$cat <- rep(col_name, 2)
    PrecipAcc$yagg$cat <- rep(col_name, length(PrecipAcc$yagg$year))

    SavePrecip_res <- rbind(SavePrecip_res, PrecipAcc$yagg)
    SeasonPrecip_res <- rbind(SeasonPrecip_res, PrecipAcc$seaonal)
    SaveAET_res[, (col_name) := AETAcc$val]
    SaveE0_res[, (col_name) := E0Acc$val]
}

save(SavePrecip_res, file = paste0(base_dir, "PrecipitationMetrics_nested_f.Rdata"))
save(SeasonPrecip_res, file = paste0(base_dir, "PrecipitationSeason_nested_f.Rdata"))
save(SaveAET_res, file = paste0(base_dir, "AETMetrics_nested_f.Rdata"))

rm(ActEvapo_res, PEvapo_res, Precipitation_res, Snowfraction_res)
gc()

# =============================================================================
# SECTION 8: Compute Derived Metrics (Aridity, Gini, Seasonality)
# =============================================================================
message("Computing derived metrics ...")

# --- Mean annual AET ---
ys1 <- year_range
SaveAET_df <- as.data.frame(SaveAET_res)
my_aet <- which(SaveAET_df$year %in% ys1)
KAET_means <- colMeans(SaveAET_df[my_aet, -1], na.rm = TRUE)

# --- Mean annual precipitation and aridity index (AET/P) ---
myp <- which(SavePrecip_res$year %in% ys1)
KeepP <- SavePrecip_res[myp, c("year", "val.sum", "cat")]

KeepPag <- stats::aggregate(
    list(pmean = KeepP$val.sum),
    by = list(catch = KeepP$cat),
    FUN = function(x) mean(x, na.rm = TRUE)
)
KeepPag <- do.call(data.frame, KeepPag)

kpe <- match(KeepPag$catch, cnames)
KeepPag$aet <- KAET_means[kpe]
KeepPag$AET.P <- KeepPag$aet / KeepPag$pmean

# --- Gini concentration index ---
KeepG <- SavePrecip_res[myp, c("year", "val.gini", "cat")]

KeepGag <- stats::aggregate(
    list(px = KeepG$val.gini),
    by = list(catch = KeepG$cat),
    FUN = function(x) mean(x, na.rm = TRUE)
)
KeepGag <- do.call(data.frame, KeepGag)

# --- Precipitation seasonality index and change ---
ls_sp <- nrow(SeasonPrecip_res)
SeasonPrecipSI_res <- SeasonPrecip_res[-seq(1, ls_sp, 2), ]
SeasonPrecipSI_res$SIchange <- SeasonPrecipSI_res$Seasonality_2 -
    SeasonPrecipSI_res$Seasonality_1

# =============================================================================
# SECTION 9: Enrich Catchments with Climate Zones and Elevation
# =============================================================================
message("Enriching catchments with climate and elevation ...")

clim_lookup <- c(
    setNames(rep("Tropical", 3), as.character(1:3)),
    setNames(rep("Arid", 4), as.character(4:7)),
    setNames(rep("Temperate", 9), as.character(8:16)),
    setNames(rep("Cold", 12), as.character(17:28)),
    setNames(rep("Polar", 2), as.character(29:30))
)

r_clim <- terra::rast(path_clim)
catchments_gpkg$clim_class <- exactextractr::exact_extract(
    r_clim, sf::st_transform(catchments_gpkg, terra::crs(r_clim)),
    fun = function(values, coverage_fractions) {
        if (all(is.na(values))) {
            return(NA_character_)
        }
        maj <- as.character(names(which.max(table(values[!is.na(values)]))))
        unname(clim_lookup[maj])
    }
)

r_dem <- terra::rast(path_dem)
catchments_gpkg$elev_m <- exactextractr::exact_extract(
    r_dem, sf::st_transform(catchments_gpkg, terra::crs(r_dem)),
    fun = "mean"
)

catchments_gpkg$clim_class <- factor(
    catchments_gpkg$clim_class,
    levels = c("Polar", "Cold", "Temperate", "Arid", "Tropical")
)
catchments_gpkg$elev_class <- cut(
    catchments_gpkg$elev_m,
    breaks = c(-Inf, 200, 500, 1000, 2000, Inf),
    labels = c("< 200 m", "200-500 m", "500-1000 m", "1000-2000 m", "> 2000 m")
)
catchments_gpkg$area_class <- cut(
    catchments_gpkg$area_km2,
    breaks = quantile(catchments_gpkg$area_km2,
        probs = c(0, .25, .50, .75, 1),
        na.rm = TRUE
    ),
    labels = c("Q1 (smallest)", "Q2", "Q3", "Q4 (largest)"),
    include.lowest = TRUE
)

# =============================================================================
# SECTION 10: Prepare Spatial Data for Plotting
# =============================================================================
message("Preparing spatial data for maps ...")

# Match catchments to our outlet set
gpkg_in_outlets <- match(cnames, as.character(catchments_gpkg$catch_id))
catchments_plot <- catchments_gpkg[gpkg_in_outlets[!is.na(gpkg_in_outlets)], ]

# Match index: position of each catchments_plot$catch_id in cnames
km_res <- match(as.character(catchments_plot$catch_id), cnames)

# Spatial points for map insets (residual centroids in EPSG:3035)
ppoints_res <- catchments_plot |>
    st_centroid() |>
    st_transform(crs = 3035)

# --- Attach metric values to plot layer ---
# Mean annual precipitation (rain + snow)
rain_dt <- copy(Rain_res)
rain_dt[, year := year(timeStampX)]
yearly_rain <- rain_dt[, lapply(.SD, sum, na.rm = TRUE), by = year, .SDcols = cnames]
mean_annual_rain <- colMeans(yearly_rain[, .SD, .SDcols = cnames], na.rm = TRUE)

snow_dt <- copy(Snow_res)
snow_dt[, year := year(timeStampX)]
yearly_snow <- snow_dt[, lapply(.SD, sum, na.rm = TRUE), by = year, .SDcols = cnames]
mean_annual_snow <- colMeans(yearly_snow[, .SD, .SDcols = cnames], na.rm = TRUE)

catchments_plot$mean_rain <- mean_annual_rain[match(
    as.character(catchments_plot$catch_id),
    names(mean_annual_rain)
)]
catchments_plot$mean_snow <- mean_annual_snow[match(
    as.character(catchments_plot$catch_id),
    names(mean_annual_snow)
)]
catchments_plot$mean_precip <- catchments_plot$mean_rain + catchments_plot$mean_snow

# Snow fraction (mean across years)
SF_means <- colMeans(SaveAgg_SF[, -1, with = FALSE], na.rm = TRUE)
catchments_plot$snowfraction <- SF_means[km_res]

# AET
catchments_plot$AET <- KAET_means[km_res]

# Aridity index (AET/P)
mai_res <- match(as.character(catchments_plot$catch_id), KeepPag$catch)
catchments_plot$AET_P <- KeepPag$AET.P[mai_res]

# Gini
mai_g <- match(as.character(catchments_plot$catch_id), KeepGag$catch)
catchments_plot$gini <- KeepGag$px[mai_g]

# Precipitation seasonality
mai_si <- match(as.character(catchments_plot$catch_id), SeasonPrecipSI_res$cat)
catchments_plot$PSI <- SeasonPrecipSI_res$Seasonality_1[mai_si]
catchments_plot$PSIch <- SeasonPrecipSI_res$SIchange[mai_si]

rm(rain_dt, snow_dt)
gc()

# =============================================================================
# SECTION 11: Map Generation
# =============================================================================
message("Generating maps ...")

# --- Shared theme for all maps ---
map_theme <- theme(
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 11, colour = "grey40"),
    axis.title       = element_text(size = tsize),
    panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
    panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
    panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
    panel.grid.minor = element_line(colour = "grey90"),
    legend.title     = element_text(size = tsize),
    legend.text      = element_text(size = osize),
    legend.position  = "right",
    legend.key       = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size  = unit(0.8, "cm")
)

# --- 11.1 Mean Annual Precipitation ---
palet_precip <- hcl.colors(11, palette = "Roma", rev = FALSE)

ggplot(basemap) +
    geom_sf(fill = "grey95", color = "grey70", linewidth = 0.2) +
    geom_sf(
        data = catchments_plot,
        aes(geometry = geom, fill = mean_precip),
        color = "grey20", stroke = 0.2, alpha = 0.9
    ) +
    scale_fill_gradientn(
        colors = palet_precip, oob = scales::squish,
        name = "Mean annual\nprecipitation (mm/y)",
        breaks = seq(0, 2000, by = 200), limits = c(0, 1500)
    ) +
    coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
    labs(
        title = "Mean annual precipitation - residual catchments",
        subtitle = "Average over 1951-2020, deaggregated to inter-catchment area",
        x = "Longitude", y = "Latitude"
    ) +
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 18, reverse = FALSE)) +
    map_theme

ggsave(paste0(plot_dir, "mean_annual_precipitation_vf.png"),
    width = 22, height = 16, units = "cm", dpi = 300
)

# --- 11.2 Snow Fraction ---
palet_snow <- hcl.colors(9, palette = "YlGnBu", rev = TRUE)

ggplot(basemap) +
    geom_sf(fill = "white", color = NA) +
    geom_sf(
        data = catchments_plot,
        aes(geometry = geom, fill = snowfraction),
        color = "black", stroke = 0
    ) +
    geom_sf(fill = NA, color = "grey20") +
    scale_fill_gradientn(
        colors = palet_snow, oob = scales::squish,
        name = "Snow fraction\n(residual catchment)",
        breaks = seq(0, 1, 0.1), limits = c(0, 0.5)
    ) +
    coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
    labs(x = "Longitude", y = "Latitude") +
    guides(
        fill = guide_colourbar(barwidth = 0.5, barheight = 20, reverse = FALSE),
        size = "none"
    ) +
    map_theme

ggsave(paste0(plot_dir, "snowfraction_map_f.png"),
    width = 20, height = 15, units = "cm", dpi = 500
)

# --- 11.3 AET ---
palet_aet <- hcl.colors(9, palette = "RdYlBu", rev = FALSE)

ggplot(basemap) +
    geom_sf(fill = "white", color = NA) +
    geom_sf(
        data = catchments_plot,
        aes(geometry = geom, fill = AET),
        color = "black", stroke = 0
    ) +
    geom_sf(fill = NA, color = "grey20") +
    scale_fill_gradientn(
        colors = palet_aet, oob = scales::squish,
        name = "AET (mm/y)\n(residual catchment)", limits = c(300, 700)
    ) +
    coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
    labs(x = "Longitude", y = "Latitude") +
    guides(
        fill = guide_colourbar(barwidth = 0.5, barheight = 20, reverse = FALSE),
        size = "none"
    ) +
    map_theme

ggsave(paste0(plot_dir, "AET_nested.png"),
    width = 20, height = 15, units = "cm", dpi = 500
)

# --- 11.4 Aridity Index (AET/P) ---
palet_arid <- hcl.colors(9, palette = "RdYlBu", rev = TRUE)

ggplot(basemap) +
    geom_sf(fill = "white", color = NA) +
    geom_sf(
        data = catchments_plot,
        aes(geometry = geom, fill = AET_P),
        color = "black", stroke = 0
    ) +
    geom_sf(fill = NA, color = "grey20") +
    scale_fill_gradientn(
        colors = palet_arid, oob = scales::squish,
        name = "AET/P (residual)", limits = c(0, 1)
    ) +
    coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
    labs(x = "Longitude", y = "Latitude") +
    guides(
        fill = guide_colourbar(barwidth = 0.5, barheight = 20, reverse = FALSE),
        size = "none"
    ) +
    map_theme

ggsave(paste0(plot_dir, "AridityIndex_nested_f.png"),
    width = 20, height = 15, units = "cm", dpi = 500
)

# --- 11.5 Gini Concentration ---
ggplot(basemap) +
    geom_sf(fill = "white", color = NA) +
    geom_sf(
        data = catchments_plot,
        aes(geometry = geom, fill = gini),
        color = "black", stroke = 0.1
    ) +
    geom_sf(fill = NA, color = "grey20") +
    scale_fill_gradientn(
        colors = palet_arid, oob = scales::squish,
        name = "Precip. GINI\n(residual catchment)", limits = c(0.5, 1)
    ) +
    coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
    labs(x = "Longitude", y = "Latitude") +
    guides(
        fill = guide_colourbar(barwidth = 0.5, barheight = 20, reverse = FALSE),
        size = "none"
    ) +
    map_theme

ggsave(paste0(plot_dir, "PGini_nested.png"),
    width = 20, height = 15, units = "cm", dpi = 500
)

# --- 11.6 Precipitation Seasonality ---
ggplot(basemap) +
    geom_sf(fill = "white", color = NA) +
    geom_sf(
        data = catchments_plot,
        aes(geometry = geom, fill = PSI),
        color = "black", stroke = 0.1
    ) +
    geom_sf(fill = NA, color = "grey20") +
    scale_fill_gradientn(
        colors = palet_arid, oob = scales::squish,
        name = "Precip. seasonality\n(residual catchment)", limits = c(0, 0.5)
    ) +
    coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
    labs(x = "Longitude", y = "Latitude") +
    guides(
        fill = guide_colourbar(barwidth = 0.5, barheight = 20, reverse = FALSE),
        size = "none"
    ) +
    map_theme

ggsave(paste0(plot_dir, "PSI_nested.png"),
    width = 20, height = 15, units = "cm", dpi = 500
)

# --- 11.7 Precipitation Seasonality Change ---
ggplot(basemap) +
    geom_sf(fill = "white", color = NA) +
    geom_sf(
        data = catchments_plot,
        aes(geometry = geom, fill = PSIch),
        color = "black", stroke = 0.1
    ) +
    geom_sf(fill = NA, color = "grey20") +
    scale_fill_gradientn(
        colors = palet_arid, oob = scales::squish,
        name = "Precip. seasonality change\n(residual catchment)",
        limits = c(-0.1, 0.1)
    ) +
    coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
    labs(x = "Longitude", y = "Latitude") +
    guides(
        fill = guide_colourbar(barwidth = 0.5, barheight = 20, reverse = FALSE),
        size = "none"
    ) +
    map_theme

ggsave(paste0(plot_dir, "PSIchange_nested.png"),
    width = 20, height = 15, units = "cm", dpi = 500
)

# =============================================================================
# SECTION 12: Batch Deaggregation of All Variables
# =============================================================================
# Deaggregates each variable and saves as CSV with "_nested" suffix.
# Discharge (disWin) uses volume-based mass balance (Q = TRUE).
message("Starting batch deaggregation ...")

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

    Vari_res <- deaggregate_to_residual(Vari, imm_children_list,
        area_lookup, resarea_lookup,
        Q = is_discharge
    )

    # Clip negative values for non-discharge variables
    if (!is_discharge) {
        Vari_res <- Vari_res[, lapply(.SD, function(x) pmax(0, x))]
    }

    # Prepend sorted timestamp column
    Vari_res <- data.table(time = time_v[order(time_v)], Vari_res)

    out_file <- paste0(base_dir, "data/", var, "_nested_1951_2020.csv")
    fwrite(Vari_res, out_file,
        sep = ",", na = "-9999", row.names = FALSE,
        quote = FALSE, showProgress = TRUE, nThread = 2
    )

    rm(Vari, Vari_res)
    gc()
}

message("deaggregate_and_metrics.R completed successfully.")
