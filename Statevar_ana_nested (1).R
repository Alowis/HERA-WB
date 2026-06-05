# =============================================================================
# Statevar_ana_nested.R
# -----------------------------------------------------------------------------
# Re-implementation of Statevar_ana.R where every climatic / hydrologic
# variable is first deaggregated to the *residual* (inter-catchment) area
# before computing derived metrics (snow fraction, aridity, Gini, seasonality).
#
# CONCEPT
# -------
# The original CSV files store spatially-aggregated mean values over the full
# upstream area of each outlet (i.e. including all nested sub-catchments).
# For a nested outlet i with immediate children j1, j2, ...:
#
#   residual_value_i = (Value_i * Area_i  -  Sum(Value_j * Area_j)) / ResArea_i
#
# where ResArea_i is already stored in the .gpkg as residual_area_km2.
#
# Headwater catchments (immediate_nested_ids == NA) keep their original values.
#
# SELF-CONTAINED: all nesting relationships are read from the .gpkg file
# produced by Nested_Catchments.R — no in-memory matrix objects needed.
# =============================================================================

# --------------------------------------------------------------------------- #
# 0.  Libraries                                                               #
# --------------------------------------------------------------------------- #
library(sf)
library(dplyr)
library(data.table)
library(lubridate)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)

source("functions_regime.R")

# --------------------------------------------------------------------------- #
# 1.  Paths                                                                   #
# --------------------------------------------------------------------------- #
hydroDir  <- "D:/tilloal/Documents/LFRuns_utils/data/"
regimeDir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
gpkg_path <- paste0(regimeDir, "data/catchments_analysis_results.gpkg")

# --------------------------------------------------------------------------- #
# 2.  Load nesting structure from the .gpkg                                   #
# --------------------------------------------------------------------------- #
message("Reading catchments_analysis_results.gpkg ...")
catchments_gpkg <- st_read(gpkg_path, quiet = TRUE)

# Expected columns (written by Nested_Catchments.R):
#   catch_id, area_km2, residual_area_km2,
#   immediate_nested_ids   <- comma-separated catch_ids of immediate children
#   all_nested_ids         <- comma-separated catch_ids of all descendants

stopifnot(all(c("catch_id", "area_km2", "residual_area_km2",
                "immediate_nested_ids") %in% names(catchments_gpkg)))

# Build a lookup: catch_id -> list of immediate child catch_ids (as character)
# immediate_nested_ids is NA for headwaters, else "id1, id2, ..." string
parse_ids <- function(x) {
  if (is.na(x) || trimws(x) == "") return(character(0))
  trimws(strsplit(as.character(x), ",")[[1]])
}

# Named list: for each catch_id, the vector of immediate child catch_ids
imm_children_list <- setNames(
  lapply(catchments_gpkg$immediate_nested_ids, parse_ids),
  as.character(catchments_gpkg$catch_id)
)

# area_km2 lookup by catch_id (character key)
area_lookup <- setNames(catchments_gpkg$area_km2,
                        as.character(catchments_gpkg$catch_id))

# residual_area_km2 lookup
resarea_lookup <- setNames(catchments_gpkg$residual_area_km2,
                           as.character(catchments_gpkg$catch_id))

message(sprintf("Loaded %d catchments (%d headwaters, %d nested).",
                nrow(catchments_gpkg),
                sum(sapply(imm_children_list, length) == 0),
                sum(sapply(imm_children_list, length) >  0)))

# --------------------------------------------------------------------------- #
# 3.  Base spatial setup  (identical to Statevar_ana.R)                      #
# --------------------------------------------------------------------------- #
world  <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"), ]

outletname   <- "outletsv8_hybas07_01min"
outhybas_raw <- outletopen(hydroDir, outletname)
outhybas_raw$latlong <- paste(round(outhybas_raw$Var1, 4),
                              round(outhybas_raw$Var2, 4), sep = " ")

# Plot bounding-box coordinates in EPSG:3035
outll    <- outletopen(hydroDir, outletname)
cord.dec <- SpatialPoints(outll[, c(2, 3)],
                          proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco      <- cord.UTM@coords

w2      <- st_transform(world, crs = 3035)
basemap <- w2
tsize   <- 12
osize   <- 12

# Upstream area file
outhybas_raw$idlalo  <- paste(outhybas_raw$idlo, outhybas_raw$idla, sep = " ")
outhybas_raw$latlong <- paste(round(outhybas_raw$Var1, 4),
                               round(outhybas_raw$Var2, 4), sep = " ")
UpArea_full <- UpAopen(hydroDir, "upArea_European_01min.nc", outhybas_raw)

out1         <- outletopen(hydroDir, "efas_rnet_100km_01min")
out1$latlong <- paste(round(out1$Var1, 4), round(out1$Var2, 4), sep = " ")
outhybas_eu  <- inner_join(out1, outhybas_raw, by = "latlong")

# Upstream-catchment filter
UpCat    <- st_read(dsn = paste0(regimeDir, "Upstreamgroups.shp"), quiet = TRUE)
upup     <- UpCat[which(UpCat$up == "Upstream Catchments"), ]
hybas2uc <- match(upup$outlets_x, outhybas_eu$outlets.y)
outhybas <- outhybas_eu[hybas2uc, ]

matcat <- match(outhybas$latlong, UpArea_full$latlong)
UpArea <- UpArea_full[matcat, ]

# Column selector: positions in CSV files matching our outlet IDs
# (read one header row to get the column names)
header_ref <- fread(paste0(hydroDir, "tss/HERA_Histo/SnowUpsX_1951_2020.csv"),
                    nrows = 0, header = TRUE)
matcol <- match(UpArea$outlets, as.numeric(colnames(header_ref)))
matcol <- matcol[!is.na(matcol)]
cnames <- as.character(UpArea$outlets)   # character IDs of our outlets

# --------------------------------------------------------------------------- #
# 4.  Core deaggregation function                                             #
# --------------------------------------------------------------------------- #
# Operates entirely from the .gpkg-derived lookup tables — no matrix needed.
#
# INPUT
#   dt_raw          : data.table [nTime x nOutlets], columns named by catch_id
#   imm_children    : named list(catch_id -> character vector of child ids)
#   area_lkp        : named numeric vector of area_km2  (key = catch_id)
#   resarea_lkp     : named numeric vector of residual_area_km2
#
# OUTPUT
#   data.table same shape as dt_raw, nested columns replaced by residual value
# --------------------------------------------------------------------------- #
deaggregate_to_residual <- function(dt_raw,
                                    imm_children,
                                    area_lkp,
                                    resarea_lkp) {
  res_dt    <- copy(dt_raw)
  col_names <- names(dt_raw)            # outlet IDs present in this CSV

  # Only process outlets that actually appear in the data
  outlets_in_data <- intersect(names(imm_children), col_names)

  for (p_id in outlets_in_data) {

    child_ids <- imm_children[[p_id]]

    # --- HEADWATER: nothing to subtract ---
    if (length(child_ids) == 0) next

    p_area   <- area_lkp[p_id]
    res_area <- resarea_lkp[p_id]

    if (is.na(res_area) || res_area <= 0) next

    # Keep only children that have a column in this data.table
    valid_children <- intersect(child_ids, col_names)
    if (length(valid_children) == 0) next

    child_areas <- area_lkp[valid_children]

    # --- MASS BALANCE ---
    # child volume matrix: each column = child_series * child_area
    child_vols <- dt_raw[, mapply(`*`, .SD, child_areas),
                         .SDcols = valid_children]

    total_child_vol <- if (length(valid_children) > 1) {
      rowSums(child_vols, na.rm = TRUE)
    } else {
      as.vector(child_vols)
    }

    res_series <- (dt_raw[[p_id]] * p_area - total_child_vol) / res_area

    set(res_dt, j = p_id, value = res_series)
  }

  return(res_dt)
}

# --------------------------------------------------------------------------- #
# 5.  Load raw Rain & Snow, deaggregate, build Precipitation & SnowFraction  #
# --------------------------------------------------------------------------- #
message("Loading Rain ...")
Rain_raw <- fread(paste0(hydroDir, "tss/HERA_SocCF/RainUpsX_1951_2020.csv"),
                  header = TRUE)
time       <- Rain_raw$V1
timeStampX <- time[order(time)]
Rain_raw   <- Rain_raw[order(time), .SD, .SDcols = matcol]

message("Loading Snow ...")
Snow_raw <- fread(paste0(hydroDir, "tss/HERA_Histo/SnowUpsX_1951_2020.csv"),
                  header = TRUE)
Snow_raw <- Snow_raw[order(time), .SD, .SDcols = matcol]

message("Deaggregating Rain ...")
Rain_res <- deaggregate_to_residual(Rain_raw, imm_children_list,
                                    area_lookup, resarea_lookup)
Rain_res <- Rain_res[, lapply(.SD, function(x) pmax(0, x))]

message("Deaggregating Snow ...")
Snow_res <- deaggregate_to_residual(Snow_raw, imm_children_list,
                                    area_lookup, resarea_lookup)
Snow_res <- Snow_res[, lapply(.SD, function(x) pmax(0, x))]

Precipitation_res <- Rain_res + Snow_res
Snowfraction_res  <- Snow_res / Precipitation_res   # NaN where P = 0 is OK

rm(Rain_raw, Snow_raw, Rain_res, Snow_res); gc()

# --------------------------------------------------------------------------- #
# 6.  Snow Fraction metric (yearly)                                           #
# --------------------------------------------------------------------------- #
message("Computing Snow Fraction ...")
dsel    <- hour(timeStampX)
dt_step <- 4                                        # 6-hourly steps

SaveAgg_SF <- c(1951:2020)
for (col_name in cnames) {
  print(col_name)
  Trun <- preprocess_frac(time = timeStampX, input_var = Snowfraction_res,
                          col_name, dsel, dt_step)
  Yagg <- process_frac(Trun)
  SaveAgg_SF <- cbind(SaveAgg_SF, Yagg$val)
}
SaveAgg_SF           <- data.frame(SaveAgg_SF)
colnames(SaveAgg_SF)[-1] <- cnames
save(SaveAgg_SF, file = paste0(regimeDir, "SnowFraction_nested.Rdata"))

# --------------------------------------------------------------------------- #
# 7.  Spatial points for plotting (centroids of residual geometries)          #
# --------------------------------------------------------------------------- #
# Match gpkg rows to our outlet list so point order is consistent
gpkg_in_outlets <- match(cnames, as.character(catchments_gpkg$catch_id))
catchments_plot  <- catchments_gpkg[gpkg_in_outlets[!is.na(gpkg_in_outlets)], ]

res_cents   <- st_centroid(st_geometry(catchments_plot))
ppoints_res <- st_as_sf(
  data.frame(
    catch_id = catchments_plot$catch_id,
    upa      = catchments_plot$area_km2,          # full upstream (for sizing)
    res_area = catchments_plot$residual_area_km2
  ),
  geometry = res_cents,
  crs      = st_crs(catchments_plot)
) |> st_transform(crs = 3035)

# Match vector: position of each ppoints_res$catch_id in cnames
km_res <- match(as.character(ppoints_res$catch_id), cnames)

# --------------------------------------------------------------------------- #
# 8.  Plot: Snow Fraction (residual)                                          #
# --------------------------------------------------------------------------- #
ys1 <- 1951:1970
my  <- which(!is.na(match(SaveAgg_SF$SaveAgg_SF, ys1)))
SF_means              <- colMeans(SaveAgg_SF[my, -1], na.rm = TRUE)
ppoints_res$snowfraction <- SF_means[km_res]

palet <- hcl.colors(9, palette = "YlGnBu", rev = TRUE)

ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = ppoints_res,
          aes(geometry = geometry, size = upa, fill = snowfraction),
          color = "black", shape = 21, stroke = 0) +
  geom_sf(fill = NA, color = "grey20") +
  scale_x_continuous(breaks = seq(-30, 40, by = 5)) +
  scale_size(range = c(1, 4), trans = "sqrt",
             name = expression(paste("Upstream area ", (km^2))),
             breaks = c(101, 1000, 10000, 100000, 500000),
             labels = c("100", "1000", "10 000", "100 000", "500 000")) +
  scale_fill_gradientn(colors = palet, oob = scales::squish,
                       name = "Snow fraction\n(residual catchment)",
                       breaks = seq(0, 1, 0.2), limits = c(0, 0.5)) +
  coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20, reverse = FALSE),
         size = "none") +
  theme(axis.title       = element_text(size = tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
        legend.title     = element_text(size = tsize),
        legend.text      = element_text(size = osize),
        legend.position  = "right",
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key       = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size  = unit(0.8, "cm"))

ggsave(paste0(regimeDir, "plots/snowfraction_nested.png"),
       width = 20, height = 15, units = "cm", dpi = 1500)

# --------------------------------------------------------------------------- #
# 9.  AET and PET: load, deaggregate, compute yearly metrics                 #
# --------------------------------------------------------------------------- #
message("Loading and deaggregating AET and PET ...")

ActEvapo_raw <- fread(paste0(hydroDir, "tss/HERA_Histo/ActEvapo_1951_2020.csv"),
                      header = TRUE)
ActEvapo_raw <- ActEvapo_raw[order(time), .SD, .SDcols = matcol]

PEvapo_raw   <- fread(paste0(hydroDir, "tss/HERA_Histo/etUpsX_1951_2020.csv"),
                      header = TRUE)
PEvapo_raw   <- PEvapo_raw[order(time), .SD, .SDcols = matcol]

ActEvapo_res <- deaggregate_to_residual(ActEvapo_raw, imm_children_list,
                                        area_lookup, resarea_lookup)
ActEvapo_res <- ActEvapo_res[, lapply(.SD, function(x) pmax(0, x))]

PEvapo_res   <- deaggregate_to_residual(PEvapo_raw, imm_children_list,
                                        area_lookup, resarea_lookup)
PEvapo_res   <- PEvapo_res[, lapply(.SD, function(x) pmax(0, x))]

rm(ActEvapo_raw, PEvapo_raw); gc()

# --------------------------------------------------------------------------- #
# 10. Precipitation, AET and E0 yearly aggregates                             #
# --------------------------------------------------------------------------- #
message("Computing yearly aggregates ...")

SaveE0_res      <- c(1951:2020)
SaveAET_res     <- c(1951:2020)
SavePrecip_res  <- NULL
SeasonPrecip_res <- NULL

for (col_name in cnames) {
  print(col_name)
  AET    <- preprocess_in(time = timeStampX, input_var = ActEvapo_res,
                          col_name, dsel, dt_step)
  Precip <- preprocess_in(time = timeStampX, input_var = Precipitation_res,
                          col_name, dsel, dt_step)
  TE0    <- preprocess_in(time = timeStampX, input_var = PEvapo_res,
                          col_name, dsel, dt_step)

  E0Acc     <- process_data(TE0)
  AETAcc    <- process_data(AET)
  PrecipAcc <- process_precip(Precip)

  PrecipAcc$seaonal$cat <- rep(col_name, 2)
  PrecipAcc$yagg$cat    <- rep(col_name, length(PrecipAcc$yagg$year))

  SavePrecip_res   <- rbind(SavePrecip_res,   PrecipAcc$yagg)
  SeasonPrecip_res <- rbind(SeasonPrecip_res, PrecipAcc$seaonal)
  SaveAET_res      <- cbind(SaveAET_res,      AETAcc$val)
  SaveE0_res       <- cbind(SaveE0_res,       E0Acc$val)
}

save(SavePrecip_res,   file = paste0(regimeDir, "PrecipitationMetrics_nested.Rdata"))
save(SeasonPrecip_res, file = paste0(regimeDir, "PrecipitationSeason_nested.Rdata"))
save(SaveAET_res,      file = paste0(regimeDir, "AETMetrics_nested.Rdata"))

rm(ActEvapo_res, PEvapo_res, Precipitation_res, Snowfraction_res); gc()

# --------------------------------------------------------------------------- #
# 11. Plot: AET (residual)                                                    #
# --------------------------------------------------------------------------- #
SaveAET_res_df <- data.frame(SaveAET_res)
my             <- which(!is.na(match(SaveAET_res_df$SaveAET_res, ys1)))
KAET_means     <- colMeans(SaveAET_res_df[my, -1], na.rm = TRUE)
ppoints_res$AET1 <- KAET_means[km_res]

palet <- hcl.colors(9, palette = "RdYlBu", rev = FALSE)

ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = ppoints_res,
          aes(geometry = geometry, size = upa, fill = AET1),
          color = "black", shape = 21, stroke = 0) +
  geom_sf(fill = NA, color = "grey20") +
  scale_x_continuous(breaks = seq(-30, 40, by = 5)) +
  scale_size(range = c(1, 4), trans = "sqrt",
             name = expression(paste("Upstream area ", (km^2))),
             breaks = c(101, 1000, 10000, 100000, 500000),
             labels = c("100", "1000", "10 000", "100 000", "500 000")) +
  scale_fill_gradientn(colors = palet, oob = scales::squish,
                       name = "AET (mm/y)\n(residual catchment)") +
  coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20, reverse = FALSE),
         size = "none") +
  theme(axis.title       = element_text(size = tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
        legend.title     = element_text(size = tsize),
        legend.text      = element_text(size = osize),
        legend.position  = "right",
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key       = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size  = unit(0.8, "cm"))

ggsave(paste0(regimeDir, "plots/AET_nested.png"),
       width = 20, height = 15, units = "cm", dpi = 1500)

# --------------------------------------------------------------------------- #
# 12. AET/P aridity index (residual)                                          #
# --------------------------------------------------------------------------- #
myp     <- which(!is.na(match(SavePrecip_res$year, ys1)))
KeepP   <- SavePrecip_res[myp, c(1, 2, 4)]

KeepPag <- stats::aggregate(list(pmean = KeepP$val.sum),
                            by  = list(catch = KeepP$cat),
                            FUN = function(x) mean(x, na.rm = TRUE))
KeepPag <- do.call(data.frame, KeepPag)

kpe           <- match(KeepPag$catch, cnames)
KeepPag$aet   <- KAET_means[kpe]
KeepPag$AET.P <- KeepPag$aet / KeepPag$pmean

mai_res            <- match(as.character(ppoints_res$catch_id), KeepPag$catch)
ppoints_res$AET.P1 <- KeepPag$AET.P[mai_res]

palet <- hcl.colors(9, palette = "RdYlBu", rev = TRUE)

ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = ppoints_res,
          aes(geometry = geometry, size = upa, fill = AET.P1),
          color = "black", shape = 21, stroke = 0) +
  geom_sf(fill = NA, color = "grey20") +
  scale_x_continuous(breaks = seq(-30, 40, by = 5)) +
  scale_size(range = c(1, 4), trans = "sqrt",
             name = expression(paste("Upstream area ", (km^2))),
             breaks = c(101, 1000, 10000, 100000, 500000),
             labels = c("100", "1000", "10 000", "100 000", "500 000")) +
  scale_fill_gradientn(colors = palet, oob = scales::squish,
                       name = "AET/P (residual)", limits = c(0, 1)) +
  coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20, reverse = FALSE),
         size = "none") +
  theme(axis.title       = element_text(size = tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
        legend.title     = element_text(size = tsize),
        legend.text      = element_text(size = osize),
        legend.position  = "right",
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key       = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size  = unit(0.8, "cm"))

ggsave(paste0(regimeDir, "plots/AridityIndex_nested.png"),
       width = 20, height = 15, units = "cm", dpi = 1500)

# --------------------------------------------------------------------------- #
# 13. Rainfall Gini concentration index (residual)                            #
# --------------------------------------------------------------------------- #
KeepG <- SavePrecip_res[myp, c(1, 3, 4)]

KeepGag <- stats::aggregate(list(px = KeepG$val.gini),
                            by  = list(catch = KeepG$cat),
                            FUN = function(x) mean(x, na.rm = TRUE))
KeepGag <- do.call(data.frame, KeepGag)

mai_g             <- match(as.character(ppoints_res$catch_id), KeepGag$catch)
ppoints_res$gini1 <- KeepGag$px[mai_g]

ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = ppoints_res,
          aes(geometry = geometry, size = upa, fill = gini1),
          color = "black", shape = 21, stroke = 0.1) +
  geom_sf(fill = NA, color = "grey20") +
  scale_x_continuous(breaks = seq(-30, 40, by = 5)) +
  scale_size(range = c(1, 4), trans = "sqrt",
             name = expression(paste("Upstream area ", (km^2))),
             breaks = c(101, 1000, 10000, 100000, 500000),
             labels = c("100", "1000", "10 000", "100 000", "500 000")) +
  scale_fill_gradientn(colors = palet, oob = scales::squish,
                       name = "Precip. concentration\n(residual catchment)",
                       limits = c(0.5, 1)) +
  coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20, reverse = FALSE),
         size = "none") +
  theme(axis.title       = element_text(size = tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
        legend.title     = element_text(size = tsize),
        legend.text      = element_text(size = osize),
        legend.position  = "right",
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key       = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size  = unit(0.8, "cm"))

ggsave(paste0(regimeDir, "plots/PGini_nested.png"),
       width = 20, height = 15, units = "cm", dpi = 500)

# --------------------------------------------------------------------------- #
# 14. Precipitation Seasonality Index (residual)                              #
# --------------------------------------------------------------------------- #
ls_sp                <- length(SeasonPrecip_res$cat)
SeasonPrecipSI_res   <- SeasonPrecip_res[-seq(1, ls_sp, 2), ]
SeasonPrecipSI_res$SIchange <- SeasonPrecipSI_res$Seasonality_2 -
                               SeasonPrecipSI_res$Seasonality_1

mai_si             <- match(as.character(ppoints_res$catch_id), SeasonPrecipSI_res$cat)
ppoints_res$PSI    <- SeasonPrecipSI_res$Seasonality_1[mai_si]
ppoints_res$PSIch  <- SeasonPrecipSI_res$SIchange[mai_si]

# --- Seasonality level ---
ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = ppoints_res,
          aes(geometry = geometry, size = upa, fill = PSI),
          color = "black", shape = 21, stroke = 0.1) +
  geom_sf(fill = NA, color = "grey20") +
  scale_x_continuous(breaks = seq(-30, 40, by = 5)) +
  scale_size(range = c(1, 4), trans = "sqrt",
             name = expression(paste("Upstream area ", (km^2))),
             breaks = c(101, 1000, 10000, 100000, 500000),
             labels = c("100", "1000", "10 000", "100 000", "500 000")) +
  scale_fill_gradientn(colors = palet, oob = scales::squish,
                       name = "Precip. seasonality\n(residual catchment)",
                       limits = c(0, 0.5)) +
  coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20, reverse = FALSE),
         size = "none") +
  theme(axis.title       = element_text(size = tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
        legend.title     = element_text(size = tsize),
        legend.text      = element_text(size = osize),
        legend.position  = "right",
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key       = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size  = unit(0.8, "cm"))

ggsave(paste0(regimeDir, "plots/PSI_nested.png"),
       width = 20, height = 15, units = "cm", dpi = 500)

# --- Seasonality change ---
ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = ppoints_res,
          aes(geometry = geometry, size = upa, fill = PSIch),
          color = "black", shape = 21, stroke = 0.1) +
  geom_sf(fill = NA, color = "grey20") +
  scale_x_continuous(breaks = seq(-30, 40, by = 5)) +
  scale_size(range = c(1, 4), trans = "sqrt",
             name = expression(paste("Upstream area ", (km^2))),
             breaks = c(101, 1000, 10000, 100000, 500000),
             labels = c("100", "1000", "10 000", "100 000", "500 000")) +
  scale_fill_gradientn(colors = palet, oob = scales::squish,
                       name = "Precip. seasonality change\n(residual catchment)",
                       limits = c(-0.1, 0.1)) +
  coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
  labs(x = "Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20, reverse = FALSE),
         size = "none") +
  theme(axis.title       = element_text(size = tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
        legend.title     = element_text(size = tsize),
        legend.text      = element_text(size = osize),
        legend.position  = "right",
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key       = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size  = unit(0.8, "cm"))

ggsave(paste0(regimeDir, "plots/PSIchange_nested.png"),
       width = 20, height = 15, units = "cm", dpi = 500)

# --------------------------------------------------------------------------- #
# 15. Batch deaggregation of other model output variables                     #
# --------------------------------------------------------------------------- #
# Deaggregates and saves every variable listed below.
# Output files go to regimeDir/data/ with a "_nested" suffix.

outputfilenames <- c(
  "rainUpsX", "snowMeltUpsX", "qUzUpsX", "qLZUpsX",
  "surfaceRunoffUpsX", "infUpsX", "ActEvapo", "disWin",
  "percUZLZUpsX", "dSubToUzUpsX", "prefFlowUpsX", "lossUpsX"
)

for (var in outputfilenames) {
  message("Batch deaggregating: ", var)

  in_file <- paste0(hydroDir, "tss/HERA_Histo/", var, "_1951_2020.csv")
  if (!file.exists(in_file)) {
    message("  File not found, skipping: ", in_file)
    next
  }

  Vari      <- fread(in_file, header = TRUE)
  Vari      <- Vari[order(Vari$V1), ]
  matcol_v  <- match(UpArea$outlets, as.numeric(colnames(Vari)))
  matcol_v  <- matcol_v[!is.na(matcol_v)]
  Vari      <- Vari[, .SD, .SDcols = matcol_v]

  Vari_res  <- deaggregate_to_residual(Vari, imm_children_list,
                                       area_lookup, resarea_lookup)
  Vari_res  <- Vari_res[, lapply(.SD, function(x) pmax(0, x))]

  fwrite(Vari_res,
         paste0(regimeDir, "data/", var, "_nested_1951_2020.csv"),
         sep = ",", na = "-9999", row.names = FALSE,
         quote = FALSE, showProgress = TRUE, nThread = 2)

  rm(Vari, Vari_res); gc()
}

message("Statevar_ana_nested.R completed successfully.")
