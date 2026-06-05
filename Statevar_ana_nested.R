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
gpkg_path <- paste0(regimeDir, "data/catchments_analysis_final_v3.gpkg")


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

# --------------------------------------------------------------------------- #
# 2.  Load nesting structure from the .gpkg                                   #
# --------------------------------------------------------------------------- #
message("Reading catchments_analysis_final.gpkg ...")
catchments_gpkg <- st_read(gpkg_path, quiet = TRUE)

# --------------------------------------------------------------------------- #
# 2b.  Diagnostic plot: headwaters vs nested catchments                       #
# --------------------------------------------------------------------------- #

# Count nested catchments per outlet (0 = headwater)
catchments_gpkg$n_nested <- sapply(
  catchments_gpkg$immediate_nested_ids,
  function(x) {
    if (is.na(x) || trimws(x) == "NA" || trimws(x) == "") return(0L)
    length(trimws(strsplit(x, ",")[[1]]))
  }
)

catchments_gpkg$is_headwater <- catchments_gpkg$n_nested == 0

# Centroids in EPSG:3035 for plotting
plot_pts <- catchments_gpkg |>
  st_centroid() |>
  st_transform(crs = 3035)

# World basemap
world   <- ne_countries(scale = "medium", returnclass = "sf")
basemap <- st_transform(world, crs = 3035)

# Bounding box from the catchment centroids themselves
coords  <- st_coordinates(plot_pts)

# Split into headwaters and nested for layering
hw_pts     <- plot_pts[plot_pts$is_headwater, ]
nested_pts <- plot_pts[!plot_pts$is_headwater, ]

hw_poly     <- catchments_gpkg[catchments_gpkg$is_headwater, ]
nested_poly <- catchments_gpkg[!catchments_gpkg$is_headwater, ]

ggplot(basemap) +
  geom_sf(fill = "grey95", color = "grey70", linewidth = 0.2) +
  
  # Headwaters: fixed grey fill, size still encodes upstream area
  geom_sf(data  = hw_pts,
          aes(size = area_km2),
          fill  = "grey60", color = "grey30",
          shape = 21, stroke = 0.2, alpha = 0.7) +
  
  # Nested outlets: color encodes number of immediate children
  geom_sf(data = nested_pts,
          aes(size = area_km2, fill = nesting_level),
          color = "grey20", shape = 21, stroke = 0.2, alpha = 0.9) +
  
  scale_fill_gradientn(
    colors = hcl.colors(9, palette = "Zissou 1", rev = FALSE),
    name   = "No. immediate\nnested catchments",
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_size(
    range  = c(0.8, 5),
    trans  = "sqrt",
    name   = expression(paste("Upstream area (", km^2, ")")),
    breaks = c(100, 1000, 10000, 100000, 500000),
    labels = c("100", "1 000", "10 000", "100 000", "500 000")
  ) +
  coord_sf(
    xlim = c(min(coords[, 1]) - 1e5, max(coords[, 1]) + 1e5),
    ylim = c(min(coords[, 2]) - 1e5, max(coords[, 2]) + 1e5)
  ) +
  labs(
    title    = "Catchment nesting structure",
    subtitle = sprintf("%d headwaters (grey)  |  %d nested outlets (coloured)",
                       sum(catchments_gpkg$is_headwater),
                       sum(!catchments_gpkg$is_headwater)),
    x = "Longitude", y = "Latitude"
  ) +
  guides(
    fill = guide_colourbar(barwidth = 0.5, barheight = 15, reverse = FALSE),
    size = guide_legend(override.aes = list(fill = "grey60", color = "grey30"))
  ) +
  theme(
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 11, colour = "grey40"),
    axis.title       = element_text(size = 12),
    panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
    panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
    panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
    panel.grid.minor = element_line(colour = "grey90"),
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    legend.position  = "right",
    legend.key       = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size  = unit(0.8, "cm")
  )



ggsave(paste0(regimeDir, "plots/nesting_structure.png"),
       width = 22, height = 16, units = "cm", dpi = 300)

#histogram of upstream area (residual)
library(scales)
pu<-ggplot(catchments_gpkg, aes(x=residual_area_km2)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=20,alpha=0.9,lwd=1)+
  scale_y_continuous(breaks=seq(0,600, by=100),name="Number of outlets")+
  scale_x_log10(name=expression(paste("Upstream area ", (km^2),sep = " ")),
                breaks=c(1,10,100,1000,10000,100000), minor_breaks = log10_minor_break(),
                labels=c("1","10","100","1 000","10 000","100 000")) +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=16),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  annotate("label", x=20000, y=400, label= paste0("n = ",length(CatUpA$upa)),size=6)

pu
ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/histo_outlets_residuals.jpg", pu, width=20, height=15, units=c("cm"),dpi=1500)


library(cowplot)

# --- Main map ---
p_map <- ggplot(basemap) +
  geom_sf(fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data  = hw_poly,
          fill  = "lightblue", color = "grey30",
          shape = 21, stroke = 0.2, alpha = 0.7) +
  geom_sf(data = nested_poly,
          aes(fill = nesting_level),
          color = "grey20", shape = 21, stroke = 0.2, alpha = 0.9) +
  scale_fill_gradientn(
    colors = hcl.colors(9, palette = "BluYl", rev = FALSE),
    name   = "No. \nnested catchments", limits=c(0,50),oob = scales::squish,
    breaks = scales::pretty_breaks(n = 5)
  ) +
  # scale_size(
  #   range  = c(0.8, 5),
  #   trans  = "sqrt",
  #   name   = expression(paste("Upstream area (", km^2, ")")),
  #   breaks = c(100, 1000, 10000, 100000, 500000),
  #   labels = c("100", "1 000", "10 000", "100 000", "500 000")
  # ) +
  coord_sf(
    xlim = c(min(coords[, 1]) - 1e5, max(coords[, 1]) + 1e5),
    ylim = c(min(coords[, 2]) - 1e5, max(coords[, 2]) + 1e5)
  ) +
  labs(
    title    = "Catchments",
    subtitle = sprintf("%d headwaters (lightblue)  |  %d nested outlets (coloured)",
                       sum(catchments_gpkg$is_headwater),
                       sum(!catchments_gpkg$is_headwater)),
    x = "Longitude", y = "Latitude"
  ) +
  guides(
    fill = guide_colourbar(barwidth = 0.5, barheight = 15, reverse = FALSE),
    size = guide_legend(override.aes = list(fill = "grey60", color = "grey30"))
  ) +
  theme(
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 11, colour = "grey40"),
    axis.title       = element_text(size = 12),
    panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
    panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
    panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
    panel.grid.minor = element_line(colour = "grey90"),
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    legend.position  = "right",
    legend.key       = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size  = unit(0.8, "cm")
  )

# --- Histogram inset ---
p_hist <- ggplot(catchments_gpkg, aes(x = residual_area_km2)) +
  geom_histogram(color = "steelblue", fill = "slategray1",
                 bins = 20, alpha = 0.9, lwd = .2) +
  scale_y_continuous(breaks = seq(0, 600, by = 100),
                     name = "Number of outlets") +
  scale_x_log10(
    name   = expression(paste("Residual area (", km^2, ")")),
    breaks = c(1, 10, 100, 1000, 10000, 100000),
    minor_breaks = log10_minor_break(),
    labels = c("1", "10", "100", "1 000", "10 000", "100 000")
  ) +
  # annotate("label", x = 20000, y = 400,
  #          label = paste0("n = ", nrow(catchments_gpkg)), size = 1.2) +
  theme(
    axis.title       = element_text(size = 4),
    axis.ticks.length = unit(.5, "pt"),
    axis.text        = element_text(size = 3),
    axis.ticks       = element_line(linewidth =.1),
    axis.title.x      = element_text(margin = margin(t = 0.2)),
    axis.title.y      = element_text(margin = margin(r = 0.2)),
    axis.text.x      = element_text(margin = margin(t = 0.2)),
    axis.text.y      = element_text(margin = margin(r = 0.2)),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black",linewidth =.1),
    panel.grid.major = element_line(colour = "grey80",linewidth =.1),
    panel.grid.minor.x = element_line(colour = "grey90", linetype = "dashed",linewidth =.1),
    plot.background  = element_rect(fill = "white", colour = "grey50",
                                    linewidth = 0.4),  # border around inset
    plot.margin      = margin(2, 3, 1, 1) 
  )

# --- Combine: histogram inset in top-right corner of the map ---
combined <- ggdraw(p_map) +
  draw_plot(p_hist,
            x      = 0.55,   # left edge of inset (fraction of full plot width)
            y      = 0.72,   # bottom edge of inset
            width  = 0.18,   # inset width  (adjust to taste)
            height = 0.18)   # inset height (adjust to taste)

# combined
ggsave(paste0(regimeDir, "plots/Figure1_catchments_v2.png"),
       combined,
       width = 22, height = 16, units = "cm", dpi = 500)


################################################################################

mean(catchments_gpkg$residual_area_km2)
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
                sum(catchments_gpkg$immediate_nested_ids == "NA", na.rm = TRUE),
                sum(catchments_gpkg$immediate_nested_ids != "NA", na.rm = TRUE)))





# Upstream area file
outhybas_raw$idlalo  <- paste(outhybas_raw$idlo, outhybas_raw$idla, sep = " ")
outhybas_raw$latlong <- paste(round(outhybas_raw$Var1, 4),
                               round(outhybas_raw$Var2, 4), sep = " ")
UpArea_full <- UpAopen(hydroDir, "upArea_European_01min.nc", outhybas_raw)

out1         <- outletopen(hydroDir, "efas_rnet_100km_01min")
out1$latlong <- paste(round(out1$Var1, 4), round(out1$Var2, 4), sep = " ")
outhybas_eu  <- inner_join(out1, outhybas_raw, by = "latlong")

# Upstream-catchment filter
# UpCat    <- st_read(dsn = paste0(regimeDir, "Upstreamgroups.shp"), quiet = TRUE)
# length((which(UpCat$up == "Upstream Catchments")))
# upup     <- UpCat[which(UpCat$up == "Upstream Catchments"), ]
# hybas2uc <- match(upup$outlets_x, outhybas_eu$outlets.y)
# outhybas <- outhybas_eu[hybas2uc, ]
# 
matcat <- match(outhybas_eu$latlong, UpArea_full$latlong)
UpArea <- UpArea_full[matcat, ]

# Column selector: positions in CSV files matching our outlet IDs
# (read one header row to get the column names)
header_ref <- fread(paste0(hydroDir, "tss/HERA_Histo/SnowUpsX_1951_2020.csv"),
                    nrows = 0, header = TRUE)
matcol <- match(UpArea_full$outlets, as.numeric(colnames(header_ref)))
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
                                    resarea_lkp,
                                    Q=FALSE) {
  
  res_dt    <- copy(dt_raw)
  col_names <- names(dt_raw)            # outlet IDs present in this CSV

  # Only process outlets that actually appear in the data
  outlets_in_data <- intersect(names(imm_children), col_names)
  names(imm_children)
  #outlets_in_data[11]
  for (p_id in outlets_in_data) {
    #p_id="307087"
    #p_id="236676"
    #print(p_id)
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
    # if the variable is discharge is m3/s transfer to m3 for the mass balance
    if (Q==T){
      #factor is 3600*24
      factor <- 86400
      child_vols <- dt_raw[, .SD * (factor), .SDcols = valid_children]
      tot_vol=dt_raw[[p_id]] * factor
    }else{
      # child volume matrix: each column = child_series * child_area
      child_vols <- dt_raw[, mapply(`*`, .SD, child_areas),
                           .SDcols = valid_children]
      # head(child_vols)
      tot_vol=dt_raw[[p_id]] * p_area
    }
    # head(tot_vol)
    total_child_vol <- if (length(valid_children) > 1) {
      rowSums(child_vols, na.rm = TRUE)
    } else {
      as.vector(child_vols)
    }
  
    # head(total_child_vol)
    # merdd=data.frame(tot_vol[1:10000],tot_vol[1:10000]-total_child_vol[1:10000],tot_vol[1:10000]/total_child_vol[1:10000])
    # # 
    # plot(tot_vol[1:10000],total_child_vol[1:10000],type="p")
    # points(total_child_vol[1:10000],col=2)
    if (Q==T){
      res_series <- (tot_vol - total_child_vol) / factor
    }else{
      res_series <- (tot_vol - total_child_vol) / res_area
    }
    # 

    # plot(dt_raw[[p_id]][1:10000],col=2)
    # plot(res_series[1:10000])
    # (res_series[[1]])
    #print(mean((res_series[[1]]),na.rm=T))
    if(mean(res_series[[1]],na.rm=T)>1000) print(pid)
    set(res_dt, j = p_id, value = res_series)
  }

  return(res_dt)
}



#check for disaggregation
catchments=catchments_gpkg
# Compare polygon area vs LISFLOOD upstream area for all catchments
area_comparison <- data.frame(
  catch_id   = catchments$catch_id,
  area_poly  = catchments$area_km2,
  area_grid  = UpArea$upa[match(catchments$catch_id, UpArea$outlets)]
)
area_comparison$ratio <- area_comparison$area_poly / area_comparison$area_grid
summary(area_comparison$ratio)
# Should be close to 1.0 for all catchments
# Large deviations (< 0.9 or > 1.1) flag problematic outlets

#If the ratio deviates significantly, 
#replace area_km2 with the grid-based upstream area (upa) in the mass balance, since that is what LISFLOOD actually used:

area_lookup    <- setNames(UpArea$upa,
                           as.character(UpArea$outlets))
resarea_lookup <- setNames(
  UpArea$upa[match(as.numeric(catchments$catch_id), UpArea$outlets)] -
    sapply(as.character(catchments$catch_id), function(id) {
      children <- imm_children_list[[id]]
      if (length(children) == 0) return(0)
      sum(UpArea$upa[match(children, UpArea$outlets)], na.rm = TRUE)
    }),
  as.character(catchments$catch_id)
)

# All from grid-based upa
resarea_lookup <- setNames(
  sapply((catchments$catch_id), function(p_id) {
    p_upa    <- UpArea$upa[match(as.numeric(p_id), (UpArea$outlets))]
    children <- imm_children_list[[p_id]]
    if (length(children) == 0) return(p_upa)
    children_upa <- sum(UpArea$upa[match(children, UpArea$outlets)], na.rm = TRUE)
    p_upa - children_upa
  }),
  as.character(catchments$catch_id)
)

# For each nested outlet, check what fraction of child area is actually covered
coverage <- sapply(names(imm_children_list), function(p_id) {
  children       <- imm_children_list[[p_id]]
  if (length(children) == 0) return(NA)
  covered        <- intersect(children, cnames)
  area_covered   <- sum(area_lookup[covered],   na.rm = TRUE)
  area_total     <- sum(area_lookup[children],  na.rm = TRUE)
  area_covered / area_total
})
# Any value < 1.0 means some child area is unaccounted for
low_coverage <- which(coverage < 1)
cat(length(low_coverage), "outlets have < 95% child area coverage\n")

catchments$residual_fraction <- resarea_lookup[as.character(catchments$catch_id)] /
  area_lookup[as.character(catchments$catch_id)]

# Outlets where the residual is < 5% of the parent area are unreliable
plot(catchments$residual_fraction)


negative_res <- names(resarea_lookup)[resarea_lookup < 0]


for (p_id in negative_res) {
  id_p     <- which(catchments$catch_id == as.numeric(p_id))
  children <- imm_children_list[[p_id]]
  
  # Union of all children geometries
  child_idx   <- which(catchments$catch_id %in% as.numeric(children))
  child_union <- st_make_valid(st_union(catchments[child_idx, ]))
  parent_geom <- st_make_valid(st_geometry(catchments[id_p, ]))
  
  # How much of the child union falls OUTSIDE the parent?
  outside      <- st_difference(child_union, parent_geom)
  outside=st_make_valid(outside)
  outside_area <- as.numeric(st_area(outside)) / 1e6
  
  cat(sprintf("Outlet %s: %.1f km2 of children lie outside parent polygon\n",
              p_id, outside_area))
}

unreliable <- catchments$catch_id[catchments$residual_fraction < 0.01]
cat(length(unreliable), "outlets have residual area < 5% of parent\n")

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

# pid=which(names(Rain_res)=="235555")
# imm_children_list["235555"]
# Rain_res <- Rain_res[, lapply(.SD, function(x) pmax(0, x))]
# sample=as.numeric(unlist(Rain_res[c(1:10000),"235555"]))
# plot(c(1:10000),sample)
# sampleR=as.numeric(unlist(Rain_raw[c(1:10000),"235555"]))
# points(c(1:10000),sampleR,col=2)


message("Deaggregating Snow ...")
Snow_res <- deaggregate_to_residual(Snow_raw, imm_children_list,
                                    area_lookup, resarea_lookup)


Snow_res <- Snow_res[, lapply(.SD, function(x) pmax(0, x))]

Precipitation_res <- Rain_res + Snow_res
Snowfraction_res  <- Snow_res / Precipitation_res   # NaN where P = 0 is OK

# rm(Rain_raw, Snow_raw, Rain_res, Snow_res); gc()





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
save(SaveAgg_SF, file = paste0(regimeDir, "SnowFraction_nested2.Rdata"))

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
# Average yearly precipitation map (residual catchments)                      #
# --------------------------------------------------------------------------- #

# -- 1. Compute average yearly precipitation per outlet --------------------
#    Rain_res is in mm per time step; sum to yearly then average across years

rain_dt <- copy(Rain_res)
rain_dt[, year := year(timeStampX)]

# Sum within each year, then average across years -> mean annual precipitation
yearly_rain <- rain_dt[, lapply(.SD, sum, na.rm = TRUE),
                       by = year,
                       .SDcols = cnames]

mean_annual_rain <- colMeans(yearly_rain[, .SD, .SDcols = cnames],
                             na.rm = TRUE)


snow_dt <- copy(Snow_res)
snow_dt[, year := year(timeStampX)]

# Sum within each year, then average across years -> mean annual precipitation
yearly_snow <-  snow_dt[, lapply(.SD, sum, na.rm = TRUE),
                       by = year,
                       .SDcols = cnames]

mean_annual_snow <- colMeans(yearly_snow[, .SD, .SDcols = cnames],
                             na.rm = TRUE)

# -- 2. Attach to spatial points --------------------------------------------
# Reuse ppoints_res built in Section 7 (residual centroids, EPSG:3035)
# ppoints_res$mean_precip <- mean_annual_rain[
#   match(as.character(ppoints_res$catch_id), names(mean_annual_rain))
# ]

catchments_plot$mean_rain<- mean_annual_rain[
  match(as.character(catchments_plot$catch_id), names(mean_annual_rain))
]

catchments_plot$mean_snow <- mean_annual_snow[
  match(as.character(catchments_plot$catch_id), names(mean_annual_snow))
]
catchments_plot$mean_precip <- catchments_plot$mean_rain+catchments_plot$mean_snow
# -- 3. Plot ----------------------------------------------------------------
palet <- hcl.colors(11, palette = "Roma", rev = F)

ggplot(basemap) +
  geom_sf(fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data  = catchments_plot,
          aes(geometry = geom,
              fill     = mean_precip),
          color = "grey20", shape = 21, stroke = 0.2, alpha = 0.9) +
  scale_fill_gradientn(
    colors = palet,
    oob    = scales::squish,
    name   = "Mean annual\nprecipitation (mm/y)",
    breaks = seq(0, 2000, by = 200),
    limits=c(0,1500)
  ) +
  # scale_size(
  #   range  = c(0.8, 5),
  #   trans  = "sqrt",
  #   name   = expression(paste("Upstream area (", km^2, ")")),
  #   breaks = c(100, 1000, 10000, 100000, 500000),
  #   labels = c("100", "1 000", "10 000", "100 000", "500 000")
  # ) +
  coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
  labs(
    title    = "Mean annual precipitation – residual catchments",
    subtitle = "Average over 1951–2020, deaggregated to inter-catchment area",
    x = "Longitude", y = "Latitude"
  ) +
  guides(
    fill = guide_colourbar(barwidth = 0.5, barheight = 18, reverse = FALSE),
    size = "none"
  ) +
  theme(
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

ggsave(paste0(regimeDir, "plots/mean_annual_precipitation_vf.png"),
       width = 22, height = 16, units = "cm", dpi = 300)

# --------------------------------------------------------------------------- #
# 8.  Plot: Snow Fraction (residual)                                          #
# --------------------------------------------------------------------------- #
ys1 <- 1951:2020
my  <- which(!is.na(match(SaveAgg_SF$SaveAgg_SF, ys1)))
SF_means              <- colMeans(SaveAgg_SF[my, -1], na.rm = TRUE)
catchments_plot$snowfraction <- SF_means[km_res]

palet <- hcl.colors(9, palette = "YlGnBu", rev = TRUE)

ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = catchments_plot,
          aes(geometry = geom, fill = snowfraction),
          color = "black", shape = 21, stroke = 0) +
  geom_sf(fill = NA, color = "grey20") +
  scale_x_continuous(breaks = seq(-30, 40, by = 5)) +
  # scale_size(range = c(1, 4), trans = "sqrt",
  #            name = expression(paste("Upstream area ", (km^2))),
  #            breaks = c(101, 1000, 10000, 100000, 500000),
  #            labels = c("100", "1000", "10 000", "100 000", "500 000")) +
  scale_fill_gradientn(colors = palet, oob = scales::squish,
                       name = "Snow fraction\n(residual catchment)",
                       breaks = seq(0, 1, 0.1), limits = c(0, 0.5)) +
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

ggsave(paste0(regimeDir, "plots/snowfraction_map_f.png"),
       width = 20, height = 15, units = "cm", dpi = 500)

rm(Rain_raw,Snow_raw,rain_dt,snow_dt)
gc()
# --------------------------------------------------------------------------- #
# 9.  AET and PET: load, deaggregate, compute yearly metrics                 #
# --------------------------------------------------------------------------- #
message("Loading and deaggregating AET and PET ...")

ActEvapo_raw <- fread(paste0(hydroDir, "tss/HERA_Histo/ActEvapo_1951_2020.csv"),
                      header = TRUE)
colnames(ActEvapo_raw)[1]="time"
time       <- ActEvapo_raw$time
timeStampX <- time[order(time)]
ActEvapo_raw <- ActEvapo_raw[order(time), .SD, .SDcols = matcol]
ActEvapo_res <- deaggregate_to_residual(dt_raw=ActEvapo_raw, imm_children=imm_children_list,
                                        area_lkp=area_lookup, resarea_lkp=resarea_lookup)
ActEvapo_res <- ActEvapo_res[, lapply(.SD, function(x) pmax(0, x))]
rm(ActEvapo_raw)
gc()

PEvapo_raw   <- fread(paste0(hydroDir, "tss/HERA_Histo/etUpsX_1951_2020.csv"),
                      header = TRUE)
PEvapo_raw   <- PEvapo_raw[order(time), .SD, .SDcols = matcol]

PEvapo_res   <- deaggregate_to_residual(PEvapo_raw, imm_children_list,
                                        area_lookup, resarea_lookup)
PEvapo_res   <- PEvapo_res[, lapply(.SD, function(x) pmax(0, x))]

rm(PEvapo_raw); gc()


#continue

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

save(SavePrecip_res,   file = paste0(regimeDir, "PrecipitationMetrics_nested_f.Rdata"))
save(SeasonPrecip_res, file = paste0(regimeDir, "PrecipitationSeason_nested_f.Rdata"))
save(SaveAET_res,      file = paste0(regimeDir, "AETMetrics_nested_f.Rdata"))

rm(ActEvapo_res, PEvapo_res, Precipitation_res, Snowfraction_res); gc()
gc()
# --------------------------------------------------------------------------- #
# 11. Plot: AET (residual)                                                    #
# --------------------------------------------------------------------------- #
SaveAET_res_df <- data.frame(SaveAET_res)
my             <- which(!is.na(match(SaveAET_res_df$SaveAET_res, ys1)))
KAET_means     <- colMeans(SaveAET_res_df[my, -1], na.rm = TRUE)
catchments_plot$AET1 <- KAET_means[km_res]

palet <- hcl.colors(9, palette = "RdYlBu", rev = FALSE)

ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = catchments_plot,
          aes(geometry = geom, fill = AET1),
          color = "black", shape = 21, stroke = 0) +
  geom_sf(fill = NA, color = "grey20") +
  scale_x_continuous(breaks = seq(-30, 40, by = 5)) +
  # scale_size(range = c(1, 4), trans = "sqrt",
  #            name = expression(paste("Upstream area ", (km^2))),
  #            breaks = c(101, 1000, 10000, 100000, 500000),
  #            labels = c("100", "1000", "10 000", "100 000", "500 000")) +
  scale_fill_gradientn(colors = palet, oob = scales::squish,
                       name = "AET (mm/y)\n(residual catchment)", limits =c(300,700)) +
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
catchments_plot$AET.P1 <- KeepPag$AET.P[mai_res]

palet <- hcl.colors(9, palette = "RdYlBu", rev = TRUE)

ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = catchments_plot,
          aes(geometry = geom, fill = AET.P1),
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

ggsave(paste0(regimeDir, "plots/AridityIndex_nested_f.png"),
       width = 20, height = 15, units = "cm", dpi = 1500)

# --------------------------------------------------------------------------- #
# 13. Rainfall Gini concentration index (residual)                            #
# --------------------------------------------------------------------------- #
KeepG <- SavePrecip_res[myp, c(1, 3, 4)]

KeepGag <- stats::aggregate(list(px = KeepG$val.gini),
                            by  = list(catch = KeepG$cat),
                            FUN = function(x) mean(x, na.rm = TRUE))
KeepGag <- do.call(data.frame, KeepGag)

mai_g             <- match(as.character(catchments_plot$catch_id), KeepGag$catch)
catchments_plot$gini1 <- KeepGag$px[mai_g]

ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = catchments_plot,
          aes(geometry = geom, fill = gini1),
          color = "black", shape = 21, stroke = 0.1) +
  geom_sf(fill = NA, color = "grey20") +
  scale_x_continuous(breaks = seq(-30, 40, by = 5)) +
  scale_size(range = c(1, 4), trans = "sqrt",
             name = expression(paste("Upstream area ", (km^2))),
             breaks = c(101, 1000, 10000, 100000, 500000),
             labels = c("100", "1000", "10 000", "100 000", "500 000")) +
  scale_fill_gradientn(colors = palet, oob = scales::squish,
                       name = "Precip. GINI\n(residual catchment)",
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

mai_si             <- match(as.character(catchments_plot$catch_id), SeasonPrecipSI_res$cat)
catchments_plot$PSI    <- SeasonPrecipSI_res$Seasonality_1[mai_si]
catchments_plot$PSIch  <- SeasonPrecipSI_res$SIchange[mai_si]

# --- Seasonality level ---
ggplot(basemap) +
  geom_sf(fill = "white", color = NA) +
  geom_sf(data = catchments_plot,
          aes(geometry = geom, fill = PSI),
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
  geom_sf(data = catchments_plot,
          aes(geometry = geom, fill = PSIch),
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
outputfilenames <- c(
 "disWin"
)
Q=FALSE
for (var in outputfilenames) {
  var=outputfilenames
  if(var=="disWin") Q=TRUE
  message("Batch deaggregating: ", var)

  in_file <- paste0(hydroDir, "tss/HERA_Histo/", var, "_1951_2020.csv")
  if (!file.exists(in_file)) {
    message("  File not found, skipping: ", in_file)
    next
  }

  Vari      <- fread(in_file, header = TRUE)
  time=Vari$V1
  Vari      <- Vari[order(Vari$V1), ]
  matcol_v  <- match(UpArea$outlets, as.numeric(colnames(Vari)))
  matcol_v  <- matcol_v[!is.na(matcol_v)]
  Vari      <- Vari[, .SD, .SDcols = matcol_v]

  Vari_res  <- deaggregate_to_residual(Vari, imm_children_list,
                                       area_lookup, resarea_lookup,Q)
  if(Q==FALSE){
    Vari_res  <- Vari_res[, lapply(.SD, function(x) pmax(0, x))]
  }
  Vari_res <- data.table(time[order(time)],Vari_res)
  nik=Vari_res[["307087"]]
  mean(nik)
  fwrite(Vari_res,
         paste0(regimeDir, "data/", var, "_nested_1951_2020.csv"),
         sep = ",", na = "-9999", row.names = FALSE,
         quote = FALSE, showProgress = TRUE, nThread = 2)

  rm(Vari, Vari_res); gc()
}

message("Statevar_ana_nested.R completed successfully.")
