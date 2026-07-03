# --------------------------------------------------------------------------- #
# process_hydro_variable()
# --------------------------------------------------------------------------- #
# Reads a nested CSV, computes mean annual and monthly aggregates,
# saves monthly aggregates to disk, and optionally plots a map.
#
# ARGUMENTS
#   var_name      : short label used in plot titles and output filenames
#   file_path     : full path to the CSV file
#   cnames        : character vector of outlet IDs (column names to keep)
#   catchments_plot: sf object with catch_id and geom columns
#   basemap       : world basemap sf object
#   nco           : coordinate matrix for coord_sf limits
#   regimeDir     : root output directory
#   palette       : hcl.colors palette name
#   rev_palette   : logical, reverse palette
#   legend_label  : string for the fill legend title
#   plot_limits   : numeric(2), c(min, max) for scale_fill_gradientn limits
#   plot_breaks   : numeric vector of legend breaks
#   time_col      : name of the timestamp column in the CSV (default "V1")
#   save_plot     : logical, whether to save the map (default TRUE)
#   time_override : optional POSIXct vector if time must be parsed externally
#
# RETURNS
#   list(mean_annual, monthly_means) — both as named numeric vectors / data.table
# --------------------------------------------------------------------------- #



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

regimeDir="D:/tilloal/Documents/01_Projects/RegimeShifts/"

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
cnames <- (UpArea$outlets)   # character IDs of our outlets

# Match gpkg rows to our outlet list so point order is consistent
gpkg_in_outlets <- match(cnames, as.numeric(catchments_gpkg$catch_id))
catchments_plot  <- catchments_gpkg[gpkg_in_outlets[!is.na(gpkg_in_outlets)], ]

message("Loading Rain ...")
Rain_raw <- fread(paste0(regimeDir, "data/tss/HERA_SocCF/rainUpsX_1951_2020.csv"),
                  header = TRUE)
time       <- Rain_raw$V1
timeStampX <- time[order(time)]


process_hydro_variable <- function(var_name,
                                   file_path,
                                   cnames,
                                   catchments_plot,
                                   time,
                                   basemap,
                                   nco,
                                   regimeDir,
                                   palette       = "Blues 3",
                                   rev_palette   = TRUE,
                                   legend_label  = NULL,
                                   plot_limits   = NULL,
                                   plot_breaks   = waiver(),
                                   time_col      = "V1",
                                   save_plot     = TRUE,
                                   sm            = FALSE,
                                   time_override = NULL) {
  
  message("\n=== Processing: ", var_name, " ===")
  
  # ---- 1. Load -----------------------------------------------------------
  message("  Loading ", basename(file_path), " ...")
  dt <- fread(file_path, header = TRUE)
  
  # Parse time — use override if provided, otherwise use passed-in time vector
  if (!is.null(time_override)) {
    ts <- time_override
  } else {
    ts <- time
  }
  dt  <- data.table(time, dt)
  dt  <- dt[order(ts), ]
  ts  <- ts[order(ts)]
  
  # Keep only outlet columns present in cnames
  cols_present <- intersect(cnames, names(dt))
  missing_cols <- setdiff(cnames, names(dt))
  if (length(missing_cols) > 0) {
    message("  WARNING: ", length(missing_cols), " outlets missing from file")
  }
  dt <- dt[, .SD, .SDcols = cols_present]
  
  
  anchor = as.POSIXct("1978-11-01")   # month-day string, default keeps current behaviour
 
  # ---- 2. Add all time index columns once --------------------------------
  dt_time <- copy(dt)
  
  # # Basic date components
  dt_time[, `:=`(
    year   = year(time),
    month  = month(time),
    date   = as.Date(time)
  )]
  
  # Fixed‑epoch indices relative to anchor_date
  day_offset <- ((as.Date(ts) - as.Date(anchor)))
  # win7   = floor(day_offset / (7))
  # plot(win7[1:100])
  dt_time[, `:=`(
    win7   = floor(day_offset / 7),                     # 7‑day window index
    win14  = floor(day_offset / 14),                    # 14‑day window index
    mon    = (year - year(anchor)) * 12L +
      (month - month(anchor))                    # month offset
  )]
  # dt_time[, `:=`(
  #   year   = year(ts),
  #   month  = month(ts),
  #   date   = as.Date(ts),
  #   week7  = paste(year(ts),
  #                  formatC(ceiling(yday(ts) / 7),  width = 2, flag = "0"),
  #                  sep = "-"),
  #   week14 = paste(year(ts),
  #                  formatC(ceiling(yday(ts) / 14), width = 2, flag = "0"),
  #                  sep = "-")
  # )]
  
  out_dir <- paste0(regimeDir, "data/aggregates/", var_name, "/")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Aggregation function: sum for fluxes, mean for state variables (sm=TRUE)
  agg_fun <- if (sm) mean else sum
  
  # ---- 3. Daily ----------------------------------------------------------
  message("  Computing daily aggregates ...")
  
  
  
  daily_totals <- dt_time[,
                          lapply(.SD, agg_fun, na.rm = TRUE),
                          by = .(year, date),
                          .SDcols = cols_present
  ]
  daily_totals[, doy := yday(date)]   # keep day-of-year for means
  
  daily_means <- daily_totals[,
                              lapply(.SD, mean, na.rm = TRUE),
                              by = doy,
                              .SDcols = cols_present
  ]
  setorder(daily_means, doy)
  

  
  fwrite(daily_means,
         paste0(out_dir, var_name, "_daily_means.csv"),
         sep = ",", na = "-9999")
  fwrite(daily_totals,
         paste0(out_dir, var_name, "_daily_all_years.csv"),
         sep = ",", na = "-9999")
  message("  Saved daily means and all-years")
  
  # ---- 4. 7-day ----------------------------------------------------------
  message("  Computing 7-day aggregates ...")
  
  
  
  
  week7_totals <- dt_time[,
                          lapply(.SD, agg_fun, na.rm = TRUE),
                          by = win7,
                          .SDcols = cols_present
  ]
  #week7_totals[, window := as.integer(sub(".*-", "", win7))]
  
  setnames(week7_totals, "win7", "window")
  setorder(week7_totals, window)
  #week7_totals$start=as.Date(anchor) + week7_totals$window*7
  
  # Add start/end dates for each 7‑day window
  week7_totals[, period_start := as.Date(anchor) + window * 7L]
  week7_totals[, period_end   := period_start + 6L]   # 7‑day window
  
  setcolorder(week7_totals, c("window", "period_start", "period_end", cols_present))
  
  
  # week7_means <- week7_totals[,
  #                             lapply(.SD, mean, na.rm = TRUE),
  #                             by = window,
  #                             .SDcols = cols_present
  # ]
  # setorder(week7_means, window)
  
  # fwrite(week7_means,
  #        paste0(out_dir, var_name, "_7day_means.csv"),
  #        sep = ",", na = "-9999")
  fwrite(week7_totals,
         paste0(out_dir, var_name, "_7day_all_years.csv"),
         sep = ",", na = "-9999")
  message("  Saved 7-day means and all-years")
  
  # ---- 5. 14-day ---------------------------------------------------------
  message("  Computing 14-day aggregates ...")
  
  week14_totals <- dt_time[,
                           lapply(.SD, agg_fun, na.rm = TRUE),
                           by =  win14,
                           .SDcols = cols_present
  ]
  
  setnames(week14_totals, "win14", "window")
  setorder(week14_totals, window)
  #week7_totals$start=as.Date(anchor) + week7_totals$window*7
  
  # Add start/end dates for each 7‑day window
  week14_totals[, period_start := as.Date(anchor) + window * 14L]
  week14_totals[, period_end   := period_start + 13L]  
  #week14_totals[, window := as.integer(sub(".*-", "", week14))]
  
  setcolorder(week14_totals, c("window", "period_start", "period_end", cols_present))
  
  # week14_means <- week14_totals[,
  #                               lapply(.SD, mean, na.rm = TRUE),
  #                               by = window,
  #                               .SDcols = cols_present
  # ]
  # setorder(week14_means, window)
  
  # fwrite(week14_means,
  #        paste0(out_dir, var_name, "_14day_means.csv"),
  #        sep = ",", na = "-9999")
  fwrite(week14_totals,
         paste0(out_dir, var_name, "_14day_all_years.csv"),
         sep = ",", na = "-9999")
  message("  Saved 14-day means and all-years")
  
  # ---- 6. Monthly --------------------------------------------------------
  message("  Computing monthly aggregates ...")
  
  monthly_totals <- dt_time[,
                            lapply(.SD, agg_fun, na.rm = TRUE),
                            by = mon,
                            .SDcols = cols_present
  ]
  setnames(monthly_totals, "mon", "month_idx")
  setorder(monthly_totals, month_idx)
  
  # Month start/end using lubridate’s safe month arithmetic
  monthly_totals[, period_start := anchor %m+% months(month_idx)]
  monthly_totals[, period_end   := (anchor %m+% months(month_idx + 1L)) - days(1)]
  
  setcolorder(monthly_totals, c("month_idx", "period_start", "period_end", cols_present))
 
  monthly_totals$period_start=as.Date(monthly_totals$period_start+3600)
  monthly_totals$period_end=as.Date(monthly_totals$period_end+3600)
  

  
  monthly_means <- monthly_totals[,
                                  lapply(.SD, mean, na.rm = TRUE),
                                  by = month_idx,
                                  .SDcols = cols_present
  ]
  setorder(monthly_means, month_idx)
  
  
  fwrite(monthly_means,
         paste0(out_dir, var_name, "_monthly_means.csv"),
         sep = ",", na = "-9999")
  fwrite(monthly_totals,
         paste0(out_dir, var_name, "_monthly_all_years.csv"),
         sep = ",", na = "-9999")
  message("  Saved monthly means and all-years")
  
  # ---- 7. Mean annual ----------------------------------------------------
  yearly_totals <- dt_time[,
                           lapply(.SD, agg_fun, na.rm = TRUE),
                           by = year,
                           .SDcols = cols_present
  ]
  mean_annual <- colMeans(yearly_totals[, .SD, .SDcols = cols_present],
                          na.rm = TRUE)
  
  # ---- 8. Attach to catchments_plot and map ------------------------------
  col_out <- paste0("mean_", var_name)
  catchments_plot[[col_out]] <- mean_annual[
    match(as.numeric(catchments_plot$catch_id), as.numeric(names(mean_annual)))
  ]
  
  if (save_plot) {
    message("  Plotting map ...")
    
    if (is.null(legend_label)) legend_label <- paste0("Mean annual\n", var_name, " (mm/y)")
    palet <- hcl.colors(11, palette = palette, rev = rev_palette)
    
    p <- ggplot(basemap) +
      geom_sf(fill = "grey95", color = "grey70", linewidth = 0.2) +
      geom_sf(data  = catchments_plot,
              aes(geometry = geom, fill = .data[[col_out]]),
              color = "grey20", shape = 21, stroke = 0.2, alpha = 0.9) +
      scale_fill_gradientn(
        colors = palet,
        oob    = scales::squish,
        name   = legend_label,
        breaks = plot_breaks,
        limits = plot_limits
      ) +
      coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
      labs(
        title    = paste("Mean annual", var_name),
        subtitle = "Average over 1951-2020, deaggregated to inter-catchment area",
        x = "Longitude", y = "Latitude"
      ) +
      guides(
        fill = guide_colourbar(barwidth = 0.5, barheight = 18, reverse = FALSE)
      ) +
      theme(
        plot.title       = element_text(size = 14, face = "bold"),
        plot.subtitle    = element_text(size = 11, colour = "grey40"),
        axis.title       = element_text(size = 12),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.title     = element_text(size = 12),
        legend.text      = element_text(size = 10),
        legend.position  = "right",
        legend.key       = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size  = unit(0.8, "cm")
      )
    
    ggsave(paste0(regimeDir, "plots/mean_", var_name, ".png"),
           p, width = 22, height = 16, units = "cm", dpi = 300)
    message("  Saved plot -> mean_", var_name, ".png")
  }
  
  return(invisible(list(
    mean_annual    = mean_annual,
    daily_means    = daily_means,
    daily_totals   = daily_totals,
   
    week7_totals   = week7_totals,

    week14_totals  = week14_totals,
    monthly_means  = monthly_means,
    monthly_totals = monthly_totals
  )))
}


# --------------------------------------------------------------------------- #
# Run for all variables
# --------------------------------------------------------------------------- #

# Upper zone outflow
quz_out <- process_hydro_variable(
  var_name      = "quz",
  file_path     = paste0(regimeDir, "/data/qUzUpsX_nested_1951_2020.csv"),
  cnames        = cnames,
  catchments_plot = catchments_plot,
  time          = timeStampX,
  basemap       = basemap,
  nco           = nco,
  regimeDir     = regimeDir,
  palette       = "Blues 3",
  rev_palette   = TRUE,
  legend_label  = "Mean annual\nUZ outflow (mm/y)",
  plot_limits   = c(0, 500)
)

# Percolation upper to lower zone
qutl_out <- process_hydro_variable(
  var_name      = "qutl",
  file_path     = paste0(regimeDir, "/data/percUZLZUpsX_nested_1951_2020.csv"),
  cnames        = cnames,
  catchments_plot = catchments_plot,
  time          = timeStampX,
  basemap       = basemap,
  nco           = nco,
  regimeDir     = regimeDir,
  palette       = "Teal",
  rev_palette   = FALSE,
  legend_label  = "Mean annual\nUZ->LZ perc. (mm/y)",
  plot_limits   = c(0, 300)
)

# Lower zone outflow
qlz_out <- process_hydro_variable(
  var_name      = "qlz",
  file_path     = paste0(regimeDir, "/data/qLzUpsX_nested_1951_2020.csv"),
  cnames        = cnames,
  catchments_plot = catchments_plot,
  time          = timeStampX,
  basemap       = basemap,
  nco           = nco,
  regimeDir     = regimeDir,
  palette       = "YlGnBu",
  rev_palette   = TRUE,
  legend_label  = "Mean annual\nLZ outflow (mm/y)",
  plot_limits   = c(0, 400)
)

# Discharge out 
Q_out <- process_hydro_variable(
  var_name      = "Q",
  file_path     = paste0(regimeDir, "/data/disWin_nested_1951_2020.csv"),
  cnames        = cnames,
  catchments_plot = catchments_plot,
  time          = timeStampX,
  basemap       = basemap,
  nco           = nco,
  regimeDir     = regimeDir,
  palette       = "YlGnBu",
  rev_palette   = TRUE,
  sm=T,
  legend_label  = "Mean annual\nQ outflow (m3/s)",
  plot_limits   = c(0, 100)
)

# # Infiltration — note different time vector
# timeI   <- as.POSIXct(fread(paste0(hydroDir, "/tss/HERA_Histo/InfUpsX_1951_2020_v2.csv"),
#                             select = "V1")$V1)

infil_out <- process_hydro_variable(
  var_name      = "infiltration",
  file_path     = paste0(regimeDir, "/data/InfUpsX_nested_1951_2020_v2.csv"),
  cnames        = cnames,
  catchments_plot = catchments_plot,
  time          = timeStampX,
  basemap       = basemap,
  nco           = nco,
  regimeDir     = regimeDir,
  palette       = "YlGnBu",
  rev_palette   = FALSE,
  legend_label  = "Mean annual\ninfiltration (mm/y)",
  plot_limits   = c(0, 600)     # pass the externally parsed time for this file
)


# Runoff — from nested residual CSV
aet_out <- process_hydro_variable(
  var_name      = "ActEvapo",
  file_path     = paste0(regimeDir, "/data/ActEvapo_nested_1951_2020.csv"),
  cnames        = cnames,
  catchments_plot = catchments_plot,
  time          = timeStampX,
  basemap       = basemap,
  nco           = nco,
  regimeDir     = regimeDir,
  palette       = "Roma",
  rev_palette   = FALSE,
  legend_label  = "Mean annual\nAET (mm/y)",
  plot_limits   = c(0, 600),
  plot_breaks   = seq(0, 600, by = 100)
)

# soil moisture upper layer
th2_out <- process_hydro_variable(
  var_name      = "root_soil_moisture",
  file_path     = paste0(regimeDir, "/data/theta2totalX_nested_1951_2020.csv"),
  cnames        = cnames,
  catchments_plot = catchments_plot,
  time          = timeStampX,
  basemap       = basemap,
  nco           = nco,
  regimeDir     = regimeDir,
  palette       = "RdYlBu",
  rev_palette   = FALSE,
  sm            = TRUE,
  legend_label  = "Mean annual soil moisture (mm3/mm3)",
  plot_limits   = c(0.1, 0.4)
)


# Runoff — from nested residual CSV
runoff_out <- process_hydro_variable(
  var_name      = "runoff",
  file_path     = paste0(regimeDir, "/data/surfaceRunoffUpsX_nested_1951_2020.csv"),
  cnames        = cnames,
  catchments_plot = catchments_plot,
  time          = timeStampX,
  basemap       = basemap,
  nco           = nco,
  regimeDir     = regimeDir,
  palette       = "Roma",
  rev_palette   = FALSE,
  legend_label  = "Mean annual\nrunoff (mm/y)",
  plot_limits   = c(0, 100),
  plot_breaks   = seq(0, 100, by = 20)
)

# soil moisture upper layer
th1_out <- process_hydro_variable(
  var_name      = "surface_soil_moisture",
  file_path     = paste0(regimeDir, "/data/theta1totalX_nested_1951_2020.csv"),
  cnames        = cnames,
  catchments_plot = catchments_plot,
  time          = timeStampX,
  basemap       = basemap,
  nco           = nco,
  regimeDir     = regimeDir,
  palette       = "RdYlBu",
  rev_palette   = FALSE,
  sm            = TRUE,
  legend_label  = "Mean annual soil moisture (mm3/mm3)",
  plot_limits   = c(0.1, 0.4)
)

#snow
snow <- process_hydro_variable(
  var_name      = "snow_water_equivalent",
  file_path     = paste0(regimeDir, "/data/scovUps_nested_1951_2020.csv"),
  cnames        = cnames,
  catchments_plot = catchments_plot,
  time          = timeStampX,
  basemap       = basemap,
  nco           = nco,
  regimeDir     = regimeDir,
  palette       = "RdYlBu",
  rev_palette   = FALSE,
  sm            = TRUE,
  legend_label  = "Mean Snow Water Equivalent (mm)",
  plot_limits   = c(0, 100)
)

