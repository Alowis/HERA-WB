# =============================================================================
# Figure 1: 4-panel evolution
#   A: Infiltration ratio — Avignon (307920)
#   B: Sealed area fraction — Avignon (307920)
#   C: Snowmelt (Jan–Jun) — Aquila (325150)
#   D: Summer root zone soil moisture — Aquila (325150)
#
# Figure 1b: Flow components decomposition — Avignon
#
# Figure 2: Continental climate stripes (2 panels)
#   A: Max monthly SWE (area-weighted)
#   B: Min monthly baseflow qlz (area-weighted)
#
# Figure 4: Continental mean AET — LISFLOOD vs GLEAM (mm/year)
# Figure 4b: Total yearly AET volume (km³)
# Figure S1: Spatial AET bias map for a given year
# Figure S2: GLEAM step-change diagnostic
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
library(cowplot)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
agg_dir <- file.path(base_dir, "data", "aggregates")
tss_dir <- file.path(base_dir, "data", "tss_postprocess")
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
landuse_dir <- "D:/tilloal/Documents/06_Floodrivers/landuse/"
out_dir <- file.path(base_dir, "output", "temporal_evolution")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Catchment choices --------------------------------------------------------
catch_avignon <- "307920"
catch_avignon <- "205374"
catch_aquila <- "325150"
catch_aquila <- "284595"
# 192468 france belgique

# --- Helper: assign dates to TSS data ----------------------------------------
assign_tss_dates <- function(dt) {
    n_rows <- nrow(dt)
    n_days <- as.integer(as.Date("2020-12-31") - as.Date("1951-01-01")) + 1
    if (n_rows == n_days) {
        dt[, date := seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")]
    } else if (n_rows == n_days * 4) {
        dt[, day_idx := rep(seq_len(n_days), each = 4)]
        date_seq <- seq.Date(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")
        dt[, date := date_seq[day_idx]]
    } else {
        warning("Unexpected row count (", n_rows, "). Assuming daily.")
        dt[, date := seq.Date(as.Date("1951-01-01"), length.out = n_rows, by = "day")]
    }
    return(dt)
}

# --- Helper: find column name in TSS header -----------------------------------
find_col <- function(catch_id, header_cols) {
    if (catch_id %in% header_cols) {
        return(catch_id)
    }
    x_catch <- paste0("X", catch_id)
    if (x_catch %in% header_cols) {
        return(x_catch)
    }
    stop("Catchment ", catch_id, " not found in CSV header.")
}

# =============================================================================
# FIGURE 1 — PANEL A: Infiltration ratio (Avignon)
# =============================================================================
cat("[Fig1-A] Infiltration ratio — Avignon...\n")

inf_path <- file.path(tss_dir, "infUpsX_nested_1951_2020.csv")
rain_path <- file.path(tss_dir, "rainUpsX_nested_1951_2020.csv")
snow_path <- file.path(tss_dir, "snowUpsX_nested_1951_2020.csv")

inf_header <- fread(inf_path, nrows = 0, header = TRUE)
col_av <- find_col(catch_avignon, names(inf_header))

inf_dt <- fread(inf_path, select = col_av, header = TRUE)
rain_dt <- fread(rain_path, select = col_av, header = TRUE)
snow_dt <- fread(snow_path, select = col_av, header = TRUE)

inf_dt[, infiltration := inf_dt[[1]]]
inf_dt[, precipitation := rain_dt[[1]] + snow_dt[[1]]]
inf_dt <- assign_tss_dates(inf_dt)

# If 6-hourly, aggregate to daily
if ("day_idx" %in% names(inf_dt)) {
    inf_dt <- inf_dt[, .(
        infiltration = mean(infiltration, na.rm = TRUE),
        precipitation = mean(precipitation, na.rm = TRUE)
    ), by = date]
}

# Monthly aggregation then ratio
inf_monthly <- inf_dt[, .(
    infil = mean(infiltration, na.rm = TRUE),
    precip = mean(precipitation, na.rm = TRUE)
), by = .(year = year(date), month = month(date))]
inf_monthly[, date := as.Date(paste(year, month, "15", sep = "-"))]
inf_monthly[, inf_ratio := fifelse(precip > 0, infil / precip, NA_real_)]
inf_monthly[, inf_ratio := frollmean(inf_ratio, 3, align = "center")]

pA <- ggplot(inf_monthly, aes(x = date, y = inf_ratio)) +
    geom_line(color = "steelblue", linewidth = 0.4) +
    geom_smooth(
        method = "lm", se = TRUE,
        color = "darkblue", linewidth = 0.8
    ) +
    labs(
        title = "A) Infiltration ratio — Avignon (307920)",
        x = NULL, y = "Infiltration / Precipitation (-)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

pA
# =============================================================================
# FIGURE 1 — PANEL B: Sealed area fraction (Avignon)
# =============================================================================
cat("[Fig1-B] Sealed fraction — Avignon...\n")

sealed_files <- list.files(landuse_dir,
    pattern = "^fracsealed_European_01min_\\d{4}\\.nc$",
    full.names = TRUE
)

sealed_years <- as.integer(sub(".*_(\\d{4})\\.nc$", "\\1", sealed_files))
sealed_order <- order(sealed_years)
sealed_files <- sealed_files[sealed_order]
sealed_years <- sealed_years[sealed_order]

catchments <- st_read(gpkg_path, quiet = TRUE)
if (!isTRUE(st_crs(catchments) == st_crs(4326))) {
    catchments <- st_transform(catchments, 4326)
}

catch_sub_av <- catchments[as.character(catchments$catch_id) == catch_avignon, ]

sealed_av <- data.table(year = sealed_years, frac_sealed = NA_real_)

for (i in seq_along(sealed_files)) {
    r <- tryCatch(terra::rast(sealed_files[i]), error = function(e) NULL)
    if (is.null(r)) next
    if (nlyr(r) > 1) r <- r[[1]]
    sealed_av$frac_sealed[i] <- exact_extract(r, catch_sub_av, "mean")
}
sealed_av <- sealed_av[!is.na(frac_sealed)]

pB <- ggplot(sealed_av, aes(x = year, y = frac_sealed)) +
    geom_line(color = "firebrick", linewidth = 0.8) +
    geom_point(color = "firebrick", size = 1.2, alpha = 0.6) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    labs(
        title = paste0("B) Sealed land fraction — Avignon (", catch_avignon, ")"),
        x = NULL, y = "Sealed fraction"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

pB

cat("[Fig1-B] Winter root zone soil moisture — Belgium ..\n")

rsm_path <- file.path(
    agg_dir, "surface_soil_moisture",
    "surface_soil_moisture_monthly_all_years.csv"
)
rsm_dt <- fread(rsm_path)
rsm_dt[, date := as.Date(period_start)]

rsm_col <- find_col(catch_avignon, names(rsm_dt))
rsm_catch <- rsm_dt[, .(date, rsm = get(rsm_col))]
# rsm_catch <- rsm_catch[, .month := month(date)]
# plot(rsm_catch,type="o")
rsm_winter <- rsm_catch[month(date) %in% c(1), ]
plot(rsm_winter$rsm)
# Winter months:
rsm_winter <- rsm_catch[month(date) %in% c(1, 2, 3, 4), .(
    summer_rsm = mean(rsm, na.rm = TRUE)
), by = .(year = year(date))]

pB1 <- ggplot(rsm_winter, aes(x = year, y = summer_rsm)) +
    geom_line(color = "sienna", linewidth = 0.6) +
    geom_smooth(
        method = "lm", se = TRUE,
        color = "saddlebrown", linewidth = 0.8
    ) +
    labs(
        title = "D) Winter surface soil moisture (Jan–March) — Aquila (325150)",
        x = NULL, y = "Root zone SM (mm)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))
pB1
# =============================================================================
# FIGURE 1 — PANEL C: Snowmelt (Aquila)
# =============================================================================
cat("[Fig1-C] Snowmelt — Aquila...\n")

melt_path <- file.path(tss_dir, "snowMeltUpsX_nested_1951_2020.csv")
melt_header <- fread(melt_path, nrows = 0, header = TRUE)
col_aq <- find_col(catch_aquila, names(melt_header))

melt_dt <- fread(melt_path, select = col_aq, header = TRUE)
setnames(melt_dt, 1, "snowmelt")
melt_dt <- assign_tss_dates(melt_dt)

# If 6-hourly, aggregate to daily
if ("day_idx" %in% names(melt_dt)) {
    melt_dt <- melt_dt[, .(snowmelt = mean(snowmelt, na.rm = TRUE)), by = date]
}

# Annual sum of snowmelt (Jan-Jun)
melt_annual <- melt_dt[month(date) <= 6, .(
    snowmelt_sum = sum(snowmelt, na.rm = TRUE)
), by = .(year = year(date))]

pC <- ggplot(melt_annual, aes(x = year, y = snowmelt_sum)) +
    geom_line(color = "darkcyan", linewidth = 0.6) +
    geom_smooth(
        method = "lm", se = TRUE,
        color = "darkblue", linewidth = 0.8
    ) +
    labs(
        title = "C) Annual snowmelt (Jan–Jun) — Aquila (325150)",
        x = NULL, y = "Snowmelt (mm)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

pC
# =============================================================================
# FIGURE 1 — PANEL D: Summer root zone soil moisture (Aquila)
# =============================================================================
cat("[Fig1-D] Summer root zone soil moisture — Aquila...\n")

rsm_path <- file.path(
    agg_dir, "root_soil_moisture",
    "root_soil_moisture_monthly_all_years.csv"
)
rsm_dt <- fread(rsm_path)
rsm_dt[, date := as.Date(period_start)]

rsm_col <- find_col(catch_aquila, names(rsm_dt))
rsm_catch <- rsm_dt[, .(date, rsm = get(rsm_col))]

# Summer months: Jun–Sep
rsm_summer <- rsm_catch[month(date) %in% 6:9, .(
    summer_rsm = mean(rsm, na.rm = TRUE)
), by = .(year = year(date))]

pD <- ggplot(rsm_summer, aes(x = year, y = summer_rsm)) +
    geom_line(color = "sienna", linewidth = 0.6) +
    geom_smooth(
        method = "lm", se = TRUE,
        color = "saddlebrown", linewidth = 0.8
    ) +
    labs(
        title = "D) Summer root zone soil moisture (Jun–Sep) — Aquila (325150)",
        x = NULL, y = "Root zone SM (mm)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

pD
# =============================================================================
# FIGURE 1 — COMPOSE AND SAVE (Avignon)
# =============================================================================
cat("[Fig1] Composing 4-panel Avignon figure...\n")

fig1 <- plot_grid(pA, pB, pC, pD, ncol = 2, align = "hv")

ggsave(file.path(out_dir, "Figure1_Avignon_evolution.png"),
    fig1,
    width = 14, height = 9, dpi = 200
)
cat("  -> Saved Figure 1 (Avignon).\n")

# =============================================================================
# FIGURE 1b — Flow components decomposition (Avignon)
#   Shows absolute annual mean of qlz, quz, and direct runoff
# =============================================================================
cat("[Fig1b] Flow components decomposition — Avignon...\n")

# Load flow components for Avignon
qlz_path_av <- file.path(agg_dir, "qlz", "qlz_monthly_all_years.csv")
quz_path_av <- file.path(agg_dir, "quz", "quz_monthly_all_years.csv")
runoff_path_av <- file.path(agg_dir, "runoff", "runoff_monthly_all_years.csv")

qlz_dt_av <- fread(qlz_path_av)
quz_dt_av <- fread(quz_path_av)
runoff_dt_av <- fread(runoff_path_av)

qlz_dt_av[, date := as.Date(period_start)]
quz_dt_av[, date := as.Date(period_start)]
runoff_dt_av[, date := as.Date(period_start)]

qlz_col_av <- find_col(catch_avignon, names(qlz_dt_av))
quz_col_av <- find_col(catch_avignon, names(quz_dt_av))
runoff_col_av <- find_col(catch_avignon, names(runoff_dt_av))

bfi_dt <- data.table(
    date = qlz_dt_av$date,
    qlz = qlz_dt_av[[qlz_col_av]],
    quz = quz_dt_av[[quz_col_av]],
    runoff = runoff_dt_av[[runoff_col_av]]
)

# Annual means of each component
flow_annual <- bfi_dt[, .(
    Baseflow = mean(qlz, na.rm = TRUE),
    Interflow = mean(quz, na.rm = TRUE),
    `Direct runoff` = mean(runoff, na.rm = TRUE)
), by = .(year = year(date))]

# Melt to long format for ggplot
flow_long <- melt(flow_annual,
    id.vars = "year",
    variable.name = "component", value.name = "flow_mm"
)

p_flow_abs <- ggplot(flow_long, aes(x = year, y = flow_mm, color = component)) +
    geom_line(linewidth = 0.6) +
    geom_smooth(method = "loess", span = 0.5, se = FALSE, linewidth = 0.9) +
    scale_color_manual(values = c(
        "Baseflow" = "#1b7837",
        "Interflow" = "#5e3c99",
        "Direct runoff" = "#e66101"
    )) +
    labs(
        title = "Flow components — Avignon (307920)",
        subtitle = "Annual mean of each component (absolute values)",
        x = NULL, y = "Flow (mm/month)",
        color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

# Stacked area version for proportions
flow_annual_long <- melt(flow_annual,
    id.vars = "year",
    variable.name = "component", value.name = "flow_mm"
)
flow_annual_long[, component := factor(component,
    levels = c("Direct runoff", "Interflow", "Baseflow")
)]

p_flow_stack <- ggplot(flow_annual_long, aes(x = year, y = flow_mm, fill = component)) +
    geom_area(alpha = 0.8) +
    scale_fill_manual(values = c(
        "Interflow" = "#5e3c99",
        "Baseflow" = "#1b7837",
        "Direct runoff" = "#e66101"
    )) +
    labs(
        title = "Flow partitioning — Avignon (307920)",
        subtitle = "Stacked annual mean components",
        x = NULL, y = "Flow (mm/month)",
        fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

fig1b <- plot_grid(p_flow_abs, p_flow_stack, ncol = 1, align = "v")

ggsave(file.path(out_dir, "Figure1b_flow_components_Avignon.png"),
    fig1b,
    width = 10, height = 9, dpi = 200
)
cat("  -> Saved Figure 1b (flow components).\n")

# =============================================================================
# FIGURE 1_NEW — 6-panel evolution (3 catchments × 2 variables)
#   a: Yearly max daily runoff — Paris (224577)
#   b: Sealed area — Paris (224577)
#   c: Mean yearly soil moisture — Rio Cega (335731)
#   d: Yearly total AET — Rio Cega (335731)
#   e: Apr–Oct snowmelt — Ticino (290666)
#   f: Apr–Oct root zone SM — Ticino (290666)
# =============================================================================
cat("[Fig1_new] Building 6-panel figure...\n")

catch_paris <- "224577"
catch_cega <- "335731"
catch_ticino <- "290666"

# Get upstream areas from catchments gpkg
if (!exists("catchments")) catchments <- st_read(gpkg_path, quiet = TRUE)
get_area <- function(cid) {
    idx <- which(as.character(catchments$catch_id) == cid)
    if (length(idx) == 0) {
        return(NA)
    }
    round(catchments$residual_area_km2[idx])
}
area_paris <- get_area(catch_paris)
area_cega <- get_area(catch_cega)
area_ticino <- get_area(catch_ticino)

# Helper: compute Sen's slope + MK significance and format annotation
trend_annotation <- function(x, y) {
    valid <- !is.na(x) & !is.na(y)
    xv <- x[valid]
    yv <- y[valid]
    n <- length(yv)
    if (n < 10) {
        return("n < 10")
    }
    # Sen's slope
    slopes <- outer(seq_len(n), seq_len(n), function(i, j) {
        (yv[i] - yv[j]) / (xv[i] - xv[j])
    })
    sen <- median(slopes[upper.tri(slopes)], na.rm = TRUE)
    # Mann-Kendall p-value
    mk <- cor.test(xv, yv, method = "kendall")
    sig <- if (mk$p.value < 0.001) "***" else if (mk$p.value < 0.01) "**" else if (mk$p.value < 0.05) "*" else "n.s."
    sprintf("Sen = %.3g/yr  %s (p=%.3g)", sen, sig, mk$p.value)
}

# --- Panel a: Yearly max daily runoff — Paris ---------------------------------
cat("  [a] Max daily runoff — Paris...\n")
runoff_daily_path <- file.path(agg_dir, "runoff", "runoff_daily_all_years.csv")
ro_daily <- fread(runoff_daily_path)
ro_daily[, date := seq.Date(as.Date("1951-01-01"), by = "day", length.out = .N)]
ro_daily[, year := year(date)]

ro_col_paris <- find_col(catch_paris, setdiff(names(ro_daily), c("period_start", "period_end", "date", "year")))
ro_paris_annual <- ro_daily[, .(max_runoff = max(get(ro_col_paris), na.rm = TRUE)), by = year]

p1a <- ggplot(ro_paris_annual, aes(x = year, y = max_runoff)) +
    geom_line(color = "darkorange3", linewidth = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
    annotate("text",
        x = 2000, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
        label = trend_annotation(ro_paris_annual$year, ro_paris_annual$max_runoff)
    ) +
    labs(
        title = sprintf("a) Yearly max daily runoff \u2014 %s (Paris, %.0f km\u00b2)", catch_paris, area_paris),
        x = NULL, y = "Max daily runoff (mm/day)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

# --- Panel b: Sealed area — Paris --------------------------------------------
cat("  [b] Sealed fraction — Paris...\n")
sealed_files <- list.files(landuse_dir,
    pattern = "^fracsealed_European_01min_\\d{4}\\.nc$", full.names = TRUE
)
sealed_years <- as.integer(sub(".*_(\\d{4})\\.nc$", "\\1", sealed_files))
sealed_order <- order(sealed_years)
sealed_files <- sealed_files[sealed_order]
sealed_years <- sealed_years[sealed_order]

if (!exists("catchments")) catchments <- st_read(gpkg_path, quiet = TRUE)
if (!isTRUE(st_crs(catchments) == st_crs(4326))) {
    catchments <- st_transform(catchments, 4326)
}
catch_sub_paris <- catchments[as.character(catchments$catch_id) == catch_paris, ]

sealed_paris <- data.table(year = sealed_years, frac_sealed = NA_real_)
for (i in seq_along(sealed_files)) {
    r <- tryCatch(terra::rast(sealed_files[i]), error = function(e) NULL)
    if (is.null(r)) next
    if (nlyr(r) > 1) r <- r[[1]]
    sealed_paris$frac_sealed[i] <- exact_extract(r, catch_sub_paris, "mean")
}
sealed_paris <- sealed_paris[!is.na(frac_sealed)]

p1b <- ggplot(sealed_paris, aes(x = year, y = frac_sealed)) +
    geom_line(color = "firebrick", linewidth = 0.8) +
    geom_point(color = "firebrick", size = 1.2, alpha = 0.6) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    annotate("text",
        x = 2000, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
        label = trend_annotation(sealed_paris$year, sealed_paris$frac_sealed)
    ) +
    labs(
        title = sprintf("b) Sealed land fraction \u2014 %s (Paris, %.0f km\u00b2)", catch_paris, area_paris),
        x = NULL, y = "Sealed fraction"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

# --- Panel c: Mean yearly soil moisture — Rio Cega ----------------------------
cat("  [c] Mean yearly SM — Rio Cega...\n")
sm_path_c <- file.path(
    agg_dir, "root_soil_moisture",
    "root_soil_moisture_monthly_all_years.csv"
)
sm_dt_c <- fread(sm_path_c)
sm_dt_c[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
sm_dt_c[, year := year(date)]
sm_dt_c[,month := month(date)]
sm_col_cega <- find_col(catch_cega, setdiff(
    names(sm_dt_c),
    c("month_idx", "period_start", "period_end", "date", "year")
))
sm_cega_annual <- sm_dt_c[month %in% 1:12, .(mean_sm = mean(get(sm_col_cega), na.rm = TRUE)), by = year]
# Remove first year (spin-up)
sm_cega_annual <- sm_cega_annual[year > min(year)]

p1c <- ggplot(sm_cega_annual, aes(x = year, y = mean_sm)) +
    geom_line(color = "sienna", linewidth = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "saddlebrown", linewidth = 0.8) +
    annotate("text",
        x = 2000, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
        label = trend_annotation(sm_cega_annual$year, sm_cega_annual$mean_sm)
    ) +
    labs(
        title = sprintf("c) Mean yearly root SM \u2014 %s (Rio Cega, %.0f km\u00b2)", catch_cega, area_cega),
        x = NULL, y = "Mean SM (mm)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

# --- Panel d: Yearly total AET — Rio Cega ------------------------------------
cat("  [d] Yearly total AET — Rio Cega...\n")
aet_path_c <- file.path(agg_dir, "ActEvapo", "ActEvapo_monthly_all_years.csv")
aet_dt_c <- fread(aet_path_c)
aet_dt_c[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
aet_dt_c[, year := year(date)]
aet_dt_c[,month := month(date)]
aet_col_cega <- find_col(catch_cega, setdiff(
    names(aet_dt_c),
    c("month_idx", "period_start", "period_end", "date", "year")
))
aet_cega_annual <- aet_dt_c[month %in% 1:12, .(total_aet = sum(get(aet_col_cega), na.rm = TRUE)), by = year]

p1d <- ggplot(aet_cega_annual, aes(x = year, y = total_aet)) +
    geom_line(color = "#d62728", linewidth = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.8) +
    annotate("text",
        x = 2000, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
        label = trend_annotation(aet_cega_annual$year, aet_cega_annual$total_aet)
    ) +
    labs(
        title = sprintf("d) Yearly total AET \u2014 %s (Rio Cega, %.0f km\u00b2)", catch_cega, area_cega),
        x = NULL, y = "Total AET (mm/year)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

# --- Panel e: Apr–Oct snowmelt — Ticino ---------------------------------------
cat("  [e] Apr-Oct snowmelt — Ticino...\n")
melt_path_e <- file.path(agg_dir, "snowmelt", "snowmelt_monthly_all_years.csv")
if (!file.exists(melt_path_e)) {
    # Fallback: use TSS
    melt_path_e <- file.path(tss_dir, "snowMeltUpsX_nested_1951_2020.csv")
    melt_dt_e <- fread(melt_path_e, header = TRUE)
    melt_dt_e <- assign_tss_dates(melt_dt_e)
    if ("day_idx" %in% names(melt_dt_e)) {
        melt_cols_e <- setdiff(names(melt_dt_e), c("day_idx", "date"))
        melt_dt_e <- melt_dt_e[, lapply(.SD, sum, na.rm = TRUE), by = date, .SDcols = melt_cols_e]
    }
    melt_dt_e[, year := year(date)]
    melt_dt_e[, month := month(date)]
    melt_col_tic <- find_col(catch_ticino, setdiff(names(melt_dt_e), c("date", "year", "month", "day_idx")))
    melt_ticino <- melt_dt_e[month %in% 5:10, .(snowmelt_sum = sum(get(melt_col_tic), na.rm = TRUE)), by = year]
} else {
    melt_dt_e <- fread(melt_path_e)
    melt_dt_e[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
    melt_dt_e[, year := year(date)]
    melt_dt_e[, month := month(date)]
    melt_col_tic <- find_col(catch_ticino, setdiff(
        names(melt_dt_e),
        c("month_idx", "period_start", "period_end", "date", "year", "month")
    ))
    melt_ticino <- melt_dt_e[month %in% 4:9, .(snowmelt_sum = sum(get(melt_col_tic), na.rm = TRUE)), by = year]
}

# Remove first year (spin-up)
melt_ticino <- melt_ticino[year > min(year)]
p1e <- ggplot(melt_ticino, aes(x = year, y = snowmelt_sum)) +
    geom_line(color = "darkcyan", linewidth = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "darkblue", linewidth = 0.8) +
    annotate("text",
        x = 2000, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
        label = trend_annotation(melt_ticino$year, melt_ticino$snowmelt_sum)
    ) +
    labs(
        title = sprintf("e) Apr\u2013Sept snowmelt \u2014 %s (Ticino, %.0f km\u00b2)", catch_ticino, area_ticino),
        x = NULL, y = "Snowmelt (mm)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

# --- Panel f: Apr–Oct root zone SM — Ticino -----------------------------------
cat("  [f] Apr-Oct baseflow — Ticino...\n")
quz_path_f <- file.path(
    agg_dir, "quz",
    "quz_monthly_all_years.csv"
)
quz_dt_f <- fread(quz_path_f)
quz_dt_f[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
quz_dt_f[, year := year(date)]
quz_dt_f[, month := month(date)]

qlz_path_f <- file.path(
  agg_dir, "qlz",
  "qlz_monthly_all_years.csv"
)
qlz_dt_f <- fread(qlz_path_f)
qlz_dt_f[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
qlz_dt_f[, year := year(date)]
qlz_dt_f[, month := month(date)]

runoff_path_f <- file.path(
  agg_dir, "quz",
  "quz_monthly_all_years.csv"
)
runoff_dt_f <- fread(runoff_path_f)
runoff_dt_f[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
runoff_dt_f[, year := year(date)]
runoff_dt_f[, month := month(date)]

# Identify catchment columns (same across all three)
meta_cols_q <- c("month_idx", "period_start", "period_end", "date", "year", "month")
catch_cols_q <- setdiff(names(quz_dt_f), meta_cols_q)

# Total Q = quz + qlz + runoff (element-wise sum across all catchments)
totalQ_dt_f <- copy(quz_dt_f)
totalQ_dt_f[, (catch_cols_q) := Map(`+`,
                                    Map(`+`, quz_dt_f[, ..catch_cols_q], qlz_dt_f[, ..catch_cols_q]),
                                    runoff_dt_f[, ..catch_cols_q]
)]


rsm_ticino <- totalQ_dt_f[month %in% 6:8, .(mean_rsm = mean(get(rsm_col_tic), na.rm = TRUE)), by = year]
# Remove first year (spin-up)
rsm_ticino <- rsm_ticino[year > min(year)]

p1f <- ggplot(rsm_ticino, aes(x = year, y = mean_rsm)) +
    geom_line(color = "sienna", linewidth = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "saddlebrown", linewidth = 0.8) +
    annotate("text",
        x = 2000, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
        label = trend_annotation(rsm_ticino$year, rsm_ticino$mean_rsm)
    ) +
    labs(
        title = sprintf("f) Jun\u2013Aug Discharge \u2014 %s (Ticino, %.0f km\u00b2)", catch_ticino, area_ticino),
        x = NULL, y = "Q (mm)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

p1f
# --- Compose 6-panel figure ---------------------------------------------------
cat("  Composing 6-panel figure...\n")
fig1_new <- plot_grid(p1a, p1b, p1c, p1d, p1e, p1f,
    ncol = 2, align = "hv"
)

ggsave(file.path(out_dir, "Figure1_new_6panels_vF.png"),
    fig1_new,
    width = 14, height = 12, dpi = 200
)
cat("  -> Saved Figure1_new_6panels.png\n")

# =============================================================================
# FIGURE 2 — CONTINENTAL CLIMATE STRIPES
#   A: Max monthly SWE (area-weighted)
#   B: Min monthly baseflow qlz (area-weighted)
# =============================================================================

# --- Helper: compute area-weighted stripes ------------------------------------
compute_stripes <- function(variable, agg_fun = "max", month_filter = NULL) {
    # Load aggregates
    var_path <- file.path(
        agg_dir, variable,
        paste0(variable, "_monthly_all_years.csv")
    )
    if (!file.exists(var_path)) stop("File not found: ", var_path)

    dt <- fread(var_path)
    dt[, date := as.Date(period_end)]

    # Optional month filter (e.g., only Jan-Jun for SWE)
    if (!is.null(month_filter)) {
        dt <- dt[month(date) %in% month_filter]
    }

    # Identify catchment columns
    meta_cols <- c("month_idx", "period_start", "period_end", "date")
    catch_cols <- setdiff(names(dt), meta_cols)

    # Load catchments for area weighting (reuse if already loaded)
    if (!exists("catchments_for_stripes", envir = .GlobalEnv)) {
        cats <- st_read(gpkg_path, quiet = TRUE)
        assign("catchments_for_stripes", cats, envir = .GlobalEnv)
    }
    cats <- get("catchments_for_stripes", envir = .GlobalEnv)

    common_ids <- intersect(catch_cols, as.character(as.numeric(cats$catch_id)))
    area_vec <- cats$residual_area_km2[match(common_ids, as.character(as.numeric(cats$catch_id)))]
    weights <- area_vec / sum(area_vec, na.rm = TRUE)

    # Area-weighted mean per timestep (handle NAs)
    mat <- as.matrix(dt[, ..common_ids])
    dt[, weighted_mean := apply(mat, 1, function(row) {
        valid <- !is.na(row)
        if (sum(valid) == 0) {
            return(NA_real_)
        }
        sum(row[valid] * weights[valid]) / sum(weights[valid])
    })]

    # Annual aggregation
    if (agg_fun == "max") {
        annual <- dt[, .(annual_val = max(weighted_mean, na.rm = TRUE)), by = .(year = year(date))]
    } else if (agg_fun == "min") {
        annual <- dt[, .(annual_val = min(weighted_mean, na.rm = TRUE)), by = .(year = year(date))]
    } else {
        annual <- dt[, .(annual_val = mean(weighted_mean, na.rm = TRUE)), by = .(year = year(date))]
    }

    # Anomaly
    long_mean <- mean(annual$annual_val, na.rm = TRUE)
    annual[, anomaly := annual_val - long_mean]

    return(annual)
}

# --- Helper: make stripe plot (bar style with color scale) --------------------
make_stripe_plot <- function(annual_dt, title_text, color_lims = NULL,
                             palette = "RdBu", rev_pal = TRUE) {
    if (is.null(color_lims)) {
        mx <- max(abs(annual_dt$anomaly), na.rm = TRUE)
        color_lims <- c(-mx, mx)
    }

    ggplot(annual_dt, aes(x = year, y = anomaly, fill = anomaly)) +
        geom_col(width = 0.8) +
        scale_fill_gradient2(
            low = "#2166AC", mid = "white", high = "#B2182B",
            midpoint = 0, limits = color_lims, oob = squish,
            name = "Anomaly\n(mm)"
        ) +
        scale_x_continuous(expand = c(0.01, 0), breaks = seq(1955, 2020, by = 10)) +
        labs(title = title_text, x = NULL, y = "Anomaly (mm)") +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
            axis.text.x = element_text(size = 9),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.position = "right",
            legend.key.height = unit(1.5, "cm"),
            plot.margin = margin(5, 10, 5, 10)
        )
}

# --- Panel A: Max monthly SWE ------------------------------------------------
cat("[Fig2-A] Stripes: max monthly SWE...\n")
swe_stripes <- compute_stripes("snow_water_equivalent",
    agg_fun = "max",
    month_filter = 1:6
)
pS_A <- make_stripe_plot(swe_stripes,
    title_text = "A) Continental max monthly SWE — annual anomaly",
    palette = "RdBu", rev_pal = FALSE
)

# --- Panel B: Min monthly baseflow (qlz) -------------------------------------
cat("[Fig2-B] Stripes: min monthly qlz...\n")
qlz_stripes <- compute_stripes("qlz", agg_fun = "min")
pS_B <- make_stripe_plot(qlz_stripes,
    title_text = "B) Continental min monthly baseflow (qlz) — annual anomaly",
    palette = "RdBu", rev_pal = FALSE
)

# =============================================================================
# FIGURE 2 — COMPOSE AND SAVE
# =============================================================================
cat("[Fig2] Composing stripes figure...\n")

fig2 <- plot_grid(pS_A, pS_B, ncol = 1, align = "v")

ggsave(file.path(out_dir, "Figure2_continental_stripes.png"),
    fig2,
    width = 14, height = 6, dpi = 200
)
cat("  -> Saved Figure 2 (stripes).\n")

# =============================================================================
# FIGURE 4 — Continental mean AET: LISFLOOD vs GLEAM
# =============================================================================
cat("[Fig4] Continental mean AET: LISFLOOD vs GLEAM...\n")

# LISFLOOD AET: use monthly aggregates, compute area-weighted annual mean
aet_lf_path <- file.path(agg_dir, "ActEvapo", "ActEvapo_monthly_all_years.csv")

aet_lf_dt <- fread(aet_lf_path)
aet_lf_dt[, date := as.Date(period_end)]
aet_lf_dt[, period := as.Date(period_end) - as.Date(period_start)]

meta_cols_aet <- c("month_idx", "period_start", "period_end", "date", "period")
catch_cols_aet <- setdiff(names(aet_lf_dt), meta_cols_aet)


if (!exists("catchments_for_stripes", envir = .GlobalEnv)) {
    cats <- st_read(gpkg_path, quiet = TRUE)
    assign("catchments_for_stripes", cats, envir = .GlobalEnv)
}
cats <- get("catchments_for_stripes", envir = .GlobalEnv)

common_ids_aet <- intersect(catch_cols_aet, as.character(as.numeric(cats$catch_id)))
area_vec_aet <- cats$residual_area_km2[match(common_ids_aet, as.character(as.numeric(cats$catch_id)))]
weights_aet <- area_vec_aet / sum(area_vec_aet, na.rm = TRUE)

mat_aet <- as.matrix(aet_lf_dt[, ..common_ids_aet])
# Convert monthly sums to daily values (divide by number of days in period)
n_days_period <- as.numeric(aet_lf_dt$period)
mat_aet <- mat_aet / n_days_period
# Store daily values back so downstream subsetting also uses daily rates
aet_lf_dt[, (common_ids_aet) := as.data.table(mat_aet)]

aet_lf_dt[, weighted_mean := apply(mat_aet, 1, function(row) {
    valid <- !is.na(row)
    if (sum(valid) == 0) {
        return(NA_real_)
    }
    sum(row[valid] * weights_aet[valid]) / sum(weights_aet[valid])
})]

lf_annual <- aet_lf_dt[, .(aet_mean = mean(weighted_mean, na.rm = TRUE) * 365.25),
    by = .(year = year(date))
]
lf_annual[, source := "LISFLOOD"]

# GLEAM AET: use homogenized monthly file
gleam_homog_path <- file.path(
    base_dir, "output", "aet_diego",
    "1.homogenized", "gleam_monthly_homog.csv"
)
gleam_dt <- fread(gleam_homog_path, header = TRUE)
gleam_dt[, date := as.Date(paste0(date, "-15"), format = "%Y-%m-%d")]

# Same area-weighted mean approach
gleam_catch_cols <- setdiff(names(gleam_dt), "date")
common_ids_gleam <- intersect(gleam_catch_cols, common_ids_aet)
weights_gleam <- area_vec_aet[match(common_ids_gleam, common_ids_aet)]
weights_gleam <- weights_gleam / sum(weights_gleam, na.rm = TRUE)

mat_gleam <- as.matrix(gleam_dt[, ..common_ids_gleam])
gleam_dt[, weighted_mean := apply(mat_gleam, 1, function(row) {
    valid <- !is.na(row)
    if (sum(valid) == 0) {
        return(NA_real_)
    }
    sum(row[valid] * weights_gleam[valid]) / sum(weights_gleam[valid])
})]

gleam_annual <- gleam_dt[, .(aet_mean = mean(weighted_mean, na.rm = TRUE) * 365.25),
    by = .(year = year(date))
]
gleam_annual[, source := "GLEAM v4.3a"]

# Combine and plot
aet_combined <- rbind(lf_annual, gleam_annual)

p_aet <- ggplot(aet_combined, aes(x = year, y = aet_mean, color = source)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 0.9) +
    scale_color_manual(values = c("LISFLOOD" = "steelblue", "GLEAM v4.3a" = "firebrick")) +
    labs(
        title = "Continental mean AET — LISFLOOD vs GLEAM (1951–2020)",
        x = NULL, y = "Annual mean AET (mm/year)",
        color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

ggsave(file.path(out_dir, "Figure4_AET_LISFLOOD_vs_GLEAM.png"),
    p_aet,
    width = 10, height = 5, dpi = 200
)
cat("  -> Saved Figure 4 (AET comparison).\n")

# =============================================================================
# FIGURE 4c — Continental mean soil moisture: ESA CCI vs LISFLOOD
# =============================================================================
cat("[Fig4c] Continental mean SM: ESA CCI vs LISFLOOD...\n")

# LISFLOOD surface soil moisture: monthly aggregates
sm_lf_path <- file.path(
    agg_dir, "surface_soil_moisture",
    "surface_soil_moisture_monthly_all_years.csv"
)
sm_lf_dt <- fread(sm_lf_path)
sm_lf_dt[, date := as.Date(period_start)]

meta_cols_sm <- c("month_idx", "period_start", "period_end", "date")
catch_cols_sm <- setdiff(names(sm_lf_dt), meta_cols_sm)

# Area-weighted mean (reuse catchments already loaded)
common_ids_sm <- intersect(catch_cols_sm, as.character(as.numeric(cats$catch_id)))
area_vec_sm <- cats$residual_area_km2[match(common_ids_sm, as.character(as.numeric(cats$catch_id)))]
weights_sm <- area_vec_sm / sum(area_vec_sm, na.rm = TRUE)

mat_sm <- as.matrix(sm_lf_dt[, ..common_ids_sm])
sm_lf_dt[, weighted_mean := apply(mat_sm, 1, function(row) {
    valid <- !is.na(row)
    if (sum(valid) == 0) {
        return(NA_real_)
    }
    sum(row[valid] * weights_sm[valid]) / sum(weights_sm[valid])
})]

lf_sm_annual <- sm_lf_dt[, .(sm_mean = mean(weighted_mean, na.rm = TRUE)),
    by = .(year = year(date))
]
lf_sm_annual[, source := "LISFLOOD"]

# ESA CCI soil moisture: homogenized monthly
cci_path <- file.path(
    base_dir, "output", "soil_moisture_diego",
    "1.Diego_Merged", "esacci_homogenized.csv"
)
cci_dt <- fread(cci_path, header = TRUE)
cci_dt[, date := as.Date(date)]

# Column names have X prefix — remove for matching
cci_catch_cols <- setdiff(names(cci_dt), "date")
clean_cci_ids <- sub("^X", "", cci_catch_cols)
setnames(cci_dt, cci_catch_cols, clean_cci_ids)

# Match to common catchments with area weights
common_ids_cci <- intersect(clean_cci_ids, common_ids_sm)
weights_cci <- area_vec_sm[match(common_ids_cci, common_ids_sm)]
weights_cci <- weights_cci / sum(weights_cci, na.rm = TRUE)

mat_cci <- as.matrix(cci_dt[, ..common_ids_cci])
cci_dt[, weighted_mean := apply(mat_cci, 1, function(row) {
    valid <- !is.na(row)
    if (sum(valid) == 0) {
        return(NA_real_)
    }
    sum(row[valid] * weights_cci[valid]) / sum(weights_cci[valid])
})]

cci_sm_annual <- cci_dt[, .(sm_mean = mean(weighted_mean, na.rm = TRUE)),
    by = .(year = year(date))
]
cci_sm_annual[, source := "ESA CCI"]

# Combine and plot
sm_combined <- rbind(lf_sm_annual, cci_sm_annual)

p_sm_comp <- ggplot(sm_combined, aes(x = year, y = sm_mean, color = source)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 0.9) +
    scale_color_manual(values = c("LISFLOOD" = "steelblue", "ESA CCI" = "darkorange")) +
    labs(
        title = "Continental mean soil moisture — LISFLOOD vs ESA CCI",
        x = NULL, y = "Annual mean SM (mm)",
        color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

ggsave(file.path(out_dir, "Figure4c_SM_LISFLOOD_vs_ESACCI.png"),
    p_sm_comp,
    width = 10, height = 5, dpi = 200
)
cat("  -> Saved Figure 4c (SM comparison).\n")

# =============================================================================
# FIGURE 4d — Continental mean root zone soil moisture (annual)
# =============================================================================
cat("[Fig4d] Continental mean root zone soil moisture...\n")

rsm_lf_path <- file.path(
    agg_dir, "root_soil_moisture",
    "root_soil_moisture_monthly_all_years.csv"
)
rsm_lf_dt <- fread(rsm_lf_path)
# Use row position for reliable dates (month_idx starts at -334 = Jan 1951)
rsm_lf_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
rsm_lf_dt[, month := month(date)]
rsm_lf_dt[, year := year(date)]

meta_rsm <- c("month_idx", "period_start", "period_end", "date", "month", "year")
catch_cols_rsm <- setdiff(names(rsm_lf_dt), meta_rsm)

# Area-weighted mean (reuse cats and weights from AET section)
common_ids_rsm <- intersect(catch_cols_rsm, as.character(as.numeric(cats$catch_id)))
area_vec_rsm <- cats$residual_area_km2[match(common_ids_rsm, as.character(as.numeric(cats$catch_id)))]
weights_rsm <- area_vec_rsm / sum(area_vec_rsm, na.rm = TRUE)

mat_rsm <- as.matrix(rsm_lf_dt[, ..common_ids_rsm])
rsm_lf_dt[, weighted_mean := apply(mat_rsm, 1, function(row) {
    valid <- !is.na(row)
    if (sum(valid) == 0) {
        return(NA_real_)
    }
    sum(row[valid] * weights_rsm[valid]) / sum(weights_rsm[valid])
})]

rsm_annual <- rsm_lf_dt[, .(rsm_mean = mean(weighted_mean, na.rm = TRUE)),
    by = .(year = year(date))
]

p_rsm <- ggplot(rsm_annual, aes(x = year, y = rsm_mean)) +
    geom_line(color = "sienna", linewidth = 0.7) +
    geom_smooth(
        method = "loess", span = 0.4, se = FALSE,
        color = "saddlebrown", linewidth = 0.9
    ) +
    labs(
        title = "Continental mean root zone soil moisture (1951–2020)",
        x = NULL, y = "Annual mean root zone SM (mm)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

ggsave(file.path(out_dir, "Figure4d_root_SM_continental.png"),
    p_rsm,
    width = 10, height = 5, dpi = 200
)
cat("  -> Saved Figure 4d (root zone SM).\n")


# =============================================================================
# FIGURE 4e — Continental mean direct runoff(annual)
# =============================================================================
cat("[Fig4d] Continental mean direct runoff...\n")

rsm_lf_path <- file.path(
    agg_dir, "runoff",
    "runoff_monthly_all_years.csv"
)
rsm_lf_dt <- fread(rsm_lf_path)
# Use row position for reliable dates (month_idx starts at -334 = Jan 1951)
rsm_lf_dt[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
rsm_lf_dt[, month := month(date)]
rsm_lf_dt[, year := year(date)]

meta_rsm <- c("month_idx", "period_start", "period_end", "date", "month", "year")
catch_cols_rsm <- setdiff(names(rsm_lf_dt), meta_rsm)

# Area-weighted mean (reuse cats and weights from AET section)
common_ids_rsm <- intersect(catch_cols_rsm, as.character(as.numeric(cats$catch_id)))
area_vec_rsm <- cats$residual_area_km2[match(common_ids_rsm, as.character(as.numeric(cats$catch_id)))]
weights_rsm <- area_vec_rsm / sum(area_vec_rsm, na.rm = TRUE)

mat_rsm <- as.matrix(rsm_lf_dt[, ..common_ids_rsm])
rsm_lf_dt[, weighted_mean := apply(mat_rsm, 1, function(row) {
    valid <- !is.na(row)
    if (sum(valid) == 0) {
        return(NA_real_)
    }
    sum(row[valid] * weights_rsm[valid]) / sum(weights_rsm[valid])
})]

rsm_annual <- rsm_lf_dt[, .(rsm_mean = mean(weighted_mean, na.rm = TRUE)),
    by = .(year = year(date))
]

p_rsm <- ggplot(rsm_annual, aes(x = year, y = rsm_mean)) +
    geom_line(color = "sienna", linewidth = 0.7) +
    geom_smooth(
        method = "loess", span = 0.4, se = FALSE,
        color = "saddlebrown", linewidth = 0.9
    ) +
    labs(
        title = "Continental direct runoff (1951–2020)",
        x = NULL, y = "Annual mean direct runoff (mm)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

ggsave(file.path(out_dir, "Figure4d_runoff_continental.png"),
    p_rsm,
    width = 10, height = 5, dpi = 200
)
cat("  -> Saved Figure 4d (root zone SM).\n")

# =============================================================================
# FIGURE 4b — Total yearly AET volume (km³) — LISFLOOD vs GLEAM
# =============================================================================
cat("[Fig4b] Total yearly AET volume...\n")

# Compute total annual AET volume per catchment then sum
# AET is in mm/day (after conversion), area in km²
# Volume = sum_over_catchments( mean_daily_AET_mm * 365.25 * area_km2 * 1e6 ) / 1e9
# = km³/year

# LISFLOOD: sum monthly values (already daily rate * stored back), need annual sum
# Actually the daily rate * days_in_year gives annual mm, then * area -> volume
# Simpler: for each year, annual mean (mm/day) * 365.25 * area_km2 * 1e-6 = km³

# Get area for common catchments (not normalized weights — actual km²)
area_km2_aet <- cats$residual_area_km2[match(common_ids_aet, as.character(as.numeric(cats$catch_id)))]

# LISFLOOD annual volume
lf_vol <- aet_lf_dt[,
    {
        year_mat <- as.matrix(.SD)
        annual_mean <- colMeans(year_mat, na.rm = TRUE) # mm/day per catchment
        # volume in km³: mm/day * 365.25 days * km² * 1e6 m²/km² * 1e-3 m/mm * 1e-9 km³/m³
        vol_km3 <- sum(annual_mean * 365.25 * area_km2_aet * 1e-6, na.rm = TRUE)
        .(aet_volume_km3 = vol_km3)
    },
    by = .(year = year(date)),
    .SDcols = common_ids_aet
]
lf_vol[, source := "LISFLOOD"]

# GLEAM annual volume (use same area vector for matching catchments)
area_km2_gleam <- area_km2_aet[match(common_ids_gleam, common_ids_aet)]

gleam_vol <- gleam_dt[,
    {
        year_mat <- as.matrix(.SD)
        annual_mean <- colMeans(year_mat, na.rm = TRUE)
        vol_km3 <- sum(annual_mean * 365.25 * area_km2_gleam * 1e-6, na.rm = TRUE)
        .(aet_volume_km3 = vol_km3)
    },
    by = .(year = year(date)),
    .SDcols = common_ids_gleam
]
gleam_vol[, source := "GLEAM v4.3a"]

vol_combined <- rbind(lf_vol, gleam_vol)

p_vol <- ggplot(vol_combined, aes(x = year, y = aet_volume_km3, color = source)) +
    geom_line(linewidth = 0.7) +
    geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 0.9) +
    scale_color_manual(values = c("LISFLOOD" = "steelblue", "GLEAM v4.3a" = "firebrick")) +
    labs(
        title = "Total yearly AET volume — LISFLOOD vs GLEAM",
        x = NULL, y = expression("AET volume (km"^3 * "/year)"),
        color = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom"
    )

ggsave(file.path(out_dir, "Figure4b_AET_volume_LISFLOOD_vs_GLEAM.png"),
    p_vol,
    width = 10, height = 5, dpi = 200
)
cat("  -> Saved Figure 4b (AET volume comparison).\n")

# =============================================================================
# FIGURE 5 — Mean AET ratio bias map: LISFLOOD / GLEAM per catchment
# =============================================================================
cat("[Fig5] AET ratio bias map (mean LISFLOOD / mean GLEAM)...\n")

# Compute long-term mean AET per catchment for overlapping period (1980–2020)
# LISFLOOD (already in mm/day from earlier conversion)
lf_overlap <- aet_lf_dt[year(date) >= 1980]
lf_mean_catch <- colMeans(as.matrix(lf_overlap[, ..common_ids_aet]), na.rm = TRUE)

# GLEAM
gleam_mean_catch <- colMeans(as.matrix(gleam_dt[, ..common_ids_gleam]), na.rm = TRUE)

# Compute ratio for common catchments
ratio_ids <- intersect(names(lf_mean_catch), names(gleam_mean_catch))
ratio_vals <- lf_mean_catch[ratio_ids] / gleam_mean_catch[ratio_ids]

# Handle Inf/NaN from division by zero
ratio_vals[!is.finite(ratio_vals)] <- NA

ratio_df <- data.frame(
    catch_id = ratio_ids,
    ratio = as.numeric(ratio_vals),
    stringsAsFactors = FALSE
)

# Merge with catchment geometries
cats_ratio <- cats[as.character(as.numeric(cats$catch_id)) %in% ratio_ids, ]
cats_ratio$ratio <- ratio_df$ratio[match(
    as.character(as.numeric(cats_ratio$catch_id)), ratio_df$catch_id
)]

cats_ratio <- st_transform(cats_ratio, 3035)

if (!exists("basemap")) {
    basemap <- ne_countries(scale = "medium", returnclass = "sf") |>
        st_transform(3035)
}
bbox_r <- st_bbox(cats_ratio)

# Color scale: ratio centered on 1 (no bias)
palet_ratio <- hcl.colors(11, palette = "RdBu", rev = TRUE)

p_ratio <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_ratio, aes(fill = ratio), color = NA) +
    scale_fill_gradientn(
        colors = palet_ratio,
        limits = c(0.5, 1.5),
        oob = squish,
        name = "LISFLOOD /\nGLEAM"
    ) +
    coord_sf(
        xlim = c(bbox_r["xmin"], bbox_r["xmax"]),
        ylim = c(bbox_r["ymin"], bbox_r["ymax"]), expand = FALSE
    ) +
    labs(
        title = "Mean AET ratio: LISFLOOD / GLEAM (1980–2020)",
        subtitle = "Values > 1: LISFLOOD overestimates. Values < 1: underestimates."
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
    )

ggsave(file.path(out_dir, "Figure5_AET_ratio_bias_map.png"),
    p_ratio,
    width = 10, height = 8, dpi = 200
)
cat("  -> Saved Figure 5 (AET ratio bias map).\n")

# =============================================================================
# FIGURE S1 — Spatial AET bias map: LISFLOOD - GLEAM for a given year
# =============================================================================
bias_year <- 2007 # Change this to plot a different year

cat("[FigS1] AET bias map (LISFLOOD - GLEAM) for year", bias_year, "...\n")

# Compute annual mean per catchment for the chosen year — LISFLOOD
lf_year <- aet_lf_dt[year(date) == bias_year, ..common_ids_aet]
lf_annual_catch <- colMeans(lf_year, na.rm = TRUE)

# Compute annual mean per catchment for the chosen year — GLEAM
gleam_year <- gleam_dt[year(date) == bias_year, ..common_ids_gleam]
gleam_annual_catch <- colMeans(gleam_year, na.rm = TRUE)

# Compute bias for catchments present in both
bias_ids <- intersect(names(lf_annual_catch), names(gleam_annual_catch))
bias_vals <- lf_annual_catch[bias_ids] - gleam_annual_catch[bias_ids]

bias_df <- data.frame(
    catch_id = bias_ids,
    bias = as.numeric(bias_vals),
    stringsAsFactors = FALSE
)

# Merge with catchment geometries for spatial plot
library(rnaturalearth)

cats_plot <- cats[as.character(as.numeric(cats$catch_id)) %in% bias_ids, ]
cats_plot$bias <- bias_df$bias[match(
    as.character(as.numeric(cats_plot$catch_id)), bias_df$catch_id
)]

# Transform to EPSG:3035 for European map
cats_plot <- st_transform(cats_plot, 3035)
basemap <- ne_countries(scale = "medium", returnclass = "sf") |>
    st_transform(3035)

# Bounding box from catchments
bbox <- st_bbox(cats_plot)

# Symmetric color limits
bias_lim <- max(abs(bias_df$bias), na.rm = TRUE) * 0.8
palet_bias <- hcl.colors(11, palette = "RdBu", rev = TRUE)

p_bias <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_plot, aes(fill = bias), color = NA) +
    scale_fill_gradientn(
        colors = palet_bias,
        limits = c(-bias_lim, bias_lim),
        oob = squish,
        name = "Bias\n(mm/day)"
    ) +
    coord_sf(
        xlim = c(bbox["xmin"], bbox["xmax"]),
        ylim = c(bbox["ymin"], bbox["ymax"]),
        expand = FALSE
    ) +
    labs(
        title = paste0("AET bias (LISFLOOD - GLEAM) — Year ", bias_year),
        subtitle = "Annual mean per catchment (mm/day)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
    )

ggsave(file.path(out_dir, paste0("FigureS1_AET_bias_map_", bias_year, ".png")),
    p_bias,
    width = 10, height = 8, dpi = 200
)
cat("  -> Saved Figure S1 (AET bias map,", bias_year, ").\n")

# =============================================================================
# FIGURE S2 — GLEAM step-change diagnostic
#   Maps the per-catchment difference in mean GLEAM AET between two periods
#   to identify which catchments drive the ~1998-2000 upward jump.
#   Also shows the same difference for LISFLOOD as a control panel.
# =============================================================================
cat("[FigS2] GLEAM step-change diagnostic...\n")

# --- User choice: periods to compare -----------------------------------------
period_before <- 1990:1999
period_after <- 2000:2009

# --- Per-catchment mean AET in each period (GLEAM) ----------------------------
gleam_before <- gleam_dt[year(date) %in% period_before, ..common_ids_gleam]
gleam_after <- gleam_dt[year(date) %in% period_after, ..common_ids_gleam]

gleam_mean_before <- colMeans(gleam_before, na.rm = TRUE)
gleam_mean_after <- colMeans(gleam_after, na.rm = TRUE)
gleam_step <- gleam_mean_after - gleam_mean_before

# --- Per-catchment mean AET in each period (LISFLOOD) -------------------------
lf_before <- aet_lf_dt[year(date) %in% period_before, ..common_ids_aet]
lf_after <- aet_lf_dt[year(date) %in% period_after, ..common_ids_aet]

lf_mean_before <- colMeans(lf_before, na.rm = TRUE)
lf_mean_after <- colMeans(lf_after, na.rm = TRUE)
lf_step <- lf_mean_after - lf_mean_before

# --- Build spatial data for common catchments ---------------------------------
step_ids <- intersect(names(gleam_step), names(lf_step))

step_df <- data.frame(
    catch_id = step_ids,
    gleam_step = as.numeric(gleam_step[step_ids]),
    lf_step = as.numeric(lf_step[step_ids]),
    stringsAsFactors = FALSE
)
# Difference of differences: how much more did GLEAM jump vs LISFLOOD?
step_df$excess_gleam <- step_df$gleam_step - step_df$lf_step

# Merge with geometry
cats_step <- cats[as.character(as.numeric(cats$catch_id)) %in% step_ids, ]
cats_step <- merge(cats_step, step_df,
    by.x = "catch_id", by.y = "catch_id", all.x = FALSE
)
cats_step <- st_transform(cats_step, 3035)

if (!exists("basemap")) {
    basemap <- ne_countries(scale = "medium", returnclass = "sf") |>
        st_transform(3035)
}
bbox_s <- st_bbox(cats_step)

# --- Panel A: GLEAM step-change map -------------------------------------------
step_lim <- max(abs(step_df$gleam_step), na.rm = TRUE) * 0.8
palet_step <- hcl.colors(11, palette = "RdBu", rev = TRUE)

pStep_A <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_step, aes(fill = gleam_step), color = NA) +
    scale_fill_gradientn(
        colors = palet_step,
        limits = c(-step_lim, step_lim),
        oob = squish,
        name = "\u0394 AET\n(mm/day)"
    ) +
    coord_sf(
        xlim = c(bbox_s["xmin"], bbox_s["xmax"]),
        ylim = c(bbox_s["ymin"], bbox_s["ymax"]), expand = FALSE
    ) +
    labs(
        title = paste0(
            "A) GLEAM AET change: ",
            min(period_after), "–", max(period_after),
            " vs ", min(period_before), "–", max(period_before)
        ),
        subtitle = "Mean AET difference per catchment"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "right"
    )

# --- Panel B: LISFLOOD step-change map (control) ------------------------------
lf_step_lim <- max(abs(step_df$lf_step), na.rm = TRUE) * 0.8

pStep_B <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_step, aes(fill = lf_step), color = NA) +
    scale_fill_gradientn(
        colors = palet_step,
        limits = c(-lf_step_lim, lf_step_lim),
        oob = squish,
        name = "\u0394 AET\n(mm/day)"
    ) +
    coord_sf(
        xlim = c(bbox_s["xmin"], bbox_s["xmax"]),
        ylim = c(bbox_s["ymin"], bbox_s["ymax"]), expand = FALSE
    ) +
    labs(
        title = paste0(
            "B) LISFLOOD AET change: ",
            min(period_after), "–", max(period_after),
            " vs ", min(period_before), "–", max(period_before)
        ),
        subtitle = "Mean AET difference per catchment (control)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "right"
    )

# --- Panel C: Excess GLEAM jump (GLEAM step - LISFLOOD step) ------------------
excess_lim <- max(abs(step_df$excess_gleam), na.rm = TRUE) * 0.8

pStep_C <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(data = cats_step, aes(fill = excess_gleam), color = NA) +
    scale_fill_gradientn(
        colors = palet_step,
        limits = c(-excess_lim, excess_lim),
        oob = squish,
        name = "Excess\n(mm/day)"
    ) +
    coord_sf(
        xlim = c(bbox_s["xmin"], bbox_s["xmax"]),
        ylim = c(bbox_s["ymin"], bbox_s["ymax"]), expand = FALSE
    ) +
    labs(
        title = "C) Excess GLEAM jump (GLEAM \u0394 - LISFLOOD \u0394)",
        subtitle = "Positive = GLEAM jumped more than LISFLOOD"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "right"
    )

# --- Panel D: Histogram of excess jump per catchment --------------------------
pStep_D <- ggplot(step_df, aes(x = excess_gleam)) +
    geom_histogram(bins = 40, fill = "firebrick", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    labs(
        title = "D) Distribution of excess GLEAM step-change",
        x = "GLEAM \u0394 - LISFLOOD \u0394 (mm/day)",
        y = "Number of catchments"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

# --- Compose and save ---------------------------------------------------------
fig_s2 <- plot_grid(pStep_A, pStep_B, pStep_C, pStep_D,
    ncol = 2, align = "hv"
)

ggsave(file.path(out_dir, "FigureS2_GLEAM_stepchange_diagnostic.png"),
    fig_s2,
    width = 16, height = 12, dpi = 200
)
cat("  -> Saved Figure S2 (GLEAM step-change diagnostic).\n")

# --- Print summary stats ------------------------------------------------------
cat("\n--- GLEAM step-change summary ---\n")
cat("  Period before:", min(period_before), "-", max(period_before), "\n")
cat("  Period after: ", min(period_after), "-", max(period_after), "\n")
cat("  Mean GLEAM step (all catchments):", round(mean(step_df$gleam_step, na.rm = TRUE), 2), "mm/day\n")
cat("  Mean LISFLOOD step:              ", round(mean(step_df$lf_step, na.rm = TRUE), 2), "mm/day\n")
cat("  Mean excess GLEAM jump:          ", round(mean(step_df$excess_gleam, na.rm = TRUE), 2), "mm/day\n")
cat(
    "  % catchments with excess > 2 mm/day:",
    round(100 * mean(step_df$excess_gleam > 2, na.rm = TRUE), 1), "%\n"
)
cat("  Top 10 catchments by excess GLEAM jump:\n")
top10 <- step_df[order(-step_df$excess_gleam), ][1:10, ]
print(top10[, c("catch_id", "gleam_step", "lf_step", "excess_gleam")])

cat("\nDone! All figures saved to:", out_dir, "\n")

# =============================================================================
# FIGURE S3 — Map of trends in winter surface soil moisture (all catchments)
# =============================================================================
cat("[FigS3] Computing winter surface SM trends per catchment...\n")

# Load monthly surface soil moisture
sm_trend_path <- file.path(
    agg_dir, "surface_soil_moisture",
    "surface_soil_moisture_monthly_all_years.csv"
)
sm_trend_dt <- fread(sm_trend_path)
sm_trend_dt[, date := as.Date(period_start)]

# Winter months: Dec, Jan, Feb, Mar
sm_trend_dt <- sm_trend_dt[month(date) %in% c(12, 1, 2, 3)]

# Identify catchment columns
meta_sm <- c("month_idx", "period_start", "period_end", "date")
catch_cols_sm <- setdiff(names(sm_trend_dt), meta_sm)

# Compute winter mean per year (assign Dec to the following winter year)
sm_trend_dt[, winter_year := fifelse(month(date) == 12, year(date) + 1L, year(date))]

# Remove incomplete first/last winter
sm_trend_dt <- sm_trend_dt[winter_year > min(winter_year) & winter_year < max(winter_year)]

# Compute Sen's slope (linear trend) for each catchment
# Using simple linear regression (year as predictor) for speed
cat("  Computing linear trends for", length(catch_cols_sm), "catchments...\n")

trend_results <- data.table(
    catch_id = catch_cols_sm,
    slope = NA_real_,
    p_value = NA_real_
)

for (i in seq_along(catch_cols_sm)) {
    cid <- catch_cols_sm[i]
    winter_annual <- sm_trend_dt[, .(
        winter_sm = mean(get(cid), na.rm = TRUE)
    ), by = .(year = winter_year)]

    # Skip if too few valid years
    valid <- winter_annual[!is.na(winter_sm)]
    if (nrow(valid) < 10) next

    fit <- lm(winter_sm ~ year, data = valid)
    trend_results[i, slope := coef(fit)[2]]
    trend_results[i, p_value := summary(fit)$coefficients[2, 4]]
}

# Remove catchments with no trend computed
trend_results <- trend_results[!is.na(slope)]

cat("  Trends computed for", nrow(trend_results), "catchments.\n")
cat(
    "  Significant (p<0.05) positive trends:",
    sum(trend_results$slope > 0 & trend_results$p_value < 0.05), "\n"
)
cat(
    "  Significant (p<0.05) negative trends:",
    sum(trend_results$slope < 0 & trend_results$p_value < 0.05), "\n"
)

# Merge with catchment geometries
if (!exists("cats")) {
    cats <- st_read(gpkg_path, quiet = TRUE)
}
cats_trend <- cats[as.character(as.numeric(cats$catch_id)) %in% trend_results$catch_id, ]
cats_trend <- merge(cats_trend, trend_results,
    by.x = "catch_id", by.y = "catch_id", all.x = FALSE
)
cats_trend <- st_transform(cats_trend, 3035)

if (!exists("basemap")) {
    basemap <- ne_countries(scale = "medium", returnclass = "sf") |>
        st_transform(3035)
}
bbox_t <- st_bbox(cats_trend)

# Symmetric color limits for slope
slope_lim <- quantile(abs(trend_results$slope), 0.95, na.rm = TRUE)
palet_trend <- hcl.colors(11, palette = "BrBG", rev = FALSE)

# Mark non-significant trends with transparency
cats_trend$significant <- cats_trend$p_value < 0.05

p_sm_trend <- ggplot() +
    geom_sf(data = basemap, fill = "grey95", color = "grey60", linewidth = 0.2) +
    geom_sf(
        data = cats_trend[!cats_trend$significant, ],
        aes(fill = slope), color = NA, alpha = 0.3
    ) +
    geom_sf(
        data = cats_trend[cats_trend$significant, ],
        aes(fill = slope), color = NA, alpha = 0.9
    ) +
    scale_fill_gradientn(
        colors = palet_trend,
        limits = c(-slope_lim, slope_lim),
        oob = squish,
        name = "Trend\n(mm/year²)"
    ) +
    coord_sf(
        xlim = c(bbox_t["xmin"], bbox_t["xmax"]),
        ylim = c(bbox_t["ymin"], bbox_t["ymax"]), expand = FALSE
    ) +
    labs(
        title = "Trends in winter surface soil moisture (DJFM)",
        subtitle = "Linear slope per catchment (1951–2020). Faded = non-significant (p=0.05)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10),
        legend.position = "right"
    )

ggsave(file.path(out_dir, "FigureS3_winter_SM_trends_map.png"),
    p_sm_trend,
    width = 10, height = 8, dpi = 200
)
cat("  -> Saved Figure S3 (winter SM trends map).\n")

cat("\nDone! All figures saved to:", out_dir, "\n")
