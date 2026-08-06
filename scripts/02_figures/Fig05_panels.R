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
catch_roer <- "189926"
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
area_roer <- get_area(catch_roer)
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

ro_col_roer <- find_col(catch_roer, setdiff(names(ro_daily), c("period_start", "period_end", "date", "year")))
ro_roer_annual <- ro_daily[, .(max_runoff = max(get(ro_col_roer), na.rm = TRUE)), by = year]

p1a <- ggplot(ro_roer_annual, aes(x = year, y = max_runoff)) +
    geom_line(color = "darkgrey", linewidth = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "grey20", linewidth = 0.8) +
    annotate("text",
        x = 2000, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
        label = trend_annotation(ro_roer_annual$year, ro_roer_annual$max_runoff)
    ) +
    labs(
        title = sprintf("a) Yearly max daily runoff \u2014 %s (Roer, %.0f km\u00b2)", catch_roer, area_roer),
        x = NULL, y = "Max daily runoff (mm/day)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))
p1a


# --- Panel a: Yearly max daily runoff — Paris ---------------------------------
cat("  [a] Max daily discharge — roer...\n")
quz_daily_path <- file.path(agg_dir, "quz", "quz_daily_all_years.csv")
qlz_daily_path <- file.path(agg_dir, "qlz", "qlz_daily_all_years.csv")
runoff_daily_path <- file.path(agg_dir, "runoff", "runoff_daily_all_years.csv")
dis_daily_path <- file.path(agg_dir, "discharge", "discharge_daily_all_years.csv")
ssm_daily_path <- file.path(agg_dir, "surface_soil_moisture", "surface_soil_moisture_daily_all_years.csv")

ro_daily <- fread(runoff_daily_path)
ro_daily[, date := seq.Date(as.Date("1951-01-01"), by = "day", length.out = .N)]
ro_daily[, year := year(date)]

quz <- fread(quz_daily_path)
quz[, date := seq.Date(as.Date("1951-01-01"), by = "day", length.out = .N)]
quz[, year := year(date)]

qlz <- fread(qlz_daily_path)
qlz[, date := seq.Date(as.Date("1951-01-01"), by = "day", length.out = .N)]
qlz[, year := year(date)]

dis_daily <- fread(dis_daily_path)
dis_daily[, date := seq.Date(as.Date("1951-01-01"), by = "day", length.out = .N)]
dis_daily[, year := year(date)]


sm_daily <- fread(ssm_daily_path)
sm_daily[, date := seq.Date(as.Date("1951-01-01"), by = "day", length.out = .N)]
sm_daily[, year := year(date)]

# Identify catchment columns (same across all three)
meta_cols_q <- c("month_idx", "period_start", "period_end", "doy", "date", "year")
catch_cols_q <- setdiff(names(quz), meta_cols_q)
totalQ_dt <- copy(quz)
totalQ_dt[, (catch_cols_q) := Map(
    `+`,
    Map(`+`, quz[, ..catch_cols_q], qlz[, ..catch_cols_q]),
    ro_daily[, ..catch_cols_q]
)]

ro_col_roer <- find_col(catch_roer, setdiff(names(ro_daily), c("period_start", "period_end", "date", "year", "doy")))

q_roer_annual <- totalQ_dt[, .(max_runoff = max(get(ro_col_roer), na.rm = TRUE)), by = year]

dis_roer_annual <- dis_daily[, .(max_runoff = max(get(ro_col_roer), na.rm = TRUE)), by = year]
dis_roer_annual$qmm <- dis_roer_annual$max_runoff * 86400 / (area_roer * 1000)

p1b <- ggplot(dis_roer_annual, aes(x = year, y = qmm)) +
    geom_line(color = "royalblue", linewidth = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "royalblue4", linewidth = 0.8) +
    annotate("text",
        x = 2000, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
        label = trend_annotation(dis_roer_annual$year, dis_roer_annual$qmm)
    ) +
    labs(
        title = sprintf("a) Yearly max daily discharge \u2014 %s (Roer, %.0f km\u00b2)", catch_roer, area_roer),
        x = NULL, y = "Max daily discharge (mm/day)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

p1b
# --- Panel supp: Sealed area — roer --------------------------------------------
cat("  [b] Sealed fraction — roer...\n")
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
catch_sub_roer <- catchments[as.character(catchments$catch_id) == catch_roer, ]

sealed_roer <- data.table(year = sealed_years, frac_sealed = NA_real_)
for (i in seq_along(sealed_files)) {
    r <- tryCatch(terra::rast(sealed_files[i]), error = function(e) NULL)
    if (is.null(r)) next
    if (nlyr(r) > 1) r <- r[[1]]
    sealed_roer$frac_sealed[i] <- exact_extract(r, catch_sub_roer, "mean")
}
sealed_roer <- sealed_roer[!is.na(frac_sealed)]

p1s <- ggplot(sealed_roer, aes(x = year, y = frac_sealed)) +
    geom_line(color = "firebrick", linewidth = 0.8) +
    geom_point(color = "firebrick", size = 1.2, alpha = 0.6) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    annotate("text",
        x = 2000, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
        label = trend_annotation(sealed_roer$year, sealed_roer$frac_sealed)
    ) +
    labs(
        title = sprintf("b) Sealed land fraction \u2014 %s (roer, %.0f km\u00b2)", catch_roer, area_roer),
        x = NULL, y = "Sealed fraction"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))
p1s
# --- Panel c: Mean yearly soil moisture — Rio Cega ----------------------------
cat("  [c] Mean yearly SM — Rio Cega...\n")
sm_path_c <- file.path(
    agg_dir, "root_soil_moisture",
    "root_soil_moisture_monthly_all_years.csv"
)
sm_dt_c <- fread(sm_path_c)
sm_dt_c[, date := seq.Date(as.Date("1951-01-01"), by = "month", length.out = .N)]
sm_dt_c[, year := year(date)]
sm_dt_c[, month := month(date)]
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
aet_dt_c[, month := month(date)]
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
    geom_line(color = "lightblue3", linewidth = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "cyan4", linewidth = 0.8) +
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

p1e
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
totalQ_dt_f[, (catch_cols_q) := Map(
    `+`,
    Map(`+`, quz_dt_f[, ..catch_cols_q], qlz_dt_f[, ..catch_cols_q]),
    runoff_dt_f[, ..catch_cols_q]
)]

q_col_tic <- find_col(catch_ticino, names(dis_daily))

dis_daily[, month := month(date)]
dis_ticino_annual <- dis_daily[month %in% 6:8, .(mean_q = mean(get(q_col_tic), na.rm = TRUE)), by = year]

# m3/s to mm/month
dis_ticino_annual$qmm <- dis_ticino_annual$mean_q * 86400 / (area_ticino * 1000) * 30.7
dis_ticino_annual <- dis_ticino_annual[year > min(year)]

ptn <- totalQ_dt_f[[q_col_tic]]
plot(ptn)
rsm_ticino <- totalQ_dt_f[month %in% 6:8, .(mean_q = mean(get(q_col_tic), na.rm = TRUE)), by = year]
# Remove first year (spin-up)
rsm_ticino <- rsm_ticino[year > min(year)]

plot(rsm_ticino$mean_q, dis_ticino_annual$qmm)
p1f <- ggplot(dis_ticino_annual, aes(x = year, y = qmm)) +
    geom_line(color = "royalblue", linewidth = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "royalblue4", linewidth = 0.8) +
    annotate("text",
        x = 2000, y = Inf, vjust = 1.5, hjust = 0, size = 3.2,
        label = trend_annotation(dis_ticino_annual$year, dis_ticino_annual$qmm)
    ) +
    labs(
        title = sprintf("f) Jun\u2013Aug Discharge \u2014 %s (Ticino, %.0f km\u00b2)", catch_ticino, area_ticino),
        x = NULL, y = "Q (mm/month)"
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
    width = 14, height = 12, dpi = 400
)
cat("  -> Saved Figure5.png\n")











# =============================================================================
# Supplementary Figure — CONTINENTAL CLIMATE STRIPES
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
    } else if (agg_fun == "sum") {
        annual <- dt[, .(annual_val = sum(weighted_mean, na.rm = TRUE)), by = .(year = year(date))]
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
                             palette = "RdBu", rev_pal = F) {
    if (is.null(color_lims)) {
        mx <- max(abs(annual_dt$anomaly), na.rm = TRUE)
        color_lims <- c(-mx, mx)
    }

    pal <- hcl.colors(11, palette = palette, rev = rev_pal)
    ggplot(annual_dt, aes(x = year, y = anomaly, fill = anomaly)) +
        geom_col(width = 0.8) +
        scale_fill_gradientn(
            colors = pal,
            limits = color_lims, oob = squish,
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
    palette = "RdBu", rev_pal = F
)

pS_A
# --- Panel B: Min monthly baseflow (qlz) -------------------------------------
cat("[Fig2-B] Stripes: total yearly AET...\n")
qlz_stripes <- compute_stripes("ActEvapo", agg_fun = "sum")
pS_B <- make_stripe_plot(qlz_stripes,
    title_text = "B) Continental yearly AET — annual anomaly",
    palette = "RdBu", rev_pal = T
)
pS_B
# =============================================================================
# FIGURE SUPPLEMENT — COMPOSE AND SAVE
# =============================================================================
cat("[Fig2] Composing stripes figure...\n")

fig2 <- plot_grid(pS_A, pS_B, ncol = 1, align = "v")

ggsave(file.path(out_dir, "Figure2_continental_stripes.png"),
    fig2,
    width = 14, height = 14, dpi = 400
)
cat("  -> Saved Figure 2 (stripes).\n")


# =============================================================================
# CLASSIC CLIMATE STRIPES — All variables (continental area-weighted average)
# =============================================================================
# Produces Ed Hawkins-style warming stripes: one colored vertical bar per year,
# no axes, no gaps. Color encodes the anomaly relative to the long-term mean.
# Each variable gets its own stripe panel, then they are assembled into a
# multi-panel figure.
# =============================================================================

cat("[Stripes] Building classic climate stripes for all variables...\n")

# --- Classic stripe plot function ---------------------------------------------
#' Create a classic Ed Hawkins-style stripe plot
#'
#' @param annual_dt   data.table with columns: year, annual_val, anomaly
#' @param title_text  Character string for panel title
#' @param palette     HCL palette name (default "RdBu")
#' @param rev_pal     Logical, reverse palette direction
#' @return ggplot object
make_classic_stripe <- function(annual_dt, title_text, palette = "RdBu",
                                rev_pal = FALSE) {
    # Normalize anomaly to [-1, 1] range so all panels share the same color scale
    mx <- max(abs(annual_dt$anomaly), na.rm = TRUE)
    if (mx == 0) mx <- 1 # avoid division by zero
    annual_dt[, norm_anomaly := anomaly / mx]

    pal <- hcl.colors(11, palette = palette, rev = rev_pal)

    ggplot(annual_dt, aes(x = year, y = 1, fill = norm_anomaly)) +
        geom_tile(width = 1, height = 1) +
        scale_fill_gradientn(
            colors = pal,
            limits = c(-1, 1), oob = squish,
            name = NULL
        ) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        labs(title = title_text) +
        theme_void(base_size = 10) +
        theme(
            plot.title = element_text(
                face = "bold", hjust = 0, size = 10,
                margin = margin(b = 2)
            ),
            legend.position = "none",
            plot.margin = margin(2, 5, 2, 5)
        )
}

# --- Define all variables to stripe -------------------------------------------
# Each entry: list(variable folder name, aggregation function, month filter,
#                  display name, palette, reverse palette)
stripe_vars <- list(
    list(
        var = "rainfall", agg = "mean", months = NULL,
        label = "Rainfall (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "snowfall", agg = "mean", months = NULL,
        label = "Snowfall (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "snowmelt", agg = "mean",
        label = "Snowmelt (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "snow_water_equivalent", agg = "mean",
        label = "Max SWE (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "ActEvapo", agg = "mean", months = NULL,
        label = "Actual ET (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "infiltration", agg = "mean", months = NULL,
        label = "Infiltration (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "runoff", agg = "mean", months = NULL,
        label = "Surface runoff (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "quz", agg = "mean", months = NULL,
        label = "Interflow (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "qlz", agg = "mean", months = NULL,
        label = "Baseflow (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "percolation", agg = "sum", months = NULL,
        label = "Percolation (total)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "root_soil_moisture", agg = "mean",
        label = "Root zone SM (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "surface_soil_moisture", agg = "mean",
        label = "Surface SM (mean)", palette = "RdBu", rev = TRUE
    ),
    list(
        var = "discharge", agg = "mean", months = NULL,
        label = "Discharge (mean)", palette = "RdBu", rev = TRUE
    )
)

# --- Compute stripes for each variable and build plot list --------------------
stripe_plots <- list()

for (sv in stripe_vars) {
    var_path <- file.path(
        agg_dir, sv$var,
        paste0(sv$var, "_monthly_all_years.csv")
    )
    if (!file.exists(var_path)) {
        message("  Skipping (not found): ", sv$var)
        next
    }

    cat(sprintf("  Computing stripes: %s ...\n", sv$label))

    tryCatch(
        {
            annual_dt <- compute_stripes(
                variable = sv$var,
                agg_fun = sv$agg,
                month_filter = sv$months
            )

            p <- make_classic_stripe(
                annual_dt,
                title_text = sv$label,
                palette = sv$palette,
                rev_pal = sv$rev
            )
            stripe_plots[[sv$var]] <- p
        },
        error = function(e) {
            message("  ERROR for ", sv$var, ": ", conditionMessage(e))
        }
    )
}

# --- Add a shared x-axis label strip at the bottom ----------------------------
# Create a minimal plot that just shows the year axis
year_axis <- ggplot(data.frame(x = 1951:2020, y = 1), aes(x = x, y = y)) +
    geom_blank() +
    scale_x_continuous(
        expand = c(0, 0),
        breaks = seq(1955, 2020, by = 5),
        limits = c(1950.5, 2020.5)
    ) +
    theme_void(base_size = 10) +
    theme(
        axis.text.x = element_text(size = 8, margin = margin(t = 2)),
        axis.ticks.x = element_line(linewidth = 0.3),
        axis.ticks.length.x = unit(2, "pt"),
        plot.margin = margin(0, 5, 5, 5)
    )

# --- Compose multi-panel stripes figure ---------------------------------------
n_panels <- length(stripe_plots)
if (n_panels > 0) {
    # Combine all stripe panels vertically with equal heights
    # Add year axis at the bottom with smaller relative height
    all_panels <- c(stripe_plots, list(year_axis = year_axis))
    rel_heights <- c(rep(1, n_panels), 0.4)

    stripes_body <- plot_grid(
        plotlist = all_panels,
        ncol = 1, align = "v",
        rel_heights = rel_heights
    )

    # Create a shared legend (colorbar) on the right side
    # Labels indicate positive/negative anomaly direction
    legend_dt <- data.table(
        x = 1, y = seq(-1, 1, length.out = 100),
        fill = seq(-1, 1, length.out = 100)
    )
    pal_legend <- hcl.colors(11, palette = "RdBu", rev = FALSE)

    p_legend <- ggplot(legend_dt, aes(x = x, y = y, fill = fill)) +
        geom_tile(width = 1, height = 0.02) +
        scale_fill_gradientn(colors = pal_legend, limits = c(-1, 1)) +
        annotate("text",
            x = 1, y = 1.12, label = "+ anomaly",
            size = 3.2, fontface = "bold", hjust = 0.5
        ) +
        annotate("text",
            x = 1, y = -1.12, label = "\u2013 anomaly",
            size = 3.2, fontface = "bold", hjust = 0.5
        ) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0.15, 0)) +
        coord_cartesian(clip = "off") +
        theme_void() +
        theme(
            legend.position = "none",
            plot.margin = margin(10, 10, 10, 10)
        )

    # Combine body + legend side by side
    fig_stripes <- plot_grid(
        stripes_body, p_legend,
        ncol = 2, rel_widths = c(1, 0.06)
    )

    ggsave(file.path(out_dir, "Figure_classic_stripes_all_variables.png"),
        fig_stripes,
        width = 15, height = n_panels * 1.1 + 0.8,
        dpi = 400
    )
    cat(sprintf("  -> Saved classic stripes figure (%d panels).\n", n_panels))
} else {
    warning("No stripe panels were generated.")
}
