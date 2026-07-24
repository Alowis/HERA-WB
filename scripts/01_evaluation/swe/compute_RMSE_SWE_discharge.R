###############################################################################
# COMPUTE RMSE FOR SWE AND DISCHARGE
#
# Computes normalised RMSE (NRMSE = RMSE / mean_obs) per catchment/station
# at four temporal aggregations (Daily / 7-Day / 14-Day / Monthly).
# Produces maps + violin plots in the same style as the Spearman rho figures.
#
# Inputs:
#   SWE: output/swe_diego/2.stats/SWE_time_series_all_scenarios.rds
#   Discharge: data/river_discharge/*.csv + output/station_catchment_matches.csv
#
# Outputs:
#   output/rmse_validation/RMSE_SWE_all_windows.csv
#   output/rmse_validation/RMSE_discharge_all_windows.csv
#   output/rmse_validation/Fig_RMSE_combined.png
###############################################################################

library(sf)
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(rnaturalearth)

base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
tsize <- 10
osize <- 9

# ===========================================================================
# 1. COMPUTE RMSE FOR SWE (from time series RDS)
# ===========================================================================
cat("[1/4] Computing RMSE for SWE...\n")

path_swe_stats <- file.path(base_dir, "output/swe_diego/2.stats")
swe_ts <- readRDS(file.path(path_swe_stats, "SWE_time_series_all_scenarios.rds"))

compute_rmse_swe <- function(obs_dt, mod_dt, label) {
    obs_dt <- as.data.table(obs_dt)
    mod_dt <- as.data.table(mod_dt)
    ids <- setdiff(names(mod_dt), "date")
    rbindlist(lapply(ids, function(id) {
        x <- mod_dt[[id]]
        y <- obs_dt[[id]]
        valid <- which(!is.na(x) & !is.na(y))
        n <- length(valid)
        if (n < 30) {
            return(list(
                catch_id = id, nrmse = NA_real_, nmbe = NA_real_,
                mean_obs = NA_real_, n_obs = n, scenario = label
            ))
        }
        xv <- x[valid]
        yv <- y[valid]
        rmse_val <- sqrt(mean((xv - yv))^2)
        lyp <-length(which(yv>0))
        mean_obs <- mean(yv)
        smean_obs <- mean(yv[which(yv>0)])
        nrmse <-  (rmse_val)/smean_obs
        nmbe <-  mean((xv - yv))/smean_obs
        # rmse_mm <- if (mean_obs != 0) rmse_val / abs(mean_obs) else rmse_val / abs(smean_obs)
        list(
            catch_id = id, nrmse = nrmse, nmbe = nmbe,
            mean_obs = smean_obs, n_obs = lyp, scenario = label
        )
    }))
}

swe_rmse <- rbindlist(list(
    compute_rmse_swe(swe_ts$daily$obs, swe_ts$daily$mod, "Daily"),
    compute_rmse_swe(swe_ts$`7d`$obs, swe_ts$`7d`$mod, "7-Day"),
    compute_rmse_swe(swe_ts$`15d`$obs, swe_ts$`15d`$mod, "14-Day"),
    compute_rmse_swe(swe_ts$monthly$obs, swe_ts$monthly$mod, "Monthly")
))
swe_rmse$scenario <- factor(swe_rmse$scenario,
    levels = c("Daily", "7-Day", "14-Day", "Monthly")
)

cat(sprintf("  SWE RMSE computed: %d rows\n", nrow(swe_rmse)))

# ===========================================================================
# 2. COMPUTE RMSE FOR DISCHARGE (station-level, same loop as validation)
# ===========================================================================
cat("[2/4] Computing RMSE for Discharge...\n")

# Load discharge data
matches <- read.csv(file.path(base_dir, "output/station_catchment_matches.csv"),
    stringsAsFactors = FALSE
)
stations_raw <- read.csv(file.path(base_dir, "Stations_ValidationF.csv"),
    stringsAsFactors = FALSE
)

# Load observed
Q_data <- read.csv(file.path(base_dir, "data/river_discharge/Q_19502020.csv"), header = FALSE)
time <- as.Date(as.numeric(Q_data$V2) - 1, origin = "0000-01-01")[-c(1, 2)]
Q_data <- Q_data[-1, -1]
Station_obs_IDs <- as.numeric(as.vector(t(Q_data[1, ])))[-1]
Q_data <- Q_data[-1, ]
Q_data <- as.data.frame(lapply(Q_data, as.numeric))
rm_obs <- which(lubridate::year(time) == 1950)

# Load simulated
Q_sim <- read.csv(file.path(base_dir, "data/river_discharge/HERA_Val4_19502020.csv"), header = FALSE)
Q_sim <- Q_sim[-1, -1]
Station_sim_IDs <- as.numeric(as.vector(t(Q_sim[1, ])))
Q_sim <- Q_sim[-1, ]
Q_sim <- as.data.frame(lapply(Q_sim, as.numeric))

# Corrected stations
Q_cor_path <- file.path(base_dir, "data/river_discharge/HERA_CorStatv3_19512020.csv")
has_corrected <- file.exists(Q_cor_path)
if (has_corrected) {
    Q_cor <- read.csv(Q_cor_path, header = FALSE)
    Q_cor <- Q_cor[-1, -1]
    Station_cor_IDs <- as.numeric(as.vector(t(Q_cor[1, ])))
    Q_cor <- Q_cor[-1, ]
    Q_cor <- as.data.frame(lapply(Q_cor, as.numeric))
} else {
    Station_cor_IDs <- numeric(0)
}

date_obs <- seq(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")
date_sim <- date_obs
date_cor <- seq(as.Date("1951-01-01"), as.Date("2020-12-31"), by = "day")

# RMSE computation function per station with optional smoothing
compute_rmse_station <- function(matches, k = 1) {
    results <- vector("list", nrow(matches))
    for (i in seq_len(nrow(matches))) {
        sid <- as.numeric(matches$station_id[i])
        cid <- as.character(matches$catch_id[i])
        upa <- as.numeric(matches$station_upa[i])
        use_corrected <- has_corrected && (sid %in% Station_cor_IDs)
        sim_col <- if (use_corrected) which(Station_cor_IDs == sid) else which(Station_sim_IDs == sid)
        obs_col <- which(Station_obs_IDs == sid)
        if (length(sim_col) == 0 || length(obs_col) == 0) {
            results[[i]] <- data.frame(
                catch_id = cid, station_id = sid,
                rmse = NA, nrmse = NA, mean_obs = NA, n_obs = 0,
                stringsAsFactors = FALSE
            )
            next
        }
        sim_vec <- if (use_corrected) Q_cor[, sim_col[1]] else Q_sim[, sim_col[1]]
        obs_vec <- Q_data[-rm_obs, obs_col[1] + 1]
        dt_sim <- data.table(date = if (use_corrected) date_cor else date_sim, sim = as.numeric(sim_vec))
        dt_obs <- data.table(date = date_obs, obs = as.numeric(obs_vec))
        paired <- merge(dt_sim, dt_obs, by = "date", all = FALSE)
        paired <- paired[!is.na(sim) & !is.na(obs)]
        if (nrow(paired) < 365) {
            results[[i]] <- data.frame(
                catch_id = cid, station_id = sid,
                nrmse = NA, nmbe = NA, mean_obs = NA, n_obs = nrow(paired),
                stringsAsFactors = FALSE
            )
            next
        }
        sv <- paired$sim
        ov <- paired$obs
        if (k > 1) {
            sv <- data.table::frollmean(sv, n = k, align = "right", na.rm = TRUE)
            ov <- data.table::frollmean(ov, n = k, align = "right", na.rm = TRUE)
        }
        valid <- which(!is.na(sv) & !is.na(ov))
        rmse_val <- sqrt(mean((sv[valid] - ov[valid])^2))
        mean_obs <- mean(ov[valid])
        nrmse <-  (rmse_val)/mean_obs
        nmbe <-  mean((sv[valid] - ov[valid]))/mean_obs
        
        results[[i]] <- data.frame(
            catch_id = cid, station_id = sid,
            nrmse = nrmse, nmbe = nmbe,
            mean_obs = mean_obs, n_obs = length(valid),
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, results)
}

dis_rmse_daily <- compute_rmse_station(matches, k = 1)
dis_rmse_daily$scenario <- "Daily"
cat("  Daily done\n")
dis_rmse_7d <- compute_rmse_station(matches, k = 7)
dis_rmse_7d$scenario <- "7-Day"
cat("  7-Day done\n")
dis_rmse_14d <- compute_rmse_station(matches, k = 14)
dis_rmse_14d$scenario <- "14-Day"
cat("  14-Day done\n")
dis_rmse_30d <- compute_rmse_station(matches, k = 30)
dis_rmse_30d$scenario <- "Monthly"
cat("  Monthly done\n")

dis_rmse <- rbind(dis_rmse_daily, dis_rmse_7d, dis_rmse_14d, dis_rmse_30d)
dis_rmse$scenario <- factor(dis_rmse$scenario,
    levels = c("Daily", "7-Day", "14-Day", "Monthly")
)

# ===========================================================================
# 3. SAVE CSVs
# ===========================================================================
cat("[3/4] Saving results...\n")
path_out <- file.path(base_dir, "output", "rmse_validation")
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

fwrite(swe_rmse, file.path(path_out, "RMSE_SWE_all_windows.csv"))
fwrite(dis_rmse, file.path(path_out, "RMSE_discharge_all_windows.csv"))

# ===========================================================================
# 4. FIGURES (same layout as Spearman rho combined figure)
# ===========================================================================
cat("[4/4] Generating figures...\n")

# Load catchments for maps
file_shp <- file.path(base_dir, "data/catchments_analysis_final_v3.gpkg")
shp <- st_read(file_shp, quiet = TRUE)
norm_id <- function(x) {
    x <- as.character(x)
    suppressWarnings(ifelse(grepl("^\\d+$", x), as.character(as.integer(x)), x))
}
shp$join_id <- norm_id(shp$catch_id)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
Europe_laea <- st_transform(world, crs = 3035)
shp_laea <- st_transform(shp, 3035)
bbox_laea <- st_bbox(shp_laea)

# NRMSE palette (red = bad, green = good, reversed from rho)
palet_rmse <- hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE)

# Shared map function for NRMSE
build_rmse_map <- function(sp_valid, ttl) {
    ggplot(Europe_laea) +
        geom_sf(fill = "grey87", color = NA) +
        geom_sf(
            data = sp_valid, aes(fill = nrmse),
            color = "grey40", linewidth = 0.03, alpha = 0.9
        ) +
        geom_sf(data = Europe_laea, fill = NA, color = "grey30", linewidth = 0.15) +
        scale_fill_gradientn(
            colors = palet_rmse, limits = c(0, 2), oob = squish,
            breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,4,5,6,7,8,9),
            name = "RMSE",
            guide = guide_colorbar(
                direction = "horizontal", title.position = "top",
                barwidth = 8, barheight = 0.5
            )
        ) +
        coord_sf(
            xlim = c(bbox_laea["xmin"], bbox_laea["xmax"]),
            ylim = c(bbox_laea["ymin"], bbox_laea["ymax"]), expand = FALSE
        ) +
        labs(subtitle = ttl) +
        theme_minimal(base_size = 8) +
        theme(
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
            panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
            panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
            plot.subtitle = element_text(face = "bold", size = tsize, hjust = 0.5),
            plot.margin = margin(1, 1, 1, 1), legend.position = "bottom"
        )
}


# NMBE palette
palet_nmbe <- hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE)

#Shared map function for NMBE
build_nmbe_map <- function(sp_valid, ttl) {
  ggplot(Europe_laea) +
    geom_sf(fill = "grey87", color = NA) +
    geom_sf(
      data = sp_valid, aes(fill = nmbe),
      color = "grey40", linewidth = 0.03, alpha = 0.9
    ) +
    geom_sf(data = Europe_laea, fill = NA, color = "grey30", linewidth = 0.15) +
    scale_fill_gradientn(
      colors = palet_nmbe, limits = c(-2, 2), oob = squish,
      breaks = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,4,5,6,7,8,9),
      name = "NMBE",
      guide = guide_colorbar(
        direction = "horizontal", title.position = "top",
        barwidth = 8, barheight = 0.5
      )
    ) +
    coord_sf(
      xlim = c(bbox_laea["xmin"], bbox_laea["xmax"]),
      ylim = c(bbox_laea["ymin"], bbox_laea["ymax"]), expand = FALSE
    ) +
    labs(subtitle = ttl) +
    theme_minimal(base_size = 8) +
    theme(
      axis.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
      panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
      plot.subtitle = element_text(face = "bold", size = tsize, hjust = 0.5),
      plot.margin = margin(1, 1, 1, 1), legend.position = "bottom"
    )
}
# --- Discharge RMSE maps (aggregate to catchment) ---
dis_rmse$upa_ratio <- matches$station_upa[match(dis_rmse$station_id, matches$station_id)] /
    matches$catchment_upa[match(dis_rmse$station_id, matches$station_id)]

make_dis_rmse_sf <- function(scenario_label) {
    sub <- dis_rmse[dis_rmse$scenario == scenario_label & !is.na(dis_rmse$nrmse), ]
    cat_agg <- sub %>%
        group_by(catch_id) %>%
        summarise(
            nrmse = weighted.mean(nrmse, w = pmax(upa_ratio, 0.01, na.rm = TRUE), na.rm = TRUE),
            .groups = "drop"
        ) %>%
        as.data.frame()
    cat_agg$join_id <- norm_id(cat_agg$catch_id)
    sp <- left_join(shp, cat_agg, by = "join_id")
    sp <- st_transform(sp, 3035)
    sp[!is.na(sp$nrmse), ]
}

dis_rmap1 <- build_rmse_map(make_dis_rmse_sf("Daily"), "Daily")
dis_rmap2 <- build_rmse_map(make_dis_rmse_sf("7-Day"), "7-Day")
dis_rmap3 <- build_rmse_map(make_dis_rmse_sf("14-Day"), "14-Day")
dis_rmap4 <- build_rmse_map(make_dis_rmse_sf("Monthly"), "Monthly")


make_dis_nmbe_sf <- function(scenario_label) {
  sub <- dis_rmse[dis_rmse$scenario == scenario_label & !is.na(dis_rmse$nmbe), ]
  cat_agg <- sub %>%
    group_by(catch_id) %>%
    summarise(
      nmbe = weighted.mean(nmbe, w = pmax(upa_ratio, 0.01, na.rm = TRUE), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame()
  cat_agg$join_id <- norm_id(cat_agg$catch_id)
  sp <- left_join(shp, cat_agg, by = "join_id")
  sp <- st_transform(sp, 3035)
  sp[!is.na(sp$nmbe), ]
}

dis_bmap1 <- build_nmbe_map(make_dis_nmbe_sf("Daily"), "Daily")
dis_bmap2 <- build_nmbe_map(make_dis_nmbe_sf("7-Day"), "7-Day")
dis_bmap3 <- build_nmbe_map(make_dis_nmbe_sf("14-Day"), "14-Day")
dis_bmap4 <- build_nmbe_map(make_dis_nmbe_sf("Monthly"), "Monthly")


# --- SWE NRMSE maps ---
make_swe_rmse_sf <- function(scenario_label) {
    sub <- swe_rmse[scenario == scenario_label & !is.na(nrmse)]
    sub$join_id <- norm_id(sub$catch_id)
    sp <- left_join(shp, as.data.frame(sub), by = "join_id")
    sp <- st_transform(sp, 3035)
    sp[!is.na(sp$nrmse), ]
}

swe_rmap1 <- build_rmse_map(make_swe_rmse_sf("Daily"), "Daily")
swe_rmap2 <- build_rmse_map(make_swe_rmse_sf("7-Day"), "7-Day")
swe_rmap3 <- build_rmse_map(make_swe_rmse_sf("14-Day"), "14-Day")
swe_rmap4 <- build_rmse_map(make_swe_rmse_sf("Monthly"), "Monthly")

# --- Violin plots for NRMSE ---
scen_pal <- c(
    "Daily" = "#bdbdbd", "7-Day" = "#74c4e4",
    "14-Day" = "#2c9e4b", "Monthly" = "#1a3f7a"
)

# Discharge NRMSE violin
dis_rmse_valid <- dis_rmse[!is.na(dis_rmse$nrmse) & abs(dis_rmse$nrmse) < 5, ]
dis_rmse_ann <- dis_rmse_valid %>%
    group_by(scenario) %>%
    summarise(med = round(median(nrmse, na.rm = TRUE), 2), .groups = "drop")

violin_dis_rmse <- ggplot(dis_rmse_valid, aes(scenario, nrmse, fill = scenario)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", colour = "grey40") +
    geom_boxplot(
        width = 0.12, outlier.shape = NA, colour = "grey25",
        fill = "white", fatten = 2
    ) +
    geom_text(
        data = dis_rmse_ann, aes(scenario, med, label = paste0("md=", med)),
        vjust = -0.6, size = 2.2, fontface = "bold", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = scen_pal, guide = "none") +
    scale_y_continuous(limits = c(0, 3), breaks = seq(-3, 3, 0.5)) +
    labs(title = NULL, x = NULL, y = "NRMSE") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), plot.margin = margin(2, 2, 2, 2))

# SWE NRMSE violin
swe_rmse_valid <- swe_rmse[!is.na(nrmse) & nrmse < 5]
swe_rmse_ann <- swe_rmse_valid[, .(med = round(median(nrmse, na.rm = TRUE), 2)), by = scenario]

violin_swe_rmse <- ggplot(swe_rmse_valid, aes(scenario, nrmse, fill = scenario)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", colour = "grey40") +
    geom_boxplot(
        width = 0.12, outlier.shape = NA, colour = "grey25",
        fill = "white", fatten = 2
    ) +
    geom_text(
        data = as.data.frame(swe_rmse_ann),
        aes(scenario, med, label = paste0("md=", med)),
        vjust = -0.6, size = 2.2, fontface = "bold", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = scen_pal, guide = "none") +
    scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
    labs(title = NULL, x = NULL, y = "NRMSE") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), plot.margin = margin(2, 2, 2, 2))


violin_swe_nmbe <- ggplot(swe_rmse_valid, aes(scenario, nmbe, fill = scenario)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_violin(alpha = 0.35, trim = FALSE, scale = "width", colour = "grey40") +
  geom_boxplot(
    width = 0.12, outlier.shape = NA, colour = "grey25",
    fill = "white", fatten = 2
  ) +
  geom_text(
    data = as.data.frame(swe_rmse_ann),
    aes(scenario, med, label = paste0("md=", med)),
    vjust = -0.6, size = 2.2, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = scen_pal, guide = "none") +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 0.5)) +
  labs(title = NULL, x = NULL, y = "NMBE") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), plot.margin = margin(2, 2, 2, 2))

violin_swe_nmbe

# --- Assemble (same layout as rho figure) ---
dis_rmap1 <- dis_rmap1 + ggtitle("(a) ")
swe_rmap1 <- swe_rmap1 + ggtitle("(b) ")

col_a <- dis_rmap1 / dis_rmap2 / dis_rmap3 / dis_rmap4
col_b <- swe_rmap1 / swe_rmap2 / swe_rmap3 / swe_rmap4
col_c <- (violin_dis_rmse + labs(title = "(c) ") +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 10))) /
    (violin_swe_rmse + labs(title = "(d) ") +
        theme(plot.title = element_text(face = "bold", hjust = 0, size = 10)))

fig_rmse <- (col_a | col_b | col_c) +
    plot_layout(widths = c(1, 1, 1.6), guides = "collect") &
    theme(
        legend.position = "bottom",
        legend.justification = "center",
        plot.margin = margin(0, -2, 0, -2)
    )

fig_rmse <- fig_rmse +
    plot_annotation(
        title = "HERA validation: Normalised RMSE (Discharge & SWE)",
        subtitle = "NRMSE = RMSE / mean(obs) | Dashed line = NRMSE = 1",
        theme = theme(
            plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
            plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0.5),
            plot.margin = margin(5, 0, 5, 0)
        )
    )

ggsave(file.path(path_out, "Fig_RMSE_combined.png"), fig_rmse,
    width = 36, height = 30, units = "cm", dpi = 300, bg = "white"
)

cat(sprintf("Done. Outputs in: %s\n", path_out))
