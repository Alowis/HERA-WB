###############################################################################
# COMBINED AET + SOIL MOISTURE FIGURE
#
# Layout: 3 columns
#   Col a: AET catchment maps (Daily / 7-Day / 14-Day / Monthly)
#   Col b: Soil Moisture catchment maps (same 4 windows)
#   Col c: Violin plots (AET top, SM bottom)
#
# Single shared colour legend at the bottom.
###############################################################################

library(sf)
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(scales)
library(rnaturalearth)

base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
sm_stats_dir <- "Z:/ClimateRun4/nahaUsers/tilloal/SoilMositure_2026_06_03/2.Diego_Analysis/0.Stats_time_windows/"
tsize <- 10
osize <- 9

# ===========================================================================
# 1. LOAD DATA
# ===========================================================================
cat("Loading data...\n")

# --- AET stats ---
path_aet_stats <- file.path(base_dir, "output/aet_diego/2.stats")
aet_daily <- as.data.frame(readRDS(file.path(path_aet_stats, "stats_daily.rds")))
aet_7d <- as.data.frame(readRDS(file.path(path_aet_stats, "stats_7d.rds")))
aet_15d <- as.data.frame(readRDS(file.path(path_aet_stats, "stats_15d.rds")))
aet_month <- as.data.frame(readRDS(file.path(path_aet_stats, "stats_month.rds")))

# --- Soil Moisture stats (from RDS on network drive) ---
sm_daily <- as.data.frame(readRDS(file.path(sm_stats_dir, "stats_daily.rds")))
sm_7d <- as.data.frame(readRDS(file.path(sm_stats_dir, "stats_7d.rds")))
sm_15d <- as.data.frame(readRDS(file.path(sm_stats_dir, "stats_15d.rds")))
sm_month <- as.data.frame(readRDS(file.path(sm_stats_dir, "stats_month.rds")))

# Fix SM catch_id (remove leading X from R column name conversion)
sm_daily$catch_id <- sub("^X", "", sm_daily$catch_id)
sm_7d$catch_id <- sub("^X", "", sm_7d$catch_id)
sm_15d$catch_id <- sub("^X", "", sm_15d$catch_id)
sm_month$catch_id <- sub("^X", "", sm_month$catch_id)

# --- Catchment polygons ---
file_shp <- file.path(base_dir, "data/catchments_analysis_final_v3.gpkg")
shp <- st_read(file_shp, quiet = TRUE)
norm_id <- function(x) {
    x <- as.character(x)
    suppressWarnings(ifelse(grepl("^\\d+$", x), as.character(as.integer(x)), x))
}
shp$join_id <- norm_id(shp$catch_id)

# ===========================================================================
# 2. BASEMAP + BOUNDING BOX
# ===========================================================================
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
Europe_laea <- st_transform(world, crs = 3035)
shp_laea <- st_transform(shp, 3035)
bbox_laea <- st_bbox(shp_laea)

palet <- hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = TRUE, fixup = TRUE)

# ===========================================================================
# 3. SHARED MAP FUNCTION
# ===========================================================================
build_map <- function(sp_valid, ttl) {
    ggplot(Europe_laea) +
        geom_sf(fill = "grey87", color = NA) +
        geom_sf(
            data = sp_valid, aes(fill = rho),
            color = "grey40", linewidth = 0.03, alpha = 0.9
        ) +
        geom_sf(data = Europe_laea, fill = NA, color = "grey30", linewidth = 0.15) +
        scale_fill_gradientn(
          colors = palet, limits = c(0, 1), oob = squish,
          breaks = c(0, 0.25, 0.5, 0.75, 1),
          name = "Spearman \u03c1",
          guide = guide_colorbar(
            direction = "vertical", title.position = "top",
            barwidth = 1.5, barheight = 20
          )
        )  +
        coord_sf(
            xlim = c(bbox_laea["xmin"], bbox_laea["xmax"]),
            ylim = c(bbox_laea["ymin"], bbox_laea["ymax"]), expand = FALSE
        ) +
        labs(subtitle = ttl) +
        theme_minimal(base_size = 10) +
        theme(
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
            panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
            legend.title = element_text(face = "bold", size = 10, hjust = 0.5),
            panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
            plot.subtitle = element_text(face = "bold", size = tsize, hjust = 0.5),
            plot.margin = margin(1, 1, 1, 1), legend.position = "bottom"
        )
}

# Helper to join stats to shapefile
make_sf <- function(stats_df) {
    stats_df <- as.data.frame(stats_df)
    stats_df$join_id <- norm_id(stats_df$catch_id)
    sp <- left_join(shp, stats_df, by = "join_id")
    sp <- st_transform(sp, 3035)
    sp[!is.na(sp$rho), ]
}

# ===========================================================================
# 4. AET MAPS
# ===========================================================================
cat("Building AET maps...\n")

aet_map1 <- build_map(make_sf(aet_daily), "Daily")
aet_map2 <- build_map(make_sf(aet_7d), "7-Day")
aet_map3 <- build_map(make_sf(aet_15d), "14-Day")
aet_map4 <- build_map(make_sf(aet_month), "Monthly")

# ===========================================================================
# 5. SOIL MOISTURE MAPS
# ===========================================================================
cat("Building Soil Moisture maps...\n")

sm_map1 <- build_map(make_sf(sm_daily), "Daily")
sm_map2 <- build_map(make_sf(sm_7d), "7-Day")
sm_map3 <- build_map(make_sf(sm_15d), "14-Day")
sm_map4 <- build_map(make_sf(sm_month), "Monthly")

# ===========================================================================
# 6. VIOLIN PLOTS
# ===========================================================================
cat("Building violin plots...\n")

scen_pal <- c(
    "Daily" = "#bdbdbd", "7-Day" = "#74c4e4",
    "14-Day" = "#2c9e4b", "Monthly" = "#1a3f7a"
)

# AET violin
aet_all <- rbind(
    cbind(aet_daily, scenario = "Daily"),
    cbind(aet_7d, scenario = "7-Day"),
    cbind(aet_15d, scenario = "14-Day"),
    cbind(aet_month, scenario = "Monthly")
)
aet_all$scenario <- factor(aet_all$scenario,
    levels = c("Daily", "7-Day", "14-Day", "Monthly")
)
aet_valid <- aet_all[!is.na(aet_all$rho), ]

aet_ann <- aet_valid %>%
    group_by(scenario) %>%
    summarise(
        med = round(median(rho, na.rm = TRUE), 2),
        pct = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 0),
        .groups = "drop"
    )

violin_aet <- ggplot(aet_valid, aes(scenario, rho, fill = scenario)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", colour = "grey40") +
    geom_boxplot(
        width = 0.12, outlier.shape = NA, colour = "grey25",
        fill = "white", fatten = 2
    ) +
    geom_text(
        data = aet_ann, aes(scenario, med, label = paste0("md=", med)),
        vjust = -0.6, size = 2.2, fontface = "bold", inherit.aes = FALSE
    ) +
    geom_text(
        data = aet_ann, aes(scenario, 0.97, label = paste0(pct, "%\u22650.5")),
        size = 2.2, colour = "grey30", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = scen_pal, guide = "none") +
    scale_y_continuous(limits = c(-0.35, 1.02), breaks = seq(-0.25, 1, 0.25)) +
    labs(title = NULL, x = NULL, y = "Spearman \u03c1") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), plot.margin = margin(2, 2, 2, 2))

# Soil Moisture violin
sm_all <- rbind(
    cbind(sm_daily, scenario = "Daily"),
    cbind(sm_7d, scenario = "7-Day"),
    cbind(sm_15d, scenario = "14-Day"),
    cbind(sm_month, scenario = "Monthly")
)
sm_all$scenario <- factor(sm_all$scenario,
    levels = c("Daily", "7-Day", "14-Day", "Monthly")
)
sm_valid <- sm_all[!is.na(sm_all$rho), ]

sm_ann <- sm_valid %>%
    group_by(scenario) %>%
    summarise(
        med = round(median(rho, na.rm = TRUE), 2),
        pct = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 0),
        .groups = "drop"
    )

violin_sm <- ggplot(sm_valid, aes(scenario, rho, fill = scenario)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", colour = "grey40") +
    geom_boxplot(
        width = 0.12, outlier.shape = NA, colour = "grey25",
        fill = "white", fatten = 2
    ) +
    geom_text(
        data = sm_ann, aes(scenario, med, label = paste0("md=", med)),
        vjust = -0.6, size = 2.2, fontface = "bold", inherit.aes = FALSE
    ) +
    geom_text(
        data = sm_ann, aes(scenario, 0.97, label = paste0(pct, "%\u22650.5")),
        size = 2.2, colour = "grey30", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = scen_pal, guide = "none") +
    scale_y_continuous(limits = c(-0.35, 1.02), breaks = seq(-0.25, 1, 0.25)) +
    labs(title = NULL, x = NULL, y = "Spearman \u03c1") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), plot.margin = margin(2, 2, 2, 2))

# ===========================================================================
# 7. ASSEMBLE COMBINED FIGURE
# ===========================================================================
cat("Assembling combined figure...\n")

# Add (a) and (b) labels to first map in each column
aet_map1 <- aet_map1 + ggtitle("(a)")+
theme(plot.title = element_text(face = "bold", hjust = 0, size = 14))
sm_map1 <- sm_map1 + ggtitle("(b)")+
theme(plot.title = element_text(face = "bold", hjust = 0, size = 14))

# Column a: AET maps stacked
col_a <- aet_map1 / aet_map2 / aet_map3 / aet_map4

# Column b: SM maps stacked
col_b <- sm_map1 / sm_map2 / sm_map3 / sm_map4

# Column c: Violins with (c) and (d) labels
col_c <- (violin_aet + labs(title = "(c) ") +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 14))) /
    (violin_sm + labs(title = "(d) ") +
        theme(plot.title = element_text(face = "bold", hjust = 0, size = 14)))

# Combine all three columns with shared legend
fig_combined <- (col_a | col_b | col_c) +
    plot_layout(widths = c(1, 1, 1.6), guides = "collect") &
    theme(
        legend.position = "left",
        legend.justification = "center",
        plot.margin = margin(2, 2, 2, 2)
    )

# fig_combined <- fig_combined +
#     plot_annotation(
#         title = "HERA validation: Actual Evapotranspiration and Soil Moisture",
#         subtitle = "Spearman \u03c1 at four temporal aggregations (LISFLOOD vs GLEAM / ESA CCI)",
#         theme = theme(
#             plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
#             plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0.5),
#             plot.margin = margin(5, 0, 5, 0)
#         )
#     )

# ===========================================================================
# 8. SAVE
# ===========================================================================
path_out <- file.path(base_dir, "output", "combined_validation")
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(path_out, "Fig_AET_SM_combined.png"), fig_combined,
    width = 36, height = 30, units = "cm", dpi = 300, bg = "white"
)

cat(sprintf(
    "Done. Figure saved: %s\n",
    file.path(path_out, "Fig_AET_SM_combined.png")
))
