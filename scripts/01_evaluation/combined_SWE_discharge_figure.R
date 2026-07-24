###############################################################################
# COMBINED SWE + DISCHARGE FIGURE
#
# Layout: 3 columns
#   Col a: Discharge catchment maps (Daily / 7-Day / 14-Day / Monthly)
#   Col b: SWE catchment maps (same 4 windows)
#   Col c: Violin plots (Discharge top, SWE bottom)
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
tsize <- 10
osize <- 9

# ===========================================================================
# 1. LOAD DATA
# ===========================================================================
cat("Loading data...\n")

# SWE stats
path_swe <- file.path(base_dir, "output/swe_diego/2.stats")
swe_daily <- as.data.frame(readRDS(file.path(path_swe, "stats_daily.rds")))
swe_7d <- as.data.frame(readRDS(file.path(path_swe, "stats_7d.rds")))
swe_15d <- as.data.frame(readRDS(file.path(path_swe, "stats_15d.rds")))
swe_month <- as.data.frame(readRDS(file.path(path_swe, "stats_month.rds")))

# Discharge stats (station-level, all windows)
dis_all <- fread(file.path(base_dir, "output/discharge_validation_all_windows.csv"))
dis_daily <- as.data.frame(dis_all[scenario == "Daily"])
dis_7d <- as.data.frame(dis_all[scenario == "7-Day"])
dis_14d <- as.data.frame(dis_all[scenario == "14-Day"])
dis_30d <- as.data.frame(dis_all[scenario == "Monthly"])

# Catchment polygons
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

# Global bounding box from ALL catchments (ensures no clipping)
shp_laea <- st_transform(shp, 3035)
bbox_laea <- st_bbox(shp_laea)

palet <- hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = TRUE, fixup = TRUE)

# ===========================================================================
# 3. SHARED MAP FUNCTION (publication style)
# ===========================================================================

build_map <- function(sp_valid, ttl) {
    ggplot(Europe_laea) +
        geom_sf(fill = "grey87", color = NA) +
        geom_sf(
            data = sp_valid, aes(fill = rho),
            color = "transparent", linewidth = 0.03, alpha = 0.9
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
        ) +
        coord_sf(
            xlim = c(bbox_laea["xmin"], bbox_laea["xmax"]),
            ylim = c(bbox_laea["ymin"], bbox_laea["ymax"]),
            expand = FALSE
        ) +
        labs(subtitle = ttl) +
        theme_minimal(base_size = 10) +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
            panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
            panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
            plot.subtitle = element_text(face = "bold", size = 10, hjust = 0.5),
            legend.title = element_text(face = "bold", size = 10, hjust = 0.5),
            plot.margin = margin(1, 1, 1, 1),
            legend.position = "left"
        )
}

# ===========================================================================
# 4. DISCHARGE MAPS (catchment-level aggregation)
# ===========================================================================
cat("Building discharge maps...\n")

dis_to_catchment <- function(sdf) {
    sdf <- as.data.frame(sdf)
    sdf$upa_ratio <- sdf$station_upa / sdf$catchment_upa
    sdf %>%
        group_by(catch_id) %>%
        summarise(rho = {
            v <- !is.na(spearman_rho)
            if (sum(v) == 0) {
                NA_real_
            } else if (sum(v) == 1) {
                spearman_rho[v]
            } else {
                weighted.mean(spearman_rho[v], w = upa_ratio[v], na.rm = TRUE)
            }
        }, .groups = "drop") %>%
        as.data.frame()
}

make_dis_sf <- function(sdf) {
    cat_df <- dis_to_catchment(sdf)
    cat_df$join_id <- norm_id(cat_df$catch_id)
    sp <- left_join(shp, cat_df, by = "join_id")
    sp <- st_transform(sp, 3035)
    sp[!is.na(sp$rho), ]
}

dis_map1 <- build_map(make_dis_sf(dis_daily), "Daily")
dis_map2 <- build_map(make_dis_sf(dis_7d), "7-Day")
dis_map3 <- build_map(make_dis_sf(dis_14d), "14-Day")
dis_map4 <- build_map(make_dis_sf(dis_30d), "Monthly")

# ===========================================================================
# 5. SWE MAPS (catchment-level)
# ===========================================================================
cat("Building SWE maps...\n")

make_swe_sf <- function(stats_df) {
    stats_df <- as.data.frame(stats_df)
    stats_df$join_id <- norm_id(stats_df$catch_id)
    sp <- left_join(shp, stats_df, by = "join_id")
    sp <- st_transform(sp, 3035)
    sp[!is.na(sp$rho), ]
}

swe_map1 <- build_map(make_swe_sf(swe_daily), "Daily")
swe_map2 <- build_map(make_swe_sf(swe_7d), "7-Day")
swe_map3 <- build_map(make_swe_sf(swe_15d), "14-Day")
swe_map4 <- build_map(make_swe_sf(swe_month), "Monthly")

# ===========================================================================
# 6. VIOLIN PLOTS
# ===========================================================================
cat("Building violin plots...\n")

scen_pal <- c(
    "Daily" = "#bdbdbd", "7-Day" = "#74c4e4",
    "14-Day" = "#2c9e4b", "Monthly" = "#1a3f7a"
)

# Discharge violin
dis_valid <- as.data.frame(dis_all[!is.na(spearman_rho)])
dis_valid$scenario <- factor(dis_valid$scenario,
    levels = c("Daily", "7-Day", "14-Day", "Monthly")
)
dis_ann <- dis_valid %>%
    group_by(scenario) %>%
    summarise(
        med = round(median(spearman_rho, na.rm = TRUE), 2),
        pct = round(mean(spearman_rho >= 0.5, na.rm = TRUE) * 100, 0),
        .groups = "drop"
    )

violin_dis <- ggplot(dis_valid, aes(scenario, spearman_rho, fill = scenario)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", colour = "grey40") +
    geom_boxplot(
        width = 0.12, outlier.shape = NA, colour = "grey25",
        fill = "white", fatten = 2
    ) +
    geom_text(
        data = dis_ann, aes(scenario, med, label = paste0("md=", med)),
        vjust = -0.6, size = 2.2, fontface = "bold", inherit.aes = FALSE
    ) +
    geom_text(
        data = dis_ann, aes(scenario, 0.97, label = paste0(pct, "%\u22650.5")),
        size = 2.2, colour = "grey30", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = scen_pal, guide = "none") +
    scale_y_continuous(limits = c(-0.35, 1.02), breaks = seq(-0.25, 1, 0.25)) +
    labs(title = NULL, x = NULL, y = "Spearman \u03c1") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), plot.margin = margin(2, 2, 2, 2))

# SWE violin
swe_all <- rbind(
    cbind(swe_daily, scenario = "Daily"),
    cbind(swe_7d, scenario = "7-Day"),
    cbind(swe_15d, scenario = "14-Day"),
    cbind(swe_month, scenario = "Monthly")
)
swe_all$scenario <- factor(swe_all$scenario,
    levels = c("Daily", "7-Day", "14-Day", "Monthly")
)
swe_valid <- swe_all[!is.na(swe_all$rho), ]
swe_ann <- swe_valid %>%
    group_by(scenario) %>%
    summarise(
        med = round(median(rho, na.rm = TRUE), 2),
        pct = round(mean(rho >= 0.5, na.rm = TRUE) * 100, 0),
        .groups = "drop"
    )

violin_swe <- ggplot(swe_valid, aes(scenario, rho, fill = scenario)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_violin(alpha = 0.35, trim = FALSE, scale = "width", colour = "grey40") +
    geom_boxplot(
        width = 0.12, outlier.shape = NA, colour = "grey25",
        fill = "white", fatten = 2
    ) +
    geom_text(
        data = swe_ann, aes(scenario, med, label = paste0("md=", med)),
        vjust = -0.6, size = 2.2, fontface = "bold", inherit.aes = FALSE
    ) +
    geom_text(
        data = swe_ann, aes(scenario, 0.97, label = paste0(pct, "%\u22650.5")),
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

dis_map1 <- dis_map1 + ggtitle("(a)") +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 14))
# Column a: Discharge maps (4 stacked)
col_a <- (dis_map1 / dis_map2 / dis_map3 / dis_map4) 


# Column b: SWE maps (4 stacked)
swe_map1 <- swe_map1 + ggtitle("(b)") +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 14))
col_b <- (swe_map1 / swe_map2 / swe_map3 / swe_map4) 


# Column c: Violins (discharge top, SWE bottom) with labels
col_c <- ((violin_dis + labs(title = "(c)") +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 14))) /
    (violin_swe + labs(title = "(d)") +
        theme(plot.title = element_text(face = "bold", hjust = 0, size = 14))))

# Combine all three columns with shared legend and tight spacing
fig_combined <- (col_a | col_b | col_c) +
    plot_layout(widths = c(1, 1, 1.6), guides = "collect") &
    theme(
        legend.position = "left",
        legend.justification = "center",
        plot.margin = margin(2, 2, 2, 2)
    )

# fig_combined <- fig_combined +
#     plot_annotation(
#         title = "HERA validation: River Discharge and Snow Water Equivalent",
#         subtitle = "Spearman \u03c1 at four temporal aggregations",
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

ggsave(file.path(path_out, "Fig_SWE_discharge_combined.png"), fig_combined,
    width = 36, height = 30, units = "cm", dpi = 300, bg = "white"
)

cat(sprintf(
    "Done. Figure saved: %s\n",
    file.path(path_out, "Fig_SWE_discharge_combined.png")
))
