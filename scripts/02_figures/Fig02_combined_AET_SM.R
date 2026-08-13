###############################################################################
# COMBINED AET + SOIL MOISTURE FIGURE
#
# Layout: 3 columns
#   Col a: AET catchment maps (Daily / Top20 / Bot20 / Regime)
#   Col b: Soil Moisture catchment maps (same 4 panels)
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
tsize <- 15
osize <- 14

# ===========================================================================
# 1. LOAD DATA
# ===========================================================================
cat("Loading data...\n")

# --- AET stats ---
path_aet_stats <- file.path(base_dir, "output/aet_diego/2.stats")
aet_daily <- as.data.frame(readRDS(file.path(path_aet_stats, "stats_daily.rds")))
aet_tex <- as.data.frame(readRDS(file.path(path_aet_stats, "pearson_extremes_daily.rds")))
aet_reg <- as.data.frame(readRDS(file.path(path_aet_stats, "stats_regime.rds")))

aet_extra <- data.frame(aet_daily, aet_tex, aet_reg)

# --- Soil Moisture stats (from RDS on network drive) ---
sm_daily <- as.data.frame(readRDS(file.path(sm_stats_dir, "stats_daily_snowmasked.rds")))
sm_tex <- as.data.frame(readRDS(file.path(sm_stats_dir, "pearson_extremes_daily.rds")))
sm_reg <- as.data.frame(readRDS(file.path(sm_stats_dir, "stats_regime.rds")))

sm_extra <- data.frame(sm_daily, sm_tex, sm_reg)

# Fix SM catch_id (remove leading X from R column name conversion)
sm_extra$catch_id <- sub("^X", "", sm_extra$catch_id)
if ("catch_id.1" %in% names(sm_extra)) sm_extra$catch_id.1 <- sub("^X", "", sm_extra$catch_id.1)

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
            data = sp_valid, aes(fill = var),
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
            ylim = c(bbox_laea["ymin"], bbox_laea["ymax"]), expand = FALSE
        ) +
        labs(subtitle = ttl) +
        theme_minimal(base_size = 14) +
        theme(
            axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
            panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
            panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
            plot.margin = margin(1, 1, 1, 1), legend.position = "left"
        )
}

# Helper to join stats to shapefile with selectable variable
make_sf <- function(stats_df, vsel) {
    stats_df <- as.data.frame(stats_df)
    mn <- which(names(stats_df) == vsel)
    names(stats_df)[mn] <- "var"
    stats_df$join_id <- norm_id(stats_df$catch_id)
    sp <- left_join(shp, stats_df, by = "join_id")
    sp <- st_transform(sp, 3035)
    sp[!is.na(sp$var), ]
}

# ===========================================================================
# 4. AET MAPS
# ===========================================================================
cat("Building AET maps...\n")

aet_map1 <- build_map(make_sf(aet_extra, "rho"), "Daily")
aet_map2 <- build_map(make_sf(aet_extra, "pearson_top20"), "Top20")
aet_map3 <- build_map(make_sf(aet_extra, "pearson_bot20"), "Bot20")
aet_map4 <- build_map(make_sf(aet_extra, "rho.1"), "Regime")

# ===========================================================================
# 5. SOIL MOISTURE MAPS
# ===========================================================================
cat("Building Soil Moisture maps...\n")

sm_map1 <- build_map(make_sf(sm_extra, "rho"), "Daily")
sm_map2 <- build_map(make_sf(sm_extra, "pearson_top20"), "Top20")
sm_map3 <- build_map(make_sf(sm_extra, "pearson_bot20"), "Bot20")
sm_map4 <- build_map(make_sf(sm_extra, "rho.1"), "Regime")

# ===========================================================================
# 6. VIOLIN PLOTS
# ===========================================================================
cat("Building violin plots...\n")

levels_v2 <- c("Daily", "Top20", "Bot20", "Regime")
scen_pal <- c(
    "Daily" = "#bdbdbd", "Top20" = "#74c4e4",
    "Bot20" = "#2c9e4b", "Regime" = "#1a3f7a"
)

# --- AET violin ---
daily_long_aet <- data.frame(
    catch_id = aet_extra$catch_id,
    rho = aet_extra$rho,
    p_val = aet_extra$p_val,
    n_eff = aet_extra$n_eff,
    status = aet_extra$status,
    scenario = "Daily"
)

regime_long_aet <- data.frame(
    catch_id = aet_extra$catch_id,
    rho = aet_extra$rho.1,
    p_val = aet_extra$p_val,
    n_eff = aet_extra$n_eff,
    status = aet_extra$status,
    scenario = "Regime"
)

top20_long_aet <- data.frame(
    catch_id = aet_extra$catch_id,
    rho = aet_extra$pearson_top20,
    p_val = NA_real_,
    n_eff = NA_real_,
    status = aet_extra$status.1,
    scenario = "Top20"
)

bot20_long_aet <- data.frame(
    catch_id = aet_extra$catch_id,
    rho = aet_extra$pearson_bot20,
    p_val = NA_real_,
    n_eff = NA_real_,
    status = aet_extra$status.1,
    scenario = "Bot20"
)

aet_all <- rbind(daily_long_aet, regime_long_aet, top20_long_aet, bot20_long_aet)
aet_all$scenario <- factor(aet_all$scenario, levels = levels_v2)
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
        vjust = 1.4, size = 2.6, fontface = "bold", inherit.aes = FALSE
    ) +
    geom_text(
        data = aet_ann, aes(scenario, 1.01, label = paste0(pct, "%\u22650.5")),
        size = 2.6, colour = "grey30", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = scen_pal, guide = "none") +
    scale_y_continuous(limits = c(-0.35, 1.02), breaks = seq(-0.25, 1, 0.25)) +
    labs(title = NULL, x = NULL, y = "Spearman \u03c1") +
    theme_bw(base_size = 15) +
    theme(panel.grid.minor = element_blank(), plot.margin = margin(2, 2, 2, 2))

# --- Soil Moisture violin ---
daily_long_sm <- data.frame(
    catch_id = sm_extra$catch_id,
    rho = sm_extra$rho,
    p_val = sm_extra$p_val,
    n_eff = sm_extra$n_eff,
    status = sm_extra$status,
    scenario = "Daily"
)

regime_long_sm <- data.frame(
    catch_id = sm_extra$catch_id,
    rho = sm_extra$rho.1,
    p_val = sm_extra$p_val,
    n_eff = sm_extra$n_eff,
    status = sm_extra$status,
    scenario = "Regime"
)

top20_long_sm <- data.frame(
    catch_id = sm_extra$catch_id,
    rho = sm_extra$pearson_top20,
    p_val = NA_real_,
    n_eff = NA_real_,
    status = sm_extra$status.1,
    scenario = "Top20"
)

bot20_long_sm <- data.frame(
    catch_id = sm_extra$catch_id,
    rho = sm_extra$pearson_bot20,
    p_val = NA_real_,
    n_eff = NA_real_,
    status = sm_extra$status.1,
    scenario = "Bot20"
)

sm_all <- rbind(daily_long_sm, regime_long_sm, top20_long_sm, bot20_long_sm)
sm_all$scenario <- factor(sm_all$scenario, levels = levels_v2)
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
        vjust = 1.4, size = 2.6, fontface = "bold", inherit.aes = FALSE
    ) +
    geom_text(
        data = sm_ann, aes(scenario, 1.01, label = paste0(pct, "%\u22650.5")),
        size = 2.6, colour = "grey20", inherit.aes = FALSE
    ) +
    scale_fill_manual(values = scen_pal, guide = "none") +
    scale_y_continuous(limits = c(-0.35, 1.02), breaks = seq(-0.25, 1, 0.25)) +
    labs(title = NULL, x = NULL, y = "Spearman \u03c1") +
    theme_bw(base_size = 15) +
    theme(panel.grid.minor = element_blank(), plot.margin = margin(2, 2, 2, 2))

# ===========================================================================
# 7. ASSEMBLE COMBINED FIGURE
# ===========================================================================
cat("Assembling combined figure...\n")

# Add (a) and (b) labels to first map in each column
aet_map1 <- aet_map1 + ggtitle("(a)") +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 18))
sm_map1 <- sm_map1 + ggtitle("(b)") +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 18))

# Column a: AET maps stacked
col_a <- aet_map1 / aet_map2 / aet_map3 / aet_map4

# Column b: SM maps stacked
col_b <- sm_map1 / sm_map2 / sm_map3 / sm_map4

# Column c: Violins with (c) and (d) labels
col_c <- (violin_aet + labs(title = "(c)") +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 18))) /
    (violin_sm + labs(title = "(d)") +
        theme(plot.title = element_text(face = "bold", hjust = 0, size = 18)))

# Combine all three columns with shared legend
fig_combined <- (col_a | col_b | col_c) +
    plot_layout(widths = c(1, 1, 1.6), guides = "collect") &
    theme(
        legend.position = "left",
        legend.justification = "center",
        plot.margin = margin(2, 2, 2, 2)
    )

# ===========================================================================
# 8. SAVE
# ===========================================================================
path_out <- file.path(base_dir, "output", "figures")
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(path_out, "Fig_2_AET_SM_combined.png"), fig_combined,
    width = 36, height = 30, units = "cm", dpi = 400, bg = "white"
)

cat(sprintf(
    "Done. Figure saved: %s\n",
    file.path(path_out, "Fig_2_AET_SM_combined.png")
))
