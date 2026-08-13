# =============================================================================
# plot_catchment_nesting_map.R
# =============================================================================
# Produces a map of all catchments showing headwaters vs nested outlets,
# with a histogram inset of the residual area distribution.
# Output: plots/Figure1_catchments_v2.png
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(sf)
library(data.table)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(cowplot)

# --- Paths --------------------------------------------------------------------
base_dir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
gpkg_path <- file.path(base_dir, "data", "catchments_analysis_final_v3.gpkg")
plot_dir <- file.path(base_dir, "output", "figures")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load catchments ----------------------------------------------------------
message("Reading catchments GeoPackage ...")
catchments_gpkg <- st_read(gpkg_path, quiet = TRUE)

# --- Classify headwaters vs nested --------------------------------------------
catchments_gpkg$n_nested <- sapply(
    catchments_gpkg$immediate_nested_ids,
    function(x) {
        if (is.na(x) || trimws(x) == "NA" || trimws(x) == "") {
            return(0L)
        }
        length(trimws(strsplit(x, ",")[[1]]))
    }
)
catchments_gpkg$is_headwater <- catchments_gpkg$n_nested == 0

message(sprintf(
    "  %d headwaters  |  %d nested outlets",
    sum(catchments_gpkg$is_headwater),
    sum(!catchments_gpkg$is_headwater)
))

# --- Prepare spatial objects in EPSG:3035 -------------------------------------
# Centroids for bounding box
plot_pts <- catchments_gpkg |>
    st_centroid() |>
    st_transform(crs = 3035)
coords <- st_coordinates(plot_pts)

# Split polygons
hw_poly <- catchments_gpkg[catchments_gpkg$is_headwater, ]
nested_poly <- catchments_gpkg[!catchments_gpkg$is_headwater, ]

# World basemap
world <- ne_countries(scale = "medium", returnclass = "sf")
basemap <- st_transform(world, crs = 3035)

# --- Main map -----------------------------------------------------------------
p_map <- ggplot(basemap) +
    geom_sf(fill = "grey95", color = "grey70", linewidth = 0.2) +
    geom_sf(
        data = hw_poly,
        fill = "lightblue", color = "grey30",
        stroke = 0.2, alpha = 0.7
    ) +
    geom_sf(
        data = nested_poly,
        aes(fill = nesting_level),
        color = "grey20", stroke = 0.2, alpha = 0.9
    ) +
    scale_fill_gradientn(
        colors = hcl.colors(9, palette = "BluYl", rev = FALSE),
        name = "No. \nnested catchments",
        limits = c(0, 50), oob = scales::squish,
        breaks = scales::pretty_breaks(n = 5)
    ) +
    coord_sf(
        xlim = c(min(coords[, 1]) - 1e5, max(coords[, 1]) + 1e5),
        ylim = c(min(coords[, 2]) - 1e5, max(coords[, 2]) + 1e5)
    ) +
    labs(
        title = "Catchments",
        subtitle = sprintf(
            "%d headwaters (lightblue)  |  %d nested outlets (coloured)",
            sum(catchments_gpkg$is_headwater),
            sum(!catchments_gpkg$is_headwater)
        ),
        x = "Longitude", y = "Latitude"
    ) +
    guides(
        fill = guide_colourbar(barwidth = 0.5, barheight = 15, reverse = FALSE)
    ) +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, colour = "grey40"),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour = "black"),
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.position = "right",
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(0.8, "cm")
    )

# --- Histogram inset ----------------------------------------------------------
p_hist <- ggplot(catchments_gpkg, aes(x = residual_area_km2)) +
    geom_histogram(
        color = "steelblue", fill = "slategray1",
        bins = 20, alpha = 0.9, lwd = 0.2
    ) +
    scale_y_continuous(
        breaks = seq(0, 600, by = 100),
        name = "Number of outlets"
    ) +
    scale_x_log10(
        name = expression(paste("Residual area (", km^2, ")")),
        breaks = c(1, 10, 100, 1000, 10000, 100000),
        minor_breaks = log10_minor_break(),
        labels = c("1", "10", "100", "1 000", "10 000", "100 000")
    ) +
    theme(
        axis.title = element_text(size = 4),
        axis.ticks.length = unit(0.5, "pt"),
        axis.text = element_text(size = 3),
        axis.ticks = element_line(linewidth = 0.1),
        axis.title.x = element_text(margin = margin(t = 0.2)),
        axis.title.y = element_text(margin = margin(r = 0.2)),
        axis.text.x = element_text(margin = margin(t = 0.2)),
        axis.text.y = element_text(margin = margin(r = 0.2)),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_rect(
            linetype = "solid", fill = NA,
            colour = "black", linewidth = 0.1
        ),
        panel.grid.major = element_line(colour = "grey80", linewidth = 0.1),
        panel.grid.minor.x = element_line(
            colour = "grey90", linetype = "dashed", linewidth = 0.1
        ),
        plot.background = element_rect(
            fill = "white", colour = "grey50", linewidth = 0.4
        ),
        plot.margin = margin(2, 3, 1, 1)
    )

# --- Combine: map + histogram inset ------------------------------------------
combined <- ggdraw(p_map) +
    draw_plot(p_hist,
        x = 0.55,
        y = 0.72,
        width = 0.18,
        height = 0.18
    )

ggsave(file.path(plot_dir, "Figure1_catchments_v2.png"),
    combined,
    width = 22, height = 16, units = "cm", dpi = 500
)

message("Saved: plots/Figure1_catchments_v2.png")
