# Match streamflow stations to analysis catchments --------------------------
# Two-pass matching: strict (outlet proximity + area ratio) then relaxed
# (point-in-polygon + area constraint)

# Library calling ----------------------------------------------------------
library(sf)
library(dplyr)
library(tidyr)
library(ncdf4)
library(data.table)
library(sp)

# Path configuration -------------------------------------------------------
source("config/paths.R")

# Input paths
stations_path <- file.path(base_dir, "Stations_ValidationF.csv")
catchments_path <- gpkg_path

# Input validation ---------------------------------------------------------
if (!file.exists(stations_path)) {
    stop("Station metadata file not found: ", stations_path)
}
if (!file.exists(catchments_path)) {
    stop("Catchment polygons file not found: ", catchments_path)
}

# --------------------------------------------------------------------------- #
# Function declarations (from project conventions)                            #
# --------------------------------------------------------------------------- #
outletopen <- function(dir, outletname, nrspace = rep(NA, 5)) {
    ncbassin <- paste0(dir, "/", outletname, ".nc")
    ncb <- nc_open(ncbassin)
    name.vb <- names(ncb[["var"]])
    namev <- name.vb[1]
    name.lon <- "lon"
    name.lat <- "lat"
    if (!is.na(nrspace[1])) {
        start <- as.numeric(nrspace[c(2, 4)])
        count <- as.numeric(nrspace[c(3, 5)]) - start + 1
    } else {
        londat <- ncvar_get(ncb, name.lon)
        llo <- length(londat)
        latdat <- ncvar_get(ncb, name.lat)
        lla <- length(latdat)
        start <- c(1, 1)
        count <- c(llo, lla)
    }
    londat <- ncvar_get(ncb, name.lon, start = start[1], count = count[1])
    llo <- length(londat)
    latdat <- ncvar_get(ncb, name.lat, start = start[2], count = count[2])
    lla <- length(latdat)
    outlets <- ncvar_get(ncb, namev, start = start, count = count)
    outlets <- as.vector(outlets)
    outll <- expand.grid(londat, latdat)
    lonlatloop <- expand.grid(c(1:llo), c(1:lla))
    outll$idlo <- lonlatloop$Var1
    outll$idla <- lonlatloop$Var2
    outll <- outll[which(!is.na(outlets)), ]
    outlets <- outlets[which(!is.na(outlets))]
    outll <- data.frame(outlets, outll)
    nc_close(ncb)
    return(outll)
}

UpAopen <- function(dir, outletname, Sloc_final) {
    ncbassin <- paste0(dir, outletname)
    ncb <- nc_open(ncbassin)
    name.vb <- names(ncb[["var"]])
    namev <- name.vb[2]
    name.lon <- "lon"
    name.lat <- "lat"
    londat <- ncvar_get(ncb, name.lon)
    llo <- length(londat)
    latdat <- ncvar_get(ncb, name.lat)
    lla <- length(latdat)
    start <- c(1, 1)
    count <- c(llo, lla)
    londat <- ncvar_get(ncb, name.lon, start = start[1], count = count[1])
    llo <- length(londat)
    latdat <- ncvar_get(ncb, name.lat, start = start[2], count = count[2])
    lla <- length(latdat)
    outlets <- ncvar_get(ncb, namev, start = start, count = count)
    outlets <- as.vector(outlets) / 1000000
    outll <- expand.grid(londat, latdat)
    lonlatloop <- expand.grid(c(1:llo), c(1:lla))
    outll$upa <- outlets
    outll$idlo <- lonlatloop$Var1
    outll$idla <- lonlatloop$Var2
    outll$latlong <- paste(round(outll$Var1, 4), round(outll$Var2, 4), sep = " ")
    outfinal <- inner_join(outll, Sloc_final, by = "latlong")
    nc_close(ncb)
    return(outfinal)
}

# [1/5] Load station metadata and catchment polygons -----------------------
cat("[1/5] Loading station metadata and catchment polygons...\n")

# --- Station metadata ---
stations_raw <- read.csv(stations_path, stringsAsFactors = FALSE)

# Filter: exclude NA in Var1/Var2, exclude out-of-bounds
valid_coords <- !is.na(stations_raw$Var1) & !is.na(stations_raw$Var2)
valid_lon <- stations_raw$Var1 >= -180 & stations_raw$Var1 <= 180
valid_lat <- stations_raw$Var2 >= -90 & stations_raw$Var2 <= 90
keep <- valid_coords & valid_lon & valid_lat

n_total <- nrow(stations_raw)
n_valid <- sum(keep)
n_excluded <- n_total - n_valid

stations_sf <- st_as_sf(
    stations_raw[keep, ],
    coords = c("Var1", "Var2"),
    crs = 4326
)

cat(sprintf(
    "  Stations: %d valid, %d excluded (total %d)\n",
    n_valid, n_excluded, n_total
))

st_write(stations_sf, paste0(base_dir, "stations.shp"))
# --- Catchment polygons ---
cat("  Loading catchment polygons...\n")
catchments_gpkg <- st_read(catchments_path, quiet = TRUE)

# Reproject to EPSG:4326 if CRS differs
if (st_crs(catchments_gpkg)$epsg != 4326) {
    ids_before <- catchments_gpkg$catch_id
    catchments_gpkg <- st_transform(catchments_gpkg, crs = 4326)
    ids_after <- catchments_gpkg$catch_id
    missing_ids <- setdiff(ids_before, ids_after)
    if (length(missing_ids) > 0) {
        stop(
            "catch_id values lost during reprojection: ",
            paste(missing_ids, collapse = ", ")
        )
    }
}

cat(sprintf("  Catchments: %d polygons loaded\n", nrow(catchments_gpkg)))

# --- Outlets and upstream area (project convention) ---
cat("  Loading outlet coordinates and upstream area...\n")

outletname <- "outletsv8_hybas07_01min"
outhybas_raw <- outletopen(hydro_dir, outletname)
outhybas_raw$idlalo <- paste(outhybas_raw$idlo, outhybas_raw$idla, sep = " ")
outhybas_raw$latlong <- paste(round(outhybas_raw$Var1, 4),
    round(outhybas_raw$Var2, 4),
    sep = " "
)

# Upstream area from NetCDF
UpArea_full <- UpAopen(hydro_dir, "upArea_European_01min.nc", outhybas_raw)

# Filter to EU domain using river network mask
out1 <- outletopen(hydro_dir, "efas_rnet_100km_01min")
out1$latlong <- paste(round(out1$Var1, 4), round(out1$Var2, 4), sep = " ")
outhybas_eu <- inner_join(out1, outhybas_raw, by = "latlong")

# Match to upstream area
matcat <- match(outhybas_eu$latlong, UpArea_full$latlong)
UpArea <- UpArea_full[matcat, ]

# Column selector: get catch_ids present in LISFLOOD discharge files
header_ref <- fread(
    paste0(base_dir, "data/tss/HERA_Histo/disWin_1951_2020.csv"),
    nrows = 0, header = TRUE
)
matcol <- match(UpArea_full$outlets, as.numeric(colnames(header_ref)))
matcol <- matcol[!is.na(matcol)]
cnames <- as.character(UpArea$outlets) # character IDs of our outlets

# Match gpkg rows to our outlet list so point order is consistent
gpkg_in_outlets <- match(cnames, as.numeric(catchments_gpkg$catch_id))
catchments_plot <- catchments_gpkg[gpkg_in_outlets[!is.na(gpkg_in_outlets)], ]

# Build outlet sf object with catch_id and upstream area for matching
outlets_df <- data.frame(
    catch_id = as.character(UpArea$outlets),
    lon = UpArea$Var1.x,
    lat = UpArea$Var2.x,
    catchment_upa = UpArea$upa,
    stringsAsFactors = FALSE
)
# Remove rows with NA coordinates
outlets_df <- outlets_df[!is.na(outlets_df$lon) & !is.na(outlets_df$lat), ]

outlets_sf <- st_as_sf(outlets_df, coords = c("lon", "lat"), crs = 4326)

cat(sprintf(
    "  Outlets: %d catchment outlet points with upstream area\n",
    nrow(outlets_sf)
))


# [2/5] Performing strict matching ------------------------------------------
cat("[2/5] Performing strict matching...\n")

# Exclude stations with NA upstream area
stations_valid <- stations_sf[!is.na(stations_sf$upa), ]
cat(sprintf(
    "  Stations with valid upstream area: %d / %d\n",
    nrow(stations_valid), nrow(stations_sf)
))

# Exclude outlets (catchments) with NA upstream area
outlets_valid <- outlets_sf[!is.na(outlets_sf$catchment_upa), ]
cat(sprintf(
    "  Outlets with valid upstream area: %d / %d\n",
    nrow(outlets_valid), nrow(outlets_sf)
))

# Compute geodesic distance matrix (meters) between all valid stations and outlets
dist_matrix <- st_distance(stations_valid, outlets_valid) # units: meters
dist_matrix_km <- as.numeric(dist_matrix) / 1000
dim(dist_matrix_km) <- dim(dist_matrix)

# Build all station-outlet pair combinations with distance and upa_ratio
n_stations <- nrow(stations_valid)
n_outlets <- nrow(outlets_valid)

# Extract upstream area vectors
station_upa_vec <- stations_valid$upa
# catchment_upa_vec <- outlets_valid$catchment_upa
ordc <- match(as.numeric(catchments_plot$catch_id), outlets_sf$catch_id)

catchment_upa_vec <- catchments_plot$residual_area_km2[ordc]

# Create data frame of all pairs
pairs_list <- expand.grid(
    station_idx = seq_len(n_stations),
    outlet_idx = seq_len(n_outlets)
)

pairs_list$station_id <- stations_valid$V1[pairs_list$station_idx]
pairs_list$catch_id <- catchments_plot$catch_id[pairs_list$outlet_idx]
pairs_list$distance_km <- dist_matrix_km[
    cbind(pairs_list$station_idx, pairs_list$outlet_idx)
]
pairs_list$station_upa <- station_upa_vec[pairs_list$station_idx]
pairs_list$catchment_upa <- catchment_upa_vec[pairs_list$outlet_idx]
pairs_list$upa_ratio <- pairs_list$station_upa / pairs_list$catchment_upa

# Filter to strict match candidates: distance <= 5 km AND ratio in [0.5, 1.5]
strict_candidates <- pairs_list[
    pairs_list$distance_km <= 5 &
        pairs_list$upa_ratio >= 0.5 &
        pairs_list$upa_ratio <= 1.5,
    c(
        "station_id", "catch_id", "distance_km",
        "station_upa", "catchment_upa", "upa_ratio"
    )
]
rownames(strict_candidates) <- NULL

cat(sprintf(
    "  Strict match candidates: %d pairs (distance <= 2 km & ratio in [0.5, 1.5])\n",
    nrow(strict_candidates)
))

# # Keep only the closest station per catchment
# strict_candidates <- strict_candidates[order(strict_candidates$distance_km), ]
# strict_candidates <- strict_candidates[!duplicated(strict_candidates$station_id), ]
#
# cat(sprintf(
#   "  Strict match candidates: %d pairs (distance <= 2 km & ratio in [0.5, 1.5])\n",
#   nrow(strict_candidates)
# ))
# length(unique(strict_candidates$catch_id))
# --- Iterative one-to-one conflict resolution ---
# Algorithm:
# 1. Each station selects its best catchment (smallest distance, tie-break: smallest |upa diff|)
# 2. Each catchment with multiple stations selects best station (smallest |upa diff|, tie-break: smallest distance)
# 3. Repeat until all assignments are one-to-one
#
if (nrow(strict_candidates) == 0) {
    strict_matches <- strict_candidates[, c(
        "catch_id", "station_id", "distance_km",
        "station_upa", "catchment_upa", "upa_ratio"
    )]
    cat("  No strict match candidates to resolve.\n")
} else {
    # Add absolute upa difference column for tie-breaking
    strict_candidates$upa_diff <- abs(
        strict_candidates$station_upa - strict_candidates$catchment_upa
    )

    pool <- strict_candidates
    resolved <- FALSE
    iteration <- 0

    while (!resolved) {
        iteration <- iteration + 1

        # Step 1: Each station picks its best catchment
        # (smallest distance, tie-break: smallest |upa_diff|)
        pool <- pool[order(pool$station_id, pool$distance_km, pool$upa_diff), ]
        station_best <- pool[!duplicated(pool$station_id), ]

        # Step 2: Each catchment with multiple stations picks the best station
        # (smallest |upa_diff|, tie-break: smallest distance)
        station_best <- station_best[
            order(
                station_best$catch_id, station_best$upa_diff,
                station_best$distance_km
            ),
        ]
        catchment_best <- station_best[!duplicated(station_best$catch_id), ]

        # Check if one-to-one: no duplicate station_id in catchment_best
        if (!any(duplicated(catchment_best$station_id))) {
            resolved <- TRUE
        } else {
            # Remove losing pairs: keep only the winning assignments
            # Stations that appear multiple times keep only their best catchment
            # (smallest upa_diff among the catchment_best rows for that station)
            catchment_best <- catchment_best[
                order(
                    catchment_best$station_id, catchment_best$upa_diff,
                    catchment_best$distance_km
                ),
            ]
            winners <- catchment_best[!duplicated(catchment_best$station_id), ]

            # Reduce the pool: remove pairs involving resolved stations/catchments
            resolved_stations <- winners$station_id
            resolved_catchments <- winners$catch_id

            # Keep winners plus unresolved candidates
            pool_remaining <- pool[
                !(pool$station_id %in% resolved_stations) &
                    !(pool$catch_id %in% resolved_catchments),
            ]

            pool <- rbind(
                winners[, names(pool)],
                pool_remaining
            )

            # Safety check: if pool hasn't changed, force resolution
            if (nrow(pool) == 0) {
                resolved <- TRUE
                catchment_best <- winners
            }
        }

        # Safety limit to prevent infinite loops
        if (iteration > 100) {
            warning("Conflict resolution exceeded 100 iterations, stopping.")
            break
        }
    }

    # Final strict matches
    strict_matches <- catchment_best[, c(
        "catch_id", "station_id", "distance_km",
        "station_upa", "catchment_upa", "upa_ratio"
    )]
    rownames(strict_matches) <- NULL

    cat(sprintf(
        "  Strict matches resolved: %d one-to-one assignments (%d iterations)\n",
        nrow(strict_matches), iteration
    ))
}
#
# strict_matches=strict_candidates
cat(sprintf("  Strict matches found: %d\n", nrow(strict_matches)))


# [3/5] Performing relaxed matching -----------------------------------------
cat("[3/5] Performing relaxed matching...\n")

# Identify catchments with no strict match assignment
all_catch_ids <- as.character(catchments_plot$catch_id)
strict_catch_ids <- as.character(strict_matches$catch_id)
unmatched_ids <- setdiff(all_catch_ids, strict_catch_ids)

cat(sprintf(
    "  Unmatched catchments for relaxed matching: %d / %d\n",
    length(unmatched_ids), length(all_catch_ids)
))

if (length(unmatched_ids) == 0) {
    # No unmatched catchments — empty relaxed_matches
    relaxed_matches <- data.frame(
        catch_id = character(0),
        station_id = character(0),
        distance_km = numeric(0),
        station_upa = numeric(0),
        catchment_upa = numeric(0),
        stringsAsFactors = FALSE
    )
    cat("  No unmatched catchments to process.\n")
} else {
    # Subset to unmatched catchment polygons
    unmatched_polys <- catchments_plot[
        as.character(catchments_plot$catch_id) %in% unmatched_ids,
    ]

    # All stations are eligible for relaxed matching (Req 3.1)
    # including those already assigned via strict match to OTHER catchments
    candidates_sf <- stations_sf

    # --- Point-in-polygon: find stations within each unmatched polygon ---
    # st_intersects returns a sparse list: for each polygon, which stations are inside
    pip_result <- st_intersects(unmatched_polys, candidates_sf)

    # Build a long-form data frame of (catch_id, station_idx) pairs
    pip_pairs <- do.call(rbind, lapply(seq_along(pip_result), function(i) {
        station_indices <- pip_result[[i]]
        if (length(station_indices) == 0) {
            return(NULL)
        } else {
            data.frame(
                poly_idx = i,
                catch_id = as.character(unmatched_polys$catch_id[i]),
                station_idx = station_indices,
                stringsAsFactors = FALSE
            )
        }
    }))
    if (is.null(pip_pairs) || nrow(pip_pairs) == 0) {
        relaxed_matches <- data.frame(
            catch_id = character(0),
            station_id = character(0),
            distance_km = numeric(0),
            station_upa = numeric(0),
            catchment_upa = numeric(0),
            stringsAsFactors = FALSE
        )
        cat("  No stations found within any unmatched catchment polygon.\n")
    } else {
        # Attach station attributes
        pip_pairs$station_id <- as.character(candidates_sf$V1[pip_pairs$station_idx])
        pip_pairs$station_upa <- candidates_sf$upa[pip_pairs$station_idx]
        pip_pairs$station_Rlen <- candidates_sf$Rlen[pip_pairs$station_idx]

        # Attach catchment upstream area from outlets_sf
        outlets_upa_lookup <- data.frame(
            catch_id = as.character(catchments_plot$catch_id),
            catchment_upa = catchments_plot$residual_area_km2,
            stringsAsFactors = FALSE
        )
        pip_pairs <- merge(pip_pairs, outlets_upa_lookup, by = "catch_id", all.x = TRUE)

        # --- Filter: station_upa <= catchment_upa (Req 3.2) ---
        pip_pairs <- pip_pairs[
            !is.na(pip_pairs$station_upa) &
                !is.na(pip_pairs$catchment_upa) &
                pip_pairs$station_upa <= pip_pairs$catchment_upa,
        ]

        cat(sprintf(
            "  Candidate station-catchment pairs after UPA filter: %d\n",
            nrow(pip_pairs)
        ))

        if (nrow(pip_pairs) == 0) {
            relaxed_matches <- data.frame(
                catch_id = character(0),
                station_id = character(0),
                distance_km = numeric(0),
                station_upa = numeric(0),
                catchment_upa = numeric(0),
                stringsAsFactors = FALSE
            )
            cat("  No candidates remain after upstream area filtering.\n")
        } else {
            # --- Multi-polygon conflict resolution (Req 3.4) ---
            # If a station falls in multiple unmatched catchments, assign to
            # the catchment with the smallest polygon area first.
            # Compute polygon areas in km2
            poly_areas_m2 <- st_area(unmatched_polys) # units: m^2
            poly_areas_km2 <- as.numeric(poly_areas_m2) / 1e6
            area_lookup <- data.frame(
                catch_id = as.character(unmatched_polys$catch_id),
                poly_area_km2 = poly_areas_km2,
                stringsAsFactors = FALSE
            )
            pip_pairs <- merge(pip_pairs, area_lookup, by = "catch_id", all.x = TRUE)

            # For stations appearing in multiple catchments, assign to smallest area
            # Sort by station_id then poly_area_km2 (ascending)
            pip_pairs <- pip_pairs[order(pip_pairs$station_id, pip_pairs$poly_area_km2), ]

            # Identify stations in multiple catchments
            dup_stations <- pip_pairs$station_id[duplicated(pip_pairs$station_id)]

            if (length(dup_stations) > 0) {
                # For each duplicated station, keep only the smallest-area catchment
                # The station is "consumed" by the smallest-area catchment
                conflict_stations <- unique(dup_stations)
                cat(sprintf(
                    "  Multi-polygon conflicts: %d stations in multiple catchments\n",
                    length(conflict_stations)
                ))

                # Keep first occurrence per station (smallest area due to sorting)
                pip_pairs_resolved <- pip_pairs[!duplicated(pip_pairs$station_id), ]

                # Also keep pairs where the station is NOT in conflict
                # (already handled by !duplicated keeping the first)
            } else {
                pip_pairs_resolved <- pip_pairs
            }

            # --- Best station selection per catchment (Req 3.3) ---
            # For each catchment with multiple valid candidates:
            # Select station with longest Rlen; tie-break by distance to centroid

            # Compute centroids of unmatched polygons
            centroids <- st_centroid(unmatched_polys)

            # Compute geodesic distance from each candidate station to its
            # catchment centroid
            pip_pairs_resolved$dist_to_centroid_km <- NA_real_

            for (i in seq_len(nrow(pip_pairs_resolved))) {
                cid <- pip_pairs_resolved$catch_id[i]
                sidx <- pip_pairs_resolved$station_idx[i]

                # Get centroid for this catchment
                centroid_row <- which(
                    as.character(unmatched_polys$catch_id) == cid
                )
                if (length(centroid_row) == 0) next

                # Compute geodesic distance between station and centroid
                d <- st_distance(
                    candidates_sf[sidx, ],
                    centroids[centroid_row, ]
                )
                pip_pairs_resolved$dist_to_centroid_km[i] <- as.numeric(d) / 1000
            }


            pip_pairs_resolved <- pip_pairs_resolved[
                order(
                    pip_pairs_resolved$catch_id,
                    -pip_pairs_resolved$station_Rlen,
                    pip_pairs_resolved$dist_to_centroid_km
                ),
            ]
            best_per_catchment <- pip_pairs_resolved

            # Build relaxed_matches output
            relaxed_matches <- data.frame(
                catch_id = best_per_catchment$catch_id,
                station_id = best_per_catchment$station_id,
                distance_km = round(best_per_catchment$dist_to_centroid_km, 2),
                station_upa = best_per_catchment$station_upa,
                catchment_upa = best_per_catchment$catchment_upa,
                stringsAsFactors = FALSE
            )
            rownames(relaxed_matches) <- NULL

            cat(sprintf("  Relaxed matches found: %d\n", nrow(relaxed_matches)))
        }
    }
}

cat(sprintf(
    "  Total matches so far: %d strict + %d relaxed = %d\n",
    nrow(strict_matches), nrow(relaxed_matches),
    nrow(strict_matches) + nrow(relaxed_matches)
))

cat(sprintf(
    "  Total individual catchment matched with stations: %d strict + %d relaxed = %d\n",
    length(unique(strict_matches$catch_id)), length(unique(relaxed_matches$catch_id)),
    length(unique(strict_matches$catch_id)) + length(unique(relaxed_matches$catch_id))
))
# [4/5] Exporting matching table -------------------------------------------
cat("[4/5] Exporting matching table...\n")

# Combine strict and relaxed matches into a single data frame
n_strict <- length(unique(strict_matches$catch_id))
n_relaxed <- length(unique(relaxed_matches$catch_id))
n_total_matches <- n_strict + n_relaxed

if (n_total_matches == 0) {
    # Zero-match case: write header-only CSV with warning
    cat("  WARNING: No station-to-catchment matches were produced.\n")

    output_df <- data.frame(
        catch_id = character(0),
        station_id = character(0),
        match_type = character(0),
        distance_km = numeric(0),
        station_upa = numeric(0),
        catchment_upa = numeric(0),
        upa_ratio = numeric(0),
        station_name = character(0),
        river_name = character(0),
        record_length = integer(0),
        stringsAsFactors = FALSE
    )
} else {
    # --- Build strict portion ---
    if (n_strict > 0) {
        strict_df <- data.frame(
            catch_id = as.character(strict_matches$catch_id),
            station_id = as.character(strict_matches$station_id),
            match_type = "strict",
            distance_km = round(strict_matches$distance_km, 2),
            station_upa = strict_matches$station_upa,
            catchment_upa = strict_matches$catchment_upa,
            upa_ratio = round(strict_matches$upa_ratio, 3),
            stringsAsFactors = FALSE
        )
    } else {
        strict_df <- data.frame(
            catch_id = character(0),
            station_id = character(0),
            match_type = character(0),
            distance_km = numeric(0),
            station_upa = numeric(0),
            catchment_upa = numeric(0),
            upa_ratio = numeric(0),
            stringsAsFactors = FALSE
        )
    }

    # --- Build relaxed portion ---
    if (n_relaxed > 0) {
        relaxed_df <- data.frame(
            catch_id = as.character(relaxed_matches$catch_id),
            station_id = as.character(relaxed_matches$station_id),
            match_type = "relaxed",
            distance_km = round(relaxed_matches$distance_km, 2),
            station_upa = relaxed_matches$station_upa,
            catchment_upa = relaxed_matches$catchment_upa,
            upa_ratio = NA_real_,
            stringsAsFactors = FALSE
        )
    } else {
        relaxed_df <- data.frame(
            catch_id = character(0),
            station_id = character(0),
            match_type = character(0),
            distance_km = numeric(0),
            station_upa = numeric(0),
            catchment_upa = numeric(0),
            upa_ratio = numeric(0),
            stringsAsFactors = FALSE
        )
    }

    # Combine strict and relaxed
    combined_df <- rbind(strict_df, relaxed_df)

    # --- Attach station metadata: station_name, river_name, record_length ---
    # Station metadata columns available: V1, Var1, Var2, upa, csource, Rlen
    # No explicit station_name or river_name columns exist in the metadata,
    # so we set them to NA. record_length comes from Rlen.
    metadata_lookup <- data.frame(
        station_id = as.character(stations_raw$V1),
        station_name = NA_character_,
        river_name = NA_character_,
        record_length = as.integer(stations_raw$Rlen),
        stringsAsFactors = FALSE
    )
    # Remove duplicate station IDs (keep first occurrence)
    metadata_lookup <- metadata_lookup[!duplicated(metadata_lookup$station_id), ]

    output_df <- merge(combined_df, metadata_lookup,
        by = "station_id", all.x = TRUE
    )

    # Reorder columns to match required format
    output_df <- output_df[, c(
        "catch_id", "station_id", "match_type", "distance_km",
        "station_upa", "catchment_upa", "upa_ratio",
        "station_name", "river_name", "record_length"
    )]

    # Sort by catch_id ascending
    output_df <- output_df[order(output_df$catch_id), ]
    rownames(output_df) <- NULL
}

# Create output directory if it does not exist
output_dir <- file.path(base_dir, "output")
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

catmatch <- catchments_gpkg[!is.na(match(catchments_gpkg$catch_id, output_df$catch_id)), ]
sum(catmatch$residual_area_km2) / sum(catchments_gpkg$residual_area_km2)
# Write CSV (UTF-8, header row)
output_path <- file.path(output_dir, "station_catchment_matches.csv")
write.csv(output_df,
    file = output_path, row.names = FALSE,
    fileEncoding = "UTF-8"
)

cat(sprintf("  Output written to: %s\n", output_path))

# Print summary
n_unmatched <- length(all_catch_ids) - n_total_matches
cat(sprintf(
    "  Summary: %d catchments matched (%d strict + %d relaxed), %d unmatched\n",
    n_total_matches, n_strict, n_relaxed, n_unmatched
))
