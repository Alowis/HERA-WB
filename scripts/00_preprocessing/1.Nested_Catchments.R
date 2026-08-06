# Install once --------------------------------------------------------------
# install.packages(c(
#   "whitebox"   # wrapper for WhiteboxTools (fast watershed)
# ))
# 
# whitebox::install_whitebox()
# Load ----------------------------------------------------------------------
library(terra)     # raster work
library(sf)        # vector work
library(dplyr)     # data manipulation
library(igraph)    # graph utilities
library(whitebox)  # WhiteboxTools (optional but recommended)
library(units)
library(tidyr)
library(dplyr)
library(rstudioapi)

# Set working directory to the script's folder
setwd(dirname(getActiveDocumentContext()$path))
WD <- getwd()
if (!is.null(WD)) setwd(WD)
source("functions_regime.R")

#setwd("D:/tilloal/Documents/LFRuns_utils")
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
outletname="outletsv8_hybas07_01min"
outhybas=outletopen(hydroDir,outletname)
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")

## 2‑a  LDD raster (NetCDF)  -------------------------------------------------
nc_path <- "D:/tilloal/Documents/06_Floodrivers/mapscal/ldd_European_01min.nc"                 # <-- change to your file
ldd_raw <- rast(nc_path)   # replace "LDD" with the exact variable name


#Hybas07
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=inner_join(Catamere07,outhybas,by= c("llcoord"="latlong"))
cst7=st_transform(Catf7,  crs=3035)

#Plot parameters
palet=c(hcl.colors(9, palette = "BuPu", alpha = NULL, rev = TRUE, fixup = TRUE))
outll=outletopen(hydroDir,outletname)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
basemap=w2
palet=c(hcl.colors(9, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE))

###load UpArea -----
#load upstream area
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data/")
outletname="upArea_European_01min.nc"
outhybas$idlalo=paste(outhybas$idlo, outhybas$idla, sep=" ")
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4), sep=" ")
UpArea=UpAopen(hydroDir,outletname,outhybas)
head(UpArea)
#remove Q obs that are outside
matcat=match(outhybas$latlong,UpArea$latlong)
UpArea=UpArea[matcat,]

UpArea_l=UpArea[,c(1,2,3,6)]
CatUpA=inner_join(Catf7,UpArea_l,by=c("llcoord"="latlong"))


#load matched hybas catchments as well
hybas_path <- "D:/tilloal/Documents/01_Projects/RegimeShifts/data/HYBAS07_pixelized.shp"
hybas_pix <- st_read(hybas_path, quiet = TRUE) %>%
  st_make_valid() %>%   
  st_transform(crs(ldd_raw)) 

hybas_pix <- hybas_pix %>%
  mutate(area_km2 = as.numeric(st_area(.)) / 1e6)

# 3. Dissolve by region ID and sum areas
regions_dissolved <- hybas_pix %>%
  group_by(PFAF_ID) %>%               # replace 'region_id' with your ID column name
  summarise(
    total_area = sum(area_km2),    # sum of original polygon areas
    geometry = st_union(geometry)       # merge geometries
  )

st_geometry(Catf7)=NULL
Catf7t=inner_join(Catf7,regions_dissolved,by= c("PFAF_ID"))

#keep only rivers in EU domain
out1=outletopen(hydroDir,"efas_rnet_100km_01min")
out1$latlong=paste(round(out1$Var1,4),round(out1$Var2,4),sep=" ")
outhybas=inner_join(out1,outhybas, by="latlong")

## 2‑b  Outlet points and shapefile ---------------------------------------------------------
outlets <- st_as_sf(UpArea, coords = c("Var1.x", "Var2.x"), crs = 4326)
# ppl <- st_transform(ppl, crs = 3035)          # make CRS identical to the raster

allpoly_path <- "D:/tilloal/Documents/01_Projects/RegimeShifts/All_catchment_raw.shp"
catch_tot <- st_read(allpoly_path, quiet = TRUE) %>%
  st_transform(crs(ldd_raw)) 



## 2‑c other catchment polygons -------------------------------------------
poly_path <- "D:/tilloal/Documents/01_Projects/RegimeShifts/data/catchments_others2.shp"
catch_pol <- st_read(poly_path, quiet = TRUE) %>%
  st_transform(crs(ldd_raw))               # same CRS as everything else


catch_pol <- st_read(poly_path, quiet = TRUE) %>%   # quiet = TRUE hides the progress messages
  st_make_valid() %>%                               # fixes self‑intersections, etc.
  st_transform(crs(ldd_raw))                          # any CRS you like; keep it consistent later


#load matched hybas catchments as well
hybas_path <- "D:/tilloal/Documents/01_Projects/RegimeShifts/data/catchments_hybas2.shp"
catch_hy <- st_read(hybas_path, quiet = TRUE) %>%
  st_make_valid() %>%   
  st_transform(crs(ldd_raw)) 

#check how HYBAS catchments compare with the catchments obtained from upstream area only.

hyvsup=match(catch_hy$otlts_y_y,as.numeric(catch_tot$Name))

catch_tot <- catch_tot %>%
  mutate(area_km2 = as.numeric(st_area(.)) / 1e6)

plot(catch_tot$area_km2[hyvsup],catch_hy$SUB_ARE)
#catch_hy<-catch_hy[,c(30,36)]
colnames(catch_hy)[1]=colnames(catch_pol)[1]
catch_pol$group="other"
catch_hy$group="hybas"
#all catchments are here for the nested analysis check
catch_ph=rbind(catch_pol,catch_hy)


# -----------------------------------------------------------------
# 2‑c  Make sure we have a column called `catch_id`
# -----------------------------------------------------------------
# If the original file already has a field that looks like an ID (e.g. "ID",
# "CatchID", "PolyID", …) we rename it to `catch_id`.  If not, we create one.
catchments=catch_tot
#catchments$Name=catchments$otlts_y_x
#catchments<-catchments[-which(is.na(catchments$Name)),]

#catchments<-catch_pol
possible_id_names <- names(catchments)[
  grepl("^id$|^Name|^poly|^feat", names(catchments), ignore.case = TRUE)
]

if (length(possible_id_names) > 0) {
  # Keep the first column that matches the pattern and rename it
  catchments <- catchments %>% dplyr::rename(catch_id = !!sym(possible_id_names[1]))
} else {
  # No ID column – create a simple 1,2,3,… identifier
  catchments <- catchments %>% mutate(catch_id = row_number())
}

# Verify that the column really exists
stopifnot("catch_id" %in% names(catchments))
# --------------------------------------------------------------
# 3. Find direct children (loop version) ----------------------
# --------------------------------------------------------------

# Fix the polygons before the spatial join
catchments <- st_make_valid(catchments)

# If you still get errors, try a tiny buffer (the "zero-buffer" trick)
# This often snaps those duplicate vertices back into place
#catchments <- st_buffer(catchments, dist = 0)
# contain_list <- st_contains(catchments, catchments, sparse = TRUE)

nest_matrix <- st_contains(catchments, catchments, sparse = FALSE)


# Disable S2 for speed and to avoid the "degenerate edge" errors
#sf_use_s2(FALSE)

# 1. Calculate centroids for all catchments
centroids <- st_centroid(st_geometry(catchments))

#instead of centroids i use outlets coordinates
outlets=st_make_valid(CatUpA)%>%   
  st_transform(crs(ldd_raw)) 

# 3. Spatial Join: Find which catchment index each outlet falls into
# This assumes one outlet per catchment
matched_data <- match(as.numeric(catchments$catch_id),outlets$outlets.y)
outlets<-outlets[matched_data,]

# 2. Check which centroids are inside which polygons
# This is an order of magnitude faster than st_contains(poly, poly)
nest_matrix1 <- st_intersects(catchments, outlets, sparse = FALSE)
diag(nest_matrix1) <- FALSE

nest_matrix2 <- st_intersects(catchments, outlets)
# 2. Initialize a clean nesting matrix
# Calculate areas in km2 once
catchments <- catchments %>%
  mutate(area_km2 = as.numeric(st_area(.)) / 1e6)



# 2. Extract geometry list once (faster than repeatedly calling st_geometry)
geom_list <- st_geometry(catchments)

n <- nrow(catchments)
refined_nest_matrix <- matrix(FALSE, nrow = n, ncol = n)
cat("Refining nesting matrix...\n")

for (i in 1:n) {
  if (i %% 50 == 0) cat("Processing catchment:", i, "/", n, "\n")
  candidates <- nest_matrix2[[i]]
  candidates <- candidates[candidates != i]
  
  if (length(candidates) == 0) next
  
  # Optimization 1: Only keep candidates smaller than the current parent (i)
  # This cuts the number of st_intersection calls roughly in half.
  valid_children <- candidates[catchments$area_km2[candidates] < catchments$area_km2[i]]
  
  if (length(valid_children) == 0) next
  
  # Optimization 2: Vectorized intersection
  # We intersect the parent geometry with ALL potential children at once
  parent_geom <- geom_list[i]
  child_geoms <- geom_list[valid_children]
  parent_buffered <- st_buffer(parent_geom, 1)
  # st_intersection here returns the parts of children inside the parent
  suppressMessages(sf_use_s2(FALSE)) # switch to GEOS — handles degenerate edges gracefully
  
  inter_areas <- suppressMessages(vapply(valid_children, function(j) {
    overlap <- tryCatch(
      st_intersection(
        st_make_valid(parent_geom),
        st_make_valid(geom_list[j])
      ),
      error   = function(e) NULL,
      warning = function(w) suppressWarnings(
        st_intersection(
          st_make_valid(parent_geom),
          st_make_valid(geom_list[j])
        )
      )
    )
    if (is.null(overlap) || length(overlap) == 0 || st_is_empty(overlap)) return(0)
    as.numeric(st_area(overlap)) / 1e6
  }, numeric(1))
  )
  suppressMessages(sf_use_s2(TRUE))  # restore s2 for everything else
  # intersections <- st_intersection(child_geoms, parent_buffered)
  # intersections <- st_make_valid(intersections)
  # # Calculate areas for all intersections at once
  # inter_areas <- as.numeric(st_area(intersections)) / 1e6
  
  # Compare against child areas
  child_areas <- catchments$area_km2[valid_children]
  overlap_ratios <- inter_areas / child_areas
  #plot(overlap_ratios)
  # Update the matrix for all j that pass the 90% threshold
  passing_indices <- valid_children[overlap_ratios > 0.90]
  if (length(passing_indices) > 0) {
    refined_nest_matrix[i, passing_indices] <- TRUE
  }
}

# Find indices where the relationship is TRUE
relations <- which(refined_nest_matrix, arr.ind = TRUE)

# Create a mapping table
nesting_df <- data.frame(
  parent_id = catchments$catch_id[relations[,1]], # Change 'ID' to your column name
  child_id  = catchments$catch_id[relations[,2]]
)

# # Find indices where the relationship is TRUE
# relations2 <- which(nest_matrix2, arr.ind = TRUE)
# 
# # Create a mapping table
# nesting_df2 <- data.frame(
#   parent_id = catchments$catch_id[relations2[,1]], # Change 'ID' to your column name
#   child_id  = catchments$catch_id[relations2[,2]]
# )

# Count how many polygons contain each specific polygon
catchments$nesting_level <- rowSums(refined_nest_matrix)


id_a <- which(catchments$catch_id == 292302)
id_b <- which(catchments$catch_id == 290666)

# Check what the matrix currently says
cat("292302 contains 290666:", refined_nest_matrix[id_a, id_b], "\n")
cat("290666 contains 292302:", refined_nest_matrix[id_b, id_a], "\n")

# Check if 290666's outlet falls inside 292302's polygon
outlet_b  <- outlets[outlets$outlets.y == 290666, ]
polygon_a <- catchments[id_a, ]

cat("Outlet 290666 intersects polygon 292302:",
    st_intersects(polygon_a, outlet_b, sparse = FALSE)[1, 1], "\n")

# Check the true geometric overlap between the two polygons
overlap <- st_intersection(st_geometry(catchments[id_a, ]),
                           st_geometry(catchments[id_b, ]))
overlap_area   <- as.numeric(st_area(overlap)) / 1e6
pct_of_290666  <- overlap_area / catchments$area_km2[id_b] * 100
cat(sprintf("Geometric overlap: %.1f km2 (%.1f%% of 290666)\n",
            overlap_area, pct_of_290666), "\n")

cat("Source of 292302:", catchments$group[id_a], "\n")
cat("Source of 290666:", catchments$group[id_b], "\n")
cat("Area 292302:",      catchments$area_km2[id_a], "km2\n")
cat("Area 290666:",      catchments$area_km2[id_b], "km2\n")
cat("UpArea 292302:",    UpArea$upa[match(292302, UpArea$outlets)], "km2\n")
cat("UpArea 290666:",    UpArea$upa[match(290666, UpArea$outlets)], "km2\n")

# Was 290666 even a candidate for 292302 in nest_matrix2?
id_b %in% nest_matrix2[[id_a]]

cat("matched_data for 290666:", matched_data[id_b], "\n")
cat("outlets row for 290666:", outlets$outlets.x[id_b], "\n")

cat("area_km2 292302:", catchments$area_km2[id_a], "\n")
cat("area_km2 290666:", catchments$area_km2[id_b], "\n")
cat("290666 passes area filter:", catchments$area_km2[id_b] < catchments$area_km2[id_a], "\n")

# Reproduce exactly what the loop does for i = id_a
candidates    <- nest_matrix2[[id_a]]
candidates    <- candidates[candidates != id_a]
valid_children <- candidates[catchments$area_km2[candidates] < catchments$area_km2[id_a]]

cat("id_b in candidates:     ", id_b %in% candidates, "\n")
cat("id_b in valid_children: ", id_b %in% valid_children, "\n")

# Now reproduce the intersection step
parent_geom   <- geom_list[id_a]
child_geoms   <- geom_list[valid_children]
parent_buffered <- st_buffer(parent_geom, 1)




plot_nested_branch <- function(sf_data, matrix, row_index) {
  # Find which columns (children) are TRUE for this row (parent)
  child_indices <- which(matrix[row_index, ])
  
  # Create a subset for plotting
  parent_geom <- sf_data[row_index, ]
  children_geom <- sf_data[child_indices, ]
  
  ggplot() +
    geom_sf(data = parent_geom, fill = "darkblue", color = "blue", size = 1) +
    geom_sf(data = children_geom, fill = "white", alpha = 0.5, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Catchment ID:", sf_data$catch_id[row_index]),
         subtitle = paste(length(child_indices), "Sub-catchments nested inside"))
}

# Example: Plot the first catchment in your file
plot_nested_branch(sf_data=catchments, matrix=refined_nest_matrix, row_index=197)
# Level 0 = The outermost "terminal" catchments (pouring into the ocean/sink)
# Higher Level = Smaller sub-catchments deep upstream

# Find the 'Roots' (Outermost catchments)
# They are not contained by anyone else (column sum is 0)
roots <- which(colSums(refined_nest_matrix) == 0)

# Find the 'Leaves' (Smallest headwaters)
# They contain no one else (row sum is 0)
leaves <- which(rowSums(refined_nest_matrix) == 0)

# Add this to your data for easy mapping
catchments$type <- "Mid-stream"
catchments$type[roots] <- "Root (Downstream)"
catchments$type[leaves] <- "Leaf (Headwater)"

# catchments2<-catchments[-which(is.na(catchments$catch_id)),]
# # Visualizing the types
# visu<-ggplot(catchments2) +
#   geom_sf(aes(fill = type)) +
#   scale_fill_viridis_d(alpha=0.2) +
#   theme_void()
# 
# ggsave(visu,file="D:/tilloal/Documents/01_Projects/RegimeShifts/plots/visucat.jpg", width=20, height=24, units=c("cm"),dpi=500)

# library(igraph)
# 
# # Create the graph object
# network <- graph_from_data_frame(nesting_df[c(1:12),], directed = TRUE)
# 
# # Plot the hierarchy
# plot(network, 
#      layout = layout_as_tree(network, mode = "out"),
#      vertex.size = 5,
#      vertex.label.cex = 0.7,
#      edge.arrow.size = 0.5,
#      main = "Catchment Nesting Hierarchy")



# 1. Function to find immediate children
# get_immediate_children <- function(parent_idx, nest_matrix) {
#   # All descendants
#   all_children <- which(nest_matrix[parent_idx, ])
#   
#   # A child is 'immediate' if it is not contained by any OTHER child
#   if (length(all_children)>1){
#     immediate <- all_children[colSums(nest_matrix[all_children, all_children]) == 0]
#   }else{
#     immediate <- all_children[sum(nest_matrix[all_children, all_children]) == 0]
#   }
#   return(immediate)
# }
# 
# # Pick a parent catchment index (e.g., the most downstream one)
# parent_idx <- 328
# child_indices <- get_immediate_children(parent_idx, nest_matrix)
# 
# # Get the geometries
# parent_geom <- catchments[parent_idx, ]
# children_geom <- catchments[child_indices, ]
# 
# # 2. Combine all children into one shape (Union)
# children_union <- st_union(children_geom)
# 
# # 1. Cast the children to a single geometry first
# children_geom_single <- st_union(children_geom)
# 
# # 2. Perform the difference
# # Using st_geometry() ensures we are working with the shapes, not the data frame
# diff_geom <- st_difference(st_geometry(parent_geom), children_geom_single)
# 
# 
# # 1. Isolate the geometries
# parent_shape <- st_geometry(parent_geom)
# children_union <- st_union(st_geometry(children_geom))
# 
# # 2. Perform the difference
# # If children cover the parent 100%, this returns an empty geometry
# diff_geom <- st_difference(parent_shape, children_union)
# 
# # 3. Clean up the result
# # st_as_sf handles the "nrow == length" check automatically by creating a new table
# inter_catchment <- st_as_sf(data.frame(ID = parent_geom$catch_id), geometry = diff_geom)
# 
# # 4. Optional: If the result is a GEOMETRYCOLLECTION (points + lines + polygons)
# # we usually only want the Polygons
# inter_catchment <- st_collection_extract(inter_catchment, "POLYGON")
# 
# ggplot() +
#   # The original full parent (faded)
#   geom_sf(data = parent_geom, fill = "grey90", color = "black", linetype = "dotted") +
#   # The sub-catchments (the holes)
#   geom_sf(data = children_geom, fill = "white", color = "red") +
#   # The Inter-catchment area (the actual contributing land)
#   geom_sf(data = inter_catchment, fill = "springgreen4", alpha = 0.7) +
#   theme_minimal() +
#   labs(title = "Inter-catchment (Residual) Area",
#        subtitle = "Green represents area draining directly to the reach between stations")
# 

# # catchments<-catch_pol
# # 1. Setup an empty list to store results
# residual_catchments_list <- list()
# 
# # 2. Start the loop
# for (i in 1:nrow(catchments)) {
#   print(i)
#   parent_row <- catchments[i, ]
#   
#   # Get immediate children indices using the function we defined earlier
#   child_idx <- get_immediate_children(i, nest_matrix)
#   
#   if (length(child_idx) == 0) {
#     # It's a headwater: the residual area is the whole catchment
#     residual_catchments_list[[i]] <- parent_row
#     
#   } else {
#     # It's a downstream catchment: we need to subtract children
#     parent_geom <- st_geometry(parent_row)
#     children_union <- st_union(st_geometry(catchments[child_idx, ]))
#     
#     # Use st_make_valid and st_buffer(0) to fix any topological ghosts
#     diff_geom <- st_difference(st_make_valid(parent_geom), 
#                                st_make_valid(children_union))
#     
#     # Check if anything is left after subtraction
#     if (length(diff_geom) == 0 || st_is_empty(diff_geom)) {
#       message(paste("Catchment", i, "is fully covered by children. Skipping."))
#       next
#     }
#     
#     # Ensure we only keep POLYGONS (removes lines/points from edges)
#     diff_geom <- st_collection_extract(diff_geom, "POLYGON") %>% st_union()
#     
#     # Create the new row by replacing the geometry in the parent row
#     # this preserves all your original attributes (ID, name, etc.)
#     new_row <- parent_row
#     st_geometry(new_row) <- diff_geom
#     
#     residual_catchments_list[[i]] <- new_row
#   }
# }
# 
# # 3. Combine everything back into a single SF object
# # Filter out NULLs in case some catchments were fully skipped
# final_residuals <- do.call(rbind, residual_catchments_list[!sapply(residual_catchments_list, is.null)])
# 
# # 4. Recalculate the actual area of these new pieces
# final_residuals$inter_area_m2 <- as.numeric(st_area(final_residuals))

# 1. Ensure we have the total area first
#catchments<-catchments[-which(is.na(catchments$catch_id)),]
catchments <- catchments %>%
  mutate(total_area_km2 = as.numeric(st_area(.)) / 1e6)

residual_results <- list()

for (i in 1:nrow(catchments)) {
  print(i)
  #i=304
  #i=199
  # A. Identify ALL nested children (the whole family tree upstream)
  all_child_indices <- which(refined_nest_matrix[i, ])
  all_child_ids <- paste(catchments$catch_id[all_child_indices], collapse = ", ")
  
  # B. Identify only IMMEDIATE children (for the subtraction logic)
  # A child is immediate if it is NOT contained by any other child in the set
  immediate_child_indices <- all_child_indices[colSums(refined_nest_matrix[all_child_indices, all_child_indices, drop=FALSE]) == 0]
 
  # A child is 'immediate' if it is not contained by any OTHER child
  if (length(all_child_indices)>1){
    immediate_child_indices <- all_child_indices[colSums(refined_nest_matrix[all_child_indices, all_child_indices]) == 0]
   }else{
    immediate_child_indices <- all_child_indices[sum(refined_nest_matrix[all_child_indices, all_child_indices]) == 0]
  }
  immediate_child_ids <- paste(catchments$catch_id[immediate_child_indices], collapse = ", ")
  parent_row <- catchments[i, ]
  
  if (length(immediate_child_indices) == 0) {
    # Headwater logic
    parent_row$residual_area_km2 <- parent_row$total_area_km2
    parent_row$all_nested_ids <- NA
    residual_results[[i]] <- parent_row
    
  } else {
    # Subtraction logic (using only immediate children to avoid topology errors)
    #sf_use_s2(TRUE)
    parent_geom <- st_geometry(parent_row)
    parent_shape <- st_geometry(parent_geom)
    children_geom <- catchments[immediate_child_indices, ]
    #children_union <- st_union(children_geom)
    #children_union <- st_union(st_geometry(catchments[immediate_child_indices, ]))
    children_union <- st_union(st_geometry(children_geom))
    #diff_geom <- st_difference(parent_shape, children_union)
    
    # 3. Clean up the result
    # st_as_sf handles the "nrow == length" check automatically by creating a new table
   # inter_catchment <- st_as_sf(data.frame(ID = parent_geom$catch_id), geometry = diff_geom)
    
    # 4. Optional: If the result is a GEOMETRYCOLLECTION (points + lines + polygons)
    # we usually only want the Polygons
   # inter_catchment <- st_collection_extract(inter_catchment, "POLYGON")
    
    # Subtract and clean up
    diff_geom <- st_difference(st_make_valid(parent_shape), st_make_valid(children_union)) %>%
      st_collection_extract("POLYGON") 
    
    # 2. Extract only the LARGEST polygon (removes slivers/islands)
    # st_cast turns a MULTIPOLYGON into individual POLYGON rows
    individual_polygons <- st_make_valid(st_cast(diff_geom, "POLYGON"))
    
    if (length(individual_polygons) > 1) {
      areas <- st_area(individual_polygons)
      diff_geom <- individual_polygons[which.max(areas)]
    } else {
      diff_geom <- individual_polygons
    }
    
    # Check if empty after cleaning
    if (length(diff_geom) == 0 || st_is_empty(diff_geom)) {
      message(paste("Catchment", parent_row$catch_id, "is fully covered. Skipping."))
      next
    }
    
    # Create the result row
    new_row <- parent_row
    st_geometry(new_row) <- diff_geom
    a=st_area(parent_geom)/1e6
    b=st_area(diff_geom)/1e6
    new_row$residual_area_km2 <- as.numeric(st_area(diff_geom)) / 1e6
    new_row$all_nested_ids <- all_child_ids
    new_row$immediate_nested_ids <- immediate_child_ids
    
    if(new_row$residual_area_km2>new_row$area_km2+1){
      print(c(a,b))
      ggplot() +
        # The original full parent (faded)
        geom_sf(data = parent_geom, fill = "grey90", color = "black", linetype = "dotted") +
        # The sub-catchments (the holes)
       geom_sf(data = children_geom, fill = "white", color = "red") +
        # The Inter-catchment area (the actual contributing land)
       geom_sf(data = new_row, fill = "green", alpha = 0.7) +
        theme_minimal() +
        labs(title = "Inter-catchment (Residual) Area",
             subtitle = "Green represents area draining directly to the reach between stations")
      
      
    } 
    
    residual_results[[i]] <- new_row
  }
}

# 3. Final Table Assembly
final_residuals <- do.call(rbind, residual_results)

which(is.na(final_residuals$catch_id))
#final_residuals<-final_residuals[-which(is.na(final_residuals$catch_id)),]

# Recommended: Save as GeoPackage (preserves long names and long strings)
st_write(final_residuals, "D:/tilloal/Documents/01_Projects/RegimeShifts/data/catchments_analysis_final_v3.gpkg", delete_dsn = TRUE)


library(ggplot2)
library(viridis)

# Count the number of nested catchments
# We split the string by comma and count the elements
final_residuals$nesting_count <- sapply(strsplit(as.character(final_residuals$all_nested_ids), ","), 
                                        function(x) {
                                          if (all(is.na(x))) return(0) 
                                          length(trimws(x))
                                        })

ggplot(data = final_residuals) +
  # Plot the residual polygons colored by how many catchments they contain
  geom_sf(aes(fill=nesting_count) ,color = "white", size = 0.1) +
  # Using the 'viridis' scale makes it readable for colorblindness and printing
 scale_fill_viridis_c(option = "plasma", trans="sqrt", name = "Number of\nNested Sub-catchments") +
  theme_minimal() +
  labs(
    title = "Catchment Nesting Complexity Map",
    subtitle = "Darker areas represent major basins; lighter areas are headwaters",
    caption = paste("Total catchments processed:", nrow(final_residuals))
  ) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

#I have to redo this and somehow add other catchments....

# Recommended: Save as GeoPackage (preserves long names and long strings)
st_write(final_residuals, "catchments_analysis_results_v3.gpkg", delete_dsn = TRUE)


