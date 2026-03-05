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

#Hybas07
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=inner_join(Catamere07,outhybas,by= c("llcoord"="latlong"))

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
cst7=st_transform(Catf7,  crs=3035)
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

#keep only rivers in EU domain
out1=outletopen(hydroDir,"efas_rnet_100km_01min")
out1$latlong=paste(round(out1$Var1,4),round(out1$Var2,4),sep=" ")
outhybas=inner_join(out1,outhybas, by="latlong")


CatUpA=inner_join(Catf7,UpArea,by=c("llcoord"="latlong"))

matcat=match(outhybas$latlong,CatUpA$llcoord)
CatUpA=CatUpA[matcat,]
ratUp=CatUpA$SUB_AREA/CatUpA$upa

plot(ratUp, log="y")
abline(h=1.5)
abline(h=0.5)

CatUpA$outlet=1
CatUpA$outlet[which(ratUp>1.2)]=2
CatUpA$outlet[which(ratUp<0.8)]=3
length(CatUpA$outlet[which(ratUp>1.5)])
CatUpA$outlet[which(CatUpA$outlet==3 & CatUpA$upa<2e4)]=2



Cathybas=CatUpA[which(CatUpA$outlet==1),]

## 2‑a  LDD raster (NetCDF)  -------------------------------------------------
nc_path <- "D:/tilloal/Documents/06_Floodrivers/mapscal/ldd_European_01min.nc"                 # <-- change to your file
ldd_raw <- rast(nc_path)   # replace "LDD" with the exact variable name
# (If you already have a raster on disk you can skip this block.)

## 2‑b  Outlet points ---------------------------------------------------------
outlets <- st_as_sf(UpArea, coords = c("Var1.x", "Var2.x"), crs = 4326)
# ppl <- st_transform(ppl, crs = 3035)          # make CRS identical to the raster


## 2‑c  Existing catchment polygons -------------------------------------------
poly_path <- "D:/tilloal/Documents/01_Projects/RegimeShifts/OtherCatchments_polygon.shp"
catch_pol <- st_read(poly_path, quiet = TRUE) %>%
  st_transform(crs(ldd_raw))               # same CRS as everything else


catch_pol <- st_read(poly_path, quiet = TRUE) %>%   # quiet = TRUE hides the progress messages
  st_make_valid() %>%                               # fixes self‑intersections, etc.
  st_transform(crs(ldd_raw))                          # any CRS you like; keep it consistent later


#load matched hybas catchments as well
hybas_path <- "D:/tilloal/Documents/01_Projects/RegimeShifts/data/catchments_hybas.shp"
catch_hy <- st_read(hybas_path, quiet = TRUE) %>%
  st_make_valid() %>%   
  st_transform(crs(ldd_raw)) 

catch_hy<-catch_hy[,c(30,36)]
colnames(catch_hy)[1]=colnames(catch_pol)[1]
catch_pol$group="other"
catch_hy$group="hybas"
#all my catchments are here for the nested analysis
catch_ph=rbind(catch_pol,catch_hy)
# -----------------------------------------------------------------
# 2‑c  Make sure we have a column called `catch_id`
# -----------------------------------------------------------------
# If the original file already has a field that looks like an ID (e.g. "ID",
# "CatchID", "PolyID", …) we rename it to `catch_id`.  If not, we create one.
catchments=catch_ph
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
catchments<-catchments[-which(is.na(catchments$catch_id)),]
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
matched_data <- match(catchments$catch_id,outlets$outlets.x)
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
  intersections <- st_intersection(child_geoms, parent_buffered)
  intersections <- st_make_valid(intersections)
  # Calculate areas for all intersections at once
  inter_areas <- as.numeric(st_area(intersections)) / 1e6
  
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

# Find indices where the relationship is TRUE
relations2 <- which(nest_matrix2, arr.ind = TRUE)

# Create a mapping table
nesting_df2 <- data.frame(
  parent_id = catchments$catch_id[relations2[,1]], # Change 'ID' to your column name
  child_id  = catchments$catch_id[relations2[,2]]
)

# Count how many polygons contain each specific polygon
catchments$nesting_level <- rowSums(refined_nest_matrix)

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
st_write(final_residuals, "D:/tilloal/Documents/01_Projects/RegimeShifts/data/catchments_analysis_final.gpkg", delete_dsn = TRUE)


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
  geom_sf(aes(fill=residual_area_km2) ,color = "white", size = 0.1) +
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
st_write(final_residuals, "catchments_analysis_results.gpkg", delete_dsn = TRUE)

#meteo(for next time)
library(dplyr)

# Let's say your variable is 'precip_mm'
# 1. Create a function to calculate residual meteorology

#upload from Q_components file
Qv="bluewater"
#Qe is the sum of Q components

if (Qv=="bluewater"){
  Qe=fread(file="D:/tilloal/Documents/01_Projects/RegimeShifts/data/BlueWater.csv",header=TRUE)
}

library(tidyverse)


library(sf)
library(dplyr)
library(data.table)
library(sf)

catchments=final_residuals
matcat=match(outhybas$latlong,UpArea$latlong)
UpArea=UpArea[matcat,]

outputfilenames = c('rainUpsX','snowMeltUpsX','qUzUpsX', 'qLZUpsX','surfaceRunoffUpsX',
                   'infUpsX','ActEvapo','disWin')
for (var in outputfilenames){
  #var=outputfilenames[1]
  print(var)
  Vari<- fread(paste0(hydroDir,"/tss/HERA_Histo/",var,"_1951_2020.csv"),header=TRUE)
  time=Vari$V1
  Vari=Vari[order(time),]
  matcol=match(UpArea$outlets,as.numeric(colnames(Vari)))
  Vari <- Vari[, .SD, .SDcols = matcol]
  
  # 2. Create a copy of the data.table to store residual results
  # This preserves the time/index columns
  Qe<-Vari
  res_met_dt <- copy(Qe)
  
  # 3. The Loop
  residual_results <- list()
  
  for (i in 1:nrow(catchments)) {
    print(i)
    #i=304
    p_id <- as.character(catchments$catch_id[i])
    p_area <- catchments$area_km2[i]
    
    # Identify children using your refined_nest_matrix
    all_child_idx <- which(refined_nest_matrix[i, ])
    imm_child_idx <- all_child_idx[colSums(refined_nest_matrix[all_child_idx, all_child_idx, drop=FALSE]) == 0]
    
    # --- CASE 1: HEADWATERS (No children) ---
    if (length(imm_child_idx) == 0) {
      # Geometry and Met stay the same
      new_geom <- st_geometry(catchments[i, ])
      res_area <- p_area
      # No changes needed to res_met_dt[ , p_id] as it's already a copy
      
    } else {
      # --- CASE 2: NESTED CATCHMENTS ---
      child_ids <- as.character(catchments$catch_id[imm_child_idx])
      child_areas <- catchments$area_km2[imm_child_idx]
      
      # A. Geometry Subtraction (Largest Piece Logic)
      children_union <- st_union(st_geometry(catchments[imm_child_idx, ]))
      diff_poly <- st_make_valid(st_difference(st_make_valid(st_geometry(catchments[i, ])), 
                                               st_make_valid(children_union)) %>%
                                   st_collection_extract("POLYGON") %>% 
                                   st_cast("POLYGON"))
      
      if(length(diff_poly) == 0) next
      new_geom <- diff_poly[which.max(st_area(diff_poly))]
      res_area <- as.numeric(st_area(new_geom)) / 1e6
      
      # B. Mass Balance using data.table syntax
      # We calculate: (ParentVal * ParentArea - Sum(ChildVal * ChildArea)) / ResArea
      
      # Calculate total volume of children for all rows (hours) at once
      # .SD allows us to perform math across the specific child columns
      child_vols <- Qe[, mapply(`*`, .SD, child_areas), .SDcols = child_ids]
      
      # If multiple children, rowSums them; if one child, it's already a vector
      if (length(child_ids) > 1) {
        total_child_vol <- rowSums(child_vols, na.rm = TRUE)
      } else {
        total_child_vol <- as.vector(child_vols)
      }
      
      # Use 'set' for high-performance column update
      # pmax(0, ...) ensures no negative precipitation/radiation
      res_series <- pmax(0, (Qe[[p_id]] * p_area - total_child_vol) / res_area)
      
      # plot(Qe[[p_id]][1:6000]*(4*p_area*1000000)/(1000*3600*24) ,
      #      col=1,type="l",ylim=c(1,6500))
      # lines(total_child_vol*(4*1000000)/(1000*3600*24),col=2)
      # lines(res_series[1:6000]*(4*res_area*1000000)/(1000*3600*24),col=3,type="l")
      # rs=sum(res_series*(4*res_area*1000000)/(1000*3600*24))/sum(Qe[[p_id]]*(4*p_area*1000000)/(1000*3600*24))
      # rs*100
      set(res_met_dt, j = p_id, value = res_series)
    }
    
    # C. Store geometry results
    res_row <- catchments[i, ]
    st_geometry(res_row) <- new_geom
    res_row$residual_area_km2 <- res_area
    residual_results[[i]] <- res_row
  }
  
  # 4. Finalize
  final_residual_sf <- do.call(rbind, residual_results)
  
  #compare for catchments for which I know there should be no difference
  ardeche="314805"
  arde1=res_met_dt[[ardeche]]
  plot(arde1[1:300], type="l")
  arde0=Qe[[ardeche]]
  lines(arde0,col=2)
  
  fwrite(res_met_dt,
         paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/data/",var,"_CC_1951_2020.csv"),
         sep = ",",           # Standard CSV separator
         na = "-9999",        # Common "No Data" flag for hydrologic models (or use "NA")
         row.names = FALSE,   # Don't include the row indices (1, 2, 3...)
         quote = FALSE,       # No quotes around numbers/IDs (keeps file size small)
         showProgress = TRUE, # Watch the progress bar for large files
         nThread = 2)         # Use 4 CPU cores to speed up the compression/writing)
  
}



calc_residual_meteo <- function(parent_idx, met_var, sf_data, n_matrix) {
  
  # Get immediate children
  child_idx <- get_immediate_children(parent_idx, n_matrix)
  
  parent_val <- sf_data[[met_var]][parent_idx]
  parent_area <- sf_data$total_area_km2[parent_idx]
  res_area <- sf_data$residual_area_km2[parent_idx]
  
  if (length(child_idx) == 0) {
    # Headwater: residual value is just the parent value
    return(parent_val)
  } else {
    # Sum of (Value * Area) for all immediate children
    children_total <- sum(sf_data[[met_var]][child_idx] * sf_data$total_area_km2[child_idx])
    
    # Apply the mass balance formula
    residual_val <- ((parent_val * parent_area) - children_total) / res_area
    return(residual_val)
  }
}

# 2. Apply to your dataset (e.g., for Precipitation)
final_residuals$precip_residual <- sapply(1:nrow(final_residuals), function(i) {
  calc_residual_meteo(i, "precip_mm", final_residuals, nest_matrix)
})




# Remove self‑containment
for (i in seq_along(contain_list)) {
  contain_list[[i]] <- setdiff(contain_list[[i]], i)
}

parent_id_vec    <- integer(nrow(catch_pol))   # 0 = no parent
children_ids_list <- vector("list", nrow(catch_pol))

for (i in seq_len(nrow(catch_pol))) {
  print(i)
  child_rows <- contain_list[[i]]
  
  if (length(child_rows) == 0) {
    children_ids_list[[i]] <- integer(0)
    next
  }
  
  child_ids <- catch_pol$catch_id[child_rows]
  children_ids_list[[i]] <- child_ids
  
  # Record immediate parent (choose the smallest‑area parent if multiple)
  for (cid in child_ids) {
    child_row <- which(catch_pol$catch_id == cid)
    
    if (parent_id_vec[child_row] == 0) {
      parent_id_vec[child_row] <- catch_pol$catch_id[i]
    } else {
      current_parent_row <- which(catch_pol$catch_id == parent_id_vec[child_row])
      if (as.numeric(st_area(catch_pol[i, ])) <
          as.numeric(st_area(catch_pol[current_parent_row, ]))) {
        parent_id_vec[child_row] <- catch_pol$catch_id[i]
      }
    }
  }
}

catch_pol <- catch_pol %>%
  mutate(parent_id   = parent_id_vec,
         children_ids = children_ids_list)

# --------------------------------------------------------------
# 4. Union children & compute residual (parent – children) ----
# --------------------------------------------------------------
# --------------------------------------------------------------
# 4‑a  Prepare containers for the results
# --------------------------------------------------------------
children_union_list <- vector("list", nrow(catch_pol))   # union of children (or NULL)
residual_geom_list  <- vector("list", nrow(catch_pol))   # parent – children



# --------------------------------------------------------------
# 4‑b  Loop over every catchment and compute the two geometries
# --------------------------------------------------------------
for (i in seq_len(nrow(catch_pol))) {
  # ----- 4‑b‑i  Union of all direct children -----------------
  print(i)
  child_ids <- catch_pol$children_ids[[i]]
  
  if (length(child_ids) == 0) {
    # No children → union is NULL
    children_union_list[[i]] <- NULL
  } else {
    # Subset the children geometries and take the union
    child_geoms <- catch_pol %>%
      filter(catch_id %in% child_ids) %>%
      st_geometry()
    children_union_list[[i]] <- st_union(child_geoms)
  }
  
}
children_union_list <- lapply(children_ids_list, function(ids) {
  if (length(ids) == 0) return(NULL)
  st_union(catch_pol %>% filter(catch_id %in% ids) %>% st_geometry())
})

children_union_list <- vector("list", nrow(catch_pol))

for (i in seq_len(nrow(catch_pol))) {
  child_ids <- catch_pol$children_ids[[i]]
  children_union_list[[i]] <- safe_union_children(child_ids)
}

# Store the result back in the catchment table (optional)
catch_pol$children_union <- children_union_list


catch_pol2 <- catch_pol %>%
  mutate(children_ids   = children_ids_list,
         children_union = children_union_list)


# --------------------------------------------------------------
# 3‑a  Create a data.frame that pairs each catch_id with its union
# --------------------------------------------------------------
union_df <- data.frame(
  catch_id = catch_pol2$catch_id,
  union_geom = I(catch_pol2$children_union)   # I() keeps the list‑column as is
)


# --------------------------------------------------------------
# Helper: turn a single element (NULL or sfc) into a valid geometry
# --------------------------------------------------------------
make_valid_geom <- function(g) {
  if (is.null(g) || length(g) == 0) {
    # No children → return an *empty* geometry collection.
    # This is recognised by sf as a valid (but empty) geometry.
    return(st_geometrycollection())
  }
  # If we already have an sfc object, just return it
  if (inherits(g, "sfc")) return(g)
  
  # Fallback – should not happen, but keep the code safe
  stop("Unexpected geometry type")
}

clean_polygon <- function(g) {
  # 1) Exact duplicate vertices
  g <- st_remove_repeated_points_R(g)
  
  # 2) Zero‑width buffer (fixes self‑intersection, tiny slivers)
  g <- st_buffer(g, 0)
  
  # 3) Robust make_valid (lwgeom version)
  g <- lwgeom::st_make_valid(g)
  
  g
}

# --------------------------------------------------------------
# A.  Helper that works on a *single* geometry (sfg)
# --------------------------------------------------------------
clean_polygon_sfg <- function(g) {
  # g is an sfg (e.g. POLYGON, MULTIPOLYGON, LINESTRING, …)
  
  # 1) Remove exact duplicate vertices
  g <- st_remove_repeated_points_R(g)   # works on an sfg as well
  
  # 2) Buffer‑zero – fixes many self‑intersections & tiny slivers
  g <- st_buffer(g, 0)
  
  # 3) Robust make_valid (lwgeom version – splits multipart if needed)
  #    If you do not have lwgeom, just use st_make_valid(g) instead.
  if (requireNamespace("lwgeom", quietly = TRUE)) {
    g <- lwgeom::st_make_valid(g)
  } else {
    g <- st_make_valid(g)
  }
  
  g
}

st_remove_repeated_points_R <- function(x) {
  if (inherits(x, "sf")) {
    st_geometry(x) <- st_remove_repeated_points_R(st_geometry(x))
    return(x)
  }
  if (!inherits(x, "sfc"))
    stop("`x` must be an sf or sfc object", call. = FALSE)
  
  # Helper that removes consecutive duplicates from a single ring (matrix)
  clean_ring <- function(ring_mat) {
    # keep the first point, then any point that differs from the previous one
    keep <- c(TRUE, !(ring_mat[-1, 1] == ring_mat[-nrow(ring_mat), 1] &
                        ring_mat[-1, 2] == ring_mat[-nrow(ring_mat), 2]))
    out <- ring_mat[keep, , drop = FALSE]
    
    # ensure closed ring (first == last)
    if (!all(out[1, ] == out[nrow(out), ]))
      out <- rbind(out, out[1, ])
    out
  }
  
  # Apply to each geometry
  lapply(x, function(g) {
    if (is.null(g)) return(g)
    
    # POLYGON ---------------------------------------------------------
    if (inherits(g, "POLYGON")) {
      rings <- lapply(g, clean_ring)
      st_polygon(rings)
    }
    # MULTIPOLYGON ----------------------------------------------------
    else if (inherits(g, "MULTIPOLYGON")) {
      polys <- lapply(g, function(poly) {
        lapply(poly, clean_ring)
      })
      st_multipolygon(polys)
    }
    # LINESTRING / MULTILINESTRING – we usually don't need to clean them,
    # but we can drop duplicates if you wish:
    else if (inherits(g, c("LINESTRING", "MULTILINESTRING"))) {
      # For lines we keep the first point of each segment only
      # (duplicate vertices are harmless for most operations)
      g
    }
    else {
      # Other geometry types (POINT, GEOMETRYCOLLECTION, etc.) are returned unchanged
      g
    }
  }) %>% st_sfc(crs = st_crs(x))
}

# -----------------------------------------------------------------
# 4‑a  Empty‑geometry placeholder (used when a catch has no children)
# -----------------------------------------------------------------
empty_geom <- st_geometrycollection()   # valid empty geometry

# -----------------------------------------------------------------
# 4‑b  Function that safely creates a union for a set of child IDs
# -----------------------------------------------------------------
safe_union_children <- function(child_ids) {
  # No children → return an empty geometry
  if (length(child_ids) == 0) return(empty_geom)
  
  # Grab the child geometries
  child_geoms <- catch_pol2 %>%
    filter(catch_id %in% child_ids) %>%
    st_geometry()
  
  # --------------------------------------------------------------
  # Try the *plain* union first – it is fast and works in most cases
  # --------------------------------------------------------------
  plain_union <- tryCatch(
    st_union(child_geoms),
    error = function(e) NULL   # will be NULL if GEOS complains
  )
  
  # --------------------------------------------------------------
  # If the plain union succeeded and is not degenerate, keep it
  # --------------------------------------------------------------
  if (!is.null(plain_union) && !st_is_empty(plain_union)) {
    # Quick validity check – if it passes we are done
    if (st_is_valid(plain_union)) return(plain_union)
  }
  
  # --------------------------------------------------------------
  # Otherwise clean each child and union again (this almost always
  # fixes duplicate‑vertex problems)
  # --------------------------------------------------------------
  cleaned_children <- map(child_geoms, clean_polygon_sfg)
  
  # Re‑union the cleaned geometries (wrap in tryCatch just in case)
  cleaned_union <- tryCatch(
    st_union(do.call(c, cleaned_children)),
    error = function(e) {
      message("  *** STILL FAILED after cleaning for children: ", paste(child_ids, collapse = ","))
      empty_geom   # fallback – keep an empty geometry so the script continues
    }
  )
  
  # Final safety net – if the cleaned union is still invalid, replace it
  if (!st_is_valid(cleaned_union) || st_is_empty(cleaned_union)) {
    message("  *** UNION remains invalid or empty for parent with children: ", paste(child_ids, collapse = ","))
    return(empty_geom)
  }
  
  cleaned_union
}

# --------------------------------------------------------------
# Build a new sfc list that has the same length as catch_pol
# --------------------------------------------------------------
union_sfc <- vector("list", nrow(catch_pol))

for (i in seq_len(nrow(catch_pol))) {
  union_sfc[[i]] <- make_valid_geom(catch_pol2$children_union[[i]])
}

# Turn the plain list into an sfc object and assign the CRS of the catchments
union_sfc <- st_sfc(union_sfc, crs = st_crs(catch_pol))

# --------------------------------------------------------------
# 3‑b  Convert to an sf object (rows with NULL become empty geometries)
# --------------------------------------------------------------

union_sf <- st_sf(
  catch_id   = catch_pol2$catch_id,   # keep the identifier
  union_geom = union_sfc,            # the geometry column we just built
  crs        = st_crs(catch_pol2)    # ensure the CRS is set
)

# Quick sanity check
print(union_sf)

# --------------------------------------------------------------
# 3‑c  Quick look at the first few rows
# --------------------------------------------------------------
head(union_sf)

is_list_column <- inherits(union_sf$union_geom, "list")
cat("union_geom is a list‑column :", is_list_column, "\n")
# -----------------------------------------------------------------
# Helper: turn NULL / empty entries into an *empty geometry collection*
# -----------------------------------------------------------------
make_valid_geom <- function(g) {
  if (is.null(g) || length(g) == 0) {
    # Empty geometry – still a valid sf object
    return(st_geometrycollection())
  }
  if (inherits(g, "sfc")) return(g)   # already a geometry list
  stop("Unexpected object type in union_geom")
}

# Build a new sfc vector that has the same length as union_sf
union_sfc <- vector("list", nrow(union_sf))

for (i in seq_len(nrow(union_sf))) {
  union_sfc[[i]] <- make_valid_geom(union_sf$union_geom[[i]])
}

# Collapse the plain list into an *sfc* object and keep the CRS of the original layer
union_sfc <- st_sfc(union_sfc, crs = st_crs(union_sf))

# Replace the old (list) column with the new geometry column
union_sf <- union_sf %>%
  st_set_geometry(NULL) %>%          # drop the old geometry column completely
  mutate(union_geom = union_sfc) %>% # add the new, proper geometry column
  st_as_sf()                         # coerce back to an sf object

union_sf <- union_sf %>%
  mutate(
    union_geom = st_cast(union_geom, "POLYGON")   # or "MULTIPOLYGON"
  )


union_sf <- union_sf %>%
  mutate(
    # -----------------------------------------------------------------
    # 4‑a  Is the union geometry valid? (TRUE / FALSE / NA)
    # -----------------------------------------------------------------
    is_valid = st_is_valid(union_geom),
    
    # -----------------------------------------------------------------
    # 4‑b  Geometry type (e.g. "POLYGON", "MULTIPOLYGON", NA)
    # -----------------------------------------------------------------
    geom_type = ifelse(st_is_empty(union_geom),
                       NA_character_,
                       as.character(st_geometry_type(union_geom))),
    # 
    # # -----------------------------------------------------------------
    # # 4‑c  Area (numeric, in the units of the CRS)
    # # -----------------------------------------------------------------
    area = ifelse(st_is_empty(union_geom),
                  NA_real_,
                  as.numeric(st_area(union_geom))),
    # 
    # # -----------------------------------------------------------------
    # # 4‑d  Number of parts (useful when you get MULTIPOLYGON)
    # # -----------------------------------------------------------------
    # n_parts = ifelse(st_is_empty(union_geom),
    #                  NA_integer_,
    #                  lengths(st_geometry(union_geom)))
  )
which(union_sf$is_valid!=T)

# Returns a logical vector: TRUE = valid, FALSE = invalid
valid_vec <- st_is_valid(union_sf$union_geom)

# If you want *why* a feature is invalid, use the detail version
invalid_detail <- st_is_valid(union_sf[!valid_vec, ])

# Print a short summary
cat("Valid features :", sum(valid_vec), "\n")
cat("Invalid features:", sum(!valid_vec), "\n")


  # ----- 4‑b‑ii  Residual (parent minus children) ----------
  parent_geom <- st_geometry(catch_pol[i, ])
  
  if (is.null(children_union_list[[i]])) {
    # No children → the residual is just the original polygon
    residual_geom_list[[i]] <- parent_geom
  } else {
    # Subtract the union of children from the parent
    # This may return an empty geometry collection if the children cover the parent completely.
    residual_geom_list[[i]] <- st_difference(parent_geom,
                                             children_union_list[[i]])
  }
}
