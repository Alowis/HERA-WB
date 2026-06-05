# Install packages if not already installed
# install.packages("raster")
# install.packages("sf")

# Load the necessary libraries
library(raster)
library(sf)
library(stars)
# Define the list of PCRaster file paths
pcraster_files1 <- list.files(path = "D:/tilloal/Documents/01_Projects/RegimeShifts/catmasks/", pattern = "\\.map$", full.names = TRUE)
pcraster_files2 <- list.files(path = "D:/tilloal/Documents/01_Projects/RegimeShifts/catmasks/hybas", pattern = "\\.map$", full.names = TRUE)

pcraster_files=c(pcraster_files1,pcraster_files2)
# Initialize an empty list to store the polygons
polygons_list <- list()
i=0
# Loop through each PCRaster file
for (pcraster_file in pcraster_files) {
  i=i+1
  if (i %% 50 == 0) cat("Processing catchment:", pcraster_file, "\n")

  r <- read_stars(pcraster_file)
  
  # Convert raster to polygons
  p_sf <- st_as_sf(r, as_points = FALSE, merge = TRUE)
  
  # # Extract the numeric part of the filename
  file_name <- basename(pcraster_file)
  numeric_part <- gsub("[^0-9]", "", file_name)  # Remove non-numeric characters

  # Add the numeric part as a new column in the sf object
  p_sf$Name <- numeric_part
  
  # Merge polygons into a single polygon
  merged_polygon <- st_union(p_sf)
  
  # Create a new sf object with the merged polygon
  merged_sf <- st_sf(Name = numeric_part, geometry = merged_polygon)
  
  # Append the merged sf object to the list
  polygons_list[[length(polygons_list) + 1]] <- merged_sf
 
  # Append the sf object to the list
  #polygons_list[[length(polygons_list) + 1]] <- p_sf
}


# Define a template with the desired column names
desired_columns <- c("geometry","Name")

# polygons_list2<-polygons_list
# # Standardize each sf object
# polygons_list2 <- lapply(polygons_list2, function(p) {
#   # Check if the necessary columns are present
#   missing_cols <- setdiff(desired_columns, names(p))
#   
#   # Add missing columns with NA values
#   for (col in missing_cols) {
#     p[[col]] <- NA
#   }
#   
#   # Ensure the order and presence of columns
#   p <- p[desired_columns]
#   
#   return(p)
# })

# Combine all the polygon sf objects into one
combined_polygons <- do.call(rbind, polygons_list)
st_crs(combined_polygons) <- 4326 
# Optionally, you might want to ensure there are no overlapping geometries
#combined_polygons <- st_union(combined_polygons)

# Save the combined polygons as a shapefile
st_write(combined_polygons, "D:/tilloal/Documents/01_Projects/RegimeShifts/All_catchment_raw.shp")


#Step 2: load upstream catchments, match them with ID and with catchments, compute total area

#load tiff file of catchments
outletopen=function(dir,outletname,nrspace=rep(NA,5)){
  ncbassin=paste0(dir,"/",outletname,".nc")
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  if (!is.na(nrspace[1])){
    start=as.numeric(nrspace[c(2,4)])
    count=as.numeric(nrspace[c(3,5)])-start+1
  }else{
    londat = ncvar_get(ncb,name.lon) 
    llo=length(londat)
    latdat = ncvar_get(ncb,name.lat)
    lla=length(latdat)
    start=c(1,1)
    count=c(llo,lla)
  }
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll=outll[which(!is.na(outlets)),]
  outlets=outlets[which(!is.na(outlets))]
  outll=data.frame(outlets,outll)
  return (outll)
}
library(stars)

setwd ("D:/tilloal/Documents/01_Projects/RegimeShifts")

hydroDir<-("D:/tilloal/Documents/01_Projects/RegimeShifts")
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
outletname="catchment_upstreamF"
outletname="fluxmerde"

outupstream=outletopen(hydroDir,outletname)
outupstream$latlong=paste(round(outupstream$Var1,4),round(outupstream$Var2,4),sep=" ")

outpustreamT<-outupstream[which(outupstream$outlets<3),]

#match with outhybas
# Prepare work environment -------

hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
outletname="outletsv8_hybas07_01min"
outhybas=outletopen(hydroDir,outletname)
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4),sep=" ")

outpustreamT<-inner_join(outpustreamT,outhybas,by="latlong")


#Hybas07
Catchmentrivers7=read.csv(paste0(hydroDir,"/Catchments/from_hybas_eu_onlyid.csv"),encoding = "UTF-8", header = T, stringsAsFactors = F)
hybas07 <- read_sf(dsn = paste0(hydroDir,"/Catchments/hydrosheds/hybas_eu_lev07_v1c.shp"))
hybasf7=fortify(hybas07) 
Catamere07=inner_join(hybasf7,Catchmentrivers7,by= "HYBAS_ID")
Catamere07$llcoord=paste(round(Catamere07$POINT_X,4),round(Catamere07$POINT_Y,4),sep=" ") 
Catf7=inner_join(Catamere07,outhybas,by= c("llcoord"="latlong"))

#Aggregate at catchment level and plot
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


#Plot Hybas points by catchment area

###load UpArea -----
#load upstream area
# main_path = 'D:/tilloal/Documents/06_Floodrivers/'
# valid_path = paste0(main_path,'DataPaper/')
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data/")
outletname="upArea_European_01min.nc"
#dir=valid_path
# outhybas$idlalo=paste(outhybas$idlo, outhybas$idla, sep=" ")
# outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4), sep=" ")
UpArea=UpAopen(hydroDir,outletname,outhybas)
head(UpArea)


#keep only rivers in EU domain
out1=outletopen(hydroDir,"efas_rnet_100km_01min")
out1$latlong=paste(round(out1$Var1,4),round(out1$Var2,4),sep=" ")
mout=which(!is.na(match(outhybas$latlong,out1$latlong)))
outhybas1=outhybas[mout,]


#save txt file
outtxt

Catf7=inner_join(Catamere07,outhybas1,by= c("llcoord"="latlong"))

#Plot parameters
cord.dec=UpArea[,c(1,2)]
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

ppl <- st_as_sf(UpArea, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,size=upa),color="transparent",alpha=.9,shape=21,stroke=0,fill="blue")+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  # scale_fill_gradientn(
  #   colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
  #   breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5,reverse=F),
         size= guide_legend(override.aes = list(fill = "grey50")))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))



CatUpA=inner_join(Catf7,UpArea,by=c("llcoord"="latlong"))

matcat=match(outhybas1$latlong,CatUpA$llcoord)
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



Catmrdik=CatUpA[which(CatUpA$outlet!=1),c(22,23,21)]
st_geometry(Catmrdik)=NULL
colnames(Catmrdik)<-c("X","Y","ID")
Catmrdik=as.matrix(Catmrdik)

Catother <- read_sf(dsn ="D:/tilloal/Documents/01_Projects/RegimeShifts/OtherCatchments_polygon.shp")
Catother$outlets.x<-as.numeric(Catother$Name)
Catother<-Catother[,-1]
#link Cathother with catmrdik
Catmrdik=CatUpA[which(CatUpA$outlet!=1),c(22,23,21,1)]
st_geometry(Catmrdik)<-NULL

#match option
mcat=match(Catother$outlets.x,Catmrdik$outlets.x)
Catother$HYBAS_ID=Catmrdik$HYBAS_ID[mcat]
# Catother=inner_join(Catmrdik,Catother, by=c("outlets.x"))
# Catother<-Catother[,-5]

Cathybas1<-CatUpA[which(CatUpA$outlet==1),c(21,1)]

cx=crs(Cathybas1)
# Catother1<-st_as_sf(Catother, crs = 4326)
# Catother1 <- st_transform(Catother, crs = cx)
CatRcompose<-rbind(Catother,Cathybas1)
CatRcompose<-CatRcompose[-which(is.na(CatRcompose$outlets.x)),]
CatRcompose<-st_as_sf(CatRcompose)


CatRcompose$upstream=0
mc=match(outpustreamT$outlets.y,CatRcompose$outlets.x)
CatRcompose$upstream[mc]=1
Catrp=CatRcompose[which(CatRcompose$upstream==1),]


library(lwgeom)  # Ensure lwgeom is installed
fixed_polygons <- st_make_valid(Catrp)
# Union the polygons to handle overlaps
unioned_polygons <- st_union(fixed_polygons)

# Calculate the total area without overlaps
fixed_union <- st_make_valid(unioned_polygons)
total_area <- st_area(fixed_union)

# Print total area in square meters
print(total_area)

# Convert to square kilometers
total_area_sq_km <- total_area / 1e6
print(total_area_sq_km)

#plot CatRcompose

ppl <- st_as_sf(CatRcompose, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=Catrp,aes(geometry=geometry),color="transparent",alpha=.9,shape=21,stroke=0, fill="royalblue")+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  # scale_fill_gradientn(
  #   colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
  #   breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5,reverse=F),
         size= guide_legend(override.aes = list(fill = "grey50")))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))



mup=match(outpustreamT$outlets.y,CatUpA$outlets.x)
CatLR=CatUpA[-mup,]
st_geometry(CatLR)=NULL
CatLR$up="Large Rivers"
CatO=CatUpA[mup,]
CatO$up="Upstream Catchments"

ppl <- st_as_sf(CatLR, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

ppc <- st_as_sf(CatO, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppc <- st_transform(ppc, crs = 3035)

# Compute centroids
CatO_centroids <- st_centroid(ppc)
alle=match(colnames(ppl),colnames(CatO_centroids))
ppoints=rbind(ppl,CatO_centroids[,alle])



ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa, fill=up),color="black",alpha=.8,shape=21,stroke=0.1)+ 
  # geom_sf(data=CatO_centroids,aes(geometry=geometry,size=upa,fill=factor(outlet)),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  # scale_fill_gradientn(
  #   colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
  #   breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  scale_fill_manual(name= "group",values=c("Large Rivers"="royalblue","Upstream Catchments"="forestgreen"))+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_legend(override.aes = list(size = 8), reverse=T), size = "none")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/upstreamgroups.png", width=20, height=15, units=c("cm"),dpi=1500)

ppointCx<-ppoints

st_write(ppointCx,dsn="D:/tilloal/Documents/01_Projects/RegimeShifts/Upstreamgroups.shp")
muc=match(outpustreamT$outlets.y,CatUpA$outlets.x)
upaup=CatUpA$upa[muc]

hist(upaup)
s1<-sum(upaup)
s1

#load the netcdf of pixel area
workDir<-("//ies.jrc.it/H07/ClimateRun4/nahaUsers/tilloal/ERA5l_x_lisflood/EFAS1rcmin_GitLFS_efas5_24032023/")
aname="/pixarea_European_01min.nc"
nccal=paste0(workDir,aname)
nca=nc_open(nccal)
name.vb=names(nca[['var']])
namev=name.vb[2]
#time <- ncvar_get(ncb,"time")

#timestamp corretion
name.lon="lon"
name.lat="lat"
start=c(1,1)
count=c(llo,lla)
londat = ncvar_get(nca,name.lon,start=start[1],count=count[1]) 
llo=length(londat)
latdat = ncvar_get(nca,name.lat,start=start[2],count=count[2])
lla=length(latdat)
idcal = ncvar_get(nca,namev,start = start, count= count) 
idcal=as.vector(idcal)
max(idcal,na.rm=T)
area=expand.grid(londat,latdat)
#lonlatloop=expand.grid(c(1:llo),c(1:lla))
area$V1=idcal/1000000
sa=sum(area$V1)
area$latlong=paste(round(area$Var1,4),round(area$Var2,4),sep=" ")

hydroDir<-("D:/tilloal/Documents/06_Floodrivers/DataPaper/GIS/")
GridHERA=raster( paste0(hydroDir,"HERA_domain_01min.tif"))
GHERA=as.data.frame(GridHERA,xy=T)
GHERA=GHERA[which(!is.na(GHERA$HERA_domain_01min)),]
#Joining area with HERA domain
GHERA$latlong=paste(round(GHERA$x,4),round(GHERA$y,4),sep=" ")

mag=match(GHERA$latlong,area$latlong)
GHERA$area=area$V1[mag]
areaH=sum(GHERA$area)
total_area_sq_km/areaH
