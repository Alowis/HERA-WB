


#Analysis of the timeseries outputs of LISFLOOD-----
#Each output corresponds to a watershed In Europe-----
#Extraction of frost season-----
# Library calling --------------------------------------------------
suppressWarnings(suppressMessages(library(ncdf4)))
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(reshape)
library(raster)
library(dplyr)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(raster)
library(RtsEva)
library(moments)
#  Function declaration ---------------------------------------------------
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
domainopen=function(dir,outletname,nrspace=rep(NA,5)){
  ncbassin=paste0(dir,"/",outletname,".nc")
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  if ("Band1"%in% name.vb)namev="Band1"
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
  
  # outll=outll[which(!is.na(outlets)),]
  # outlets=outlets[which(!is.na(outlets))]
  outll=data.frame(outlets,outll)
  return (outll)
}
seasony=function(x){
  theta=x*(2*pi/365.25)
  # plot(theta)
  
  xi=1/(length(theta))*sum(cos(theta))
  yi=1/(length(theta))*sum(sin(theta))
  if (xi<=0){
    Di=(atan(yi/xi)+pi)*(365.25/(2*pi))
  }else if(xi>0 & yi>=0){
    Di=(atan(yi/xi))*(365.25/(2*pi))
  }else if(xi>0 & yi<0){
    Di=(atan(yi/xi)+2*pi)*(365.25/(2*pi))
  }
  R=sqrt(xi^2+yi^2)
  return(c(Di,R))
}
UpAopen=function(dir,outletname,Sloc_final){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[2]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncb,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat)
  lla=length(latdat)
  start=c(1,1)
  count=c(llo,lla)
  
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)/1000000
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$upa=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  #outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outll$latlong=paste(round(outll$Var1,4),round(outll$Var2,4),sep=" ")
  outfinal=inner_join(outll, Sloc_final, by="latlong")
  return (outfinal)
}

log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}
# Prepare work environment -------
setwd("D:/tilloal/Documents/LFRuns_utils")
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
write.table(Catmrdik, file = "D:/tilloal/Documents/01_Projects/RegimeShifts/catextra.txt", sep = "\t", row.names = FALSE, col.names = FALSE)   
Catmrdtest=Catmrdik[c(1:10),]
write.table(Catmrdtest, file = "D:/tilloal/Documents/01_Projects/RegimeShifts/catest.txt", sep = "\t", row.names = FALSE, col.names = FALSE)   


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



CatLR=CatUpA[which(CatUpA$upa>2e4),]
st_geometry(CatLR)=NULL
CatO=CatUpA[-which(CatUpA$upa>2e4),]
ppl <- st_as_sf(CatLR, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

ppc <- st_as_sf(CatO, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppc <- st_transform(ppc, crs = 3035)

# Compute centroids
CatO_centroids <- st_centroid(ppc)
alle=match(colnames(ppl),colnames(CatO_centroids))
ppoints=rbind(ppl,CatO_centroids[,alle])

#Land use
library(exactextractr)

#in 1951
rastforest=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracforest_1951.tif")
rastsealed=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracsealed_1951.tif")
rastirrigated=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracirrigated_1951.tif")
rastother=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracother_1951.tif")
rastrice=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracrice_1951.tif")
rastwater=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracwater_1951.tif")
palet=c(hcl.colors(9, palette = "BuPu", alpha = NULL, rev = TRUE, fixup = TRUE))


LandUse=CatRcompose
forestchange<- exact_extract(rastforest, LandUse, 'mean')
LandUse$forestchange=forestchange

ziz=match(ppoints$HYBAS_ID,LandUse$HYBAS_ID)
ppoints$forest=LandUse$forestchange[ziz]

LandUse$sealed <- exact_extract(rastsealed, LandUse, 'mean')[ziz]
LandUse$irrigated <- exact_extract(rastirrigated, LandUse, 'mean')[ziz]
LandUse$other <- exact_extract(rastother, LandUse, 'mean')[ziz]
LandUse$rice <- exact_extract(rastrice, LandUse, 'mean')[ziz]
LandUse$water <- exact_extract(rastwater, LandUse, 'mean')[ziz]

save(LandUse,file="D:/tilloal/Documents/01_Projects/RegimeShifts/Landuse_init.Rdata")


#Change
rastforest=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracforest_ch20201951.tif")
rastsealed=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracsealed_ch20201951.tif")
rastirrigated=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracirrigated_ch20201951.tif")
rastother=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracother_ch20201951.tif")
rastrice=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracrice_ch20201951.tif")
rastwater=raster( "D:/tilloal/Documents/06_Floodrivers/landuse/fracwater_ch20201951.tif")
palet=c(hcl.colors(9, palette = "Spectral", alpha = NULL, rev = FALSE, fixup = TRUE))


name_legend<-" fraction (%)"
name_legend<-" fraction change 2020-1951 (%)"

GHshpp=CatRcompose
forestchange<- exact_extract(rastforest, GHshpp, 'mean')
GHshpp$forestchange=forestchange

ziz=match(ppoints$HYBAS_ID,GHshpp$HYBAS_ID)
ppoints$forest=GHshpp$forestchange[ziz]
  
ppoints$sealedchange <- exact_extract(rastsealed, GHshpp, 'mean')[ziz]
ppoints$irrigatedchange <- exact_extract(rastirrigated, GHshpp, 'mean')[ziz]
ppoints$otherchange <- exact_extract(rastother, GHshpp, 'mean')[ziz]
ppoints$ricechange <- exact_extract(rastrice, GHshpp, 'mean')[ziz]
ppoints$waterchange <- exact_extract(rastwater, GHshpp, 'mean')[ziz]
mhh=match(Catf7$HYBAS_ID,GHshpp$HYBAS_ID)
GHshppH=ppoints[,-c(2:24,26:34)]
df_GHshppH=(GHshppH)
st_geometry(df_GHshppH)<-NULL
df_GHshppH=data.frame(df_GHshppH)

lucmap=list()
luclass=c("forest","sealed","irrigated","other","rice","water")

for (li in 1:length(luclass))
{
  lu=luclass[li]
  print(lu)
  GHshppH$fill=as.numeric(df_GHshppH[,2+li])
  
  flims=(quantile(GHshppH$fill,c(0.01,0.95),na.rm=T))
  #lims=c(0,round(max(abs(flims)),1))
  lims=c(-round(max(abs(flims)),1),round(max(abs(flims)),1))
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),2),round(max(abs(flims)),2))
  }
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),4),round(max(abs(flims)),4))
  }

  #lims=c(-25,25)
  lucmap<-ggplot(basemap) +
    geom_sf(fill="gray95",color="transparent",size=0.5)+
    geom_sf(data=GHshppH,aes(geometry=geometry,size=upa,fill=fill*100),color="black",alpha=1,shape=21,stroke=0.1)+
    geom_sf(fill="transparent",color="gray30",size=0.5)+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_fill_gradientn(
      colors=palet,
      limits=lims*100,oob = scales::squish,
      name=paste0(lu,name_legend))   +
    labs(x="Longitude", y = "Latitude")+
    guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
           size= "none")+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=16),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))+
    ggtitle(lu)
  lucmap
  ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/lu_change_",lu,"_1951-2020_upa.jpg"), lucmap,width=30, height=20, units=c("cm"),dpi=300)
  
}


#Find a way to load nectdf and convert them to rasters


#I need to use cutmaps to extract my exact upstream area for all the guys


#WATER DEMAND
#Same for water demand change



#divide by area
workDir<-("D:/tilloal/Documents/06_Floodrivers/")
ncpix= paste0(workDir,"/mapscal/pixarea_European_01min.nc")
ncp=nc_open(ncpix)
t=ncp$var[[2]]
tsize<-t$varsize
pixarea=ncvar_get(ncp,names(ncp[['var']])[2]) 
plon=ncvar_get(ncp,"lon") 
plat=ncvar_get(ncp,"lat") 
pixarea=as.matrix(t(pixarea))

rast_totwd<-raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/wateruse_sums/all_demands_2020.tif")
rast_totwd<-raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/wateruse_sums/all_demands_1951.tif")
#Convert RasterLayer to matrix
rast_tmat <- (as.matrix(rast_totwd))
#multiply by 30.4 to have the real sum values
rast_tmat=rast_tmat*30.4
#m3/m2 to m3
rast_tmat=rast_tmat*pixarea
#convert mm to m3
rast_tmat=rast_tmat/1e3
#m3 to km3
rast_tmat=rast_tmat/1e9
rast_tmout <- raster(nrows=nrow(rast_totwd), ncols=ncol(rast_totwd), ext=extent(rast_totwd))
crs(rast_tmout) <- crs(rast_tmat)
values(rast_tmout) <- rast_tmat
rast_totwd=rast_tmout

ppoints$totalwd2020 <- exact_extract(rast_totwd, GHshpp, 'sum')[ziz]
sum(ppoints$totalwd2020)
ppoints$wdpkm2=ppoints$totalwd2020/ppoints$upa*1000*1000




# mhh=match(UnHY,GHshpp$HYBAS_ID)
# hybas07$totalwd2020 <- exact_extract(rast_totwd, hybas07, 'sum')
# hybas07H=hybas07[mhh,]
# sum(hybas07H$totalwd2020,na.rm=T)
# hybas07H$wdpkm2=hybas07H$totalwd2020/hybas07H$UP_AREA*1000*1000
palet2=c(hcl.colors(9, palette = "BuPu", alpha = NULL, rev = TRUE, fixup = TRUE))
#unit is now in mm
wd="total water demand in 1951"
tsize=12
flims=(quantile(ppoints$wdpkm2,c(0.1,0.95),na.rm=T))
lims=c(0,round(max(abs(flims)),1))
lims=c(0,120)
wdmap<-ggplot(basemap) +
  geom_sf(fill="gray95",color="transparent",size=0.5)+
  # geom_sf(data=ppoints,aes(fill=wdpkm2,geometry=geometry),alpha=1,color="transparent")+
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=wdpkm2),color="black",alpha=1,shape=21,stroke=0)+
  geom_sf(fill="transparent",color="gray30",size=0.5)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet2,
    limits=lims,oob = scales::squish, trans="sqrt", breaks=c(1,5,20,50,100,200),
    name=paste0("(mm/year)"))   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=16),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  ggtitle(wd)

wdmap
ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/TotalWaterDemand_1951.jpg"), wdmap,width=30, height=20, units=c("cm"),dpi=300)



rast_ene=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/ene_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_enemat <- (as.matrix(rast_ene))
#multiply by 30.4 to have the real sum values
rast_enemat=rast_enemat*30.4
#m3/m2 to m3
rast_enemat=rast_enemat*pixarea
#convert mm to m3
rast_enemat=rast_enemat/1e3
#m3 to km3
rast_enemat=rast_enemat/1e9


rast_eneout <- raster(nrows=nrow(rast_ene), ncols=ncol(rast_ene), ext=extent(rast_ene))
crs(rast_eneout) <- crs(rast_ene)
values(rast_eneout) <- rast_enemat
rast_ene=rast_eneout

rast_dom=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/dom_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_dommat <- (as.matrix(rast_dom))
#multiply by 30.4 to have the real sum values
rast_dommat=rast_dommat*30.4
#m3/m2 to m3
rast_dommat=rast_dommat*pixarea
#convert mm to m3
rast_dommat=rast_dommat/1e3
#m3 to km3
rast_dommat=rast_dommat/1e9
rast_domout <- raster(nrows=nrow(rast_dom), ncols=ncol(rast_dom), ext=extent(rast_dom))
crs(rast_domout) <- crs(rast_dom)
values(rast_domout) <- rast_dommat
rast_dom=rast_domout
rast_liv=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/liv_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_livmat <- (as.matrix(rast_liv))
#multiply by 30.4 to have the real sum values
rast_livmat=rast_livmat*30.4
#m3/m2 to m3
rast_livmat=rast_livmat*pixarea
#convert mm to m3
rast_livmat=rast_livmat/1e3
#m3 to km3
rast_livmat=rast_livmat/1e9
rast_livout <- raster(nrows=nrow(rast_liv), ncols=ncol(rast_liv), ext=extent(rast_liv))
crs(rast_livout) <- crs(rast_liv)
values(rast_livout) <- rast_livmat
rast_liv=rast_livout

rast_ind=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/ind_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_indmat <- (as.matrix(rast_ind))
#multiply by 30.4 to have the real sum values
rast_indmat=rast_indmat*30.4
#m3/m2 to m3
rast_indmat=rast_indmat*pixarea
#convert mm to m3
rast_indmat=rast_indmat/1e3
#m3 to km3
rast_indmat=rast_indmat/1e9
rast_indout <- raster(nrows=nrow(rast_ind), ncols=ncol(rast_ind), ext=extent(rast_ind))
crs(rast_indout) <- crs(rast_ind)
values(rast_indout) <- rast_indmat
rast_ind=rast_indout

rast_total=raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/all_ysum_ch20201951.tif")
#Convert RasterLayer to matrix
rast_totalmat <- (as.matrix(rast_total))
#multiply by 30.4 to have the real sum values
rast_totalmat=rast_totalmat*30.4
#m3/m2 to m3
rast_totalmat=rast_totalmat*pixarea
#convert mm to m3
rast_totalmat=rast_totalmat/1e3
#m3 to km3
rast_totalmat=rast_totalmat/1e9
rast_totalout <- raster(nrows=nrow(rast_total), ncols=ncol(rast_total), ext=extent(rast_total))
crs(rast_totalout) <- crs(rast_total)
values(rast_totalout) <- rast_totalmat
rast_total=rast_totalout


ppoints$enechange <- exact_extract(rast_ene, GHshpp, 'sum')[ziz]
ppoints$domchange <- exact_extract(rast_dom, GHshpp, 'sum')[ziz]
ppoints$livchange <- exact_extract(rast_liv, GHshpp, 'sum')[ziz]
ppoints$indchange <- exact_extract(rast_ind, GHshpp, 'sum')[ziz]
ppoints$totalchange <- exact_extract(rast_total, GHshpp, 'sum')[ziz]


GHshppH=ppoints[,-c(2:24,26:34)]
df_GHshppH=(GHshppH)
st_geometry(df_GHshppH)<-NULL
df_GHshppH=data.frame(df_GHshppH)

# mhh=match(UnHY,hybas07$HYBAS_ID)
# hybas07H=hybas07[mhh,]
# df_hybas07H=data.frame(hybas07H)
# st_geometry(df_hybas07H)<-NULL

wdclass=c("ene","dom","liv","ind","total")
wdmap=list()
for (li in 1:length(wdclass))
{
  wd=wdclass[li]
  print(wd)
  GHshppH$fill=as.numeric(df_GHshppH[,10+li])
  GHshppH$fill=GHshppH$fill/GHshppH$upa*1000*1000
  
  flims=(quantile(GHshppH$fill,c(0.05,0.95),na.rm=T))
  lims=c(-round(max(abs(flims)),1),round(max(abs(flims)),1))
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),2),round(max(abs(flims)),2))
  }
  if (diff(lims)==0){
    lims=c(-round(max(abs(flims)),4),round(max(abs(flims)),4))
  }
  tsize=12
  wdmap[[li]]<-ggplot(basemap) +
    geom_sf(fill="gray95",color="transparent",size=0.5)+
    #geom_sf(data=HRwgs84h,aes(fill=fill,geometry=geometry),alpha=1,color="transparent")+
    geom_sf(data=GHshppH,aes(geometry=geometry,size=upa,fill=fill),color="black",alpha=1,shape=21,stroke=0.1)+
    geom_sf(fill="transparent",color="gray30",size=0.5)+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_fill_gradientn(
      colors=palet,
      limits=lims,oob = scales::squish,
      name=paste0("Change (mm/year)"))   +
    labs(x="Longitude", y = "Latitude")+
    guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
           size= "none")+
    theme(axis.title=element_text(size=tsize),
          title = element_text(size=16),
          axis.text=element_text(size=osize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "right",
          panel.grid.major = element_line(colour = "grey70"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(1, "cm"))+
    ggtitle(wd)
  ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/wdchangeCat_",wd,"_20201951.jpg"),wdmap[[li]],width=30, height=20, units=c("cm"),dpi=300)
  
}




#RESERVOIRS



#Volume

rast_totVolume1<-raster( "D:/tilloal/Documents/LFRuns_utils/data/reservoirs/reservoirs_volumes_1951.tif")
#rast_totwd<-raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/wateruse_sums/all_demands_1951.tif")
#Convert RasterLayer to matrix
rast_totV1mat <- (as.matrix(rast_totVolume1))
# #multiply by 30.4 to have the real sum values
# rast_tmat=rast_tmat*30.4
# #m3/m2 to m3
# rast_tmat=rast_tmat*pixarea
# #convert mm to m3
# rast_tmat=rast_tmat/1e3
#m3 to km3
rast_totV1mat=rast_totV1mat/1e9
rast_tvout <- raster(nrows=nrow(rast_totV1mat), ncols=ncol(rast_totV1mat), ext=extent(rast_totVolume1))
crs(rast_tvout) <- crs(rast_totVolume1)
values(rast_tvout) <- rast_totV1mat
rast_totV1=rast_tvout

ppoints$totalvolumeR1 <- exact_extract(rast_totV1, GHshpp, 'sum')[ziz]

ppoints$tv1pkm2=ppoints$totalvolumeR1/ppoints$upa*1000*1000




# mhh=match(UnHY,GHshpp$HYBAS_ID)
# hybas07$totalwd2020 <- exact_extract(rast_totwd, hybas07, 'sum')
# hybas07H=hybas07[mhh,]
# sum(hybas07H$totalwd2020,na.rm=T)
# hybas07H$wdpkm2=hybas07H$totalwd2020/hybas07H$UP_AREA*1000*1000
palet2=c(hcl.colors(9, palette = "BuPu", alpha = NULL, rev = TRUE, fixup = TRUE))
#unit is now in mm
wd="total reservoir storage capacity in 1951"
tsize=12
flims=(quantile(ppoints$totalvolumeR1,c(0.1,0.95),na.rm=T))
lims=c(0,round(max(abs(flims)),1))
#lims=c(0,120)
rstormap<-ggplot(basemap) +
  geom_sf(fill="gray95",color="transparent",size=0.5)+
  # geom_sf(data=ppoints,aes(fill=wdpkm2,geometry=geometry),alpha=1,color="transparent")+
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=totalvolumeR1),color="black",alpha=1,shape=21,stroke=0.1)+
  geom_sf(fill="transparent",color="gray30",size=0.5)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet2,
    limits=lims,oob = scales::squish,
    name=paste0("(km3)"))   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=16),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  ggtitle(wd)

rstormap
ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/TotalReservoirCapacity_1951.jpg"), rstormap,width=30, height=20, units=c("cm"),dpi=300)



#Reservoir ratio


resOpen=function(dir,outletname,ValidSY){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncb,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat)
  lla=length(latdat)
  start=c(1,1)
  count=c(llo,lla)
  
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(lla:1))
  outll$res.ratio=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outlp=outll[which(!is.na(outll$res.ratio)),]
  mf=match(ValidSY$idlalo,outlp$idlalo)
  ValidSY$res.ratio=outlp$res.ratio[mf]
  return (ValidSY)
}
#Initial

#1951
#load reservoir ratio map
res_path<-("D:/tilloal/Documents/LFRuns_utils/data/reservoirs/")
outletname="res_ratio_European_01min_1951.nc"
dir=res_path

outhybas1$idlalo=paste(outhybas1$idlo,outhybas1$idla,sep=" ")
ResData=resOpen(res_path,outletname,outhybas1)  
rsx=ResData
length(rsx$res.ratio[which(rsx$res.ratio<0.5)])/length(rsx$res.ratio)
rsx$res.group=1
rsx$res.group[which(rsx$res.ratio>0 & rsx$res.ratio<=0.2)]=2
rsx$res.group[which(rsx$res.ratio>0.2 & rsx$res.ratio<=0.5)]=3
rsx$res.group[which(rsx$res.ratio>=0.5 & rsx$res.ratio<=1)]=3
rsx$res.group[which(rsx$res.ratio>=1 & rsx$res.ratio<=2)]=4
rsx$res.group[which(rsx$res.ratio>2)]=5


colorz = c('#d73027','orange','#fee090','lightblue','royalblue',"darkblue")
kgelabs=c("< -0.41","-0.41 - 0","0 - 0.2", "0.2 - 0.5","0.5 - 0.75", ">0.75")

matloc=match(ppoints$outlets.x,outhybas1$outlets)
ppoints$Ratio1 <-rsx$res.ratio[matloc]

#ppoints$tv1pkm2=ppoints$totalvolumeR1/ppoints$upa*1000*1000

# mhh=match(UnHY,GHshpp$HYBAS_ID)
# hybas07$totalwd2020 <- exact_extract(rast_totwd, hybas07, 'sum')
# hybas07H=hybas07[mhh,]
# sum(hybas07H$totalwd2020,na.rm=T)
# hybas07H$wdpkm2=hybas07H$totalwd2020/hybas07H$UP_AREA*1000*1000
palet2=c(hcl.colors(9, palette = "BuPu", alpha = NULL, rev = TRUE, fixup = TRUE))
#unit is now in mm
wd="reservoir capacity ratio (C) in 1951"
tsize=12
max(ppoints$Ratio1)
flims=(quantile(ppoints$Ratio1,c(0.01,0.99),na.rm=T))
lims=c(0,round(max(abs(flims)),1))
#lims=c(0,2)
rratmap<-ggplot(basemap) +
  geom_sf(fill="gray95",color="transparent",size=0.5)+
  # geom_sf(data=ppoints,aes(fill=wdpkm2,geometry=geometry),alpha=1,color="transparent")+
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=Ratio1),color="black",alpha=1,shape=21,stroke=0.1)+
  geom_sf(fill="transparent",color="gray30",size=0.5)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet2,
    limits=lims,oob = scales::squish,
    name=paste0("C"))   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=16),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  ggtitle(wd)

rratmap
ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/ReservoirRatio_1951.jpg"), rratmap,width=30, height=20, units=c("cm"),dpi=300)




#2020

#load reservoir ratio map
res_path<-("D:/tilloal/Documents/LFRuns_utils/data/reservoirs/")
outletname="res_ratio_European_01min_2020.nc"
dir=res_path

outhybas1$idlalo=paste(outhybas1$idlo,outhybas1$idla,sep=" ")
ResData=resOpen(res_path,outletname,outhybas1)  
rsx=ResData
length(rsx$res.ratio[which(rsx$res.ratio<0.5)])/length(rsx$res.ratio)
rsx$res.group=1
rsx$res.group[which(rsx$res.ratio>0 & rsx$res.ratio<=0.2)]=2
rsx$res.group[which(rsx$res.ratio>0.2 & rsx$res.ratio<=0.5)]=3
rsx$res.group[which(rsx$res.ratio>=0.5 & rsx$res.ratio<=1)]=3
rsx$res.group[which(rsx$res.ratio>=1 & rsx$res.ratio<=2)]=4
rsx$res.group[which(rsx$res.ratio>2)]=5


colorz = c('#d73027','orange','#fee090','lightblue','royalblue',"darkblue")
kgelabs=c("< -0.41","-0.41 - 0","0 - 0.2", "0.2 - 0.5","0.5 - 0.75", ">0.75")

matloc=match(ppoints$outlets.x,outhybas1$outlets)
ppoints$Ratio2 <-rsx$res.ratio[matloc]

#ppoints$tv1pkm2=ppoints$totalvolumeR1/ppoints$upa*1000*1000

# mhh=match(UnHY,GHshpp$HYBAS_ID)
# hybas07$totalwd2020 <- exact_extract(rast_totwd, hybas07, 'sum')
# hybas07H=hybas07[mhh,]
# sum(hybas07H$totalwd2020,na.rm=T)
# hybas07H$wdpkm2=hybas07H$totalwd2020/hybas07H$UP_AREA*1000*1000
palet2=c(hcl.colors(9, palette = "BuPu", alpha = NULL, rev = TRUE, fixup = TRUE))
#unit is now in mm
wd="reservoir capacity ratio (C) in 2020"
tsize=12
max(ppoints$Ratio2)
flims=(quantile(ppoints$Ratio1,c(0.01,0.99),na.rm=T))
lims=c(0,round(max(abs(flims)),1))
lims=c(0,.5)
rratmap<-ggplot(basemap) +
  geom_sf(fill="gray95",color="transparent",size=0.5)+
  # geom_sf(data=ppoints,aes(fill=wdpkm2,geometry=geometry),alpha=1,color="transparent")+
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=Ratio2),color="black",alpha=1,shape=21,stroke=0.1)+
  geom_sf(fill="transparent",color="gray30",size=0.5)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet2,
    limits=lims,oob = scales::squish,
    name=paste0("C"))   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=16),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  ggtitle(wd)

rratmap
ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/ReservoirRatio_2020.jpg"), rratmap,width=30, height=20, units=c("cm"),dpi=300)



#Difference

#load reservoir ratio map
res_path<-("D:/tilloal/Documents/LFRuns_utils/data/reservoirs/")
outletname="res_ratio_diff_2020-1951.nc"
dir=res_path

outhybas1$idlalo=paste(outhybas1$idlo,outhybas1$idla,sep=" ")
ResData=resOpen(res_path,outletname,outhybas1)  
rsx=ResData
length(rsx$res.ratio[which(rsx$res.ratio<0.5)])/length(rsx$res.ratio)
rsx$res.group=1
rsx$res.group[which(rsx$res.ratio>0 & rsx$res.ratio<=0.2)]=2
rsx$res.group[which(rsx$res.ratio>0.2 & rsx$res.ratio<=0.5)]=3
rsx$res.group[which(rsx$res.ratio>=0.5 & rsx$res.ratio<=1)]=3
rsx$res.group[which(rsx$res.ratio>=1 & rsx$res.ratio<=2)]=4
rsx$res.group[which(rsx$res.ratio>2)]=5


colorz = c('#d73027','orange','#fee090','lightblue','royalblue',"darkblue")
kgelabs=c("< -0.41","-0.41 - 0","0 - 0.2", "0.2 - 0.5","0.5 - 0.75", ">0.75")

matloc=match(ppoints$outlets.x,outhybas1$outlets)
ppoints$RatioD <-rsx$res.ratio[matloc]

#ppoints$tv1pkm2=ppoints$totalvolumeR1/ppoints$upa*1000*1000

# mhh=match(UnHY,GHshpp$HYBAS_ID)
# hybas07$totalwd2020 <- exact_extract(rast_totwd, hybas07, 'sum')
# hybas07H=hybas07[mhh,]
# sum(hybas07H$totalwd2020,na.rm=T)
# hybas07H$wdpkm2=hybas07H$totalwd2020/hybas07H$UP_AREA*1000*1000
palet2=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
#unit is now in mm
wd="change in reservoir capacity ratio in (2020-1951)"
tsize=12
max(ppoints$Ratio1)
flims=(quantile(ppoints$RatioD,c(0.01,0.99),na.rm=T))
lims=c(0,round(max(abs(flims)),1))
lims=c(-0.5,.5)
rratmap<-ggplot(basemap) +
  geom_sf(fill="gray95",color="transparent",size=0.5)+
  # geom_sf(data=ppoints,aes(fill=wdpkm2,geometry=geometry),alpha=1,color="transparent")+
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=RatioD),color="black",alpha=1,shape=21,stroke=0.1)+
  geom_sf(fill="transparent",color="gray30",size=0.5)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet2,
    limits=lims,oob = scales::squish,
    name=paste0("C"))   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=16),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  ggtitle(wd)

rratmap
ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/ReservoirRatio_change.jpg"), rratmap,width=30, height=20, units=c("cm"),dpi=300)



rast_totVolume1<-raster( "D:/tilloal/Documents/LFRuns_utils/data/reservoirs/reservoirs_volumes_1951.tif")
#rast_totwd<-raster( "D:/tilloal/Documents/06_Floodrivers/wateruse/wateruse_sums/all_demands_1951.tif")
#Convert RasterLayer to matrix
rast_totV1mat <- (as.matrix(rast_totVolume1))
# #multiply by 30.4 to have the real sum values
# rast_tmat=rast_tmat*30.4
# #m3/m2 to m3
# rast_tmat=rast_tmat*pixarea
# #convert mm to m3
# rast_tmat=rast_tmat/1e3
#m3 to km3
rast_totV1mat=rast_totV1mat/1e9
rast_tvout <- raster(nrows=nrow(rast_totV1mat), ncols=ncol(rast_totV1mat), ext=extent(rast_totVolume1))
crs(rast_tvout) <- crs(rast_totVolume1)
values(rast_tvout) <- rast_totV1mat
rast_totV1=rast_tvout

ppoints$totalvolumeR1 <- exact_extract(rast_totV1, GHshpp, 'sum')[ziz]

ppoints$tv1pkm2=ppoints$totalvolumeR1/ppoints$upa*1000*1000




# mhh=match(UnHY,GHshpp$HYBAS_ID)
# hybas07$totalwd2020 <- exact_extract(rast_totwd, hybas07, 'sum')
# hybas07H=hybas07[mhh,]
# sum(hybas07H$totalwd2020,na.rm=T)
# hybas07H$wdpkm2=hybas07H$totalwd2020/hybas07H$UP_AREA*1000*1000
palet2=c(hcl.colors(9, palette = "BuPu", alpha = NULL, rev = TRUE, fixup = TRUE))
#unit is now in mm
wd="total reservoir storage capacity in 1951"
tsize=12
flims=(quantile(ppoints$totalvolumeR1,c(0.1,0.95),na.rm=T))
lims=c(0,round(max(abs(flims)),1))
#lims=c(0,120)
rstormap<-ggplot(basemap) +
  geom_sf(fill="gray95",color="transparent",size=0.5)+
  # geom_sf(data=ppoints,aes(fill=wdpkm2,geometry=geometry),alpha=1,color="transparent")+
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=totalvolumeR1),color="black",alpha=1,shape=21,stroke=0)+
  geom_sf(fill="transparent",color="gray30",size=0.5)+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet2,
    limits=lims,oob = scales::squish,
    name=paste0("(km3)"))   +
  labs(x="Longitude", y = "Latitude")+
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
  theme(axis.title=element_text(size=tsize),
        title = element_text(size=16),
        axis.text=element_text(size=osize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(1, "cm"))+
  ggtitle(wd)

rstormap
ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/TotalReservoirCapacity_1951.jpg"), rstormap,width=30, height=20, units=c("cm"),dpi=300)


#Save a data.frame with all the vairables>

p2g=match(ppoints$HYBAS_ID,GHshppH$HYBAS_ID)

cn=colnames(ppoints)
SaveVars=ppoints[,c(1,20,25,36:53)]

save(SaveVars,file="D:/tilloal/Documents/01_Projects/RegimeShifts/cathcment_variables_change.RData")
st_geometry(SaveVars)<-NULL
colnames(SaveVars)
Affinvar=SaveVars[,-c(1:3,7,10:15,17:20)]
library(corrplot)
res <- cor(Affinvar)
round(res, 2)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
