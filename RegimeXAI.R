
#load functions

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
outhybas$idlalo=paste(outhybas$idlo, outhybas$idla, sep=" ")
outhybas$latlong=paste(round(outhybas$Var1,4),round(outhybas$Var2,4), sep=" ")
UpArea=UpAopen(hydroDir,outletname,outhybas)
head(UpArea)


#keep only rivers in EU domain
out1=outletopen(hydroDir,"efas_rnet_100km_01min")
out1$latlong=paste(round(out1$Var1,4),round(out1$Var2,4),sep=" ")
outhybas=inner_join(out1,outhybas, by="latlong")


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


#load data



#load catchment and clusters
ppoints<-st_read("D:/tilloal/Documents/01_Projects/RegimeShifts/Regime_shifts_hclust.shp")

#load hydrometeo inputs
load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/PrecipitationMetrics.Rdata")
load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/AETMetrics.Rdata")
load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/E0Metrics.Rdata")
load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/runofffraction.Rdata")
load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/snowfraction.Rdata")


#preprocess meteo inputs
SaveAET=data.frame(SaveAET)
colnames(SaveAET)=colnames(SaveAgg)
SaveE0=data.frame(SaveE0)
colnames(SaveE0)=colnames(SaveAgg)

ys1=c(1951:1980)
my=which(!is.na(match(SaveAgg$SaveAgg,ys1)))
Keep=(SaveAgg[my,-1])
mean_snowfracion=colMeans(Keep)

cnames<-names(mean_snowfracion)
#match ppoints with different tables
ppoints$otlts_x=ppoints$outlets_x
sf2pp=match(ppoints$otlts_x,names(mean_snowfracion))
ppoints$snowfraction=mean_snowfracion[sf2pp]


palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=snowfraction),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="snow fraction (%)",
    breaks=seq(0,1,.2), limits=c(0,0.5)) +
  #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
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


#boxplot by cluster

#snowfraction

agSF=aggregate(list(sfrac=ppoints$snowfraction),
                by = list(rclass=ppoints$C0name),
                FUN = function(x) c(mean(x)))


# merdecol=match(ValidSY$UpAgroup,names(meds))
# ppoints$col=meds[merdecol]
# class_labels <- c("Nival", "Nivo-pluvial","Pluvio-nival", "Persistent Pluvial",  "Erratic Pluvial")
# palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))
# color_regimesn5<-c("Nival"="red4","Nivo-pluvial"="orchid","Pluvio-nival"="turquoise4","Persistent Pluvial"="olivedrab3","Erratic Pluvial"="orangered")
ppoints$C0name<- factor(ppoints$C0name , levels=class_labels)
p1<-ggplot(ppoints, aes(x=factor(C0name), y=snowfraction)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=C0name),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(0,1),name="snow fraction",breaks = seq(-1,1,by=0.2),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(breaks=class_labels, name="Regime class")+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_regimesn5,name="Regimes")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p1

#Precipitation

myp=which(!is.na(match(SavePrecip$year,ys1)))

KeepP=SavePrecip[myp,c(1,2,5)]

KeepPag<- stats::aggregate(list(pmean=KeepP$val.sum), by = list(catch = KeepP$cat), 
                           FUN = function(x) c(pmean=mean(x, na.rm = T)))
KeepPag <- do.call(data.frame, KeepPag)

#Add ET0
my=which(!is.na(match(SaveE0$SaveAgg,ys1)))
KeepE0=(SaveE0[my,-1])
Ke0_means=colMeans(KeepE0)

kpe=match(KeepPag$catch,cnames)

KeepPag$et=Ke0_means[kpe]

KeepPag$AI=KeepPag$pmean/KeepPag$et

#Plot of AI
mai=match(ppoints$otlts_x,KeepPag$catch)
ppoints$AI1=KeepPag$AI[mai]
ppoints$ET01=KeepPag$et[mai]
ppoints$Precip1=KeepPag$pmean[mai]

p2<-ggplot(ppoints, aes(x=factor(C0name), y=AI1)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=C0name),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(0,6),name="Aridity Index",breaks = seq(-1,8,by=1),minor_breaks = seq(-1,8,0.5))+
  scale_x_discrete(breaks=class_labels, name="Regime class")+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_regimesn5,name="Regimes")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p2


p3<-ggplot(ppoints, aes(x=factor(C0name), y=Precip1)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=C0name),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(0,2500),name="Yearly precipitation (mm)",breaks = seq(0,3000,by=500),minor_breaks = seq(0,3000,100))+
  scale_x_discrete(breaks=class_labels, name="Regime class")+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_regimesn5,name="Regimes")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p3





#Plot of rainfall concentration
myp=which(!is.na(match(SavePrecip$year,ys1)))

KeepP=SavePrecip[myp,c(1,4,5)]

KeepPag<- stats::aggregate(list(px=KeepP$val.xtr), by = list(catch = KeepP$cat), 
                           FUN = function(x) c(pmean=mean(x, na.rm = T)))
KeepPag <- do.call(data.frame, KeepPag)

mai=match(ppoints$otlts_x,KeepPag$catch)
ppoints$px1=KeepPag$px[mai]

p4<-ggplot(ppoints, aes(x=factor(C0name), y=px1)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=C0name),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(0,1),name="Precipitation concentration",breaks = seq(0,1,by=.2),minor_breaks = seq(0,1,.1))+
  scale_x_discrete(breaks=class_labels, name="Regime class")+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_regimesn5,name="Regimes")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p4



my=which(!is.na(match(SaveAET$SaveAgg,ys1)))
KeepAET=(SaveAET[my,-1])
Ket_means=colMeans(KeepAET)


KeepPag$aet1=Ket_means[kpe]



palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=AI1),color="black",alpha=1,shape=21,stroke=.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="AET/P", limits = c(0,1)) +
  #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
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


#load landuse input
load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/landuse_init.Rdata")

mai=match(ppoints$otlts_x,LandUse$outlets.x)
ppoints$forest=LandUse$forestchange[mai]
ppoints$sealed=LandUse$sealed[mai]
ppoints$water=LandUse$water[mai]
ppoints$irrigated=LandUse$irrigated[mai]

p5<-ggplot(ppoints, aes(x=factor(Cshift), y=forest)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(0,1),name="forest fraction",breaks = seq(0,1,by=.2),minor_breaks = seq(0,1,.1))+
  scale_x_discrete(breaks=class_labels, name="Regime class")+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  #scale_fill_manual(values=color_regimesn5,name="Regimes")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p5

p6<-ggplot(ppoints, aes(x=factor(C0name), y=sealed)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=C0name),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(0,.3),name="sealed fraction",breaks = seq(0,1,by=.1),minor_breaks = seq(0,1,.05))+
  scale_x_discrete(breaks=class_labels, name="Regime class")+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_regimesn5,name="Regimes")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p6


p7<-ggplot(ppoints, aes(x=factor(C0name), y=irrigated)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=C0name),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(0,1),name="irrigated fraction",breaks = seq(0,1,by=.1),minor_breaks = seq(0,1,.05), trans="sqrt")+
  scale_x_discrete(breaks=class_labels, name="Regime class")+
  coord_cartesian(ylim=c(0,0.5))+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_regimesn5,name="Regimes")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p7

load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/cathcment_variables_change.RData")



mai=match(ppoints$otlts_x,SaveVars$outlets.x)
ppoints$Rratio=SaveVars$Ratio1[mai]



p8<-ggplot(ppoints, aes(x=factor(C0name), y=Rratio)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=C0name),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(0,5),name="reservoir ratio",breaks = seq(0,1,by=.2),minor_breaks = seq(0,1,.1), trans="sqrt")+
  scale_x_discrete(breaks=class_labels, name="Regime class")+
  coord_cartesian(ylim=c(0,0.5))+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_regimesn5,name="Regimes")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p8



#Now focus on changes

class_labels <- c( "Earlier snomelt", "Strong wetting - winter peak", "Strong drying - winter peak","Drying - summer peak","Drying - spring and summer" ,
                   "Wetter winter - dryer summer","Drying - winter peak" )
color_class_labels <- c(
  "Earlier snomelt" = "skyblue",
  "Strong wetting - winter peak" = "forestgreen",
  "Strong drying - winter peak" = "red4",
  "Drying - summer peak" = "orangered",
  "Drying - spring and summer" = "hotpink",
  "Wetter winter - dryer summer" = "mediumpurple",
  "Drying - winter peak" = "darkorange"
)

color_class_plot <- c(
  "1" = "skyblue",
  "2" = "forestgreen",
  "3" = "red4",
  "4" = "orangered",
  "5" = "hotpink",
  "6" = "mediumpurple",
  "7" = "darkorange"
)
#hydrometeo changes

ys2=c(1991:2020)
my=which(!is.na(match(SaveAgg$SaveAgg,ys2)))
Keep=(SaveAgg[my,-1])

mean_snowfracion2=colMeans(Keep)

#match ppoints with different tables

sf2pp=match(ppoints$otlts_x,names(mean_snowfracion2))
ppoints$snowfraction2=mean_snowfracion2[sf2pp]
ppoints$sfchange=(ppoints$snowfraction2-ppoints$snowfraction)/ppoints$snowfraction

plot(ppoints$sfchange)
palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
psf<-ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=sfchange*100),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="relative \nsnow fraction \nchange (%)",
    breaks=seq(-100,100,10), limits=c(-40,40)) +
  #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
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
        legend.key.size = unit(.8, "cm"))+
  ggtitle ("Snow fraction change between 1951-1980 and 1991-2020")

psf

#Precipitation

myp=which(!is.na(match(SavePrecip$year,ys2)))

KeepP=SavePrecip[myp,c(1,2,5)]

KeepPag<- stats::aggregate(list(pmean=KeepP$val.sum), by = list(catch = KeepP$cat), 
                           FUN = function(x) c(pmean=mean(x, na.rm = T)))
KeepPag <- do.call(data.frame, KeepPag)

#Add ET0
my=which(!is.na(match(SaveE0$SaveAgg,ys2)))
KeepE0=(SaveE0[my,-1])
Ke0_means2=colMeans(KeepE0)

kpe=match(KeepPag$catch,cnames)

KeepPag$et2=Ke0_means2[kpe]

KeepPag$AI2=KeepPag$pmean/KeepPag$et2

#Plot of AI
mai=match(ppoints$otlts_x,KeepPag$catch)
ppoints$AI2=KeepPag$AI2[mai]
ppoints$Precip2=KeepPag$pmean[mai]
ppoints$ET02<-KeepPag$et2[mai]

ppoints$chAI=(ppoints$AI2-ppoints$AI1)/ppoints$AI1
ppoints$chET0=(ppoints$ET02-ppoints$ET01)/ppoints$ET01

ppoints$chP=(ppoints$Precip2-ppoints$Precip1)/ppoints$Precip1

plot(ppoints$chP)
plot(ppoints$chET0)
palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
ppr<-ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=chP*100),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="relative \nyearly precipitation \nchange (%)",
    breaks=seq(-100,100,10), limits=c(-20,20)) +
  #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
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
        legend.key.size = unit(.8, "cm"))+
  ggtitle ("Yearly precipitation change between 1951-1980 and 1991-2020")

ppr

pET0<-ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=chET0*100),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="relative ETO \nchange (%)",
    breaks=seq(-100,100,10), limits=c(-20,20)) +
  #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
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
        legend.key.size = unit(.8, "cm"))+
  ggtitle ("ET0 change between 1951-1980 and 1991-2020")


myp=which(!is.na(match(SavePrecip$year,ys2)))

KeepP=SavePrecip[myp,c(1,4,5)]

KeepPag<- stats::aggregate(list(px=KeepP$val.xtr), by = list(catch = KeepP$cat), 
                           FUN = function(x) c(pmean=mean(x, na.rm = T)))
KeepPag <- do.call(data.frame, KeepPag)

mai=match(ppoints$otlts_x,KeepPag$catch)
ppoints$px2=KeepPag$px[mai]

ppoints$chpx<-(ppoints$px2-ppoints$px1)/ppoints$px1
plot(ppoints$chpx)

palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
ppx<-ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=chpx*100),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="relative PrC \nchange (%)",
    breaks=seq(-100,100,10), limits=c(-20,20)) +
  #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
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
        legend.key.size = unit(.8, "cm"))+
  ggtitle ("Precipitation concentration change between 1951-1980 and 1991-2020")


#AET

my=which(!is.na(match(SaveAET$SaveAgg,ys1)))
KeepAET=(SaveAET[my,-1])
KAET_means=colMeans(KeepAET)
km=match(ppoints$otlts_x,cnames)
ppoints$AET1=KAET_means[km]

my=which(!is.na(match(SaveAET$SaveAgg,ys2)))
KeepAET=(SaveAET[my,-1])
KAET_means=colMeans(KeepAET)
ppoints$AET2=KAET_means[km]

ppoints$AETch=(ppoints$AET2-ppoints$AET1)/ppoints$AET1
  
  
palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=AETch*100),color="black",alpha=1,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="relative AET \nchange (%)",
    breaks=seq(-100,100,10), limits=c(-20,20)) +
  #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
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


#runoff fraction


ys1=c(1951:1980)
my=which(!is.na(match(SaveRF$SaveRF,ys1)))
Keep=(SaveRF[my,-1])
mean_runfracion=colMeans(Keep)

#match ppoints with different tables

rf2pp=match(ppoints$otlts_x,names(mean_runfracion))
ppoints$runfraction1=mean_runfracion[rf2pp]


my=which(!is.na(match(SaveRF$SaveRF,ys2)))
Keep=(SaveRF[my,-1])
mean_runfracion=colMeans(Keep)

#match ppoints with different tables

ppoints$runfraction2=mean_runfracion[rf2pp]
ppoints$rfchange=(ppoints$runfraction2-ppoints$runfraction1)/ppoints$runfraction1
plot(ppoints$rfchange)


palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
prf<-ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=rfchange*100),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="relative Runoff fraction \nchange (%)",
    breaks=seq(-100,100,10), limits=c(-20,20)) +
  #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = .5, barheight = 20,reverse=F),
         size= "none")+
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
        legend.key.size = unit(.8, "cm"))+
  ggtitle ("Runoff fraction change between 1951-1980 and 1991-2020")


library(gridExtra)
combined_plot <- grid.arrange(psf, ppr, pET0, ppx, prf, ncol = 3, nrow = 2)

# Save the arranged plots to a file
ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/Shifts_combined_map2.png", plot = combined_plot, width = 30, height = 20)


#now bring all my variables together

#variables to keep
cpp=colnames(ppoints)
vkeep=c("HYBAS_ID","outlets_x","upa","Cluster","Cshift","chP","chET0","chpx","rfchange","sfchange")
mvk=match(vkeep,cpp)
Xvars=ppoints[,mvk]

mai=match(ppoints$otlts_x,SaveVars$outlets.x)
oflof<-SaveVars[mai,c(4,5,6,9,16,21)]
colnames(oflof)[c(1,5,6)]<-c("forestchange","water_demandchange","reservoir_ratiochange")
Xvars=cbind(Xvars,oflof)
st_geometry(Xvars)<-NULL

ppoints2<-cbind(ppoints,oflof)


XpVars=Xvars[,c(6:16)]

plot(XpVars$chP,XpVars$chpx)
library(corrplot)
res <- cor(XpVars,method="pearson")
rs <- cor.mtest(XpVars, method="pearson", conf.level=0.95)
round(res, 2)
corrplot(res, type = "upper", method="ellipse",
         tl.col = "black", tl.srt = 45, pch.cex = 0.9)$corrPos->p1

text(p1$x, p1$y, round(p1$corr, 2))


#variable to be predicted


# RegimeChange=rep(0,length(Xvars$HYBAS_I))
# RegimeChange[which(Xvars$C1name!=as.character(Xvars$C0name))]=1


XpVars$RegimeChange=factor(ppoints$Cshift)


#Histogram of explaining variables by class

color_class_labels <- c(
  "Earlier snomelt" = "skyblue",
  "Strong wetting - winter peak" = "forestgreen",
  "Strong drying - winter peak" = "red4",
  "Drying - summer peak" = "orangered",
  "Drying - spring and summer" = "hotpink",
  "Wetter winter - dryer summer" = "mediumpurple",
  "Drying - winter peak" = "darkorange"
)

p1<-ggplot(ppoints2, aes(y=factor(Cshift), x=sfchange)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_x_continuous(limits = c(-1,1),name="snow fraction change",breaks = seq(-1,1,by=0.2),minor_breaks = seq(-1,1,0.1))+
  scale_y_discrete(breaks=class_labels, name="Regime class")+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_class_labels,name="Shifts")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p1

p2<-ggplot(ppoints2, aes(y=factor(Cshift), x=chP)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_x_continuous(limits = c(-0.5,0.5),name="total precipitation change",breaks = seq(-1,1,by=0.2),minor_breaks = seq(-1,1,0.1))+
  scale_y_discrete(breaks=class_labels, name="Regime class")+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_class_labels,name="Shifts")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p2


p3<-ggplot(ppoints2, aes(y=factor(Cshift), x=chET0)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_x_continuous(limits = c(-0.5,0.5),name="ET0 change",breaks = seq(-1,1,by=0.2),minor_breaks = seq(-1,1,0.1))+
  scale_y_discrete(breaks=class_labels, name="Regime class")+
  coord_cartesian(xlim=c(-.2,.2))+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_class_labels,name="Shifts")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p3

p4<-ggplot(ppoints2, aes(y=factor(Cshift), x=chpx)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_x_continuous(limits = c(-0.5,.5),name="precipitation concentration change",breaks = seq(-1,1,by=0.2),minor_breaks = seq(-1,1,0.1))+
  scale_y_discrete(breaks=class_labels, name="Regime class")+
  coord_cartesian(xlim=c(-.3,.3))+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_class_labels,name="Shifts")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p4

p5<-ggplot(ppoints2, aes(y=factor(Cshift), x=rfchange)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_x_continuous(limits = c(-1,1),name="direct runoff fraction change",breaks = seq(-1,1,by=0.2),minor_breaks = seq(-1,1,0.1))+
  scale_y_discrete(breaks=class_labels, name="Regime class")+
  coord_cartesian(xlim=c(-.5,.5))+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_class_labels,name="Shifts")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p5



ppoints2$reservoir_ratiochange[which(ppoints2$reservoir_ratiochange==0)]=1e-3
p6<-ggplot(ppoints2, aes(y=factor(Cshift), x=reservoir_ratiochange)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_x_continuous(limits = c(1e-3,5),name="reservoir concentration change",breaks = c(0.01,seq(0,.1,by=.05),.2,.5,1,2,5),minor_breaks = seq(-1,10,0.1), trans="log")+
  scale_y_discrete(breaks=class_labels, name="Regime class")+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_class_labels,name="Shifts")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p6

p7<-ggplot(ppoints2, aes(y=factor(Cshift), x=forestchange)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_x_continuous(limits = c(-.5,.5),name="forest fraction change",breaks = seq(-1,10,by=.1),minor_breaks = seq(-1,10,0.1))+
  scale_y_discrete(breaks=class_labels, name="Regime class")+
  coord_cartesian(xlim=c(-.1,.3))+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_class_labels,name="Shifts")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p7


p8<-ggplot(ppoints2, aes(y=factor(Cshift), x=sealedchange)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_x_continuous(limits = c(-1,1),name="sealed fraction change",breaks = seq(-1,10,by=.05),minor_breaks = seq(-1,10,0.1))+
  scale_y_discrete(breaks=class_labels, name="Regime class")+
  coord_cartesian(xlim=c(-0,.1))+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_class_labels,name="Shifts")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p8

p9<-ggplot(ppoints2, aes(y=factor(Cshift), x=water_demandchange)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_x_continuous(limits = c(-5,5),name="water demand change",breaks = seq(-1,10,by=.2),minor_breaks = seq(-1,10,0.1))+
  scale_y_discrete(breaks=class_labels, name="Regime class")+
  coord_cartesian(xlim=c(-.5,.5))+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_class_labels,name="Shifts")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p9


p10<-ggplot(ppoints2, aes(y=factor(Cshift), x=irrigatedchange)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=Cshift),linewidth=0.8,outlier.alpha = 0.4)+
  scale_x_continuous(limits = c(-5,5),name="irrigated change",breaks = seq(-1,10,by=.05),minor_breaks = seq(-1,10,0.1))+
  scale_y_discrete(breaks=class_labels, name="Regime class")+
  coord_cartesian(xlim=c(-.1,.1))+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  scale_fill_manual(values=color_class_labels,name="Shifts")+
  # geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p10
# Arrange the plots into a grid layout
library(gridExtra)
combined_plot <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,p10, ncol = 3, nrow = 4)

# Save the arranged plots to a file
ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/Shifts_combined_plot4.png", plot = combined_plot, width = 30, height = 30)





#Random forest model ----



# random forest
library(randomForest)


rf_data=XpVars
# Assume df is your data frame and 'class' is your target variable
# Split data into training and test sets
set.seed(123)  # For reproducibility
train_indices <- sample(1:nrow(rf_data), size = 0.7 * nrow(rf_data))
train_data <- rf_data[train_indices, ]
test_data <- rf_data[-train_indices, ]

# Train a random forest model
# Assuming 'class' is the column name of the target variable
rf_model <- randomForest(RegimeChange ~ ., data = train_data, ntree = 400, mtry = 5, importance = TRUE)

# Print the model summary
print(rf_model)

# Predict class membership for the test data
predictions <- predict(rf_model, newdata = test_data)

# Evaluate model performance
# Create a confusion matrix
confusion_matrix <- table(test_data$RegimeChange, predictions)
print(confusion_matrix)

# Calculate accuracy
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(paste("Accuracy:", accuracy))

# Predict class membership for all the data
predictionsAll <- predict(rf_model, newdata = rf_data)




#Improve thr algo, later
library(caret)

# Set up cross-validation parameters
control <- trainControl(method = "cv", number = 10)

param_grid=expand.grid(mtry = c(1:4),
                       splitrule = c( "variance", "extratrees"),
                       min.node.size = c(1, 5))
# Train the model using caret's train function
# This uses randomForest as the underlying method
set.seed(123)  # For reproducibility
rf_model_cv <- train(RegimeChange ~ ., 
                     data = train_data, 
                     method = "ranger",
                     trControl = control,
                     tuneGrid = expand.grid(mtry = 5),  # Adjust mtry as needed
                     ntree = 500)

# Print the results of cross-validation
print(rf_model_cv)

# Access the best model's accuracy
best_accuracy <- max(rf_model_cv$results$Accuracy)
print(paste("Best Accuracy from 10-fold CV:", best_accuracy))

predictiontest <- predict(rf_model_cv, newdata = test_data)
print(predictiontest)

confusion_matrixT <- table(test_data$RegimeChange, predictiontest)
print(confusion_matrixT)

# Calculate accuracy
accuracy <- sum(diag(confusion_matrixT)) / sum(confusion_matrixT)
print(paste("Accuracy:", accuracy))

predictionFinal<-predict(rf_model_cv, newdata = rf_data)
print(predictionFinal)


ppoints$predchange<-predictionFinal
ppoints$Cshift
class_labels <- c( "Earlier snomelt", "Strong wetting - winter peak", "Strong drying - winter peak","Drying - summer peak","Drying - spring and summer" ,
                   "Wetter winter - dryer summer","Drying - winter peak" )


color_class_labels <- c(
  "Earlier snomelt" = "skyblue",
  "Strong wetting - winter peak" = "forestgreen",
  "Strong drying - winter peak" = "red4",
  "Drying - summer peak" = "orangered",
  "Drying - spring and summer" = "hotpink",
  "Wetter winter - dryer summer" = "mediumpurple",
  "Drying - winter peak" = "darkorange"
)
color_vector <- sapply(class_labels, function(x) color_class_labels[x])
print(color_vector)

color_class_plot <- c(
  "1" = "skyblue",
  "2" = "forestgreen",
  "3" = "red4",
  "4" = "orangered",
  "5" = "hotpink",
  "6" = "mediumpurple",
  "7" = "darkorange"
)

ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=factor(predchange)),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  # scale_fill_gradientn(
  #   colors=palet,oob = scales::squish, name="relative Runoff fraction \nchange (%)",
  #   breaks=seq(-100,100,10), limits=c(-20,20)) +
  scale_fill_manual(values=color_class_labels)+
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
        legend.key.size = unit(.8, "cm"))+
  ggtitle ("")

#ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/Shifts_classes_30y1d.jpg"), width=30, height=20, units=c("cm"),dpi=300) 

vi1<-varImp(rf_model_cv, scale = T)
plot(vi1)



library(ggplot2)
library(caret)
library(kernelshap)
library(shapviz)
library(xgboost)
#library(randomForest)

# library(DALEX)
library(shapper)



library(reticulate)
conda_create(envname = "r-shap", packages = c("shap", "python=3.6.8"))   
use_condaenv("r-shap", required = TRUE)

# # Take a subsample of the data for bg_X if needed
# bg_X <- rf_data[, -which(names(rf_data) == "RegimeChange")]
# if (nrow(bg_X) > 50) {
#   set.seed(123)  # For reproducibility
#   bg_X <- bg_X[sample(nrow(bg_X), 50), ]
# }
# 
# 
# 
# #needs shapper and python
# ive_rf <- individual_variable_effect(rf_model_cv, data = bg_X, predict_function = predict_function,
#                                      new_observation = bg_X[1:30,], nsamples = 50)
# ive_rf
# 
# 
# pl3 <- plot(ive_rf, bar_width = 4, show_predicted = T,
#             cols = c("id","ylevel"), rows = "label")
# 
# print(pl3)
# 
# 
# plot(ive_rf)


# 
# # Verify Python configuration
py_config()
sys <- import("sys", convert = TRUE); sys$path
paste0(system.file(package = "reticulate"),"/python")
data(iris)

# Prepare the data
X_train <- rf_data[, -which(names(rf_data) == "RegimeChange")]
y_train <- rf_data[, which(names(rf_data) == "RegimeChange")]

# Train a Random Forest model
rf_model <- randomForest(x = X_train, y = y_train, ntree = 100)

# Create a DALEX explainer for the Random Forest model
explainer_rf <- explain(rf_model_cv, data = X_train, y = y_train, label = "Random Forest model")

Xsample<-X_train[1,]
# Calculate SHAP values with the shapper package
waziz=shap.sample(X_train, K=10)
shap_values_rf <- shap(explainer_rf, Xsample)
# Plot SHAP values for the first observation
plot(shap_values_rf)



library(ranger)
set.seed(1)


param_grid=expand.grid(mtry = c(5,10),
                       splitrule = c( "gini"),
                       min.node.size = c(10))
# Train the model using caret's train function
# This uses randomForest as the underlying method
set.seed(123)  # For reproducibility
rf_model_cv <- train(RegimeChange ~ ., 
                     data = train_data, 
                     method = "ranger",
                     trControl = control,
                     importance = "impurity",
                     tuneGrid = param_grid,  # Adjust mtry as needed
                     ntree = 500)

plot(rf_model_cv)


fit <- ranger(RegimeChange ~ ., 
              data = rf_data, 
              mtry = 5,  # Adjust mtry as needed
              num.trees = 500,
              importance = "impurity",
              
              probability = TRUE)
fit$prediction.error
bg_X <- rf_data
if (nrow(bg_X) > 500) {
  set.seed(123)  # For reproducibility
  bg_X <- bg_X[sample(nrow(bg_X), 500), ]
}

ps_prob <- permshap(fit, X = bg_X[,-12],bg_n=50) |> 
  shapviz()



sv_importance(ps_prob)



# s <- kernelshap(fit, bg_X[,-12], bg_X = bg_X)
# 
# # Compute SHAP values using kernelshap
# s <- kernelshap(
#   object = rf_model_cv$finalModel,  # Access the final model from caret
#   X = bg_X,
#   pred_fun = predict_function,
#   bg_N = 1
# )
# # Step 2: Turn them into a shapviz object
# sv <- shapviz(s)

# Step 3: Gain insights...


color_class_plot <- c(
  "Class_1" = "skyblue",
  "Class_2" = "forestgreen",
  "Class_3" = "red4",
  "Class_4" = "orangered",
  "Class_5" = "hotpink",
  "Class_6" = "mediumpurple",
  "Class_7" = "darkorange"
)
# For the 'mshapviz' object (all classes), sv_importance can also plot:
# This will show importance for all three classes, typically in faceted plots.
sv_importance(ps_prob)+
  scale_fill_manual(values=(color_class_labels))+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.x = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("shap values from random forest model")
ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/SHAP_prediction_rf_30y1d_v2.jpg"), width=30, height=20, units=c("cm"),dpi=300) 



sv_dependence(ps_prob, "sealedchange")



#xgboost ----

library(xgboost)

# Example assuming 'your_data' is your dataframe and 'target' is the target variable
rf_data<-rf_data[order(as.numeric(rf_data$RegimeChange)),]
target <- as.integer(rf_data$RegimeChange) - 1  # Convert factor to integer starting at 0
data_matrix <- model.matrix(~ . + 0, data = rf_data[, -which(names(rf_data) == "RegimeChange")])

dtrain <- xgb.DMatrix(data = data_matrix, label = target)

num_classes <- length(unique(target))
params <- list(
  objective = "multi:softprob",  # Use "multi:softprob" for probability outputs
  eval_metric = "mlogloss",      # Multi-class log loss
  target=target,
  data=data_matrix,
  num_class = num_classes,       # Number of classes
  max_depth = 6,                 # Maximum depth of a tree
  eta = 0.3,                     # Learning rate
  nthread = 2                    # Number of parallel threads
)


nrounds <- 100  # Number of boosting rounds
xgb_model <- xgb.train(params = params, data = dtrain, nrounds = nrounds)


# Create test data matrix in the same way as the training data
dtest <- xgb.DMatrix(data = data_matrix)  # Use your actual test data matrix if available

# Predict
predictions <- predict(xgb_model, newdata = dtest)

# Reshape predictions to a matrix with num_classes columns
pred_matrix <- matrix(predictions, ncol = num_classes, byrow = TRUE)

# Get the predicted class for each instance
predicted_classes <- max.col(pred_matrix) - 1

# Example of confusion matrix
actual <- target  # Replace with your actual target variable for the test data if available
table(Predicted = predicted_classes, Actual = actual)

# Assuming `xgb_model` is your trained model
importance_matrix <- xgb.importance(feature_names = colnames(data_matrix), model = xgb_model)

# Print the feature importance
print(importance_matrix)

# Visualize the feature importance
xgb.plot.importance(importance_matrix)



data1 = rf_data[, -which(names(rf_data) == "RegimeChange")]

sv_xgb_all_classes <- shapviz(xgb_model, X_pred = data_matrix, X = data1)
# To get a shapviz object for a specific class (e.g., the 3rd class, 'virginica' in iris):
sv_xgb_class_1 <- shapviz(xgb_model,  X_pred = data_matrix, X = data1, which_class = 1)


# Determine which class corresponds to the index
target_class_index <- which(class_levels == "1") - 1 

sv_importance(sv_xgb_class_1,"beeswarm",show_numbers = TRUE)

#ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/SHAP_prediction_xgboost_merde_30y1d.jpg"), width=30, height=20, units=c("cm"),dpi=300) 


color_class_plot <- c(
  "Class_1" = "skyblue",
  "Class_2" = "forestgreen",
  "Class_3" = "red4",
  "Class_4" = "orangered",
  "Class_5" = "hotpink",
  "Class_6" = "mediumpurple",
  "Class_7" = "darkorange"
)
# For the 'mshapviz' object (all classes), sv_importance can also plot:
# This will show importance for all three classes, typically in faceted plots.

sv_xgb_all<-sv_xgb_all_classes


rlevels= c("Class_5", "Class_4", "Class_7","Class_1", "Class_3", "Class_2","Class_6")
sv_importance(sv_xgb_all)+
  scale_fill_manual(values=(color_class_plot),labels=rev(class_labels))+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.x = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("shap values from xgboost model")


z=get_shap_values(sv_xgb_all)
.get_imp <- function(z, sort_features = TRUE) {
  if (is.matrix(z)) {
    imp <- colMeans(abs(z))
    if (sort_features) {
      imp <- sort(imp, decreasing = TRUE)
    }
    return(imp)
  }
  # list/mshapviz
  imp <- sapply(z, function(x) colMeans(abs(x)))
  if (sort_features) {
    imp <- imp[order(-rowSums(imp)), ]
  }
  return(imp)
}

imp <- .get_imp(get_shap_values(sv_xgb_all), sort_features = T)
ord <- rownames(imp)

imp_df <- data.frame(feature = factor(ord, rev(ord)), value = imp)
colnames(imp_df) = class_labels

long <- melt(imp_df, id.vars = c("feature"), variable.name = "value")

p <- ggplot2::ggplot(long, ggplot2::aes(x = value, y = feature, fill=variable)) +
  ggplot2::geom_bar( width = 0.4, stat = "identity", position="dodge") +
  ggplot2::labs(x = "mean(|SHAP value|)", y = ggplot2::element_blank())

p

ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/SHAP_prediction_xgboost_30y1d_v2.jpg"), width=30, height=20, units=c("cm"),dpi=300) 

# B. Local Explanation: SHAP Dependence Plot using sv_dependence()
# Shows SHAP values of a feature against its feature values, can be colored by another feature.
# This is for the 'virginica' class.
sv_dependence(sv_xgb_iris_1, v = "sfchange")


# For an 'mshapviz' object, sv_dependence can create multiple plots:
# This will show dependence for 'Petal.Length' across all three species.
sv_dependence(sv_xgb_iris_all_classes, v = "sfchange")

# C. Local Explanation: SHAP Waterfall Plot using sv_waterfall()
# Explains a single prediction incrementally.
# Let's explain the first observation for the 'virginica' class (which_class = 3)
sv_waterfall(sv_xgb_iris_1, row_id = 1)

# D. Local Explanation: SHAP Force Plot using sv_force()
# Creates a force plot for one observation.
sv_force(sv_xgb_class_1)



shap_values <- shap.values(xgb_model = model, X_train = data1)
shap_values$mean_shap_score
shap_values_class <- shap_values$shap_score

shap_long <- shap.prep(xgb_model = model, X_train = data_matrix)
# To prepare the long-format data:
# **SHAP summary plot**
shap.plot.summary(shap_long, scientific = T, kind="bar")
shap.plot.summary(shap_long, x_bound  = 1.5, dilute = 10)


# set 6 clustering groups of observations.  
plot_data <- shap.prep.stack.data(shap_contrib = shap_values$shap_score)
# you may choose to zoom in at a location, and set y-axis limit using `y_parent_limit`  
shap.plot.force_plot(plot_data,zoom_in_location = 500)

shap.plot.force_plot_bygroup(plot_data)


sv_obj_class3 <- shapviz(explainer_rf, which_class = 3)
# Compute SHAP values
mer=as.data.frame(data_matrix)
x.interest=mer[1000,]
shapley <- Shapley$new(predictor,x.interest = x.interest)

# Plot SHAP values for an instance, e.g., the first instance
plot(shapley, instance = 1)  # Change instance number as needed

# SHAP summary plot for global feature impact
plot(shapley)

# Dependence plot for a specific feature
plot(shapley, feature = "feature_name")


# Alternatively, compute SHAP values for all instances
shapley_values <- shapley$results

# Assuming you have a SHAP matrix and corresponding feature names
shap_values <- predict(xgb_model, newdata = dtest, predcontrib = TRUE)


# Compute SHAP values using shapviz
sv_rf <- shapviz(rf_model_cv, X = data1)  # Exclude the target variable from X


ive_rf <- individual_variable_effect(xgb_model, data = X_train, predict_function = p_function,
                                     new_observation = X_train[1:10,], nsamples = 50)
ive_rf


pl3 <- plot(ive_rf, bar_width = 4, show_attributions = T,show_predicted=F,
            cols = c("id","ylevel"), rows = "label")

pl3

# Define a prediction function

# Define a custom prediction function for the XGBoost model
# Custom prediction function for probabilities
predict_function <- function(model, newdata) {
  newdata_matrix <- xgb.DMatrix(data = as.matrix(newdata))
  preds <- predict(model, newdata = newdata_matrix, reshape = TRUE)
  # Convert probabilities to class labels if needed
  class_preds <- max.col(preds) - 1
  return(class_preds)
}







# Compute SHAP values
X_explain <- data_matrix[1:10, ]
shap_values <- kernelshap(pred_function=predict_function, X_explain , X = data_matrix, model = xgb_model)
# Compute SHAP values using kernelshap
s <- kernelshap(
  object = rf_model_cv$finalModel,  # Access the final model from caret
  X = rf_data[, -which(names(rf_data) == "RegimeChange")],
  pred_fun = predict_function,
  bg_X = bg_X
)

# Check the SHAP values
print(s)

sv <- shapviz(s)

# Check the SHAP values
print(s)

sv <- shapviz(s)
sv_importance(sv)




merd=rf_data[ , which(names(rf_data) != "RegimeChange")]
# Create a predictor object using the iml package
predictor <- Predictor$new(
  model = rf_model_cv, 
  data = merd, 
  y = rf_data$RegimeChange,
  type="prob"
)
predictor$model$t

if (!("sfchange" %in% names(rf_data))) {
  stop("The variable 'sfchange' is not found in your dataset.")
}

imp <- FeatureImp$new(predictor, loss = "ce")
plot(imp)
# Compute SHAP values

interact <- Interaction$new(predictor)
plot(interact)




# Predict class probabilities on the same dataset or a new dataset
probabilities <- predict(rf_model_cv, newdata = rf_data, type = "prob")

# View the probabilities
head(probabilities)


x.interest=merd[1000,]
shap_values <- Shapley$new(predictor,x.interest = x.interest)


# Plot SHAP values for the first instance in your dataset
plot(shap_values, instance = 1)

# To view SHAP values for all instances
shap_values_data <- shap_values$results
print(shap_values_data)




