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

#load average temperature at catchment level -----
tavg<-read.csv(paste0(hydroDir,"/tss/tAvgUpsv2_1950_2020.csv"))
#reoder tavg according to time
time=as.POSIXct(tavg$X)
tavg=tavg[order(time),]
timeStampx=time[order(time)]
tstamp=data.frame(timeStampx)

cname=colnames(tavg)
cname[1]="time"
cname=gsub("X","",cname)
colnames(tavg)=cname
tavg$tday=as.Date(tavg$time)

#loop on each pixel to create a 30-daysrolling verag temperature file -----
tid=seq(1,length(timeStampx))
dsel=hour(timeStampx)
tsday=as.Date(timeStampx[which(dsel==12 | dsel==13)])

zdays=matrix(-9999,ncol=length(cname[-1]),nrow=length(timeStampx))
#matrix for data aggregated to daily time scale
mavg=matrix(-9999,ncol=3,nrow=length(tsday)*length(cname[-1]))

id=0
for (cn in cname[-1])
{
  id=id+1
  col=which(colnames(tavg)==cn)
  print(col)
  tavg_out=tavg[,col]

  #30-days running average
  tavg_mo2=round(tsEvaNanRunningMean(tavg_out,120),2)
  points(tavg_mo2[1:4800],col=4,type="o")
  
  #extract daily values
  dataD=tavg_mo2[which(dsel==12 | dsel==13)]
  zdt=tid[which(dsel==12 | dsel==13)]
  lo=length(dataD)
  zdn=as.numeric(rep(cn,length(zdt)))
  td=tid
  zddf=cbind(zdt,zdn,dataD)
  mavg[c(((id-1)*lo+1):(id*lo)),]=zddf
  zdays[,id]=tavg_mo2
}

#zdays in the 30-days averaged temperature in the considered catchments
zdays=data.frame(timeStampx, zdays)
zdays=zdays[,-1]
names(zdays) = cname

mavg=as.data.frame(mavg)
mavg$time=as.Date(timeStampx[mavg$V1])
mavg$year=year(mavg$time)

AvgFrostAgg=aggregate(list(temp=mavg$V3),
                      by = list(outlet=mavg$V2,yr=mavg$year),
                      FUN = function(x) c(frost=length(which(x<0)),len=length(x),mean=mean(x)))
AvgFrostAgg <- do.call(data.frame, AvgFrostAgg)

#trend on number of catchments experiencin frost per year
AvgFrostAgg_trend=aggregate(list(frost=AvgFrostAgg$temp.frost,Temp=AvgFrostAgg$temp.mean),
                        by = list(outlets=AvgFrostAgg$yr),
                        FUN = function(x) c(len=length(which(x>0)),mean=mean(x)))
AvgFrostAgg_trend <- do.call(data.frame, AvgFrostAgg_trend)

#trend plot
plot(AvgFrostAgg_trend$outlets,AvgFrostAgg_trend$frost.len,type="l")
abline(lm(frost.len ~ outlets, data=AvgFrostAgg_trend),col="red")


AvgFAgg=aggregate(list(tf=AvgFrostAgg$temp.frost),
                      by = list(outlet=AvgFrostAgg$outlet),
                      FUN = function(x) c(frost=length(which(x>0)),len=length(x),mean=mean(x)))
AvgFAgg <- do.call(data.frame, AvgFAgg)

names(AvgFAgg)[1]="outlets"
AvgFAgg$frost=AvgFAgg$tf.mean
#AvgFAgg$frost[which(AvgFAgg$tf.frost<50)]=0



#Now I have to link catchment which have not enough frost seasons 
noFrost=AvgFAgg$outlets[which(AvgFAgg$frost==0)]
colrep=match(noFrost,cname)

#remove the columns without frost
frostcat=zdays[,-colrep]
#save 30-days averaged temperature at 6-hourly timestep resolution
#save(frostcat, file=paste0(hydroDir,"/Drought/catchment_frost.Rdata"))

#save 30-days averaged temperature at daily timestep resolution
#save(mavg, file=paste0(hydroDir,"/Drought/catchment_frost_daily.Rdata"))



### [Plot] - Figure S5 - Average number of frost days per year over the period 1951-2020 ----
countries <- read_sf(dsn = paste0(hydroDir,"/Countries/Countries_EFAS.shp"))
catchEFAS <- st_intersection(cst7, countries)
Frostplot=inner_join(AvgFAgg, catchEFAS,by="outlets")

ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=frost,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet, limits=c(0,182),
    oob = scales::squish,na.value="transparent", name="days")   +
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("Average number of frost days per year")


ggsave("D:/tilloal/Documents/LFRuns_utils/Figures/FrostdaysPerYr.jpg", width=30, height=21, units=c("cm"),dpi=1000)

outlets=unique(AvgFrostAgg$outlet)
plot(AvgFrostAgg$yr[which(AvgFrostAgg$outlet==1029)],AvgFrostAgg$temp.frost[which(AvgFrostAgg$outlet==1029)],type="l")
#

mavg=as.data.frame(mavg)
mavg$time=timeStampx[mavg$V1]
mavg$year=year(mavg$time)
AvgFrostAgg=aggregate(list(temp=mavg$V3),
                   by = list(outlet=mavg$V2),
                   FUN = function(x) c(frost=length(which(x<0)),len=length(x),mean=mean(x)))
AvgFrostAgg <- do.call(data.frame, AvgFrostAgg)

AvgFrostAgg$ratpix=(AvgFrostAgg$temp.frost/AvgFrostAgg$temp.len)*100


names(AvgFrostAgg)[1]="outlets"
Frostplot=inner_join(AvgFrostAgg, catchEFAS,by="outlets")

### [Plot] - Figure S5 supplement - number of frost days (%) the period 1951-2020 ----
ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=ratpix,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradientn(
    colors=palet, limits=c(0,50),
    oob = scales::squish,na.value="grey", name="frost days (%)")   +
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


### [Plot] - Figure S5 - Mean daily temperature for the period 1951-2020 ----
palet=c(hcl.colors(9, palette = "YlOrRd", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="gray85")+
  geom_sf(data=Frostplot,aes(fill=temp.mean,geometry=geometry),color="transparent")+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  scale_fill_gradient2(
    low = "royalblue",
    mid = "white",
    high = "darkorange",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    transform = "identity",
    guide = "colourbar",
    aesthetics = "fill",
    limits=c(-5,20),
    oob = scales::squish, 
    name="mean 30 days \naveraged temperature"
  )
  labs(x="Longitude", y = "Latitude")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))