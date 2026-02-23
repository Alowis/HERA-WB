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

#Upstream catchments
UpCat<-st_read(dsn="D:/tilloal/Documents/01_Projects/RegimeShifts/Upstreamgroups.shp")
upup<-UpCat[which(UpCat$up=="Upstream Catchments"),]
hybas2uc<-match(upup$outlets_x,outhybas$outlets.y)
outhybas<-outhybas[hybas2uc,]

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


#Now compare the upstream area with the hybas shapefile area

CatUpA=inner_join(Catf7,UpArea,by=c("llcoord"="latlong"))
matcat=match(outhybas$latlong,CatUpA$llcoord)
CatUpA=CatUpA[matcat,]
plot(CatUpA$SUB_AREA,CatUpA$upa)

median(CatUpA$upa)
hist(CatUpA$upa)
library(scales)
pu<-ggplot(CatUpA, aes(x=upa)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=15,alpha=0.9,lwd=1)+
  scale_y_continuous(breaks=seq(0,600, by=100),name="Number of stations")+
  scale_x_log10(name=expression(paste("Upstream area ", (km^2),sep = " ")),
                breaks=c(100,1000,10000,100000), minor_breaks = log10_minor_break(),
                labels=c("100","1 000","10 000","100 000")) +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=16),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  annotate("label", x=450000, y=500, label= paste0("n = ",length(CatUpA$upa)),size=6)

pu
#ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/histo_outlets.jpg", pu, width=20, height=15, units=c("cm"),dpi=1500)


#separate catchment outlets from big rivers

ratUp=CatUpA$SUB_AREA/CatUpA$upa

plot(ratUp, log="y")
abline(h=1.5)
abline(h=0.5)

CatUpA$outlet=1
CatUpA$outlet[which(ratUp>1.5)]=2
CatUpA$outlet[which(ratUp<0.5)]=3
CatUpA$outlet[which(CatUpA$outlet==3 & CatUpA$upa<1e4)]=2

CatLR=CatUpA[which(CatUpA$outlet==3),]
st_geometry(CatLR)=NULL
CatO=CatUpA[-which(CatUpA$outlet==3),]

# 
# CatLR=CatUpA[which(CatUpA$upa>2e4),]
# st_geometry(CatLR)=NULL
# CatO=CatUpA[-which(CatUpA$upa>2e4),]

ppl <- st_as_sf(CatLR, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)


ppc <- st_as_sf(CatO, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppc <- st_transform(ppc, crs = 3035)

# Compute centroids
CatO_centroids <- st_centroid(ppc)

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,size=upa),color="transparent",alpha=.8,shape=21,stroke=0, fill="darkblue")+ 
  geom_sf(data=CatO_centroids,aes(geometry=geometry,size=upa,fill=factor(outlet)),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  # scale_fill_gradientn(
  #   colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
  #   breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  # guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5,reverse=F),
  #        size= guide_legend(override.aes = list(fill = "grey50")))+
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


#ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/catclass3.jpg"), width=20, height=20, units=c("cm"),dpi=1000) 
#load data for each catchmemt at catchment level -----
#tavg<-read.csv(paste0(hydroDir,"/tss/tAvgUpsv2_1950_2020.csv"))
#do this for every variable
#the idea is to open several variable and disaggregate the flow

library(data.table)

Q <- fread(paste0(hydroDir,"/tss/HERA_Histo/disWin_1951_2020.csv"),header=TRUE)

#convert Q to specific dicharge 

#reoder according to time
time=as.POSIXct(Q$V1)
Q=Q[order(time),]
timeStampx=time[order(time)]
tstamp=data.frame(timeStampx)
#convert to daily values
dt = difftime(timeStampx[2],timeStampx[1],units="days")
dt= as.numeric(dt)

#remove Q obs that are outside
matcat=match(outhybas$latlong,UpArea$latlong)
UpArea=UpArea[matcat,]

matcol=match(UpArea$outlets,as.numeric(colnames(Q)))
cnames=names(Q)[matcol]
plot(matcol)

Q2 <- Q[, .SD, .SDcols = matcol]


rm(Q)
gc()
#convert to daily values

Ticino=Q[["290666"]]
Ticino=data.frame(timeStampx,Ticino)
cupa=UpArea$upa[which(UpArea$outlets==290666)]
#Use specific discharge
#Ticino$Qs=Ticino$Ticino/(cupa*1000000)*1000*3600*24
Ticino$Qs=Ticino$Ticino
Ticido=Ticino[,c(1,3)]
# Ziz=plotRegime(Ticido,"Ticino")




RegimeFast=function(data){
  names(data)=c("date","Q")
  deb=data$date[4]
  mois.deb <-  seq(as.Date(deb), by="month", length=12)
  
  period=paste0(format(range(data$date)[1],"%b %Y"),"-",format(range(data$date)[2],"%b %Y"))
  jours <- as.numeric(format(data$date,"%j"))
  ## Iddinces des lignes de meme categorie
  ind.j <- tapply(seq(length(jours)), jours, c)
  ind.j <- ind.j[-366]
  ## Calcul
  Qc <- data.frame(date=as.numeric(names(ind.j)),
                   mean=sapply(ind.j, function(x) mean(data$Q[x], na.rm=TRUE)))
  return(Qc)
}

library(data.table)
library(lubridate)
library(fda) # For Fourier basis

process_column <- function(column_name, Q,tstamp, UpArea, dt, nbasis=6, yearstart, wsize) {
  
  
  #ly <- length(years)
  yr=yearstart
  yrwindow <- yr:(yr + wsize)
  subs <- which(!is.na(match(year(as.Date(tstamp)), yrwindow)))
  
  #upstream area
  cupa=UpArea$upa[which(UpArea$outlets==as.numeric(column_name))]
  
  # Compute Qs
  data=Q[[column_name]][subs]
  data=data.frame(tstamp[subs],data)
  data$Qs <- data[,2]
  
  #specific discharge
  # data$Qs <- data[,2]/ (cupa * 1000000) * 1000 * 3600 * 24

  
  # Apply running mean
  data$daily <- tsEvaNanRunningMean(data$Qs, windowSize = 1 / dt)
  
  # Rename columns
  names(data)[1] <- "date"
  
  # Filter for specific hours
  dsel <- hour(data$date)
  datD <- data[(dsel == 12 | dsel == 13), c(1,4)]
  
  # Regime over the first 20 years
  # years <- unique(year(dataD$date))
  # wsize <- 19
  # ly <- length(years)
  # yrwindow <- yr:(yr + wsize)
  # subs <- which(!is.na(match(year(dataD$date), yrwindow)))
  # datD <- dataD[subs, ]
  #sd specific discharge
  
  qsd=sd(datD[,2])/ (cupa * 1000000) * 1000 * 3600 * 24
  
  #mean specific discharge
  qm <- mean(datD[,2])/ (cupa * 1000000) * 1000 * 3600 * 24
  
  # Apply running mean again
  datD[,2]=tsEvaNanRunningMean(datD[,2],31)
  
  #coefficient mensuel de debit
  datD[,2]=datD[,2]/(mean(datD[,2],na.rm=T))
  
  #specific discharge
  #datD[,2]=datD[,2]/(cupa * 1000000) * 1000 * 3600 * 24
  
  Regi <- RegimeFast(datD)
  #normalisation to really cluster the shape and not the values
  #Regi$mean=Regi$mean/max(Regi$mean)
    
  # Fit Fourier
  value_data <- c(Regi$mean)
  rangeval <- c(min(Regi$date), max(Regi$date))
  basis <- create.fourier.basis(rangeval, nbasis)
  
  fd_object <- Data2fd(argvals = Regi$date, y = value_data, basisobj = basis)
  # summary(fd_object)
  #plot(fd_object, main = paste("Functional Data Representation -", column_name))
  
  cv=qsd/qm
  # Save parameters
  params <- c(as.numeric(column_name), fd_object$coefs,qm,qsd)
  regimeS<-Regi
  rlist=list(params,regimeS)
  return(rlist)
}


#different clistering option:
#nbasis=4,5,6,7
#adding cv
#adding mean specific discharge and sd qsp

# results <- lapply(names(Q2)[-1], function(col_name) {
#   process_column(col_name, Q2,timeStampx, UpArea, dt, nbasis=5)
# })


#classic loop to avoid rsession crashes

# Initialize a list to store results
results <- list()

# Get the column names except the first one
column_names <- names(Q2)
i=0
total_iterations=length(column_names)
last_printed_percent <- 0
# Iterate over each column name


for (col_name in column_names) {
  # Process the column and store the result
  # col_name=column_names[1]
  result <- process_column(col_name, Q2, timeStampx, UpArea, dt, nbasis=8, yearstart=1951,wsize=29)
  results[[col_name]] <- result
  
  # Calculate the current percentage of completion
  i=i+1
  current_percent <- floor((i / total_iterations) * 100)
  
  # Check if we've reached a new percentage threshold
  if (current_percent > last_printed_percent) {
    cat(sprintf("Progress: %d%% complete\n", current_percent))
    last_printed_percent <- current_percent
  }
}


#save results to avoid redoing this step
save(results,file="D:/tilloal/Documents/01_Projects/RegimeShifts/Regime_fourier_30y_30d_b8_meansd.Rdata")

# Now `results` contains the parameters for each column
# merda=process_column(column_name = "290666", Q,timeStampx, UpArea, dt, yearstart=1951)
# 
# merda[[1]]
# plot(merda[[2]])
# sum(merda[[2]]$mean)
# Assume these are the saved parameters
# Replace these with the actual saved values
saved_coefs <- results[[2000]][[1]][-1]
regimo=results[[2000]][[2]]

#Now clustering of regimes

# Assuming results is already created from the lapply function



#load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Regime_fourier_30y_7d.Rdata")
# Initialize an empty list to store the extracted parameters
extracted_params <- list()

# Iterate over the results list to extract parameters
for (i in seq_along(results)) {
  # Access the parameters from each list element
  params <- results[[i]][[1]]  # params is the first element of the sublist
  
  # params=params[-length(params)]
  # Append to the extracted parameters list
  extracted_params[[i]] <- params
}

# Convert the list of parameters to a data frame
params_df <- do.call(rbind, extracted_params)
params_df <- as.data.frame(params_df)

# Assign column names to the data frame
# Assuming the first column is the column name (or ID) and the rest are coefficients
colnames(params_df) <- c("Column_ID", paste0("Coef_", seq(1, ncol(params_df) - 1)))


# Assuming params_df is your data frame with the parameters
# Remove any identifier or non-numeric columns if present
# For example, if the first column is an identifier:
params_numeric <- params_df[, -1]
# params_numero=params_numeric
# params_numero[,length(params_numeric[1,])]=params_numero[,length(params_numeric[1,])]^2
#params_numeric=as.data.frame(scale(params_numeric))
colnames(params_numeric) <- c( paste0("Coef_", seq(1, ncol(params_numeric))))
lx=length(params_numeric[1,])
params_numero=params_numeric[,-c(lx-1,lx)]
params_numeric<-params_numero
# Compute the distance matrix
dist_matrix <- dist(params_numeric, method = "euclidian")

# Perform hierarchical clustering
#hc <- hclust(dist_matrix, method = "ward.D2")  # You can choose other methods like "complete", "average", etc.
hc <-eclust(params_numeric, "hclust", hc_metric = "euclidian", graph = F, nboot=10)
# Plot the dendrogram
plot(hc, main = "Dendrogram of Parameter Sets", xlab = "Parameter Set", sub = "", cex = 0.9)
abline(h=110,col=2)


library(cluster)   # For silhouette
library(factoextra) # For visualizations
library(stats)     # For hierarchical clustering
library(ggplot2)   # For plotting




# Assuming params_numeric is your data frame with numeric parameters
set.seed(123) # For reproducibility

# K-means clustering
max_clusters <- 10
# Hierarchical clustering
# Silhouette method for hierarchical clustering
avg_silhouette_hc <- numeric(max_clusters - 1)

for (k in 2:max_clusters) {
  # Obtain cluster assignments
  hc <-eclust(params_numeric, "hclust", k=k, hc_metric = "euclidian", graph = F)
  cluster_assignments <- cutree(hc, k = k)
  # Calculate silhouette
  silhouette_values <- silhouette(cluster_assignments, dist(params_numeric))
  # Store average silhouette width
  avg_silhouette_hc[k - 1] <- mean(silhouette_values[, 3])
}

# Plot silhouette scores for hierarchical clustering
plot(2:max_clusters, avg_silhouette_hc, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Average silhouette width",
     main = "Silhouette Method for Hierarchical Clustering")

# Silhouette method for k-means clustering
silhouette_kmeans <- sapply(2:max_clusters, function(k) {
  km <- eclust(params_numeric, "kmean", k=k, hc_metric = "euclidian", graph = F)
  ss <- silhouette(km$cluster, dist(params_numeric))
  #plot(ss[,1],ss[,3])
  mean(ss[, 3])
})

# Plot silhouette scores for k-means
plot(2:max_clusters, silhouette_kmeans, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Average silhouette width",
     main = "Silhouette Method for K-Means Clustering")

# Kmeans clustering
#km.res <- eclust(params_numeric, "kmean", k = 5, hc_metric = "euclidian", graph = F)
# Optionally, cut the tree into a desired number of clusters
hc <-eclust(coef_num, "hclust",k=9, hc_metric = "euclidian", graph = F, min.cluster.size=100)
hc$gap_stat$Tab
num_clusters <- 5# Specify the number of clusters you want
clusters <- cutree(hc, k = num_clusters)

fviz_nbclust(params_numeric, hcut, method = "gap_stat",k.max=15, nboot=2) 
fviz_nbclust(params_numeric, kmeans, method = "gap_stat",k.max=15, nboot=10) 
fviz_nbclust(params_numeric, hcut, method = "silhouette",k.max=8) 
fviz_nbclust(params_numeric, hcut, method = "wss",k.max=20) 
fviz_nbclust(params_numeric, kmeans, method = "wss",k.max=20) 


num_clusters <- 5
kmeans_result <- eclust(params_numero, "kmean", k = num_clusters, hc_metric = "euclidian", graph = F)
clusters=kmeans_result$cluster
# Add cluster information to the original data frame

params_df$Cluster <- clusters

# Print the data frame with cluster assignments
print(params_df)

fanny_result <- eclust(params_numero, "fanny", k = 5, hc_metric = "euclidian", graph = F, memb.exp=1.5)
unique(fanny_result$clustering)
fanny_result$diss



## Regime archetypes ----

library(ggplot2)
library(fda)

# Assuming params_df is your data frame with the parameters and cluster assignments
# params_numeric contains the numeric parameters

params_numeric <- params_df[, -1]
# Compute archetype (mean) for each cluster
archetypes <- aggregate(params_numeric, by = list(Cluster = params_df$Cluster), FUN = mean)

# Remove the first column which is the cluster identifier used for grouping
archetypes <- archetypes[,-1]

# Define the range and basis for the fd object
# Replace with the actual range and number of basis functions used in your original fd objects
rangeval <- c(1, 365)  # Example range
nbasis <- 8  # Example number of basis functions
basis <- create.fourier.basis(rangeval, nbasis)

# Determine the number of clusters
num_clusters <- nrow(archetypes)

# Arrange plots in a single column
par(mfrow = c(1, 1))  # Set mfrow to number of clusters by 1 for a single column layout

# Create and plot fd objects for each archetype
for (i in 1:num_clusters) {
  # Create fd object for the current cluster archetype
  
  lx=length(archetypes[1,])
  coefs <- as.numeric(archetypes[i,-(c(lx-2,lx-1,lx)) ])
  cv <- as.numeric(archetypes[i,(lx-1)])/as.numeric(archetypes[i,(lx-2)])
  fd_object <- fd(coef = coefs, basisobj = basis)
  
  # Plot the fd object
  plot(fd_object, main = paste("Archetype of Cluster",i," mean CV = ",round(cv,2)), ylab = "Value", xlab = "Time")
}

#Plot where are the archetypes

alle=match(colnames(ppl),colnames(CatO_centroids))
ppoints=rbind(ppl,CatO_centroids[,alle])
matplot=match(ppoints$outlets.x,params_df$Column_ID)

ppoints$Cluster=params_df$Cluster[matplot]


#7 regimes

#regime 1: Full Nival 
#regime 2: Nivo-pluvial persistent
#regime 3: Pluvio-nival
#regime 4: Pluvial ocenaic?
#regime 5: Nivo-Pluvial 
#regime 6: Nival dry
#regime 7: Erratic


#5 regimes (with cv as variable or not)
#regime 1: Mountain Nival 
#regime 2: Nivo-pluvial
#regime 3: Pluvial
#regime 4: Transitional Nival
#regime 5: Erratic Pluvial



#6 regimes
#regime 1: Full Nival 
#regime 2: Nivo-pluvial 
#regime 3: Pluvial
#regime 4: Transitional Nival 1
#regime 5: Transitional Nival 2
#regime 6: Erratic Pluvial

color_regimes5<-c("2"="red4","4"="orchid","1"="turquoise4","5"="olivedrab3","3"="orangered")
color_regimes5bold<-c("2"="tomato4","4"="violetred4","1"="royalblue4","5"="olivedrab4","3"="orange4")


# color_regimes5<-c("1"="red4","2"="orchid","4"="turquoise4","3"="greenyellow","5"="orangered")
# color_regimes5bold<-c("1"="tomato4","2"="violetred4","4"="royalblue4","3"="olivedrab4","5"="orange4")

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,fill=factor(Cluster),size=upa),color="black",alpha=.8,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  # scale_fill_gradientn(
  #   colors=palet,oob = scales::squish, name="record length (years)",
  #   breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  scale_fill_manual(values=color_regimes5)+
  #scale_fill_brewer(palette = "Set1")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides( fill = guide_legend(override.aes = list(size = 10)),
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




# Load necessary libraries
library(fda)
library(ggplot2)


# plot all the regimes belonging to one archetype along with the mean regime and archetype

# Initialize an empty list to store the extracted parameters
extracted_regimes <- matrix(ncol=length(results),nrow=365)

# Iterate over the results list to extract parameters
for (i in seq_along(results)) {
  # Access the parameters from each list element
  regime <- results[[i]][[2]]  # params is the first element of the sublist
  
  # Append to the extracted parameters list
  extracted_regimes[,i] <- regime[,2]
}


extracted_regimes=data.frame(extracted_regimes)
colnames(extracted_regimes)=params_df$Column_ID

rm(results)
gc()
#plot for one cluster
library(fda.usc)
#do a loop with this

i=5
# Assume fd_object is your functional data object
coefs <- as.numeric(archetypes[i,-c(nbasis+1,nbasis+2) ])
cv <- as.numeric(archetypes[i,(nbasis+1)])/as.numeric(archetypes[i,(nbasis+2)])
fd_object <- fd(coef = coefs, basisobj = basis)
# Define a grid for evaluation
argvals_grid <- seq(from = min(fd_object$basis$rangeval),
                    to = max(fd_object$basis$rangeval),
                    length.out = 1000)  # Choose a suitable number of points

# Evaluate the fd object on the grid
evaluated_values <- eval.fd(argvals_grid, fd_object)

# Convert the evaluated values to a data frame
df_plot <- data.frame(
  argvals = rep(argvals_grid, ncol(evaluated_values)),
  value = as.vector(evaluated_values),
  curve = rep(1:ncol(evaluated_values), each = length(argvals_grid))
)

c1=params_df$Column_ID[which(params_df$Cluster==i)]
 
cr1=match(c1,colnames(extracted_regimes))
Regc1=extracted_regimes[,cr1]

Regc1$day <- 1:nrow(Regc1)

regime_means <- data.frame(day=Regc1$day,rm=rowMeans(Regc1[,-length(Regc1[1,])], na.rm = TRUE))
fdataobj<-fdata(t(Regc1[,-length(Regc1[1,])]),Regc1$day)
out.mode=depth.mode(fdataobj)
regime_hmodmed <- data.frame(day=Regc1$day,rm=as.numeric(out.mode$median$data))
# Assume your data frame is named `df`
# Reshape the data from wide to long format
df_long <- pivot_longer(Regc1, cols = -day, names_to = "variable", values_to = "value")

# Plot all variables on the same plot
ggplot() +
  geom_line(data=df_long, aes(x = day, y = (value), group=variable),color=color_regimes5[5],alpha=0.1) +
  geom_line(data=regime_hmodmed, aes(x= day, y=(rm)),lwd=2, col=color_regimes5bold[5])+
  #geom_line(data=regime_means, aes(x= day, y=sqrt(rm)),lwd=2, col="black")+
  geom_line(data=df_plot, aes(x = argvals, y = (value), group = curve),col="black")+
  theme_minimal() +
  labs(x = "Observation", y = "Value", title = "Plot of All Variables") +
  theme(legend.position = "none") +
scale_y_continuous(trans = "sqrt")
 

# Initialize an empty list to store plots


## beam plot of all regimes ----

class_labels <- c("Transitional Nival", "Nival", "Erratic Pluvial", "Nivo-pluvial", "Persistent Pluvial")

class_labels <- c("Nival", "Nivo-pluvial","Transitional Nival", "Pluvio-nival",  "Erratic Pluvial",  "Persistent Pluvial")

color_regimesn5<-c("Nival"="red4","Nivo-pluvial"="orchid","Transitional Nival"="turquoise4","Persistent Pluvial"="olivedrab3","Erratic Pluvial"="orangered")
color_regimes5<-c("1"="turquoise4","2"="red4","3"="orangered","4"="orchid","5"="olivedrab3")
color_regimes5bold<-c("2"="tomato4","4"="violetred4","1"="royalblue4","5"="olivedrab4","3"="orange4")

color_regimesn6<-c("Nival"="red4","Nivo-pluvial"="orchid","Transitional Nival"="slateblue1", "Pluvio-nival"="turquoise4","Persistent Pluvial"="olivedrab3","Erratic Pluvial"="orangered")
color_regimes6<-c("1"="red4","2"="orchid","3"="slateblue1","4"="turquoise4","5"="orangered","6"="olivedrab3")
color_regimes6bold<-c("1"="tomato4","2"="violetred4","3"="slateblue4","4"="royalblue4","5"="orange4","6"="olivedrab4")
# color_regimes5=color_regimes5[sort(names(color_regimes5bold))]
# color_regimesn5<-color_regimesn5[(as.numeric(names(color_regimes5bold)))]
# class_labels<-class_labels[(as.numeric(names(color_regimes5bold)))]
#color_regimes5bold=color_regimes5bold[order(as.numeric(names(color_regimes5bold)))]
# color_regimes5<-c("1"="red4","2"="orchid","3"="olivedrab3","4"="turquoise4","5"="orangered")
# color_regimes5bold<-c("1"="tomato4","2"="violetred4","3"="olivedrab4","4"="royalblue4","5"="orange4")
plot_list <- list()
# Loop from 1 to 7
for (i in 1:6) {
  # Assume fd_object is your functional data object
  ri=class_labels[i]
  lx=length(archetypes[1,])
  coefs <- as.numeric(archetypes[i,-(c(lx-2,lx-1,lx)) ])
  cv <- as.numeric(archetypes[i,(lx-1)])/as.numeric(archetypes[i,(lx-2)])
  fd_object <- fd(coef = coefs, basisobj = basis)
  
  # Define a grid for evaluation
  argvals_grid <- seq(from = min(fd_object$basis$rangeval),
                      to = max(fd_object$basis$rangeval),
                      length.out = 1000)  # Choose a suitable number of points
  
  # Evaluate the fd object on the grid
  evaluated_values <- eval.fd(argvals_grid, fd_object)
  
  # Convert the evaluated values to a data frame
  df_plot <- data.frame(
    argvals = rep(argvals_grid, ncol(evaluated_values)),
    value = as.vector(evaluated_values),
    curve = rep(1:ncol(evaluated_values), each = length(argvals_grid))
  )
  
  c1 <- params_df$Column_ID[which(params_df$Cluster == i)]
  cr1 <- match(c1, colnames(extracted_regimes))
  Regc1 <- extracted_regimes[, cr1]
  
  Regc1$day <- 1:nrow(Regc1)
  
  regime_means <- data.frame(day = Regc1$day, rm = rowMeans(Regc1[, -length(Regc1[1, ])], na.rm = TRUE))
  fdataobj <- fdata(t(Regc1[, -length(Regc1[1, ])]), Regc1$day)
  out.mode <- depth.mode(fdataobj)
  regime_hmodmed <- data.frame(day = Regc1$day, rm = as.numeric(out.mode$median$data))
  
  # Reshape the data from wide to long format
  df_long <- pivot_longer(Regc1, cols = -day, names_to = "variable", values_to = "value")
  
  ylim=c(quantile(df_long$value,.0001),quantile(df_long$value,.9999))
  # Create the plot
  p <- ggplot() +
    geom_line(data = df_long, aes(x = day, y = (value), group = variable), col = color_regimes6[i], alpha = 0.1) +
    geom_line(data = regime_hmodmed, aes(x = day, y = (rm)), lwd = 2, col = color_regimes6bold[i]) +
    geom_line(data = df_plot, aes(x = argvals, y = (value), group = curve), col = "black") +
    theme_minimal() +
    labs(x = "DOY", y = "Q/Qmean", title = paste("All Rivers for regime", ri)) +
    theme(legend.position = "none") +
    scale_y_continuous(trans = "sqrt") +
    coord_cartesian(ylim=ylim) 
  # p
  # ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/Regime_",ri,"_30y30dkmeans.jpg"),p, width=30, height=20, units=c("cm"),dpi=300) 
  # Store the plot into the list
  plot_list[[i]] <- p
}

plot_list[[1]]
# Access individual plots using plot_list[[1]], plot_list[[2]], ..., plot_list[[7]]





#ML part ----

connard=match("362439",colnames(Regc1))
 
plot(Regc1[,connard], type="o")


# random forest
library(randomForest)

rclass=params_df$Cluster[match(colnames(extracted_regimes),params_df$Column_ID)]
regimes_rf=data.frame(t(extracted_regimes))
regimes_rf$class=factor(rclass)


params_rf=params_numeric
params_rf$class <- factor(clusters)

rf_data=regimes_rf
# Assume df is your data frame and 'class' is your target variable
# Split data into training and test sets
set.seed(123)  # For reproducibility
train_indices <- sample(1:nrow(rf_data), size = 0.7 * nrow(rf_data))
train_data <- rf_data[train_indices, ]
test_data <- rf_data[-train_indices, ]

# Train a random forest model
# Assuming 'class' is the column name of the target variable
rf_model <- randomForest(class ~ ., data = train_data, ntree = 200, mtry = 5, importance = TRUE)

# Print the model summary
print(rf_model)

# Predict class membership for the test data
predictions <- predict(rf_model, newdata = test_data)

# Evaluate model performance
# Create a confusion matrix
confusion_matrix <- table(test_data$class, predictions)
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

# Train the model using caret's train function
# This uses randomForest as the underlying method
set.seed(123)  # For reproducibility
rf_model_cv <- train(class ~ ., 
                     data = train_data, 
                     method = "rf",
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

confusion_matrixT <- table(test_data$class, predictiontest)
print(confusion_matrixT)

# Calculate accuracy
accuracy <- sum(diag(confusion_matrixT)) / sum(confusion_matrixT)
print(paste("Accuracy:", accuracy))



# Set up cross-validation parameters
control <- trainControl(method = "cv", number = 10)

# Define a grid of hyperparameters to try
tune_grid <- expand.grid(mtry = c( 4,5,6))  

tune_grid=data.frame(tune_grid)
# Train the model using caret's train function
set.seed(123)  # For reproducibility
rf_model_tuned <- train(class ~ ., 
                        data = rf_data, 
                        method = "rf",
                        trControl = control,
                        tuneGrid = tune_grid,
                        ntree = 200)  # You can adjust ntree here

# Print the results of the tuning
print(rf_model_tuned)

# Get the best hyperparameters
best_params <- rf_model_tuned$bestTune
print(paste("Best parameters:", best_params))

best_accuracy <- max(rf_model_tuned$results$Accuracy)
print(paste("Best Accuracy from 10-fold CV:", best_accuracy))

#model to keep:
best_rf=rf_model_tuned$finalModel

predtune=predict(rf_model_tuned,rf_data)
# Create a confusion matrix
confusion_matrix_tuned <- table(rf_data$class, predtune)
print(confusion_matrix_tuned)

rf_data$rfclass=predtune

# Calculate accuracy
accuracy <- sum(diag(confusion_matrix_tuned)) / sum(confusion_matrix_tuned)
print(paste("Accuracy:", accuracy))



#POST ML -----
#Now I compute the regime on the last 20 years and check the differences
tstamp=timeStampx
process_regime <- function(column_name, Q,tstamp, UpArea, dt, yearstart) {
  
  wsize <- 29
  #ly <- length(years)
  yr=yearstart
  yrwindow <- yr:(yr + wsize)
  subs <- which(!is.na(match(year(as.Date(tstamp)), yrwindow)))
  
  #upstream area
  cupa=UpArea$upa[which(UpArea$outlets==as.numeric(column_name))]
  
  # Compute Qs
  data=Q[[column_name]][subs]
  data=data.frame(tstamp[subs],data)
  data$Qs <- data[,2]
  #specific discharge
  # data$Qs <- data[,2]/ (cupa * 1000000) * 1000 * 3600 * 24
  
  
  # Apply running mean
  data$daily <- tsEvaNanRunningMean(data$Qs, windowSize = 1 / dt)
  
  # Rename columns
  names(data)[1] <- "date"
  
  # Filter for specific hours
  dsel <- hour(data$date)
  datD <- data[(dsel == 12 | dsel == 13), c(1,4)]
  
  
  
  
  # Regime over the first 20 years
  #years <- unique(year(dataD$date))

  # subs <- which(!is.na(match(year(dataD$date), yrwindow)))
  # datD <- dataD[subs, ]
  #sd specific discharge
  qsd=sd(datD[,2])/ (cupa * 1000000) * 1000 * 3600 * 24
  
  #mean specific discharge
  qm <- mean(datD[,2])/ (cupa * 1000000) * 1000 * 3600 * 24
  
  # Apply running mean again
  #dataD[,2]=tsEvaNanRunningMean(dataD[,2],7)
  
  #coefficient mensuel de debit
  datD[,2]=datD[,2]/(mean(datD[,2],na.rm=T))
  Regi <- RegimeFast(datD)
  #normalisation to really cluster the shape and not the values
  #Regi$mean=Regi$mean/max(Regi$mean)
  
  # 
  # # Fit Fourier
  # value_data <- c(Regi$mean)
  # rangeval <- c(min(Regi$date), max(Regi$date))
  # basis <- create.fourier.basis(rangeval, nbasis)
  # 
  # fd_object <- Data2fd(argvals = Regi$date, y = value_data, basisobj = basis)
  # # summary(fd_object)
  # #plot(fd_object, main = paste("Functional Data Representation -", column_name))
  # 
  cv=qsd/qm
  # Save parameters
  params <- c(as.numeric(column_name), qm,qsd)
  regimeS<-Regi
  rlist=list(params,regimeS)
  return(rlist)
}



# Get the column names except the first one
column_names <- names(Q2)
i=0
total_iterations=length(column_names)
last_printed_percent <- 0
# Iterate over each column name

reregim=list()
for (col_name in column_names) {
  # Process the column and store the result
  result <- process_regime(col_name, Q2, timeStampx, UpArea, dt,1991)
  reregim[[col_name]] <- result
  
  # Calculate the current percentage of completion
  i=i+1
  current_percent <- floor((i / total_iterations) * 100)
  
  # Check if we've reached a new percentage threshold
  if (current_percent > last_printed_percent) {
    cat(sprintf("Progress: %d%% complete\n", current_percent))
    last_printed_percent <- current_percent
  }
}



# Load necessary package
library(cluster)


# Function to assign new observation to a cluster
assign_to_cluster <- function(new_obs, archetypes, clusters, num_clusters) {
  # Compute the distance from the new observation to each cluster center

  
  cluster_centers<-t(archetypes[,-length(archetypes[1,])])
  # Calculate the distance from the new observation to each cluster center
  distances <- apply(cluster_centers, 2, function(center) {
    sum((new_obs - center)^2)
  })
  # distances <- sqrt(distances)
  # Calculate the probability of the observation belonging to each cluster
  probabilities <- exp(-distances) / sum(exp(-distances))
  
  dist2=distances[-which.min(distances)]
  wm=min(dist2)
  m2=match(wm,distances)
  # Assign the observation to the nearest cluster
  assigned_cluster <- c(which.max(probabilities),min(distances),m2,wm)
  
  return(assigned_cluster)
}

# Example new observation (ensure it's preprocessed like df)
new_obs <- c(1.5, 2.3, 3.1)  # Replace with actual values





column_names <- names(Q2)
i=0
total_iterations=length(column_names)
last_printed_percent <- 0
result2<-list()
for (col_name in column_names) {
  # Process the column and store the result
  #col_name=column_names[1]
  #col_name="1029"
  result <- process_column(col_name, Q2, timeStampx, UpArea, dt, nbasis=8, yearstart=1991,wsize=29)

  new_obs<-result[[1]][-1]
  assigned_cluster <- assign_to_cluster(new_obs, archetypes, clusters, num_clusters)
  assigned_cluster
  result2[[col_name]] <- list(assigned_cluster,result[[2]])
  # Calculate the current percentage of completion
  i=i+1
  current_percent <- floor((i / total_iterations) * 100)
  
  # Check if we've reached a new percentage threshold
  if (current_percent > last_printed_percent) {
    cat(sprintf("Progress: %d%% complete\n", current_percent))
    last_printed_percent <- current_percent
  }
}

save(result2,file="D:/tilloal/Documents/01_Projects/RegimeShifts/Regime_fourier_30y_30d_b8_meansd.Rdata")
#load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Regime_fourier_30y2_1d_b5_meansd.Rdata")
# Initialize an empty list to store the extracted parameters
new_clusters <- c()
cpar=c()
# Iterate over the results list to extract parameters
for (i in seq_along(result2)) {
  # Access the parameters from each list element
  params <- result2[[i]][[1]]  # params is the first element of the sublist
  
  # params=params[-length(params)]
  # Append to the extracted parameters list
  cpar=rbind(cpar,params)
  new_clusters[i] <- params[1]
}

# Convert the list of parameters to a data frame
new_clusters <-  data.frame(ColID=as.numeric(names(result2)),newclust=new_clusters)


# Assign the new observation to a cluster
assigned_cluster <- assign_to_cluster(params_numeric[100,], params_numeric, clusters, num_clusters)
print(paste("The new observation is assigned to cluster:", assigned_cluster))

params_df$Cluster[100]




# Initialize an empty list to store the extracted parameters
new_regimes <- matrix(ncol=length(reregim),nrow=365)

# Iterate over the results list to extract parameters
for (i in seq_along(reregim)) {
  # Access the parameters from each list element
  regime <- reregim[[i]][[2]]  # params is the first element of the sublist
  
  # Append to the extracted parameters list
  new_regimes[,i] <- regime[,2]
}

#check if therre are changes

# new_regimes=data.frame(new_regimes)
# colnames(new_regimes)=params_df$Column_ID
# newregimes_rf=data.frame(t(new_regimes))
# 
# predictions <- predict(rf_model_tuned, newdata = newregimes_rf)

# Print the predictions
print(predictions)

confusion_matrix <- table(params_df$Cluster, new_clusters$newclust)
print(confusion_matrix)

# Optionally, if you want to add the predictions to the new_data dataframe:
# newregimes_rf$predicted_class <- predictions
# head(newregimes_rf)
# 
# rownames(newregimes_rf)=params_df$Column_ID
# 
# matplot=match(ppoints$outlets.x,rownames(newregimes_rf))
# 
# ppoints$Cluster2=newregimes_rf$predicted_class[matplot]
# 
# mp2=match(ppoints$outlets.x,rownames(rf_data))
# 
# ppoints$ClusterPred=rf_data$rfclass[mp2]



mp3=match(ppoints$outlets.x,new_clusters$ColID)

ppoints$ClusterR=new_clusters$newclust[mp3]

#regime 1: Full Nival 
#regime 2: Nivo-pluvial persistent
#regime 3: Pluvio-nival
#regime 4: Pluvial ocenaic?
#regime 5: Nivo-Pluvial erratic
#regime 6: Nival dry
#regime 7: Erratic


class_labels <- c("Nival", "Nivo-pluvial", "Persistent Pluvial","Transitional Nival",  "Erratic Pluvial")

class_labels <- c("Pluvio-nival", "Nival", "Erratic Pluvial", "Nivo-pluvial", "Persistent Pluvial")

class_labels <- c("Nival", "Nivo-pluvial","Transitional Nival", "Pluvio-nival",  "Erratic Pluvial",  "Persistent Pluvial")
predictions1=na.omit(ppoints$ClusterPred)
predictions2=na.omit(ppoints$Cluster2)

ppoints$C0name <- class_labels[ppoints$Cluster]
ppoints$C0pname <- class_labels[ppoints$ClusterPred]
ppoints$C1pname <- class_labels[ppoints$Cluster2]

ppoints$C1name <- class_labels[ppoints$ClusterR]

color_regimesn5<-c("Nival"="red4","Nivo-pluvial"="orchid","Pluvio-nival"="turquoise4","Persistent Pluvial"="olivedrab3","Erratic Pluvial"="orangered")

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,fill=factor(C0name),size=upa),color="black",alpha=.8,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  # scale_fill_gradientn(
  #   colors=palet,oob = scales::squish, name="record length (years)",
  #   breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
  #scale_fill_brewer(palette = "Set1")+
  scale_fill_manual(values=color_regimesn6,name="Regimes")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides( fill = guide_legend(override.aes = list(size = 10)),
          size= guide_legend(override.aes = list(fill = "grey50")))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        #legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("Regimes in 1995")

ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/Regime_map_1955_30y1dkmeans.jpg"), width=30, height=20, units=c("cm"),dpi=300) 

st_write(ppoints, "D:/tilloal/Documents/01_Projects/RegimeShifts/Regime_map_1995_30y1d_kmeans.shp")   


# Load necessary package
library(ggplot2)

# Assume predictions is a vector of predicted class labels
# If you haven't already added predictions to a data frame, you can create one:
prediction_data <- data.frame(predicted_class = ppoints$C1name)

# Plot the frequency of each class in a histogram
ggplot(prediction_data, aes(x = predicted_class)) +
  geom_bar(aes(fill=predicted_class), color = "black") +
  scale_fill_manual(values=color_regimesn6,name="Regimes")+
  theme_minimal() +
  labs(x = "Predicted Class", y = "Frequency", title = "Frequency of Predicted Classes")



class_labels <- c("Nival", "Nivo-pluvial","Transitional Nival", "Pluvio-nival",  "Erratic Pluvial",  "Persistent Pluvial")


# Install and load necessary packages
# install.packages("networkD3")
library(networkD3)

#class_labels <- c("Nival", "Nivo-pluvial", "Persistent Pluvial","Transitional Nival",  "Erratic Pluvial")
#class_labels <- c("Pluvio-nival", "Nival", "Erratic Pluvial", "Nivo-pluvial", "Persistent Pluvial")
predictions1=na.omit(ppoints$Cluster)
predictions2=na.omit(ppoints$ClusterR)

predictions1_char <- class_labels[predictions1]
predictions2_char <- class_labels[predictions2]

# Install and load necessary packages
# install.packages("ggsankey")
# install.packages("ggplot2")
library(ggsankey)
library(ggplot2)
library(dplyr)

# Assume predictions1 and predictions2 are vectors of predicted class labels from two datasets
# Create a data frame for each prediction
d <- data.frame(cbind(predictions1_char, predictions2_char))
names(d) <- c('Prediction1', 'Prediction2')

# Convert data into a long format suitable for ggsankey
df <- d %>%
  make_long(Prediction1, Prediction2)

# Optional: Define the order of nodes on the y-axis
brk <- unique(c(as.character(predictions1_char), as.character(predictions2_char)))
brk <- c("Nival", "Nivo-pluvial", "Transitional Nival","Pluvio-nival", "Persistent Pluvial",  "Erratic Pluvial")
df$node <- factor(df$node, levels = brk)
df$next_node <- factor(df$next_node, levels = brk)

# Calculate proportions for labeling
ltot <- nrow(d)
reagg <- df %>%
  dplyr::group_by(node, x) %>%
  tally()

df2 <- merge(df, reagg, by.x = 'node', by.y = 'node', all.x = FALSE)
df2 <- df2[which(df2$x.x == df2$x.y), ]
df2$np <- round(df2$n / ltot * 100)

# Create the Sankey plot
pl <- ggplot(df2, aes(x = x.x,
                      next_x = next_x,
                      node = node,
                      next_node = next_node,
                      fill = factor(node),
                      label = paste0(node, "\n", np, "%"))) +
  geom_sankey(flow.alpha = 0.5,
              node.color = "black",
              show.legend = TRUE) +
  geom_sankey_label(size = 3,
                    color = "black",
                    fill = "white") +
  theme_sankey(base_size = 18) +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = color_regimesn6) +  # Adjust colors as needed
  ggtitle("Change in Hydrological regimes \nbetween 1951-1980 and 1991-2020")

# Print the plot
print(pl)

ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/Regime_Shifts_30y_1d_kmeans.jpg"), pl, width=30, height=20, units=c("cm"),dpi=300) 
  




#difference between old and new regime ----

#mean ratios
load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Regime_fourier_30y_1d_b5_meansd.Rdata")
load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Regime_fourier_30y2_1d_b5_meansd.Rdata")


# Load necessary libraries
library(fda)
library(ggplot2)


# plot all the regimes belonging to one archetype along with the mean regime and archetype

# Initialize an empty list to store the extracted parameters
extracted_regimes <- matrix(ncol=length(results),nrow=365)

# Iterate over the results list to extract parameters
for (i in seq_along(results)) {
  # Access the parameters from each list element
  regime <- results[[i]][[2]]  # params is the first element of the sublist
  regime_rm=tsEvaNanRunningMean(regime[,2],30)
  # Append to the extracted parameters list
  extracted_regimes[,i] <- regime_rm
}

cn=names(results)
extracted_regimes=data.frame(extracted_regimes)
colnames(extracted_regimes)<-cn


rm(results)
gc()



# Initialize an empty list to store the extracted parameters
new_regimes <- matrix(ncol=length(result2),nrow=365)

# Iterate over the results list to extract parameters
for (i in seq_along(result2)) {
  # Access the parameters from each list element
  regime <- result2[[i]][[2]]  # params is the first element of the sublist
  
  # Append to the extracted parameters list
  #averaging the regime diff
  regime_rm=tsEvaNanRunningMean(regime[,2],30)
  new_regimes[,i] <- regime_rm
}


new_regimes=data.frame(new_regimes)
colnames(new_regimes)<-cn



regimes_time=list(start=extracted_regimes,end=new_regimes)
# 
save(regimes_time,file="D:/tilloal/Documents/01_Projects/RegimeShifts/Regimes_periods_v2.Rdata")

#load("D:/tilloal/Documents/01_Projects/RegimeShifts/Regimes_periods_v1.Rdata")
process_Qmean <- function(column_name, Q,tstamp, UpArea, dt, yearstart, wsize) {
  
  #ly <- length(years)
  yr=yearstart
  yrwindow <- yr:(yr + wsize)
  subs <- which(!is.na(match(year(as.Date(tstamp)), yrwindow)))
  
  #upstream area
  cupa=UpArea$upa[which(UpArea$outlets==as.numeric(column_name))]
  
  # Compute Qs
  data=Q[[column_name]][subs]
  data=data.frame(tstamp[subs],data)
  data$Qs <- data[,2]
  #specific discharge
  # data$Qs <- data[,2]/ (cupa * 1000000) * 1000 * 3600 * 24
  
  
  # Apply running mean
  data$daily <- tsEvaNanRunningMean(data$Qs, windowSize = 1 / dt)
  
  # Rename columns
  names(data)[1] <- "date"
  
  # Filter for specific hours
  dsel <- hour(data$date)
  datD <- data[(dsel == 12 | dsel == 13), c(1,4)]

  #mean specific discharge
  qm <- mean(datD[,2])/ (cupa * 1000000) * 1000 * 3600 * 24


  return(qm)
}

column_names <- names(Q2)
i=0
total_iterations=length(column_names)
last_printed_percent <- 0
resqm<- c()
for (col_name in column_names) {
  # Process the column and store the result
  # col_name=column_names[1]
  qm1 <- process_Qmean(col_name, Q2, timeStampx, UpArea, dt, yearstart=1951,wsize=29)
  qm2 <- process_Qmean(col_name, Q2, timeStampx, UpArea, dt, yearstart=1991,wsize=29)
  qrat <- qm2/qm1
  
  rq<- as.numeric(c(col_name,qm1,qm2,qrat))
  resqm=rbind(resqm,rq)
  # Calculate the current percentage of completion
  i=i+1
  current_percent <- floor((i / total_iterations) * 100)
  
  # Check if we've reached a new percentage threshold
  if (current_percent > last_printed_percent) {
    cat(sprintf("Progress: %d%% complete\n", current_percent))
    last_printed_percent <- current_percent
  }
}

load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/meanQratio_periods.Rdata")
mean(resqm[,2])
hist(resqm[,3])
library(dplyr)

new_regimes <- regimes_time[[2]]
extracted_regimes<- regimes_time[[1]]

SI=c()
SIx<-c()
for (i in 1:length(resqm[,1])){
  test<-sum(abs(new_regimes[,i]*resqm[i,3]-resqm[i,3]))/(resqm[i,3]*365)
  SI<-c(SI,test)
  #pinaise=length(which((new_regimes[,i]*resqm[i,3]-resqm[i,3])<0))/365
  pinaise=length(which((new_regimes[,i]<quantile(extracted_regimes[,i],.1))))/365
  SIx<-c(SIx,pinaise)
}

plot(SIx)
  
SIi=c()
SIix<-c()
for (i in 1:length(resqm[,1])){
  test<-sum(abs(extracted_regimes[,i]*resqm[i,2]-resqm[i,2]))/(resqm[i,2]*365)
  SIi<-c(SIi,test)
  pinaise=length(which((extracted_regimes[,i]<quantile(extracted_regimes[,i],.1))))/365
  SIix<-c(SIix,pinaise)
}

mean(SIi)

msi=match(ppoints$outlets.x,as.numeric(colnames(extracted_regimes)))

ppoints$SI1=SIi[msi]
ppoints$SI2=SI[msi]
ppoints$SIx2=SIx[msi]
ppoints$SIx1=SIix[msi]
palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=SI1),color="black",alpha=1,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="Seasonality index",
    breaks=seq(-1,1,.2), limits=c(0,.5)) +
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

ppoints$SIch=ppoints$SI2-ppoints$SI1
ppoints$SIxch=ppoints$SIx2-ppoints$SIx1
min(ppoints$SIch)
palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=SIch),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="Change",
    breaks=seq(-.5,.5,.05), limits=c(-.2,.2)) +
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
  ggtitle("Change in flow seasonality index")

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/FSI_change.png", width=20, height=15, units=c("cm"),dpi=1500)


#Mean flow
resqm=data.frame(resqm)
mfm=match(ppoints$outlets.x,resqm$X1)
ppoints$qmean_ratio=resqm$X4[mfm]
palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=qmean_ratio),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="Change",
    breaks=seq(-.5,2,.05), limits=c(0.8,1.2)) +
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
  ggtitle("Change in Qmean")

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/qmean_change.png", width=20, height=15, units=c("cm"),dpi=500)




#Clustering-----

#new_regimes_m <- new_regimes/resqm[,4]

#new_regimes_m <- new_regimes
new_regimes_m <- sweep(new_regimes, 2, resqm[,4], "*")

new_regimes_m[,1]/new_regimes[,1]

regime_diff=new_regimes_m-extracted_regimes



i=0
total_iterations=length(cn)
last_printed_percent <- 0
coef_df<- c()
library(fda)
for (i in 1:length(resqm[,1])) {
  # Process the column and store the result
  # col_name=column_names[1]
  value_data <- c(regime_diff[,i])
  #plot(value_data)
  rangeval <- c(1,365)
  nbasis=6
  basis <- create.fourier.basis(rangeval, nbasis)
  #basis <- create.bspline.basis(rangeval, nbasis)
  
  fd_object <- Data2fd(argvals = c(1:365), y = value_data, basisobj = basis)
  # summary(fd_object)
  # plot(fd_object, main = paste("Functional Data Representation"))
  # plot(regime_diff[,i])
  # plot(extracted_regimes[,i],type="l",lwd=3)
  # lines(new_regimes_m[,i],lwd=3, col=2)
  
  coefs<-c(resqm[i,1],fd_object$coefs)
  coef_df<-rbind(coef_df,coefs)
  # Calculate the current percentage of completion
  current_percent <- floor(((i / total_iterations) * 100)/10)
  
  # Check if we've reached a new percentage threshold
  if (current_percent > last_printed_percent) {
    cat(sprintf("Progress: %d%% complete\n", current_percent*10))
    last_printed_percent <- current_percent
  }
}

colnames(coef_df) <- c("Column_ID", paste0("Coef_", seq(1, ncol(coef_df) - 1)))
coef_df=data.frame(coef_df)

# Fit Fourier
i=1900
value_data <- c(regime_diff[,i])
rangeval <- c(1,365)
nbasis=5
basis <- create.fourier.basis(rangeval, nbasis)
#basis <- create.bspline.basis(rangeval, nbasis)

fd_object <- Data2fd(argvals = c(1:365), y = value_data, basisobj = basis)
summary(fd_object)
plot(fd_object, main = paste("Functional Data Representation"))
regimean=tsEvaNanRunningMean(regime_diff[,i],30)
lines(regimean, type="o")
abline(h=0)
plot(extracted_regimes[,i],type="l",lwd=3)
lines(new_regimes_m[,i],lwd=3, col=2)


###cluster fourier----
library(cluster)   # For silhouette
library(factoextra) # For visualizations
library(stats)     # For hierarchical clustering
library(ggplot2)   # For plotting

rownames(coef_df)<-cn
coef_num=data.frame(coef_df[,-1])
# Perform hierarchical clustering
#hc <- hclust(dist_matrix, method = "ward.D2")  # You can choose other methods like "complete", "average", etc.
k=5
hc <-eclust(coef_num, "hclust", k=k, hc_metric = "euclidian", graph = F)
# Plot the dendrogram
plot(hc, main = "Dendrogram of Parameter Sets", xlab = "Parameter Set", sub = "", cex = 0.9)
abline(h=110,col=2)





# Assuming params_numeric is your data frame with numeric parameters
set.seed(123) # For reproducibility

# K-means clustering
max_clusters <- 10
# Hierarchical clustering
# Silhouette method for hierarchical clustering
avg_silhouette_hc <- numeric(max_clusters - 1)

for (k in 2:max_clusters) {
  # Obtain cluster assignments
  hc <-eclust(coef_num, "hclust", k=k, hc_metric = "euclidian", graph = F)
  cluster_assignments <- cutree(hc, k = k)
  # Calculate silhouette
  silhouette_values <- silhouette(cluster_assignments, dist(coef_num))
  # Store average silhouette width
  avg_silhouette_hc[k - 1] <- mean(silhouette_values[, 3])
}

# Plot silhouette scores for hierarchical clustering
plot(2:max_clusters, avg_silhouette_hc, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Average silhouette width",
     main = "Silhouette Method for Hierarchical Clustering")

# Silhouette method for k-means clustering
silhouette_kmeans <- sapply(2:max_clusters, function(k) {
  km <- eclust(coef_num, "kmean", k=k, hc_metric = "euclidian", graph = F)
  ss <- silhouette(km$cluster, dist(coef_num))
  #plot(ss[,1],ss[,3])
  mean(ss[, 3])
})

# Plot silhouette scores for k-means
plot(2:max_clusters, silhouette_kmeans, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K",
     ylab = "Average silhouette width",
     main = "Silhouette Method for K-Means Clustering")

# Kmeans clustering
#km.res <- eclust(params_numeric, "kmean", k = 5, hc_metric = "euclidian", graph = F)
# Optionally, cut the tree into a desired number of clusters
hc <-eclust(coef_num, "hclust",k=9, hc_metric = "euclidian", graph = F, min.cluster.size=100)
hc$gap_stat$Tab
num_clusters <- 5# Specify the number of clusters you want
clusters <- cutree(hc, k = num_clusters)

fviz_nbclust(coef_num, hcut, method = "gap_stat",k.max=15, nboot=2) 
fviz_nbclust(coef_num, kmeans, method = "gap_stat",k.max=15, nboot=10) 
fviz_nbclust(coef_num, hcut, method = "silhouette",k.max=8) 
fviz_nbclust(coef_num, hcut, method = "wss",k.max=20) 


num_clusters <- 5


# fviz_nbclust(coef_num, hcut, method = "gap_stat",k.max=8, nboot=2) 
# fviz_nbclust(coef_num, kmeans, method = "gap_stat",k.max=8, nboot=10) 
# fviz_nbclust(coef_num, kmeans, method = "silhouette",k.max=8) 
# fviz_nbclust(coef_num, kmeans, method = "wss",k.max=20) 
# 
kmeans_result <- eclust(coef_num, "kmean", k=num_clusters, nboot=20,hc_metric = "euclidian", graph = F)
kmeans_result
clusters=kmeans_result$cluster
# # Add cluster information to the original data frame

coef_df$Cluster <- clusters


library(ggplot2)
library(fda)

# Assuming params_df is your data frame with the parameters and cluster assignments
# params_numeric contains the numeric parameters

# Compute archetype (mean) for each cluster
archetypes <- aggregate(coef_num, by = list(Cluster = coef_df$Cluster), FUN = mean)

# Remove the first column which is the cluster identifier used for grouping
archetypes <- archetypes[,-1]

# Define the range and basis for the fd object
# Replace with the actual range and number of basis functions used in your original fd objects
rangeval <- c(1, 365)  # Example range
nbasis <- 5  # Example number of basis functions
basis <- create.fourier.basis(rangeval, nbasis)

# Determine the number of clusters
num_clusters <- nrow(archetypes)

# Arrange plots in a single column
par(mfrow = c(1, 1))  # Set mfrow to number of clusters by 1 for a single column layout

resqm=data.frame(resqm)
resqm$cluster=clusters
# Create and plot fd objects for each archetype
for (i in 1:num_clusters) {
  # Create fd object for the current cluster archetype
  
  coefs <- as.numeric(archetypes[i, ])
  # cv <- as.numeric(archetypes[i,(nbasis+2)])
  fd_object <- fd(coef = coefs, basisobj = basis)
  selr=which(resqm$cluster==i)
  mqc=(mean(resqm[selr,4])-1)*100
  # Plot the fd object
  plot(fd_object, main = paste("Archetype of Cluster",i, " | mean Q change = ", round(mqc)," %"), ylab = "Value", xlab = "Time")
}

#Plot where are the archetypes


# alle=match(colnames(ppl),colnames(CatO_centroids))
# ppoints=rbind(ppl,CatO_centroids[,alle])
matplot=match(ppoints$outlets.x,coef_df$Column_ID)

ppoints$Cluster=coef_df$Cluster[matplot]


#7 regimes

#regime 1: Full Nival 
#regime 2: Nivo-pluvial persistent
#regime 3: Pluvio-nival
#regime 4: Pluvial ocenaic?
#regime 5: Nivo-Pluvial 
#regime 6: Nival dry
#regime 7: Erratic


#5 regimes (with cv as variable or not)
#regime 1: Mountain Nival 
#regime 2: Nivo-pluvial
#regime 3: Pluvial
#regime 4: Transitional Nival
#regime 5: Erratic Pluvial



#6 regimes
#regime 1: Full Nival 
#regime 2: Nivo-pluvial 
#regime 3: Pluvial
#regime 4: Transitional Nival 1
#regime 5: Transitional Nival 2
#regime 6: Erratic Pluvial

color_regimes5<-c("2"="red4","4"="orchid","1"="turquoise4","5"="olivedrab3","3"="orangered")
color_regimes5bold<-c("2"="tomato4","4"="violetred4","1"="royalblue4","5"="olivedrab4","3"="orange4")


# color_regimes5<-c("1"="red4","2"="orchid","4"="turquoise4","3"="greenyellow","5"="orangered")
# color_regimes5bold<-c("1"="tomato4","2"="violetred4","4"="royalblue4","3"="olivedrab4","5"="orange4")
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
color_vector <- sapply(class_labels, function(x) color_class_labels[x])
print(color_vector)


ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,fill=factor(Cluster),size=upa),color="black",alpha=.8,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  # scale_fill_gradientn(
  #   colors=palet,oob = scales::squish, name="record length (years)",
  #   breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  scale_fill_manual(values=color_class_plot)+
  #scale_fill_brewer(palette = "Set1")+
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides( fill = guide_legend(override.aes = list(size = 10)),
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


#Archetype 1: winter wetting and snowmelt shift
#Archetype 2: wetter winter and dryer summer
#Archetype 3: earlier snowmelt
#Archetype 4: general drying peaking in spring
#Archetype 5: severe general drying peaking in winter


#Archetype 1: earlier snowmelt
#Archetype 2: strong winter wetting 
#Archetype 3: severe general drying peaking in winter
#Archetype 4: general drying peaking in summer
#Archetype 5: spring and summer drying
#Archetype 6: wetter winter and dryer summer
#Archetype 7: moderate general drying peaking in winter

## Change archetypes beam plot ----

class_labels <- c( "Earlier snomelt", "Strong wetting - winter peak", "Strong drying - winter peak","Drying - summer peak","Drying - spring and summer" ,
                   "Wetter winter - dryer summer","Drying - winter peak" )

color_vector <- sapply(class_labels, function(x) color_class_labels[x])
print(color_vector)

color_vectorbold<-c("turquoise4","darkgreen","brown4","orange4","coral4","darkorchid4","darkgoldenrod4")
label_to_number <- c(
  "Earlier snomelt" = 1,
  "Strong wetting - winter peak" = 2,
  "Strong drying - winter peak" = 3,
  "Drying - summer peak" = 4,
  "Drying - spring and summer" = 5,
  "Wetter winter - dryer summer" = 6,
  "Drying - winter peak" = 7
)


color_regimes5<-c("1"="turquoise4","2"="red4","3"="orangered","4"="orchid","5"="olivedrab3")
color_regimes5bold<-c("2"="tomato4","4"="violetred4","1"="royalblue4","5"="olivedrab4","3"="orange4")
# color_regimes5=color_regimes5[sort(names(color_regimes5bold))]
# color_regimesn5<-color_regimesn5[(as.numeric(names(color_regimes5bold)))]
# class_labels<-class_labels[(as.numeric(names(color_regimes5bold)))]
color_regimes5bold=color_regimes5bold[order(as.numeric(names(color_regimes5bold)))]
# color_regimes5<-c("1"="red4","2"="orchid","3"="olivedrab3","4"="turquoise4","5"="orangered")
# color_regimes5bold<-c("1"="tomato4","2"="violetred4","3"="olivedrab4","4"="royalblue4","5"="orange4")
plot_list <- list()
library(fda.usc)
# Loop from 1 to 7
for (i in 1:5) {
  i=1
  # Assume fd_object is your functional data object
  ri=class_labels[i]
  coefs <- as.numeric(archetypes[i,])
  selr=which(resqm$cluster==i)
  mqc=(mean(resqm[selr,4])-1)*100
  fd_object <- fd(coef = coefs, basisobj = basis)
  
  # Define a grid for evaluation
  argvals_grid <- seq(from = min(fd_object$basis$rangeval),
                      to = max(fd_object$basis$rangeval),
                      length.out = 1000)  # Choose a suitable number of points
  
  # Evaluate the fd object on the grid
  evaluated_values <- eval.fd(argvals_grid, fd_object)
  
  # Convert the evaluated values to a data frame
  df_plot <- data.frame(
    argvals = rep(argvals_grid, ncol(evaluated_values)),
    value = as.vector(evaluated_values),
    curve = rep(1:ncol(evaluated_values), each = length(argvals_grid))
  )
  
  c1 <- coef_df$Column_ID[which(coef_df$Cluster == i)]
  cr1 <- match(c1, colnames(regime_diff))
  Regc1 <- regime_diff[, cr1]
  
  Regc1$day <- 1:nrow(Regc1)
  
  regime_means <- data.frame(day = Regc1$day, rm = rowMeans(Regc1[, -length(Regc1[1, ])], na.rm = TRUE))
  fdataobj <- fdata(t(Regc1), Regc1$day)
  out.mode <- depth.mode(fdataobj)
  regime_hmodmed <- data.frame(day = Regc1$day, rm = as.numeric(out.mode$median$data))
  
  # Reshape the data from wide to long format
  df_long <- pivot_longer(Regc1, cols = -day, names_to = "variable", values_to = "value")
  
  ylim=c(quantile(df_long$value,.005),quantile(df_long$value,.995))
  # Create the plot
  p <- ggplot() +
    geom_line(data = df_long, aes(x = day, y = (value), group = variable), col = color_vector[i], alpha = 0.1) +
    geom_line(data = regime_hmodmed, aes(x = day, y = (rm)), lwd = 2, col = color_vectorbold[i]) +
    geom_line(data = df_plot, aes(x = argvals, y = (value), group = curve), col = "black") +
    theme_minimal() +
    labs(x = "DOY", y = "R1-R2", title = paste("All Rivers for shift", ri)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim=ylim) 
  p
  
  ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/Shifts_",ri,"_30y1d.jpg"),p, width=30, height=20, units=c("cm"),dpi=300) 
  # Store the plot into the list
  plot_list[[i]] <- p
}

plot_list[[6]]

ppoints$Cshift <- class_labels[ppoints$Cluster]

#st_write(ppoints, "D:/tilloal/Documents/01_Projects/RegimeShifts/Regime_shifts_hclust.shp")  




#plot change in SI and change in qmean by cluster -----

msi=match(ppoints$outlets.x,as.numeric(colnames(extracted_regimes)))


ppoints$SI1=SIi[msi]
ppoints$SI2=SI[msi]
ppoints$SIch=ppoints$SI2-ppoints$SI1

mfm=match(ppoints$outlets.x,resqm$X1)
ppoints$qmean_ratio=resqm$X4[mfm]


ggplot()+
  geom_point(data=ppoints, aes(x=SIch,y=qmean_ratio,fill=factor(Cluster),size=upa),color="black",alpha=.5,shape=21,stroke=0.1)+
  scale_fill_manual(values=color_class_plot)+
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  guides( fill = guide_legend(override.aes = list(size = 10)),
          size= guide_legend(override.aes = list(fill = "grey50")))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "white", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
  

######################### OLD STUFF ############################################












#Sliding window regime

#10 years subsample of data
years=unique(year(Ticido$timeStampx))
wsize=19
ly=length(years)
Regsave=c()
for(yr in years[-c((ly-wsize+1):ly)]){
  print(yr)
  yrwindow=c(yr:(yr+wsize))
  subs=which(!is.na(match(year(Ticido$timeStampx),yrwindow)))
  Ticid=Ticino[subs,]
  Regi=RegimeFast(Ticid)
  Regsave=rbind(Regsave,Regi$mean)
  
}

plot(Regsave[1,])


# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

dcad=seq(1951,2001,by=10)
df=data.frame(t(Regsave))
colnames(df)=years[-c((ly-wsize+1):ly)]
keep=match(dcad,years[-c((ly-wsize+1):ly)])
df=df[,keep]
# Reshape the data from wide to long format
df_long <- df %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = -row_id, names_to = "variable", values_to = "value")

# Plot using ggplot2
ggplot(df_long, aes(x = row_id, y = value, color = variable)) +
  geom_line(lwd=3) +
  scale_color_viridis_d(option = "viridis", begin = 0, end = 1) +
  theme_minimal() +
  labs(title = "Line Plot of 60 Variables", x = "Row ID", y = "Value", color = "Variable")



# Extract the column you want to convert
value_data <- c(df[,1],df[1,1])

# Define the range of the data
rangeval <- c(1,366)

# Define the basis system (e.g., B-spline basis)
# Here, we choose 10 basis functions
nbasis <- 5
basis <- create.fourier.basis(rangeval, nbasis)

fdata_obj <- fdata(value_data, argvals = c(1:365), rangeval = rangeval)

# Convert the data to a functional data object
fd_object <- Data2fd(argvals = c(1:(366)), y = value_data, basisobj = basis)

# Print summary
summary(fd_object)

# Plot the functional data object
plot(fd_object, main = "Functional Data Representation")


Ticino=Q[["290666"]]
plot(Ticino)
Ticino=data.frame(timeStampx,Ticino)
cupa=UpArea$upa[which(UpArea$outlets==290666)]
Ticino$Qs=Ticino$Ticino/(cupa*1000000)*1000*3600*24
Ticido=Ticino[,c(1,3)]

Ticido$daily=tsEvaNanRunningMean(Ticido[,2], windowSize =1/dt)
# plot(timeAndSeries$data[1:700])
# lines(timeAndSeries$rollmean[1:700],lwd=2,col="red")
data=Ticido
names(data)[c(2,3)]=c("Qs","Qs")
names(data)[1]="date"
dsel=hour(data$date)
dataD=data[which(dsel==12 | dsel==13),c(1,3)]
plot(dataD)

dataD[,2]=tsEvaNanRunningMean(dataD[,2],30)


RegimeQ=RegimeFast(Ticido)
#regime over the first 20 years
years=unique(year(dataD$date))
wsize=19
ly=length(years)
Regsave=c()
yr=years[1]
yrwindow=c(yr:(yr+wsize))
subs=which(!is.na(match(year(dataD$date),yrwindow)))
datD=dataD[subs,]
Regi=RegimeFast(datD)
plot(Regi)
Regsave=rbind(Regsave,Regi$mean)

#fit Fourier
# Extract the column you want to convert
value_data <- c(Regi$mean)
# Define the range of the data
rangeval <- c(min(Regi$date),max(Regi$date))
# Define the basis system (e.g., B-spline basis)
# Here, we choose 10 basis functions
nbasis <- 5
basis <- create.fourier.basis(rangeval, nbasis)

#fdata_obj <- fdata(value_data, argvals = Regi$date, rangeval = rangeval)

# Convert the data to a functional data object
fd_object <- Data2fd(argvals =Regi$date, y = value_data, basisobj = basis)

# Print summary
summary(fd_object)

# Plot the functional data object
plot(fd_object, main = "Functional Data Representation")

#save parameters
params=c(idsta,fd_object$coefs)












tsEstimateAverageSeasonality= function (timeStamps, seasonalitySeries, timeWindow) 
{
  avgYearLength <- 365.2425
  nMonthInYear <- 12
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") {
    dt <- dt/24
    firsthour = as.integer(format(timeStamps[1], "%H"))
    while (firsthour >= 19) {
      timeStamps = timeStamps + 3600
      firsthour = as.integer(format(timeStamps[1], "%H"))
    }
  }
  avgYearLength <- avgYearLength/dt
  timstamp <- unique(as.integer(format(timeStamps, "%Y")))
  nYear <- round(length(timeStamps)/avgYearLength)
  tw <- round(timeWindow/avgYearLength)
  avgMonthLength <- avgYearLength/nMonthInYear
  firstTmStmp <- timeStamps[1]
  lastTmStmp <- timeStamps[length(timeStamps)]
  mond <- as.numeric(format(timeStamps, "%m"))
  mony <- format(timeStamps, "%Y-%m")
  monthTmStampStart <- c(timeStamps[1], timeStamps[which(diff(mond) != 
                                                           0) + 1])
  monthTmStampEnd <- c(timeStamps[which(diff(mond) != 0)], 
                       timeStamps[length(timeStamps)])
  seasonalitySeries <- data.frame(time = timeStamps, series = seasonalitySeries, 
                                  month = mony, mond = mond)
  grpdSsn_ <- stats::aggregate(seasonalitySeries$series, by = list(month = seasonalitySeries$month), 
                               FUN = function(x) mean(x, na.rm = T))
  nYears <- ceiling(length(grpdSsn_$month)/nMonthInYear)
  grpdSsn <- matrix(NaN, nYears * nMonthInYear, 1)
  grpdSsn[1:length(grpdSsn_$x)] <- grpdSsn_$x
  grpdSsnMtx <- matrix(grpdSsn, nrow = 12, byrow = F)
  mnSsn_ <- rowMeans(grpdSsnMtx)
  mnSsn_t <- c()
  for (y in (1 + tw/2):(nYears - tw/2)) {
    grpdSsnMti <- grpdSsnMtx[, c((y - tw/2):(y + tw/2))]
    twz <- round(timeWindow/avgYearLength)
    twv <- max(2, twz)
    mnSsn_i <- rowMeans(grpdSsnMti[, c(1:twv)])
    it <- (1:nMonthInYear)
    x <- it/6 * pi
    dx <- pi/6
    a0 <- mean(mnSsn_)
    a1 <- 1/pi * sum(cos(x) * mnSsn_i) * dx
    b1 <- 1/pi * sum(sin(x) * mnSsn_i) * dx
    a2 <- 1/pi * sum(cos(2 * x) * mnSsn_i) * dx
    b2 <- 1/pi * sum(sin(2 * x) * mnSsn_i) * dx
    a3 <- 1/pi * sum(cos(3 * x) * mnSsn_i) * dx
    b3 <- 1/pi * sum(sin(3 * x) * mnSsn_i) * dx
    mnSsni <- a0 + (a1 * cos(x) + b1 * sin(x)) + (a2 * cos(2 * x) + b2 * sin(2 * x)) + (a3 * cos(3 * x) + b3 * sin(3 * x))
    mnSsn_t <- c(mnSsn_t, mnSsni)
  }
  it <- (1:nMonthInYear)
  x <- it/6 * pi
  dx <- pi/6
  a0 <- mean(mnSsn_)
  a1 <- 1/pi * sum(cos(x) * mnSsn_) * dx
  b1 <- 1/pi * sum(sin(x) * mnSsn_) * dx
  a2 <- 1/pi * sum(cos(2 * x) * mnSsn_) * dx
  b2 <- 1/pi * sum(sin(2 * x) * mnSsn_) * dx
  a3 <- 1/pi * sum(cos(3 * x) * mnSsn_) * dx
  b3 <- 1/pi * sum(sin(3 * x) * mnSsn_) * dx
  mnSsn <- a0 + (a1 * cos(x) + b1 * sin(x)) + (a2 * cos(2 * x) + b2 * sin(2 * x)) 
  #+ (a3 * cos(3 * x) + b3 * sin(3 * x))
  pt <- matrix(1, nYears, 1)
  monthAvgVec <- rep(mnSsn, nYears)
  imnth <- c(0:(length(monthAvgVec) - 1))
  avgTmStamp <- as.Date(firstTmStmp) + (avgMonthLength * dt)/2 + 
    imnth * (avgMonthLength * dt)
  imReg <- c(0:(length(mnSsn) - 1))
  regimeTmStamp <- 1 + (avgMonthLength * dt)/2 + imReg * (avgMonthLength * 
                                                            dt)
  monthAvgVec <- c(monthAvgVec[1], monthAvgVec, monthAvgVec[length(monthAvgVec)])
  avgTmStamp <- c(firstTmStmp, avgTmStamp, lastTmStmp)
  regimeVec <- c(mnSsn[1], mnSsn, mnSsn[length(mnSsn)])
  regimeTmStamp <- c(1, regimeTmStamp, 365)
  sidefill = (length(avgTmStamp) - length(mnSsn_t))/2
  monthAvgVex <- c(mnSsn_t[1:(sidefill)], mnSsn_t, mnSsn_t[(length(mnSsn_t) - 
                                                              (sidefill - 1)):length(mnSsn_t)])
  avgTmStamp <- as.numeric(avgTmStamp)
  timeStampsN <- as.numeric(timeStamps)
  regime <- pracma::interp1(regimeTmStamp, regimeVec, c(1:365), 
                            method = "spline")
  averageSeasonalitySeries <- pracma::interp1(avgTmStamp, monthAvgVec, 
                                              timeStampsN, method = "spline")
  varyingSeasonalitySeries <- pracma::interp1(avgTmStamp, monthAvgVex, 
                                              timeStampsN, method = "spline")
  return(list(regime = regime, Seasonality = data.frame(averageSeasonalitySeries = averageSeasonalitySeries, 
                                                        varyingSeasonalitySeries = varyingSeasonalitySeries)))
}
rm(tsEstimateAverageSeasonality)
library(RtsEva)

haaa=tsEstimateAverageSeasonality(Ticido$timeStampx,Ticido$Qs,timeWindow = 365*30*4)
plot(haaa$regime)

Ticisno=data.frame(timeStampx,Ticisno)
Ziz=plotRegime(Ticisno,"Ticino")

cname=colnames(tavg)
cname[1]="time"
cname=gsub("X","",cname)
colnames(tavg)=cname
tavg$tday=as.Date(tavg$time)