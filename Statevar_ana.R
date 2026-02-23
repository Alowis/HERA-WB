source("functions_regime.R")

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

#keep only rivers in EU domain
out1=outletopen(hydroDir,"efas_rnet_100km_01min")
out1$latlong=paste(round(out1$Var1,4),round(out1$Var2,4),sep=" ")
outhybas=inner_join(out1,outhybas, by="latlong")

#Upstream catchments
UpCat<-st_read(dsn="D:/tilloal/Documents/01_Projects/RegimeShifts/Upstreamgroups.shp")
upup<-UpCat[which(UpCat$up=="Upstream Catchments"),]
hybas2uc<-match(upup$outlets_x,outhybas$outlets.y)
outhybas<-outhybas[hybas2uc,]

ppl <- st_as_sf(UpArea, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

#remove Q obs that are outside
matcat=match(outhybas$latlong,UpArea$latlong)
UpArea=UpArea[matcat,]

#Plot of important variables ----


## Variable plotting Section ----
#to be saved for every catchment at yearly timescale:
#direct runoff fraction
#snow fraction
#concentration of precipitation
#total precipitation
#aridity index
#AEP
#w parameter budyko
#ET0

###Rain and snow -----
Rain<- fread(paste0(hydroDir,"/tss/HERA_SocCF/RainUpsX_1951_2020.csv"),header=TRUE)
time=Rain$V1
timeStampX=time[order(time)]
Rain=Rain[order(time),]

Snow<- fread(paste0(hydroDir,"/tss/HERA_Histo/SnowUpsX_1951_2020.csv"),header=TRUE)
Snow=Snow[order(time),]

matcol=match(UpArea$outlets,as.numeric(colnames(Rain)))
cnames=names(Rain)[matcol]

Snow <- Snow[, .SD, .SDcols = matcol]
Rain <- Rain[, .SD, .SDcols = matcol]

Precipitation <- Rain + Snow
Snowfraction<-Snow/Precipitation

dsel <- hour(timeStampX)
dt=4
SaveAgg=c(1951:2020)
for (col_name in cnames) {
  # col_name="290666"
  print(col_name)
  Trun <- preprocess_frac(time=timeStampX, input_var=Snowfraction, col_name, dsel,dt)
  Yagg <- process_frac(Trun)
  #plot(Yagg)
  SaveAgg=cbind(SaveAgg,Yagg$val)
}

SaveAgg=data.frame(SaveAgg)
colnames(SaveAgg)[-1]=cnames

save(SaveAgg,file="D:/tilloal/Documents/01_Projects/RegimeShifts/SnowFraction.Rdata")

#loading snowfraction
load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/SnowFraction.Rdata")
Keep=t(SaveAgg[which(SaveAgg$SaveAgg==1955),])[-1]

ys1=c(1951:1970)
my=which(!is.na(match(SaveAgg$SaveAgg,ys1)))
Keep=(SaveAgg[my,-1])
Keep_means=colMeans(Keep)

#Plot of snow fraction per catchment

CatUpA=inner_join(Catf7,UpArea,by=c("llcoord"="latlong"))
matcat=match(outhybas$latlong,CatUpA$llcoord)
CatUpA=CatUpA[matcat,]
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

ppl <- st_as_sf(CatLR, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

ppc <- st_as_sf(CatO, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppc <- st_transform(ppc, crs = 3035)

# Compute centroids
CatO_centroids <- st_centroid(ppc)
alle=match(colnames(ppl),colnames(CatO_centroids))
ppoints=rbind(ppl,CatO_centroids[,alle])

km=match(ppoints$outlets.y,cnames)
ppoints$snowfraction=Keep_means[km]

### Snow fraction ----
palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=snowfraction),color="black",alpha=1,shape=21,stroke=0)+ 
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

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/snowfraction.png", width=20, height=15, units=c("cm"),dpi=1500)

rm(Rain,Snow)
gc()

### Aridity, Evaportation, w budyko coefficient ----
ActEvapo<- fread(paste0(hydroDir,"/tss/HERA_Histo/ActEvapo_1951_2020.csv"),header=TRUE)
ActEvapo=ActEvapo[order(time),]

PEvapo<- fread(paste0(hydroDir,"/tss/HERA_Histo/etUpsX_1951_2020.csv"),header=TRUE)
PEvapo=PEvapo[order(time),]

ActEvapo <- ActEvapo[, .SD, .SDcols = matcol]
PEvapo <- PEvapo[, .SD, .SDcols = matcol]

SaveE0=c(1951:2020)
SaveAET=c(1951:2020)
SavePrecip=c()
SeasonPrecip=c()
for (col_name in cnames) {
  #col_name=cnames[100]
  col_name="290666"
  print(col_name)
  #TE0 <- preprocess_in(time=timeStampX, input_var=PEvapo, col_name, dsel,dt)
  AET <- preprocess_in(time=timeStampX, input_var=ActEvapo, col_name, dsel,dt)
  Precip <- preprocess_in(time=timeStampX, input_var=Precipitation, col_name, dsel,dt)
  E0Acc <- process_data(TE0)
  AETAcc <- process_data(AET)
  PrecipAcc <- process_precip(Precip)
  PrecipAcc$seaonal$cat=rep(col_name,2)
  PrecipAcc$yagg$cat=rep(col_name,length(PrecipAcc$yagg$year))
  #plot(Yagg)
  SavePrecip=rbind(SavePrecip,PrecipAcc$yagg)
  SeasonPrecip=rbind(SeasonPrecip,PrecipAcc$seaonal)
  SaveAET=cbind(SaveAET,AETAcc$val)
  SaveE0=cbind(SaveE0,E0Acc$val)
}

save(SavePrecip,file="D:/tilloal/Documents/01_Projects/RegimeShifts/PrecipitationMetrics_up.Rdata")
save(SeasonPrecip,file="D:/tilloal/Documents/01_Projects/RegimeShifts/PrecipitationSeason_up.Rdata")
save(SaveAET,file="D:/tilloal/Documents/01_Projects/RegimeShifts/AETMetrics_up.Rdata")
#save(SaveE0,file="D:/tilloal/Documents/01_Projects/RegimeShifts/E0Metrics_up.Rdata")



SaveAET=data.frame(SaveAET)

#Select a specific range of years

ys1=c(1951:1970)
my=which(!is.na(match(SaveAET$SaveAET,ys1)))
KeepAET=(SaveAET[my,-1])
KAET_means=colMeans(KeepAET)


km=match(ppoints$outlets.y,cnames)
ppoints$AET1=KAET_means[km]

palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=AET1),color="black",alpha=1,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="AET (mm/y)") +
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

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/snowfraction.png", width=20, height=15, units=c("cm"),dpi=1500)


### Precipitation ----

myp=which(!is.na(match(SavePrecip $year,ys1)))

KeepP=SavePrecip[myp,c(1,2,4)]

KeepPag<- stats::aggregate(list(pmean=KeepP$val.sum), by = list(catch = KeepP$cat), 
                             FUN = function(x) c(pmean=mean(x, na.rm = T)))
KeepPag <- do.call(data.frame, KeepPag)

#Add ET0
SaveE0=data.frame(SaveE0)
my=which(!is.na(match(SaveE0$SaveE0,ys1)))
KeepE0=(SaveE0[my,-1])
Ke0_means=colMeans(KeepE0)

kpe=match(KeepPag$catch,cnames)

KeepPag$aet=KAET_means[kpe]

KeepPag$AET.P=KeepPag$aet/KeepPag$pmean

#Plot of AET/P
mai=match(ppoints$outlets.y,KeepPag$catch)
ppoints$AET.P1=KeepPag$AET.P[mai]

palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=AET.P1),color="black",alpha=1,shape=21,stroke=0)+ 
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

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/AridityIndex1.png", width=20, height=15, units=c("cm"),dpi=1500)



### ActEvapo/P ----

SaveAET=data.frame(SaveAET)
my=which(!is.na(match(SaveAET$SaveAET,ys1)))
KeepAET=(SaveAET[my,-1])
Ket_means=colMeans(KeepAET)


KeepPag$aet=Ket_means[kpe]

KeepPag$AI=KeepPag$aet/KeepPag$pmean
mai=match(ppoints$outlets.y,KeepPag$catch)
ppoints$AI1=KeepPag$AI[mai]

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

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/AET.P_P1.png", width=20, height=15, units=c("cm"),dpi=1500)


### rainfall concentration ----
myp=which(!is.na(match(SavePrecip$year,ys1)))

KeepP=SavePrecip[myp,c(1,3,4)]

KeepPag<- stats::aggregate(list(px=KeepP$val.gini), by = list(catch = KeepP$cat), 
                           FUN = function(x) c(pmean=mean(x, na.rm = T)))
KeepPag <- do.call(data.frame, KeepPag)

min(KeepPag$px)
mai=match(ppoints$outlets.y,KeepPag$catch)
ppoints$gini1=KeepPag$px[mai]

palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=gini1),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="Precipitation concentration", limits = c(0.5,1)) +
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

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/PGini1.png", width=20, height=15, units=c("cm"),dpi=500)


ls<-length(SeasonPrecip$cat)
SeasonPrecipmean=SeasonPrecip[seq(1,ls,2),]
SeasonPrecipSI=SeasonPrecip[-seq(1,ls,2),]
SeasonPrecipSI$SIchange=SeasonPrecipSI$Seasonality_2-SeasonPrecipSI$Seasonality_1
mai=match(ppoints$outlets.x, SeasonPrecipSI$cat)
ppoints$PSI=SeasonPrecipSI$Seasonality_1[mai]

palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=PSI),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="Precipitation seasonality", limits = c(0,.5)) +
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

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/PSI1.png", width=20, height=15, units=c("cm"),dpi=500)


ppoints$PSIch=SeasonPrecipSI$SIchange[mai]

palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=PSIch),color="black",alpha=1,shape=21,stroke=0.1)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="Precipitation seasonality change", limits = c(-0.1,.1)) +
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

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/PSI1.png", width=20, height=15, units=c("cm"),dpi=500)


