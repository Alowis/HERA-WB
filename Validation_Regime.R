source("functions_regime.R")
library(xts)
library(hydroGOF)
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
#outhybas<-outhybas[hybas2uc,]

ppl <- st_as_sf(UpArea, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

#remove Q obs that are outside
matcat=match(outhybas$latlong,UpArea$latlong)
UpArea=UpArea[matcat,]



ValidSf=read.csv(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Stations_ValidationF.csv")[,-1]

ValidSY=ValidSf[which(ValidSf$removal!="YES"),]
ValidSY=ValidSY[,c(1:16)]

#ValidSY is the final set of stations used
vEFAS=ValidSY[which(ValidSY$csource=="EFAS"),]
length(which(ValidSY$Rlen<(30*365)))


length(ValidSY$Var1[which(ValidSY$upa<200)])

#keep only records that match my outlets, or at least are on the same river nearby

### Distance between "official gauges" and EFAS points -----------------------
points=ValidSY[,c(3,5,1,2)]
Vsfloc=st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
outloc=st_transform(UpCat, crs = 4326)
#outloc=st_as_sf(UpCat, coords = c("Var1", "Var2"), crs = 4326)
#Sfloc=inner_join(points,outhybas,by=c("V1"="StationID"))
mOutSt=c()
for (i in 1:length(Vsfloc$upa)){
  #i=8
  r=Vsfloc$V1[i]
  cat(paste0(r,"\n"))
  v1=Vsfloc[which(Vsfloc$V1==r),]
  v2=outloc
  oula=st_distance(v2,v1)
  oula=oula/1000
  oula=as.numeric(oula)
  oulala<-oula[which(oula<50)]
  oulaID=v2[which(oula<50),]
  rup=oulaID$upa/Vsfloc$upa[i]
  rkup=which(rup<1.1 & rup>0.9)
  length(rkup)
  if (length(rkup>0)){
    print(rkup)
    keeploc=oulaID[rkup,]
    keeploc$station=r
    keeploc$dist=oulala[rkup]
    mOutSt=rbind(mOutSt,keeploc)
  } 

}


#cleaning of Moust
idrm=c()
for (is in unique(mOutSt$station)){
  mob<-which(mOutSt$station==is)
  if (length(mob)>1){
    print(mOutSt$dist[mob])
    tk<-which.min(mOutSt$dist[mob])
    irm=mob[-tk]
    idrm=c(idrm,irm)
  }
}
mOutSt1=mOutSt[-idrm,]
#plot the matching locations
#Plot parameters
cord.dec=ValidSY[,c(1,2)]
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

ppl <- mOutSt1
ppl <- st_transform(ppl, crs = 3035)
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,size=upa, fill=dist),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="distance to gauge",
    breaks=c(0,5,10), labels=c(0,5,10)) +
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


#ggsave("Plots/ValidStations_revisions.jpg", width=20, height=24, units=c("cm"),dpi=1500)

#load Q at outlets
Q <- fread(paste0(hydroDir,"/tss/HERA_Histo/disWin_1951_2020.csv"),header=TRUE)
time1=Q$V1
timeStampX=time1[order(time1)]
Q=Q[order(time1),]


cloc="290666"
River="Ticino"
cloc=as.character(mOutSt$outlets_y[1])
cupa=mOutSt$upa[which(mOutSt$outlets_y==as.numeric(cloc))]
Regimerun<-RegimeFunctionM(dvar=Q, cloc, yearstart=1951, wsize=70, timeStampX,Q=TRUE, cupa)
plot(Regimerun)

findout<-colnames(Q)
pptn=which(!is.na(match(colnames(Q),as.character(UpArea$outlets))))
mval=which(!is.na(match(colnames(Q),(as.character(mOutSt1$outlets_x)))))
Q_sim=Q[,..mval]
Q_sim=as.data.frame(Q_sim)

#match Qs with Moust
Qidm=mOutSt1[,c(28,35)]

main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
Q_data <- read.csv(paste0(valid_path,'out/Q_19502020.csv'), header = F)  # CSVs with observations
time=as.Date(as.numeric(Q_data$V2)-1,origin="0000-01-01")[-c(1,2)]
Q_data=Q_data[-1,-1]
Station_data_IDs <- unique(Qidm$station)
Outlets_data_IDs <- Qidm$outlets_y
Station_obs_IDs <- as.numeric(as.vector(t(Q_data[1, -1])))

Qsim_ids=as.numeric(colnames(Q_sim))
skills_cor=c()

mean_daily_value<- function (timeseries) 
{
  daily_timeseries <- apply.daily(timeseries, mean)
  tday = unique(as.Date(as.character(timeseries[, 1])))
  return(data.frame(date = tday, Qmd = daily_timeseries))
}

RegimeQ<-function(data){
  # Create a new column 'Qs' by multiplying the second column by 4*nbdays
  dt=30
  dt=dt1*dt2
  data$Qs <- tsEvaNanRunningMean(data[, 2], dt)
  names(data)[1]="datetime"
  
  
  # 3.1  Identify the first year in the data
  
  first_year <- year(min(data$datetime))
  
  # 3.2  Build a table of the 12 centre dates of that year
  
  centres_year1 <- tibble(
    centre = seq(
      from = as.POSIXct(paste0(first_year, "-01-15 00:00:00"), tz = "UTC"),
      by   = "1 month",
      length.out = 12
    )
  ) %>%
    mutate(
      month = month(centre),
      day   = yday(centre)          # always 15, but we keep it for clarity
    )
  
  # 3.3  Pull the same month‑day from *all* years
  
  df_center_all_years <- data %>%
    mutate(
      yr   = year(datetime),
      dy   = yday(datetime),
      mon  = month(datetime),
    ) %>%
    inner_join(centres_year1, by = c("dy" = "day")) %>%
    # At this point we have every observation that falls on the 15th of a month.
    # If you need the *exact* 00:00 time, keep only the nearest 6‑hour slot:
    group_by(yr, dy, mon) %>% 
    mutate(
      dist = abs(as.numeric(difftime(datetime, 
                                     as.POSIXct(paste0(yr, "-",
                                                       sprintf("%02d", mon), "-15 00:00:00"),
                                                tz = "UTC"),
                                     units = "secs")))
    ) %>%
    slice_min(order_by = dist, n = 1, with_ties = FALSE) %>%  # one group per year‑month
    ungroup()
  
  
  
  
  # Select only the first and third columns
  data <- data.frame(df_center_all_years[, c(1, 3)])
  # Apply the RegimeFast function to the processed data
  Regime <- RegimeFast(data)
  results=list(data=data,Regime=Regime)
  return(results)
  
}
#continue with a simple look of kge on monthly regime to show perfomances
for (s in 1:length(Station_data_IDs)){
  Station=Station_data_IDs[s]
  #Station="6457010"
  # s=which(Station_data_IDs==Station)
  Outlet<-Outlets_data_IDs[s]
  id_obs=match(Station,Station_obs_IDs)
  id_sim=match(Outlet,Qsim_ids)
  Station_obs_IDs[id_obs]
  uparea=ValidSY$upa[which(ValidSY$V1==Station)]
  if (length(uparea)>0){
    print(s)
    colnames(Q_sim)[id_sim]
    HERA_loc=Q_sim[,id_sim]
    
    HERA_Q<-data.frame(timeStampX,HERA_loc)
    HERA_qd<-mean_daily_value(HERA_Q)
    
    Q_loc=as.numeric(Q_data[-c(1:366),id_obs+1])
    #aggregation to monthly Q
    #regime
    timeStampD<-unique(as.Date(timeStampX))
    Qobs=data.frame(timeStampD,Q_loc)
    Robs=RegimeQ(Qobs)
    # plot(HERA_qd)
    #remove days where obs are NAs for a fair comparison
    HERA_qd$HERA_loc[which(is.na(Qobs$Q_loc))]=NA
    Rsim=RegimeQ(HERA_qd)
    Robs2=RegimeFast(Qobs)
    Rsim2=RegimeFast(HERA_qd)
    # plot(Rsim$data$Qs[c(240:280)], type="o")
    # lines(Robs$data$Qs[240:280], col=2)
    # 
    # plot(Rsim$Regim,type="o")
    # lines(Robs$Regime,col=2)
    #same for Qsim
    kge_Qmon=KGE(Rsim$data$Qs,Robs$data$Qs, na.rm=TRUE, method="2012",out.type="full")
    kge_Reg=KGE(Rsim$Regime$mean,Robs$Regime$mean, na.rm=TRUE, method="2012",out.type="full")
    kge_Reg2=KGE(Rsim2$mean,Robs2$mean, na.rm=TRUE, method="2012",out.type="full")
    skills_cor=rbind(skills_cor,c(Station,kge_Qmon$KGE.value,kge_Qmon$KGE.elements,
                                  kge_Reg$KGE.value,kge_Reg$KGE.elements,kge_Reg2$KGE.value,
                                  kge_Reg2$KGE.elements))
  }
}

hist(skills_cor[,7], breaks=50)
mean(skills_cor[,7])
length(which(skills_cor[,7]>0.5))/length(skills_cor[,7])
plot(skills_cor[,4],skills_cor[,6],xlim=c(-0,2))


#look at locs with shitty KGE
strangishit<-skills_cor[,1][which(skills_cor[,7]<0.2)]

dcheck<-ValidSY[which(!is.na(match(ValidSY$V1,strangishit))),]

ddcheck<-mOutSt1[which(!is.na(match(mOutSt1$station,strangishit))),]

ppl <- ddcheck
ppl <- st_transform(ppl, crs = 3035)
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,size=upa, fill=dist),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="distance to gauge",
    breaks=c(0,5,10,15,20), labels=c(0,5,10,15,20)) +
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



#Use the matched stations for show off
mOutStSave=mOutSt1[,c(1,19:27,34:36)]
mOutStSave=inner_join(mOutStSave,ValidSY, by=c("station"="V1"))

#add performances
skills_cor=data.frame(skills_cor)

plot(skills_cor$r.1,skills_cor$r.2)
mean(skills_cor$V10)
mean(skills_cor$V6)

skillsave<-skills_cor[,c(1,2,3,6,7,10,11)]
colnames(skillsave)=c("station","KGE.qmon","r.qmon","KGE.Regmon","r.Regmon","KGE.Regd","r.Regd")
mOutStSave=inner_join(mOutStSave,skillsave, by=c("station"))
write.csv(mOutStSave,file="D:/tilloal/Documents/01_Projects/RegimeShifts/outlets_stations.csv")

## Figure S2 Upstream area--------------------------------------

### Figure S2.a Map of upstream area--------------------------------------

min(mOutStSave$upa.y)
length(which(mOutStSave$upa.y<=25000))/length(mOutStSave$HYBAS_ID)
#Plot parameters
cord.dec=ValidSYl[,c(1,2)]
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

ppl <- st_as_sf(ValidSYl, coords = c("Var1", "Var2"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,fill=Rlen,size=upa),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
    breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
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


ggsave("Plots/ValidStations_revisions.jpg", width=20, height=24, units=c("cm"),dpi=1500)
