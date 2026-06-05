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


# --------------------------------------------------------------------------- #
# 1.  Paths                                                                   #
# --------------------------------------------------------------------------- #
hydroDir  <- "D:/tilloal/Documents/LFRuns_utils/data/"
regimeDir <- "D:/tilloal/Documents/01_Projects/RegimeShifts/"
gpkg_path <- paste0(regimeDir, "data/catchments_analysis_final_v3.gpkg")

###load UpArea -----
#load upstream area
hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data/")
outletname="upArea_European_01min.nc"
UpArea=UpAopen(hydroDir,outletname,outhybas)
head(UpArea)

#keep only rivers in EU domain
out1=outletopen(hydroDir,"efas_rnet_100km_01min")
out1$latlong=paste(round(out1$Var1,4),round(out1$Var2,4),sep=" ")
outhybas=inner_join(out1,outhybas, by="latlong")
outhybas$idlalo=paste(outhybas$idlo, outhybas$idla, sep=" ")

#remove Q obs that are outside
matcat=match(outhybas$latlong,UpArea$latlong)
UpArea=UpArea[matcat,]

which(is.na(matcat))
# Column selector: positions in CSV files matching our outlet IDs
# (read one header row to get the column names)
header_ref <- fread(paste0(hydroDir, "tss/HERA_Histo/SnowUpsX_1951_2020.csv"),
                    nrows = 0, header = TRUE)
matcol <- match(UpArea_full$outlets, as.numeric(colnames(header_ref)))
matcol <- matcol[!is.na(matcol)]
cnames <- as.character(UpArea$outlets) 

which(is.na(matcol))
message("Reading catchments_analysis_final.gpkg ...")
catchments_gpkg <- st_read(gpkg_path, quiet = TRUE)

catchments_plot=catchments_gpkg
#1 Components of Q ----

## Atmospheric ----
ActEvapo<- fread(paste0(hydroDir,"/tss/HERA_Histo/ActEvapo_1951_2020.csv"),header=TRUE)
time=ActEvapo$V1
timeStampX=time[order(time)]
ActEvapo=ActEvapo[order(time),]

matcol=match(UpArea$outlets,as.numeric(colnames(ActEvapo)))
cnames=names(ActEvapo)[matcol]

ActEvapo <- ActEvapo[, .SD, .SDcols = matcol]

Rain<- fread(paste0(hydroDir,"/tss/HERA_SocCF/RainUpsX_1951_2020.csv"),header=TRUE)
Rain=Rain[order(time),]

Snowmelt<- fread(paste0(hydroDir,"/tss/HERA_Histo/snowMeltUpsX_1951_2020.csv"),header=TRUE)
Snowmelt=Snowmelt[order(time),]

Rain <- Rain[, .SD, .SDcols = matcol]
Snowmelt <- Snowmelt[, .SD, .SDcols = matcol]

## Runoff and Q ----

#Runoff
Runoff<- fread(paste0(regimeDir,"/data/surfaceRunoffUpsX_nested_1951_2020.csv"),header=TRUE)
time=Runoff$V1
timeStampX=time[order(time)]
Runoff=Runoff[order(time),]


runoff_dt <- copy(Runoff)
runoff_dt[, year := year(timeStampX)]

# Sum within each year, then average across years -> mean annual precipitation
yearly_runoff <- runoff_dt[, lapply(.SD, sum, na.rm = TRUE),
                       by = year,
                       .SDcols = cnames]

mean_annual_runoff <- colMeans(yearly_runoff[, .SD, .SDcols = cnames],
                             na.rm = TRUE)

#map of runoff

catchments_plot$mean_runoff<- mean_annual_runoff[
  match(as.numeric(catchments_plot$catch_id), as.numeric(names(mean_annual_runoff)))
]

# -- 3. Plot ----------------------------------------------------------------
palet <- hcl.colors(11, palette = "Roma", rev = F)

ggplot(basemap) +
  geom_sf(fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data  = catchments_plot,
          aes(geometry = geom,
              fill     = mean_runoff),
          color = "grey20", shape = 21, stroke = 0.2, alpha = 0.9) +
  scale_fill_gradientn(
    colors = palet,
    oob    = scales::squish,
    name   = "Mean annual\nrunoff (mm/y)",
    breaks = seq(0, 2000, by = 10),
    limits=c(0,100)
  ) +
  # scale_size(
  #   range  = c(0.8, 5),
  #   trans  = "sqrt",
  #   name   = expression(paste("Upstream area (", km^2, ")")),
  #   breaks = c(100, 1000, 10000, 100000, 500000),
  #   labels = c("100", "1 000", "10 000", "100 000", "500 000")
  # ) +
  coord_sf(xlim = range(nco[, 1]), ylim = range(nco[, 2])) +
  labs(
    title    = "Mean annual direct runoff",
    subtitle = "Average over 1951–2020, deaggregated to inter-catchment area",
    x = "Longitude", y = "Latitude"
  ) +
  guides(
    fill = guide_colourbar(barwidth = 0.5, barheight = 18, reverse = FALSE),
    size = "none"
  ) +
  theme(
    plot.title       = element_text(size = 14, face = "bold"),
    plot.subtitle    = element_text(size = 11, colour = "grey40"),
    axis.title       = element_text(size = tsize),
    panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
    panel.border     = element_rect(linetype = "solid", fill = NA, colour = "black"),
    panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
    panel.grid.minor = element_line(colour = "grey90"),
    legend.title     = element_text(size = tsize),
    legend.text      = element_text(size = osize),
    legend.position  = "right",
    legend.key       = element_rect(fill = "transparent", colour = "transparent"),
    legend.key.size  = unit(0.8, "cm")
  )

ggsave(paste0(regimeDir, "plots/mean_runoff_vf.png"),
       width = 22, height = 16, units = "cm", dpi = 300)



#Q out of upper zone
quz<- fread(paste0(hydroDir,"/tss/HERA_Histo/qUzUpsX_1951_2020.csv"),header=TRUE)
quz=quz[order(time),]

#Q from upper to lowe zone
qutl<- fread(paste0(hydroDir,"/tss/HERA_Histo/percUZLZUpsX_1951_2020.csv"),header=TRUE)
qutl=qutl[order(time),]

#Q out of lower zone
qlz<- fread(paste0(hydroDir,"/tss/HERA_Histo/qLzUpsX_1951_2020.csv"),header=TRUE)
qlz=qlz[order(time),]

#Q from rain and snowmelt to soil
Infil<- fread(paste0(hydroDir,"/tss/HERA_Histo/InfUpsX_1951_2020_v2.csv"),header=TRUE)
timeI=as.POSIXct(Infil$V1)
Infil=Infil[order(timeI),]
timeStampI=timeI[order(timeI)]

Runoff <- Runoff[, .SD, .SDcols = matcol]
quz <- quz[, .SD, .SDcols = matcol]
qlz <- qlz[, .SD, .SDcols = matcol]
qutl <- qutl[, .SD, .SDcols = matcol]
Infil <- Infil[, .SD, .SDcols = matcol]
#Q from soil to upper zone

qStGW<-fread(paste0(hydroDir,"/tss/HERA_Histo/dSubToUzUpsX_1951_2020.csv"),header=TRUE)
qStGW=qStGW[order(time),]

#Q from rain and snowmet to upper zone
Pflow<-fread(paste0(hydroDir,"/tss/HERA_Histo/prefFlowUpsX_1951_2020.csv"),header=TRUE)
Pflow=Pflow[order(time),]

qStGW <- qStGW[, .SD, .SDcols = matcol]
Pflow <- Pflow[, .SD, .SDcols = matcol]

#GW loss
GWloss<-fread(paste0(hydroDir,"/tss/HERA_Histo/lossUpsX_1951_2020.csv"),header=TRUE)
GWloss=GWloss[order(time),]

GWloss <- GWloss[, .SD, .SDcols = matcol]
#Q
Q <- fread(paste0(hydroDir,"/tss/HERA_Histo/disWin_1951_2020.csv"),header=TRUE)
Q=Q[order(time),]
# 
# #remove Q obs that are outside
# matcat=match(outhybas$latlong,UpArea$latlong)
# UpArea=UpArea[matcat,]
# 
# matcol=match(UpArea$outlets,as.numeric(colnames(Q)))
# cnames=names(Q)[matcol]
# plot(matcol)
# 
Q2 <- Q[, .SD, .SDcols = matcol]
rm(Q)
gc()

#Qe creation ----

#Blue water
Qe=Runoff+quz+qlz

#save blue water:
write.csv(Qe,file="D:/tilloal/Documents/01_Projects/RegimeShifts/data/BlueWater.csv")

#Green water
#water which does not go into the soil or rivers is intercepted and available for evaporation
EvAv=Rain+Snowmelt-Pflow-Runoff-Infil
#evaporation from soil is act evapo minus intercepted water evapo
S2At=ActEvapo-EvAv
#Soil balace in infiltration - soil evapo - soil to GW
Se= Infil - S2At - qStGW

#save green water:
write.csv(Se,file="D:/tilloal/Documents/01_Projects/RegimeShifts/data/GreenWater.csv")

#look at a small catchment

cloc="303662"
Qloc=Q2[[cloc]]
Qloc=data.frame(timeStampX,Qloc)
cupa=UpArea$upa[which(UpArea$outlets==303662)]
Qloc$Qs=Qloc$Qloc/(cupa*1000000)*1000*3600*24
Qloc=Qloc[,c(1,3)]

Qloc$daily=tsEvaNanRunningMean(Qloc[,2], windowSize =1/dt)

data=Qloc
names(data)[c(2,3)]=c("Qs","Qs")
names(data)[1]="date"
dsel=hour(data$date)
dataD=data[which(dsel==12 | dsel==13),c(1,3)]
plot(dataD)

dataD[,2]=tsEvaNanRunningMean(dataD[,2],30)

RegimeQd=RegimeFast(Qloc)
plot(RegimeQd)


#Same shit for sum of Q component
Qeloc=Qe[[cloc]]
Qeloc=data.frame(timeStampX,Qeloc)
Qeloc$Qs=Qeloc$Qeloc*4
Qeloc=Qeloc[,c(1,3)]

RegimeQe=RegimeFast(Qeloc)
plot(RegimeQe, type="l")
lines(RegimeQd, col=2)

cloc="303662"

Qlzloc=qlz[[cloc]]
Qlzloc=data.frame(timeStampX,Qlzloc)
Qlzloc$Qs=Qlzloc$Qlzloc*4
Qlzloc=Qlzloc[,c(1,3)]

Rqlz=RegimeFast(Qlzloc)
plot(Rqlz)
Quzloc=quz[[cloc]]
Quzloc=data.frame(timeStampX,Quzloc)
Quzloc$Qs=Quzloc$Quzloc*4
Quzloc=Quzloc[,c(1,3)]

Qulzloc=Qlzloc
Qulzloc[,2]=Qlzloc[,2]+Quzloc[,2]

Rqz=RegimeFast(Qulzloc)

Qrzloc=Runoff[[cloc]]
Qrzloc=data.frame(timeStampX,Qrzloc)
Qrzloc$Qs=Qrzloc$Qrzloc*4
Qrzloc=Qrzloc[,c(1,3)]

Qrloc=Qrzloc
Qrloc[,2]=Qulzloc[,2]+Qrzloc[,2]
Rqzr=RegimeFast(Qrloc)

RegimeQd$qzr=Rqzr$mean
RegimeQd$d=RegimeQd$mean-RegimeQd$qzr
qlim=c(0,1.2*max(RegimeQ$mean))
plot(Rqlz$date, Rqlz$mean, type="n", axes=FALSE,ylim=qlim,xaxs="i",yaxs="i",
     xlab = NA, ylab = expression(paste("Q (mm/d)")))

polygon(c(Rqlz$date,rev(Rqlz$date)),c(rep(0,length(Rqlz$mean)),rev(Rqlz$mean)),
        col=alpha("lightblue",.7),border="transparent")
polygon(c(Rqlz$date,rev(Rqlz$date)),c(Rqlz$mean,rev(Rqz$mean)),
        col=alpha("royalblue",.7),border="transparent")
polygon(c(Rqz$date,rev(Rqz$date)),c(Rqz$mean,rev(Rqzr$mean)),
        col=alpha("purple",.7),border="transparent")
lines(RegimeQ$date, RegimeQ$mean, col="black",lwd=1)
abline(v = format(mois.deb,"%j"), col="lightgrey", lty=3)
axis(2)
axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=1)
box()
legend("topleft", leg=c("baseflow","sub-surface flow","runoff","Q"),
       lwd=c(4,4,4,2), col=c("lightblue","royalblue","purple","black"),
       cex=1, lty=c(1,1,1,1),bty="n")

#Ok it matches perfectly on small catchments

#comparison of regimes between two periods

#regime of each component of Q ----

#Ardeche
cloc="303662"
River="Ardeche"

#Ticino river
cloc="290666"
River="Ticino"
#Sieg
cloc="196109"
River="Sieg"
#Var
cloc="311057"
River="Var"

#thames @ Reading
cloc="184737"
River="Upper Thames"

#raising issue with handling of GW losses
#Rio Zujar
cloc="362849"
River="Rio Zujar"
#Upper Arno
cloc="313396"
River="Upper Arno"

#Vistula
cloc="204769"
River="Upper Vistula"


#load rivers from validation


#inputs---- 

##rain----
Regimerain<-RegimeFunctionM(dvar=Rain, cloc, yearstart=1951, wsize=70, timeStampX)
plot(Regimerain)

##snowmelt----
Regimesm<-RegimeFunctionM(dvar=Snowmelt, cloc, yearstart=1951, wsize=70, timeStampX)
plot(Regimesm)

#outputs to atmosphere----

##Actual evapotranpiration----
Regimeaet<-RegimeFunctionM(dvar=ActEvapo, cloc, yearstart=1951, wsize=70, timeStampX)
plot(Regimeaet)

# soil state variables----

## Q from runoff ----
Regimerun<-RegimeFunctionM(dvar=Runoff, cloc, yearstart=1951, wsize=70, timeStampX)
plot(Regimerun)

## Q from lower zone ----
Regimeqlz<-RegimeFunctionM(dvar=qlz, cloc, yearstart=1951, wsize=70, timeStampX)
plot(Regimeqlz)

## Q from upper zone ----
Regimequz<-RegimeFunctionM(dvar=quz, cloc, yearstart=1951, wsize=70, timeStampX)
plot(Regimequz)

## Q from upper to lower zone ----
Regimequtl<-RegimeFunctionM(dvar=qutl, cloc, yearstart=1951, wsize=70, timeStampX)
plot(Regimequtl)

## Q from soil to GW ----
RegimeqStGW<-RegimeFunctionM(dvar=qStGW, cloc, yearstart=1951, wsize=70, timeStampX)
plot(RegimeqStGW)

## Q from rain+sm to GW----
RegimePflow<-RegimeFunctionM(dvar=Pflow, cloc, yearstart=1951, wsize=70, timeStampX)
plot(RegimePflow)

##Q from rain+SM to soil
RegimeInf<-RegimeFunctionM(dvar=Infil, cloc, yearstart=1951, wsize=70, timeStampX)
plot(RegimeInf)

## Q out ----
RegimeQ<-RegimeFunctionM(dvar=Q2, cloc, yearstart=1951, wsize=70, timeStampX,cupa=cupa,Q=TRUE)
#RegimeQ$mean=RegimeQ$mean/(cupa*1000000)*1000*3600*24
plot(RegimeQ$mean)
RegimeQr=Regimeqlz$mean+Regimequz$mean+Regimerun$mean
lines(RegimeQr)
##GWloss
Regimeloss<-RegimeFunctionM(dvar=GWloss, cloc, yearstart=1951, wsize=70, timeStampX)
plot(Regimeloss)
sum(Regimeloss$mean)
#now I add up different components to reproduce discharge
#rain + snowelt
Regimersm=Regimerain
Regimersm$mean=Regimerain$mean+Regimesm$mean
Regimersm$sm2p=Regimesm$mean/Regimersm$mean
Regimersm$r2p=Regimerain$mean/Regimersm$mean
plot(Regimersm)

#I remove evapotranspiration

RegimeAll=Regimersm
RegimeAll$atmo=Regimeaet$mean
plot(RegimeAll$atmo)

RegimeAll$pflow=RegimePflow$mean
RegimeAll$infil=RegimeInf$mean
RegimeAll$runoff=Regimerun$mean
RegimeAll$gwre=RegimeqStGW$mean
RegimeAll$utol=Regimequtl$mean
RegimeAll$qlz=Regimeqlz$mean
RegimeAll$quz=Regimequz$mean
RegimeAll$qav=RegimeAll$mean-RegimeAll$infil
RegimeAll$qavsm=RegimeAll$qav*Regimersm$sm2p
RegimeAll$qavr=RegimeAll$qav*Regimersm$r2p
RegimeAll$loss=Regimeloss$mean
sum(RegimeAll$lzb)


RegimeAll$evapavail=RegimeAll$mean-RegimeAll$pflow-RegimeAll$runoff-RegimeAll$infil
plot(RegimeAll$evapavail)
RegimeAll$atmob=RegimeAll$atmo-RegimeAll$evapavail
plot(RegimeAll$atmo)
lines(RegimeAll$atmob)

#soil balance= infiltration-evapotranspiration
RegimeAll$soilQ=RegimeAll$atmob
sum(RegimeAll$soilQ)
sum(RegimeAll$infil)
sum(RegimeAll$gwre)
#Evap from interception 
RegimeAll$intev=RegimeAll$evapavail
plot(RegimeAll$intev)
#upper zone balance= pflow+gwr-quz-utl
RegimeAll$uzb=RegimeAll$pflow+RegimeAll$gwre-RegimeAll$utol-RegimeAll$quz
#lower zone balance= utl-qlz
RegimeAll$lzb=RegimeAll$utol-RegimeAll$qlz-RegimeAll$loss

#what is really available for exit
RegimeAll$lzb2=RegimeAll$utol-RegimeAll$loss

#correction of GWperc with regard to GWloss
#RegimeAll$utol[which(RegimeAll$utol<RegimeAll$loss)]=RegimeAll$loss[which(RegimeAll$utol<RegimeAll$loss)]

RegimeAll$lzb=RegimeAll$utol-RegimeAll$qlz-RegimeAll$loss

RegimeAll$gwb=RegimeAll$pflow+RegimeAll$gwre-RegimeAll$qlz-RegimeAll$quz-RegimeAll$loss
sum(RegimeAll$lzb)
#upper zone out= quz-qutl
RegimeAll$uzqout=RegimeAll$quz-RegimeAll$utol


#1 plot for evapotranpiration balance and 1 for Q?


#GGplot

month=c(1:12)

Qcompo2=data.frame(month,pr=RegimeAll$mean,
                   rain=RegimeAll$qavr,
                   snow=RegimeAll$qavsm,
                   intev=RegimeAll$intev,
                   gwre=-RegimeAll$gwre,
                   sevap=RegimeAll$soilQ,
                   run=RegimeAll$runoff,
                   bf=RegimeAll$qlz,
                   flow=RegimeAll$quz,
                   infil=RegimeAll$infil,
                   loss=RegimeAll$loss)

Qcompo2$sout= rowSums(Qcompo2[ , c(6,8:12)], na.rm = TRUE)

Qcompo2$intev=-Qcompo2$intev
Qcompo2$sevap=-Qcompo2$sevap
library(tidyverse)   # pulls in ggplot2, tidyr, dplyr, etc.

df_long <- Qcompo2 %>%
  pivot_longer(cols = c(infil,snow,rain, intev, sevap, gwre),
               names_to   = "variable",
               values_to  = "value")



df_total <- Qcompo2 %>%
  mutate(total = infil + sevap + gwre) %>%
  select(month, total)

sum(df_total$total)
desired_order <- (c("snow", "rain","infil","intev","sevap", "gwre"))

df_long <- df_long %>% 
  mutate(
    variable = factor(variable, levels = desired_order),
    month    = factor(month, levels = sort(unique(month)))
  )

head(df_long)

my_fill_pal <- c(rain = "steelblue2",snow="purple",
                 infil="springgreen3", intev="palevioletred",sevap = "orangered", gwre="brown")

legend_labels <- c(rain = "rainfall", snow = "snowmelt", 
                   infil = "infiltration", intev ="evaporation from interception",
                   sevap="evapotranspiration", gwre="groundwater recharge")



Qcompo2$rain+Qcompo2$snow+Qcompo2$intev-Qcompo2$bf-Qcompo2$run-Qcompo2$flow-Qcompo2$loss

Qcompo2$lzb=-RegimeAll$lzb
#Remove storage moments?
Qcompo2$lzb[which(Qcompo2$lzb>Qcompo2$bf)]=Qcompo2$bf[which(Qcompo2$lzb>Qcompo2$bf)]
Qcompo2$bfgw=Qcompo2$lzb
Qcompo2$bfp=Qcompo2$bf-Qcompo2$lzb

#Only show positive contributions, GW storage is unreliable
Qcompo2$bfgw[which(Qcompo2$bfgw<0)]=0
Qcompo2$bfp[which(Qcompo2$lzb<0)]=Qcompo2$bf[which(Qcompo2$lzb<0)]

# 
Qcompo2$sflow=-Qcompo2$gwre

Qcompo2$flowp=Qcompo2$flow-Qcompo2$sflow

#if flow from soil is higher than quz, flow from soil is going to lower zone
Qcompo2$bfp[which(Qcompo2$flowp<0)]=Qcompo2$bfp[which(Qcompo2$flowp<0)]+Qcompo2$flowp[which(Qcompo2$flowp<0)]
Qcompo2$flowp[which(Qcompo2$flowp<0)]=0

desired_order2 <- rev(c("bfgw", "bfp","sflow","flowp","run"))

df_long2 <- Qcompo2 %>%
  pivot_longer(cols = c(bfgw,bfp,sflow,flowp, run),
               names_to   = "variable",
               values_to  = "value")

df_long2 <- df_long2 %>% 
  mutate(
    variable = factor(variable, levels = desired_order2),
    month    = factor(month, levels = sort(unique(month)))
  )
head(df_long2)

stopifnot(identical(levels(df_long2$variable), desired_order))
stopifnot(identical(levels(df_long2$month), levels(df_total2$month)))

df_total2 <- Qcompo2 %>%
  mutate(total = bfgw + bfp + sflow + flowp + run ) %>%
  select(month, total)



ggplot(df_long,
       aes(x = factor(month),          # treat month as a categorical axis
           y = value,
           fill = factor(variable))) +        # colour by variable (bf, rf, frf)
  geom_col() +    
  geom_line(data = df_total,
            aes(x = factor(month), y = total, group = 1, colour = "soil_balance"),
            size = 1.2,
            inherit.aes = FALSE) +            # line
  geom_point(data = df_total,
             aes(x = factor(month), y = total, colour = "soil_balance"),
             size = 3,
             inherit.aes = FALSE) +  
  scale_fill_manual(values = my_fill_pal, 
                    breaks = desired_order, labels=legend_labels) +
  scale_colour_manual(name = "",                # legend title for the line
                      values = c(soil_balance = "black")) +
  labs(x = "Months",
       y = "Q (mm)",
       fill = " ",
       title = paste0("Green water regime for river ",River)) +
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





my_fill_pal2 <- c(bfgw = "brown",
                  bfp="brown2", sflow= "springgreen4",flowp="steelblue1",run = "purple")

legend_labels2 <- c(bfgw = "baseflow from GW", bfp = "baseflow from pr", sflow="flow from soil",
                    flowp = "sub-surface flow", run ="runoff")


ggplot(df_long2,
       aes(x = factor(month),          # treat month as a categorical axis
           y = value,
           fill = variable)) +        # colour by variable (bf, rf, frf)
  geom_col() +    
  geom_line(data = df_total2,
            aes(x = factor(month), y = total, group = 1, colour = "Total"),
            size = 1.2,
            inherit.aes = FALSE) +            # line
  geom_point(data = df_total2,
             aes(x = factor(month), y = total, colour = "Total"),
             size = 3,
             inherit.aes = FALSE) +            # points on the line
  scale_fill_manual(values = my_fill_pal2, 
                    breaks = desired_order2, labels=legend_labels2) +
  scale_colour_manual(name = "",                # legend title for the line
                      values = c(Total = "black")) +
  labs(x = "Months",
       y = "Q (mm)",
       fill = " ",
       title = paste0("Blue water regime for river ", River)) +
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


#Ok, save both plots for 2-3 rivers and show to Luc


#Now compare both results and real Q

plot(df_total$total,type="l")
plot(df_total2$total,col=2)
lines(RegimeQ$mean,col=3)


# OLD CODE TO be reviewed -----

#I now have all components of Q
#GW is contributing to support the discharge in summer month by adding to baseflow and uzflow

#changes in regimes 
Regimeqlz1<-RegimeFunctionM(dvar=qlz, cloc, yearstart=1951, wsize=70, timeStampX)
plot(Regimeqlz1)

Regimeqlz2<-RegimeFunctionM(dvar=qlz, cloc, yearstart=1991, wsize=29, timeStampX)
plot(Regimeqlz2)
Regimeqlz2$dq=Regimeqlz2$mean-Regimeqlz1$mean
plot(Regimeqlz2$date,Regimeqlz2$dq,type="o")
abline(h=0,lwd=3)


Regimequz1<-RegimeFunctionM(dvar=quz, cloc, yearstart=1951, wsize=29, timeStampX)
plot(Regimequz1)

Regimequz2<-RegimeFunctionM(dvar=quz, cloc, yearstart=1991, wsize=29, timeStampX)
plot(Regimequz2)
Regimequz2$dq=Regimequz2$mean-Regimequz1$mean
plot(Regimequz2$date,Regimequz2$dq, type="o")
abline(h=0,lwd=3)

Regimeqr1<-RegimeFunctionM(dvar=Runoff, cloc, yearstart=1951, wsize=29, timeStampX)
plot(Regimeqr1)

Regimeqr2<-RegimeFunctionM(dvar=Runoff, cloc, yearstart=1991, wsize=29, timeStampX)
plot(Regimeqr2)
Regimeqr2$dq=Regimeqr2$mean-Regimeqr1$mean
plot(Regimeqr2$date,Regimeqr2$dq, type="o")
abline(h=0,lwd=3)

#GGplot

month=c(1:12)

CompoRchange=data.frame(month,bf=Regimeqlz2$dq,rf=Regimequz2$dq,frf=Regimeqr2$dq)

CompoR1=data.frame(month,bf=Regimeqlz1$mean,rf=Regimequz1$mean,frf=Regimeqr1$mean)

library(tidyverse)   # pulls in ggplot2, tidyr, dplyr, etc.

df_long <- CompoRchange %>%
  pivot_longer(cols = c(bf, rf, frf),
               names_to   = "variable",
               values_to  = "value")

head(df_long)

df_total <- CompoRchange %>%
  mutate(total = bf + rf + frf) %>%
  select(month, total)

ggplot(df_long,
       aes(x = factor(month),          # treat month as a categorical axis
           y = value,
           fill = variable)) +        # colour by variable (bf, rf, frf)
  geom_col() +    
  geom_line(data = df_total,
            aes(x = factor(month), y = total, group = 1, colour = "Total"),
            size = 1.2,
            inherit.aes = FALSE) +            # line
  geom_point(data = df_total,
             aes(x = factor(month), y = total, colour = "Total"),
             size = 3,
             inherit.aes = FALSE) +            # points on the line
  scale_fill_brewer(palette = "Set2") +
  scale_colour_manual(name = "",                # legend title for the line
                      values = c(Total = "black")) +
  labs(x = "Month",
       y = "Value",
       fill = "Series",
       title = "Change in Regime") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Same plot with rain and snow

ActEvapo<- fread(paste0(hydroDir,"/tss/HERA_Histo/ActEvapo_1951_2020.csv"),header=TRUE)
ActEvapo=ActEvapo[order(time),]

# Rain<- fread(paste0(hydroDir,"/tss/HERA_SocCF/RainUpsX_1951_2020.csv"),header=TRUE)
# Rain=Rain[order(time),]
# 
# Snow<- fread(paste0(hydroDir,"/tss/HERA_Histo/SnowUpsX_1951_2020.csv"),header=TRUE)
# Snow=Snow[order(time),]

Snowmelt<- fread(paste0(hydroDir,"/tss/HERA_Histo/snowMeltUpsX_1951_2020.csv"),header=TRUE)
Snowmelt=Snowmelt[order(time),]

Infil<- fread(paste0(hydroDir,"/tss/HERA_Histo/InfUpsX_1951_2020.csv"),header=TRUE)
timeI=as.POSIXct(Infil$V1)
Infil=Infil[order(timeI),]
timeStampI=timeI[order(timeI)]

ActEvapo <- ActEvapo[, .SD, .SDcols = matcol]
#Rain <- Rain[, .SD, .SDcols = matcol]
Snowmelt <- Snowmelt[, .SD, .SDcols = matcol]
Infil <- Infil[, .SD, .SDcols = matcol]




RegimeSM1<-RegimeFunctionM(dvar=Snowmelt, cloc, yearstart=1951, wsize=29, timeStampX)
plot(RegimeSM1)

RegimeSM2<-RegimeFunctionM(dvar=Snowmelt, cloc, yearstart=1991, wsize=29, timeStampX)
plot(RegimeSM2)
RegimeSM2$dq=RegimeSM2$mean-RegimeSM1$mean
plot(RegimeSM2$date,RegimeSM2$dq,type="o")
abline(h=0,lwd=3)


RegimeRain1<-RegimeFunctionM(dvar=Rain, cloc, yearstart=1951, wsize=29, timeStampX)


RegimeRain2<-RegimeFunctionM(dvar=Rain, cloc, yearstart=1991, wsize=29, timeStampX)

RegimeRain2$dq=RegimeRain2$mean-RegimeRain1$mean
plot(RegimeRain2$date,RegimeRain2$dq, type="o")
abline(h=0,lwd=3)

RegimeAEP1<-RegimeFunctionM(dvar=ActEvapo, cloc, yearstart=1951, wsize=29, timeStampX)
plot(RegimeAEP1)

RegimeAEP2<-RegimeFunctionM(dvar=ActEvapo, cloc, yearstart=1991, wsize=29, timeStampX)
plot(RegimeAEP2)
RegimeAEP2$dq=RegimeAEP2$mean-RegimeAEP1$mean
plot(RegimeAEP2$date,RegimeAEP2$dq, type="o")
abline(h=0,lwd=3)


RegimeInf1<-RegimeFunctionM(dvar=Infil, cloc, yearstart=1951, wsize=29, timeStampX)
plot(RegimeInf1)

RegimeInf2<-RegimeFunctionM(dvar=Infil, cloc, yearstart=1991, wsize=29, timeStampX)
plot(RegimeInf2)
RegimeInf2$dq=RegimeInf2$mean-RegimeInf1$mean
plot(RegimeInf2$date,RegimeInf2$dq, type="o")
abline(h=0,lwd=3)




#Q vs variables

Rmx<-data.frame(month,value=RegimeSM1$mean+RegimeRain1$mean-(RegimeAEP1$mean),variable=rep("met",12))
Rq<-data.frame(month,value=Regimequz1$mean+Regimeqlz1$mean+Regimeqr1$mean,variable=rep("Q",12))
sres=sum(abs(Rq$value-Rmx$value))
crinf<-(Rq$value-Rmx$value)/sres
sinf=sum(RegimeInf1$mean)
InfBf<-crinf*sinf
sum(abs(InfBf))
sum(RegimeInf1$mean)
outflow1=(RegimeInf1$mean+InfBf)

Rmx<-data.frame(month,value=RegimeSM2$mean+RegimeRain2$mean-(RegimeAEP2$mean),variable=rep("met",12))
Rq<-data.frame(month,value=Regimequz2$mean+Regimeqlz2$mean+Regimeqr2$mean,variable=rep("Q",12))
sres=sum(abs(Rq$value-Rmx$value))
crinf=(Rq$value-Rmx$value)/sres
sinf=sum(RegimeInf2$mean)
InfBf<-crinf*sinf
sum(abs(InfBf))
sum(RegimeInf2$mean)
outflow2=(RegimeInf2$mean+InfBf)

outflowd=outflow2-outflow1
#GGplot
sum(outflowd)
-sum(RegimeInf2$dq)
sum(RegimeRain2$dq)
sum(RegimeSM2$dq)
-sum(RegimeAEP2$dq)
month=c(1:12)

CompoRchange=data.frame(month,sm=RegimeSM2$dq,rain=RegimeRain2$dq,aep=-RegimeAEP2$dq,inf=-RegimeInf2$dq,outf=outflowd)

CompoR1=data.frame(month,sm=RegimeSM1$mean,rain=RegimeRain1$mean,aep=-RegimeAEP1$mean,inf=-RegimeInf1$mean)


library(tidyverse)   # pulls in ggplot2, tidyr, dplyr, etc.

# 2.1  Put the data in “long” format

df_long <- CompoRchange %>%
  pivot_longer(cols = c(sm, rain, aep,inf,outf),
               names_to   = "variable",
               values_to  = "value")

head(df_long)

df_total <- CompoRchange %>%
  mutate(total = sm + rain + inf+ aep+outf) %>%
  select(month, total)

ggplot(df_long,
       aes(x = factor(month),          # treat month as a categorical axis
           y = value,
           fill = variable)) +        # colour by variable (bf, rf, frf)
  geom_col() +     
  geom_line(data = df_total,
            aes(x = factor(month), y = total, group = 1, colour = "Total"),
            size = 1.2,
            inherit.aes = FALSE) +            # line
  geom_point(data = df_total,
             aes(x = factor(month), y = total, colour = "Total"),
             size = 3,
             inherit.aes = FALSE) + 
  scale_colour_manual(name = "",                # legend title for the line
                      values = c(Total = "black")) +# geom_col = bar with height = value
  scale_fill_brewer(palette = "Set2") +   # nice colour palette (optional)
  labs(x = "Month",
       y = "Value",
       fill = "Series",
       title = "Change in Regime") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





#Q vs variables

Rmx<-data.frame(month,value=RegimeSM1$mean+RegimeRain1$mean-(RegimeAEP1$mean),variable=rep("met",12))

Rq<-data.frame(month,value=Regimequz1$mean+Regimeqlz1$mean+Regimeqr1$mean,variable=rep("Q",12))

sres=sum(abs(Rq$value-Rmx$value))
crinf<-(Rq$value-Rmx$value)/sres
sinf=sum(RegimeInf1$mean)
InfBf<-crinf*sinf
sum(InfBf)
outflow=(RegimeInf1$mean+InfBf)

sum(-RegimeInf1$mean)
plot(-RegimeInf1$mean, ylim=c(-50,50))
points(InfBf,col=2)
Rmc<-data.frame(month,value=RegimeSM1$mean+RegimeRain1$mean-(RegimeAEP1$mean)+InfBf,variable=rep("corrmet",12))


Rqq=rbind(Rmx,Rmc,Rq)
ggplot(Rqq,
       aes(x = factor(month),          # treat month as a categorical axis
           y = value,fill=variable)) +        # colour by variable (bf, rf, frf)
  geom_col(position = "dodge", width = 0.7) +                       # geom_col = bar with height = value
  scale_fill_brewer(palette = "Set2") +   # nice colour palette (optional)
  labs(x = "Month",
       y = "Value",
       fill = "Series",
       title = "Change in Regime") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



CompoR1=data.frame(month,sm=RegimeSM1$mean,rain=RegimeRain1$mean,aep=-RegimeAEP1$mean,inf=-RegimeInf1$mean,outflow=outflow)



library(tidyverse)   # pulls in ggplot2, tidyr, dplyr, etc.

df_long <- CompoR1 %>%
  pivot_longer(cols = c(sm, rain, aep,inf,outflow),
               names_to   = "variable",
               values_to  = "value")

head(df_long)

ggplot(df_long,
       aes(x = factor(month),          # treat month as a categorical axis
           y = value,
           fill = variable)) +        # colour by variable (bf, rf, frf)
  geom_col() +                       # geom_col = bar with height = value
  scale_fill_brewer(palette = "Set2") +   # nice colour palette (optional)
  labs(x = "Month",
       y = "Value",
       fill = "Series",
       title = "Change in Regime") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Q vs varia

AEloc=ActEvapo[[cloc]]
AEloc=data.frame(timeStampX,AEloc)
AEloc$Qs=AEloc[,2]*4
AEloc=AEloc[,c(1,3)]
RegimeAE=RegimeFast(AEloc)
plot(RegimeAE, type="l")



RegimeSM<-RegimeFunction(Snowmelt, cloc, timeStampX)
plot(RegimeSM, type="l")

RegimeR<-RegimeFunction(Rain, cloc, timeStampX)
plot(RegimeR, type="l")

RegimeInf<-RegimeFunction(Infil, cloc, timeStampI)
lines(RegimeInf, type="l")

plot(RegimeInf$mean/RegimeR$mean)
Rtest=data.frame(RegimeSM$date, RegimeSM$mean+RegimeR$mean-RegimeInf$mean)
RegimeQ$WB=RegimeSM$mean+RegimeR$mean-RegimeAE$mean

plot(RegimeQ$WB,type="o")
lines(RegimeQ$mean,col=2)
dinf=(RegimeQ$mean-RegimeQ$WB)

plot(dinf)
lines(absd,col=2)
RQ=RegimeQ
RQ$mean=RegimeQ$mean+RegimeInf$mean-RegimeAE$mean
plot(Rtest, type="l")
lines(RegimeQ,col=2)
lines(RQ,col=3)


cor(Rtest$RegimeSM.mean...RegimeR.mean...RegimeInf.mean,RegimeQ$mean)
Rsnow=RegimeFast(TSnowf)
Rrain=RegimeFast(TRain)
RActE=RegimeFast(TAct)
RActE$mean=-RActE$mean
plot(RActE)
plot(Rrain)
Rprecip=Rrain
Rprecip$mean=Rprecip$mean+Rsnow$mean
Rbalance=Rprecip
Rbalance$mean=Rprecip$mean+RActE$mean


Cairo::Cairo(
  20, #length
  15, #width
  file = "D:/tilloal/Documents/01_Projects/RegimeShifts/plots/TicinoWB.jpg",
  type = "png", #tiff
  bg = "transparent", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)

qlim=c(-3,1.2*max(Rprecip$mean))
plot(Rprecip$date, Rprecip$mean, type="n", axes=FALSE,ylim=qlim,xaxs="i",yaxs="i",
     xlab = NA, ylab = expression(paste("Water balance(mm/d)")))
abline(v = format(mois.deb,"%j"), col="lightgrey", lty=3)
polygon(c(Rrain$date,rev(Rrain$date)),c(rep(0,length(Rrain$mean)),rev(Rrain$mean)),
        col=alpha("darkblue",.7),border="transparent")
polygon(c(Rrain$date,rev(Rrain$date)),c(rep(0,length(RActE$mean)),rev(RActE$mean)),
        col=alpha("orange",.7),border="transparent")
polygon(c(Rprecip$date,rev(Rprecip$date)),c(Rprecip$mean,rev(Rrain$mean)),
        col=alpha("grey",.7),border="transparent")
lines(Rbalance$date, Rbalance$mean, col="black",lwd=3)
axis(2)
axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=1)
box()
legend("topleft", leg=c("snowfall","rainfall","Actual Evapotranspiration","P-ET"),
       lwd=c(6,6,6,2), col=c("grey","darkblue","orange","black"),
       cex=1, lty=c(1,1,1,1),bty="n")

dev.off()


snowinfluence=RegimeSM$mean/RegimeR$mean
mean(snowinfluence)



#Gini coefficient
install.packages("ineq")
library(ineq)
x=Rain[[cloc]]
x=RegimeQ$mean
gini<-Gini(x)
gini
dt=4
SaveRF=c(1951:2020)
for (col_name in cnames) {
  # col_name=cnames[100]
  # col_name="290666"
  print(col_name)
  Trun <- preprocess_frac(time=timeStampX, input_var=Runoff_fraction, col_name, dsel,dt)
  Yagg <- process_frac(Trun)
  #plot(Yagg)
  SaveRF=cbind(SaveRF,Yagg$val)
}

SaveRF=data.frame(SaveRF)
colnames(SaveRF)[-1]=cnames

save(SaveRF,file="D:/tilloal/Documents/01_Projects/RegimeShifts/RunoffFraction.Rdata")

ys1=c(1951:1970)
my=which(!is.na(match(SaveRF$SaveRF,ys1)))
Keep=(SaveRF[my,-1])
Keep_RFmeans=colMeans(Keep)


km=match(ppoints$outlets.y,cnames)
ppoints$runoffraction=Keep_RFmeans[km]

palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = TRUE, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=runoffraction),color="black",alpha=1,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="runoff fraction",
    breaks=seq(0,1,.05), limits=c(0,0.2)) +
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

ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/runoffraction.png", width=20, height=15, units=c("cm"),dpi=1500)






#PART 2 Simple experimet with 2 key metrics----

#Save the green and blue water from here as well as normal Q