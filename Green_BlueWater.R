
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



#function returning the mean and seasonality index
Q1_mean_s <- function(data, yr) {
  
  
  yrwindow <- yr
  subs <- which(!is.na(match(year(as.Date(data[,1])), yrwindow)))
  data<-data[subs,]
  
  
  # Compute Qs
  data$Qs <- data[,2]
  #specific discharge
  # data$Qs <- data[,2]/ (cupa * 1000000) * 1000 * 3600 * 24
  
  
  #mean
  qm <- mean(data[,2])
  
  # # Apply running mean again
  # data[,2]=tsEvaNanRunningMean(data[,2],30)
  # 
  # Regi <- RegimeFast(data)
  # 
  # #Seasonality index
  # SI<-sum(abs(Regi[,2]-qm))/(qm*365)
  # Save parameters
  # params <- c(mean=qm,SI=SI)
  
  return(qm)
}

Q1_sum_s <- function(data, yr) {
  
  
  yrwindow <- yr
  subs <- which(!is.na(match(year(as.Date(data[,1])), yrwindow)))
  data<-data[subs,]
  
  
  # Compute Qs
  data$Qs <- data[,2]
  #specific discharge
  # data$Qs <- data[,2]/ (cupa * 1000000) * 1000 * 3600 * 24
  
  
  #mean
  #plot(data$Qs)
  qm <- sum(data[,2])
  
  # # Apply running mean again
  # data[,2]=tsEvaNanRunningMean(data[,2],30)
  # 
  # Regi <- RegimeFast(data)
  # 
  # #Seasonality index
  # SI<-sum(abs(Regi[,2]-qm))/(qm*365)
  # Save parameters
  # params <- c(mean=qm,SI=SI)
  
  return(qm)
}

#Step 1> prepocess

#I need to create the green and blue water from components

ActEvapo<- fread(paste0(hydroDir,"/tss/HERA_Histo/ActEvapo_1951_2020.csv"),header=TRUE)
time=ActEvapo$V1
timeStampX=time[order(time)]
matcol=match(UpArea$outlets,as.numeric(colnames(ActEvapo)))
cnames=names(ActEvapo)[matcol]
rm(ActEvapo)
gc()

#upload from Q_components file
Qv="greenwater"
#Qe is the sum of Q components

if (Qv=="bluewater"){
  Qe=fread(file="D:/tilloal/Documents/01_Projects/RegimeShifts/data/BlueWater.csv",header=TRUE)
}
#Qe can also be the soil balance
if (Qv=="greenwater"){
  Qe=fread(file="D:/tilloal/Documents/01_Projects/RegimeShifts/data/GreenWater.csv",header=TRUE)
}


RegimeDev=function(data){
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
                   q05=sapply(ind.j, function(x) quantile(data$Q[x],.05, na.rm=TRUE)),
                   q95=sapply(ind.j, function(x) quantile(data$Q[x],.95, na.rm=TRUE)))
  return(Qc)
}

yearstart=1951
wsize=30
dvar=Qe
RegimeFunctionDev <- function(dvar, cloc, yearstart, wsize,timeStampX,cupa=NA,Q=FALSE) {
  # Convert data into a data frame with timeStampX
  
  #ly <- length(years)
  yr=yearstart
  yrwindow <- yr:(yr + wsize)
  subs <- which(!is.na(match(year(as.Date(timeStampX)), yrwindow)))
  
  data<-dvar[[cloc]][subs]
  data <- data.frame(timeStampX[subs], data)
  # Create a new column 'Qs' by multiplying the second column by 4*nbdays
  dt1 = min(diff(timeStampX), na.rm = T)
  dt = as.numeric(dt1)
  tdim = attributes(dt1)$units
  if (tdim == "hours") 
    dt = dt/24
  if (tdim == "seconds") 
    dt = dt/3600
  dt1=1/dt
  
  if (Q==TRUE){
    data[, 2]=data[, 2]/(4*cupa*1000000)*1000*3600*24
    # data$Qs <- tsEvaNanRunningMean(data[, 2], dt)
  }
  dt2=30
  dt=dt1*dt2
  data$Qs <- RunningSum(data[, 2], dt)
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
  Regime <- RegimeDev(data)
  Rmean <- RegimeFast(data)
  
  # Return the result
  return(Regime)
}







dsel <- hour(timeStampX)
dt=4
yrlist=c(1951:2020)
Regysave=c()
Regy_d_save=c()
Regy_w_save=c()
resave=c()
metrics=c()

###MODIFICATIONS for soil moisture balance
for (cloc in cnames){
  #cloc="311057"
  print(cloc)
  #Qep=preprocess_in(timeStampX,Qe,cloc,dsel,dt)
  #plot(Qep)
  #Qy1=Q1_sum_s(Qep,2020)
  data=preprocess_in(timeStampX,Qe,cloc,dsel,dt)
  #Rqyi<-RegimeFunctionM(dvar=Qe, cloc, yearstart=1951, wsize=1, timeStampX)
  #quantiles 5% and 95% of regimes
  Rqyi <- RegimeFunctionDev(dvar = Qe,
                          cloc = cloc,
                          yearstart = 1951,
                          wsize = 30 ,
                          timeStampX = timeStampX)
  Regyt=Rqyi$date
  center_year1=data.frame(day=Rqyi$date)
  n_years   <- length(yrlist)
  # `res` will hold mean and SI for each year → 2 columns
  res_mat   <- matrix(NA_real_, nrow = n_years, ncol = 4,
                      dimnames = list(NULL, c("mean", "SI","wet-dev","dry-dev")))
  # `Regy_means` will hold the column of `Rqy1$mean` for each year
  # (plus the original date column that you already have)
  Regy_means <- matrix(NA_real_, nrow = length(Regyt), ncol = n_years)
  
  Regy_wet <- matrix(NA_real_, nrow = length(Regyt), ncol = n_years)
  Regy_dry <- matrix(NA_real_, nrow = length(Regyt), ncol = n_years)
  # Optional: keep the years for later reference
  year_vec  <- integer(n_years)
  res=c()
  pb <- txtProgressBar(min = 0, max = n_years, style = 3)   # lightweight progress bar
  
  for (i in seq_len(n_years)) {
    year <- yrlist[i]
    year_vec[i] <- year
    
    # ---- Year‑specific calculations ---------------------------------
    Rqy1 <- RegimeFunctionM(dvar = Qe,
                            cloc = cloc,
                            yearstart = year,
                            wsize = 0,
                            timeStampX = timeStampX)
    
    yDev1=Rqy1$mean-Rqyi$q95
    yDev2=Rqy1$mean-Rqyi$q05
    mw=which(yDev1>0)
    wdev=length(mw)
    md=which(yDev2<0)
    ddev=length(md)
    
    # mean for this year (using the pre‑computed `data_const`)
    if (Qv=="bluewater"){
      mean_year <- Q1_mean_s(data, year)
    }
    
    if (Qv=="greenwater"){
      #for soil sum and not mean is used
      mean_year <- Q1_sum_s(data, year)
      #for soil, bring all values to positive by adding minim value
      Rqy1$mean=Rqy1$mean-min(Rqy1$mean)
      m1=sum(abs(Rqy1$mean))
    }
    Rqy1$wdev=0
    Rqy1$ddev=0
    Rqy1$wdev[mw]=1
    Rqy1$ddev[md]=1
    #for Soil m1=0
    #m1=0
    # Seasonality index (SI) – note the use of the pre‑computed `m1_const`
    #SI <- sum(abs(Rqy1[, 2] - m1/12)) / (m1)
    
    SI <- sum(abs(Rqy1[, 2] - m1/12))/ m1
    
    # if (m1>0) SI=SI
    # if (m1<0) SI=-SI
    #print(SI)
    # ---- Store results ------------------------------------------------
    res_mat[i, ] <- c(mean_year, SI, wdev, ddev)  # two‑column result matrix
    Regy_means[, i] <- Rqy1$mean              # column of means for this year
    Regy_wet[, i] <- Rqy1$wdev   
    Regy_dry[, i] <- Rqy1$ddev  
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # ---- Assemble final objects -----------------------------------------
  # 1) `res` as a data.frame (or tibble) with a `year` column
  res <- data.frame(year = year_vec, res_mat)
  #plot(res$mean)
  # 2) `Regy` – original dates + one column per year of `Rqy1$mean`
  Regy <- cbind(mday = Regyt, Regy_means)
  Regy=data.frame(Regy)
  
  Regy_w <- cbind(mday = Regyt, Regy_wet)
  Regy_w=data.frame(Regy_w)
  
  Regy_d <- cbind(mday = Regyt, Regy_dry)
  Regy_d=data.frame(Regy_d)
  # colnames(Regy)[-1] <- paste0("mean_", yrlist)   # nice column names
  # 
  # for (year in yrlist){
  #   print(year)
  #   Rqy1<-RegimeFunctionMF(dvar=Qe, cloc, yearstart=year, wsize=1, timeStampX, center_year1)
  #   mean=Q1_mean_s(data,year)
  #   m1=mean(Rqy1$mean)
  #   #Seasonality index
  #   SI<-sum(abs(Rqy1[,2]-m1))/(m1*12) 
  #   results<-c(mean,SI)
  #   res<-rbind(res,results)
  #   # Rqy1$year=rep(year,12)
  #   Regy=cbind(Regy,Rqy1$mean)
  # }
  # Regy=data.frame(Regy)
  colnames(Regy)[c(2:71)]=yrlist
  colnames(Regy_d)[c(2:71)]=colnames(Regy_w)[c(2:71)]=yrlist
  Regy$cloc=rep(cloc,length(Regy[,1]))
  Regy_d$cloc=rep(cloc,length(Regy[,1]))
  Regy_w$cloc=rep(cloc,length(Regy[,1]))
  res=data.frame(res)
  o=ecdf(res$mean)
  perc=o(res$mean)
  res$dm=c(0,abs(diff(perc)))
  
  resp1=res[c(1:30),]
  plot=F
  if (plot==T){ 
    n=4
    mean(resp1[,n])
    qh1=quantile(resp1[,n],.75)
    ql1=quantile(resp1[,n],.25)
    qh2=quantile(resp1[,n],.90)
    ql2=quantile(resp1[,n],.10)
    plot(res[,n],type="o",lwd=2,pch=19)
    abline(h=mean(resp1[,n]),lwd=2)
    abline(h=qh1,lty=2,col="blue",lwd=2)
    abline(h=qh2,lty=2,col="darkblue",lwd=2)
    abline(h=ql1,lty=2,col="orange",lwd=2)
    abline(h=ql2,lty=2,col="red",lwd=2)
    abline(v=30,lwd=2,col="grey")
    abline(v=41,lwd=2,col="grey")
    
    
    mean(resp1[,3])
    qh1=quantile(resp1[,3],.75)
    ql1=quantile(resp1[,3],.25)
    qh2=quantile(resp1[,3],.90)
    ql2=quantile(resp1[,3],.10)
    plot(res[,3],type="o",lwd=2,pch=19)
    abline(h=mean(resp1[,3]),lwd=2)
    abline(h=qh1,lty=2,col="blue",lwd=2)
    abline(h=qh2,lty=2,col="darkblue",lwd=2)
    abline(h=ql1,lty=2,col="orange",lwd=2)
    abline(h=ql2,lty=2,col="red",lwd=2)
    abline(v=30,lwd=2,col="grey")
    abline(v=41,lwd=2,col="grey")
  }
  
  #compare periods:
  
  resp2<-res[c(41:70),]
  
  
  #kge style
  
  #availability
  b=mean(resp2$mean)/mean(resp1$mean)
  #variability
  v=mean(resp2$dm)/mean(resp1$dm)
  #seasonality
  s=mean(resp2$SI)/mean(resp1$SI)
  metric=1-((b-1)^2+(v-1)^2+(s-1)^2)
  
  colnames(res)=c("Qmean","SI","Variability")
  res$year=yrlist
  res$loc=rep(cloc,length(res$Qmean))
  metrix=data.frame(b,v,s,metric,cloc)
  
  Regysave=rbind(Regysave,Regy)
  Regy_d_save=rbind(Regy_d_save,Regy_d)
  Regy_w_save=rbind(Regy_w_save,Regy_w)
  metrics=rbind(metrics,metrix)
  resave<-rbind(resave,res)
  
}





mean(metrics$v)
names(resave)=c("Year","Avail","SI","wet-dev","dry-dev","Variability","Year2","loc")

#save the stuff
Savelist=list(Regimes=Regysave,Regimes_dry=Regy_d_save,Regimes_wet=Regy_w_save,yearly_res=resave,metrics=metrics)

if (Qv=="bluewater"){
  save(Savelist,file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_blue1.Rdata")
}
if (Qv=="greenwater"){
  save(Savelist,file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_green1.Rdata")
}



## load Biogeographic regions ----
biogeo <- read_sf(dsn = paste0(hydroDir,"/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=right_join(biogeomatch,UpArea, by="latlong")


#starting from outputs

if (Qv=="bluewater"){
 load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_blue1.Rdata")
}
if (Qv=="greenwater"){
 load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_green1.Rdata")
}

resave=Savelist$yearly_res

studyvar=c("Avail","SI")
for (st in studyvar){
  st=studyvar[1]
  mname=match(st,names(resave))
  names(resave)[2]="Q"
  names(resave)[3]="Avail"
}
names(resave)[2]="Q"
names(resave)[3]="Avail"
#match resave with upstream area
mup=match(as.numeric(resave$loc),UpArea$outlets)
resave$upa=UpArea$upa[mup]

regios=match(as.numeric(resave$loc),biogeo_rivers$outlets)
resave$region=biogeo_rivers$code[regios]

resave=resave[which(resave$region=="Atlantic"),]

#aggregate by year
Qmean_yAg<- aggregate(resave$Avail, 
                      by = list(year = resave$Year),
                      function(x) c(mean= mean(x,na.rm=T),q1=quantile(x,0.1,na.rm=T),q2=quantile(x,0.9,na.rm=T)))
Qmean_yAg<-do.call(data.frame, Qmean_yAg)

plot(Qmean_yAg$year,Qmean_yAg[,3],type="o")
QMAv=tsEvaNanRunningMean(Qmean_yAg[,3],10)
lines(Qmean_yAg$year,QMAv, col=2, lwd=2)
abline(h=mean(Qmean_yAg[,3]))


SI_yAg<- aggregate(resave$SI, 
                   by = list(year = resave$Year),
                   function(x) c(mean= mean(x,na.rm=T),q1=quantile(x,0.1,na.rm=T),q2=quantile(x,0.9,na.rm=T)))
SI_yAg<-do.call(data.frame,SI_yAg)


V_yAg<- aggregate(resave$Variability, 
                  by = list(year = resave$Year),
                  function(x) c(mean= mean(x,na.rm=T),q1=quantile(x,0.1,na.rm=T),q2=quantile(x,0.9,na.rm=T)))
V_yAg<-do.call(data.frame,V_yAg)

plot(SI_yAg$year,SI_yAg$x.mean,type="o")
head(resave)


cloc="290666"
River="Ticino"

rc1=resave[which(resave$loc==cloc),]
#do a map?
mean(rc1$`dry-dev`/12)*100
plot(rc1$Avail, type="l")
#30 years averages


#Catchment stuff

CatUpA=inner_join(Catf7,UpArea,by=c("llcoord"="latlong"))
matcat=match(outhybas$latlong,CatUpA$llcoord)
CatUpA=CatUpA[matcat,]
ratUp=CatUpA$SUB_AREA/CatUpA$upa

plot(ratUp, log="y")
abline(h=1.5)
abline(h=0.5)

CatUpA$outlet=1
CatUpA$outlet[which(ratUp>1.2)]=3
CatUpA$outlet[which(ratUp<0.8)]=3
CatUpA$outlet[which(CatUpA$outlet==3 & CatUpA$upa<1e4)]=2

CatLR=CatUpA[which(CatUpA$outlet==3),]
st_geometry(CatLR)=NULL
CatO=CatUpA[-which(CatUpA$outlet==3),]




#map of each metric at the beginning of the period
### --------------------------------------------------------------
##  ONE‑STOP CODE: quantiles, centering, binning, counting & plot
## --------------------------------------------------------------



# 1. Input data ---------------------------------------------------

resaveS=resave
#resave=resaveS
resave=resave[-which(resave_med$Year==1951),]
length(unique(CatO$outlets.y))
length(unique(resaveS$loc))
BigR<-which(!is.na(match(as.numeric(resave$loc),CatO$outlets.y)))
resave=resave[BigR,]
# 2. Define the 30‑year baseline ---------------------------------
baseline_start <- min(resave$Year, na.rm = TRUE)         # earliest year in the series
baseline_years <- baseline_start:(baseline_start + 29)   # first 30 calendar years

studyVar="Blue water availability"
# 3. Quantiles on the baseline (per location) --------------------
loc_q <- resave %>%
  filter(Year %in% baseline_years) %>%                  # reference period only
  group_by(loc) %>%
  summarise(
    med = quantile(Avail, 0.50, na.rm = TRUE),           # 50 % → will become 0
    q05 = quantile(Avail, 0.05, na.rm = TRUE),
    q25 = quantile(Avail, 0.25, na.rm = TRUE),
    q75 = quantile(Avail, 0.75, na.rm = TRUE),
    q95 = quantile(Avail, 0.95, na.rm = TRUE),
    .groups = "drop"
  )


# --------------------------------------------------------------
# 3. Quantiles on the baseline (per location) – weighted by upa
# --------------------------------------------------------------
# loc_q <- resave %>%
#   filter(Year %in% baseline_years) %>%          # keep only the 30‑yr reference period
#   group_by(loc) %>%                             # one set of quantiles per location
#   summarise(
#     # ---- Hmisc version (most readable) -----------------------
#     med = wtd.quantile(SI,          weights = upa, probs = 0.50, na.rm = TRUE),
#     q10 = wtd.quantile(SI,          weights = upa, probs = 0.10, na.rm = TRUE),
#     q25 = wtd.quantile(SI,          weights = upa, probs = 0.25, na.rm = TRUE),
#     q75 = wtd.quantile(SI,          weights = upa, probs = 0.75, na.rm = TRUE),
#     q90 = wtd.quantile(SI,          weights = upa, probs = 0.90, na.rm = TRUE),
#     
#     # ---- (optional) matrixStats version -----------------------
#     # med = weightedQuantile(SI, w = upa, probs = 0.50, na.rm = TRUE),
#     # q10 = weightedQuantile(SI, w = upa, probs = 0.10, na.rm = TRUE),
#     # q25 = weightedQuantile(SI, w = upa, probs = 0.25, na.rm = TRUE),
#     # q75 = weightedQuantile(SI, w = upa, probs = 0.75, na.rm = TRUE),
#     # q90 = weightedQuantile(SI, w = upa, probs = 0.90, na.rm = TRUE),
#     
#     .groups = "drop"
#   )
# 
# # Look at the first few rows of the weighted‑quantile table
# head(loc_q)

# 4. Attach thresholds to the full data set ----------------------
df <- resave %>%
  left_join(loc_q, by = "loc")

# 5. Centre Qmean on its median (median → 0) --------------------
df <- df %>%
  mutate(Q_centered = Avail - med)               # negative = below median

# 6. Express the other quantiles as deviations from the median ---
df <- df %>%
  mutate(
    d05 = q05 - med,   # negative
    d25 = q25 - med,   # negative
    d75 = q75 - med,   # positive
    d95 = q95 - med    # positive
  )

# 7. Bin each observation into the six categories ---------------
df <- df %>%
  mutate(
    bin = case_when(
      Q_centered > d95                     ~ "above_95",   # most positive
      Q_centered > d75 & Q_centered <= d95 ~ "75_95",
      Q_centered > 0   & Q_centered <= d75 ~ "med_75",
      Q_centered > d25 & Q_centered <= 0   ~ "25_med",
      Q_centered > d05 & Q_centered <= d25 ~ "05_25",
      Q_centered <= d05                    ~ "below_05",
      TRUE                                 ~ NA_character_
    )
  )

# # 8. Percentage of locations per bin for each year --------------
# perc_by_year <- df %>%
#   group_by(Year, bin) %>%
#   summarise(n_loc = n(), .groups = "drop") %>%               # raw counts
#   complete(Year, bin, fill = list(n_loc = 0)) %>%            # make sure every bin exists
#   group_by(Year) %>%
#   mutate(
#     total_loc = sum(n_loc),
#     perc      = 100 * n_loc / total_loc                     # % of locations
#   ) %>%
#   ungroup() %>%
#   select(Year, bin, perc)


#frequency of local deviation per catchment and per year (averaged over 10 year)
resave$dd10=tsEvaNanRunningMean(resave$Avail,10)
qp=unique(df$med[which(df$loc==cloc)])
plot(resave$Avail[which(resave$loc==cloc)], type="l")
abline(h=qp)
lines(resave$dd10[which(resave$loc==cloc)],col=2,lwd=3)

head (df)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)

# 1. Fit the models and tidy the output
# Make sure to install and load the library

library(segmented)
lin_mod <- lm(Avail ~ Year, data = df)
seg_mod <- segmented(lin_mod, seg.Z = ~Year)
summary(seg_mod) # Tells you the "Estimated Breakpoint"

library(quantreg)
fit <- rq(Avail ~ Year, tau = 0.9, data = df) # Trend of the 90th percentile

library(quantreg)

slope_results_quant <- df %>%
  group_by(loc) %>%
  nest() %>%
  mutate(
    # Fit Quantile Regression (tau = 0.5 is the median)
    model = map(data, ~ rq(Avail ~ Year, tau = 0.5, data = .x)),
    
    # Use broom::tidy to get the slope (estimate for 'Year')
    tidied = map(model, broom::tidy),
    
    # Calculate Mean for the Relative Slope conversion
    mean_avail = map_dbl(data, ~ mean(.x$Avail, na.rm = TRUE))
  ) %>%
  unnest(tidied) %>%
  filter(term == "Year") %>%
  mutate(
    abs_slope = estimate,
    slope_pct = (abs_slope / mean_avail) * 100
  ) %>%
  dplyr::select(loc, abs_slope, slope_pct)


# install.packages("trend")
library(trend)
length(unique(df$loc))
slope_results <- df %>%
  group_by(loc) %>%
  nest() %>%
  mutate(
    # 1. Calculate the Mean of Avail for each group
    mean_avail = map_dbl(data, ~ mean(.x$Avail, na.rm = TRUE)),
    
    # 2. Compute the Sen's Slope
    model = map(data, ~ sens.slope(.x$Avail)),
    abs_slope = map_dbl(model, ~ as.numeric(.x$estimates)),
    
    # 3. Convert to percentage of the mean
    slope_pct = (abs_slope / mean_avail) * 100,
    
    # 4. Extract p-value
    p.value = map_dbl(model, ~ as.numeric(.x$p.value))
  ) %>%
  dplyr::select(loc, mean_avail, abs_slope, slope_pct, p.value)
# 4. View results


ppl <- st_as_sf(CatLR, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)

ppc <- st_as_sf(CatO, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppc <- st_transform(ppc, crs = 3035)

# Compute centroids
CatO_centroids <- st_centroid(ppc)
alle=match(colnames(ppl),colnames(CatO_centroids))
ppoints=CatO_centroids


print(slope_results)

plot(slope_results$slope_pct)
slope_results$outlets=as.numeric(slope_results$loc)

km=match(ppoints$outlets.y,slope_results$outlets)
ppoints$slope=slope_results$slope_pct[km]

# map_df=inner_join(UpArea,slope_results,by=c("outlets"))


#do this for all catchments
# ppl <- st_as_sf(map_df, coords = c("Var1.x", "Var2.x"), crs = 4326)
# ppl <- st_transform(ppl, crs = 3035)



## Snow fraction ----
palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = FALSE, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=slope),color="black",alpha=1,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 6), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="water balance (%)",
    limits=c(-1,1)) +
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

#keep digging 


# ----------------------------------------------------------------
# 8. Sum of upa per year / bin ----------------------------------
# ----------------------------------------------------------------
area_pct_by_year <- df %>%
  group_by(Year, bin) %>%
  summarise(
    upa_sum = sum(upa, na.rm = TRUE),        # total area in this bin & year
    .groups = "drop"
  ) %>%
  # make sure every year has a row for every bin (missing combos → 0)
  complete(Year, bin, fill = list(upa_sum = 0)) %>%
  group_by(Year) %>%
  mutate(
    total_upa = sum(upa_sum),                # total area for the year
    pct       = 100 * upa_sum / total_upa    # **percentage of total area**
  ) %>%
  ungroup() %>%
  dplyr::select(Year, bin, pct)


extreme_locs<-area_pct_by_year[which(area_pct_by_year$bin=="above_95" |area_pct_by_year$bin=="below_05"),]

extreme_y<- stats::aggregate(list(frac=extreme_locs$pct), by = list(year = extreme_locs$Year), 
                             FUN = function(x) c(sum=sum(x, na.rm = T)))
extreme_y <- do.call(data.frame, extreme_y)

plot(extreme_y, type="l")
# 9. Give a negative sign to the three “below‑median” bins ----------
# plot_df <- perc_by_year %>%
#   mutate(
#     signed_perc = ifelse(bin %in% c("below_10", "10_25", "25_med"),
#                          -perc,        # negative side
#                          perc)        # positive side
#   )

plot_df <- area_pct_by_year %>%
  mutate(
    signed_upa = ifelse(bin %in% c("below_05", "05_25", "25_med"),
                        -pct,                # negative side
                        pct)
  )# positive side


# ----------------------------------------------------------------
# 8. 30‑yr centred moving‑average for **each** bin -------------
# ----------------------------------------------------------------
trend_df <- plot_df %>%
  arrange(Year) %>%
  group_by(bin) %>%
  mutate(
    trend = rollapply(signed_upa,
                      width = 30,               # 30‑year window
                      FUN = mean,
                      align = "right",
                      fill = NA,
                      na.rm = TRUE)
  ) %>%
  ungroup()

# ----------------------------------------------------------------
# 8. Fit a *linear* trend (OLS) for each bin over the whole series
# ----------------------------------------------------------------
trend_lm <- plot_df %>%
  group_by(bin) %>%                     # one model per bin
  nest() %>%                            # create a list‑column called `data`
  mutate(
    model  = map(data, ~ lm(signed_upa ~ Year, data = .x)),   # fit the lm
    # get fitted values for every year belonging to this bin
    fitted = map2(model, bin, ~ {
      yrs   <- area_pct_by_year$Year[area_pct_by_year$bin == .y]
      # `augment()` adds a column `.fitted` with the predictions
      broom::augment(.x, newdata = tibble(Year = yrs))
    })
  ) %>%
  dplyr::select(fitted) %>%                    # keep only the fitted data
  unnest(fitted) %>%                    # back to a tidy data frame
  dplyr::rename(trend = .fitted)               # rename the column with the predicted line

trend_stats <- plot_df %>%
  group_by(bin) %>%
  do(broom::tidy(lm(signed_upa ~ Year, data = .))) %>%   # one tidy table per bin
  filter(term == "Year") %>%                             # keep only the slope row
  mutate(
    signif = ifelse(p.value < 0.05, "significant", "ns")
  ) %>%
  dplyr::select(bin, estimate, std.error, p.value, signif)

trend_fitted <- plot_df %>%
  left_join(trend_stats, by = "bin") %>%          # bring the slope‑info onto the data
  group_by(bin) %>%
  mutate(trend = predict(lm(signed_upa ~ Year, data = cur_data()))) %>%
  ungroup()


# ----------------------------------------------------------------
# 9. Combine the two pairs of bins
#    (below_10 + 10_25)   &   (above_90 + 75_90)
# ----------------------------------------------------------------
# trend_wide <- trend_fitted %>%
#   select(Year, bin, trend, signif) %>%
#   pivot_wider(names_from = bin, values_from = trend)

trend_wide <- trend_fitted %>%
  # bring *both* the numeric trend and the flag into the wide table
  dplyr::select(Year, bin, trend, signif) %>%
  pivot_wider(
    names_from = bin,
    values_from = c(trend, signif),
    names_sep = "_"          # creates columns like trend_25_med, signif_25_med
  )

# trend_plot1 <- trend_wide %>%
#   # keep only the columns we need for plotting
#   select(Year, `25_med`, med_75) %>%
#   pivot_longer(cols = c(`25_med`, med_75),
#                names_to = "bin",
#                values_to = "trend")
#  
trend_plot1 <- trend_wide %>%
  # keep only the two middle bins and their significance flags
 dplyr::select(Year,
         trend_25_med, signif_25_med,
         trend_med_75, signif_med_75) %>%
  pivot_longer(
    cols = c(trend_25_med, trend_med_75),
    names_to = "bin",
    values_to = "trend",
    names_prefix = "trend_"
  ) %>%
  # bring the matching significance flag back
  left_join(
    trend_wide %>%
      dplyr::select(Year,
             signif_25_med, signif_med_75) %>%
      pivot_longer(
        cols = c(signif_25_med, signif_med_75),
        names_to = "bin_sig",
        values_to = "signif",
        names_prefix = "signif_"
      ),
    by = c("Year", "bin" = "bin_sig")
  )


# trend_combined <- trend_wide %>%
#   mutate(
#     # combine the lower pair
#     below_comb = `25_med` + `10_25`,
#     # combine the upper pair
#     above_comb = med_75 + `75_90`
#   ) %>%
#   # keep only the columns we need for plotting
#   select(Year, below_comb, above_comb, signif) %>%
#   pivot_longer(cols = c(below_comb, above_comb),
#                names_to = "combo_bin",
#                values_to = "combo_trend") 

trend_combined <- trend_wide %>%
  mutate(
    # ----- numeric combined trends ---------------------------------
    below_comb = trend_25_med + trend_05_25,
    above_comb = trend_med_75 + trend_75_95,
    
    # ----- significance of the combined trend ----------------------
    # “significant” if *either* component is significant
    signif_below_comb = ifelse(signif_25_med == "significant" |
                                 signif_05_25 == "significant",
                               "significant", "ns"),
    signif_above_comb = ifelse(signif_med_75 == "significant" |
                                 signif_75_95 == "significant",
                               "significant", "ns")
  ) %>%
  # keep only the columns we need for the plot
  dplyr::select(Year,
         below_comb, signif_below_comb,
         above_comb, signif_above_comb) %>%
  # ----- back to long format (trend values only) -----------------
pivot_longer(
  cols = c(below_comb, above_comb),
  names_to = "combo_bin",
  values_to = "trend"
) %>%
  # ----- add the matching significance flag ----------------------
mutate(
  signif = ifelse(combo_bin == "below_comb",
                  signif_below_comb,
                  signif_above_comb)
) %>%
  # drop the helper columns that are no longer needed
  dplyr::select(Year, combo_bin, trend, signif)


#10y moving average ----
library(dplyr)
library(tidyr)
library(Rtseva)

# ----------------------------------------------------------------
# 8. Compute 10-year Running Mean (instead of Linear Trend)
# ----------------------------------------------------------------
trend_fitted <- plot_df %>%
  arrange(bin, Year) %>%            # Crucial for running means
  group_by(bin) %>%
  mutate(
    # Use the Rtseva function for the 10y window
    trend = tsEvaNanRunningMean(signed_upa, windowSize = 10),
    # Since running means don't have a single p-value, 
    # we mark them as "MA" (Moving Average) or keep a placeholder
    signif = "MA" 
  ) %>%
  ungroup()

# ----------------------------------------------------------------
# 9. Pivot and Combine Bins
# ----------------------------------------------------------------
# We pivot to wide format to easily sum the "combined" categories
trend_wide <- trend_fitted %>%
  dplyr::select(Year, bin, trend) %>%
  pivot_wider(
    names_from = bin,
    values_from = trend,
    names_prefix = "trend_"
  )

# --- 9a. Middle Bins Plotting Data ---
trend_plot1 <- trend_wide %>%
  dplyr::select(Year, trend_25_med, trend_med_75) %>%
  pivot_longer(
    cols = starts_with("trend_"),
    names_to = "bin",
    values_to = "trend",
    names_prefix = "trend_"
  ) %>%
  mutate(signif = "MA") # Placeholder for the label

# --- 9b. Combined Bins (Aggregated Extremes) ---
trend_combined <- trend_wide %>%
  mutate(
    # Combine the numeric running means
    below_comb = trend_25_med + trend_05_25 ,
    above_comb = trend_med_75 + trend_75_95 ,
    below_comb2 = trend_25_med + trend_05_25 + trend_below_05,
    above_comb2 = trend_med_75 + trend_75_95 + trend_above_95
  ) %>%
  dplyr::select(Year, below_comb, above_comb, below_comb2,above_comb2) %>%
  pivot_longer(
    cols = c(below_comb, above_comb, below_comb2, above_comb2),
    names_to = "combo_bin",
    values_to = "trend"
  ) %>%
  mutate(signif = "MA")



# ----------------------------------------------------------------
# 10. **Re‑order the factor for the *stacking* you want**
# ----------------------------------------------------------------
##   <- change the vector below to whatever stack order you need.
##      The first element will be the *bottom* slice.
plot_df$bin <- factor(plot_df$bin,
                      levels = c(
                        "below_05",   # 1️⃣ bottom
                        "05_25",      # 2️⃣
                        "25_med", 
                        "above_95", # 3️⃣
                        "75_95",      # 5️⃣
                        "med_75"   # 4️⃣
                      ))





# 10. RED‑to‑BLUE colour palette ------------------------------------
# A smooth gradient from deep red (most negative) → orange → yellow → light‑blue → deep blue (most positive)
pal_fun <- colorRampPalette(c("#D73027", "#F46D43", "#FEE090",
                              "#FEE090", "#4575B4", "#313695"))
my_palette <- setNames(pal_fun(6),
                       c("below_05", "05_25", "25_med",
                         "med_75", "75_95", "above_95"))

my_palette2 <- c( "above_comb"="lightblue", "above_comb2"= "darkblue","below_comb"="gold", "below_comb2"="darkred")

# 11. Plot – stacked bar using geom_bar(stat = "identity") ------
ggplot(plot_df,
       aes(x = Year,           # discrete x‑axis
           y = signed_upa,            # height (negative = below median)
           fill = bin)) +
  geom_bar(stat = "identity",        # classic stacked‑bar call
           width = 0.8,
           position = "stack",alpha=.8) +
  # Force the legend to follow the same order as the factor levels
  scale_fill_manual(values = my_palette,
                    breaks = c("below_05", "05_25", "25_med",
                               "med_75", "75_95", "above_95"),
                    name = "Category") +
  labs(
    title    = "Distribution of area under different flow conditions",
    subtitle = paste0(studyVar),
    x        = "Year",
    y        = "Distribution of area (%)"
  ) +
  # geom_line(data = plot_df,
  #           aes(x = as.numeric(Year), y = signed_upa,
  #               colour = bin),
  #           size = 1.2, inherit.aes = FALSE, position="stack")+

  geom_line(data = trend_combined,
            aes(x = as.numeric(Year), y = trend,
                colour = combo_bin),
            size = 1.2, inherit.aes = FALSE) +
  # geom_line(data = trend_plot1,
  #           aes(x = as.numeric(Year), y = trend,colour=bin),
  #           size = 1.2, inherit.aes = FALSE) +
  scale_color_manual(name= "10y moving averages", values = my_palette2, 
                     labels=c("0-05","05-25","75-95","95-100")) +
  # ---- linetype legend (optional) -----------------------------------
scale_linetype_manual(values = c(significant = "solid", ns = "dashed"),
                      name = "Trend significance") +
  # geom_hline(yintercept = 0, colour = "grey30", linetype = "dashed") +
  # geom_hline(yintercept = 25, colour = "grey30", linetype = "dashed") +
  # geom_hline(yintercept = -25, colour = "grey30", linetype = "dashed") +
  scale_y_continuous(
    breaks = c(-50, -25, 0, 25, 50),   
    expand = expansion(mult = c(0, .05)),
    labels = function(x) paste0(abs(x), "%"),
    #breaks = seq(-100,100,20)# show absolute % on axis
  ) +
  scale_x_continuous( breaks= seq(1950,2020,5) # show absolute % on axis
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.major.y =  element_line(color = "black", size = 0.5, linetype = "solid"),
    panel.grid.major.x =  element_line(color = "grey", size = 0.5, linetype = "solid"),
    panel.grid.minor   = element_blank()
  )



library(dplyr)
library(tidyr)
library(purrr)
library(rsample)
library(dplyr)

# 1. Define the quantile function (keep it simple)
calc_qs <- function(x) {
  # Sample with replacement
  resample <- sample(x, size = length(x), replace = TRUE)
  # Return as a tiny data frame/tibble
  data.frame(
    q05 = quantile(resample, 0.05, na.rm = TRUE),
    q95 = quantile(resample, 0.95, na.rm = TRUE)
  )
}

# 2. Run it across all locations
# We use reframe() to handle the 1000 rows generated per location
n_boots <- 100

boot_results_fast <- resave %>%
  filter(Year %in% baseline_years) %>%
  group_by(loc) %>%
  reframe(
    # replicate() runs our function n times and binds the results
    replicate(n_boots, calc_qs(Avail), simplify = FALSE) %>% bind_rows()
  )

# 3. Now compute your CI
loc_q_conf <- boot_results_fast %>%
  group_by(loc) %>%
  summarise(across(everything(), 
                   list(low = ~quantile(.x, 0.025), 
                        high = ~quantile(.x, 0.975)), 
                   .names = "{.col}_{.fn}"))




# 4. Attach thresholds to the full data set ----------------------
df2 <- resave %>%
  left_join(loc_q_conf, by = "loc")

df2 <- df2 %>%
  left_join(loc_q, by = "loc")

# 7. Bin each observation into the six categories ---------------
df2 <- df2 %>%
  mutate(
    bin0 = case_when(
      Avail > q95                     ~ "above_95",   # most positive
      Avail < q05                     ~ "below_05",
      TRUE                                 ~ NA_character_
    )
  )

df2 <- df2 %>%
  mutate(
    bin1 = case_when(
      Avail > q95_high                     ~ "above_95h",   # most positive
      Avail < q05_high                    ~ "below_05h",
      TRUE                                 ~ NA_character_
    )
  )

df2 <- df2 %>%
  mutate(
    bin2 = case_when(
      Avail > q95_low                     ~ "above_95l",   # most positive
      Avail < q05_low                    ~ "below_05l",
      TRUE                                 ~ NA_character_
    )
  )

area_pct_by_year0 <- df2 %>%
  group_by(Year, bin0) %>%
  summarise(
    upa_sum = sum(upa, na.rm = TRUE),        # total area in this bin & year
    .groups = "drop"
  ) %>%
  # make sure every year has a row for every bin (missing combos → 0)
  complete(Year, bin0, fill = list(upa_sum = 0)) %>%
  group_by(Year) %>%
  mutate(
    total_upa = sum(upa_sum),                # total area for the year
    pct       = 100 * upa_sum / total_upa    # **percentage of total area**
  ) %>%
  ungroup() %>%
  dplyr::select(Year, bin0, pct)


area_pct_by_year1 <- df2 %>%
  group_by(Year, bin1) %>%
  summarise(
    upa_sum = sum(upa, na.rm = TRUE),        # total area in this bin & year
    .groups = "drop"
  ) %>%
  # make sure every year has a row for every bin (missing combos → 0)
  complete(Year, bin1, fill = list(upa_sum = 0)) %>%
  group_by(Year) %>%
  mutate(
    total_upa = sum(upa_sum),                # total area for the year
    pct       = 100 * upa_sum / total_upa    # **percentage of total area**
  ) %>%
  ungroup() %>%
  dplyr::select(Year, bin1, pct)


area_pct_by_year2 <- df2 %>%
  group_by(Year, bin2) %>%
  summarise(
    upa_sum = sum(upa, na.rm = TRUE),        # total area in this bin & year
    .groups = "drop"
  ) %>%
  # make sure every year has a row for every bin (missing combos → 0)
  complete(Year, bin2, fill = list(upa_sum = 0)) %>%
  group_by(Year) %>%
  mutate(
    total_upa = sum(upa_sum),                # total area for the year
    pct       = 100 * upa_sum / total_upa    # **percentage of total area**
  ) %>%
  ungroup() %>%
  dplyr::select(Year, bin2, pct)


plot(area_pct_by_year0$pct[which(area_pct_by_year0$bin0=="above_95")], type="o")
lines(area_pct_by_year1$pct[which(area_pct_by_year1$bin1=="above_95h")], type="l")
lines(area_pct_by_year2$pct[which(area_pct_by_year2$bin2=="above_95l")], type="l")



regimewet=Savelist[[3]]
mupa=(((match(regimewet$cloc,UpArea$outlets))))
regimewet$upa=UpArea$upa[mupa]

#remove large catchments
BigR<-which(!is.na(match(as.numeric(regimewet$cloc),CatO$outlets.y)))
regimewet=regimewet[BigR,]

head(regimewet)

library(dplyr)
library(tidyr)

regime_long <- regimewet %>%
  # 1. Pivot all columns that are years (1951 to 2020)
  #    We use matches() with a regex to grab all 4-digit numeric columns
  pivot_longer(
    cols = matches("^[0-9]{4}$"), 
    names_to = "Year", 
    values_to = "Value"
  ) %>%
  # 2. Convert Year to an integer for better sorting/math
  mutate(Year = as.integer(Year)) %>%
  # 3. Reorder columns for clarity
  dplyr::select(cloc, Year, mday, upa, Value)

# View the result
head(regime_long)

library(dplyr)
library(tidyr)

weighted_regime <- regimewet %>%
  # 1. Pivot to long format to get Year as a variable
  pivot_longer(
    cols = matches("^[0-9]{4}$"), 
    names_to = "Year", 
    values_to = "Value"
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  # 2. Group by the time variables (Year and Day)
  group_by(Year, mday) %>%
  summarise(
    # Total UPA of all locations for this year/day
    total_area = sum(upa, na.rm = TRUE),
    # UPA of only the locations where Value is 1
    active_area = sum(upa[Value == 1], na.rm = TRUE),
    # Calculate weighted percentage
    pct_active = (active_area / total_area) * 100,
    .groups = "drop"
  )

# View the results
head(weighted_regime)

mean(weighted_regime$pct_active)

annual_regime_mean_w <- weighted_regime %>%
  group_by(Year) %>%
  summarise(
    mean_pct_active = mean(pct_active, na.rm = TRUE),
    # Optional: also calculate the max activity seen in a single day that year
    max_pct_active = max(pct_active, na.rm = TRUE),
    .groups = "drop"
  )

# View results
head(annual_regime_mean_w)

annual_regime_mean_w <- annual_regime_mean_w %>%
  arrange(Year) %>%
  mutate(
    trend_10y = tsEvaNanRunningMean(mean_pct_active, 10)
  )

# Quick plot to see the trend
library(ggplot2)

ggplot(annual_regime_mean_w, aes(x = Year, y = mean_pct_active)) +
  geom_line(alpha = 0.3) + # Raw annual values
  geom_line(aes(y = trend_10y), color = "blue", size = 1) + # 10y Moving Average
  theme_minimal() +
  labs(title = "Annual Mean Percentage of Area with Wet Events",
       subtitle = "Blue line: 10-year running mean (Rtseva)",
       y = "Mean % Area Active")



regimedry=Savelist[[2]]
mupa=(((match(regimedry$cloc,UpArea$outlets))))
regimedry$upa=UpArea$upa[mupa]

#remove large catchments
BigR<-which(!is.na(match(as.numeric(regimedry$cloc),CatO$outlets.y)))
regimedry=regimedry[BigR,]

head(regimedry)

library(dplyr)
library(tidyr)

regime_long <- regimedry %>%
  # 1. Pivot all columns that are years (1951 to 2020)
  #    We use matches() with a regex to grab all 4-digit numeric columns
  pivot_longer(
    cols = matches("^[0-9]{4}$"), 
    names_to = "Year", 
    values_to = "Value"
  ) %>%
  # 2. Convert Year to an integer for better sorting/math
  mutate(Year = as.integer(Year)) %>%
  # 3. Reorder columns for clarity
  dplyr::select(cloc, Year, mday, upa, Value)

# View the result
head(regime_long)

library(dplyr)
library(tidyr)

weighted_regime <- regimedry %>%
  # 1. Pivot to long format to get Year as a variable
  pivot_longer(
    cols = matches("^[0-9]{4}$"), 
    names_to = "Year", 
    values_to = "Value"
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  # 2. Group by the time variables (Year and Day)
  group_by(Year, mday) %>%
  summarise(
    # Total UPA of all locations for this year/day
    total_area = sum(upa, na.rm = TRUE),
    # UPA of only the locations where Value is 1
    active_area = sum(upa[Value == 1], na.rm = TRUE),
    # Calculate weighted percentage
    pct_active = (active_area / total_area) * 100,
    .groups = "drop"
  )

# View the results
head(weighted_regime)

mean(weighted_regime$pct_active)

annual_regime_mean_d <- weighted_regime %>%
  group_by(Year) %>%
  summarise(
    mean_pct_active = mean(pct_active, na.rm = TRUE),
    # Optional: also calculate the max activity seen in a single day that year
    max_pct_active = max(pct_active, na.rm = TRUE),
    .groups = "drop"
  )

# View results
head(annual_regime_mean_d)

annual_regime_mean_d <- annual_regime_mean_d %>%
  arrange(Year) %>%
  mutate(
    trend_10y = tsEvaNanRunningMean(mean_pct_active, 10)
  )

# Quick plot to see the trend
library(ggplot2)

ggplot(annual_regime_mean_d, aes(x = Year, y = mean_pct_active)) +
  geom_line(alpha = 0.3) + # Raw annual values
  geom_line(aes(y = trend_10y), color = "red", size = 1) + # 10y Moving Average
  theme_minimal() +
  labs(title = "Annual Mean Percentage of Area with dry Events",
       subtitle = "Blue line: 10-year running mean (Rtseva)",
       y = "Mean % Area Active")

#load the european hydroregions


## load Biogeographic regions ----
biogeo <- read_sf(dsn = paste0(hydroDir,"/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=right_join(biogeomatch,UpArea, by="latlong")

regios=match(as.numeric(regimedry$cloc),biogeo_rivers$outlets)
regimedry$region=biogeo_rivers$code[regios]

regimedry_med=regimedry[which(regimedry$region=="Mediterranean"),]



weighted_regime <- regimedry_med %>%
  # 1. Pivot to long format to get Year as a variable
  pivot_longer(
    cols = matches("^[0-9]{4}$"), 
    names_to = "Year", 
    values_to = "Value"
  ) %>%
  mutate(Year = as.integer(Year)) %>%
  # 2. Group by the time variables (Year and Day)
  group_by(Year, mday) %>%
  summarise(
    # Total UPA of all locations for this year/day
    total_area = sum(upa, na.rm = TRUE),
    # UPA of only the locations where Value is 1
    active_area = sum(upa[Value == 1], na.rm = TRUE),
    # Calculate weighted percentage
    pct_active = (active_area / total_area) * 100,
    .groups = "drop"
  )

# View the results
head(weighted_regime)

mean(weighted_regime$pct_active)

annual_regime_mean_wm <- weighted_regime %>%
  group_by(Year) %>%
  summarise(
    mean_pct_active = mean(pct_active, na.rm = TRUE),
    # Optional: also calculate the max activity seen in a single day that year
    max_pct_active = max(pct_active, na.rm = TRUE),
    .groups = "drop"
  )

# View results
head(annual_regime_mean_wm)

annual_regime_mean_wm <- annual_regime_mean_wm %>%
  arrange(Year) %>%
  mutate(
    trend_10y = tsEvaNanRunningMean(mean_pct_active, 10)
  )

# Quick plot to see the trend
library(ggplot2)

ggplot(annual_regime_mean_wm, aes(x = Year, y = mean_pct_active)) +
  geom_line(alpha = 0.3) + # Raw annual values
  geom_line(aes(y = trend_10y), color = "blue", size = 1) + # 10y Moving Average
  theme_minimal() +
  labs(title = "Annual Mean Percentage of Area with Wet Events",
       subtitle = "Blue line: 10-year running mean (Rtseva)",
       y = "Mean % Area Active")


# Map the results ----
df$outlets=as.numeric(df$loc)
df1=df[which(df$Year==1977),]
map_df=inner_join(UpArea,df1,by=c("outlets"))


#do this for all catchments
ppl <- st_as_sf(map_df, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)


# Compute centroids
# CatO_centroids <- st_centroid(ppc)
# alle=match(colnames(ppl),colnames(CatO_centroids))
# ppoints=rbind(ppl,CatO_centroids[,alle])
# 
# km=match(ppoints$outlets.y,cnames)
# ppoints$snowfraction=Keep_means[km]

## Snow fraction ----
palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = FALSE, fixup = TRUE))
ggplot(basemap) +
  geom_sf(fill="white", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,size=upa.x,fill=Avail),color="black",alpha=1,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="water balance (%)",
     limits=c(-100,100)) +
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

#ggsave("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/snowfraction.png", width=20, height=15, units=c("cm"),dpi=1500)



#trend for every location?
#10 year moving average




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
#resave is the starting point

#fit sen slope for each location
library (Kendall)
library(trend)
year=resave$Year[which(resave$loc==cloc)]
y=resave$Avail[which(resave$loc==cloc)]
plot(y)
sen_res <- sens.slope(y)   # years are optional – if omitted they are taken as 1:n
sen_res

#do this for all catchments
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

## Snow fraction ----
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


#Here I am whith my plot, kind of final. I can check trend significance but there is not much to be seen


Myvar=(ppoints[,c(36:40)])
st_geometry(Myvar)<-NULL
Myvar=data.frame(Myvar)
library(corrplot)
res <- cor(Myvar)
round(res, 2)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)






Tsno=data.frame(timeStampx,Snowmelt[["290666"]])
#30 days running mean for a smoother regime and multiplication by 4 for mm/d
dt=4
Tsno[,2]=tsEvaNanRunningMean(Tsno[,2]*4,dt)

Trun=data.frame(timeStampx,Runoff[["290666"]])
Trun[,2]=tsEvaNanRunningMean(Trun[,2]*4,dt)


Tquz=data.frame(timeStampx,quz[["290666"]])
Tquz[,2]=tsEvaNanRunningMean(Tquz[,2]*4,dt)

Tqlz=data.frame(timeStampx,qlz[["290666"]])
Tqlz[,2]=tsEvaNanRunningMean(Tqlz[,2]*4,dt)


TAct=data.frame(timeStampx,ActEvapo[["290666"]])
TAct[,2]=tsEvaNanRunningMean(TAct[,2]*4,dt)

TInf=data.frame(timeStampI,Infil[["290666"]])
TInf[,2]=tsEvaNanRunningMean(TInf[,2]*4,dt)


TSnowf=data.frame(timeStampx,Snow[["290666"]])
TSnowf[,2]=tsEvaNanRunningMean(TSnowf[,2]*4,dt)

TRain=data.frame(timeStampx,Rain[["290666"]])
TRain[,2]=RunningSum(TRain[,2],dt)


TPET=data.frame(timeStampx,PEvapo[["290666"]])
TPET[,2]=tsEvaNanRunningMean(TPET[,2]*4,dt)
#snow
#regime
Rsnow=RegimeFast(Tsno)
plot(Rsnow)
#yearly aggregation

#step1: Location extraction 
Trun=data.frame(timeStampx,Runoff[["290666"]])
Trun[,2]=tsEvaNanRunningMean(Trun[,2]*4,dt)
#keep daily data for speed
names(Trun)[1] <- "date"
# Filter for specific hours
dsel <- hour(Trun$date)
datD <- Trun[(dsel == 12 | dsel == 13),] 
Yagg <- process_data(datD)

wparam=c()
ksv=c()
qra=c()
windowSize=20
dx <- floor(windowSize/2)
l <- length(testdata$year)

for (ii in c(1:l)){
  minindx <- max(ii - dx, 1)
  maxindx <- min(ii + dx, l)
  fdata=testdata[minindx:maxindx,]
  years=unique((fdata$year))
  fit1=budykofit(data=fdata,method="Fu",silent = TRUE)
  wparam=rbind(psave,fit1$Fu[1])
  sely=which(!is.na(match(year(datD$date),years)))
  datYr=datD[sely,2]
  qx=quantile(datYr,.95)
  q95r=sum(datYr[which(datYr>qx)])/sum(datYr)
  kurt<-kurtosis(datYr)
  ksv<-c(ksv,kurt)
  qra=c(qra,q95r)
}


plot(TRain[c(1:10),])
RainT <- process_data(TRain)
RainT$rnsum<-tsEvaNanRunningMean(RainT$val.sum,10)
plot(RainT$year,RainT$val.sum, type="o")
lines(RainT$year,RainT$rnsum)

SnowT <- process_data(TSnowf)
SnowT$rnsum<-tsEvaNanRunningMean(SnowT$val.sum,10)
plot(SnowT$year,SnowT$val.sum+RainT$val.sum, type="o")
lines(SnowT$year,SnowT$rnsum+RainT$rnsum)

ActT <- process_data(TAct)
ActT$rnsum<-tsEvaNanRunningMean(ActT$val.sum,10)
plot(ActT$year,ActT$rnsum, type="o")

PetT <- process_data(TPET)
PetT$rnsum<-tsEvaNanRunningMean(PetT$val.sum,10)
plot(PetT$year,PetT$rnsum, type="o")

LzT <- process_data(Tqlz)
LzT$rnsum<-tsEvaNanRunningMean(LzT$val.sum,10)
plot(LzT$year,LzT$rnsum, type="o")

UzT <- process_data(Tquz)
UzT$rnsum<-tsEvaNanRunningMean(UzT$val.sum,10)
plot(UzT$year,UzT$rnsum, type="o")


RofT <- process_data(Trun)
RofT$rnsum<-tsEvaNanRunningMean(RofT$val.sum,10)
plot(RofT$year,RofT$rnsum, type="o")

Ticido[,2]=tsEvaNanRunningMean(Ticido[,2],120)
RegimeQ=RegimeFast(Ticido)

PetT$val.sum/(RainT$val.sum+SnowT$val.sum)


#Intra annual concentration of P
length(which(TRain$Rain...290666...>1))
hist(TRain[c(1:3000),2])
# Rename columns
names(TRain)[1] <- "date"

# Filter for specific hours
dsel <- hour(TRain$date)
datD <- TRain[(dsel == 12 | dsel == 13),] 

datD[,2]=datD[,2]+TSnowf[(dsel == 12 | dsel == 13),2]
datD[,2][which(datD[,2]<1e-5)]=1e-5


#Budyko curve
ET.P=PetT$rnsum/(RainT$rnsum+SnowT$rnsum)

AET.P=ActT$rnsum/(RainT$rnsum+SnowT$rnsum)

plot(ET.P,AET.P)

#devtools::install_github("tylerbhampton/budykoR")

testdata=data.frame(year=PetT$year,PET.P=ET.P,AET.P=AET.P)

library(budyko)

blankBC
blankBC+coord_cartesian(xlim=c(0,3))

blankBC+geom_point(data=testdata)+coord_cartesian(xlim=c(0,5))

ogbudyko=budyko_sim(method = "fu",silent = T)
blankBC+geom_line(data=ogbudyko)+coord_cartesian(xlim=c(0,5))

blankBC+
  geom_line(data=ogbudyko)+
  geom_point(data=testdata)+
  coord_cartesian(xlim=c(0,5))







wparam=c()
ksv=c()
qra=c()
windowSize=20
dx <- floor(windowSize/2)
l <- length(testdata$year)

for (ii in c(1:l)){
  minindx <- max(ii - dx, 1)
  maxindx <- min(ii + dx, l)
  fdata=testdata[minindx:maxindx,]
  years=unique((fdata$year))
  fit1=budykofit(data=fdata,method="Fu",silent = TRUE)
  wparam=rbind(psave,fit1$Fu[1])
  sely=which(!is.na(match(year(datD$date),years)))
  datYr=datD[sely,2]
  qx=quantile(datYr,.95)
  q95r=sum(datYr[which(datYr>qx)])/sum(datYr)
  kurt<-kurtosis(datYr)
  ksv<-c(ksv,kurt)
  qra=c(qra,q95r)
}

plot(psave[,1])
plot(ksv,type="o")
plot(qra,type="o")
plot(testdata$AET.P)
snowfraction=SnowT$rnsum/(RainT$rnsum+SnowT$rnsum)

plot(snowfraction)

rofffraction=RofT$rnsum/(RofT$rnsum+UzT$rnsum+LzT$rnsum)

plot(rofffraction)

#to be saved for every catchment at yearly timescale:
#direct runoff fraction
#snow fraction
#concentration of precipitation
#total precipitation
#aridity index
#AEP
#w parameter budyko
#ET0



#fit a lognormal distribution
library(fitdistrplus)
descdist(datD$Rain...290666..., boot = 10)
sd(datD$Rain...290666...)/mean(datD$Rain...290666...)
library(moments)
#my measure of concentration can be the kurtosis
kurtosis(datD$Rain...290666...)
walla=fitdist(datD$Rain...290666..., "gamma")
walla$estimate

plot(walla)
RegimeRain=RegimeFast(TRain)
plot(RegimeRain, type="o")
cv=sd(RegimeRain$mean)/mean(RegimeRain$mean)
#fit a random forest to explain changes in the w parameter:

#>           Fu
#> param 2.0500
#> mae   0.0698
#> rsq   0.6139
#> hs    0.0000
sim1=budyko_sim(fit=fit1)
blankBC+
  geom_line(data=ogbudyko)+
  geom_line(data=sim1)+
  geom_point(data=testdata,aes(col=year))+
  coord_cartesian(xlim=c(0,5))

plot(Rqlz)
plot(Rsnow)
deb=Tsno$timeStampx[4]
mois.deb <-  seq(as.Date(deb), by="month", length=12)


Rqz=Rqlz
Rqz$mean=Rqz$mean+Rquz$mean

Rqzr=Rqz
Rqzr$mean=Rqzr$mean+Rrun$mean


Cairo::Cairo(
  20, #length
  15, #width
  file = "D:/tilloal/Documents/01_Projects/RegimeShifts/plots/TicinoRegime.jpg",
  type = "png", #tiff
  bg = "transparent", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)

qlim=c(0,1.2*max(Ziz$mean))
plot(Rqlz$date, Rqlz$mean, type="n", axes=FALSE,ylim=qlim,xaxs="i",yaxs="i",
     xlab = NA, ylab = expression(paste("Q (mm/d)")))

polygon(c(Rqlz$date,rev(Rqlz$date)),c(rep(0,length(Rqlz$mean)),rev(Rqlz$mean)),
        col=alpha("lightblue",.7),border="transparent")
polygon(c(Rqlz$date,rev(Rqlz$date)),c(Rqlz$mean,rev(Rqz$mean)),
        col=alpha("royalblue",.7),border="transparent")
polygon(c(Rqz$date,rev(Rqz$date)),c(Rqz$mean,rev(Rqzr$mean)),
        col=alpha("purple",.7),border="transparent")
# lines(Rqz$date, Rqz$mean*4, col="blue",lwd=2)
# lines(Rqzr$date, Rqzr$mean*4, col="orange",lwd=2)
#lines(Rqzrs$date, Rqzrs$mean*4, col="darkred",lwd=2)
# lines(Rqzrse$date, Rqzrse$mean*4, col="purple",lwd=2)
lines(RegimeQ$date, RegimeQ$mean, col="black",lwd=3)
abline(v = format(mois.deb,"%j"), col="lightgrey", lty=3)
axis(2)
axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=1)
box()
legend("topleft", leg=c("baseflow","sub-surface flow","runoff","Q"),
       lwd=c(4,4,4,2), col=c("lightblue","royalblue","purple","black"),
       cex=1, lty=c(1,1,1,1),bty="n")
dev.off()

#Same plot with rain and snow

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
#Same plot with rain and snow

