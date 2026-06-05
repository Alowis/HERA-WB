
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
wsize=40
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
# Regysave=c()
# Regy_d_save=c()
# Regy_w_save=c()
# resave=c()
# metrics=c()
# 
# ###MODIFICATIONS for soil moisture balance
# for (cloc in cnames){
#   #cloc="311057"
#   print(cloc)
#   #Qep=preprocess_in(timeStampX,Qe,cloc,dsel,dt)
#   #plot(Qep)
#   #Qy1=Q1_sum_s(Qep,2020)
#   data=preprocess_in(timeStampX,Qe,cloc,dsel,dt)
#   #Rqyi<-RegimeFunctionM(dvar=Qe, cloc, yearstart=1951, wsize=1, timeStampX)
#   #quantiles 5% and 95% of regimes
#   Rqyi <- RegimeFunctionDev(dvar = Qe,
#                           cloc = cloc,
#                           yearstart = 1951,
#                           wsize = 40 ,
#                           timeStampX = timeStampX)
#   Regyt=Rqyi$date
#   center_year1=data.frame(day=Rqyi$date)
#   n_years   <- length(yrlist)
#   # `res` will hold mean and SI for each year → 2 columns
#   res_mat   <- matrix(NA_real_, nrow = n_years, ncol = 4,
#                       dimnames = list(NULL, c("mean", "SI","wet-dev","dry-dev")))
#   # `Regy_means` will hold the column of `Rqy1$mean` for each year
#   # (plus the original date column that you already have)
#   Regy_means <- matrix(NA_real_, nrow = length(Regyt), ncol = n_years)
#   
#   Regy_wet <- matrix(NA_real_, nrow = length(Regyt), ncol = n_years)
#   Regy_dry <- matrix(NA_real_, nrow = length(Regyt), ncol = n_years)
#   # Optional: keep the years for later reference
#   year_vec  <- integer(n_years)
#   
#   res=c()
#   pb <- txtProgressBar(min = 0, max = n_years, style = 3)   # lightweight progress bar
#   # t0 <- proc.time()
#   for (i in seq_len(n_years)) {
#     year <- yrlist[i]
#     year_vec[i] <- year
#     
#     # # ---- Year‑specific calculations ---------------------------------
#     dt1=4
#     dt2=30
#     dt=dt1*dt2
#     # Before your year loop, precompute once per cloc:
#     Qs_cache     <- RunningSum(Qe[[cloc]],dt)
#     ts_sec_cache <- as.integer(timeStampX)
#     yr_cache     <- year(timeStampX)
#     mo_cache     <- month(timeStampX)
# 
#     Rqy1 <- RegimeFunctionM_fast(Qs_cache,
#                             ts_sec_cache,
#                             yr_cache,
#                             mo_cache,
#                             yearstart=year,wsize=0)
# 
#     
#     # Rqy1 <- RegimeFunctionM(dvar=Qe, 
#     #                              cloc = cloc,
#     #                              yearstart = 1951,
#     #                              wsize = 40 ,
#     #                              timeStampX = timeStampX)
#     # 
#     
#     yDev1=Rqy1$mean-Rqyi$q95
#     yDev2=Rqy1$mean-Rqyi$q05
#     mw=which(yDev1>0)
#     wdev=length(mw)
#     md=which(yDev2<0)
#     ddev=length(md)
#     
#     Rqy1$wdev=0
#     Rqy1$ddev=0
#     Rqy1$wdev[mw]=1
#     Rqy1$ddev[md]=1
#     
#     # mean for this year (using the pre‑computed `data_const`)
#     if (Qv=="bluewater"){
#       mean_year <- Q1_mean_s(data, year)
#       m1=sum(Rqy1$mean)
#     }
#     
#     if (Qv=="greenwater"){
#       #for soil sum and not mean is used
#       mean_year <- Q1_sum_s(data, year)
#       #for soil, bring all values to positive by adding minim value
#       Rqy1$mean=Rqy1$mean-min(Rqy1$mean)
#       m1=sum(abs(Rqy1$mean))
#     }
# 
#     #for Soil m1=0
#     #m1=0
#     # Seasonality index (SI) – note the use of the pre‑computed `m1_const`
#     #SI <- sum(abs(Rqy1[, 2] - m1/12)) / (m1)
#     
#     SI <- sum(abs(Rqy1[, 2] - m1/12))/ m1
#     
#     # if (m1>0) SI=SI
#     # if (m1<0) SI=-SI
#     #print(SI)
#     # ---- Store results ------------------------------------------------
#     res_mat[i, ] <- c(mean_year, SI, wdev, ddev)  # two‑column result matrix
#     Regy_means[, i] <- Rqy1$mean              # column of means for this year
#     Regy_wet[, i] <- Rqy1$wdev   
#     Regy_dry[, i] <- Rqy1$ddev  
#     setTxtProgressBar(pb, i)
#   }
#   close(pb)
#   # t_original <- proc.time() - t0
#   # cat("Original: ", t_original["elapsed"], "s for", n_years, "years\n")
#   
#   # ---- Assemble final objects -----------------------------------------
#   # 1) `res` as a data.frame (or tibble) with a `year` column
#   res <- data.frame(year = year_vec, res_mat)
#   #plot(res$mean)
#   # 2) `Regy` – original dates + one column per year of `Rqy1$mean`
#   Regy <- cbind(mday = Regyt, Regy_means)
#   Regy=data.frame(Regy)
#   
#   Regy_w <- cbind(mday = Regyt, Regy_wet)
#   Regy_w=data.frame(Regy_w)
#   
#   Regy_d <- cbind(mday = Regyt, Regy_dry)
#   Regy_d=data.frame(Regy_d)
# 
#   colnames(Regy)[c(2:71)]=yrlist
#   colnames(Regy_d)[c(2:71)]=colnames(Regy_w)[c(2:71)]=yrlist
#   Regy$cloc=rep(cloc,length(Regy[,1]))
#   Regy_d$cloc=rep(cloc,length(Regy[,1]))
#   Regy_w$cloc=rep(cloc,length(Regy[,1]))
#   
#   
#   res$tdev=(res$dry.dev+res$wet.dev)/12
#   res$dry.dev=res$dry.dev/12
#   res$wet.dev=res$wet.dev/12
#   runningsd=sqrt(tsEvaNanRunningVariance(res$tdev,30))
#   res=data.frame(res)
#   # o=ecdf(res$mean)
#   # perc=o(res$mean)
#   res$dm=runningsd
#   #statistics
#   baseline_start <- min(res$year, na.rm = TRUE)         # earliest year in the series
#   baseline_years <- baseline_start:(baseline_start + 39) 
#   modern_years <- seq(max(baseline_years)+1,max(res$year))
#   mbase=match(baseline_years,res$year)
#   mmod=match(modern_years,res$year)
#   valvar=c()
#   for (idv in c(2:7)){
#     #idv=4
#     zbas=res[mbase,idv]
#     zmod=res[mmod,idv]
#     # plot(zmod)
#     # # Assuming 'period1_values' and 'period2_values' are your data vectors
#     # wilcox_result <- wilcox.test(zbas, zmod, conf.int = TRUE)
#     # 
#     # # View results
#     # print(wilcox_result$p.value) # p < 0.05 means they are significantly different
#     # print(wilcox_result$estimate) # This is the 'Hodges-Lehmann' estimator of the difference
#     # vals=wilcox_result$p.value
#     
#     # 1. Calculate the observed difference in means (or medians)
#     obs_diff <- mean(zmod) - mean(zbas)
#     
#     # 2. Pool the data and shuffle it 10,000 times
#     pooled_data <- c(zmod, zbas)
#     n1 <- length(zbas)
#     
#     perm_diffs <- replicate(1000, {
#       shuffled <- sample(pooled_data)
#       # Split the shuffled data back into two fake groups
#       mean(shuffled[1:n1]) - mean(shuffled[(n1+1):length(shuffled)])
#     })
#     
#     # 3. Calculate the p-value (how many fake diffs are larger than our real one?)
#     p_val <- mean(abs(perm_diffs) >= abs(obs_diff))
#     
#     
#     valvar=c(valvar,p_val)
#   }
#   names(valvar)=colnames(res)[2:7]
#   
#   res=data.frame(res)
#   # o=ecdf(res$mean)
#   # perc=o(res$mean)
#   # res$dm=c(0,abs(diff(perc)))
#   
#   resp1=res[mbase,]
#   plot=F
#   if (plot==T){ 
#     n=4
#     mean(resp1[,n])
#     qh1=quantile(resp1[,n],.75)
#     ql1=quantile(resp1[,n],.25)
#     qh2=quantile(resp1[,n],.90)
#     ql2=quantile(resp1[,n],.10)
#     plot(res[,n],type="o",lwd=2,pch=19)
#     abline(h=mean(resp1[,n]),lwd=2)
#     abline(h=qh1,lty=2,col="blue",lwd=2)
#     abline(h=qh2,lty=2,col="darkblue",lwd=2)
#     abline(h=ql1,lty=2,col="orange",lwd=2)
#     abline(h=ql2,lty=2,col="red",lwd=2)
#     abline(v=30,lwd=2,col="grey")
#     abline(v=41,lwd=2,col="grey")
#     
#     
#     mean(resp1[,3])
#     qh1=quantile(resp1[,3],.75)
#     ql1=quantile(resp1[,3],.25)
#     qh2=quantile(resp1[,3],.90)
#     ql2=quantile(resp1[,3],.10)
#     plot(res[,3],type="o",lwd=2,pch=19)
#     abline(h=mean(resp1[,3]),lwd=2)
#     abline(h=qh1,lty=2,col="blue",lwd=2)
#     abline(h=qh2,lty=2,col="darkblue",lwd=2)
#     abline(h=ql1,lty=2,col="orange",lwd=2)
#     abline(h=ql2,lty=2,col="red",lwd=2)
#     abline(v=30,lwd=2,col="grey")
#     abline(v=41,lwd=2,col="grey")
#   }
#   
#   #compare periods:
#   
#   resp2<-res[mmod,]
#   
#   
#   #kge style
#   
#   #availability
#   #for green water
#   b=2*mean(resp2$mean)/(mean(resp1$mean)+mean(resp2$mean))
#   #variability
#   v=sd(resp2$mean)/sd(resp1$mean)
#   #seasonality
#   s=mean(resp2$SI)/mean(resp1$SI)
#   #deviation-wet
#   dw=mean(resp2$wet.dev)/mean(resp1$wet.dev)
#   #deviation-dry
#   dd=mean(resp2$dry.dev)/mean(resp1$dry.dev)
#   #deviation-all
#   dad=mean(resp2$tdev)/mean(resp1$tdev)
#  
#   
#   #colnames(res)=c("Qmean","SI","Variability")
#   res$year=yrlist
#   res$loc=rep(cloc,length(res$mean))
#   metrix=data.frame(t(c(avail1=b,var1=v,si1=s,dw=dw,dd=dd,dad=dad,valvar,cloc=cloc)))
#   
#   Regysave=rbind(Regysave,Regy)
#   Regy_d_save=rbind(Regy_d_save,Regy_d)
#   Regy_w_save=rbind(Regy_w_save,Regy_w)
#   metrics=rbind(metrics,metrix)
#   resave<-rbind(resave,res)
#   
# }
Regysave_list <- vector("list", length(cnames))
Regy_d_list   <- vector("list", length(cnames))
Regy_w_list   <- vector("list", length(cnames))
metrics_list  <- vector("list", length(cnames))
resave_list   <- vector("list", length(cnames))

for (ci in seq_along(cnames)) {
  cloc <- cnames[ci]
  print(cloc)
  
  data  <- preprocess_in(timeStampX, Qe, cloc, dsel, dt)
  Rqyi  <- RegimeFunctionDev(dvar = Qe, cloc = cloc, yearstart = 1951,
                             wsize = 40, timeStampX = timeStampX)
  Regyt   <- Rqyi$date
  n_years <- length(yrlist)
  
  res_mat    <- matrix(NA_real_, nrow = n_years, ncol = 4,
                       dimnames = list(NULL, c("mean", "SI", "wet.dev", "dry.dev")))
  Regy_means <- matrix(NA_real_, nrow = length(Regyt), ncol = n_years)
  Regy_wet   <- matrix(NA_real_, nrow = length(Regyt), ncol = n_years)
  Regy_dry   <- matrix(NA_real_, nrow = length(Regyt), ncol = n_years)
  year_vec   <- integer(n_years)
  
  pb <- txtProgressBar(min = 0, max = n_years, style = 3)
  
  for (i in seq_len(n_years)) {
    year       <- yrlist[i]
    year_vec[i] <- year
    
    Rqy1  <- RegimeFunctionM(dvar = Qe, cloc = cloc, yearstart = year,
                             wsize = 0, timeStampX = timeStampX)
    
    yDev1 <- Rqy1$mean - Rqyi$q95
    yDev2 <- Rqy1$mean - Rqyi$q05
    mw    <- which(yDev1 > 0)
    md    <- which(yDev2 < 0)
    wdev  <- length(mw)
    ddev  <- length(md)
    
    Rqy1$wdev     <- 0L
    Rqy1$ddev     <- 0L
    Rqy1$wdev[mw] <- 1L
    Rqy1$ddev[md] <- 1L
    
    if (Qv == "bluewater") {
      mean_year <- Q1_mean_s(data, year)
      m1        <- sum(Rqy1$mean)
    }
    if (Qv == "greenwater") {
      mean_year  <- Q1_sum_s(data, year)
      Rqy1$mean  <- Rqy1$mean - min(Rqy1$mean)
      m1         <- sum(abs(Rqy1$mean))
    }
    
    SI <- sum(abs(Rqy1[, 2] - m1 / 12)) / m1
    
    res_mat[i, ]    <- c(mean_year, SI, wdev, ddev)
    Regy_means[, i] <- Rqy1$mean
    Regy_wet[, i]   <- Rqy1$wdev
    Regy_dry[, i]   <- Rqy1$ddev
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  res         <- data.frame(year = year_vec, res_mat)
  res$tdev    <- (res$dry.dev + res$wet.dev) / 12
  res$dry.dev <-  res$dry.dev / 12
  res$wet.dev <-  res$wet.dev / 12
  res$dm      <- sqrt(tsEvaNanRunningVariance(res$mean_year, 30))
  
  # ── Vectorised permutation tests ───────────────────────────────────
  baseline_start <- min(res$year, na.rm = TRUE)
  baseline_years <- baseline_start:(baseline_start + 39)
  modern_years   <- seq(max(baseline_years) + 1, max(res$year))
  mbase <- match(baseline_years, res$year)
  mmod  <- match(modern_years,   res$year)
  
  valvar <- sapply(2:7, function(idv) {
    zbas        <- res[mbase, idv]
    zmod        <- res[mmod,  idv]
    obs_diff    <- mean(zmod) - mean(zbas)
    pooled_data <- c(zbas, zmod)
    n1          <- length(zbas)
    n_total     <- length(pooled_data)
    
    # all 1000 shuffles in one matrix op instead of a replicate() loop
    idx  <- matrix(sample(rep(seq_len(n_total), 1000)), nrow = n_total)
    grp1 <- colMeans(matrix(pooled_data[idx[1:n1, ]], nrow = n1))
    grp2 <- colMeans(matrix(pooled_data[idx[(n1+1):n_total, ]],
                            nrow = n_total - n1))
    mean(abs(grp1 - grp2) >= abs(obs_diff))
  })
  names(valvar) <- colnames(res)[2:7]
  
  # ── Assemble Regy objects ──────────────────────────────────────────
  Regy   <- data.frame(cbind(mday = Regyt, Regy_means))
  Regy_w <- data.frame(cbind(mday = Regyt, Regy_wet))
  Regy_d <- data.frame(cbind(mday = Regyt, Regy_dry))
  
  colnames(Regy)[2:71]   <- yrlist
  colnames(Regy_w)[2:71] <- colnames(Regy_d)[2:71] <- yrlist
  Regy$cloc   <- cloc
  Regy_w$cloc <- cloc
  Regy_d$cloc <- cloc
  
  # ── KGE-style metrics ──────────────────────────────────────────────
  resp1 <- res[mbase, ]
  resp2 <- res[mmod,  ]
  
  b   <- 2 * mean(resp2$mean)  / (mean(resp1$mean)  + mean(resp2$mean))
  v   <- sd(resp2$mean)        /  sd(resp1$mean)
  s   <- mean(resp2$SI)        /  mean(resp1$SI)
  dw  <- mean(resp2$wet.dev)   /  mean(resp1$wet.dev)
  dd  <- mean(resp2$dry.dev)   /  mean(resp1$dry.dev)
  dad <- mean(resp2$tdev)      /  mean(resp1$tdev)
  
  res$year <- yrlist
  res$loc  <- cloc
  
  metrix <- data.frame(t(c(avail1 = b, var1 = v, si1 = s,
                           dw = dw, dd = dd, dad = dad, valvar, cloc = cloc)))
  
  # ── Accumulate into lists ──────────────────────────────────────────
  Regysave_list[[ci]] <- Regy
  Regy_d_list[[ci]]   <- Regy_d
  Regy_w_list[[ci]]   <- Regy_w
  metrics_list[[ci]]  <- metrix
  resave_list[[ci]]   <- res
}

# ── Single rbind at the end ────────────────────────────────────────────
Regysave    <- do.call(rbind, Regysave_list)
Regy_d_save <- do.call(rbind, Regy_d_list)
Regy_w_save <- do.call(rbind, Regy_w_list)
metrics     <- do.call(rbind, metrics_list)
resave      <- do.call(rbind, resave_list)


mean(as.numeric(metrics$var1))
#names(resave)=c("Year","Avail","SI","wet-dev","dry-dev","Variability","Year2","loc")

#save the stuff
Savelist=list(Regimes=Regysave,Regimes_dry=Regy_d_save,Regimes_wet=Regy_w_save,yearly_res=resave,metrics=metrics)

if (Qv=="bluewater"){
  save(Savelist,file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_blue3.Rdata")
}
if (Qv=="greenwater"){
  save(Savelist,file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_green3.Rdata")
}



# Source - https://stackoverflow.com/a/68442175
# Posted by teunbrand, modified by community. See post 'Timeline' for change history
# Retrieved 2026-03-19, License - CC BY-SA 4.0

library(ggplot2)
library(gtable)
library(grid)






## load Biogeographic regions ----
biogeo <- read_sf(dsn = paste0(hydroDir,"/eea_3035_biogeo-regions_2016/BiogeoRegions2016_wag84.shp"))
biogeof=fortify(biogeo)
st_geometry(biogeof)<-NULL
biogeoregions=raster( paste0(hydroDir,"/eea_3035_biogeo-regions_2016/Biogeo_rasterized_wsg84.tif"))
Gbiogeoregions=as.data.frame(biogeoregions,xy=T)
biogeomatch=inner_join(biogeof,Gbiogeoregions,by= c("PK_UID"="Biogeo_rasterized_wsg84"))
biogeomatch$latlong=paste(round(biogeomatch$x,4),round(biogeomatch$y,4),sep=" ")
biogeo_rivers=right_join(biogeomatch,UpArea, by="latlong")


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

ppc <- st_as_sf(CatO, coords = c("Var1.x", "Var2.x"), crs = 4326)
ppc <- st_transform(ppc, crs = 3035)
# Compute centroids
CatO_centroids <- st_centroid(ppc)
#starting from outputs
Qv="bluewater"
if (Qv=="bluewater"){
 load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_blue3.Rdata")
}
if (Qv=="greenwater"){
 load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_green3.Rdata")
 cwater="Green"
  }

resave=Savelist$yearly_res
metricsave=c()
pb <- txtProgressBar(min = 0, max = length(cnames), style = 3)
i=0
for (cloc in cnames){
i=i+1
#print(cloc)
res=resave[which(resave$loc==cloc),]
# ── Vectorised permutation tests ───────────────────────────────────
baseline_start <- min(res$year, na.rm = TRUE)
baseline_years <- baseline_start:(baseline_start + 39)
modern_years   <- seq(max(baseline_years) + 1, max(res$year))
mbase <- match(baseline_years, res$year)
mmod  <- match(modern_years,   res$year)

valvar <- sapply(2:7, function(idv) {
  zbas        <- res[mbase, idv]
  zmod        <- res[mmod,  idv]
  obs_diff    <- mean(zmod) - mean(zbas)
  pooled_data <- c(zbas, zmod)
  n1          <- length(zbas)
  n_total     <- length(pooled_data)
  
  # all 1000 shuffles in one matrix op instead of a replicate() loop
  idx  <- matrix(sample(rep(seq_len(n_total), 1000)), nrow = n_total)
  grp1 <- colMeans(matrix(pooled_data[idx[1:n1, ]], nrow = n1))
  grp2 <- colMeans(matrix(pooled_data[idx[(n1+1):n_total, ]],
                          nrow = n_total - n1))
  mean(abs(grp1 - grp2) >= abs(obs_diff))
})
# 
# valvar=c()
#   for (idv in c(2:7)){
#     #idv=4
#     zbas=res[mbase,idv]
#     zmod=res[mmod,idv]
#     # 1. Calculate the observed difference in means (or medians)
#     obs_diff <- mean(zmod) - mean(zbas)
# 
#     # 2. Pool the data and shuffle it 10,000 times
#     pooled_data <- c(zmod, zbas)
#     n1 <- length(zbas)
# 
#     perm_diffs <- replicate(1000, {
#       shuffled <- sample(pooled_data)
#       # Split the shuffled data back into two fake groups
#       mean(shuffled[1:n1]) - mean(shuffled[(n1+1):length(shuffled)])
#     })
# 
#     # 3. Calculate the p-value (how many fake diffs are larger than our real one?)
#     p_val <- mean(abs(perm_diffs) >= abs(obs_diff))
# 
# 
#     valvar=c(valvar,p_val)
#   }
names(valvar) <- colnames(res)[2:7]

resp1 <- res[mbase, ]
resp2 <- res[mmod,  ]

b   <- (mean(resp2$mean)  - mean(resp1$mean))

v   <- sd(resp2$mean)        /  sd(resp1$mean)
s   <- mean(resp2$SI)        -  mean(resp1$SI)
dw  <- mean(resp2$wet.dev)   -  mean(resp1$wet.dev)
dd  <- mean(resp2$dry.dev)   -  mean(resp1$dry.dev)
dad <- mean(resp2$tdev)      -  mean(resp1$tdev)


metrix <- (data.frame(avail1 = b, var1 = v, si1 = s,
                         dw = dw, dd = dd, dad = dad, cloc = cloc))
metrix=cbind(metrix,as.data.frame(t(valvar)))

metricsave=rbind(metricsave,metrix)
setTxtProgressBar(pb, i)
}
close(pb)

region=NA
studyvar=c("Avail","SI","Variability","dw","dd","dad")
for (st in studyvar){
  print(st)
  resave=Savelist$yearly_res
  mname=match(st,names(resave))
  names(resave)[mname]="var"
  #match resave with upstream area
  mup=match(as.numeric(resave$loc),UpArea$outlets)
  resave$upa=UpArea$upa[mup]
  
  #match resave with regions
  regios=match(as.numeric(resave$loc),biogeo_rivers$outlets)
  resave$region=biogeo_rivers$code[regios]
  
  #optional - region choice
  if (!is.na(region)){
    resave=resave[which(resave$region=="Atlantic"),]
  }
  
  
  #STEP 1: MAPS OF RESULTS SAVED IN METRICS
  
  
  #plot a map of change in SI
  
  metric=metricsave
  ppoints=CatO_centroids
  metric$outlets=as.numeric(metric$cloc)
  km=match(ppoints$outlets.y,metric$outlets)
  ppoints$sig=NA
  limi=c(-20,20)
  leg="Change"
  if (st=="SI"){
    limi=c(-.1,.1)
    leg="Change (-)"
    ppoints$trend=(as.numeric(metric$si1[km]))
    #plot((as.numeric(metric$si1[km])))
    ppoints$sig=as.numeric(metric$SI[km])
    namevar="SI"
  }
  if (st=="Avail"){
    leg="Change (mm)"
    limi=c(-100,100)
    ppoints$trend=as.numeric(metric$avail1[km])*365.25
    plot(ppoints$trend)
    ppoints$sig=as.numeric(metric$mean[km])
    namevar="Availability"
  }
  if (st=="Variability"){
    limi=c(0,2)
    leg="Change (-)"
    ppoints$trend=as.numeric(metric$var1[km])
    ppoints$sig=as.numeric(metric$dm[km])
    namevar="Variability"
  }
  if (st=="dw"){
    limi=c(-5,5)
    ppoints$trend=as.numeric(metric$dw[km])*100
    ppoints$sig=as.numeric(metric$wet.dev[km])
    namevar="Wet deviations"
    leg="Change (%)"
  }
  if (st=="dd"){
    limi=c(-5,5)
    ppoints$trend=as.numeric(metric$dd[km])*100
    ppoints$sig=as.numeric(metric$dry.dev[km])
    namevar="Dry deviations"
    leg="Change (%)"
  }
  if (st=="dad"){
    limi=c(-10,10)
    ppoints$trend=as.numeric(metric$dad[km])*100
    plot((as.numeric(metric$dw[km])))
    ppoints$sig=as.numeric(metric$tdev[km])
    namevar="Dry + Wet deviations"
    leg="Change (%)"
  }
  
  ppoints_sig=ppoints[which(ppoints$sig<=0.05),]
  library(legendry)
  
  ## slope map ----
  palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
  pbg=ggplot(basemap) +
    geom_sf(fill="white", color=NA) +
    geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=trend),color="transparent",alpha=1,shape=21,stroke=0)+ 
    geom_sf(data=ppoints_sig,aes(geometry=geometry,size=upa),color="black",fill="transparent",alpha=1,shape=1,stroke=0.2)+ 
    geom_sf(fill=NA, color="grey20") +
    scale_x_continuous(breaks=seq(-30,40, by=5)) +
    scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                    sep = " ")),
               breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
    scale_fill_gradientn(
      colors=palet,oob = scales::oob_squish, name=leg, limits=limi) +
   
    #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    labs(x="Longitude", y = "Latitude",
         title = paste0(" Change in mean ", namevar," between baseline and \nrecent period")
    )+
    guides(fill=
               guide_colourbar(barheight = 20, barwidth = .5,reverse=F,
                                  # This ensures the bar looks like a continuous scale with limits
    
                                  # These parameters control the arrow look in some themes
                                  frame.colour = "black",
                                  ticks.colour = NA),
           size= "none")+
    theme_minimal()+
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
    
  ggsave(pbg,file=paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/Map01_",Qv,st,".png"), width=20, height=15, units=c("cm"),dpi=500)
  

  
}
  
region=NA
studyvar=c("mean","SI")
for (st in studyvar){
  print(st)
  resave=Savelist$yearly_res
  resave$Year=resave$year
  mname=match(st,names(resave))
  names(resave)[mname]="var"
  #match resave with upstream area
  mup=match(as.numeric(resave$loc),UpArea$outlets)
  resave$upa=UpArea$upa[mup]
  
  #match resave with regions
  regios=match(as.numeric(resave$loc),biogeo_rivers$outlets)
  resave$region=biogeo_rivers$code[regios]
  
  #optional - region choice
  if (!is.na(region)){
    resave=resave[which(resave$region=="Atlantic"),]
  }
  
#aggregate by year
  var_yAg<- aggregate(resave$var, 
                        by = list(year = resave$year),
                        function(x) c(mean= mean(x,na.rm=T),q1=quantile(x,0.1,na.rm=T),q2=quantile(x,0.9,na.rm=T)))
  var_yAg<-do.call(data.frame, var_yAg)
  
  plot(var_yAg$year,var_yAg[,2],type="o")
  VAv=tsEvaNanRunningMean(var_yAg[,2],10)
  lines(var_yAg$year,VAv, col=2, lwd=2)
  abline(h=mean(var_yAg[,3]))

  # 1. filter the input data --------------------------------------
  #resave=resave[-which(resave_med$Year==1951),]
  BigR<-which(!is.na(match(as.numeric(resave$loc),CatO$outlets.y)))
  resave=resave[BigR,]
  # 2. Define the 30‑year baseline ---------------------------------
  baseline_start <- min(resave$Year, na.rm = TRUE)         # earliest year in the series
  baseline_years <- baseline_start:(baseline_start + 39)   # first 30 calendar years

  # 3. Quantiles on the baseline (per location) --------------------
  loc_q <- resave %>%
    filter(Year %in% baseline_years) %>%                  # reference period only
    group_by(loc) %>%
    summarise(
      med = quantile(var, 0.50, na.rm = TRUE),           # 50 % → will become 0
      q05 = quantile(var, 0.05, na.rm = TRUE),
      q25 = quantile(var, 0.25, na.rm = TRUE),
      q75 = quantile(var, 0.75, na.rm = TRUE),
      q95 = quantile(var, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 4. Attach thresholds to the full data set ----------------------
  df <- resave %>%
    left_join(loc_q, by = "loc")
  
  # 5. Centre Qmean on its median (median → 0) --------------------
  df <- df %>%
    mutate(Q_centered = var - med)               # negative = below median
  
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

  library(dplyr)
  library(tidyr)
  library(purrr)
  library(broom)
  library(trend)
  slope_results <- df %>%
    group_by(loc) %>%
    nest() %>%
    mutate(
      # 1. Calculate the Mean of Avail for each group
      mean_avail = map_dbl(data, ~ mean(.x$var, na.rm = TRUE)),
      
      # 2. Compute the Sen's Slope
      model = map(data, ~ sens.slope(.x$var)),
      abs_slope = map_dbl(model, ~ as.numeric(.x$estimates)),
      
      # 3. Convert to percentage of the mean
      slope_pct = (abs_slope / mean_avail) * 100,
      
      # 4. Extract p-value
      p.value = map_dbl(model, ~ as.numeric(.x$p.value))
    ) %>%
    dplyr::select(loc, mean_avail, abs_slope, slope_pct, p.value)
  
  # ppl <- st_as_sf(CatLR, coords = c("Var1.x", "Var2.x"), crs = 4326)
  # ppl <- st_transform(ppl, crs = 3035)
  ppc <- st_as_sf(CatO, coords = c("Var1.x", "Var2.x"), crs = 4326)
  ppc <- st_transform(ppc, crs = 3035)
  # Compute centroids
  CatO_centroids <- st_centroid(ppc)
  alle=match(colnames(ppl),colnames(CatO_centroids))
  ppoints=CatO_centroids
  slope_results$outlets=as.numeric(slope_results$loc)
  km=match(ppoints$outlets.y,slope_results$outlets)
  ppoints$slope=slope_results$slope_pct[km]
  
  
  maxs=round(max(abs(ppoints$slope)),1)
  ## slope map ----
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
      colors=palet,oob = scales::squish, name="water balance (%)", limits=c(-maxs,maxs)) +
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


  # 8. Sum of upa per year / bin ----------------------------------
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

  area_pct_xtrem=area_pct_by_year %>%
    filter(bin %in% c("above_95", "below_05"))
  
  # 8. 30‑yr centred moving‑average for **each** bin -------------
  trend_area <- area_pct_xtrem%>%
    arrange(Year) %>%
    group_by(bin) %>%
    mutate(
      trend = tsEvaNanRunningMean(pct, 30)
      

    ) %>%
    ungroup()

  # ── Sum of the two bins per year, then smooth it ──────────────────────
  trend_sum <- area_pct_xtrem %>%
    group_by(Year) %>%
    summarise(pct = sum(pct), .groups = "drop") %>%
    arrange(Year) %>%
    mutate(
      bin   = "total",
      trend = tsEvaNanRunningMean(pct, 30)
      

    )
  
  # ── Bind together ─────────────────────────────────────────────────────
  trend_all <-  trend_sum
  
  # ── Baseline quantiles (first 40 years, fixed) ────────────────────────
  baseline_q <- trend_all %>%
    group_by(bin) %>%
    filter(Year <= min(Year) + 39) %>%
    summarise(
      lo = quantile(pct, 0.10, na.rm = TRUE),
      hi = quantile(pct, 0.9, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ── Plot ──────────────────────────────────────────────────────────────
  cols <- c("above_95" = "#2166ac", "below_05" = "#d6604d", "total" = "purple")
  

  # ── Baseline quantile bands (first 40 years) ──────────────────────────
  quantile_pairs <- list(
    list(lo = 0.05, hi = 0.95, alpha = 0.10),
    list(lo = 0.10, hi = 0.90, alpha = 0.15),
    list(lo = 0.25, hi = 0.75, alpha = 0.20),
    list(lo = 0.35, hi = 0.65, alpha = 0.25),
    list(lo = 0.45, hi = 0.55, alpha = 0.30)
  )
  
  # broadcast baseline quantiles across all years explicitly
  baseline_bands <- purrr::map_dfr(quantile_pairs, function(qp) {
    trend_all %>%
      group_by(bin) %>%
      filter(Year <= min(Year) + 39) %>%
      summarise(
        lo    = quantile(pct, qp$lo, na.rm = TRUE),
        hi    = quantile(pct, qp$hi, na.rm = TRUE),
        alpha = qp$alpha,
        .groups = "drop"
      ) %>%
      # join Years explicitly so ribbon has an x axis to draw along
      left_join(distinct(trend_all, Year, bin), by = "bin",
                relationship = "many-to-many")
  })
  
  # baseline median per bin
  baseline_med <- trend_all %>%
    group_by(bin) %>%
    filter(Year <= min(Year) + 39) %>%
    summarise(med = mean(pct, na.rm = TRUE), .groups = "drop")
  
  # ── Plot ───────────────────────────────────────────────────────────────
  cols <- c("above_95" = "#2166ac", "below_05" = "#d6604d", "total" = "black")
  tsize=16
  p <- ggplot(trend_all, aes(x = Year, color = bin, fill = bin)) +
    # outermost band first so inner bands paint on top
    purrr::map(quantile_pairs, function(qp) {
      geom_ribbon(
        data = baseline_bands %>% filter(alpha == qp$alpha),
        aes(x = Year, ymin = lo, ymax = hi, fill = bin),  # x explicit here
        alpha = qp$alpha,
        color = NA,
        inherit.aes = FALSE
      )
    }) +
    geom_hline(
      data = baseline_med,
      aes(yintercept = med, color = bin),
      linetype = "dashed", linewidth = 0.5
    ) +

    geom_point(aes(y = pct),    alpha = 0.6) +
    geom_line(aes(y = pct) ,linewidth = 0.3, alpha = 0.5) +
    geom_line(aes(y = trend), col="red",linewidth = 0.9) +
    scale_color_manual(values = "black") +
    scale_fill_manual( values = "black") +
    scale_y_continuous(
      breaks       = seq(0, 100, by = 5),
      minor_breaks = seq(0, 100, by = 1)
    ) +
    scale_x_continuous(
      breaks       = seq(min(trend_all$Year)+9, max(trend_all$Year)-10, by = 10),
      expand=c(0,0)
 
    ) +
    facet_wrap(~ bin, ncol = 1, scales = "free_y",labeller = as_labeller(c(
      "above_95" = "Above 95th percentile",
      "below_05" = "Below 5th percentile",
      "total"    = paste0("High and low deviations ",st," ",Qv)))) +
    guides(y = guide_axis_minor(), x = guide_axis_minor()) +
    labs(x = "Year", y = "% of total area", color = NULL, fill = NULL) +
    theme_minimal() +
    theme(
      axis.title=element_text(size=tsize),
      axis.text=element_text(size=tsize),
      panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
      legend.title = element_text(size=tsize),
      legend.text = element_text(size=tsize),
      strip.text  = element_text(size = tsize),
      legend.position = "none",
      axis.ticks.length = unit(4,  "pt"),        # major tick length
      axis.minor.ticks.length = rel(0.5)         # minor = half of major (base R ≥ 4.3)
    )
  
  print(p)
  ggsave(p,filename=paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/Fig1new_v3_",Qv,st,".png"), width=20, height=15, units=c("cm"),dpi=500)
  
  # 9. Give a negative sign to the three “below‑median” bins ----------
  
  plot_df <- area_pct_by_year %>%
    mutate(
      signed_upa = ifelse(bin %in% c("below_05", "05_25", "25_med"),
                          -pct,                # negative side
                          pct)
    )# positive side
  
  # 8. 30‑yr centred moving‑average for **each** bin -------------
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


  # 10. Compute 10-year Running Mean (instead of Linear Trend)
  
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


  # 11. Pivot and Combine Bins -----------------------------------

  # We pivot to wide format to easily sum the "combined" categories
  trend_wide <- trend_fitted %>%
    dplyr::select(Year, bin, trend) %>%
    pivot_wider(
      names_from = bin,
      values_from = trend,
      names_prefix = "trend_"
    )

  # --- 11a. Middle Bins Plotting Data ---
  trend_plot1 <- trend_wide %>%
    dplyr::select(Year, trend_25_med, trend_med_75) %>%
    pivot_longer(
      cols = starts_with("trend_"),
      names_to = "bin",
      values_to = "trend",
      names_prefix = "trend_"
    ) %>%
    mutate(signif = "MA") # Placeholder for the label
  
  # --- 11b. Combined Bins (Aggregated Extremes) ---
  trend_combined <- trend_wide %>%
    mutate(
      # Combine the numeric running means
      below_comb = trend_25_med ,
      above_comb = trend_med_75 ,
      below_comb2 = trend_25_med + trend_05_25 ,
      above_comb2 = trend_med_75 + trend_75_95 
    ) %>%
    dplyr::select(Year, below_comb, above_comb, below_comb2,above_comb2) %>%
    pivot_longer(
      cols = c(below_comb, above_comb, below_comb2, above_comb2),
      names_to = "combo_bin",
      values_to = "trend"
    ) %>%
    mutate(signif = "MA")
  


  # 12. **Re‑order the factor** -----------------
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
  
  


  
  # 13. colour palette ------------------------------------
  # A smooth gradient from deep red (most negative) → orange → yellow → light‑blue → deep blue (most positive)
  pal_fun <- colorRampPalette(c("#D73027", "#F46D43", "#FEE090",
                                "#FEE090", "#4575B4", "#313695"))
  my_palette <- setNames(pal_fun(6),
                         c("below_05", "05_25", "25_med",
                           "med_75", "75_95", "above_95"))
  
  my_palette2 <- c( "above_comb"="lightblue", "above_comb2"= "darkblue","below_comb"="gold", "below_comb2"="darkred")
  
  # 14. Plot – stacked bar using geom_bar(stat = "identity") ------
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
      subtitle = paste0(Qv," ",st),
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
                       labels=c("05-25","25_med",
                                "med_75","75-95")) +
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
 
  ggsave(paste0("D:/tilloal/Documents/01_Projects/RegimeShifts/plots/Fig1_v5_",Qv,st,".png"), width=30, height=15, units=c("cm"),dpi=500)
  

}
  #loop for deviation of water avail
domain_evolist=list()
map1<-list()
plt1<-list()
plt2<-list()
dname=c("Dry","Wet")
for (d in c(1,2)){
    dev=d+1
    devname=dname[d]
    regimedev=Savelist[[dev]]
    mupa=(((match(regimedev$cloc,UpArea$outlets))))
    regimedev$upa=UpArea$upa[mupa]
    
    #remove large catchments
    BigR<-which(!is.na(match(as.numeric(regimedev$cloc),CatO$outlets.y)))
    regimedev=regimedev[BigR,]
    
    regime_long <- regimedev %>%
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
    
    
    reg_location_yearly <- regime_long %>%
      group_by(cloc, Year) %>%
      summarise(
        # mean(Value) calculates (count of 1s) / (total months recorded)
        pct_months_active = mean(Value, na.rm = TRUE) * 100,
        
        # Keep the upa (area) for each location
        upa = first(upa),
        
        .groups = "drop"
      )
    
    # 1. Calculate the 30-year running mean for every location
    reg_location_trends <- reg_location_yearly %>%
      arrange(cloc, Year) %>%
      group_by(cloc) %>%
      mutate(
        # Smoothing the 'percentage of active months' over 30 years
        trend_30y = tsEvaNanRunningMean(pct_months_active, 30)
      ) %>%
      ungroup()
    
    global_annual_trend <- reg_location_yearly %>%
      group_by(Year) %>%
      summarise(
        # Average of the percentages, weighted by the size of the location
        avg_pct_active = weighted.mean(pct_months_active, w = upa, na.rm = TRUE)
      )
    global_conf <- global_annual_trend %>%
      filter(Year %in% baseline_years) %>%
      summarise(
        base_low  = quantile(avg_pct_active, 0.05),
        #base_mean  = me(boot_means),
        base_med  = quantile(avg_pct_active, 0.500),
        base_high = quantile(avg_pct_active, 0.95),
        .groups = "drop"
      )
    
    

    # 2. Identify which locations have the highest increase in the trend
    # We compare the end of the time series to the beginning (of the smoothed trend)
    trend_summary <- reg_location_trends %>%
      filter(!is.na(trend_30y)) %>%
      group_by(cloc) %>%
      summarise(
        start_trend = first(trend_30y),
        end_trend   = last(trend_30y),
        total_change = end_trend - start_trend,
        upa = first(upa),
        .groups = "drop"
      ) %>%
      arrange(desc(total_change))
    
    lmax=trend_summary$cloc[which.max(trend_summary$total_change)]

    cat("bootstrapping")
    set.seed(123)
    n_boots <- 100 # Increased for more stable Confidence Intervals

    # 1. Filter for baseline years across all locations
    # 2. Group by cloc to keep the bootstrapping separate
    all_loc_conf <- reg_location_yearly %>%
      filter(Year %in% baseline_years) %>%
      group_by(cloc) %>%
      reframe(
        # For each cloc, generate 100 bootstrap means
        boot_means = replicate(n_boots, {
          resample <- sample(pct_months_active, replace = TRUE)
          mean(resample, na.rm = TRUE)
        })
      ) %>%
      # 3. Calculate the CI for each cloc based on its specific bootstrap distribution
      group_by(cloc) %>%
      summarise(
        base_low  = quantile(boot_means, 0.05),
        #base_mean  = me(boot_means),
        base_med  = quantile(boot_means, 0.500),
        base_high = quantile(boot_means, 0.95),
        .groups = "drop"
      )
    
    # 1. Join the baseline CI limits to your actual trend data
    trend_analysis <- reg_location_trends %>%
      left_join(all_loc_conf, by = "cloc") %>%
      mutate(
        # Check if the current 10y trend is higher than the historical 97.5th percentile
        is_beyond_high = trend_30y > base_high,
        # Check if it is lower than the 2.5th percentile
        is_beyond_low  = trend_30y < base_low
      )
    
    # 1. Define your baseline end year
    baseline_end <- 1990
    # 2. Identify the FIRST year AFTER the baseline that the trend broke out
    breakout_years_post_baseline <- trend_analysis %>%
      # Filter for years strictly AFTER the reference period
      filter(Year > baseline_end) %>% 
      # Filter only for years where the trend was outside the CI
      filter(is_beyond_high | is_beyond_low) %>%
      group_by(cloc) %>%
      summarise(
        year_of_breakout = min(Year),
        type_of_breakout = ifelse(first(is_beyond_high), "Significantly more", "Significantly less"),
        value_at_breakout = first(trend_30y),
        .groups = "drop"
      )
    
    
    # Parameters
    min_years <- 30
    study_end <- 2020
    
    sustained_breakouts <- trend_analysis %>%
      filter(Year > baseline_end) %>%
      arrange(cloc, Year) %>%
      group_by(cloc) %>%
      mutate(
        # Create a status string to identify the state of each year
        status = case_when(
          is_beyond_high ~ "Above",
          is_beyond_low  ~ "Below",
          TRUE           ~ "Inside"
        ),
        # Identify consecutive runs of the SAME status
        run_id = rleid(status)
      ) %>%
      group_by(cloc, run_id) %>%
      mutate(
        run_length = n(),
        is_active_at_end = any(Year == study_end)
      ) %>%
      # Filter: Must be outside (Above or Below) AND (long enough OR active at end)
      filter(
        status != "Inside",
        (run_length >= min_years | is_active_at_end == TRUE)
      ) %>%
      # 2. Final Aggregation per Location
      group_by(cloc) %>%
      summarise(
        breakout_year = min(Year),
        direction     = first(status), # "Above" or "Below"
        is_active     = any(Year == study_end),
        run_length    = first(run_length),
        final_trend   = last(trend_30y),
        upa           = first(upa),
        .groups = "drop"
      )
  
    comparison <- breakout_years_post_baseline %>%
      dplyr::select(cloc, first_ever = year_of_breakout) %>%
      left_join(sustained_breakouts %>%dplyr::select(cloc, sustained = breakout_year), by = "cloc") %>%
      mutate(
        delay = sustained - first_ever,
        is_never_sustained = is.na(sustained)
      )
    
    breakout_plots=full_join(trend_summary,sustained_breakouts,by="cloc")
    breakout_unik=right_join(trend_summary,sustained_breakouts,by="cloc")
    
    regime_shift_summary <- breakout_unik %>%
      filter(is_active == TRUE) %>%
      summarise(
        count_locs = n(),
        total_area_km2 = sum(upa.x, na.rm = TRUE),
        # Assuming you have the total domain area calculated elsewhere
        pct_of_total_domain = (sum(upa.x) / sum(trend_summary$upa)) * 100
      )
    
    
    ppoints=CatO_centroids
    breakout_plots$outlets=as.numeric(breakout_plots$cloc)
    km=match(ppoints$outlets.y,breakout_plots$outlets)
    ppoints$trend=breakout_plots$total_change[km]
    ppoints$sig=NA
    ppoints$sig=breakout_plots$is_active[km]
    ppoints_sig=ppoints[which(ppoints$sig==TRUE),]
    
    ## slope map ----
    palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = T, fixup = TRUE))
    map1[[d]]<-ggplot(basemap) +
      geom_sf(fill="white", color=NA) +
      geom_sf(data=ppoints,aes(geometry=geometry,size=upa,fill=trend),color="transparent",alpha=1,shape=21,stroke=0)+ 
      geom_sf(data=ppoints_sig,aes(geometry=geometry,size=upa),color="black",fill="transparent",alpha=1,shape=1,stroke=0)+ 
      geom_sf(fill=NA, color="grey20") +
      scale_x_continuous(breaks=seq(-30,40, by=5)) +
      scale_size(range = c(1, 6), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                      sep = " ")),
                 breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
      scale_fill_gradientn(
        colors=palet,oob = scales::squish, name="deviation", limits=c(-10,10)) +
      #scale_fill_manual(values=c("1"="royalblue","2"="lightblue"))+
      coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
      labs(x="Longitude", y = "Latitude",
        title = paste0(" Deviation from 1951-1990",devname," Baseline")
      )+
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
    
    # 1. Filter for post-baseline period and categorize every year
    area_evolution_simple <- trend_analysis %>%
      filter(Year > 1990) %>%
      mutate(
        status = case_when(
          is_beyond_high ~ "Above CI",
          is_beyond_low  ~ "Below CI",
          TRUE           ~ "Stable (Within CI)"
        )
      ) %>%
      # 2. Sum the area (upa) for each status per year
      group_by(Year, status) %>%
      summarise(total_area = sum(upa, na.rm = TRUE), .groups = "drop")
    
    
    
    plt1[[d]]<-ggplot(area_evolution_simple, aes(x = Year, y = total_area, fill = status)) +
      geom_area(alpha = 0.85, size = 0.1, color = "white") +
      scale_fill_manual(values = c(
        "Above CI"           = "red", # Blue for wetter
        "Below CI"           = "blue", # Red for dryer
        "Stable (Within CI)" = "#bdbdbd"  # Grey for reference state
      )) +
      theme_minimal() +
      labs(
        title = paste0("Total Domain Area Deviating from 1951-1990",devname," Baseline"),
        subtitle = "Evolution of spatial footprint outside of historical Confidence Intervals",
        x = "Year",
        y = "Total Area (upa)",
        fill = "Status"
      ) +
      scale_y_continuous(labels = scales::comma) +
      theme(legend.position = "bottom")
    
    # 1. Flag sustained years for every location
    area_evolution <- trend_analysis %>%
      filter(Year > baseline_end) %>%
      arrange(cloc, Year) %>%
      group_by(cloc) %>%
      mutate(
        status = case_when(
          is_beyond_high ~ "Sustained Above",
          is_beyond_low  ~ "Sustained Below",
          TRUE           ~ "Stable"
        ),
        run_id = rleid(status)
      ) %>%
      group_by(cloc, run_id) %>%
      mutate(
        run_length = n(),
        is_active_at_end = any(Year == study_end),
        # A year is only counted as 'Sustained' if it meets your criteria
        is_sustained = (status != "Stable") & (run_length >= min_years | is_active_at_end)
      ) %>%
      # 2. Sum the area (upa) for each status per year
      group_by(Year, status, is_sustained) %>%
      summarise(total_area = sum(upa, na.rm = TRUE), .groups = "drop") %>%
      # 3. Consolidate: if not sustained, it counts as 'Stable'
      mutate(display_status = ifelse(is_sustained, status, "Stable")) %>%
      group_by(Year, display_status) %>%
      summarise(area_km2 = sum(total_area), .groups = "drop")
    
    
    #save all these plots for all variables
    
    domain_evolution <- reg_location_trends %>%
      # left_join(all_loc_conf, by = "cloc") %>%
      group_by(Year) %>%
      summarise(
        domain_yval     = weighted.mean(pct_months_active, w = upa, na.rm = TRUE),
        domain_trend     = weighted.mean(trend_30y, w = upa, na.rm = TRUE),
        # domain_base_low  = weighted.mean(base_low, w = upa, na.rm = TRUE),
        # domain_base_high = weighted.mean(base_high, w = upa, na.rm = TRUE),
        # # This represents the "average of the averages" from the baseline
        # domain_base_med  = weighted.mean(base_med, w = upa, na.rm = TRUE), 
        .groups = "drop"
      )
    domain_evolution$domain_base_med=global_conf$base_med
    domain_evolution$domain_base_low=global_conf$base_low
    domain_evolution$domain_base_high=global_conf$base_high
    
    # Calculate the single static mean from the baseline period to use as a constant line
    baseline_mean_val <- domain_evolution %>% 
      filter(Year %in% baseline_years) %>% 
      summarise(m = mean(domain_base_med)) %>% 
      pull(m)
    
    lpl=c(0,ceiling(max(domain_evolution$domain_yval)))
    plt2[[d]]<-ggplot(domain_evolution, aes(x = Year)) +
      # 1. The Baseline Confidence Ribbon
      geom_ribbon(aes(ymin = domain_base_low, ymax = domain_base_high), 
                   fill = "grey80", alpha = 0.4) +
      
      # 2. Horizontal Baseline Mean (The "Normal" reference)
      geom_hline(yintercept = baseline_mean_val, color = "gray", 
                 linetype = "solid", size = 0.8) +
      
      # 3. Horizontal CI Bounds
      geom_hline(aes(yintercept = domain_base_low), color = "gray20", 
                 linetype = "dashed", alpha = 0.6) +
      geom_hline(aes(yintercept = domain_base_high), color = "gray20", 
                 linetype = "dashed", alpha = 0.6) +
      scale_y_continuous(limits = lpl, expand=c(0,0))+
      scale_x_continuous(expand=c(0,0))+
      
      # 4. The actual domain-wide evolution
      geom_line(aes(y = domain_trend), color = "red", size = 1.2) +
      geom_line(aes(y = domain_yval), color = "black", size = 1) +
      
      theme_minimal() +
      labs(
        title = paste0("Percentage of area with,",devname,"deviation"),
        subtitle = "Solid Black: Yearly value | Red: 30y Trend | Dashed grey: 95% CI | Solid grey: 1951-1980 mean",
        y = "Percentage of area",
        x = "Year"
      )
    
    domain_evolist=c(domain_evolist,list(domain_evolution))
    
    
    
}


plt2[[2]]
map1[[2]]

all_lists <- list(map1, plt1, plt2)
# --- File 1: All plots with ID 1 ---
pdf("Plots_dry.pdf", width = 11, height = 8.5) # Landscape orientation
for (i in 1:length(all_lists)) {
  # Print the 1st element of the i-th list
  print(all_lists[[i]][[1]]) 
}
dev.off()

# --- File 2: All plots with ID 2 ---
pdf("Plots_wet.pdf", width = 11, height = 8.5)
for (i in 1:length(all_lists)) {
  # Print the 2nd element of the i-th list
  print(all_lists[[i]][[2]])
}
dev.off()

domain_evomix=domain_evolist[[1]]+domain_evolist[[2]]
domain_evomix[,1]=domain_evolist[[1]][,1]

# Calculate the single static mean from the baseline period to use as a constant line
baseline_mean_val <- domain_evomix %>% 
  filter(Year %in% baseline_years) %>% 
  summarise(
  base_low  = quantile(domain_yval, 0.05),
#base_mean  = me(boot_means),
  base_med  = quantile(domain_yval, 0.500),
  base_high = quantile(domain_yval, 0.95),
  .groups = "drop"
) 

domain_evomix$domain_base_high=baseline_mean_val$base_high
domain_evomix$domain_base_low=baseline_mean_val$base_low
domain_evomix$domain_base_med=baseline_mean_val$base_med

lpl=c(0,ceiling(max(domain_evomix$domain_yval)))
ggplot(domain_evomix, aes(x = Year)) +
  # 1. The Baseline Confidence Ribbon
  geom_ribbon(aes(ymin = domain_base_low, ymax = domain_base_high), 
              fill = "grey80", alpha = 0.4) +
  
  # 2. Horizontal Baseline Mean (The "Normal" reference)
  geom_hline(aes(yintercept = domain_base_med), color = "gray", 
             linetype = "solid", size = 0.8) +
  
  # 3. Horizontal CI Bounds
  geom_hline(aes(yintercept = domain_base_low), color = "gray20", 
             linetype = "dashed", alpha = 0.6) +
  geom_hline(aes(yintercept = domain_base_high), color = "gray20", 
             linetype = "dashed", alpha = 0.6) +
  scale_y_continuous(limits = lpl, expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  
  # 4. The actual domain-wide evolution
  geom_line(aes(y = domain_trend), color = "red", size = 1.2) +
  geom_line(aes(y = domain_yval), color = "black", size = 1) +
  
  theme_minimal() +
  labs(
    title = "Percentage of area with wet and dry deviation",
    subtitle = "Solid Black: Yearly value | Red: 30y Trend | Dashed grey: 95% CI | Solid grey: 1951-1980 mean",
    y = "Percentage of area",
    x = "Year"
  )






# 2. Flag extremes by matching cloc AND Month
monthly_extremes <- regime_long %>%
  group_by(Year, mday) %>%
  summarise(
    monthly_impact = weighted.mean(Value, w = upa, na.rm = TRUE),
    .groups = "drop"
  )

# Assuming your data has: cloc, Year, Month, is_extreme (0 or 1)
location_si_evolution <- regime_long %>%
  group_by(cloc, Year) %>%
  summarise(
    # R is the total number of extreme months in that year for that location
    R = sum(Value),
    # Walsh and Lawler SI Formula
    # We filter for R > 0 to avoid dividing by zero
    SI = if(R > 0) (1/R) * sum(abs(Value - (R/12))) else NA,
    
    # Keep the upa (area) for each location
    upa = first(upa),
    .groups = "drop"
  )

# 1. Calculate the 30-year running mean for every location
SI_location_trends <- location_si_evolution %>%
  arrange(cloc, Year) %>%
  group_by(cloc) %>%
  mutate(
    # Smoothing the 'percentage of active months' over 30 years
    trend_30y = tsEvaNanRunningMean(SI, 30)
  ) %>%
  ungroup()


# Calculate SI for every location and year
# This is fast because it's vectorized
location_annual_si <- regime_long %>%
  group_by(cloc, Year) %>%
  summarise(
    R = sum(Value),
    # If R=0, SI is undefined. We handle this by making it NA
    SI = if(R > 0) (1/R) * sum(abs(Value - (R/12))) else NA,
    .groups = "drop"
  ) %>%
  mutate(Period = ifelse(Year <= 1990, "Baseline", "Modern"))

library(dplyr)
library(tidyr)

si_confidence_intervals <- location_annual_si %>%
  filter(!is.na(SI)) %>% 
  group_by(cloc, Period) %>%
  summarise(
    mean_si = mean(SI),
    sd_si = sd(SI),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Standard Error
    se = sd_si / sqrt(n),
    # 95% Confidence Interval Bounds
    ci_lower = mean_si - (1.96 * se),
    ci_upper = mean_si + (1.96 * se)
  ) %>%
  # Optional: Ensure CI bounds don't go below 0 (SI cannot be negative)
  mutate(ci_lower = pmax(0, ci_lower))


#now combine the two anomalies


  # 
  # regimedry=Savelist[[2]]
  # mupa=(((match(regimedry$cloc,UpArea$outlets))))
  # regimedry$upa=UpArea$upa[mupa]
  # 
  # #remove large catchments
  # BigR<-which(!is.na(match(as.numeric(regimedry$cloc),CatO$outlets.y)))
  # regimedry=regimedry[BigR,]
  # 
  # head(regimedry)
  # 
  # library(dplyr)
  # library(tidyr)
  # 
  # regime_long <- regimedry %>%
  #   # 1. Pivot all columns that are years (1951 to 2020)
  #   #    We use matches() with a regex to grab all 4-digit numeric columns
  #   pivot_longer(
  #     cols = matches("^[0-9]{4}$"), 
  #     names_to = "Year", 
  #     values_to = "Value"
  #   ) %>%
  #   # 2. Convert Year to an integer for better sorting/math
  #   mutate(Year = as.integer(Year)) %>%
  #   # 3. Reorder columns for clarity
  #   dplyr::select(cloc, Year, mday, upa, Value)
  # 
  # # View the result
  # head(regime_long)
  # 
  # library(dplyr)
  # library(tidyr)
  # 
  # weighted_regime <- regimedry %>%
  #   # 1. Pivot to long format to get Year as a variable
  #   pivot_longer(
  #     cols = matches("^[0-9]{4}$"), 
  #     names_to = "Year", 
  #     values_to = "Value"
  #   ) %>%
  #   mutate(Year = as.integer(Year)) %>%
  #   # 2. Group by the time variables (Year and Day)
  #   group_by(Year, mday) %>%
  #   summarise(
  #     # Total UPA of all locations for this year/day
  #     total_area = sum(upa, na.rm = TRUE),
  #     # UPA of only the locations where Value is 1
  #     active_area = sum(upa[Value == 1], na.rm = TRUE),
  #     # Calculate weighted percentage
  #     pct_active = (active_area / total_area) * 100,
  #     .groups = "drop"
  #   )
  # 
  # # View the results
  # head(weighted_regime)
  # 
  # mean(weighted_regime$pct_active)
  # 
  # annual_regime_mean_d <- weighted_regime %>%
  #   group_by(Year) %>%
  #   summarise(
  #     mean_pct_active = mean(pct_active, na.rm = TRUE),
  #     # Optional: also calculate the max activity seen in a single day that year
  #     max_pct_active = max(pct_active, na.rm = TRUE),
  #     .groups = "drop"
  #   )
  # 
  # # View results
  # head(annual_regime_mean_d)
  # 
  # annual_regime_mean_d <- annual_regime_mean_d %>%
  #   arrange(Year) %>%
  #   mutate(
  #     trend_10y = tsEvaNanRunningMean(mean_pct_active, 10)
  #   )
  # 
  # # Quick plot to see the trend
  # library(ggplot2)
  # 
  # ggplot(annual_regime_mean_d, aes(x = Year, y = mean_pct_active)) +
  #   geom_line(alpha = 0.3) + # Raw annual values
  #   geom_line(aes(y = trend_10y), color = "red", size = 1) + # 10y Moving Average
  #   theme_minimal() +
  #   labs(title = "Annual Mean Percentage of Area with dry Events",
  #        subtitle = "Blue line: 10-year running mean (Rtseva)",
  #        y = "Mean % Area Active")
  # 

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
    replicate(n_boots, calc_qs(var), simplify = FALSE) %>% bind_rows()
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




#trend stuff


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

# 4. View results








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

