# Function script

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
library(data.table)
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

XCon<-function(data){
  qx=quantile(data,.95, na.rm=T)
  q95r=sum(data[which(data>qx)])/sum(data)
}

plotRegime <- function(data, catch, seuil=NULL){
  #~ data <- as.simu(data)
  names(data)=c("date","Q")
  if( is.null(seuil) ) seuil <- -Inf
  ## Suppression des lignes sans couples Qobs/Qsim
  data <- subset( data, !is.na(Q))
  ## Preparation de la sequence de mois
  deb=data$date[4]
  mois.deb <-  seq(as.Date(deb), by="month", length=12)
  
  period=paste0(format(range(data$date)[1],"%b %Y"),"-",format(range(data$date)[2],"%b %Y"))
  main = bquote(.(catch)~ " Regime | Period: "~ .(period))
  jours <- as.numeric(format(data$date,"%j"))
  ## Iddinces des lignes de meme categorie
  ind.j <- tapply(seq(length(jours)), jours, c)
  ind.j <- ind.j[-366]
  ## Calcul
  Qc <- data.frame(date=as.numeric(names(ind.j)),
                   mean=sapply(ind.j, function(x) mean(data$Q[x], na.rm=TRUE)),
                   q10=sapply(ind.j, function(x) quantile(data$Q[x],0.1, na.rm=TRUE)),
                   q90=sapply(ind.j, function(x) quantile(data$Q[x],.9, na.rm=TRUE)))
  
  qlim=c(0,1.1*max(Qc[,4]))
  plot(Qc$date, Qc$mean, type="n", axes=FALSE, ylim=qlim,
       xlab = NA, ylab = expression(paste("Debit (",m^3/s,")")))
  mtext(main,3,font = 2,line = 0.5,cex = .75)
  lines(Qc$date, Qc$mean, col="darkblue",lwd=2)
  lines(Qc$date, Qc$q10, col="blue",lwd=1,lty=2)
  lines(Qc$date, Qc$q90, col="blue",lwd=1,lty=2)
  abline(v = format(mois.deb,"%j"), col="lightgrey", lty=3)
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q10,rev(Qc$q90)),
          col=alpha("lightblue",.4),border="transparent")
  axis(2)
  axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=.8)
  box()
  legend("topleft", leg=c("mean","Q90","Q10"),
         lwd=c(2,2,1,2), col=c("darkblue","blue","blue"),
         cex=.8, lty=c(1,2,2),bty="n")
  return(Qc)
}

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

process_data <- function(data_frame) {
  # Ensure necessary libraries are loaded
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    install.packages("lubridate")
  }
  
  library(dplyr)
  library(lubridate)
  
  # Create a 'year' column from the timestamp
  data_frame$year <- year(data_frame[,1])
  
  # Aggregate data by year
  Sagg <- aggregate(list(val = data_frame[,2]),
                    by = list(year = data_frame$year),
                    FUN = function(x) c(sum = sum(x, na.rm = TRUE)))
  
  # Convert the list columns to data frame columns
  Sagg <- do.call(data.frame, Sagg)
  
  return(Sagg)
}

p_regime <- function(data, yearstart, wsize=29) {
  
  #ly <- length(years)
  yr=yearstart
  yrwindow <- yr:(yr + wsize)
  subs <- which(!is.na(match(year(as.Date(data[,1])), yrwindow)))
  data<-data[subs,]
  
  
  # Compute Qs
  data$Qs <- data[,2]
  #specific discharge
  # data$Qs <- data[,2]/ (cupa * 1000000) * 1000 * 3600 * 24
  
  
  #mean
  qm <- mean(data[,2])
  
  # Apply running mean again
  data[,2]=tsEvaNanRunningMean(data[,2],30)
  
  Regi <- RegimeFast(data)
  
  #Seasonality index
  SI<-sum(abs(Regi[,2]-qm))/(qm*365)
  # Save parameters
  params <- c(mean=qm*365,SI=SI)
  
  return(params)
}

process_precip <- function(data_frame) {
  # Ensure necessary libraries are loaded
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    install.packages("lubridate")
  }
  
  library(dplyr)
  library(lubridate)
  library(ineq)
  
  # Create a 'year' column from the timestamp
  data_frame$year <- year(data_frame[,1])
  
  # Aggregate data by year
  Sagg <- aggregate(list(val = data_frame[,2]),
                    by = list(year = data_frame$year),
                    FUN = function(x) c(sum = sum(x, na.rm = TRUE),gini=Gini(x)))
  
  # Convert the list columns to data frame columns
  Sagg <- do.call(data.frame, Sagg)
  
  Seasonality_1<-p_regime(data_frame, 1951)
  Seasonality_2<-p_regime(data_frame, 1991)
  
  
  return(list(yagg=Sagg,seaonal=data.frame(Seasonality_1,Seasonality_2)))
}

process_frac <- function(data_frame) {
  # Ensure necessary libraries are loaded
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    install.packages("lubridate")
  }
  
  library(dplyr)
  library(lubridate)
  
  # Create a 'year' column from the timestamp
  data_frame$year <- year(data_frame[,1])
  
  # Aggregate data by year
  Sagg <- aggregate(list(val = data_frame[,2]),
                    by = list(year = data_frame$year),
                    FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
  
  # Convert the list columns to data frame columns
  Sagg <- do.call(data.frame, Sagg)
  
  return(Sagg)
}

budykofit<-function (data, method, dif = "nls", res = NULL, hshift = FALSE, 
                     hs = NULL, silent = FALSE) 
{
  if (!(tolower(method) %in% c("wang-tang", "fu", 
                               "turc-pike", "zhang")) & silent == FALSE) {
    print("Error: unrecognized method")
  }
  else {
    if (silent == FALSE) {
      print(paste0("Budyko Fit: method = ", method))
    }
    if (is.null(res)) {
      res = 0.01
    }
    data = data[, c("AET.P", "PET.P")]
    data = subset(data, !is.na(AET.P))
    data = subset(data, PET.P <= 20)
    data = subset(data, !(AET.P > PET.P))
    data = data[order(data$PET.P), ]
    if (dif == "nls") {
      formula = list(fu = "AET.P~1+PET.P-(1+(PET.P)^p)^(1/p)", 
                     `turc-pike` = "AET.P~(1+(PET.P)^(-p))^(-1/p)", 
                     `wang-tang` = "AET.P~(1+PET.P-((1+PET.P)^2-4*p*(2-p)*PET.P)^(1/2))/(2*p*(2-p))", 
                     zhang = "AET.P~(1+p*PET.P)/(1+p*PET.P+(PET.P^(-1)))")[[tolower(method)]]
      startval = list(fu = list(p = 2), `turc-pike` = list(p = 1.9), 
                      `wang-tang` = list(p = 0.59), zhang = list(p = 1.4))[[tolower(method)]]
      budykoNLS = stats::nls(formula = formula, start = startval, 
                             data = data)
      paramval = summary(budykoNLS)$coefficients[[1]]
      fitmae = mean(abs(stats::residuals(budykoNLS)))
      fitrsq = summary(stats::lm(data$AET.P ~ fitted(budykoNLS)))$r.squared
      fitdat = as.data.frame(c(param = paramval, mae = round(fitmae, 
                                                             4), rsq = round(fitrsq, 4), hs = 0))
    }
    if (dif %in% c("rsq", "mae")) {
      if (tolower(method) == "fu") {
        testval = c(seq(1, 3, res), seq(3, 10, min(0.1, 
                                                   res * 10)))
      }
      if (tolower(method) == "turc-pike") {
        testval = c(seq(0, 2, res), seq(2, 10, min(0.1, 
                                                   res * 10)))
      }
      if (tolower(method) == "wang-tang") {
        testval = c(seq(-10, (-res), min(0.1, res * 10)), 
                    seq(res, 1, res))
      }
      if (tolower(method) == "zhang") {
        testval = c(seq(-0.05, 0.99, 0.01), seq(1, 5, 
                                                0.1))
      }
      if (hshift == FALSE) {
        fiterr = do.call("rbind", lapply(testval, 
                                         function(p) {
                                           fit = as.data.frame(c(param = p, mae = 0, 
                                                                 rsq = 0, hshift = 0))
                                           test = budyko_sim(fit = fit, method = method, 
                                                             res = res)
                                           AETPest = sapply((data$PET.P), function(z) {
                                             test$AET.P[round(test$PET.P, -log10(res)) == 
                                                          round(z, -log10(res))]
                                           })
                                           return(data.frame(rsq = 1 - summary(stats::lm(AETPest ~ 
                                                                                           data$AET.P))$r.squared, mae = mean(abs(data$AET.P - 
                                                                                                                                    AETPest))))
                                         }))
        wch = which.min(fiterr[, dif])[1]
        fitdat = as.data.frame(c(param = testval[wch], 
                                 mae = round(fiterr$mae[wch], 4), rsq = 1 - 
                                   round(fiterr$rsq[wch], 4), hs = 0))
      }
      else {
        if (is.null(hs)) {
          hs = c(seq(0, 1, 0.1), seq(1, ifelse(mean(data$PET.P, 
                                                    na.rm = TRUE) <= 1, 1, ceiling(mean(data$PET.P, 
                                                                                        na.rm = TRUE))), 0.1))
        }
        else {
          hs = hs
        }
        difM = matrix(data = NA, nrow = length(hs), ncol = length(testval))
        for (r in 1:nrow(difM)) {
          difM[r, ] = sapply(testval, function(p) {
            fit = as.data.frame(c(param = p, err = 0, 
                                  hshift = hs[r]))
            test = budyko_sim(fit = fit, method = method, 
                              res = res, hshift = TRUE)
            difma = mean(abs(data$AET.P - sapply((data$PET.P), 
                                                 function(z) {
                                                   test$AET.P[round(test$PET.P, -log10(res)) == 
                                                                round(z, -log10(res))]
                                                 })))
            return(difma)
          })
        }
        fitdat = as.data.frame(c(param = testval[which(difM == 
                                                         min(difM, na.rm = TRUE), arr.ind = TRUE)[1, 
                                                         ][2]], err = round(min(difM, na.rm = TRUE)[1], 
                                                                            4), hshift = hs[which(difM == min(difM, na.rm = TRUE), 
                                                                                                  arr.ind = TRUE)[1, ][1]]))
      }
    }
    names(fitdat) = method
    return(fitdat)
  }
}

RunningSum<-function (series, windowSize) 
{
  minNThreshold <- 1
  rnmn <- matrix(nrow = length(series), ncol = 1)
  dx <- floor(windowSize/2)
  l <- length(series)
  sm <- 0
  n <- 0
  smm <- c()
  snne <- c()
  spp <- c()
  for (ii in c(1:l)) {
    minindx <- max(ii - dx, 1)
    maxindx <- min(ii + dx, l)
    if (ii == 1) {
      subsrs <- series[minindx:maxindx]
      sm <- sum(subsrs, na.rm = T)
      n <- sum(!is.na(subsrs))
    }
    else {
      if (minindx > 1) {
        sprev <- series[minindx - 1]
        if (!is.na(sprev)) {
          sm <- sm - sprev
          n <- n - 1
        }
      }
      if (maxindx < l) {
        snext <- series[maxindx]
        if (!is.na(snext)) {
          sm <- sm + snext
          n <- n + 1
        }
      }
    }
    if (n > minNThreshold) {
      rnmn[ii] <- sm
    }
    else {
      rnmn[ii] <- NaN
    }
  }
  return(as.vector(rnmn))
}

preprocess_frac <- function(timeStampx, input_var, col_name, dsel, dt) {
  input_df=data.frame(timeStampx ,input_var[[col_name]])
  # Assume the second column needs processing
  #input_df[, 2]<- input_df[, 2] * 4
  input_df[, 2] <- tsEvaNanRunningMean(input_df[, 2], dt)
  # Filter the data frame based on the condition on dsel
  filtered_df <- input_df[(dsel == 12 | dsel == 13), ]
  # Return the processed result
  return(filtered_df)
}

preprocess_in <- function(timeStampX, input_var, col_name, dsel, dt) {
  input_df=data.frame(timeStampX ,input_var[[col_name]])
  # Assume the second column needs processing
  #input_df[, 2]<- input_df[, 2] * 4
  input_df[, 2] <- RunningSum(input_df[, 2], dt)
  # Filter the data frame based on the condition on dsel
  filtered_df <- input_df[(dsel == 12 | dsel == 13), ]
  # Return the processed result
  return(filtered_df)
}

RegimeFast2=function(data){
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

RegimeFunctionM <- function(dvar, cloc, yearstart, wsize,timeStampX,cupa=NA,Q=FALSE) {
  # Convert data into a data frame with timeStampX
  
  #ly <- length(years)
  yr=yearstart
  yrwindow <- yr:(yr + wsize)
  subs <- which(!is.na(match(year(as.Date(timeStampX)), yrwindow)))
  
  data<-dvar[[cloc]][subs]
  data <- data.frame(timeStampX[subs], data)
  # Create a new column 'Qs' by multiplying the second column by 4*nbdays
  dt1=4
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
  Regime <- RegimeFast(data)
  
  # Return the result
  return(Regime)
}

RegimeFunctionMF <- function(dvar, cloc, yearstart, wsize,timeStampX, centres_year1) {
  # Convert data into a data frame with timeStampX
  
  #ly <- length(years)
  yr=yearstart
  yrwindow <- yr:(yr + wsize)
  subs <- which(!is.na(match(year(as.Date(timeStampX)), yrwindow)))
  
  data<-dvar[[cloc]][subs]
  data <- data.frame(timeStampX[subs], data)
  mean(data$data)
  # Create a new column 'Qs' by multiplying the second column by 4*nbdays
  dt1=4
  dt2=29
  dt=dt1*dt2
  data$Qs <- RunningSum(data[, 2], dt)
  names(data)[1]="datetime"
  
  
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
  
  # Return the result
  return(Regime)
}
