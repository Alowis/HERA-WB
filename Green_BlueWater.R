
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

dsel <- hour(timeStampX)
dt=4
yrlist=c(1951:2020)
Regysave=c()
resave=c()
metrics=c()

###MODIFICATIONS for soil moisture balance
for (cloc in cnames){
  print(cloc)
  Qep=preprocess_in(timeStampX,Qe,cloc,dsel,dt)
  Qy1=Q1_mean_s(Qep,1951)
  data=preprocess_in(timeStampX,Qe,cloc,dsel,dt)
  Rqyi<-RegimeFunctionM(dvar=Qe, cloc, yearstart=1951, wsize=1, timeStampX)
  Regy=Rqyi$date
  center_year1=data.frame(day=Rqyi$date)
  n_years   <- length(yrlist)
  # `res` will hold mean and SI for each year → 2 columns
  res_mat   <- matrix(NA_real_, nrow = n_years, ncol = 2,
                      dimnames = list(NULL, c("mean", "SI")))
  # `Regy_means` will hold the column of `Rqy1$mean` for each year
  # (plus the original date column that you already have)
  Regy_means <- matrix(NA_real_, nrow = length(Regy), ncol = n_years)
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
    
    # mean for this year (using the pre‑computed `data_const`)
    mean_year <- Q1_mean_s(data, year)
    
    #for soil sum and not mean is used
    mean_year <- Q1_sum_s(data, year)
    
    #for soil, bring all values to positive by adding minim value
    Rqy1$mean=Rqy1$mean-min(Rqy1$mean)
    m1=sum(abs(Rqy1$mean))
    #for Soil m1=0
    #m1=0
    # Seasonality index (SI) – note the use of the pre‑computed `m1_const`
    SI <- sum(abs(Rqy1[, 2] - m1/12)) / (m1)
    
    SI <- sum(abs(Rqy1[, 2] - m1/12))/ m1
    
    # if (m1>0) SI=SI
    # if (m1<0) SI=-SI
    print(SI)
    # ---- Store results ------------------------------------------------
    res_mat[i, ] <- c(mean_year, SI)          # two‑column result matrix
    Regy_means[, i] <- Rqy1$mean              # column of means for this year
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # ---- Assemble final objects -----------------------------------------
  # 1) `res` as a data.frame (or tibble) with a `year` column
  res <- data.frame(year = year_vec, res_mat)
  
  # 2) `Regy` – original dates + one column per year of `Rqy1$mean`
  Regy <- cbind(Regy = Regy, Regy_means)
  Regy=data.frame(Regy)
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
  Regy$cloc=rep(cloc,length(Regy[,1]))
  res=data.frame(res)
  o=ecdf(res$mean)
  perc=o(res$mean)
  res$dm=c(0,abs(diff(perc)))
  
  resp1=res[c(1:30),]
  plot=F
  if (plot==T){ 
    
    mean(resp1[,2])
    qh1=quantile(resp1[,2],.75)
    ql1=quantile(resp1[,2],.25)
    qh2=quantile(resp1[,2],.90)
    ql2=quantile(resp1[,2],.10)
    plot(res[,2],type="o",lwd=2,pch=19)
    abline(h=mean(resp1[,2]),lwd=2)
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
  metrics=rbind(metrics,metrix)
  resave<-rbind(resave,res)
  
}





mean(metrics$v)
names(resave)=c("Year","Avail","SI","Variability","Year2","loc")

#save the stuff
Savelist=list(Regimes=Regysave,yearly_res=resave,metrics=metrics)

if (Qv=="bluewater"){
  save(Savelist,file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_v1.Rdata")
}
if (Qv=="greenwater"){
  save(Savelist,file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_green.Rdata")
}


#starting from outputs

load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_v1.Rdata")

load(file="D:/tilloal/Documents/01_Projects/RegimeShifts/Results_BiasVarSeason_green.Rdata")

Savelist$metrics
resave=Savelist$yearly_res
#aggregate by year
Qmean_yAg<- aggregate(resave$Avail, 
                      by = list(year = resave$Year),
                      function(x) c(mean= mean(x,na.rm=T),q1=quantile(x,0.1,na.rm=T),q2=quantile(x,0.9,na.rm=T)))
Qmean_yAg<-do.call(data.frame, Qmean_yAg)

plot(Qmean_yAg$year,Qmean_yAg$x.mean,type="o")


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

#match resave with upstream area

mup=match(as.numeric(resave$loc),UpArea$outlets)
resave$upa=UpArea$upa[mup]
#do a map?

### --------------------------------------------------------------
##  ONE‑STOP CODE: quantiles, centering, binning, counting & plot
## --------------------------------------------------------------



# 1. Input data ---------------------------------------------------

resaveS=resave
resave=resave[-which(resave$Year==1951),]

# 2. Define the 30‑year baseline ---------------------------------
baseline_start <- min(resave$Year, na.rm = TRUE)         # earliest year in the series
baseline_years <- baseline_start:(baseline_start + 29)   # first 30 calendar years

studyVar="Green water availability"
# 3. Quantiles on the baseline (per location) --------------------
loc_q <- resave %>%
  filter(Year %in% baseline_years) %>%                  # reference period only
  group_by(loc) %>%
  summarise(
    med = quantile(Avail, 0.50, na.rm = TRUE),           # 50 % → will become 0
    q10 = quantile(Avail, 0.10, na.rm = TRUE),
    q25 = quantile(Avail, 0.25, na.rm = TRUE),
    q75 = quantile(Avail, 0.75, na.rm = TRUE),
    q90 = quantile(Avail, 0.90, na.rm = TRUE),
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
    d10 = q10 - med,   # negative
    d25 = q25 - med,   # negative
    d75 = q75 - med,   # positive
    d90 = q90 - med    # positive
  )

# 7. Bin each observation into the six categories ---------------
df <- df %>%
  mutate(
    bin = case_when(
      Q_centered > d90                     ~ "above_90",   # most positive
      Q_centered > d75 & Q_centered <= d90 ~ "75_90",
      Q_centered > 0   & Q_centered <= d75 ~ "med_75",
      Q_centered > d25 & Q_centered <= 0   ~ "25_med",
      Q_centered > d10 & Q_centered <= d25 ~ "10_25",
      Q_centered <= d10                    ~ "below_10",
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
  select(Year, bin, pct)

# 9. Give a negative sign to the three “below‑median” bins ----------
# plot_df <- perc_by_year %>%
#   mutate(
#     signed_perc = ifelse(bin %in% c("below_10", "10_25", "25_med"),
#                          -perc,        # negative side
#                          perc)        # positive side
#   )

plot_df <- area_pct_by_year %>%
  mutate(
    signed_upa = ifelse(bin %in% c("below_10", "10_25", "25_med"),
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
  select(fitted) %>%                    # keep only the fitted data
  unnest(fitted) %>%                    # back to a tidy data frame
  rename(trend = .fitted)               # rename the column with the predicted line

trend_stats <- plot_df %>%
  group_by(bin) %>%
  do(broom::tidy(lm(signed_upa ~ Year, data = .))) %>%   # one tidy table per bin
  filter(term == "Year") %>%                             # keep only the slope row
  mutate(
    signif = ifelse(p.value < 0.05, "significant", "ns")
  ) %>%
  select(bin, estimate, std.error, p.value, signif)

trend_fitted <- plot_df %>%
  left_join(trend_stats, by = "bin") %>%          # bring the slope‑info onto the data
  group_by(bin) %>%
  mutate(trend = predict(lm(signed_upa ~ Year, data = cur_data()))) %>%
  ungroup()

# ## 1) Raw percentages for a single year (e.g. 1980)
# raw_1980 <- trend_df %>% filter(Year == 1980)
# 
# ## 2) How far off are they from the “theoretical” values?
# raw_1980 %>%
#   mutate(
#     target = case_when(
#       bin == "below_10" ~ 10,
#       bin == "10_25"    ~ 15,   # (25‑10) = 15 % of the distribution
#       bin == "25_med"   ~ 25,
#       bin == "med_75"   ~ 25,
#       bin == "75_90"    ~ 15,
#       bin == "above_90" ~ 10
#     ),
#     diff = trend - target
#   )


# ----------------------------------------------------------------
# 9. Combine the two pairs of bins
#    (below_10 + 10_25)   &   (above_90 + 75_90)
# ----------------------------------------------------------------
# trend_wide <- trend_fitted %>%
#   select(Year, bin, trend, signif) %>%
#   pivot_wider(names_from = bin, values_from = trend)

trend_wide <- trend_fitted %>%
  # bring *both* the numeric trend and the flag into the wide table
  select(Year, bin, trend, signif) %>%
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
  select(Year,
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
      select(Year,
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
    below_comb = trend_25_med + trend_10_25,
    above_comb = trend_med_75 + trend_75_90,
    
    # ----- significance of the combined trend ----------------------
    # “significant” if *either* component is significant
    signif_below_comb = ifelse(signif_25_med == "significant" |
                                 signif_10_25 == "significant",
                               "significant", "ns"),
    signif_above_comb = ifelse(signif_med_75 == "significant" |
                                 signif_75_90 == "significant",
                               "significant", "ns")
  ) %>%
  # keep only the columns we need for the plot
  select(Year,
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
  select(Year, combo_bin, trend, signif)


# ----------------------------------------------------------------
# 10. **Re‑order the factor for the *stacking* you want**
# ----------------------------------------------------------------
##   <- change the vector below to whatever stack order you need.
##      The first element will be the *bottom* slice.
plot_df$bin <- factor(plot_df$bin,
                      levels = c(
                        "below_10",   # 1️⃣ bottom
                        "10_25",      # 2️⃣
                        "25_med", 
                        "above_90", # 3️⃣
                        "75_90",      # 5️⃣
                        "med_75"   # 4️⃣
                      ))




# 10. RED‑to‑BLUE colour palette ------------------------------------
# A smooth gradient from deep red (most negative) → orange → yellow → light‑blue → deep blue (most positive)
pal_fun <- colorRampPalette(c("#D73027", "#F46D43", "#FEE090",
                              "#FEE090", "#4575B4", "#313695"))
my_palette <- setNames(pal_fun(6),
                       c("below_10", "10_25", "25_med",
                         "med_75", "75_90", "above_90"))

my_palette2 <- c( "med_75"="lightblue", "above_comb"= "darkblue","25_med"="gold", "below_comb"="darkred")

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
                    breaks = c("below_10", "10_25", "25_med",
                               "med_75", "75_90", "above_90"),
                    name = "Category") +
  labs(
    title    = "Distribution of area under different flow conditions",
    subtitle = paste0(studyVar),
    x        = "Year",
    y        = "Distribution of area (%)"
  ) +
  geom_line(data = trend_combined,
            aes(x = as.numeric(Year), y = trend,
                colour = combo_bin, linetype=signif),
            size = 1.2, inherit.aes = FALSE) +
  geom_line(data = trend_plot1,
            aes(x = as.numeric(Year), y = trend,colour=bin,linetype=signif),
            size = 1.2, inherit.aes = FALSE) +
  scale_color_manual(values = my_palette2) +
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
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )



#Now do the plot for different regions


# Map the results ----

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

