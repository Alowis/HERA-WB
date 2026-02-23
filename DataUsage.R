##########################################################################################
############   THIS SCRPT IS FOR RUNNING 1 SQUARE OF THE DOMAIN ON THE HPC  ##############
##########################################################################################
#Import functions from TSEVA
source("~/LFRuns_utils/TS-EVA/functions.R")

#Library importation
suppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
suppressWarnings(suppressMessages(library(rnaturalearthdata)))
suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(lubridate)))


#Functions
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
  
  mierda=which(!is.na(outlets))
  outll=outll[mierda,]
  outlets=outlets[which(!is.na(outlets))]
  outll=data.frame(outlets,outll)
  return (outll)
}
disNcopen=function(fname,dir,outloc){
  ncdis=paste0(dir,"/",fname,".nc")
  ncd=nc_open(ncdis)
  name.vb=names(ncd[['var']])
  namev=name.vb[1]
  time <- ncvar_get(ncd,"time")
  lt=length(time)
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncd,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncd,name.lat) 
  lla=length(latdat)
  
  outllplus=matrix(-9999, nrow = lt*length(outloc[,1]), ncol = 5)
  outllplus=as.data.frame(outllplus)
  for ( idc in 1:length(outloc[,1])){
    print(idc)
    idm=1+(idc-1)*lt
    start=c(outloc$idlo[idc],outloc$idla[idc],1)
    count=c(1,1,lt)
    outlets = ncvar_get(ncd,namev,start = start, count= count)
    outlets=as.vector(outlets)
    outid=rep(outloc[idc,1],length(time))
    #lonlatloop=expand.grid(c(1:lla),c(1:lt))
    lon=rep(londat[idc],length(time))
    lat=rep(latdat[idc],length(time))
    outll=data.frame(outlets,outid,lon,lat,time)
    names(outllplus)=names(outll)
    fck=length(c(idm:(idm+lt-1)))
    outllplus[c(idm:(idm+lt-1)),]= outll
  }
  
  
  # outllplus=c()
  # for ( idc in 1:llo){
  #   print(idc)
  #   start=c(idc,1,1)
  #   count=c(1,lla,1)
  #   outlets = ncvar_get(ncd,namev,start = start, count= count)
  #   lna=length(which(!is.na(outlets)))
  #   if (lna>0){
  #     print("data")
  #     exact=(which(!is.na(outlets)))
  #     for (st in exact){
  #       start=c(idc,st,1)
  #       count=c(1,1,lt)
  #       outlets = ncvar_get(ncd,namev,start = start, count= count)
  #       outlets=as.vector(outlets)
  #       #lonlatloop=expand.grid(c(1:lla),c(1:lt))
  #       lon=rep(londat[idc],length(time))
  #       lat=rep(latdat[st],length(time))
  #       outll=data.frame(outlets,lon,lat,time)
  #     }
  #     outllplus=rbind(outllplus,outll)
  #   }
  
  #}
  
  return (outllplus)
}
disNcopenloc=function(fname,dir,outloc,idc){
  ncdis=paste0(dir,"/",fname,".nc")
  ncd=nc_open(ncdis)
  name.vb=names(ncd[['var']])
  namev=name.vb[1]
  time <- ncvar_get(ncd,"time")
  lt=length(time)
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncd,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncd,name.lat) 
  lla=length(latdat)
  
  idm=1+(idc-1)*lt
  start=c(outloc$idlo[idc],outloc$idla[idc],1)
  count=c(1,1,lt)
  outlets = ncvar_get(ncd,namev,start = start, count= count)
  outlets=as.vector(outlets)
  outid=rep(outloc[idc,1],length(time))
  #lonlatloop=expand.grid(c(1:lla),c(1:lt))
  lon=rep(londat[start[1]],length(time))
  lat=rep(latdat[start[2]],length(time))
  outll=data.frame(outlets,outid,lon,lat,time)
  
  return (outll)
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
plotQj <- function(data,catch,UpAloc){
  #~ data <- as.simu(data)
  names(data)=c("date1","Q1","date2","Q2")
  ## Suppression des lignes sans couples Qobs/Qsim
  data <- subset( data, !is.na(Q1) & !is.na(Q2))
  ## Preparation de la sequence de mois
  mois.deb <-  seq(as.Date("1950-01-01"), by="month", length=12)
  ## Vecteur de catÃ©gories
  period1=paste0(format(range(data$date1)[1],"%Y"),"-",format(range(data$date1)[2],"%Y"))
  period2=paste0(format(range(data$date2)[1],"%Y"),"-",format(range(data$date2)[2],"%Y"))
  print(period2)
  n.area=UpAloc
  name.riv=catch
  main = bquote(.(name.riv)~ " Regime (Area="~.(n.area)~ km^2~") | Periods: "~ .(period1) ~" vs" ~ .(period2))
  ## Suppression des lignes sans couples Qobs/Qsim
  ## Preparation de la sequence de mois
  jours <- as.numeric(format(data$date1,"%j"))
  ## Iddinces des lignes de mÃ©me catÃ©gorie
  ind.j <- tapply(seq(length(jours)), jours, c)
  ind.j <- ind.j[-366]
  
  Qc <- data.frame(date=as.numeric(names(ind.j)),
                   obs=sapply(ind.j, function(x) mean(data$Q1[x], na.rm=TRUE)),
                   sim=sapply(ind.j, function(x) mean(data$Q2[x], na.rm=TRUE)),
                   q25o=sapply(ind.j, function(x) quantile(data$Q1[x],0.25, na.rm=TRUE)),
                   q75o=sapply(ind.j, function(x) quantile(data$Q1[x],.75, na.rm=TRUE)),
                   q25s=sapply(ind.j, function(x) quantile(data$Q2[x],0.25, na.rm=TRUE)),
                   q75s=sapply(ind.j, function(x) quantile(data$Q2[x],.75, na.rm=TRUE)))
  
  md=mean(data$Q1)
  qlim=c(0, max(unlist(Qc[,c(5,7)]),2*md))
  
  # qlim <- c(0, max(unlist(Qc[,-1]),unlist(qqo),unlist(qqs)) )
  # titre.date=paste("DÃ©bits moyens interannuels de",
  # format(range(data$date)[1],"%b %Y"),"Ã©",
  # format(range(data$date)[2],"%b %Y"),"- NASH =",
  # round(100 * nash(Qc$sim, Qc$obs),1),"%" )
  
  ## ParamÃ©tres graphiques
  # par(mar=c(2,3,2,1),mgp=c(1.7,0.5,0))
  ## PrÃ©paration du cadre graphique des Qj observÃ©s
  plot(Qc$date, Qc$obs, type="n", axes=FALSE, ylim=qlim,xaxs="i",yaxs="i",
       xlab = NA, ylab = expression(paste("Q (",m^3/s,")")))
  mtext(main,3,font = 2,line = 0.5,cex = .75)
  abline(v = format(mois.deb,"%j"), col="lightgrey", lty=2)
  lines(Qc$date, Qc$obs, col="blue",lwd=2)
  lines(Qc$date, Qc$sim, col="red",lwd=2)
  axis(2)
  axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=.8)
  box()
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q25o,rev(Qc$q75o)),
          col=alpha("lightblue",.3),border="transparent")
  
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q25s,rev(Qc$q75s)),
          col=alpha("indianred",.3),border="transparent")
  ## Trace des erreurs absolues
  polygon(c(0,Qc$date,365),c(0,abs(Qc$sim-Qc$obs),0),
          col=alpha("lightgrey",.4),border="grey")
  ## Trace de la grille mensuelle
  
  ## Trace des courbes de dÃ©bits journaliers simules et observes
  
  ## Definition de la legende
  legend("topleft", leg=c("1990-2020","1951-1981","deviation"),
         lwd=c(2,2,1,2), col=c("red","blue","grey"),
         cex=.8, lty=c(1,1,2,1), bg=alpha("white",.6))
  #return(Qc)
} 
Regime_comp=function(rspace,Nsq,river,hydroDir,catch){
  #load specific outlet file for this square:
  nrspace=rspace[Nsq,]
  outsq=outletopen(hydroDir,outletname,nrspace)
  Idstart=as.numeric(Nsq)*10000
  if (length(outsq$outlets)>0){
    outsq$outlets=seq((Idstart+1),(Idstart+length(outsq$outlets)))
  }
  
  filename=paste0("timeseries/validations/dis_",Nsq,"_1950_2020_cf")
  
  riv=st_as_sf(river, coords = c("lon", "lat"), crs = 4326)
  outloc=st_as_sf(outsq, coords = c("Var1", "Var2"), crs = 4326)
  oula=st_distance(riv,outloc)
  oula=oula/1000
  md=min(oula)
  if (as.numeric(md)<1){
    rloc=which(oula==md)
  }else{
    print("No river found")
    break
  }
  
  riv=as.matrix(round(outsq[rloc,c(2,3)],3))
  outsq$V1r=round(outsq$Var1,3)
  outsq$V2r=round(outsq$Var2,3)
  
  dists=disNcopenloc(filename,hydroDir,outsq,rloc)
  
  rloc2=which(round(outEFAS$Var1,3)==riv[1] & round(outEFAS$Var2,3)==riv[2])
  outR=outEFAS[rloc2,]
  UpAloc=round(upArea[outR$idlo,outR$idla])
  print(paste0("Upstream area: ", UpAloc))
  df.dis=dists 
  timeStamps=unique(as.Date(df.dis$time,origin="1979-01-01"))
  timeStamps=as.POSIXct(timeStamps-1/24)
  txx=timeStamps
  df.dis$timeStamps=timeStamps
  
  names(df.dis)[c(1,2)]=c("dis","outlets")
  ts1=df.dis[c(tb1:tb2),]
  
  #Now the regime
  dts1=ts1[,c(6,1)]
  
  #plot regime allyears:
  dts=df.dis[,c(6,1)]
  
  ts2=df.dis[c(tb3:tb4),]
  
  #Now the regime
  dts2=ts2[,c(6,1)]
  
  data=data.frame(dts1,dts2)
  #30 days running mean for a smoother regime
  data[,2]=tsEvaNanRunningMean(data[,2],120)
  data[,4]=tsEvaNanRunningMean(data[,4],120)
  
  plotQj(data,catch,UpAloc)
  
}
main_path = 'D:/tilloal/Documents/06_Floodrivers/' ### CHANGE THIS PATH
valid_path = paste0(main_path,'DataPaper/')
#Create the upArea df that will be used to match with station data
nca=nc_open(paste0(valid_path, 'GIS/upArea_European_01min.nc'))
name.lon="lon"
name.lat="lat"
nav=names(nca[['var']])
#Band1 is the second variable
t=nca$var[[2]]
name.var=names(nca$var)[2]
tsize<-t$varsize
tdims<-t$ndims
nt1<-tsize[tdims]
lon=ncvar_get(nca,name.lon)
lat=ncvar_get(nca,name.lat)
start <- rep(1,tdims) # begin with start=(1,1,1,...,1)
count <- tsize # begin w/count=(nx,ny,nz,...,nt), reads entire var
#Here I need to extract only some values as this is the full EFAS domain
upArea   = ncvar_get(nca,name.var,start = start, count= count) 

#convert to sqkm
upArea=upArea/1000000



#Step 1: load data
#Squares that I have on this machine
Nsq = 42
outlets="RNetwork"

# haz="flood"
var = "dis"

# workDir = "/BGFS/CLIMEX/tilloal/HydroMeteo/"
# setwd(workDir)
# hydroDir<-"/BGFS/CLIMEX/tilloal/HydroMeteo/Timeseries/dis6_Uncal2"

hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
#workDir<-("D:/tilloal/Documents/06_Floodrivers/dis")
rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]



nrspace=rspace[Nsq,]
outletname="efas_rnet_100km_01min"
nameout="UCRnet"
outEFAS=outletopen(hydroDir,outletname)


unikout=outEFAS$outlets
outEFAS$latlong=paste(round(outEFAS$Var1,4),round(outEFAS$Var2,4),sep=" ")

lagaronne=outEFAS[which(outEFAS$Var1<=-0.7 & outEFAS$Var2<45.4),]
lagaronne=lagaronne[which(lagaronne$Var1>=-0.8),]
lagaronne=lagaronne[which(lagaronne$Var2>45.32),]
outlag=lagaronne$outlets[1]

laloire=outEFAS[which(outEFAS$Var1<=4.1 & outEFAS$Var1>4.09),]
laloire=laloire[which(laloire$Var2>=46.1 & laloire$Var2<46.15),]


#Load the file
#loading the files as netcdf (needs to be checked offline)
Nsq=42
#load specific outlet file for this square:
outletname="efas_rnet_100km_01min"
nameout="UCRnet"
rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]
nrspace=rspace[Nsq,]
outsq=outletopen(hydroDir,outletname,nrspace)
Idstart=as.numeric(Nsq)*10000
if (length(outsq$outlets)>0){
  outsq$outlets=seq((Idstart+1),(Idstart+length(outsq$outlets)))
}


filename=paste0("timeseries/validations/dis_",Nsq,"_1950_2020_cf")

#Extract two times eries: 1951-1981 & 1990-2020

timebound1=c(as.POSIXct(("1951-01-01 00:00:00")))
timebound2=c(as.POSIXct(("1981-12-31 00:00:00")))

tb1=match(timebound1,df.dis$timeStamps)
tb2=match(timebound2,df.dis$timeStamps)

# Selection of pixels to be checked


Ebro=list(Nsq=42,
river=data.frame("lon"=-0.825,"lat"=41.608),
catch="Ebro @ Zaragoza")

Ardeche=list(Nsq=42,
river=data.frame("lon"=4.658,"lat"=44.258),
catch="Ardeche @ Saint-Martin")

Rhone=list(Nsq=52,
river=data.frame("lon"=4.891,"lat"=45.772),
catch="Rhone @ Lyon")

Po=list(Nsq=52,
river=data.frame("lon"=11.60,"lat"=44.89),
catch="Po @ Ferrara")

Danube=list(Nsq=63,
river=data.frame("lon"=16.64,"lat"=48.13),
catch="Danube @ Vienna")

Warta=list(river=data.frame("lon"=16.941,"lat"=52.408),
Nsq=61,
catch="Warta @ Poznan")

Vistula=list(river=data.frame("lon"=21.03,"lat"=52.24),
Nsq=72,
catch="Vistula @ Warsaw")


bug=list(river=data.frame("lon"=25.42,"lat"=57.54),
Nsq=71,
catch="bug @ bugland")
Nsq=bug$Nsq
river=bug$river
nrspace=rspace[Nsq,]
outsq=outletopen(hydroDir,outletname,nrspace)
Idstart=as.numeric(Nsq)*10000
if (length(outsq$outlets)>0){
  outsq$outlets=seq((Idstart+1),(Idstart+length(outsq$outlets)))
}

filename=paste0("timeseries/validations/dis_",Nsq,"_1950_2020_cf")

riv=st_as_sf(river, coords = c("lon", "lat"), crs = 4326)
outloc=st_as_sf(outsq, coords = c("Var1", "Var2"), crs = 4326)
oula=st_distance(riv,outloc)
oula=oula/1000
md=min(oula)
if (as.numeric(md)<1){
  rloc=which(oula==md)
}else{
  print("No river found")
  break
}

riv=as.matrix(round(outsq[rloc,c(2,3)],3))
outsq$V1r=round(outsq$Var1,3)
outsq$V2r=round(outsq$Var2,3)

dists=disNcopenloc(filename,hydroDir,outsq,rloc)
min(dists$outlets)
timeStamps=unique(as.Date(dists$time,origin="1979-01-01"))
timeStamps=as.POSIXct(timeStamps-1/24)
txx=timeStamps
dists$timeStamps=timeStamps
plot(dists$timeStamps[53000:54500],dists$outlets[53000:54500])

# river=data.frame("lon"=1.825,"lat"=50.09)
# Nsq=40
# catch="Somme @ Abbeville"


# river=data.frame("lon"=2.009,"lat"=50.374)
# Nsq=40
# catch="Canche @ Hesdin"

# river=data.frame("lon"=13.405,"lat"=52.5229)
# Nsq=61
# catch="Spree @ Berlin"

Schelde=list(river=data.frame("lon"=3.774,"lat"=51.058),
Nsq=40,
catch="Schelde @ Gent")

Allriver=list(Ebro,Ardeche,Rhone,Po,Danube,Warta,Vistula,Schelde)

Regime_comp(rspace,Nsq,river,hydroDir,catch)

for (ir in 1:length(Allriver)){
  Riviere=Allriver[[ir]]
  Nsq=Riviere$Nsq
  river=Riviere$river
  catch=Riviere$catch
  Cairo::Cairo(
    20, #length
    15, #width
    file = paste("Regimes/",catch, ".png", sep = ""),
    type = "png", #tiff
    bg = "transparent", #white or transparent depending on your requirement 
    dpi = 300,
    units = "cm" #you can change to pixels etc 
  )
  Regime_comp(rspace,Nsq,river,hydroDir,catch)
  dev.off()
}



#Sicily catchments
#Load model
#Also strange catchment in france

Imera=list(river=data.frame("lon"=13.945,"lat"=37.174),
         Nsq=66,
         catch="Imera @ Drasi")

Platani=list(river=data.frame("lon"=13.657,"lat"=37.476),
           Nsq=66,
           catch="Platani @ Passofonduto",
           sname="PASSOFONDUTO")

Simeto=list(river=data.frame("lon"=14.7933,"lat"=37.6585),
           Nsq=66,
           catch="Simeto @ Ponte Maccarrone",
           sname="PONTE MACCARRONE")


Francheville=list(river=data.frame("lon"=4.57,"lat"=47.758),
            Nsq=41,
            catch="Seine @ Francheville",
            sname="H0100010")


VEP=list(river=data.frame("lon"=4.62,"lat"=48.74),
                  Nsq=40,
                  catch="Saulx @ Vitry",
                  sname="H5172010")


Loranca=list(river=data.frame("lon"=-3.09,"lat"=40.457),
         Nsq=32,
         catch="loranca",
         sname="2003003")


#new technique
Roya=list(   catch="loranca",
             sname="2003003")

#Find Nsq with my station
sID=6139140
Roya=kgeout[which(kgeout$V1==sID),]
#load nrspace

rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]
#ido
rs1=rspace[which(rspace$stlon<=Roya$idlo & rspace$endlon>=Roya$idlo),]
rs2=rs1[which(rs1$stlat<=Roya$idla & rs1$endlat>=Roya$idla),]
Nsq=rs2$spID
nrspace=rs2



outall=outletopen(hydroDir,outletname)
outall$idlalo=paste(outall$idlo,outall$idla,sep=" ")
matcho=match(Roya$idlalo,outall$idlalo)
pic=outall[matcho,]
outsq=outletopen(hydroDir,outletname,nrspace)
Idstart=as.numeric(Nsq)*10000
if (length(outsq$outlets)>0){
  outsq$outlets=seq((Idstart+1),(Idstart+length(outsq$outlets)))
}

filename=paste0("timeseries/validations/dis_",Nsq,"_1950_2020_cf")

riv=st_as_sf(river, coords = c("lon", "lat"), crs = 4326)
outloc=st_as_sf(outsq, coords = c("Var1", "Var2"), crs = 4326)
oula=st_distance(riv,outloc)
oula=oula/1000
md=min(oula)
if (as.numeric(md)<1){
  rloc=which(oula==md)
}else{
  print("No river found")
  break
}

riv=as.matrix(round(outsq[rloc,c(2,3)],3))
outsq$V1r=round(outsq$Var1,3)
outsq$V2r=round(outsq$Var2,3)

pic$vv=paste(pic$Var1,pic$Var2,sep=" ")
outsq$vv=paste(outsq$Var1,outsq$Var2,sep=" ")
rloc=match(pic$vv,outsq$vv)
rloc


dists=disNcopenloc(filename,hydroDir,outsq,rloc)
min(dists$outlets)
timeStamps=unique(as.Date(dists$time,origin="1979-01-01"))
timeStamps=as.POSIXct(timeStamps-1/24)
txx=timeStamps
dists$timeStamps=timeStamps
EFAS_flow=dists
EFAS_flow$dQ=tsEvaNanRunningMean(EFAS_flow$outlets,4)
selectix=EFAS_flow$time-floor(EFAS_flow$time)
EFAS_flday=EFAS_flow[which(selectix==0.5),]

#Now load the observations
#What is the ID of the station?

#load efas matches
efas_stations=read.csv("efas_flooddriver_match.csv")
mystat=efas_stations[which(efas_stations$National_Station_Identifier==Loranca$sname),]

#load domis stations
Q_stations <- st_read(paste0(valid_path,'Europe_daily_combined_2022_v2_WGS84_NUTS.shp'))  # SHP with gauge locations
mystat=Q_stations[which(Q_stations$StationID==sID),]

#loop on all years to load the timeserie of that station:
y=1950
s=mystat$StationID
Q_data <- read.csv(paste0('Q_', y, '.csv'), header = F)  # CSVs with observations
Station_data_IDs <- as.vector(t(Q_data[1, ]))
ix <- which(!is.na(match(Station_data_IDs,s)))
Q_obs <- Q_data[-1, ix]
print(Q_obs)
yrlist=c(1951:2020)
for (y in yrlist){
  print(y)
  Q_data <- read.csv(paste0('Q_', y, '.csv'), header = F)  # CSVs with observations
  Station_data_IDs <- as.vector(t(Q_data[1, ]))
  ix <- which(!is.na(match(Station_data_IDs,s)))
  Q_s <- Q_data[-1, ix]
  Q_obs=c(Q_obs,Q_s)
}
plot(Q_obs)

Q_obs=Q_obs[-c(1,2)]

EFAS_flday$Qobs=Q_obs

Q_comp=EFAS_flday[which(!is.na(EFAS_flday$Qobs)),]

b1=2000
b2=3200
plot(Q_comp$timeStamps[b1:b2],Q_comp$Qobs[b1:b2],type="l",col="blue",xlab="Time",ylab="Q (m3/s)",main = paste0("hydrograph @ ",mystat$StationName,". River: ",mystat$River))
points(Q_comp$timeStamps[b1:b2],Q_comp$dQ[b1:b2], type="l",col="red")
legend(Q_comp$timeStamps[900],30, legend=c("simulated", "observed"),
       col=c("red", "blue"),lty=1, cex=0.8)


Cairo::Cairo(
  20, #length
  15, #width
  file = paste("Regimes/",Francheville$catch, "flows.png", sep = ""),
  type = "png", #tiff
  bg = "transparent", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
plot(Q_comp$timeStamps,Q_comp$Qobs,type="l",col="blue",xlab="Time",ylab="Q (m3/s)",main = paste0("hydrograph @ ",mystat$StationName,". River: ",mystat$River))
points(Q_comp$timeStamps,Q_comp$dQ, type="l",col="red")
legend(Q_comp$timeStamps[900],30, legend=c("simulated", "observed"),
       col=c("red", "blue"),lty=1, cex=0.8)
dev.off()

#now the regimes
plotQjm <- function(data,catch,UpAloc){
  #~ data <- as.simu(data)
  names(data)=c("date","Q1","Q2")
  ## Suppression des lignes sans couples Qobs/Qsim
  data <- subset( data, !is.na(Q1) & !is.na(Q2))
  ## Preparation de la sequence de mois
  mois.deb <-  seq(as.Date("1950-01-01"), by="month", length=12)
  ## Vecteur de catÃ©gories
  period1=paste0(format(range(data$date)[1],"%Y"),"-",format(range(data$date)[2],"%Y"))
  #period2=paste0(format(range(data$date2)[1],"%Y"),"-",format(range(data$date2)[2],"%Y"))
  n.area=UpAloc
  name.riv=catch
  main = bquote(.(name.riv)~ " Regime (Area="~.(n.area)~ km^2~") | Periods: "~ .(period1))
  ## Suppression des lignes sans couples Qobs/Qsim
  ## Preparation de la sequence de mois
  jours <- as.numeric(format(data$date,"%j"))
  ## Iddinces des lignes de mÃ©me catÃ©gorie
  ind.j <- tapply(seq(length(jours)), jours, c)
  ind.j <- ind.j[-366]
  
  Qc <- data.frame(date=as.numeric(names(ind.j)),
                   obs=sapply(ind.j, function(x) mean(data$Q1[x], na.rm=TRUE)),
                   sim=sapply(ind.j, function(x) mean(data$Q2[x], na.rm=TRUE)),
                   q25o=sapply(ind.j, function(x) quantile(data$Q1[x],0.25, na.rm=TRUE)),
                   q75o=sapply(ind.j, function(x) quantile(data$Q1[x],.75, na.rm=TRUE)),
                   q25s=sapply(ind.j, function(x) quantile(data$Q2[x],0.25, na.rm=TRUE)),
                   q75s=sapply(ind.j, function(x) quantile(data$Q2[x],.75, na.rm=TRUE)))
  
  md=mean(data$Q1)
  qlim=c(0, max(unlist(Qc[,c(5,7)]),2*md))
  
  # qlim <- c(0, max(unlist(Qc[,-1]),unlist(qqo),unlist(qqs)) )
  # titre.date=paste("DÃ©bits moyens interannuels de",
  # format(range(data$date)[1],"%b %Y"),"Ã©",
  # format(range(data$date)[2],"%b %Y"),"- NASH =",
  # round(100 * nash(Qc$sim, Qc$obs),1),"%" )
  
  ## ParamÃ©tres graphiques
  # par(mar=c(2,3,2,1),mgp=c(1.7,0.5,0))
  ## PrÃ©paration du cadre graphique des Qj observÃ©s
  plot(Qc$date, Qc$obs, type="n", axes=FALSE, ylim=qlim,xaxs="i",yaxs="i",
       xlab = NA, ylab = expression(paste("Q (",m^3/s,")")))
  mtext(main,3,font = 2,line = 0.5,cex = .75)
  abline(v = format(mois.deb,"%j"), col="lightgrey", lty=2)
  lines(Qc$date, Qc$obs, col="blue",lwd=2)
  lines(Qc$date, Qc$sim, col="red",lwd=2)
  axis(2)
  axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=.8)
  box()
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q25o,rev(Qc$q75o)),
          col=alpha("lightblue",.3),border="transparent")
  
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q25s,rev(Qc$q75s)),
          col=alpha("indianred",.3),border="transparent")
  ## Trace des erreurs absolues
  polygon(c(0,Qc$date,365),c(0,abs(Qc$sim-Qc$obs),0),
          col=alpha("lightgrey",.4),border="grey")
  ## Trace de la grille mensuelle
  
  ## Trace des courbes de dÃ©bits journaliers simules et observes
  
  ## Definition de la legende
  legend("topleft", leg=c("simulated","observed","deviation"),
         lwd=c(2,2,1,2), col=c("red","blue","grey"),
         cex=.8, lty=c(1,1,2,1), bg=alpha("white",.6))
  #return(Qc)
} 

data=Q_comp[,c(6,8,7)]
data[,2]=tsEvaNanRunningMean(data[,2],120)
data[,3]=tsEvaNanRunningMean(data[,3],120)
Cairo::Cairo(
  20, #length
  15, #width
  file = paste("Regimes/",Simeto$catch, ".png", sep = ""),
  type = "png", #tiff
  bg = "transparent", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
test1=plotQjm (data,catch="Roya",UpAloc=mystat$DrainingArea.km2.Provider) 
dev.off()

Q_out=Q_comp[,c(6,8,7,3,4)]
names(Q_out)[3] = "Qsim"
write.csv(Q_out,file=paste0("out/Q",mystat$StationName,".csv"))


#Q loranca from Spanish platforms
loranca=read.csv("3003_Loranca_de_Tajuña.csv",sep=";")

Q_comp$timeD=as.Date(Q_comp$timeStamps)
loranca$timeD=as.Date(loranca$fecha,format="%d/%m/%Y")

Qcomp=inner_join(Q_comp,loranca, by="timeD")


plot(Qcomp$timeD,Qcomp$Qobs,type="l",col="blue",xlab="Time",ylab="Q (m3/s)",main = paste0("hydrograph @ ",mystat$StationName,". River: ",mystat$River))
points(Qcomp$timeD,Qcomp$dQ, type="l",col="red")
points(Qcomp$timeD,Qcomp$caudal, type="l",col="green")
legend(Qcomp$timeStamps[900],30, legend=c("simulated", "observed"),
       col=c("red", "blue"),lty=1, cex=0.8)


