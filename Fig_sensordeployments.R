#### figs describing sensor deployment and success
# Fall 2021, Jordan Stark

library(ggplot2)
library(lubridate)
library(tidyr)
library(patchwork)
library(rgdal)
library(raster)
library(rasterVis)

SensorData <- read.csv(paste0(intermediate_path,"model_data.csv"))
SensorData$date <- as_date(SensorData$date)
SensorData <- SensorData[-which(is.na(SensorData$SensorID)),]

SensorData <- SensorData[,!names(SensorData) %in% c("delsurf","deldeep","fracdem_deep","fracdem_surf")]


SensorData_long <- pivot_longer(SensorData,c("vmc_Deep","vmc_Surf","soiltemp"),
                                names_to="sensortype",
                                values_to="value")
SensorData_long <- SensorData_long[!is.na(SensorData_long$value),]



CountByDay <- aggregate(SensorID ~ date + sensortype,SensorData_long,FUN=function(x) length(unique(x)))
names(CountByDay) <- c("day","sensortype","count")


ggplot(CountByDay, aes(x=day,y=count)) +
  geom_area() +
  theme_bw() +
  facet_wrap(~sensortype,ncol=1,strip.position="right")


SensorLength <- aggregate(value ~ sensortype + SiteID,SensorData_long,length,na.action=na.pass)
SensorLength <- pivot_wider(SensorLength,names_from=sensortype,values_from=value)

SiteData <- SensorData[,c("SiteID","elev","log_tci","strdist","totrad","slope","tpi","maxEVI")]
#SiteData <- SiteData[-which(is.na(SiteData$elev)),]
SiteData <- unique(SiteData)

SensorLength <- merge(SensorLength,SiteData)
SensorLength$days <- apply(SensorLength[,c("vmc_Deep","vmc_Surf")],1,max,na.rm=T)


# sample park for background distribution of tci, elev, totrad

parkbound <-readOGR(paste0(gis_path,"GRSM_DATA/GRSM_BOUNDARY_POLYGON"))    
parkbound <- parkbound[parkbound$OBJECTID<18,]

set.seed(1238)

parksample <- spsample(parkbound,10000,"regular")

# import GIS data
elev <- raster(paste0(gis_path,"gsmnp_ascii/elev.txt"))
tci <- raster(paste0(gis_path,"gsmnp_ascii/tci_cor.asc"))
totrad <- raster(paste0(gis_path,"gsmnp_ascii/totrad.txt"))

# set coordinate system
crs(totrad) <- crs(elev) <- crs(tci) <- CRS("+proj=utm +zone=17 +datum=NAD27")

log_tci <- log(tci)   


## extract data across park
parksample$elev <- extract(elev,parksample)
parksample$log_tci <- extract(log_tci,parksample)
parksample$totrad <- extract(totrad,parksample)

# plot!

e.tci <- ggplot(parksample@data, aes(x=elev,y=log_tci)) +
          geom_point(alpha=0.1,size=1) +
          geom_point(data=SensorLength,aes(x=elev,y=log_tci,size=days),color="deeppink3",alpha=0.7) +
          scale_size_binned(breaks=seq(0,700,by=100)) +
          theme_bw() +
          theme(text=element_text(size=11),
                legend.position="bottom") +
          labs(x="elevation (m)",y=("log(TCI)"))

e.rad <- ggplot(parksample@data, aes(x=elev,y=totrad)) +
            geom_point(alpha=0.1,size=1) +
            geom_point(data=SensorLength,aes(x=elev,y=totrad,size=days),color="deeppink3",alpha=0.7) +
            scale_size_binned(breaks=seq(0,700,by=100)) +
            theme_bw() +
            theme(text=element_text(size=11),
                  legend.position="bottom") +
            labs(x="elevation (m)",y=("Annual radiation"))

# figure version for paper
((e.tci + theme(axis.text.x=element_blank(),axis.title.x=element_blank())) / e.rad) + plot_layout(guides="collect") & theme(legend.position="bottom")

# figure version for presentation
(e.tci | e.rad) + plot_layout(guides="collect") & theme(text=element_text(size=20),legend.position="bottom")


#### coverage of deployed dates ####
sensordata <- read.csv(paste0(intermediate_path,"cleaned_sensordata.csv"))
sensordata$Deploy_Date <- ymd(sensordata$Deploy_Date)
sensordata$Remove_Date <- ymd(sensordata$Remove_Date)

# get min and max deployment dates for sites which had multiple sensors
min_deploy <- aggregate(Deploy_Date ~ SiteID,sensordata,min,na.rm=T)
max_remove <- aggregate(Remove_Date ~ SiteID,sensordata,max,na.rm=T)

sensordays <- merge(SensorLength,merge(min_deploy,max_remove,all=T))
sensordays <- sensordays[-which(is.na(sensordays$Remove_Date)),] #take out sensors which were not found

sensordays$max_ndays <- as.numeric(sensordays$Remove_Date - sensordays$Deploy_Date) -3
  # cleaning script removes 3 days (2 at beginning, 1 at end)
sensordays$fracdays <- sensordays$days/sensordays$max_ndays

# total fraction of days measured 
sum(sensordays$days,na.rm=T)/sum(sensordays$max_ndays,na.rm=T)


# test for survival differences along gradients
survtest <- lm(fracdays ~ elev*log_tci*totrad,sensordays)

summary(survtest)

ggplot(sensordays, aes(color=elev)) +
  geom_histogram(aes(x=fracdays),bins=10) +
  theme_classic() +
  labs(x="fraction of days observed",y="number of sensors")


#### map of sensor deployments ####
med_sensordat <- aggregate(cbind(vmc_Surf,vmc_Deep) ~ SiteID + X + Y,
                           sensordata,median,na.rm=T,na.action=na.pass)
spdat <- SpatialPointsDataFrame(coords=med_sensordat[,c("X","Y")],
                                proj4string=CRS("+init=EPSG:32617"),
                                data=med_sensordat)
spdat <- spTransform(spdat,crs(elev))
spdat@data$X.tr <- spdat@coords[,1]
spdat@data$Y.tr <- spdat@coords[,2]

minE <- cellStats(elev,min)
maxE <- cellStats(elev,max)

minVMC <- min(spdat$vmc_Deep,na.rm=T)
maxVMC <- max(spdat$vmc_Deep,na.rm=T)
scalebar_data <- data.frame(x=300000,xend=310000,y=3927000)

dev.off()


fullpark <- gplot(elev,maxpixels=1e8) +
  geom_tile(aes(fill=value),alpha=0.6) +
  geom_point(data=spdat@data,aes(x=X.tr,y=Y.tr,color=vmc_Deep)) +
  theme_void() +
  scale_fill_distiller(name="elevation (m)",
                       na.value="white",
                       palette="Greys",
                       limits=c(minE,maxE)) +
  scale_color_distiller(name="median\n10-15 cm deep\nmoisture (v/v)",
                        na.value="orange",
                        palette="YlGnBu",
                        direction=1,
                        limits=c(minVMC,maxVMC)) +
  geom_segment(data=scalebar_data,aes(x=x,xend=xend,y=y,yend=y),
               arrow=arrow(angle=90,ends="both",length=unit(0.01,"npc")))+
  geom_text(aes(x=mean(c(scalebar_data$x,scalebar_data$xend)),
                y=scalebar_data$y+2000,label="10 km"),size=3) +
  coord_fixed(expand=F) 


ne.ext <- extent(296000,314000,3952000,3964000)
elev.ne <- crop(elev,ne.ext)
spdat.ne <- crop(spdat,ne.ext)

ne.fig <- gplot(elev.ne,maxpixels=1e8) +
  geom_tile(aes(fill=value),alpha=0.6) +
  geom_point(data=spdat.ne@data,aes(x=X.tr,y=Y.tr,color=vmc_Deep),size=2.5,alpha=0.8) +
  theme_void() +
  scale_fill_distiller(name="elevation (m)",
                       na.value="white",
                       palette="Greys",
                       limits=c(minE,maxE)) +
  scale_color_distiller(name="median\n10-15 cm deep\nmoisture (v/v)",
                        na.value="orange",
                        palette="YlGnBu",
                        direction=1,
                        limits=c(minVMC,maxVMC)) +
  coord_fixed(expand=F) +
  theme(legend.position="none")

tiff(filename=paste0(fig_path,"park_points.tif"),
     width=24,height=14,units="cm",res=600,compression="lzw")
print(fullpark)
dev.off()
tiff(filename=paste0(fig_path,"ne_park_points.tif"),
     width=12,height=6,units="cm",res=600,compression="lzw")
print(ne.fig)
dev.off()

