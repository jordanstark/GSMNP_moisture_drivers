#### Script to extract MODIS ET estimates
## Jordan Stark, Fall 2022

#### setup ####
  ## packages
    library(raster)
    library(rgdal)
    library(tidyr)
    library(lubridate)
    library(sp)
    library(rasterVis)
    library(ggplot2)

    
    
  ## import data
    # sensor data
    sensordata <- read.csv(paste0(intermediate_path,"cleaned_sensordata.csv"))
    sensordata$val <- ifelse(sensordata$SRS=="Validation",T,F)
    sensordata$date <- as_date(sensordata$timestamp)
    
    sensorsites <- unique(sensordata[,c("SiteID","X","Y","val")])
    row.names(sensorsites) <- NULL
    
    spsites <- SpatialPointsDataFrame(coords=sensorsites[,c("X","Y")],
                                      proj4string=CRS("+init=EPSG:32617"),
                                      data=sensorsites)
    

    # ET data
    ET <- stack(list.files(paste0(gis_path,"MODIS_ET/"),
                              pattern="MOD16A3GF.006_ET_500m",
                              full.names=T))
    
    
    # QA/QC
    qa <- stack(list.files(paste0(gis_path,"MODIS_ET/"),
                           pattern="MOD16A3GF.006_ET_QC",
                           full.names=T))

    
     
    
    
    
  ## clean stacks
    MaskQA <- function(valuename,stack){

      qamask <- qa

      qamask[qamask > 100 ] <- -1 # values >100 are fill values where ET not calculated

      clean <- mask(stack,qamask,maskvalue=-1)

      return(clean)
    }

    ET <- MaskQA("ET",ET) * 0.1 # data scale factor
                                  # output units are kg/m2/yr
    
#### extract and summarize data ####
  # transform spatial points
    spsites <- spTransform(spsites,crs(ET))
    
  # extract data
    getYear <- function(names){
      substr(simplify2array(strsplit(names,"_doy"))[2,],
             start=0,stop=4)
    }
    
    extractPts <- function(stack,shortname){
      df <- data.frame(raster::extract(stack,spsites))
      names(df) <- paste0(shortname,getYear(names(df)))
      df$SiteID <- spsites$SiteID
      
      longdf <- data.frame(pivot_longer(df,cols=-SiteID,
                                        names_to="year",
                                        names_prefix=shortname,
                                        values_to=shortname))
      
      longdf$year <- as.numeric(longdf$year)
      
      return(longdf)
    }

    ET_df <- extractPts(ET,"ET")
    qual_df <- extractPts(qa,"QA")
 

    mean_ET <- aggregate(ET ~ SiteID, ET_df, mean,na.rm=T,na.action=na.pass)  
  
     
    
    write.csv(mean_ET, paste0(intermediate_path,"site_mean_ET.csv"),row.names=F)
    
    
#### generate meanET raster and make a quick figure
    
    meanET <- mean(ET)
    
    
    parkbound <- readOGR(paste0(gis_path,"GRSM_data/GRSM_BOUNDARY_POLYGON/GRSM_BOUNDARY_POLYGON.shp"))
    parkbound <- parkbound[parkbound$OBJECTID<18,]
    parkbound_tr <- spTransform(parkbound,crs(meanET))

    template <- raster(paste0(pred_path,"deep_vmc_summer_med.tif"))
    
    meanET_park <- crop(meanET,parkbound_tr)
    names(meanET_park) <- "ET"
    
    meanET_tr <- projectRaster(meanET_park,template,method="ngb")
    meanET_tr <- mask(meanET_tr,spTransform(parkbound,crs(template)))
    
    
    p <- gplot(meanET_tr,maxpixels=1e8) +
      geom_tile(aes(fill=value)) +
      theme_void() +
      theme(legend.key.height=unit(0.02,"npc"),
            legend.key.width=unit(0.01,"npc"),
            text=element_text(size=10)) +
      labs(x="",y="",title="Mean Annual Evapotranspiration (2010-2020)") +
      coord_fixed(expand=c(0.001)) +
      scale_fill_distiller(name="ET",
                           palette="RdYlBu",
                           direction=-1,
                           na.value="white") 
    
    
ggsave(paste0(fig_path,"mean_ET.tiff"),plot=p,width=7,height=7)    
   