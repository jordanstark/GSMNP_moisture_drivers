#### Script to extract average canopy seasonality statistics
## Jordan Stark, Fall 2021

#### setup ####
  ## packages
    library(raster)
    library(rgdal)
    library(tidyr)
    library(lubridate)
    library(sp)
    
    
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
    ET <- stack(list.files(paste0(gis_path,"MODIS_ET/8D/"),
                           pattern="MOD16A2GF.006_ET_500m",
                           full.names=T))
    
    
    # QA/QC
    # qa <- stack(list.files(paste0(gis_path,"MODIS_ET/8D/"),
    #                        pattern="MOD16A2GF.006_ET_QC",
    #                        full.names=T))
    # 
    # 
    # 
    # 
    # 
    # 
    # qa_lookup <- read.csv(paste0(gis_path,
    #                              "MODIS_ET/8D/MOD16A2GF-006-ET-QC-500m-lookup.csv"))
    
    
    ## I did not remove poor QC cells because "For the improved and reprocessed MOD16A2, users may ignore QC data layer
    #because cloud-contaminated LAI/FPAR gaps have been temporally filled before calculating ET
    #(Mu e al., 2007, also see previous section 3.2.1 and following 6.3). QC just denotes if filled
    #LAI/FPAR were used as inputs." (from user manual; https://lpdaac.usgs.gov/documents/494/MOD16_User_Guide_V6.pdf)
    
    
  ## clean stacks
    # MaskQA <- function(valuename,stack){
    # 
    #   #badvals <- qa_lookup$Value[which(qa_lookup[,valuename] !="Best")]
    # 
    #   qamask <- ET
    # 
    #   ET[ET > 32] <- -1
    # 
    #   clean <- mask(stack,qamask,maskvalue=-1)
    # 
    #   return(clean)
    # }

    
    # remove values >32700 (fill) and multiply by 0.1 (scale)
    ET[ET>32700] <- NA
    ET <- ET * 0.1

#### extract and summarize data ####
  # transform spatial points
    spsites <- spTransform(spsites,crs(midup))
    
  # extract data
    getYear <- function(names){
      substr(simplify2array(strsplit(names,"_doy"))[2,],
             start=0,stop=4)
    }
    getDOY <- function(names){
      substr(simplify2array(strsplit(names,"_doy"))[2,],
             start=5,stop=7)
    }
    
    extractPts <- function(stack,shortname){
      df <- data.frame(raster::extract(stack,spsites))
      names(df) <- paste0(shortname,"_",getYear(names(df)),"_",getDOY(names(df)))
      df$SiteID <- spsites$SiteID
      
      longdf <- data.frame(pivot_longer(df,cols=-SiteID,
                                        names_to=c("year","DOY"),
                                        names_sep="_",
                                        names_prefix=paste0(shortname,"_"),
                                        values_to=shortname))
      
      longdf$year <- as.numeric(longdf$year)
      
      return(longdf)
    }
    
    ET_df <- extractPts(ET,"ET")

    

   write.csv(ET_df, paste0(intermediate_path,"8D_ET_Sites.csv"),row.names=F)
    