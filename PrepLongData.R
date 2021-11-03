## Prepare data for models of soil moisture
## Jordan Stark
## Summer 2021

#### setup ####
  # packages
    library(lubridate)
    library(tidyr)
    library(raster)
    library(sp)
    library(rgdal)



  # paths
    data_path <- intermediate_path
    gis_path <- gis_path

#### import and clean data from each source ####
    
  ## sensor data
    # import
      sensordata <- read.csv(paste0(data_path,"cleaned_sensordata.csv"))
      sensordata$val <- ifelse(sensordata$SRS=="Validation",T,F)
      sensordata$date <- as_date(sensordata$timestamp)
      
    # extract only data from midnight
      daily_midn <- sensordata[hour(sensordata$timestamp)==0,
                               c("date","SensorID","SiteID","vmc_Surf","vmc_Deep","soiltemp")]
      
    # calculate daily change in vmc
      daily_midn$delsurf <- NA
      daily_midn$deldeep <- NA
      
      daily_midn <- daily_midn[order(daily_midn$date),]
      
      
      bysite <- split(daily_midn,daily_midn$SiteID)
      
      
      for(j in 1:length(bysite)){
        # fill in NAs
        startday <- min(bysite[[j]]$date,na.rm=T)
        endday <- max(bysite[[j]]$date,na.rm=T)
        all_days <- data.frame(date=seq(startday,endday,by="1 day"))
        
        bysite[[j]] <- merge(bysite[[j]],all_days,all=T)
        
        # calc difference from yesterday
        for(i in 1:(length(bysite[[j]][,1]-1))){
          bysite[[j]]$delsurf[i] <- bysite[[j]]$vmc_Surf[i+1] - bysite[[j]]$vmc_Surf[i]
          bysite[[j]]$deldeep[i] <- bysite[[j]]$vmc_Deep[i+1] - bysite[[j]]$vmc_Deep[i]
          
        } 
      }
      
      daily_vmc <- do.call(rbind,bysite)
      daily_vmc$doy <- yday(daily_vmc$date)
        # for merging later with rad data
      
      # make sites spatial
      sensorsites <- unique(sensordata[,c("SiteID","X","Y","val")])
      row.names(sensorsites) <- NULL
      
      # add correction to get AT1.1 radiation value (it is 1px outside map edge)
      sensorsites[length(sensorsites[,1])+1,] <- sensorsites[sensorsites$SiteID=="AT1.1",]
      sensorsites[length(sensorsites[,1]),"SiteID"] <- "AT1.1rad"
      sensorsites[length(sensorsites[,1]),"X"] <- sensorsites[length(sensorsites[,1]),"X"] - 30
      
      spsites <- SpatialPointsDataFrame(coords=sensorsites[,c("X","Y")],
                                        proj4string=CRS("+init=EPSG:32617"),
                                        data=sensorsites)

  ## function for extracting stacks to sites
      ExtractStack <- function(rstack,varname){
        sitedat <- data.frame(SiteID=spsites$SiteID,
                              raster::extract(rstack,spsites))
        longdat <- pivot_longer(sitedat,-SiteID,names_to="rastername",values_to=varname)
      }
      
  ## landform, radiation, and elevation data
      # import
        elev <- raster(paste0(gis_path,"gsmnp_ascii/elev.txt"))
        tci <- raster(paste0(gis_path,"gsmnp_ascii/tci_cor.asc"))
        rad <- stack(list.files(paste0(gis_path,"gsmnp_ascii/rad/"),
                                full.names=T),quick=T)
        strdist <- raster(paste0(gis_path,"gsmnp_ascii/logsd.txt"))
        totrad <- raster(paste0(gis_path,"gsmnp_ascii/totrad.txt"))
      
      # set coordinate system
        crs(rad) <- crs(strdist) <- crs(totrad) <- crs(elev) <- crs(tci) <- CRS("+proj=utm +zone=17 +datum=NAD27")
        
        log_tci <- log(tci)   
      
      # calculate slope and TPI
        slope <- terrain(elev,opt="slope",unit="degrees")
        tpi <- terrain(elev,opt="tpi")
        
      # import PRISM elev and calculat difference from 30m DEM
        PRISM_elev  <- projectRaster(raster(list.files(paste0(gis_path,"PRISM/Elev/"),
                                                       pattern="*.bil$",full.names=T)),elev,method="ngb")
        delta_elev <- elev - PRISM_elev
          
      # extract at sites
        spsites$elev <- raster::extract(elev,spsites)
        spsites$log_tci <- raster::extract(log_tci,spsites)
        spsites$strdist <- raster::extract(strdist,spsites)
        spsites$totrad <- raster::extract(totrad,spsites)
        spsites$totrad[spsites$SiteID=="AT1.1"] <- spsites$totrad[spsites$SiteID=="AT1.1rad"]
        spsites$slope <- raster::extract(slope,spsites)
        spsites$slope[spsites$SiteID=="AT1.1"] <- spsites$slope[spsites$SiteID=="AT1.1rad"]
        spsites$tpi <- raster::extract(tpi,spsites)
        spsites$tpi[spsites$SiteID=="AT1.1"] <- spsites$tpi[spsites$SiteID=="AT1.1rad"]
        spsites$delta.elev <- raster::extract(delta_elev,spsites)
        
        rad_long <- ExtractStack(rad,"rad")
        rad_long$rastername[rad_long$rastername=="rad073x2"] <- "rad074"
        rad_long$doy <- as.numeric(substring(rad_long$rastername,4))
        rad_long$rastername <- NULL
        
        rad_long$rad[rad_long$SiteID=="AT1.1"] <- rad_long$rad[rad_long$SiteID=="AT1.1rad"]
        
        rad_long$rad[rad_long$rad==0] <- NA
        
      

  ## PRISM VPD (interpolated with elevation)
      vpd_long <- ExtractStack(stack(list.files(paste0(gis_path,"PRISM/daily_vpd"),
                                                full.names=T),quick=T),"vpd")
      vpd_long$date <- ymd(simplify2array(strsplit(vpd_long$rastername,"vpd_"))[2,])
      vpd_long$rastername <- NULL
        
  ## PRISM precipitation (4km)
      prec_long <- ExtractStack(stack(c(list.files(paste0(gis_path,"PRISM/Precip/2019/"),
                                                   pattern="*.bil$",full.names=T),
                                        list.files(paste0(gis_path,"PRISM/Precip/2020/"),
                                                   pattern="*.bil$",full.names=T),
                                        list.files(paste0(gis_path,"PRISM/Precip/2021s/"),
                                                   pattern="*.bil$",full.names=T),
                                        list.files(paste0(gis_path,"PRISM/Precip/2021p/"),
                                                   pattern="*.bil$",full.names=T)),quick=T),
                                "prec")
      prec_long$date <- ymd(simplify2array(strsplit(prec_long$rastername,"_"))[5,])
      prec_long$rastername <- NULL
      
      # add a column for yesterday's precip
      prec_long_dm1 <- prec_long
      prec_long_dm1$date <- prec_long_dm1$date + 1
      names(prec_long_dm1) <- c("SiteID","prec_dm1","date")
      
      prec_long <- merge(prec_long,prec_long_dm1,all=T)
      

      
              
  ## microclimate temperature
      # isolate data from 2019-present
        meant_names <- list.files(paste0(gis_path,"Microclimate/MeanTs/"),
                                 full.names=T)
        meant_names <- meant_names[which(as.numeric(simplify2array(strsplit(meant_names,
                                                                          "_"))[3,]) >= 2019)]

      # extract data
        meant <- ExtractStack(stack(meant_names,quick=T),"meant")
        meant$date <- as_date(meant$rastername,format="y_%Y_d_%j")
        meant$rastername <- NULL
        meant$doy <- yday(meant$date)

  ## seasonality and peak EVI data
    # import EVI data
      midup <- stack(list.files(paste0(gis_path,"Seasonality/"),
                                pattern="MidGreenup_0",
                                full.names=T))
      maturity <- stack(list.files(paste0(gis_path,"Seasonality/"),
                                   pattern="Maturity_0",
                                   full.names=T))
      senescence <- stack(list.files(paste0(gis_path,"Seasonality/"),
                                     pattern="Senescence_0",
                                     full.names=T))
      middown <- stack(list.files(paste0(gis_path,"Seasonality/"),
                                  pattern="MidGreendown_0",
                                  full.names=T))
      minEVI <- stack(list.files(paste0(gis_path,"Seasonality/"),
                               pattern="EVI_Minimum_0",
                               full.names=T)) 
      EVIampl <- stack(list.files(paste0(gis_path,"Seasonality/"),
                                  pattern="EVI_Amplitude_0",
                                  full.names=T))
        
    # import QA/QC
      qa <- stack(list.files(paste0(gis_path,"Seasonality/"),
                             pattern="QA_Detailed_0",
                             full.names=T))
      qa_lookup <- read.csv(paste0(gis_path,
                                   "Seasonality/MCD12Q2-006-QA-Detailed-0-lookup.csv"))
      
        
        
    # clean stacks
      MaskQA <- function(valuename,stack){
        
        badvals <- qa_lookup$Value[which(qa_lookup[,valuename] !="Best")]
        
        qamask <- qa
        
        qamask[qamask %in% badvals ] <- -1
        
        if(!identical(extent(stack),extent(qamask))) stack <- crop(stack,qamask,)
        if(!identical(extent(stack),extent(qamask))) qamask <- crop(qamask,stack)
        
        clean <- mask(stack,qamask,maskvalue=-1)
        
        return(clean)
      }
      
      midup <- MaskQA("MidGreenup",midup)
      maturity <- MaskQA("Maturity",maturity)
      senescence <- MaskQA("Senescence",senescence)
      middown <- MaskQA("MidGreendown",middown)

      
      maxEVI <- minEVI+EVIampl
      
      
    # extract and summarize data
      # transform spatial points
        spsites_evi <- spsites
        spsites_evi <- spTransform(spsites_evi,crs(midup))
        spsites_evi@data[,5:10] <- list(NULL)
        
      # extract data
        getYear <- function(names){
          substr(simplify2array(strsplit(names,"_doy"))[2,],
                 start=0,stop=4)
        }
        
        extractPts <- function(stack,shortname,doy=T){
          df <- data.frame(raster::extract(stack,spsites_evi))
          names(df) <- paste0(shortname,getYear(names(df)))
          df$SiteID <- spsites_evi$SiteID
          
          longdf <- data.frame(pivot_longer(df,cols=-SiteID,
                                            names_to="year",
                                            names_prefix=shortname,
                                            values_to=shortname))
          
          longdf$year <- as.numeric(longdf$year)
          
          if(doy){
          longdf[,shortname] <- yday(as_date(longdf[,shortname],
                                             origin=ymd("1970-01-01")))
          }
          
          return(longdf)
        }
        
        midup_df <- extractPts(midup,"midup")
        mat_df <- extractPts(maturity,"mat")
        sen_df <- extractPts(senescence,"sen")    
        middown_df <- extractPts(middown,"middown")    
        maxEVI_df <- extractPts(maxEVI,"maxEVI",doy=F)
        
        allevi <- merge(midup_df,
                        merge(mat_df,
                              merge(sen_df,
                                    merge(middown_df,maxEVI_df,all=T),all=T),all=T),all=T)   
        
        allevi_summary <- aggregate(cbind(midup,mat,sen,middown,maxEVI) ~ SiteID, allevi, mean,na.rm=T,na.action=na.pass)  
        
        
        
#### merge and finalize datasets ####
  ## merge all data
      alldat <- merge(daily_vmc,
                      merge(spsites@data[,c("SiteID","elev","log_tci","strdist","totrad","slope","tpi","val")],
                            merge(prec_long,
                                  merge(vpd_long,
                                        merge(rad_long,
                                              merge(allevi_summary,meant,all=T),
                                              all=T),
                                        all=T),
                                  all=T),
                            all=T),
                      all=T)

  ## remove dummy site for extracting AT1.1 radiation   
      alldat <- alldat[alldat$SiteID!="AT1.1rad",]

      alldat$fracdem_deep <- -alldat$deldeep/alldat$vmc_Deep
      alldat$fracdem_surf <- -alldat$delsurf/alldat$vmc_Surf
      
      alldat <- alldat[-which(is.na(alldat$SiteID)),]

      
  ## isolate days when soil moisture was decreasing
      demand <- alldat[which(alldat$deldeep < -0.001 & alldat$delsurf < -0.001),]
      demand_frac <- alldat[which(alldat$fracdem_deep > 0.01 &
                                    alldat$fracdem_surf >0.01),]
      
#### save all data ####
  write.csv(alldat,paste0(data_path,"model_data.csv"),row.names=F)
  write.csv(demand,paste0(data_path,"demand_model_data.csv"),row.names=F)
  write.csv(demand_frac,paste0(data_path,"demand_fracmodel_data.csv"),row.names=F)
  
