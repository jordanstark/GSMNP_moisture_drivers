#### run Fridley microclimate model; adjust climate predictions for change in lapse

#### setup ####
  # packages
    library(raster)
    library(rgdal)
    library(lubridate)
    library(parallel)

  # paths
    fridley_GIS <- paste0(gis_path,"gsmnp_ascii/")
    climate_path <- paste0(gis_path,"Microclimate/")  
    
    
  
  # import GIS data
    tci <- log(raster(paste(fridley_GIS,"tci_cor.asc",sep="")))
    totrad <- raster(paste(fridley_GIS,"totrad.txt",sep=""))
    elev <- raster(paste(fridley_GIS,"elev.txt",sep=""))
    strdist <- raster(paste(fridley_GIS,"logsd.txt",sep=""))
    
    radnames <- list.files(paste(fridley_GIS,"rad/",sep=""),full.names=T)
    
    radstack <- stack(radnames,quick=T)

    
    
    crs(radstack) <- crs(elev) <- crs(strdist) <- crs(totrad) <- crs(tci) <- CRS("+proj=utm +zone=17 +datum=NAD27")
    
    
    #strdist <- log(raster(paste0(gis_path,"strdist_new.tif")))
    
    
    tci <- crop(tci,radstack)
    totrad <- crop(totrad,radstack)
    elev <- crop(elev,radstack)  
    strdist <- crop(strdist,radstack)
    
    names(tci) <- "tci"
    names(totrad) <- "totrad"
    names(elev) <- "elev"
    names(strdist) <- "strdist"
        
    #stack <- stack(tci,totrad,elev,strdist,radstack)
    
    
  # import lapse rates
    
    maxt_lapse <- read.csv(paste(climate_path,"maxtemp_synoptic_2011_2021.csv",sep=""))
    mint_lapse <- read.csv(paste(climate_path,"mintemp_synoptic_2011_2021.csv",sep=""))
    
    maxt_lapse$date <- ymd(maxt_lapse$date)
    mint_lapse$date <- ymd(mint_lapse$date)        
    #maxt_lapse <- maxt_lapse[maxt_lapse$date < ymd("2000-01-01") & maxt_lapse$date > ymd("1969-12-31"),]
    #mint_lapse <- mint_lapse[mint_lapse$date < ymd("2000-01-01") & mint_lapse$date > ymd("1969-12-31"),] 
    
    maxt_lapse$day <- yday(maxt_lapse$date)
    mint_lapse$day <- yday(mint_lapse$date)
    
    maxt_lapse$month <- month(maxt_lapse$date)
    mint_lapse$month <- month(mint_lapse$date)
    
    maxt_lapse$year <- year(maxt_lapse$date)
    mint_lapse$year <- year(mint_lapse$date)
    
  # import model coefficients
    maxtmod <- read.csv(paste(climate_path,"micro_mod/maxcoef_out.csv",sep=""),row.names=1)
    mintmod <- read.csv(paste(climate_path,"micro_mod/mincoef_out.csv",sep=""),row.names=1)
             
  # import files created by Micro_annualrasters.R           
  
    min_tmp_names <- list.files(paste(climate_path,"tmpMin/",sep=""),full.names=T)
    yearminstack <- stack(min_tmp_names,quick=T)

    
    max_tmp_names <- list.files(paste(climate_path,"tmpMax/",sep=""),full.names=T)
    yearmaxstack <- stack(max_tmp_names,quick=T)

    
  # function to calculate historical microclimate from files
    
    CalcDailyClim <- function(year) {
    ## historical 
      ## min temp
      min_lapse <- mint_lapse[mint_lapse$year==year,]
      min_lapse <- min_lapse[order(min_lapse$date),]
      ## max temp
      max_lapse <- maxt_lapse[maxt_lapse$year==year,]
      max_lapse <- max_lapse[order(max_lapse$date),]
      
      
      for(j in 1:length(min_lapse[,1])){
        dminlapse <- min_lapse[min_lapse$day == j,]
        
        
        minSYN  <- dminlapse$intercept + (dminlapse$slope * elev)
        
        # dealing with leap years
        stacknum <- ifelse(length(min_lapse[,1])==365,j,
                           ifelse(j<60,j,j-1))
          # this repeats feb 29 as feb 28 and then continues with correct dates
        
        rad <- radstack[[stacknum]]
        
        mint <- yearminstack[[stacknum]] +
          (mintmod["minSYN",1] * minSYN) +
          (mintmod["minSYN:RAD",1] * minSYN * rad) +
          (mintmod["minSYN:ELEV",1] * minSYN * elev) +
          (mintmod["minSYN:LOG.STRDST",1] * minSYN * strdist) +
          (mintmod["minSYN:I(log(TCI))",1] * minSYN * tci)
        
        names(mint) <- paste("y",year,"d",sprintf("%03d",j),sep="_")
        
        writeRaster(mint,
                    paste(climate_path,"minTs/y_",year,"_d_",sprintf("%03d",j),sep=""),
                    format="GTiff")
        
        
        
        dmaxlapse <- max_lapse[max_lapse$day == j,]
        
        
        maxSYN  <- dmaxlapse$intercept + (dmaxlapse$slope * elev)
        
        
        maxt <- yearmaxstack[[stacknum]] +
          (maxtmod["maxSYN",1] * maxSYN) +
          (maxtmod["maxSYN:RAD",1] * maxSYN * rad) +
          (maxtmod["maxSYN:ELEV",1] * maxSYN * elev) +
          (maxtmod["maxSYN:TOTRAD",1] * maxSYN * totrad ) +
          (maxtmod["maxSYN:LOG.STRDST",1] * maxSYN * strdist)
        
        names(maxt) <- paste("y",year,"d",sprintf("%03d",j),sep="_")
        
        writeRaster(maxt,
                    paste(climate_path,"maxTs/y_",year,"_d_",sprintf("%03d",j),sep=""),
                    format="GTiff")
        
        
        meant <- mean(mint,maxt)
        names(meant) <- paste("y",year,"d",sprintf("%03d",j),sep="_")
        
        meant_filled <- focal(meant, w=matrix(1,3,3), fun=mean, NAonly=T,na.rm=T)
        
        writeRaster(meant_filled,
                    paste(climate_path,"meanTs/y_",year,"_d_",sprintf("%03d",j),sep=""),
                    format="GTiff")
        
      }
      
    }
    
      
    clus <- makePSOCKcluster(11)
    clusterExport(cl=clus,varlist=ls())
    clusterEvalQ(cl=clus,library(raster))
    clusterApply(cl=clus,x=2011:2021,fun=CalcDailyClim)
    stopCluster(clus)
    
  # avoid writing GIS variables back to global env
    allvars <- ls()
    
    rmvars <- allvars[!allvars %in% mainvars]
    
    rm(list=rmvars)
    