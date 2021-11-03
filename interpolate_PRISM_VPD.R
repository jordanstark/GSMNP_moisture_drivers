## script to process raw PRISM files 
# Jordan Stark
# Summer 2021


#### setup ####

  # pathways
  #gis_path <- "E:/Smokies_Moisture/GIS/"
  data_path <- intermediate_path
  f_gis_path <- paste0(gis_path,"gsmnp_ascii/")
  
  out_path <- paste0(gis_path,"PRISM/daily_vpd/")
  if(!dir.exists(out_path)) dir.create(out_path)
  
  # packages
  library(raster)
  library(rgdal)
  library(sp)
  library(lubridate)
  library(tidyr)
  
    
  # import PRISM data
    maxvpd <- stack(c(list.files(paste0(gis_path,"PRISM/Max_VPD/2019/"),
                                 pattern="*.bil$",full.names=T),
                    list.files(paste0(gis_path,"PRISM/Max_VPD/2020/"),
                               pattern="*.bil$",full.names=T),
                    list.files(paste0(gis_path,"PRISM/Max_VPD/2021s/"),
                               pattern="*.bil$",full.names=T),
                    list.files(paste0(gis_path,"PRISM/Max_VPD/2021p/"),
                               pattern="*.bil$",full.names=T)),quick=T)

    
    elev  <- raster(list.files(paste0(gis_path,"PRISM/Elev/"),
                               pattern="*.bil$",full.names=T))
    
    
    
    GRSM <- readOGR(paste0(gis_path,"GRSM_data/GRSM_BOUNDARY_POLYGON"))
    GRSM <- GRSM[GRSM$OBJECTID<18,]
    GRSM <- spTransform(GRSM,crs(elev))
    
    fine_elev <- raster(paste0(f_gis_path,"elev.txt"))
    crs(fine_elev) <- CRS("+proj=utm +zone=17 +datum=NAD27")
    
    # crop raster stacks to park
    # maxvpd <- crop(maxvpd,GRSM,snap="out")
    # elev <- crop(elev,GRSM,snap="out")

    
#### function to extract date from PRISM files ####
    
    ExtractDate <- function(rastername){
      codes <- rastername |>
               strsplit("_") |>
               simplify2array()
      date <- ymd(codes[5,])
      return(date)
    }

    
#### interpolate based on elevation ####

  # randomly sample 10,000 points
    pt_samp <- spsample(GRSM,10000,"regular")

  # extract data at points 
    
    
    pt_elev <- data.frame(raster::extract(elev,pt_samp))
    names(pt_elev) <- "elev"
    pt_elev$ptnum <- 1:length(pt_elev$elev)
    
    daily_lms <- data.frame(date=as_date(NA),
                           intercept=NA,
                           slope=NA,
                           se.slope=NA,
                           r2=NA,
                           data_sd=NA)
    

    pt_data <- data.frame(raster::extract(maxvpd,pt_samp)) |>
                cbind(pt_elev) |>
                pivot_longer(cols=-c(elev,ptnum),
                   values_to="value",names_to="rastername")

    pt_data$date <- ExtractDate(pt_data$rastername)
    
    split_data <- split(pt_data,pt_data$date)
    
    
    for(j in 1:length(split_data)) {
      model <- summary(lm(value~elev,split_data[[j]]))
      daily_lms[j,"date"] <- split_data[[j]]$date[1]
      daily_lms[j,"intercept"] <- model$coefficients["(Intercept)","Estimate"]
      daily_lms[j,"slope"] <- model$coefficients["elev","Estimate"]
      daily_lms[j,"se.slope"] <- model$coefficients["elev","Std. Error"]
      daily_lms[j,"r2"] <- model$r.squared
      daily_lms[j,"data_sd"] <- sd(split_data[[j]]$value,na.rm=T)
      
    }
  
  
    hist(daily_lms$r2)
    plot(r2 ~ slope, daily_lms)
      # the low r squared days are all days with little elevational variation in vpd
    plot(r2 ~ data_sd, daily_lms)
      # they also have little variance in the data overall
      # so the 'line' (essentially mean) should still describe sites well

    
    
#### predict maxvpd based on elevation across GRSM ####

  for(i in 1:length(daily_lms$date)){
    daily_vpd <- fine_elev*daily_lms[i,"slope"] + daily_lms[i,"intercept"]
    writeRaster(daily_vpd,paste0(out_path,"vpd_",daily_lms[i,"date"],".tiff"),overwrite=T)
  }
    
    
