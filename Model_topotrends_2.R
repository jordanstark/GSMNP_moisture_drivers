#### modelling broad topographic trends in sensor data ####
## edited for manuscript, Jordan Stark, Jan 2022

#### setup ####
  ## packages
    library(lubridate)
    library(lme4)
    library(MuMIn)
    library(lmerTest)
    library(car)

  ## logit functions
    logit <- function(x) log(x/(1-x))
    ilogit <- function(x) exp(x)/(1+exp(x))
    
  ## import data
    # full dataset
    fulldf <- read.csv(paste0(intermediate_path,"model_data.csv"))
    fulldf$SiteID <- as.factor(fulldf$SiteID)
    fulldf$date <- ymd(fulldf$date)
   
# #### select sensor readings with matching deployment times ####
  ## select data from year with most sensors deployed
    main_deployment <- fulldf[fulldf$date>"2020-07-10" &
                                fulldf$date<"2021-07-11",]
    sensorlen <- aggregate(cbind(vmc_Surf,vmc_Deep) ~ SiteID + SensorID,
                           data=main_deployment,FUN=function(x)sum(!is.na(x)),
                           na.action=na.pass)
    sensorlen$min_len <- apply(sensorlen[,c("vmc_Surf","vmc_Deep")],
                               1,FUN=function(x) min(x)/365)
    
    good_id <- sensorlen[sensorlen$min_len>0.9,c("SensorID","SiteID")]
    deep_id <- sensorlen[sensorlen$vmc_Deep/365 > 0.9,c("SensorID","SiteID") ]
    
    SelectSensors <- function(data,IDs){
      data$fullid <- paste(data$SensorID,data$SiteID,sep="_")
      IDs$fullid <- paste(IDs$SensorID,IDs$SiteID,sep="_")
      
      data <- data[which(data$fullid %in% IDs$fullid),]
      data$fullid <- NULL
      
      return(data)
    }
    
    deep_all <- SelectSensors(main_deployment,deep_id)
    deep_all[,c("vmc_Surf","delsurf","fracdem_surf")] <- list(NULL)
    names(deep_all)[names(deep_all) == "vmc_Deep"] <- "vmc"

    #write.csv(deep_all,paste0(intermediate_path,"deep_timeseries.csv"),row.names=F)
    
#### effect sizes of drivers on summer data####
  scaledat <- fulldf[fulldf$doy>fulldf$mat &
                         fulldf$doy<fulldf$sen,
                       c("vmc_Deep","SiteID","elev","slope","log_tci",
                         "tpi","prec","vpd","rad","meant","maxEVI","prec_dm1")]
  names(scaledat)[names(scaledat)=="vmc_Deep"] <- "vmc"
    

    
  scaledat$elev <- scale(scaledat$elev)
  scaledat$slope <- scale(scaledat$slope)
  scaledat$log_tci <- scale(scaledat$log_tci)
  scaledat$tpi <- scale(scaledat$tpi)
  scaledat$prec <- scale(scaledat$prec)
  scaledat$vpd <- scale(scaledat$vpd)
  scaledat$rad <- scale(scaledat$rad)
  scaledat$meant <- scale(scaledat$meant)
  scaledat$maxEVI <- scale(sqrt(scaledat$maxEVI))
  
  scaledat$logit_vmc <- logit(scaledat$vmc)
  
  fullmod <- lmer(logit_vmc ~ elev + log_tci + slope + tpi + prec + rad + meant  + (1|SiteID), 
                   scaledat)
  r.squaredGLMM(fullmod)
  
  vif(fullmod)
    
  # save data and model for figures
  save(fullmod,file=paste0(model_out_path,"summer_drivers_lmer_2.RData"))
  write.csv(scaledat,paste0(intermediate_path,"scaled_summer_vmc_drivers_2.csv"))
  
#### summarize data annually ####
  ## variables of interest
    topo_vars <- c("elev","log_tci","slope","tpi","maxEVI")
    daily_vars <- c("prec","vpd","meant")
    
  ## annual medians
    deep_ann <- aggregate(as.formula(paste0("vmc ~ SiteID +", paste(topo_vars,collapse=" + "))),
                           deep_all,median,na.rm=T,na.action=na.pass)
    deep_ann$logit_vmc <- logit(deep_ann$vmc)
    
    write.csv(deep_ann, paste0(intermediate_path,"annual_vmc_drivers.csv"),row..names=F)
    
  ## annual coefficient of variation
    deep_var <- aggregate(as.formula(paste0("vmc ~ SiteID +", paste(topo_vars,collapse=" + "))),
                          deep_all,FUN=function(x) sd(x,na.rm=T)/mean(x,na.rm=T),na.action=na.pass)
    names(deep_var)[names(deep_var)=="vmc"] <- "cv_vmc"
    
  ## scaled effects of drivers on moisture variability
    summary(lm(cv_vmc ~ scale(elev) + scale(log_tci) + scale(slope) + scale(tpi) + scale(maxEVI),deep_var))

  ## effect of elevation on moisture
    elev.mod <- lm(logit_vmc ~ elev, deep_ann)
    
    # soil moisture increases by 0.00011 per masl
    
    ilogit(predict(elev.mod,data.frame(elev=267)))
    # low elevation mean vmc = 0.053
    
    ilogit(predict(elev.mod,data.frame(elev=2025)))
    # high elevation mean vmc = 0.25
  
    0.25/0.053
    # 4.7x wetter at high elevation than low elevation, compared to 1.5x more rainfall!
    
    
  ## full annual model
    annmod <- lm(logit_vmc ~ scale(elev) + scale(log_tci) + 
                   scale(slope) +
                   scale(tpi) + scale(maxEVI),deep_ann,na.action=na.fail)
      # slope has largest effect after elevation

    slope.elev.mod <- lm(logit_vmc ~ elev+slope,deep_ann)
    slope.elev.sc <- lm(logit_vmc ~ scale(elev)+scale(slope),deep_ann)
    
    slope.elev.intxn <-  lm(logit_vmc ~ scale(elev)*scale(slope),deep_ann)
    
    AIC(slope.elev.sc)
    AIC(slope.elev.intxn)
      # adding interaction makes model worse
    
    save(slope.elev.mod,file=paste0(model_out_path,"med_slope_elev.RData"))
    
    
#### topo effect over time fig ####
    moddat <- deep_all
    moddat$logit_vmc <- logit(moddat$vmc)

    
    moddat$elev <- scale(moddat$elev)
    moddat$slope <- scale(moddat$slope)
    moddat$log_tci <- scale(moddat$log_tci)

    deepmod <- lmer(logit_vmc ~ (cos(doy*0.0172)+sin(doy*0.0172))*(elev + slope) + (1|SiteID),moddat)
    r.squaredGLMM(deepmod)
    summary(deepmod)
    vif(deepmod)
    
    save(deepmod,file=paste0(model_out_path,"topo_over_time.RData"))
    