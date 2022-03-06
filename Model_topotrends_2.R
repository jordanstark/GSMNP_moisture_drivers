#### modelling broad topographic trends in sensor data ####
## edited for manuscript, Jordan Stark, Jan 2022, additional edits by JDF

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
  scale_summer_nv <- fulldf[fulldf$doy>fulldf$mat &
                        fulldf$doy<fulldf$sen &
                          fulldf$val==F,
                      c("vmc_Deep","SiteID","elev","slope","log_tci",
                        "tpi","prec","vpd","rad","meant","maxEVI","prec_dm1","doy")]
  names(scale_summer_nv)[names(scale_summer_nv)=="vmc_Deep"] <- "vmc"
  
  #change dataset to full year, and remove validation sensors
  scale_all_nv <- fulldf[fulldf$val==F,
                     c("vmc_Deep","SiteID","elev","slope","log_tci",
                       "tpi","prec","vpd","rad","meant","maxEVI","prec_dm1","doy")]
  names(scale_all_nv)[names(scale_all_nv)=="vmc_Deep"] <- "vmc"
    
  # scale and save parameters
  summernv_scales <- data.frame(var=c("elev","slope","log_tci","tpi","prec",
                                        "rad","meant","maxEVI"),
                                means=NA,
                                sds=NA)
  allnv_scales <- summernv_scales
  
  
  for(i in 1:length(summernv_scales$var)){
    var <- summernv_scales$var[i]
    
    summermean <- mean(scale_summer_nv[,var],na.rm=T)
    summersd <- sd(scale_summer_nv[,var],na.rm=T)
    
    summernv_scales$means[i] <- summermean
    summernv_scales$sds[i] <- summersd
    
    scale_summer_nv[,var] <- (scale_summer_nv[,var] - summermean)/summersd
    
    
    allmean <- mean(scale_all_nv[,var],na.rm=T)
    allsd <- sd(scale_all_nv[,var],na.rm=T)
    
    allnv_scales$means[i] <- allmean
    allnv_scales$sds[i] <- allsd
    
    scale_all_nv[,var] <- (scale_all_nv[,var] - allmean)/allsd
    
    
  }
  
  
  scale_all_nv$sindoy <- scale(sin(scale_all_nv$doy*0.0172))
  scale_all_nv$cosdoy <- scale(cos(scale_all_nv$doy*0.0172))
  
  scale_all_nv$logit_vmc <- logit(scale_all_nv$vmc)
  scale_summer_nv$logit_vmc <- logit(scale_summer_nv$vmc)
  
  fullmod <- lmer(logit_vmc ~ elev + log_tci + slope + tpi + prec + rad + meant  + (1|SiteID), 
                   scale_summer_nv)
  r.squaredGLMM(fullmod)
  
  vif(fullmod)
    
  # save data and model for figures
  save(fullmod,file=paste0(model_out_path,"summer_drivers_lmer_2.RData"))
  write.csv(scale_summer_nv,paste0(intermediate_path,"scaled_summer_vmc_drivers_2.csv"),row.names=F)
  write.csv(summernv_scales,paste0(intermediate_path,"scaling_values_summer_nonval.csv"),row.names=F)
  
  # JDF: model with seasonal effects
  
  fullmod2 <- lmer(logit_vmc ~ (elev + log_tci + slope + tpi + prec + rad + meant)*(sindoy+cosdoy)  + (1|SiteID), 
                   scale_all_nv)
  
  summary(fullmod2)
  r.squaredGLMM(fullmod2)
  #anova(fullmod,fullmod2)
  
  # save data and model for figures
  save(fullmod2,file=paste0(model_out_path,"annual_drivers_lmer_2.RData"))
  write.csv(scale_all_nv,paste0(intermediate_path,"scaled_vmc_drivers_2.csv"),row.names=F)
  write.csv(allnv_scales,paste0(intermediate_path,"scaling_values_nonval.csv"),row.names=F)
  
  
  
  #BRT: not working
  # library(dismo)
  # library(gbm)
  # gbm1 = gbm.step(data=na.omit(scaledat),gbm.x = c(3,4,6,7,9,10,14,15),gbm.y=16,family="gaussian",tree.complexity=2,learning.rate=.01,bag.fraction=.7)
  # summary(gbm1)
  
  ########## Predict validation sensors
  
  # function to scale based on original dataset
  ScaleVar <- function(value, varname, scaledat){
    return((value - scaledat$means[which(scaledat$var==varname)])/scaledat$sds[which(scaledat$var==varname)])
  }
  
  # get mean and sd for transformations
  
  scaledatV <- fulldf[fulldf$val==T,
                     c("vmc_Deep","SiteID","elev","slope","log_tci",
                       "tpi","prec","vpd","rad","meant","maxEVI","prec_dm1","doy")]
  names(scaledatV)[names(scaledatV)=="vmc_Deep"] <- "vmc"
  
  scaledatV$elev <- ScaleVar(scaledatV$elev,"elev",allnv_scales)
  scaledatV$slope <-  ScaleVar(scaledatV$slope,"slope",allnv_scales)
  scaledatV$log_tci <-  ScaleVar(scaledatV$log_tci,"log_tci",allnv_scales)
  scaledatV$tpi <- ScaleVar(scaledatV$tpi,"tpi",allnv_scales)
  scaledatV$prec <- ScaleVar(scaledatV$prec,"prec",allnv_scales)
  scaledatV$rad <- ScaleVar(scaledatV$rad,"rad",allnv_scales)
  scaledatV$meant <- ScaleVar(scaledatV$meant,"meant",allnv_scales)
  
  scaledatV$sindoy <- sin(scaledatV$doy*0.0172)
  scaledatV$cosdoy <- cos(scaledatV$doy*0.0172)
  scaledatV$logit_vmc <- logit(scaledatV$vmc)
  
  #predict with lmer model
  scaledatV$logit_vmc.pred <- predict(fullmod2,newdata=scaledatV,re.form=NA)

  #back transform predicted values
  scaledatV$vmc.pred = ilogit(scaledatV$logit_vmc.pred)
  
  #examine predicted-observed  
  plot(scaledatV$logit_vmc.pred,scaledatV$logit_vmc)
  plot(scaledatV$vmc.pred,scaledatV$vmc,xlim=c(0,.4),ylim=c(0,.4))
  abline(0,1)  
  
  #model bias
  accuracy = abs(scaledatV$vmc.pred-scaledatV$vmc)
  bias = (scaledatV$vmc.pred-scaledatV$vmc)
  plot(scaledatV$vmc,bias); abline(h=0)
  median(accuracy,na.rm=T) #0.03
  
  #validation stats by month
  scaledatV$month = month(as.Date(scaledatV$doy,origin="2019-01-01"))
  par(mfrow=c(3,4),mar=c(3,3,1,1),oma=c(3,3,1,1))
    for(i in 1:12)  {
      plot(scaledatV$vmc[scaledatV$month==i],bias[scaledatV$month==i], 
           col=as.numeric(as.factor(scaledatV$SiteID)),
           xlim=c(0,.4),ylim=c(-.15,.15),main=month.abb[i],
           xlab="VMC",ylab="Absolute Error")
      abline(h=0,col="gray",lty=2)
      text(.35,.1,  round(  median( accuracy[scaledatV$month==i],na.rm=T),3)    )
      text(.35,.13,"MAE=")
    }
  #par(mfrow=c(1,1),oma=c(3,3,1,1),new=T)
  mtext("Observed VMC",side=1,at=-.7,line=3.5,cex=1.8)
  mtext("Absolute Error",side=2,at=.45,line=51,cex=1.8)
  
  #this model seems way too accurate. Let's try a model with only elevation and season:
  fullmod3 <- lmer(logit_vmc ~ (elev)*(sindoy+cosdoy)  + (1|SiteID), scale_all_nv)
  summary(fullmod3)
  r.squaredGLMM(fullmod3) #fixed effect R2 is 23.6
  #predict with lmer model
  scaledatV$logit_vmc.pred.elev <- predict(fullmod3,newdata=scaledatV,re.form=NA)
  #back transform predicted values
  scaledatV$vmc.pred.elev = ilogit(scaledatV$logit_vmc.pred.elev)
  par(mfrow=c(1,1))
  plot(scaledatV$vmc.pred.elev,scaledatV$vmc,col=as.numeric(as.factor(scaledatV$SiteID))) ; abline(0,1)
  accuracy.elev = abs(scaledatV$vmc.pred.elev-scaledatV$vmc)
  bias.elev = (scaledatV$vmc.pred.elev-scaledatV$vmc)
  plot(scaledatV$vmc,bias.elev); abline(h=0)
  median(accuracy.elev,na.rm=T) #0.03
  plot(scaledatV$vmc.pred.elev,scaledatV$vmc.pred,col=as.numeric(as.factor(scaledatV$SiteID))); abline(0,1)
  plot(scaledatV$vmc.pred.elev,scaledatV$vmc,col="red2"); abline(0,1)
  points(scaledatV$vmc.pred,scaledatV$vmc,col="blue"); abline(0,1)
  mean(accuracy,na.rm=T)
  mean(accuracy.elev,na.rm=T) #elevation and season are sufficient to make a good model
  plot(scaledatV$elev,scaledatV$vmc,col=as.numeric(as.factor(scaledatV$SiteID)))
  
  par(mfrow=c(3,4),mar=c(3,3,1,1),oma=c(3,3,1,1))
  for(i in 1:12)  {
    plot(scaledatV$vmc[scaledatV$month==i],bias.elev[scaledatV$month==i], 
         col=as.numeric(as.factor(scaledatV$SiteID)),
         xlim=c(0,.4),ylim=c(-.15,.15),main=month.abb[i],
         xlab="VMC",ylab="Absolute Error")
    abline(h=0,col="gray",lty=2)
    text(.35,.1,  round(  median( accuracy.elev[scaledatV$month==i],na.rm=T),3)    )
    text(.35,.13,"MAE=")
  }
  #par(mfrow=c(1,1),oma=c(3,3,1,1),new=T)
  mtext("Observed VMC",side=1,at=-.7,line=3.5,cex=1.8)
  mtext("Absolute Error",side=2,at=.45,line=51,cex=1.8)
  
  
#### summarize data annually ####
  ## variables of interest
    topo_vars <- c("elev","log_tci","slope","tpi","maxEVI")
    daily_vars <- c("prec","vpd","meant")
    
  ## annual medians
    deep_ann <- aggregate(as.formula(paste0("vmc ~ SiteID +", paste(topo_vars,collapse=" + "))),
                           deep_all,median,na.rm=T,na.action=na.pass)
    deep_ann$logit_vmc <- logit(deep_ann$vmc)
    
    write.csv(deep_ann, paste0(intermediate_path,"annual_vmc_drivers.csv"),row.names=F)
    
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
    