#### modelling broad topographic trends in sensor data ####
## edited for manuscript, Jordan Stark, Jan 2022
## modified for revised manuscript, JDF Sept-Oct 2022

#setwd("C:\\Users\\fridley\\OneDrive - Clemson University\\Desktop\\SM_analyses")

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
    #fulldf <- read.csv("model_data.csv")
    fulldf$SiteID <- as.factor(fulldf$SiteID)
    fulldf$date <- ymd(fulldf$date)
   
  ## JDF: merge API dataset
    # from script "API_tests.R" on 9-13-22
    api = read.csv(paste0(intermediate_path,"APIdataset.csv"))
    #api = read.csv("APIdataset.csv")
    fulldf = merge(fulldf,api)

  ## JDF: merge ET dataset (from Jordan, email 9-20-22)  
    etdat = read.csv(paste0(intermediate_path,"8D_ET_Sites.csv"))
    #etdat = read.csv("8D_ET_Sites.csv")  
      #create date column from year and DOY
    etdat$date = as.Date(substr(strptime(paste(etdat$year, etdat$DOY), format="%Y %j"),1,10))
    fulldf = merge(fulldf,etdat,all.x=T)
    library(wql) #for interpolation function
    fulldf$ET2 = NULL
    for(i in 1:length(unique(fulldf$SiteID))) {
      true = fulldf$SiteID==unique(fulldf$SiteID)[i]
    fulldf$ET2[true] = interpTs(fulldf$ET[true],gap=NULL) } #linear interpolate to daily from 8D
    
    fulldf$ET2 = fulldf$ET2/8 #correct units
  
    fulldf$logitVMCSurf = logit(fulldf$vmc_Surf)
    fulldf$logitVMCDeep = logit(fulldf$vmc_Deep)
    
        
#################################    
#Model 1: VMC spatial patterns across seasons. We regressed VMC against elevation, a local topographic factor found to significantly moderate elevation effects (slope angle), and a ‘season’ effect that was estimated as a continuous sine-cosine transformation of day-of-year (DOY) with a period of 365 d (Fridley, 2019). We included interactions of elevation and slope with the sine-cosine seasonal effect. VMC was logit-transformed to conform to residual normality assumptions. We modeled 0-5 cm and 10-15 cm depth VMC values in separate models; with those for 0-5 cm depth presented in the supplement.
    
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
    surf_id <- sensorlen[sensorlen$vmc_Surf/365 > 0.9,c("SensorID","SiteID") ]
    
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
    
    surf_all <- SelectSensors(main_deployment,surf_id)
    surf_all[,c("vmc_Deep","deldeep","fracdem_deep")] <- list(NULL)
    names(surf_all)[names(surf_all) == "vmc_Surf"] <- "vmc"

    ##Deep model
        
    moddat <- deep_all
    moddat$logit_vmc <- logit(moddat$vmc)
    moddat$elev <- scale(moddat$elev)
    moddat$slope <- scale(moddat$slope)

    deepmod <- lmer(logit_vmc ~ (cos(doy*0.0172)+sin(doy*0.0172))*(elev + slope) + (1|SiteID),moddat)
    r.squaredGLMM(deepmod)
    summary(deepmod) #this is table for publication
    vif(deepmod)
    anova(deepmod)
    
    ##Shallow model
    
    moddat <- surf_all
    moddat$logit_vmc <- logit(moddat$vmc)
    moddat$elev <- scale(moddat$elev)
    moddat$slope <- scale(moddat$slope)
    moddat$log_tci <- scale(moddat$log_tci)
    
    surfmod <- lmer(logit_vmc ~ (cos(doy*0.0172)+sin(doy*0.0172))*(elev + slope) + (1|SiteID),moddat)
    r.squaredGLMM(surfmod)
    summary(surfmod) #this is table for publication (supplement)
    vif(surfmod)
    anova(surfmod)
    
    #summaries by site: separate for deep and shallow
    deep.range = tapply(deep_all$vmc,deep_all$SiteID,function(x)median(x,na.rm=T))
    surf.range = tapply(surf_all$vmc,surf_all$SiteID,function(x)median(x,na.rm=T))
    range(deep.range,na.rm=T)
    range(surf.range,na.rm=T)
    
#######################################
#Model 2 (a and b). 
    
    ###
    #Model 2a, daily scale (same as first WWR submission, but adding API)
    
    scaledat <- fulldf[fulldf$doy>fulldf$mat &
                         fulldf$doy<fulldf$sen,
                       c("vmc_Deep","vmc_Surf","SiteID","elev","slope","log_tci",
                         "tpi","prec","vpd","rad","meant","maxEVI","prec_dm1","APIdeep","APIsurf","doy")]
    
    
    # save scaling values
    scaledf_2 <- data.frame(var=c("elev","log_tci","slope","tpi","prec","APIdeep","rad","meant"),
                             means=NA,
                             sds=NA)
    
    for(i in 1:length(scaledf_2$var)){
      scaledf_2$means[i] <- mean(scaledat[,scaledf_2$var[i]],na.rm=T)
      scaledf_2$sds[i] <- sd(scaledat[,scaledf_2$var[i]],na.rm=T)
    }
    write.csv(scaledf_2,paste0(intermediate_path,"scaling_values_summer_mod2_Nov2022.csv"),row.names=F)
    

    scaledat$elev <- scale(scaledat$elev)
    scaledat$slope <- scale(scaledat$slope)
    scaledat$log_tci <- scale(scaledat$log_tci)
    scaledat$tpi <- scale(scaledat$tpi)
    scaledat$prec <- scale(scaledat$prec)
    scaledat$vpd <- scale(scaledat$vpd)
    scaledat$rad <- scale(scaledat$rad)
    scaledat$meant <- scale(scaledat$meant)
    scaledat$maxEVI <- scale(sqrt(scaledat$maxEVI))
    scaledat$APIdeep <- scale((scaledat$APIdeep))
    scaledat$APIsurf <- scale((scaledat$APIsurf))
    scaledat$sindoy <- scale(sin(scaledat$doy*0.0172))
    scaledat$cosdoy <- scale(cos(scaledat$doy*0.0172))
    scaledat$prec_dm1 <- scale(scaledat$prec_dm1)
    scaledat$logit_vmc_deep <- logit(scaledat$vmc_Deep)
    scaledat$logit_vmc_surf <- logit(scaledat$vmc_Surf)
    
    #daily, 10-15 cm deep
    fullmodD <- lmer(logit_vmc_deep ~ elev + log_tci + slope + tpi + prec + APIdeep + rad + meant  + (1|SiteID), scaledat)
    summary(fullmodD) #table for publication
    anova(fullmodD)
    r.squaredGLMM(fullmodD)
    vif(fullmodD)
    
    #daily, 0-5 cm deep
    fullmodS <- lmer(logit_vmc_surf ~ elev + log_tci + slope + tpi + prec + APIsurf + rad + meant  + (1|SiteID),  scaledat)
    summary(fullmodS)
    anova(fullmodS)
    r.squaredGLMM(fullmodS)
    vif(fullmodS)
    
    ## save scaled data and model for figures
    save(fullmodD,file=paste0(model_out_path,"summer_drivers_lmer_3_mod2a.RData"))
    write.csv(scaledat,paste0(intermediate_path,"scaled_summer_vmc_drivers_3_mod2a.csv"))
    
    
    
    
    ###
    #Model 2b, 8-d scale
    
    #8-day ET from MODIS is the sum of ET values of the 8 days before the listed DOY
    #could also sum precip for those 8 days, and calculate VMC as the median value over those days
    #model at 8-d resolution would then be:
    #logit(median VMC) ~ elev + log(TCI) + slope + TPI + precip + ET + (1|site)
    #or could also do logit(median VMC) ~ (elev + log(TCI) + slope + TPI) *(precip + ET) + (1|site)
    
    #create 8-d res dataset
    library(RcppRoll)
    fulldf$Pd8 = NULL
    for(i in 1:length(unique(fulldf$SiteID))) {
      true = fulldf$SiteID==unique(fulldf$SiteID)[i]
      pvec = fulldf$prec[true] 
      fulldf$Pd8[true] = roll_sum(pvec,8,fill=NA,align="right") #sum 8 days, including focal day
    }    
    fulldf8 = fulldf[!is.na(fulldf$ET),]
    fulldf8$logit_vmc_deep = logit(fulldf8$vmc_Deep)
    fulldf8$logit_vmc_surf = logit(fulldf8$vmc_Surf)
    
    mod8d <- lmer(logit_vmc_deep ~ scale(elev) + scale(log_tci) + scale(slope) + scale(tpi) + scale(Pd8) + scale(ET) + (1|SiteID), data=fulldf8)
    summary(mod8d)    
    anova(mod8d)
    r.squaredGLMM(mod8d)
    vif(mod8d)
    
    mod8d.s <- lmer(logit_vmc_surf ~ scale(elev) + scale(log_tci) + scale(slope) + scale(tpi) + scale(Pd8) + scale(ET) + (1|SiteID), data=fulldf8)
    summary(mod8d.s)    
    anova(mod8d.s)
    r.squaredGLMM(mod8d.s)
    #basically same model for surface VMC, lower R2, no TPI effect
    
    
#######################################
#Model 3: Predictive model
    
    #validation sensors only
    preddat = fulldf[fulldf$val==F,]
    
    
    # save scaling values
    scaledf_pred <- data.frame(var=c("elev","log_tci","slope","tpi","prec","APIdeep","rad","meant"),
                            means=NA,
                            sds=NA)
    
    for(i in 1:length(scaledf_pred$var)){
      scaledf_pred$means[i] <- mean(preddat[,scaledf_pred$var[i]],na.rm=T)
      scaledf_pred$sds[i] <- sd(preddat[,scaledf_pred$var[i]],na.rm=T)
    }
    write.csv(scaledf_pred,paste0(intermediate_path,"scaling_values_predmod_Nov2022.csv"),row.names=F)
    
    
    
    
    #scale predictors
    preddat$elev = scale(preddat$elev)
    preddat$slope = scale(preddat$slope)
    preddat$tpi = scale(preddat$tpi)
    preddat$APIdeep = scale(preddat$APIdeep)
    preddat$prec = scale(preddat$prec)
    preddat$rad = scale(preddat$rad)
    preddat$meant = scale(preddat$meant)
    preddat$prec_dm1 = scale(preddat$prec_dm1)
    
    #models
    predmod = lmer(logitVMCDeep ~ elev + slope + tpi + APIdeep + (prec*rad*meant) + (1|SiteID), data=preddat)
    predmod2 = lmer(logitVMCDeep ~ elev + slope + tpi + APIdeep + (rad*meant) + (1|SiteID), data=preddat)
    anova(predmod2,predmod)
      #model with pred is much better
    
    predmod3 = lmer(logitVMCDeep ~ elev + APIdeep + (prec*rad*meant) + (1|SiteID), data=preddat)
    anova(predmod,predmod3)
      #model without spatial predictors is better
    
    predmod4 = lmer(logitVMCDeep ~ APIdeep + (elev+slope+tpi)*(prec*rad*meant) + (1|SiteID), data=preddat)
    anova(predmod4,predmod)
      #predmod4 is way way better
    
    predmod5 = lmer(logitVMCDeep ~ (elev+slope+tpi)*(prec*rad*meant*APIdeep) + (1|SiteID), data=preddat)
    anova(predmod4,predmod5)
    
    predmod6 = lmer(logitVMCDeep ~ (elev+slope+tpi)*(prec+rad+meant+APIdeep) + (1|SiteID), data=preddat)
    anova(predmod6,predmod5)
      #predmod5 is the champ - but high potential for overfitting due to 5-way interaction
    
    summary(predmod6)
    r.squaredGLMM(predmod6)
    anova(predmod6)
    
    # save model 6 for figures
    save(predmod6,file=paste0(model_out_path,"predmod6.RData"))
    
    
    #Predict validation sensors
        # get mean and sd for transformations
    scaledatF <- fulldf[fulldf$val==F,
                        c("vmc_Deep","SiteID","elev","slope","log_tci",
                          "tpi","prec","vpd","rad","meant","maxEVI","prec_dm1","doy","APIdeep")]
    names(scaledatF)[names(scaledatF)=="vmc_Deep"] <- "vmc"
    
    scaledatV <- fulldf[fulldf$val==T,
                        c("vmc_Deep","SiteID","elev","slope","log_tci",
                          "tpi","prec","vpd","rad","meant","maxEVI","prec_dm1","doy","APIdeep")]
    names(scaledatV)[names(scaledatV)=="vmc_Deep"] <- "vmc"
    
    scaledatV$elev <- ( scaledatV$elev-mean(scaledatF$elev,na.rm=T) ) / sd(scaledatF$elev,na.rm=T)
    scaledatV$slope <- ( scaledatV$slope-mean(scaledatF$slope,na.rm=T) ) / sd(scaledatF$slope,na.rm=T)
    scaledatV$log_tci <- ( scaledatV$log_tci-mean(scaledatF$log_tci,na.rm=T) ) / sd(scaledatF$log_tci,na.rm=T)
    scaledatV$tpi <- ( scaledatV$tpi-mean(scaledatF$tpi,na.rm=T) ) / sd(scaledatF$tpi,na.rm=T)
    scaledatV$prec <- ( scaledatV$prec-mean(scaledatF$prec,na.rm=T) ) / sd(scaledatF$prec,na.rm=T)
    scaledatV$prec_dm1 <- ( scaledatV$prec_dm1-mean(scaledatF$prec_dm1,na.rm=T) ) / sd(scaledatF$prec_dm1,na.rm=T)
    scaledatV$rad <- ( scaledatV$rad-mean(scaledatF$rad,na.rm=T) ) / sd(scaledatF$rad,na.rm=T)
    scaledatV$meant <- ( scaledatV$meant-mean(scaledatF$meant,na.rm=T) ) / sd(scaledatF$meant,na.rm=T)
    scaledatV$sindoy <- ( sin(scaledatV$doy*0.0172)-mean(sin(scaledatF$doy*0.0172),na.rm=T) ) / sd(sin(scaledatF$doy*0.0172),na.rm=T)
    scaledatV$cosdoy <- ( cos(scaledatV$doy*0.0172)-mean(cos(scaledatF$doy*0.0172),na.rm=T) ) / sd(cos(scaledatF$doy*0.0172),na.rm=T)
    scaledatV$logit_vmc <- logit(scaledatV$vmc)
    scaledatV$APIdeep <- ( scaledatV$APIdeep-mean(scaledatF$APIdeep,na.rm=T) ) / sd(scaledatF$APIdeep,na.rm=T)
    
    #predict with lmer model
    predmod = predmod6 #insert best model
    scaledatV$logit_vmc.pred <- predict(predmod,newdata=scaledatV,re.form=NA)
    
    #back transform predicted values
    scaledatV$vmc.pred = ilogit(scaledatV$logit_vmc.pred)
    
    #examine predicted-observed  
    plot(scaledatV$logit_vmc.pred,scaledatV$logit_vmc)
    plot(scaledatV$vmc.pred,scaledatV$vmc,xlim=c(0,.4),ylim=c(0,.4))
    abline(0,1)  
    
    #model accuracy and bias
    accuracy = abs(scaledatV$vmc.pred-scaledatV$vmc)
    diff2 = accuracy^2; diff2 = diff2[!is.na(diff2)] #squared differences, non-NA values only
    RMSE = sqrt(sum(diff2)/length(diff2)) #root mean square error
    bias = (scaledatV$vmc.pred-scaledatV$vmc)
    plot(scaledatV$vmc,bias); abline(h=0)
    median(accuracy,na.rm=T) #0.027
    mean(accuracy,na.rm=T) #.038
    RMSE #.048
    
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
    

    
    
    