#### script to model moisture pattern based on summary stats ####
## Jordan Stark, Jan 2021

#### setup ####
  ## packages
    library(lubridate)
    library(ggplot2)
    library(lme4)
  

  ## data import
    # models 
    load(paste0(model_out_path,"prec_freq_mod.Rdata"))
    load(paste0(model_out_path,"prec_amt_mod.Rdata"))
    load(paste0(model_out_path,"dem_mod.Rdata"))
    load(paste0(model_out_path,"drain_mod.Rdata"))

    
    # data scaling
    prec_freq_scale <- read.csv(paste0(model_out_path,"prec_freq_scale.csv"))
    prec_amt_scale <- read.csv(paste0(model_out_path,"amt_scale.csv"))
    dem_scale <- read.csv(paste0(model_out_path,"dem_scale.csv"))
    drain_scale <- read.csv(paste0(model_out_path,"drain_scale.csv"))   

    # site characteristics
    sitechars <- read.csv(paste0(intermediate_path,"site_met_topo.csv"))
    sitechars$date <- ymd(sitechars$date)
    
    # actual moisture data
    rawdat <- read.csv(paste0(intermediate_path,"cleaned_sensordata.csv"))
    rawdat$timestamp <- ymd_hms(rawdat$timestamp)
    
    daily_dat <- rawdat[which(hour(rawdat$timestamp)==0),]
    daily_dat$date <- as_date(daily_dat$timestamp)

    
#### what to model? ####
  daterng <- seq(ymd("2021-02-01"),ymd("2021-05-01"),by="1 day")

  site <- "AT1.1"
  depth <- "Surf"

  # check that the data exist
  ggplot(daily_dat[which(daily_dat$SiteID==site),], aes(x=date)) +
    geom_line(aes(y=vmc_Surf),color="blue") +
    geom_line(aes(y=vmc_Deep),color="black") +
    theme_classic() +
    scale_x_date(limits=c(daterng[1],daterng[length(daterng)]))
  
  
  site_data <- data.frame(SiteID=site,
                          depth=depth,
                          date=daterng,
                          sindoy=sin(yday(daterng)*0.0172),
                          cosdoy=cos(yday(daterng)*0.0172),
                          elev=NA,
                          totrad=NA,
                          slope=NA,
                          rad=NA,
                          meant=NA,
                          vpd=NA,
                          tpi=NA)
  
  vmc_col <- paste0("vmc_",depth)
  

  for(i in 1:length(site_data[,1])){
    for(j in 6:length(site_data)){
      site_data[i,j] <- sitechars[which(sitechars$SiteID==site &
                                          sitechars$date==daterng[i]),
                                  names(sitechars)==names(site_data)[j]]
    }
  }
  
#### function to predict model outputs ####
  
  predict.sitemod <- function(mod, scaledat, sitedata, start_vmc=0.2,prec_start=NA,prec_rng=NA,trans=NA){
    
    sitedata$start_vmc <- start_vmc
    
    if(!is.na(prec_start)) sitedata$prec_start <- prec_start
    if(!is.na(prec_rng)) sitedata$prec_rng <- prec_rng
    
    for(i in 6:length(sitedata)){
      sc_row <- which(scaledat$var==names(sitedata)[i])
      sitedata[,i] <- (sitedata[,i] - scaledat$means[sc_row])/scaledat$sds[sc_row]
    }
    
    sitedata$pred <- predict(mod,sitedata,re.form=NA) # type = response is not working here, so manually back-transforming
    
    if(!is.na(trans)){
      if(trans=="sq"){
        sitedata$pred <- sitedata$pred^2
      } else if(trans=="exp"){
        sitedata$pred <- exp(sitedata$pred)
      } else if(trans=="ilogit"){
        sitedata$pred <- exp(sitedata$pred)/(1+exp(sitedata$pred))
      }
    }
    
    
    return(sitedata$pred)
  }
 
  
#### simulate soil moisture ####  
  start_vmc <- 0.15

  site_data$prob_prec <- NA
  site_data$prec_occ <- NA
  site_data$prec_amt <- NA
  site_data$dem <- NA
  site_data$drain <- NA
  
  site_data$pred_vmc <- NA
  site_data$pred_vmc[1] <- start_vmc
  
  for(i in 1:(length(site_data[,1])-1)) {
    site_data$prob_prec[i] <- predict.sitemod(mod=freq_mod,
                                           scaledat=prec_freq_scale,
                                           sitedata=site_data[i,],
                                           start_vmc=site_data$pred_vmc[i],
                                           trans="ilogit")
    site_data$prec_occ[i] <- rbinom(1,1,site_data$prob_prec[i])
    
    if(site_data$prec_occ[i]==1){
      site_data$prec_amt[i] <- predict.sitemod(mod=amt_mod,
                                               scaledat=prec_amt_scale,
                                               sitedata=site_data[i,],
                                               start_vmc=site_data$pred_vmc[i],
                                               trans="exp")
    } else{
      site_data$prec_amt[i] <- 0
    }
    
    if(i==1){
      site_data$drain[i] <- 0
    } else  if(site_data$prec_occ[i-1]==0){
      site_data$drain[i] <- 0
    } else{
      site_data$drain[i] <- predict.sitemod(mod=drain_mod,
                                            scaledat=drain_scale,
                                            sitedata=site_data[i,],
                                            start_vmc=site_data$pred_vmc[i],
                                            prec_rng=site_data$prec_amt[i-1],
                                            prec_start=site_data$pred_vmc[i-1],
                                            trans="exp") *14 # assuming 24h drainage
    }
    
    site_data$dem[i] <- predict.sitemod(mod=dem_mod,
                                        scaledat=dem_scale,
                                        sitedata=site_data[i,],
                                        start_vmc=site_data$pred_vmc[i],
                                        trans="sq")
    
    
    site_data$pred_vmc[i+1] <- site_data$pred_vmc[i] + site_data$prec_amt[i] - site_data$dem[i] - site_data$drain[i]
  }
  
  
  
    
    ggplot(site_data, aes(x=date)) +
      geom_line(aes(y=pred_vmc),color="red") +
      geom_line(data=daily_dat[which(daily_dat$SiteID==site),], aes_string(y=vmc_col)) +
      theme_classic() +
      scale_x_date(limits=c(min(daterng),max(daterng))) +
      ggtitle(paste0(site, ": ",depth))
        