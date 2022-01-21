#### script to model moisture pattern based on summary stats ####
## Jordan Stark, Jan 2021

#### setup ####
  # packages
    library(lubridate)
    library(ggplot2)
    library(lme4)
    library(MuMIn)
    library(patchwork)
    library(tidyr)
    library(lmerTest)
  
  # import data
    all_events <- read.csv(paste0(intermediate_path,"all_events.csv"))
    all_events$date <- ymd(all_events$date)
    all_events$vmc_rate <- all_events$rng_vmc/(all_events$n_hours)
    all_events$sindoy <- sin(yday(all_events$date)*0.0172)
    all_events$cosdoy <- cos(yday(all_events$date)*0.0172)
    all_events$m <- month(all_events$date)
    
    all_events[which(all_events$event=="dem" & all_events$vmc_rate >= 0),c("vmc_rate","rng_vmc")] <- NA
    
    all_events <- all_events[which(!is.na(all_events$rng_vmc)),]

    
    rawdat <- read.csv(paste0(intermediate_path,"cleaned_sensordata.csv"))
    rawdat$timestamp <- ymd_hms(rawdat$timestamp)
    rawdat$m <- month(rawdat$timestamp)
    rawdat$d <- day(rawdat$timestamp)
    rawdat$y <- year(rawdat$timestamp)
    
    

  # create vector describing whether each day had precip
    prec_freq_df <- rawdat[which(hour(rawdat$timestamp)==0),c("SiteID","timestamp","vmc_Deep","vmc_Surf")]
    prec_freq_df <- pivot_longer(prec_freq_df,
                                 cols=3:4,
                                 names_to="depth",
                                 names_prefix="vmc_",
                                 values_to="vmc")
    prec_freq_df <- prec_freq_df[!is.na(prec_freq_df$vmc),]
    names(prec_freq_df) <- c("SiteID","date","depth","vmc_midn")
    
    prec_freq_df$prec_occ <- F

    for(i in 1:length(all_events[,1])){
      if(all_events$event[i]=="prec"){
        prec_freq_df$prec_occ[which(prec_freq_df$SiteID == all_events$SiteID[i] &
                                      prec_freq_df$date == all_events$date[i] &
                                      prec_freq_df$depth == all_events$depth[i])] <- T
      }
    }
    
    prec_freq_df <- merge(prec_freq_df,unique(all_events[,c("SiteID","elev","aws0_150","totrad","maxEVI","midup","middown")]))
    prec_freq_df$doy <- yday(prec_freq_df$date)
    prec_freq_df$sindoy <- sin(prec_freq_df$doy*0.0172)
    prec_freq_df$cosdoy <- cos(prec_freq_df$doy*0.0172)
    prec_freq_df$doy_frac <- prec_freq_df$doy/365
  
#### prep data for models ####
  PrepModels <- function(data,variables){
    data <- data[complete.cases(data[,names(data) %in% variables]),]
    
    scale_df <- data.frame(var=variables,
                           means=NA,
                           sds=NA)
    for(i in 1:length(variables)){
      var_mean <- mean(data[,variables[i]])
      var_sd <- sd(data[,variables[i]])
      
      scale_df$means[i] <- var_mean
      scale_df$sds[i] <- var_sd
      
      data[,variables[i]] <- (data[,variables[i]] - var_mean)/var_sd
      
    }
    
    data$SiteID <- factor(data$SiteID)
    
    return(list(model_data=data,scaled_vars=scale_df))
    
  }
    

  
#### functions ####
  plot_sc_model <- function(model,depth="Surf",colorvar=NA,panelvar=NA,trans=NA,xvar=NA,ann=T){
    xvars <- names(model@frame[2:(length(model@frame)-1)])
    nreps <- ifelse(!is.na(colorvar),
                    ifelse(!is.na(panelvar), 
                           ifelse(!is.na(xvar),27,9),3),
                    ifelse(!is.na(xvar),3,1))
    
    if(ann){
      pred.df <- data.frame(matrix(ncol=length(xvars),
                                   nrow=365*nreps))
      names(pred.df) <- xvars
      
      pred.df$doy <- NA
      
      pred.df$doy <- rep(1:365,each=nreps)
      pred.df$sindoy <- sin(pred.df$doy*0.0172)
      pred.df$cosdoy <- cos(pred.df$doy*0.0172)
      
    } else if(!is.na(xvar)){
      pred.df <- data.frame(matrix(ncol=length(xvars),
                                   nrow=nreps))
      names(pred.df) <- xvars
      
      
    } else stop("xvar must be provided when ann=F")
    
    
   
    pred.df$depth <- depth
    
    
      
    pred.df[,!names(pred.df) %in% c("sindoy","cosdoy","doy","depth")] <- 0
    
    if(!is.na(colorvar)){
      pred.df[,colorvar] <- c(-2,0,2)
    }
    if(!is.na(panelvar)){
      pred.df[,panelvar] <- rep(c(-2,0,2),each=3)
    }
    if(!is.na(xvar)){
      pred.df[,xvar] <- rep(c(-2,0,2),each=nreps/3)
    }
    
    pred.df$pred <- predict(model,pred.df,re.form=NA,type="response")
    

    if(ann){
      p <- ggplot(pred.df, aes_string(x="doy",y="pred",color=colorvar)) +
        geom_point() +
        theme_bw()
    } else if(!is.na(xvar)){
      pred.df[,colorvar] <- factor(pred.df[,colorvar])
      
      p <- ggplot(pred.df, aes_string(x=xvar,y="pred",color=colorvar)) +
        geom_line() +
        theme_bw()
    } else stop("xvar must be provided when ann=F")
   
    
    if(!is.na(panelvar)){
      p <- p + facet_wrap(facets=panelvar,ncol=1) 
    }
    
    print(p)
  }
  
#### precip frequency model ####
    names(prec_freq_df)[names(prec_freq_df)=="vmc_midn"] <- "start_vmc"
    
    freq_dat <- PrepModels(prec_freq_df,c("start_vmc","totrad","elev"))
    freq_dat[[1]]$seas <- factor(ifelse(freq_dat[[1]]$doy > freq_dat[[1]]$midup &
                                          freq_dat[[1]]$doy < freq_dat[[1]]$middown, 
                                        "growing","winter"))
    
    freq_mod <- glmer(prec_occ ~ (elev + totrad)*(sindoy + cosdoy) + depth +
                        start_vmc + I(start_vmc^2) + (1|SiteID),
                      family=binomial,
                      data=freq_dat[[1]])
    
    
#### demand model ####
  dem_dat_raw <- all_events[which(all_events$event=="dem"),]
  dem_dat_raw <- dem_dat_raw[which(dem_dat_raw$rng_vmc > -0.02), ]
  
  dem_dat <- PrepModels(dem_dat_raw, 
                        c("elev","rad","meant","vpd","start_vmc"))

  dem_mod <- lmer(sqrt(-1*rng_vmc) ~ rad*meant + start_vmc + I(start_vmc^2) +
                    depth + (1|SiteID),
                  dem_dat[[1]])

  plot_sc_model(dem_mod,colorvar="meant",panelvar="elev",trans="exp",ann=F,xvar="rad")  
  
  plotdat <- dem_dat[[1]]
  plotdat$pred_dem <- predict(dem_mod,re.form=NA)
  
  plotdat$elev_cat <- cut(plotdat$elev,3)
  plotdat$meant_cat <- cut(plotdat$meant,3)
  
  ggplot(plotdat, aes(x=date,y=pred_dem,color=elev)) +
    geom_line(aes(group=SiteID)) +
    theme_classic() 

#### drainage model ####
  drain_dat_raw <- all_events[which(all_events$event=="drain"),]
  drain_dat_raw <- drain_dat_raw[which(drain_dat_raw$rng_vmc <0),]
  
  drain_dat <- PrepModels(drain_dat_raw,
                          c("elev","slope","tpi","start_vmc","prec_rng","prec_start"))
  
  drain_dat[[1]]$vmc_slope <- drain_dat[[1]]$rng_vmc / drain_dat[[1]]$n_hours

  
  drain_mod <- lmer(log(-1*vmc_slope) ~ elev + slope + tpi + prec_rng + prec_start  + I(prec_start^2) +
                      depth + (1|SiteID),
                    drain_dat[[1]])
  r.squaredGLMM(drain_mod)
  

#### precip amount model ####
  amt_dat_raw <- all_events[which(all_events$event=="prec"),]
  
  amt_dat <- PrepModels(amt_dat_raw,
                        c("elev","slope","start_vmc","totrad"))
  

  
  amt_mod <- lmer(log(rng_vmc) ~ (sindoy + cosdoy)*(elev + totrad) + slope + start_vmc + I(start_vmc^2) + depth + (1|SiteID),
                   amt_dat[[1]])
  
  
  plot_sc_model(amt_mod,colorvar="totrad",trans="exp")
  
  plot_sc_model(amt_mod,colorvar="totrad",panelvar="elev",trans="exp")  
  
  
#### save data ####
  
  # model coefs
  write.csv(summary(freq_mod)$coefficients,paste0(model_out_path,"prec_freq_coefs.csv"))
  write.csv(summary(dem_mod)$coefficients,paste0(model_out_path,"dem_coefs.csv"))
  write.csv(summary(drain_mod)$coefficients,paste0(model_out_path,"drain_coefs.csv"))
  write.csv(summary(amt_mod)$coefficients,paste0(model_out_path,"prec_amt_coefs.csv"))
  
  # scaling
  write.csv(freq_dat[[2]],paste0(model_out_path,"prec_freq_scale.csv"),row.names=F)
  write.csv(dem_dat[[2]],paste0(model_out_path,"dem_scale.csv"),row.names=F)
  write.csv(drain_dat[[2]],paste0(model_out_path,"drain_scale.csv"),row.names=F)
  write.csv(amt_dat[[2]],paste0(model_out_path,"amt_scale.csv"),row.names=F)
  
  # model objects
  save(freq_mod,file=paste0(model_out_path,"prec_freq_mod.Rdata"))
  save(amt_mod,file=paste0(model_out_path,"prec_amt_mod.Rdata"))
  save(dem_mod,file=paste0(model_out_path,"dem_mod.Rdata"))
  save(drain_mod,file=paste0(model_out_path,"drain_mod.Rdata"))
  