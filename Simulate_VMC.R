#### script to model moisture pattern based on summary stats ####
## Jordan Stark, Jan 2021

#### setup ####
  ## packages
    library(lubridate)
    library(ggplot2)
    library(ggrepel)
    library(lme4)
    library(patchwork)
  

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

#### simulation functions ####
    predict.sitemod <- function(mod, scaledat, sitedata, start_vmc=0.2,prec_start=NA,prec_rng=NA,trans=NA){
      
      sitedata$start_vmc <- start_vmc
      
      if(!is.na(prec_start)) sitedata$prec_start <- prec_start
      if(!is.na(prec_rng)) sitedata$prec_rng <- prec_rng
      
      for(j in 6:length(sitedata)){
        sc_row <- which(scaledat$var==names(sitedata)[j])
        sitedata[,j] <- (sitedata[,j] - scaledat$means[sc_row])/scaledat$sds[sc_row]
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
    
    
    
    SimulateSiteVMC <- function(site,depth,minday=NA,maxday=NA,start_vmc){
     vmc_col <- paste0("vmc_",depth)
      
     if(is.na(minday)) minday <- min(daily_dat$date[which(daily_dat$SiteID==site &
                                                            !is.na(daily_dat[,vmc_col]))],na.rm=T) 
     if(is.na(maxday)) maxday <- min(max(daily_dat$date[which(daily_dat$SiteID==site &
                                                           !is.na(daily_dat[,vmc_col]))],na.rm=T),ymd("2021-07-01"))
     
     daterng <- seq(ymd(minday),ymd(maxday),by="1 day")
     
     
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
     
     for(i in 1:length(site_data[,1])){
       for(j in 6:length(site_data)){
         site_data[i,j] <- sitechars[which(sitechars$SiteID==site &
                                             sitechars$date==daterng[i]),
                                     names(sitechars)==names(site_data)[j]]
       }
     }
     
     site_data$prob_prec <- NA
     site_data$prec_occ <- NA
     site_data$prec_amt <- NA
     site_data$dem <- NA
     site_data$drain <- NA
     
     site_data$pred_vmc <- NA
     site_data$pred_vmc[1] <- start_vmc
     
     
     for(i in 1:(length(site_data[,1])-1)) {
       # fill in missing rad values due to leap years
       if(is.na(site_data[i,"rad"])) site_data[i,"rad"] <- site_data[(i-1),"rad"]
       
       
       site_data$prob_prec[i] <- predict.sitemod(mod=freq_mod,
                                                 scaledat=prec_freq_scale,
                                                 sitedata=site_data[i,],
                                                 start_vmc=site_data$pred_vmc[i],
                                                 trans="ilogit")
       site_data$prec_occ[i] <- rbinom(1,1,site_data$prob_prec[i])
       
       if(is.na(site_data$prec_occ[i])) stop(paste0("error: NA precip on line ", i))
       
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
                                               trans="exp") *14 # assuming 14h drainage (median)
       }
       
       if(site_data$prec_amt[i]==0 & site_data$drain[i]==0){
         # demand should occur every day, but we did not take it out in the predictive models so shouldn't be included here
         site_data$dem[i] <- predict.sitemod(mod=dem_mod,
                                             scaledat=dem_scale,
                                             sitedata=site_data[i,],
                                             start_vmc=site_data$pred_vmc[i],
                                             trans="sq")
       } else {
         site_data$dem[i] <- 0
       }
       
       
       
       site_data$pred_vmc[i+1] <- site_data$pred_vmc[i] + site_data$prec_amt[i] - site_data$dem[i] - site_data$drain[i]
     }
     
     return(site_data)
     
    }
    
     

    
#### run and plot model example ####
  ATE02_mod <- SimulateSiteVMC(site="ATE02",depth="Surf",start_vmc=0.2)
    
  vmc_col <- paste0("vmc_",ATE02_mod$depth[1])
  
    
    ggplot(ATE02_mod, aes(x=date)) +
      geom_line(aes(y=pred_vmc),color="red") +
      geom_line(data=daily_dat[which(daily_dat$SiteID=="ATE02"),], aes_string(y=vmc_col)) +
      theme_classic() +
      ggtitle(paste0("ATE02", ": ","Surf"))
    
    
#### choose sites to run simulation ####
    topo_data <- unique(sitechars[,c("SiteID","elev","slope","totrad")])
    
    ggplot(topo_data, aes(x=elev,y=slope,color=totrad,label=SiteID)) +
      geom_point() +
      geom_label_repel() +
      theme_classic()

    # high elev, low slope, high rad
      # R6
    # low elev, low slope, high rad
      # BC01 ##or AT1.5
    # high elev, high slope, high rad
      # MtSt6 (ATE04 doesn't have enough data)
    # low elev, high slope, high rad
      # GM2
    
    # nearest points to each with low rad
      # R6 with MtSt8 or R9
      # BC01/AT1.5 with LM3, GM1 or R3
      # MtSt6 with ATE02 or R4
      # GM2 with AT4 (LM6 doesn't have enough data)
    
    
    sitelist <- c("R6","BC01","MtSt6","GM2","MtSt8","LM3","R4","AT4")

    sim_sites <- data.frame(SiteID=sitelist)    
    
    sim_sites <- merge(sim_sites, topo_data)
    
    sim_sites$slope_cat <- ifelse(sim_sites$slope > 20, "steep","flat")
    sim_sites$elev_cat <- ifelse(sim_sites$elev > 1400, "high","low")    
    sim_sites$totrad_cat <- ifelse(sim_sites$totrad > 2400000, "exposed","shaded")    
    
#### run simulation at all 8 sites ####
    
    nsim <- 50
    
    out_list <- vector(mode="list",length=nsim)
    
    d <- "Surf"
    
    #~ 50 seconds per rep
    
    for(r in 1:nsim){
      
      sim_out <- vector(mode="list", length=length(sim_sites[,1]))
      
      
      for(site in 1:length(sim_out)){
        startval <- sample(daily_dat[which(daily_dat$SiteID==sim_sites$SiteID[[site]]),
                                     paste0("vmc_",d)],1)
        while(is.na(startval)){
          startval <- sample(daily_dat[which(daily_dat$SiteID==sim_sites$SiteID[[site]]),
                                       paste0("vmc_",d)],1)
        }
        
        sim_out[[site]] <- SimulateSiteVMC(site=sim_sites$SiteID[site],
                                           minday=ymd("2020-07-10"),
                                           depth=d,
                                           start_vmc=startval)
        
        sim_out[[site]]$slope_cat <- sim_sites$slope_cat[site]
        sim_out[[site]]$elev_cat <- sim_sites$elev_cat[site]
        sim_out[[site]]$totrad_cat <- sim_sites$totrad_cat[site]
        
        sim_out[[site]] <- merge(sim_out[[site]],daily_dat[,c("SiteID","date",paste0("vmc_",d))],all.x=T)
      }        
      
      out_list[[r]] <- Reduce(rbind,sim_out)
      out_list[[r]]$simnum <- r
    }

    out_df <- Reduce(rbind,out_list)
    
#### plot outputs ####
    minval <- min(c(out_df$pred_vmc,out_df$vmc_Surf),na.rm=T)
    maxval <- max(c(out_df$pred_vmc,out_df$vmc_Surf),na.rm=T)

    low <- ggplot(out_df[out_df$elev_cat=="low",], aes(x=date,y=pred_vmc,group=simnum)) +
      geom_line(aes(y=vmc_Surf),color="black") +
      geom_line(alpha=0.2,color="red") +
      theme_bw() +
      facet_grid(totrad_cat ~ slope_cat) +
      geom_text(data=sim_sites[which(sim_sites$elev_cat=="low"),],
                aes(label=SiteID,x=ymd("2021-06-20"),y=0.39), inherit.aes=F) +
      scale_y_continuous(limits=c(minval,maxval)) +
      ggtitle("low elevation") +
      theme(legend.position="none")
    
    high <- ggplot(out_df[out_df$elev_cat=="high",], aes(x=date,y=pred_vmc,group=simnum)) +
      geom_line(aes(y=vmc_Surf),color="black") +
      geom_line(alpha=0.2,color="red") +
      theme_bw() +
      facet_grid(totrad_cat ~ slope_cat) +
      geom_text(data=sim_sites[which(sim_sites$elev_cat=="high"),],
                aes(label=SiteID,x=ymd("2021-06-20"),y=0.39), inherit.aes=F) +
      scale_y_continuous(limits=c(minval,maxval)) +
      ggtitle("high elevation") +
      theme(legend.position="none")

    high / low   
    
    