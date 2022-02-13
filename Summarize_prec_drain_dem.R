#### algorithms for summarizing precip, drainage, demand by site and time
## Jan 2022
## Jordan Stark

#### setup ####
  # packages
    library(lubridate) # dates and times
    library(raster) # raster processing
    library(sp) # spatial points  
    library(rgdal) # spatial transformations
    library(ggplot2) # figs
    library(segmented)  # breakpoint regression for demand vs drainage
    
    
    
  # paths
    #intermediate_path
    #gis_path

#### import and clean data from each source ####

## sensor data
  # import
    sensordata <- read.csv(paste0(intermediate_path,"cleaned_sensordata.csv"))
    sensordata$val <- ifelse(sensordata$SRS=="Validation",T,F)
    sensordata$timestamp <- ymd_hms(sensordata$timestamp)
    
## site data (topo, met etc.)
    sitedata <- read.csv(paste0(intermediate_path,"site_met_topo.csv"))
    sitedata$date <- ymd(sitedata$date)
    
    
    site_fixed <- unique(sitedata[,c("SiteID","maxEVI","elev","log_tci","strdist",
                              "totrad","slope","tpi","delta.elev","aws0_50","aws0_150")])
    
#### function to identify precip, drainage, and demand in data ####
    
    IDevents <- function(vmc,times,summary=F,dem_raw=F) {
      df <- data.frame(vmc=vmc,
                       smooth=as.numeric(stats::filter(vmc,rep(1,24),sides=1)/24),
                       time=times,
                       day=as_date(times),
                       event=NA)
      df$del <- c(vmc[-1],NA) - vmc
      
    #### stop if no data ####
      if(sum(df$vmc,na.rm=T)==0){
        return(NA)
      }
      
    #### precip ####
      # ID precip if at least 6 consecutive hours all above 90th quantile
        q90 <- quantile(df$del,p=0.9,na.rm=T)
      
        prec_q90 <- logical(length=length(df$del))
        
        for(i in 1:(length(df$del)-6)){
          prec_q90[i] <- ifelse(min(df$del[i:(i+6)])>q90,T,F)
        }
        
        prec_q90[is.na(prec_q90)] <- F


      # ID precip if ever above 98th quantile
        q98 <- quantile(df$del,p=0.98,na.rm=T)
      
        prec <- ifelse(df$del>q98,T,prec_q90)
        prec[is.na(prec)] <- F
      
      
      # fill in gaps of <2 hours between precip events
        runs <- rle(prec)

        runs$values[runs$lengths < 2 &
                      runs$values==F] <- T
        runs <- inverse.rle(runs)
        
      # interrupt event if dramatic decline in soil moisture
        q01 <- quantile(df$del,0.01,na.rm=T)
        runs[which(df$del<q01)] <- F
      
      # assign unique ID to each event
        n_up <- rle(prec)
        n_up$values <- cumsum(n_up$values) * n_up$values
      
        prec_event <- inverse.rle(n_up)
      
        
      # create summary of each ID'd event
        prec.df <- data.frame(event=n_up$values,
                               starttime=ymd_hms(NA),
                               endtime=ymd_hms(NA),
                               start_vmc=NA,
                               end_vmc=NA,
                               valid=T)
      
        prec.df <- prec.df[which(prec.df$event > 0),]
        prec.df$valid[1] <- F
        prec.df$valid[length(prec.df$valid)] <- F
      
        row.names(prec.df) <- NULL
      
      for(i in 1:(length(prec.df$event))){
        
        startindex <- which(prec_event==i)[1]  
        
        if(startindex<5) next # can't test if too close to start
        
        while(startindex > 5 &
              (!is.na(df$del[startindex-1]) &
              df$del[startindex - 1] > 0) |
              (!is.na(sum(df$del[(startindex-1):(startindex-4)])) &
               max(df$del[(startindex-1):(startindex-4)]) > q98 &
               df$del[startindex-1] > q01)){
          # precip started earlier if vmc was already going up, or past vmc went up a lot
          
          startindex <- startindex - 1
        }
        starttime <- df$time[startindex]
        
        
        endindex <- which(prec_event==i)[sum(prec_event==i)]
        
        if(endindex > (length(df$vmc) -5)) next    # can't test if too close to end    

        while(endindex < (length(df$vmc) -5) &
              (!is.na(df$del[endindex+1]) &
               df$del[endindex + 1] > 0) |
              (!is.na(sum(df$del[(endindex+1):(endindex+4)])) &
               max(df$del[(endindex+1):(endindex+4)]) > q98 &
               df$del[endindex] > q01 &
               df$del[endindex+1] > q01)){
          # precip continued later if vmc kept going up, or future vmc went up a lot
          
          endindex <- endindex +1
        }
        
        
        endtime <- df$time[endindex] 
        
        prec.df$starttime[i] <- starttime
        prec.df$endtime[i] <- endtime
        prec.df$start_vmc[i] <- df$vmc[startindex]
        prec.df$end_vmc[i] <- df$vmc[endindex + 1] # this is needed to catch the peak because of how del is calculated
      }
      
      # remove events where total change in vmc is less than q98
        prec.df$rng_vmc <- prec.df$end_vmc - prec.df$start_vmc
        prec.df$valid[prec.df$rng_vmc < q98] <- F
      
      # remove duplicates
        prec.df$event <- NULL
        prec.df <- unique(prec.df)
        
      # list as invalid if not detected
        prec.df$valid[is.na(prec.df$start_vmc)|
                        is.na(prec.df$end_vmc)] <- F

      # list as invalid if partially detected
        for(i in 2:(length(prec.df$valid)-1)){
          if(is.na(prec.df$starttime[i+1]) | is.na(prec.df$endtime[i-1])){
            prec.df$valid[i] <- F
          } else if(prec.df$starttime[i]==prec.df$starttime[i+1]){
            prec.df$valid[i] <- F
          }else if(prec.df$endtime[i]==prec.df$endtime[i-1]){
            prec.df$valid[i] <- F
          }
        }        
        
      # remove invalid events
        prec.df <- prec.df[which(prec.df$valid==T),]
        prec.df$valid <- NULL

      # final listing of precip
        df$prec <- F
        
        for(i in 1:length(prec.df[,1])){
          df$prec[which(df$time>=prec.df$starttime[i] &
                                df$time<=prec.df$endtime[i])] <- T 
    
        }
        
        df$event[which(df$prec==T)] <- "prec"


    #### drainage amd demand ####
    if(length(prec.df[,1])>1){ # don't run if no precip events
      
      drain.df <- data.frame(event=1:length(prec.df[,1]),
                             prec_rng=prec.df$rng_vmc,
                             prec_start=prec.df$start_vmc,
                             starttime=ymd_hms(NA),
                             endtime=ymd_hms(NA),
                             start_vmc=NA,
                             end_vmc=NA,
                             breakpt=NA,
                             slope=NA,
                             r2=NA)
      dem.df.raw <- data.frame(event=1:length(prec.df[,1]),
                           starttime=ymd_hms(NA),
                           endtime=ymd_hms(NA),
                           start_vmc=NA,
                           end_vmc=NA,
                           slope=NA,
                           r2=NA)
      
      df$drain <- NA
      df$dem <- NA

      for(i in 1:(length(drain.df[,1])-1)){
        startindex <- which(df$time==prec.df$endtime[i]) + 1
          # start 1h after precip
        nextstart <- which(df$time==prec.df$starttime[i+1]) - 4
          # end 4h before next ID'd precip event to avoid early increase
        
        
      # isolate vmc data until next rain event
        postvmc <- df$vmc[startindex:nextstart]
        
      # correct start index if max was after the end of ID'd precip
        check <- postvmc[1:6]
        
        startindex <- startindex + which(check == max(check,na.rm=T))[length(which(check == max(check,na.rm=T)))] -1
        
        postvmc <- df$vmc[startindex:nextstart]
      
      # if startindex is the max vmc and at least 20 h of points, continue
        if(!any(is.na(postvmc)) &
           max(postvmc)==postvmc[1] &
           length(postvmc) >= 20){
          
          index <- 1:length(postvmc)
          seg <- tryCatch(segmented(lm(postvmc ~ index)), error=function(e) NA)
            # the next part takes care of errors due to not findng a breakpoint so set this to not stop the loop
          
          if(!is.na(seg)){ # regression worked
          if(is.numeric(summary(seg)$psi[2])) { #breakpoint ID'd
          if(summary(seg)$r.squared > 0.2 & # descriptive fit
             summary(seg)$Ttable[2,1] <0 & #soil moisture declining in drainage event
             summary(seg)$Ttable[2,1] + summary(seg)$Ttable[3,1] < 0) { #soil moisture was falling during demand segment
            
            drain.df$breakpt[i] <- floor(summary(seg)$psi[2])
            drain.df$slope[i] <- summary(seg)$Ttable[2,1]
            drain.df$r2[i] <- summary(seg)$r.squared
            
            drain.df$start_vmc[i] <- postvmc[1]
            drain.df$starttime[i] <- df$time[startindex]
            drain.df$end_vmc[i] <- postvmc[drain.df$breakpt[i]]
            drain.df$endtime[i] <- df$time[startindex+drain.df$breakpt[i]]
            
            dem.df.raw$slope[i] <- drain.df$slope[i] + summary(seg)$Ttable[3,1]
            dem.df.raw$r2[i] <- drain.df$r2[i]
            
            dem.df.raw$starttime[i] <- drain.df$endtime[i] + hours(1)
            dem.df.raw$start_vmc[i] <- postvmc[drain.df$breakpt[i]+1]
            dem.df.raw$endtime[i] <- df$time[nextstart]
            dem.df.raw$end_vmc[i] <- postvmc[length(postvmc)]
             
            # this plots the fit of the segmented regression
            #plot(postvmc); lines(predict(seg),col="red2")
          }}}
          

          
        }
        

      }

     
    # clean drainage summary
     drain.df$rng_vmc <- drain.df$end_vmc - drain.df$start_vmc
     drain.df$event <- NULL
     drain.df <- drain.df[which(!is.na(drain.df$start_vmc) & !is.na(drain.df$end_vmc)),]
     
    # clean demand summary
     # remove events <24h long or where start or end vmc is missing
     dem.df.raw <- dem.df.raw[which(!is.na(dem.df.raw$start_vmc) & !is.na(dem.df.raw$end_vmc)),]
     dem.df.raw <- dem.df.raw[which(as.numeric(difftime(dem.df.raw$endtime,dem.df.raw$starttime,units="hours")) > 24),]
     
     dem.df.raw$rng_vmc <- dem.df.raw$start_vmc - dem.df.raw$end_vmc
     
     # calculate all days within range of demand events
     dem_hours <- ymd_hms(NA)

     for(i in 1:length(dem.df.raw[,1])) {
       dem_hours <- c(dem_hours,seq(dem.df.raw$starttime[i],dem.df.raw$endtime[i],by="1 hour"))
     }
     dem_days <- dem_hours[which(hour(dem_hours)==0)]


     # calculate start and end vmc for each day
     dem.df <- data.frame(starttime=dem_days,
                          start_vmc=NA,
                          endtime=dem_days + days(1),
                          end_vmc=NA,
                          max_smooth24=NA)

     for(i in 1:length(dem.df[,1])){
       startindex <- which(df$time==dem.df$starttime[i])
       dem.df$start_vmc[i] <- df$vmc[startindex]
       
       if(dem.df$endtime[i] %in% dem.df$starttime){
         # only get end vmc if demand was still happening
         endindex <- which(df$time==dem.df$endtime[i])
         dem.df$end_vmc[i] <- df$vmc[endindex]
         dem.df$max_smooth24[i] <- max(df$smooth[startindex:endindex] - df$smooth[(startindex-1):(endindex-1)])
         
       }
     }

     dem.df <- dem.df[which(!is.na(dem.df$end_vmc)),]
     dem.df <- dem.df[which(dem.df$max_smooth24 < 0),]
     
     dem.df$rng_vmc <- dem.df$end_vmc - dem.df$start_vmc
     
     
     
     ## apply calculated times to full timeseries 
     
     for(i in 1:length(drain.df[,1])){
       df$drain[which(df$time>=drain.df$starttime[i] &
                        df$time<=drain.df$endtime[i])] <- T
     }
     for(i in 1:length(dem.df[,1])){
       if(!dem_raw){
         df$dem[which(df$time>=dem.df$starttime[i] &
                        df$time<=dem.df$endtime[i])] <- T
       }
       if(dem_raw){
         df$dem[which(df$time>=dem.df.raw$starttime[i] &
                        df$time<=dem.df.raw$endtime[i])] <- T
       }
     }

     
     # apply events to data
     df$event[which(df$drain==T & is.na(df$event))] <- "drain"
     df$event[which(df$dem==T & is.na(df$event))] <- "dem"
     

    } else { # create empty data frames if no events, to avoid an error when combining

        drain.df <- data.frame(prec_rng=NA,
                               starttime=ymd_hms(NA),
                               endtime=ymd_hms(NA),
                               start_vmc=NA,
                               end_vmc=NA,
                               slope=NA,
                               r2=NA,
                               rng_vmc=NA)
        dem.df <- data.frame(starttime=ymd_hms(NA),
                                endtime=ymd_hms(NA),
                                start_vmc=NA,
                                end_vmc=NA,
                                rng_vmc=NA)
        dem.df.raw <- data.frame(starttime=ymd_hms(NA),
                                 endtime=ymd_hms(NA),
                                 start_vmc=NA,
                                 slope=NA,
                                 end_vmc=NA,
                                 rng_vmc=NA)
        prec.df <- data.frame(starttime=ymd_hms(NA),
                              endtime=ymd_hms(NA),
                              start_vmc=NA,
                              end_vmc=NA,
                              rng_vmc=NA)
      }

    # plot a portion of the data as a diagnostic      
      # ggplot(df[df$time>ymd("2021-02-01") &
      #             df$time<ymd("2021-05-01"),],
      #        aes(x=time,y=vmc,color=event,group=NA)) +
      #   geom_line(size=1) +
      #   theme_classic()

     
      
    #### combine summaries ####
      df$event[which(df$event=="FALSE")] <- NA
      
      dem.df$max_smooth24 <- NULL
        
      dem.df$starttime <- as_datetime(dem.df$starttime)
      dem.df$endtime <- as_datetime(dem.df$endtime)
      
      dem.df$prec_rng <- NA
      dem.df$prec_start <- NA
      
      dem.df.raw$prec_rng <- NA
      dem.df.raw$prec_start <- NA
      
      prec.df$prec_rng <- NA
      prec.df$prec_start <- NA
      prec.df$slope <- NA
      prec.df$r2 <- NA

      


      
      dem.df$event <- "dem"
      dem.df.raw$event <- "dem"
      drain.df$event <- "drain"
      prec.df$event <- "prec"
      
      if(!dem_raw){
        all.df <- rbind(dem.df,prec.df[,names(dem.df)],drain.df[,names(dem.df)])
      }
      if(dem_raw){
        all.df <- rbind(dem.df.raw,prec.df[,names(dem.df.raw)],drain.df[,names(dem.df.raw)])
      }
      
      if(summary==T) {return(list(df$event,all.df))} else {return(df$event)}
    }
    

    
    
    
#### ID precip and calculate hourly change in vmc ####
    
    sensordata <- sensordata[order(sensordata$timestamp),]
    

    
    bysite <- split(sensordata,sensordata$SiteID)
    summary_list <- vector(mode="list",length=length(bysite))
    names(summary_list) <- names(bysite)
    
    
    for(j in 1:length(bysite)){
      # fill in NAs
      start <- as.POSIXct(min(bysite[[j]]$timestamp,na.rm=T))
      end <- as.POSIXct(max(bysite[[j]]$timestamp,na.rm=T))
      all_hours <- data.frame(timestamp=seq(start,end,by="hour"))
      
      bysite[[j]] <- merge(bysite[[j]],all_hours,all=T)
      
      out_Surf <- IDevents(bysite[[j]]$vmc_Surf,bysite[[j]]$timestamp,summary=T,dem_raw=F)
      out_Deep <- IDevents(bysite[[j]]$vmc_Deep,bysite[[j]]$timestamp,summary=T,dem_raw=F)
      
      bysite[[j]]$event_Surf <- out_Surf[[1]]
      bysite[[j]]$event_Deep <- out_Deep[[1]]
      

      
      if(length(out_Surf)>1 & length(out_Deep)>1){
        out_Surf[[2]]$depth <- "Surf"
        out_Deep[[2]]$depth <- "Deep"
        
        summary_list[[j]] <- rbind(out_Surf[[2]],out_Deep[[2]])
        summary_list[[j]]$SiteID <- bysite[[j]]$SiteID[1]
      } else if(length(out_Surf)>1){
        out_Surf[[2]]$depth <- "Surf"
        
        summary_list[[j]] <- out_Surf[[2]]
        summary_list[[j]]$SiteID <- bysite[[j]]$SiteID[1]
      } else if(length(out_Deep)>1){
        out_Deep[[2]]$depth <- "Deep"
        
        summary_list[[j]] <- out_Deep[[2]]
        summary_list[[j]]$SiteID <- bysite[[j]]$SiteID[1]
      }


      
    }
    
    
    hour_vmc <- do.call(rbind,bysite)
    event_summary <- do.call(rbind,summary_list)
    
#### visual checks ####
    
    site <- "AT1.1"
    mintime <- ymd("2021-02-01")#min(hour_vmc$timestamp) #ymd("2021-02-01")
    maxtime <- ymd("2021-05-01")#max(hour_vmc$timestamp) #ymd("2021-03-01")
    
    ggplot(hour_vmc[hour_vmc$SiteID==site &
                      hour_vmc$timestamp > mintime &
                      hour_vmc$timestamp < maxtime,],
           aes(x=timestamp,y=vmc_Surf,color=event_Surf,group=SiteID)) +
      geom_line(size=1) +
      theme_classic() +
      theme(panel.grid.minor.x=element_line()) +
      scale_x_datetime(date_minor_breaks="1 day")
    
    ggplot(hour_vmc[hour_vmc$SiteID==site &
                      hour_vmc$timestamp > mintime &
                      hour_vmc$timestamp < maxtime,],
           aes(x=timestamp,y=vmc_Deep,color=event_Deep,group=SiteID)) +
      geom_line(size=1) +
      theme_classic()
    
    ggplot(event_summary[event_summary$event=="prec",], 
           aes(x=start_vmc,y=rng_vmc,color=SiteID)) + 
      geom_point(alpha=0.1) + 
      geom_line(stat="smooth",method="lm") + 
      theme_classic() + 
      theme(legend.position="none")
    
    ggplot(event_summary[event_summary$event=="dem",], 
           aes(x=start_vmc,y=rng_vmc,color=SiteID)) + 
      geom_point(alpha=0.1) + 
      geom_line(stat="smooth",method="lm") + 
      theme_classic() + 
      theme(legend.position="none")
    
  # combine summary dfs to compare drainage and demand data - this only works with dem_raw=T on lines 440 and 441
  #   compare_slope_df <- event_summary[event_summary$event=="drain",
  #                                     c("starttime","endtime","rng_vmc","slope","SiteID","depth")]
  #   names(compare_slope_df) <- c("drain_starttime","drain_endtime","drain_rng","drain_slope","SiteID","depth")
  #   
  #   compare_dem <- event_summary[event_summary$event=="dem",
  #                                c("starttime","endtime","rng_vmc","slope","SiteID","depth")]
  #   names(compare_dem) <- c("dem_starttime","dem_endtime","dem_rng","dem_slope","SiteID","depth")
  #   
  #   compare_dem$drain_endtime <- compare_dem$dem_starttime - hours(1)
  #   
  #   compare_slope_df <- merge(compare_slope_df,compare_dem, all=T)
  #   
  # # plot drainage vs demand
  #   ggplot(compare_slope_df,aes(x=drain_slope*-1,y=dem_slope*-1,color=SiteID)) +
  #     geom_point(alpha=0.4) +
  #     geom_line(stat="smooth",method="lm",se=F,alpha=0.7) +
  #     theme_classic() +
  #     theme(legend.position="none") 
  #   
  #   ggplot(compare_slope_df,aes(x=drain_slope*-1,y=dem_slope*-1,color=SiteID)) +
  #     geom_point(alpha=0.2) +
  #     geom_line(stat="smooth",method="lm",se=F,alpha=0.7) +
  #     theme_classic() +
  #     theme(legend.position="none") +
  #     scale_y_continuous(trans="log",breaks=c(1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1)) +
  #     scale_x_continuous(trans="log",breaks=c(1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1))
  #   
  #   ggplot(compare_slope_df,aes(x=drain_slope*-1,y=dem_slope*-1,color=SiteID)) +
  #     geom_point(alpha=0.2) +
  #     geom_line(stat="smooth",method="lm",se=F,alpha=0.7) +
  #     geom_abline(slope=1,intercept=0) +
  #     theme_classic() +
  #     theme(legend.position="none") +
  #     scale_y_continuous(trans="log",breaks=c(1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1)) +
  #     scale_x_continuous(trans="log",breaks=c(1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1))
  # 
  # 
  #   
  #       
#### summarize data for modelling ####
    
    event_summary$date <- as_date(event_summary$starttime)
    event_summary$n_hours <- as.numeric(difftime(event_summary$endtime,event_summary$starttime,units="hours")) + 1
    event_summary$n_hours[which(event_summary$event=="dem")] <- 24
    
    names(event_summary)[which(names(event_summary)=="slope")] <- "del_slope"
    
    all_dat <- merge(event_summary, sitedata)

# save data
    write.csv(all_dat, paste0(intermediate_path,"all_events.csv"),row.names=F)

  
