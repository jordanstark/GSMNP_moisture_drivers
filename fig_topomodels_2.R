##### script for making figures using topographic drivers ####
## outputs from Model_topotrends
## Jordan Stark, edited for manuscript, Jan 2022

#### Setup ####
  ## packages
    library(ggplot2)
    library(patchwork)
    library(lubridate)

  ## logit functions
    logit <- function(x) log(x/(1-x))
    ilogit <- function(x) exp(x)/(1+exp(x))
    
  ## import data
    # all data
    fulldf <- read.csv(paste0(intermediate_path,"model_data.csv"))
    fulldf$SiteID <- as.factor(fulldf$SiteID)
    fulldf$date <- ymd(fulldf$date)
    
    # scaled summer data
    scaledat <- read.csv(paste0(intermediate_path,"scaled_summer_vmc_drivers_mod2.csv"))
    scaledat$SiteID <- as.factor(scaledat$SiteID)
    
    # scaling values
    scalevals <- read.csv(paste0(intermediate_path,"scaling_values_summer_mod2.csv"))
    
    # data from deep sensors that worked summer 2020-summer 2021
    deep_all <- read.csv(paste0(intermediate_path,"deep_timeseries.csv"))
    deep_all$date <- ymd(deep_all$date)
    deep_all$SiteID <- factor(deep_all$SiteID)
    
    # scaled model of topo, met and veg drivers on summer moisture
    load(paste0(model_out_path,"summer_drivers_lmer_mod2.RData"))
    
    # model of drivers of annual medians
    load(paste0(model_out_path,"med_slope_elev.RData"))
    
    # model of seasonal changes in topo effects
    load(paste0(model_out_path,"topo_over_time.RData"))
    
#### prep data
    # calculate annual medians
    topo_vars <- c("elev","log_tci","slope","tpi")
    daily_vars <- c("rad","prec","meant")
    
    ## annual medians
    deep_ann <- aggregate(as.formula(paste0("vmc ~ SiteID +", paste(topo_vars,collapse=" + "))),
                          deep_all,median,na.rm=T,na.action=na.pass)
    deep_ann$logit_vmc <- logit(deep_ann$vmc)
    

    
    
#### figure of effect of drivers on summer moisture ####
    vars <- c("elev","log_tci","slope","tpi","prec","rad","meant")
    labels <- c("elevation (masl)","log(TCI)","slope (degrees)","TPI","precipitation (mm)","radiation (W/m2)","below-canopy temperature (degrees)")
    
    plots <- list()
    
    for(i in 1:length(vars)){
      var <- vars[i]
      
      plotdat <- scaledat
      plotdat$xval <- (plotdat[,var] * scalevals$sds[scalevals$var==var]) + scalevals$means[scalevals$var==var]
      plotdat$elev <- (plotdat[,"elev"] * scalevals$sds[scalevals$var=="elev"]) + scalevals$means[scalevals$var=="elev"]
      
      fig.df <- data.frame(elev=rep(0,41))
      fig.df$slope <- 0
      fig.df$log_tci <- 0    
      fig.df$tpi <- 0
      fig.df$prec <- 0
      fig.df$rad <- 0
      fig.df$meant <- 0

      fig.df[,var] <- seq(min(scaledat[,var],na.rm=T),max(scaledat[,var],na.rm=T),length.out=41)
      
      fig.df$vmc <- ilogit(predict(fullmod,fig.df,re.form=NA))
      fig.df$xval <- (fig.df[,var] * scalevals$sds[scalevals$var==var]) + scalevals$means[scalevals$var==var]
      
      
      
      plots[[i]] <- ggplot(fig.df, aes(x=xval,y=vmc)) +
        geom_point(data=plotdat,aes(color=elev),alpha=0.1) +
        geom_line() +
        theme_bw() + 
        theme(legend.position="bottom",
              text=element_text(size=10)) +
        labs(x=labels[i],color="elevation (masl)") +
        guides(color=guide_colorbar(title.position="top"))
      
      if(!i %in% c(1,3,5,7)){
        plots[[i]] <- plots[[i]] +
          theme(axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank())
      }
    }  
    
    wrap_plots(plots,guides="collect") + guide_area() + plot_layout(ncol=2)
    
    ggsave(paste0(fig_path,"fig6_revised_elev.png"),width=4,height=6)
    

#### comparison of deep and surface moisture across sites ####
    ggplot(fulldf, aes(x=logit(vmc_Deep),y=logit(vmc_Surf),color=SiteID)) +
      geom_point(alpha=0.05) +
      geom_smooth(method="lm",se=F,alpha=0.3,size=0.4) +
      geom_abline(slope=1,intercept=0) +
      theme_classic() + 
      coord_fixed() +
      geom_density(inherit.aes=F,aes(y=logit(vmc_Surf)),position=position_nudge(x=-0.38),bw=0.05) +
      geom_density(inherit.aes=F,aes(x=logit(vmc_Deep)),position=position_nudge(y=-0.38),bw=0.05) +
      geom_hline(yintercept=-0.38) +
      geom_vline(xintercept=-0.38) +
      theme(legend.position="none") +
      labs(x="logit(deep soil moisture) (v/v)",y="logit(surface soil moisture) (v/v)")
    
#### effect of slope and elevation on median annual moisture ####
    test.df <- data.frame(elev=rep(seq(267,2025,length.out=50),5),
                          slope=rep(seq(0,40,length.out=5),each=50))
    
    test.df$pred <- ilogit(predict(slope.elev.mod,test.df))
    
    # check elevation effect
    test.df$pred[which(test.df$elev==2025 & test.df$slope==0)]/
      test.df$pred[which(test.df$elev==267 & test.df$slope==0)]
    
    # check slope effect
    test.df$pred[which(test.df$elev==267 & test.df$slope==0)]/
      test.df$pred[which(test.df$elev==267 & test.df$slope==40)]
    
    # high elevation site with 40 deg slope is equivalent to flat site at
    predict(slope.elev.mod,data.frame(elev=2025,slope=40))
    #-1.422 = coef(slope.elev.mod)["elev"] * x + coef(slope.elev.mod)["slope"]*0
    1.422/coef(slope.elev.mod)["elev"]
    predict(slope.elev.mod,data.frame(elev=1325,slope=0))
    2025-1325
    
    
    plot_summarydat <- deep_ann
    plot_summarydat$slope.cat <- round(plot_summarydat$slope,digits=-1)
    
    ggplot(test.df, aes(x=elev,y=pred,color=factor(slope))) +
      geom_line() +
      geom_point(data=plot_summarydat,aes(x=elev,y=vmc,color=factor(slope.cat)),size=2) +
      scale_color_viridis_d("slope",option="plasma") +
      theme_classic() +
      labs(x="elevation (m a.s.l.)",y="median vmc (v/v)")
    
    
    ggplot(test.df, aes(x=elev,y=pred,color=factor(slope))) +
      geom_line(size=2) +
      geom_point(data=plot_summarydat,
                 aes(x=elev,y=vmc,color=factor(slope.cat)),
                 size=3,alpha=0.8) +
      scale_color_viridis_d("slope",option="plasma") +
      theme_classic() +
      labs(x="elevation (m a.s.l.)",y="median vmc (v/v)") +
      theme(text=element_text(size=20))
    
#### seasonal effects of slope and elevation ####
    moddat <- deep_all
    moddat$logit_vmc <- logit(moddat$vmc)
    
    slope.m <- mean(moddat$slope)
    slope.sd <- sd(moddat$slope)
    
    moddat$elev <- scale(moddat$elev)
    moddat$slope <- scale(moddat$slope)
    moddat$log_tci <- scale(moddat$log_tci)
    
    
    test.df3 <- data.frame(doy=rep(1:365,each=15))
    test.df3$elev <- rep(c(-1,0,1),5)
    test.df3$slope <- rep(c(-2,-1,0,1,2),each=3)
    
    test.df3$vmc <- ilogit(predict(deepmod,test.df3,re.form=NA))
    
    test.df3$slope <- factor(test.df3$slope,levels=seq(-2,2,by=1),
                             labels=round(seq(slope.m-(2*slope.sd),slope.m+(2*slope.sd),slope.sd),0),
                             ordered=T)
    test.df3$elev <- factor(test.df3$elev,levels=c(-1,0,1),
                            labels=c("low elevation","mid elevation","high elevation"),ordered=T)
    
    
    plotdat <- moddat
    plotdat$slope <- factor(round(plotdat$slope),levels=seq(-2,2,by=1),
                            labels=round(seq(slope.m-(2*slope.sd),slope.m+(2*slope.sd),slope.sd),0),
                            ordered=T)
    plotdat$elev <- round(plotdat$elev)
    plotdat$elev <- factor(plotdat$elev,levels=c(-1,0,1),labels=c("low elevation","mid elevation","high elevation"),ordered=T)
    
    ggplot(test.df3, aes(x=doy,y=vmc,color=slope)) +
      geom_line(data=plotdat[which(!is.na(plotdat$elev)),],
                aes(group=interaction(SiteID,year(date))),alpha=0.5) +
      geom_line(size=2) +
      scale_color_viridis_d("slope",option="plasma") +
      theme_classic() +
      facet_wrap(~elev,nrow=1) +
      labs(x="day of year",y="vmc")
    
    ggplot(test.df3, aes(x=doy,y=vmc,color=slope)) +
      geom_line(data=plotdat[which(!is.na(plotdat$elev)),],
                aes(group=interaction(SiteID,year(date))),alpha=0.5) +
      geom_line(size=2) +
      scale_color_viridis_d("slope",option="plasma") +
      theme_classic() +
      facet_wrap(~elev,nrow=1) +
      labs(x="day of year",y="vmc") +
      theme(text=element_text(size=24))
    
    