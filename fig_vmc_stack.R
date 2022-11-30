## script to make figures from VMC stack

#### setup ####
  ## packages
  library(raster)
  library(rasterVis)
  library(ggplot2)
  library(patchwork)
  library(rgdal)

  ## logit functions
  logit <- function(x) log(x/(1-x))
  ilogit <- function(x) exp(x)/(1+exp(x))
  
  ## import data
    # moisture predictions
    medr <- raster(paste0(pred_path,"deep_vmc_summer_med.tif"))
    maxr <- raster(paste0(pred_path,"deep_vmc_summer_max.tif"))
    minr <- raster(paste0(pred_path,"deep_vmc_summer_min.tif"))
    q025r <- raster(paste0(pred_path,"deep_vmc_summer_q025.tif"))
    q975r <- raster(paste0(pred_path,"deep_vmc_summer_q975.tif"))
    sdr <- raster(paste0(pred_path,"deep_vmc_summer_sd.tif"))
    meanr <- raster(paste0(pred_path,"deep_vmc_summer_mean.tif"))
    
    rng <- maxr - minr
    rng95 <- q975r - q025r
    cv <- sdr/meanr
    
    # park boundary
    # parkbound <- readOGR(paste0(gis_path,"GRSM_data/GRSM_BOUNDARY_POLYGON/GRSM_BOUNDARY_POLYGON.shp"))
    # parkbound <- spTransform(parkbound,crs(medval))
    
    
    summary_stack <- stack(medr,maxr,minr,q025r,q975r,sdr,meanr,rng,rng95,cv)
    names(summary_stack) <- c("median","maximum","minimum","quantile_025","quantile_975","SD","mean","range","range_95","CV")

#### plot results ####
pix <- 1e8 # max number of pixels to plot for each map, set to 1e7 or 1e8 for final figs or lower to test formatting

    
scalebar_data <- data.frame(x=230000,xend=240000,y=3955000)


ne.ext <- extent(290000,314000,3949000,3964000)
scalebar_data_zoom <- data.frame(x=292000,xend=293000,y=3960500)


ne.med <- crop(medr,ne.ext)
ne.rng95 <- crop(medr,ne.ext)
ne.cv <- crop(cv,ne.ext)

PlotSingle <- function(r,title="",pal="RdYlBu",logittrans=T,dir=-1,scalebar=F,ne=F,legend=F){
  minval <- cellStats(r,min)
  maxval <- cellStats(r,max)
  
  if(ne){
    r <- crop(r,ne.ext)
  }
  
    plot <-  gplot(r,maxpixels=pix)+
      geom_tile(aes(fill=value)) +
      theme_void() +
      theme(legend.key.height=unit(0.02,"npc"),
            legend.key.width=unit(0.01,"npc"),
            text=element_text(size=14)) +
      labs(x="",y="",title=title) +
      coord_fixed(expand=F) 
    
  if(logittrans){
    plot <- plot +
      scale_fill_distiller(name="VMC (%)",
                           palette=pal,
                           direction=-dir,
                           na.value="white",
                           trans="logit",
                           breaks=c(0.05,0.25,0.45,0.65),
                           limits=c(minval,maxval))
  } else {
    plot <- plot +
      scale_fill_distiller(name="VMC (%)",
                           palette=pal,
                           direction=-dir,
                           na.value="white",
                           breaks=c(0,0.1,0.2,0.3,0.4,0.5),
                           limits=c(minval,maxval)) 
  }
    
    
  if(scalebar){
    if(!ne){
      plot <- plot +
        geom_segment(data=scalebar_data,aes(x=x,xend=xend,y=y,yend=y),
                     arrow=arrow(angle=90,ends="both",length=unit(0.02,"npc")))+
        geom_text(aes(x=mean(c(scalebar_data$x,scalebar_data$xend)),
                      y=scalebar_data$y+3000,label="10 km"),size=5) 
    } else if(ne){
      
     plot <- plot +
      geom_segment(data=scalebar_data_zoom,aes(x=x,xend=xend,y=y,yend=y),
                   arrow=arrow(angle=90,ends="both",length=unit(0.02,"npc")))+
      geom_text(aes(x=mean(c(scalebar_data_zoom$x,scalebar_data_zoom$xend)),
                    y=scalebar_data_zoom$y+1000,label="1 km"),size=5) 
    }

  }
    
  if(!legend){
    plot <- plot +
              theme(legend.position="none")
  }
  
  return(plot)
  
}

medplot <- PlotSingle(medr,title="median",scalebar=T)
medplot_ne <- PlotSingle(ne.med,ne=T,scalebar=T,legend=T)
rngplot <- PlotSingle(rng95,title="range of middle 95%",pal="PiYG",dir=1)
rngplot_ne <- PlotSingle(rng95,pal="PiYG",dir=1,ne=T,legend=T)
cvplot <- PlotSingle(cv,title="coefficient of variation",pal="BrBG",dir=-1)
cvplot_ne <- PlotSingle(cv,pal="BrBG",dir=-1,ne=T,legend=T)

full_plot <- (medplot|medplot_ne)/(rngplot|rngplot_ne)/(cvplot|cvplot_ne)


ggsave(plot=full_plot,filename=paste0(fig_path,"summer_vmc_summary_fig2.tiff"),width=7,height=7)
