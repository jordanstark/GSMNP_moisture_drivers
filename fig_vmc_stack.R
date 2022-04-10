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
    
    rng <- maxr - minr
    
    # park boundary
    # parkbound <- readOGR(paste0(gis_path,"GRSM_data/GRSM_BOUNDARY_POLYGON/GRSM_BOUNDARY_POLYGON.shp"))
    # parkbound <- spTransform(parkbound,crs(medval))
    
    
    summary_stack <- stack(medr,maxr,minr)
    names(summary_stack) <- c("median","maximum","minimum")

#### plot results ####
pix <- 1e7 # max number of pixels to plot for each map, set to 1e7 or 1e8 for final figs or lower to test formatting

    
scalebar_data <- data.frame(x=230000,xend=240000,y=3955000)
    
PlotFunc <- function(rstack){
  minval <- min(cellStats(rstack,min))
  maxval <- max(cellStats(rstack,max))
  
  plots <- list()
  for(i in 1:nlayers(rstack)){
    
  
    plots[[i]] <-  gplot(rstack[[i]],maxpixels=pix)+
      geom_tile(aes(fill=value)) +

      theme_void() +
      scale_fill_distiller(name="",
                           palette="RdYlBu",
                           direction=-1,
                           na.value="white",
                           trans="logit",
                           limits=c(minval,maxval),
                           breaks=c(0.05,0.25,0.45,0.65)) +
      theme(legend.key.height=unit(0.02,"npc"),
            legend.key.width=unit(0.06,"npc"),
            legend.position="bottom",
            text=element_text(size=10)) +
      labs(x="",y="",title=names(rstack)[i]) +
      coord_fixed(expand=F) 

  }
  
  plots[[1]] <- plots[[1]] +
    geom_segment(data=scalebar_data,aes(x=x,xend=xend,y=y,yend=y),
                 arrow=arrow(angle=90,ends="both",length=unit(0.02,"npc")))+
    geom_text(aes(x=mean(c(scalebar_data$x,scalebar_data$xend)),
                  y=scalebar_data$y+3000,label="10 km"),size=3) 
  
  plots <- wrap_plots(plots,guides="collect",ncol=1) &
            theme(legend.position="bottom")
  
  return(plots)
  
}

#PlotFunc(summary_stack)


ggsave(plot=PlotFunc(summary_stack),filename=paste0(fig_path,"summer_vmc_summary.tiff"))


pix <- 1e7

PlotSingle <- function(r,title="",pal="RdYlBu",logittrans=T,dir=-1,scalebar=F){

    plot <-  gplot(r,maxpixels=pix)+
      geom_tile(aes(fill=value)) +
      theme_void() +
      theme(legend.key.height=unit(0.02,"npc"),
            legend.key.width=unit(0.01,"npc"),
            text=element_text(size=10)) +
      labs(x="",y="",title=title) +
      coord_fixed(expand=F) 
    
  if(logittrans){
    plot <- plot +
      scale_fill_distiller(name="",
                           palette=pal,
                           direction=-dir,
                           na.value="white",
                           trans="logit",
                           breaks=c(0.05,0.25,0.45,0.65)) 
  } else {
    plot <- plot +
      scale_fill_distiller(name="",
                           palette=pal,
                           direction=-dir,
                           na.value="white",
                           breaks=c(0,0.1,0.2,0.3,0.4,0.5)) 
  }
    
    
  if(scalebar){
    plot <- plot +
      geom_segment(data=scalebar_data,aes(x=x,xend=xend,y=y,yend=y),
                   arrow=arrow(angle=90,ends="both",length=unit(0.02,"npc")))+
      geom_text(aes(x=mean(c(scalebar_data$x,scalebar_data$xend)),
                    y=scalebar_data$y+3000,label="10 km"),size=3) 
  }
  
  return(plot)
  
}

medplot <- PlotSingle(medr,title="median",scalebar=T)
rngplot <- PlotSingle(rng,title="range",pal="PiYG",dir=1,logittrans=T)

#medplot/rngplot


ggsave(plot=(medplot/rngplot),filename=paste0(fig_path,"summer_vmc_summary2.tiff"))
