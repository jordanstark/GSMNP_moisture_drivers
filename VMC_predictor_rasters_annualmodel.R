## apply annual lmer model to make predictive rastesr
## Jordan Stark, Mar 2022

#### setup ####
  ## packages
    library(raster)
    library(rgdal)
    library(lubridate)
    library(ggplot2)
    library(rasterVis)
    library(patchwork)
    library(lme4)

  ## logit functions
    logit <- function(x) log(x/(1-x))
    ilogit <- function(x) exp(x)/(1+exp(x))

  ## data import
    # model 
      load(paste0(model_out_path,"annual_drivers_lmer_2.RData"))
      coefs <- data.frame(summary(fullmod2)$coefficients[,"Estimate"])
      names(coefs) <- "Estimate"
        

    # scaling
      scales <- read.csv(paste0(intermediate_path,"scaling_values_nonval.csv"))
      
    # dates to predict
      dates <- c("winter_dry"=ymd("2020-01-01"),
                 "winter_wet"=ymd("2020-01-15"),
                 "spring_wet"=ymd("2020-04-01"),
                 "spring_dry"=ymd("2020-04-15"),
                 "summer_wet"=ymd("2020-06-15"),
                 "summer_dry"=ymd("2020-07-15"),
                 "fall_dry"=ymd("2020-10-15"))
      doys <- yday(dates)
      
      sindoys <- sin(doys*0.0172)
      cosdoys <- cos(doys*0.0172)
      
    # rasters
      elev <- raster(paste0(gis_path,"gsmnp_ascii/elev.txt"))
      tci <- raster(paste0(gis_path,"gsmnp_ascii/tci.txt"))
      log_tci <- log(tci)
      rad <- stack(list.files(paste0(gis_path,"gsmnp_ascii/rad/"),
                              full.names=T),quick=T)
      names(rad)[names(rad)=="rad073x2"] <- "rad074"
      
      crs(rad) <- crs(tci) <- crs(elev) <- CRS("+proj=utm +zone=17 +datum=NAD27")
      
      rad <- rad[[which(names(rad) %in% paste0("rad",sprintf("%03d",doys)))]]
      
      slope <- terrain(elev,opt="slope",unit="degrees")
      tpi <- terrain(elev,opt="tpi")

      
      prec_files <- c(list.files(paste0(gis_path,"PRISM/Precip/2019/"),
                                 pattern="*.bil$",full.names=T),
                      list.files(paste0(gis_path,"PRISM/Precip/2020/"),
                                 pattern="*.bil$",full.names=T),
                      list.files(paste0(gis_path,"PRISM/Precip/2021s/"),
                                 pattern="*.bil$",full.names=T),
                      list.files(paste0(gis_path,"PRISM/Precip/2021p/"),
                                 pattern="*.bil$",full.names=T))
      prec_dates <- ymd(simplify2array(strsplit(prec_files,"_"))[6,])
      
      prec <- stack(prec_files[which(prec_dates %in% dates)])
      names(prec) <- names(dates)
      
      # transform precip rasters to crs of other rasters
      prec_tr <- projectRaster(prec,elev,method="ngb")
      
      meant_files <- list.files(paste0(gis_path,"Microclimate/MeanTs/"),
                                full.names=T)     
      meant_dates <- as_date(simplify2array(strsplit(meant_files,"y_|.tif"))[2,],format="%Y_d_%j")
      
      meant <- stack(meant_files[which(meant_dates %in% dates)])
      names(meant) <- names(dates)
      
#### apply scaling to rasters ####
      ScaleVar <- function(value, varname, scaledat){
        return((value - scaledat$means[which(scaledat$var==varname)])/scaledat$sds[which(scaledat$var==varname)])
      }
      
      
      elev_sc <- ScaleVar(elev,"elev",scales)
      logtci_sc <- ScaleVar(log_tci,"log_tci",scales)
      slope_sc <- ScaleVar(slope,"slope",scales)
      tpi_sc <- ScaleVar(tpi,"tpi",scales)
      prec_sc <- ScaleVar(prec_tr,"prec",scales)
      rad_sc <- ScaleVar(rad,"rad",scales)
      meant_sc <- ScaleVar(meant,"meant",scales)
      
    # save space by removing unscaled rasters
      rm(elev,log_tci,slope,tpi,prec,prec_tr,rad,meant)
      
#### apply model to rasters ####
  # values that do not change seasonally
    fixed_pred <- coefs["(Intercept)","Estimate"] + 
                      coefs["elev","Estimate"] * elev_sc +
                      coefs["log_tci","Estimate"] * logtci_sc +
                      coefs["slope","Estimate"] * slope_sc +
                      coefs["tpi","Estimate"] * tpi_sc 
  
  # final prediction including seasonal effects
    pred <- fixed_pred + 
              coefs["prec","Estimate"] * prec_sc +
              coefs["rad","Estimate"] * rad_sc +
              coefs["meant","Estimate"] * meant_sc +
              coefs["sindoy","Estimate"] * sindoys +
              coefs["cosdoy","Estimate"] * cosdoys +
              coefs["elev:sindoy","Estimate"] * sindoys * elev_sc +
              coefs["elev:cosdoy","Estimate"] * cosdoys * elev_sc +
              coefs["log_tci:sindoy","Estimate"] * sindoys * logtci_sc +
              coefs["log_tci:cosdoy","Estimate"] * cosdoys * logtci_sc +
              coefs["slope:sindoy","Estimate"] * sindoys * slope_sc +
              coefs["slope:cosdoy","Estimate"] * cosdoys * slope_sc +
              coefs["tpi:sindoy","Estimate"] * sindoys * tpi_sc +
              coefs["tpi:cosdoy","Estimate"] * cosdoys * tpi_sc +
              coefs["prec:sindoy","Estimate"] * sindoys * prec_sc +
              coefs["prec:cosdoy","Estimate"] * cosdoys * prec_sc +
              coefs["rad:sindoy","Estimate"] * sindoys * rad_sc +
              coefs["rad:cosdoy","Estimate"] * cosdoys * rad_sc +
              coefs["meant:sindoy","Estimate"] * sindoys * meant_sc +
              coefs["meant:cosdoy","Estimate"] * cosdoys * meant_sc
    
    pred <- ilogit(pred)
    
    names(pred) <- paste0(names(dates),dates)
  
#### plot results ####
  pix <- 1e7 # max number of pixels to plot for each map, set to 1e7 or 1e8 for final figs or lower to test formatting
  
  PlotFunc <- function(stack,dir=-1,title="",print=T,file=NA){
    minval <- min(cellStats(stack,min))
    maxval <- max(cellStats(stack,max))
    
    plots <- list()
    for(i in 1:nlayers(stack)){
      plots[[i]] <-  gplot(stack[[i]],maxpixels=pix)+
        geom_tile(aes(fill=value)) +
        theme_void() +
        scale_fill_distiller(name="",
                             limits=c(minval,maxval),
                             palette="RdYlGn",
                             direction=dir,
                             na.value="white") +
        theme(legend.key.height=unit(0.5,"cm")) +
        labs(x="",y="",title=names(pred)[i]) +
        coord_fixed(expand=F) 
    }
    
    out <- wrap_plots(plots) + plot_layout(guides="collect") + plot_annotation(title=title)
    
    if(!is.na(file)){
      tiff(filename=file,
           width=24,height=10,units="cm",res=400,compression="lzw")
      print(out)
      dev.off()
    }
   
    if(print) {print(out)}
    
  }

  PlotFunc(pred,dir=1,
           title="Predicted VMC",
           print=F,file=paste0(fig_path,"pred_vmc_fullmod2.tif"))

  PlotFunc(logit(pred),dir=1,
           title="logit(Predicted VMC)",
           print=F,file=paste0(fig_path,"pred_logit_vmc_fullmod2.tif"))  
  