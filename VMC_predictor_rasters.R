## apply vmc driver models to rasters in GMSNP
## Jordan Stark, Jan 2022

#### setup ####
  ## packages
    library(raster)
    library(rgdal)
    library(lubridate)
    library(ggplot2)
    library(rasterVis)
    library(patchwork)

  ## data import
    # model coefs
      prec_freq_coefs <- read.csv(paste0(model_out_path,"prec_freq_coefs.csv"))
      prec_amt_coefs <- read.csv(paste0(model_out_path,"prec_amt_coefs.csv"))
      dem_coefs <- read.csv(paste0(model_out_path,"dem_coefs.csv"))      
      drain_coefs <- read.csv(paste0(model_out_path,"drain_coefs.csv"))
      
    # data scaling
      prec_freq_scale <- read.csv(paste0(model_out_path,"prec_freq_scale.csv"))
      prec_amt_scale <- read.csv(paste0(model_out_path,"amt_scale.csv"))
      dem_scale <- read.csv(paste0(model_out_path,"dem_scale.csv"))
      drain_scale <- read.csv(paste0(model_out_path,"drain_scale.csv")) 
      
    # start vmc 
      prec_start <- start_vmc <- 0.2
      
    # amount of rainfall (for drainage model)
      prec_rng <- 0.08
      
    # pred_dates
      dates <- c("winter"=ymd("2020-01-01"),
                 "spring"=ymd("2020-04-01"),
                 "summer"=ymd("2020-07-01"),
                 "fall"=ymd("2020-10-01"))
      doys <- yday(dates)
      
    # rasters
      elev <- raster(paste0(gis_path,"gsmnp_ascii/elev.txt"))
      totrad <- raster(paste0(gis_path,"gsmnp_ascii/totrad.txt"))
      rad <- stack(list.files(paste0(gis_path,"gsmnp_ascii/rad/"),
                              full.names=T),quick=T)
      names(rad)[names(rad)=="rad073x2"] <- "rad074"
      
      crs(rad) <- crs(totrad) <- crs(elev) <- CRS("+proj=utm +zone=17 +datum=NAD27")
      
      rad <- rad[[which(names(rad) %in% paste0("rad",sprintf("%03d",doys)))]]
      
      slope <- terrain(elev,opt="slope",unit="degrees")
      tpi <- terrain(elev,opt="tpi")

      
      vpd_files <- list.files(paste0(gis_path,"PRISM/daily_vpd"),
                              full.names=T)
      vpd_dates <- ymd(simplify2array(strsplit(vpd_files,"vpd_"))[2,])
      
      vpd <- stack(vpd_files[which(vpd_dates %in% dates)])
      names(vpd) <- names(dates)
      
      meant_files <- list.files(paste0(gis_path,"Microclimate/MeanTs/"),
                                full.names=T)     
      meant_dates <- as_date(simplify2array(strsplit(meant_files,"y_|.tif"))[2,],format="%Y_d_%j")
      
      meant <- stack(meant_files[which(meant_dates %in% dates)])
      names(meant) <- names(dates)
      
#### function to apply models to rasters ####
      ScaleVar <- function(value, varname, scaledat){
        return((value - scaledat$means[which(scaledat$var==varname)])/scaledat$sds[which(scaledat$var==varname)])
      }
      
      
      
      MapModel <- function(coefs,scaledat,depth="Deep",trans=NA){
        coef_names <- coefs$X
        
        all_fixed_coefs <- c("(Intercept)","start_vmc","I(start_vmc^2)","depthSurf","prec_rng","prec_start","I(prec_start^2)","sindoy","cosdoy")
      
        fixed_index <- which(coef_names %in% all_fixed_coefs)
        fixed_coefs <- coef_names[fixed_index]
        
        raster_coefs <- coef_names[-fixed_index]
        raster_index <- which(coef_names %in% raster_coefs)
        
        rlist <- list(winter=list(),
                      spring=list(),
                      summer=list(),
                      fall=list())
        
        for(r in 1:length(raster_coefs)){
          varname <- raster_coefs[r]
          
          if(length(simplify2array(strsplit(varname,":"))) > 1){
            both_names <- simplify2array(strsplit(varname,":"))
            
            rasterind <- which(both_names %in% raster_coefs)
            fixedind <- ifelse(rasterind==1,2,1)
            
            if(both_names[fixedind] == "sindoy"){
              mult <- sin(doys * 0.0172)
            } else if(both_names[fixedind] == "cosdoy"){
              mult <- cos(doys * 0.0172)
            } else stop(paste0("interaction ",both_names," must be added to function"))
                        
            varname_clean <- both_names[rasterind]            
            
          } else {
            mult <- 1
            varname_clean <- varname
          }
          
          dat <- calc(ScaleVar(get(varname_clean),varname_clean,scaledat),
                      function(x) x * coefs$Estimate[which(coefs$X==varname)] * mult,
                      forceapply=T)

          if(nlayers(dat)==1){
            for(seas in 1:length(rlist)){
              rlist[[seas]][[r]] <- dat
            }
          } else if(nlayers(dat)==length(rlist)){
            for(seas in 1:length(rlist)){
              rlist[[seas]][[r]] <- dat[[seas]]
            }
          } else stop("nlayers does not equal ndates")
        }
        
        outlist <- list(winter=list(),
                        spring=list(),
                        summer=list(),
                        fall=list())
        
        for(seas in 1:length(rlist)){
         nlayers <- length(rlist[[seas]])
         
         outlist[[seas]] <- rlist[[seas]][[1]]
         
         for(r in 2:nlayers){
           outlist[[seas]] <- outlist[[seas]] + rlist[[seas]][[r]]
           
         }
         
         
         for(i in 1:length(fixed_coefs)){
           varname <- fixed_coefs[i]
           
           if(varname=="(Intercept)"){
             outlist[[seas]] <- outlist[[seas]] + coefs$Estimate[which(coefs$X==varname)]
             
           } else if(varname %in% c("start_vmc","prec_start","prec_rng")){
             outlist[[seas]] <- outlist[[seas]] + coefs$Estimate[which(coefs$X==varname)]*ScaleVar(get(varname),varname,scaledat)
             
           } else if(varname %in% c("I(start_vmc^2)","I(prec_start^2)")){
             cleaned_name <- substr(varname,3,nchar(varname)-3)
             
             outlist[[seas]] <- outlist[[seas]] + coefs$Estimate[which(coefs$X==varname)]*(ScaleVar(get(cleaned_name),cleaned_name,scaledat)^2)
             
           } else if(varname=="depthSurf"){
             if(depth=="Surf") {
               outlist[[seas]] <- outlist[[seas]] + coefs$Estimate[which(coefs$X==varname)]
             } else if(depth!="Deep") stop("Depth must be 'Surf' or 'Deep'")
             
           } else if(varname=="sindoy"){
             for(d in 1:length(doys)){
               outlist[[seas]][d] <- outlist[[seas]][d] + coefs$Estimate[which(coefs$X==varname)]*sin(doys[d] * 0.0172)
               
             }
           } else if(varname=="cosdoy"){
             for(d in 1:length(doys)){
               outlist[[seas]][d] <- outlist[[seas]][d] + coefs$Estimate[which(coefs$X==varname)]*cos(doys[d] * 0.0172)
               
             }
           } else stop(paste0("fixed var ",varname," must be added to function"))

           
           
         }

        }
        
        out <- stack(outlist)
        
        if(!is.na(trans)){
          if(trans=="exp"){
            out <- exp(out)
          } else if(trans=="sq"){
            out <- out^2
          } else if(trans=="ilogit"){
            out <-  exp(out)/(1+exp(out)) 
          } else stop("trans must be NA, 'exp', 'sq', or 'ilogit'")
          
        }
        
        return(out)
        
      }

#### calculate results ####
  dem_stack <- MapModel(dem_coefs,dem_scale,depth="Surf",trans="sq")
  drain_stack <- MapModel(drain_coefs,drain_scale,depth="Surf",trans="exp")
  prec_amt_stack <- MapModel(prec_amt_coefs,prec_amt_scale,depth="Surf",trans="exp")
  prec_freq_stack <- MapModel(prec_freq_coefs,prec_freq_scale,depth="Surf",trans="ilogit")
  
  
#### plot results ####
  pix <- 1e5 # max number of pixels to plot for each map, set to 1e7 or 1e8 for final figs or lower to test formatting
  
  PlotFunc <- function(stack,dir=-1,title=""){
    minval <- min(cellStats(stack,min))
    maxval <- max(cellStats(stack,max))
    
    plots <- list()
    for(i in 1:nlayers(stack)){
      plots[[i]] <-  gplot(stack[[i]],maxpixels=pix)+
        geom_tile(aes(fill=value)) +
        theme_void() +
        scale_fill_distiller(name="",
                             limits=c(minval,maxval),
                             palette="RdYlBu",
                             direction=dir,
                             na.value="white") +
        theme(legend.key.height=unit(0.5,"cm")) +
        labs(x="",y="",title=dates[i]) +
        coord_fixed(expand=F) 
    }
    
    print(wrap_plots(plots) + plot_layout(guides="collect") + plot_annotation(title=title))
   
    
  }

  PlotFunc(prec_freq_stack,dir=1,title="probability of rainfall leading to increase in surface vmc")
  PlotFunc(prec_amt_stack,dir=1,title="Surface vmc increase on days with precip")  
  PlotFunc(drain_stack,dir=-1,title=paste0("Surface vmc decrease following vmc increase of ",prec_rng, " due to precip"))
  PlotFunc(dem_stack,dir=-1,title="Surface vmc demand per day")
