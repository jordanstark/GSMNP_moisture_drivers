## apply summer lmer model to make predictive raster stack
## Jordan Stark, Mar 2022

#### setup ####
  ## packages
    library(raster)
    library(rgdal)
    library(lubridate)
    

  ## logit functions
    logit <- function(x) log(x/(1-x))
    ilogit <- function(x) exp(x)/(1+exp(x))

  ## data import
    # model 
      load(paste0(model_out_path,"summer_drivers_lmer_2.RData"))
      coefs <- data.frame(summary(fullmod)$coefficients[,"Estimate"])
      names(coefs) <- "Estimate"
        

    # scaling
      scales <- read.csv(paste0(intermediate_path,"scaling_values_nonval_2.csv"))
      
    # dates to predict
      dates <- c(seq(ymd("2020-05-01"),ymd("2020-09-01"),by="1 day"),
                 seq(ymd("2021-05-01"),ymd("2021-07-10"),by="1 day")) # last day of microclim pred
      doys <- yday(dates)
      

    # rasters
      elev <- raster(paste0(gis_path,"gsmnp_ascii/elev.txt"))
      rad <- stack(list.files(paste0(gis_path,"gsmnp_ascii/rad/"),
                              full.names=T),quick=T)
      names(rad)[names(rad)=="rad073x2"] <- "rad074"
      
      crs(rad) <- crs(elev) <- CRS("+proj=utm +zone=17 +datum=NAD27")
      
      rad <- rad[[paste0("rad",sprintf("%03d",doys))]]
      
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

      # get yesterday's precip
      prec_dates_dm1 <- prec_dates + days(1)
      
      prec_dm1 <- stack(prec_files[which(prec_dates_dm1 %in% dates)])
      names(prec_dm1) <- names(dates)
      
      
      # transform precip rasters to crs of other rasters
      prec_tr <- projectRaster(prec,elev,method="ngb",
                               filename=paste0(intermediate_path,"prec_tmp"),
                               overwrite=T)
      prec_dm1_tr <- projectRaster(prec_dm1,elev,method="ngb",
                                   filename=paste0(intermediate_path,"precdm1_tmp"),
                                   overwrite=T)
      prec_tr <- stack(paste0(intermediate_path,"prec_tmp.gri"))
      prec_dm1_tr <- stack(paste0(intermediate_path,"precdm1_tmp.gri"))
      
      
      
      meant_files <- list.files(paste0(gis_path,"Microclimate/MeanTs/"),
                                full.names=T)     
      meant_dates <- as_date(simplify2array(strsplit(meant_files,"y_|.tif"))[2,],format="%Y_d_%j")
      
      meant <- stack(meant_files[which(meant_dates %in% dates)])
      names(meant) <- names(dates)
      
    # EVI rasters for masking
      maturity <- stack(list.files(paste0(gis_path,"Seasonality/"),
                                   pattern="Maturity_0",
                                   full.names=T))
      senescence <- stack(list.files(paste0(gis_path,"Seasonality/"),
                                     pattern="Senescence_0",
                                     full.names=T))
      
      # import QA/QC
      qa <- stack(list.files(paste0(gis_path,"Seasonality/"),
                             pattern="QA_Detailed_0",
                             full.names=T))
      qa_lookup <- read.csv(paste0(gis_path,
                                   "Seasonality/MCD12Q2-006-QA-Detailed-0-lookup.csv"))
      
      
      
      # clean stacks
      MaskQA <- function(valuename,stack){
        
        badvals <- qa_lookup$Value[which(qa_lookup[,valuename] !="Best")]
        
        qamask <- qa
        
        qamask[qamask %in% badvals ] <- -1
        
        if(!identical(extent(stack),extent(qamask))) stack <- crop(stack,qamask,)
        if(!identical(extent(stack),extent(qamask))) qamask <- crop(qamask,stack)
        
        clean <- mask(stack,qamask,maskvalue=-1)
        
        return(clean)
      }
      
      maturity <- MaskQA("Maturity",maturity)
      senescence <- MaskQA("Senescence",senescence)

      getYear <- function(names){
        substr(simplify2array(strsplit(names,"_doy"))[2,],
               start=0,stop=4)
      }
      
      makeDoy <- function(num){
        yday(as_date(num,origin="1970-01-01"))
      }
      
      
      maturity <- MaskQA("Maturity",maturity)
      senescence <- MaskQA("Senescence",senescence)
      
      values(maturity) <- makeDoy(values(maturity))
      values(senescence) <- makeDoy(values(senescence))
      
      mean_mat <- calc(maturity,mean,na.rm=T)
      mean_sen <- calc(senescence,mean,na.rm=T)
      
      mean_mat_tr <- projectRaster(mean_mat,rad[[1]])
      mean_sen_tr <- projectRaster(mean_sen,rad[[1]])
      
      rm(maturity,senescence,mean_mat,mean_sen)
      
      CreateMaskStack <- function(dates,mat,sen){
        outstack <- list()
        for(i in 1:length(dates)){
          doy <- yday(dates[i])
          
          mask_mat <- mat
          mask_mat[mask_mat>doy] <- NA
          mask_sen <- sen
          mask_sen[mask_sen<doy] <- NA
          
          outstack[[i]] <- sum(mask_mat,mask_sen)
          
        }
        outstack <- stack(outstack)
        
        return(outstack)
      }
      
      seasmask <- CreateMaskStack(dates,mean_mat_tr,mean_sen_tr)
      
      
      
#### apply scaling to rasters ####
      ScaleVar <- function(value, varname, scaledat){
        return((value - scaledat$means[which(scaledat$var==varname)])/scaledat$sds[which(scaledat$var==varname)])
      }
      
      
      elev_sc <- ScaleVar(elev,"elev",scales)
      slope_sc <- ScaleVar(slope,"slope",scales)
      tpi_sc <- ScaleVar(tpi,"tpi",scales)
      prec_sc <- ScaleVar(prec_tr,"prec",scales)
      prec_dm1_sc <- ScaleVar(prec_dm1_tr,"prec_dm1",scales)
      rad_sc <- ScaleVar(rad,"rad",scales)
      meant_sc <- ScaleVar(meant,"meant",scales)
    
        
    # save space by removing unscaled rasters
      rm(elev,slope,tpi,prec,prec_tr,prec_dm1,prec_dm1_tr,rad,meant)
      
#### apply model to rasters ####
  # values that do not change seasonally
    fixed_pred <- coefs["(Intercept)",] + 
                      coefs["elev",] * elev_sc +
                      coefs["slope",] * slope_sc +
                      coefs["tpi",] * tpi_sc 
  
 for(i in 1:nlayers(prec_sc)){

   pred <- fixed_pred + 
     coefs["prec",] * prec_sc[[i]] +
     coefs["prec_dm1",] * prec_dm1_sc[[i]] +
     coefs["meant",] * meant_sc[[i]] +
     coefs["rad",] * rad_sc[[i]] +
     coefs["prec:meant",] * prec_sc[[i]] * meant_sc[[i]] +
     coefs["prec:rad",] * prec_sc[[i]] * rad_sc[[i]] +
     coefs["meant:rad",] * meant_sc[[i]] * rad_sc[[i]] +
     coefs["prec:meant:rad",] * prec_sc[[i]] * meant_sc[[i]] * rad_sc[[i]]
   
   
   pred <- ilogit(pred)

   # mask based on leafout
   pred <- mask(pred,seasmask[[i]])
   
   writeRaster(pred,paste0(pred_path,"/deep_vmc_summer/pred_",dates[i],".tif"))
   
 }
      
    
  # final prediction 
      
  preds <- stack(list.files(paste0(pred_path,"/deep_vmc_summer/"),full.names=T))
   
  
  # summaries
  med_pred <- calc(preds,median,na.rm=T,filename=paste0(pred_path,"deep_vmc_summer_med.tif"))
  max_pred <- calc(preds,max,na.rm=T,filename=paste0(pred_path,"deep_vmc_summer_max.tif"))
  min_pred <- calc(preds,min,na.rm=T,filename=paste0(pred_path,"deep_vmc_summer_min.tif"))
  sd_pred <- calc(preds,sd,na.rm=T,filename=paste0(pred_path,"deep_vmc_summer_sd.tif"))
  mean_pred <- calc(preds,mean,na.rm=T,filename=paste0(pred_path,"deep_vmc_summer_mean.tif"))
  