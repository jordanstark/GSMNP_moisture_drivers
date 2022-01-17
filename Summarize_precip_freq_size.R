## Test pattern of precip frequency and size over time and elev

# Jan 2022
# Jordan Stark

#### setup ####
  # load packages
    library(lubridate)
    library(tidyr)
    library(lme4)
    library(MuMIn)
    library(ggplot2)

  # set paths
    dataset <- weather_path
    micro_path <- paste0(gis_path,"Microclimate/")

  # import data
    AndrewsMurphy <- read.csv(paste0(dataset,"AndrewsMurphy_weather.csv"),
                              na.strings=c("QCF","MV"))
    NOAAdata<- read.csv(paste0(dataset,"NOAA_weather.csv"),
                        na.strings=c("QCF","MV"))
    
    Elevation <- read.csv(paste0(dataset,"Weather Station Elevations.csv"))

  # clean up date formats
    NOAAdata$DATE <- ymd(NOAAdata$DATE)
    AndrewsMurphy$Date <-mdy(AndrewsMurphy$Date)
    
#### get data in the right format ####
  # make names match
    AndrewsMurphy$NAME <- "AndrewsMurphy"
    names(Elevation) <- c("NAME","Elev")
    names(AndrewsMurphy) <- c("DATE","TMAX","TMIN","PRCP","NAME")
    
#### combine data for all stations #### 
    
    NOAAdata$LATITUDE<- NULL
    NOAAdata$LONGITUDE<- NULL
    NOAAdata$ELEVATION<- NULL
    NOAAdata$STATION<- NULL
    
  #CHANGE ORDER
    AndrewsMurphy <- AndrewsMurphy[,c("NAME","DATE", "PRCP","TMAX","TMIN")] 
    
  #combine 
    all_stations <- rbind(NOAAdata,AndrewsMurphy)
    all_data <- merge(all_stations,Elevation)
    
  # summarize and prep columns for models
    prcp_thresh <- 4 # how much rain do we need to count a day of precip (cm)?

    all_data$prcp_through <- ifelse(all_data$PRCP > prcp_thresh, all_data$PRCP - prcp_thresh, 0)
    
    all_data$prec_occ <- all_data$prcp_through > 0
    all_data$sindoy <- sin(yday(all_data$DATE)*0.0172)
    all_data$cosdoy <- cos(yday(all_data$DATE)*0.0172)    

  # to look only at data after first sensor install month
    #all_data <- all_data[which(all_data$DATE > ymd("2019-10-01")),]
    
    all_data$elev.sc <- scale(all_data$Elev)

      
    
#### model precip frequency ####
    freqmod <- glmer(prec_occ ~ elev.sc*(sindoy + cosdoy) + (1|NAME), 
                     family=binomial,
                     data=all_data)
    r.squaredGLMM(freqmod)
    

    sizemod <- lmer(log(prcp_through) ~ elev.sc*(sindoy + cosdoy) + (1|NAME),
                     all_data[all_data$prec_occ,])
    qqnorm(residuals(sizemod)) 
    qqline(residuals(sizemod)) # not great
    
    r.squaredGLMM(sizemod)
    
    
#### plot model effects ####
    elevs <- c(500,875,1250,1625,2000)
    elevs_sc <- (elevs - attributes(all_data$elev.sc)$"scaled:center")/attributes(all_data$elev.sc)$"scaled:scale"
    
    
    pred.df <- data.frame(doy=rep(1:365,each=5),
                          elev.sc=rep(elevs_sc,365))  
    pred.df$sindoy <- sin(pred.df$doy * 0.0172)    
    pred.df$cosdoy <- cos(pred.df$doy * 0.0172)    

    pred.df$prob_prec <- predict(freqmod,pred.df,re.form=NA,type="response")
    pred.df$prec_cm <- exp(predict(sizemod,pred.df,re.form=NA,type="response"))
    
    ggplot(pred.df, aes(x=doy,y=prob_prec,color=elev.sc)) +
      geom_point() + 
      theme_classic()
    
    
    
    ggplot(pred.df, aes(x=doy,y=prec_cm,color=elev.sc)) +
      geom_point() + 
      theme_classic()
        