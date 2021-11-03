## Processing weather data from 2011-2021
## Repeating process from Lesser & Fridley 2015 paper

# Aug 2021
# Sophie Cohen and Jordan Stark

#### setup ####
  # load packages
    library(lubridate)
    library(tidyr)

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
    
#### split data by date ####
    data_list <- split(all_data,all_data$DATE)
    
    data_list <- data_list[1:(length(data_list))]
    
#### for loop to do linear models ####
    max_models <- data.frame(date=NA, slope=NA, intercept=NA)
    min_models <- data.frame(date=NA, slope=NA, intercept=NA)
    
    
    for(i in 1:length(data_list)){
      daily_lm_max <- lm(TMAX ~ Elev, data=data_list[[i]])
      summary_lm_max <- summary(daily_lm_max)
      
      slope_max <- summary_lm_max$coefficients["Elev","Estimate"]
      intercept_max <- summary_lm_max$coefficients["(Intercept)","Estimate"]
      date_max <- as.character(data_list[[i]]$DATE[1])
      
      max_models[i,] <- c(date_max,slope_max,intercept_max)
      
      daily_lm_min <- lm(TMIN ~ Elev, data=data_list[[i]])
      summary_lm_min <- summary(daily_lm_min)
      
      slope_min <- summary_lm_min$coefficients["Elev","Estimate"]
      intercept_min <- summary_lm_min$coefficients["(Intercept)","Estimate"]
      date_min <- as.character(data_list[[i]]$DATE[1])
      
      min_models[i,] <- c(date_min,slope_min,intercept_min)
    }
    

#### save data ####
  
  write.csv(max_models,paste0(micro_path,"maxtemp_synoptic_2011_2021.csv"),
            row.names=F)
  write.csv(min_models,paste0(micro_path,"mintemp_synoptic_2011_2021.csv"),
            row.names=F)
  

  