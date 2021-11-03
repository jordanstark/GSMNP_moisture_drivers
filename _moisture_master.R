#### Master script for analysis and prediction of moisture availability in GRSM ####
## Jordan Stark
## Summer/Fall 2021

###########################################################
############## setup and directory structure ##############
###########################################################

## pathways
  sensor_script_path <- "E:/GithubRepos/Soil_temp_moisture_EMU/DataCleaning/"
    # scripts for first section on processing raw sensor code
  
  script_path <- "E:/GithubRepos/GSMNP_moisture_drivers/"
    # all other scripts
  
  sensordata_path <- "E:/Smokies_Moisture/soil_data/raw/BSensors"
    # all raw sensor data files
  
  sensormetadata_path <- "E:/Smokies_Moisture/soil_data/raw/"
    # metadata: SensorLocationHistory.csv, BSensor_BadDataList.csv, calib_coefs.csv
  
  intermediate_path <- "E:/Smokies_Moisture/soil_data/"
    # cleaned sensor data, site data, and PRISM data
  
  gis_path <- "E:/Smokies_Moisture/GIS/"
    # rasters of PRISM data, SMAP data, and Park data
    # /gsmnp_ascii/ has files from Fridley 2009
    # /Microclimate/ has synoptic temperature predictions from 2011-2021,
      # Fridley (2009) model coefficients, and predicted below-canopy temperatures
  
  weather_path <- "E:/Smokies_Moisture/weather_history/"
  
  model_out_path <- "E:/Smokies_Moisture/model_runs/"
  
  pred_path <- "E:/Smokies_Moisture/preds/"
  
  fig_path <- "E:/Smokies_Moisture/figs/"
  
  
## required packages
  library(rstudioapi) # to run scripts in background

  # all other packages are loaded in their own scripts
  
  # convenience functions
    # lubridate (date parsing)
    # stringr (text parsing)
    # tidyr (conversion between long and wide format data)
    # parallel (parallel processing)
  # GIS data processing
    # raster (raster processing)
    # rgdal (spatial formats and transformations)
    # sp (spatial point processing)
  # modelling
    # R2jags (Bayesian models interface with JAGS)
    # lme4 (linear mixed models - for checking overall patterns)
    # MuMIn (dredging for AIC and marginal/conditional r2 values for mixed models)
  # visualization
    # ggplot2
    # patchwork (multipanel figs)
    # rasterVis (spatial figs)
  
## raw inputs and data sources
  # sensor readings
    # folder of raw data in sensordata_path
    # metadata (location history) in sensormetadata_path
    # calibration parameters in sensormetadata_path
    # records of sensor failures in sensormetadata_path
  # precipitation and VPD data from PRISM
    # daily max VPD in gis_path/PRISM/MAX_VPD/[year]/
      # 2021 is split - 2021s is stable data and 2021p is provisional (more recent)
    # 4km DEM matching PRISM data in gis_path/PRISM/Elev/
  # microclimate conditions
    # weather station data following stations used in Lesser & Fridley 2015
      # in weather_path
    # microclimate model parameters from Fridley 2009
      # in gis_path/Microclimate/micro_mod
  # GIS data
    # for microclimate model:
      # TCI - ("tci.txt") in gis_path/gsmnp_ascii/
      # totrad ("totrad.txt") in gis_path/gsmnp_ascii/
      # elev ("elev.txt") in gis_path/gsmnp_ascii/
      # strdist ("logsd.txt") in gis_path/gsmnp_ascii/
      # 365 daily solar radiation rasters (folder /rad/) in gis_path/gsmnp_ascii/
    # GRSM boundary polygon in gis_path

###########################################################
################ process raw sensor data ##################
###########################################################
## clean and apply calibration
  jobRunScript(paste0(sensor_script_path,"ReadBSensors.R"),importEnv=T)
    # inputs: entire folder of raw sensor data
    #         calib.coefs.csv
    # outputs: "sensordata.csv" in intermediate_path
              # this has all calibrated sensor data
  jobRunScript(paste0(sensor_script_path,"SaveFieldData.R"),importEnv=T)
    # inputs: "sensordata.csv" from above
    #         SensorLocationHistory.csv in sensormetadata_path
    # outputs: "field_sensordata.csv" in intermediate_path
              # this has all calibrated data when sensors were deployed

## remove sensor errors and save data
  remove_list <- c("PK1","PK04b","PK04a")
  
  jobRunScript(paste0(sensor_script_path,"RemoveSensorErrors.R"),importEnv=T)
    # inputs: "field_sensordata.csv" from above
    #         "Bsensor_BadDataList.csv" in sensormetadata_path
                # note these errors were ID'd using 'CheckSensorData' script
                # which is also in sensor_script_path
    # outputs: fully cleaned, QA/QC'd sensor data in "cleaned_sensordata.csv" in intermediate_path
  
    # warnings:
      # "PK1" <- "not under canopy" 
      # "GM8" <- "very high vmcs"
      # "R4" <- "low variation in vmcs"
      # "ATC3" <- "very high vmcd"
      # "ATC4.5" <- "low vmcs"
      # "ATE02" <- "odd dips in vmc"
      # "PK04a" <- "st was 10 cm deep"
      # "PK04b" <- "st and vmcs were 10 cm deep"
      # "SD6" <- "vmcd very high, possibly need to remove starting June 2020"
      # 
###########################################################
##### prep site, climate, and sensor data for models ######
###########################################################
## prep GIS data
  jobRunScript(paste0(script_path,"Correct_TCI_raster.R"),importEnv=T)
    # inputs: tci.txt
    # outputs: tci_cor.asc
    # changes high NA values from Fridley 2009 coded as tci=25 to 95th percentile TCI in park
  
## predict relationship between above-canopy temp and elevation
  jobRunScript(paste0(script_path,"Process_weather.R"),importEnv=T)
    # Script by Sophie Cohen with minor modifications by Jordan
    # inputs: weather data ('AndrewsMurphy_weather.csv' and 'NOAA_weather.csv')
    #           from stations described in Fridley 2009 and LEsser&Fridley 2015
    #         and weather station elevations ('Weather Station Elevations.csv')
    #         in weather_path
    # outputs: 'maxtemp_synoptic_2011_2021.csv' and 'mintemp_synoptic_2011_2021.csv'
    #          these contain daily linear models for effect of elevation on temp
    #          in gis_path/Microclimate/
  
## predict relationship between VPD and elevation
  jobRunScript(paste0(script_path,"interpolate_PRISM_VPD.R"),importEnv=T)
    # takes ~ 1 hour
    # inputs: maxvpd rasters in gis_path/PRISM/Max_VPD/
    #         PRISM 4km elevation raster in gis_path/PRISM/Elev/
    #         30m DEM in gis_path/gsmnp_ascii/elev.txt
    #         GRSM boundary
    # outputs: daily rasters of maximum VPD in GRSM at 30m scale in 
                  # gis_path/PRISM/daily_vpd ; named as vpd_[date ymd].tiff

## predict below-canopy microclimate
  
  
  mainvars <- ls()
  #these are saved so that they aren't wiped out by the following scripts
  mainvars <- c(mainvars, "mainvars")
  
  source(paste0(script_path,"Micro_annual_invar.R"),local=T)
  # parallel fails with jobRunScript
  # so using source. Local=T means that the paths will copy
  # this takes ~ 15 mins to run
  # inputs: TCI ("tci_cor.asc") in gis_path/gsmnp_ascii/
  #         totrad ("totrad.txt") in gis_path/gsmnp_ascii/
  #         elev ("elev.txt") in gis_path/gsmnp_ascii/
  #         strdist ("logsd.txt") in gis_path/gsmnp_ascii/
  #         365 daily solar radiation rasters (folder /rad/) in gis_path/gsmnp_ascii/
  #         maxt model coefficients ("micro_mod/maxcoef_out.txt") in climate_path
  #         mint model coefficients ("micro_mod/mincoef_out.txt") in climate_path
  # outputs: 365 rasters each in  gis_path/Microclimate/ '/tmpMin/' '/tmpMax/'
  # files are named as 'd_001.tif' ... 'd_365.tif' 
  
  source(paste0(script_path,"Micro_alldailyrasters_present.R"),local=T)
  # this takes ~ 1 hr
  # inputs: all of the above inputs from Micro_annual_invar.R 
  #         PLUS:
  #         maxt lapse rates ("maxtemp_synoptic_1900_2011.csv") in climate_path
  #         mint lapse rates ("mintemp_synoptic_1900_2011.csv") in climate_path
  #         365 daily rasters in climate_path generate as output of Micro_annualrasters.R
  # outputs: daily rasters with temperature records
  #         climate_path "/maxTs/" has 2011-2021 max temp
  #         climate_path "/minTs/" has 2011-2021 min temp
  #         climate_path "/meanTs/" has 2011-2021 mean temp (mean of max and min)
  #         named as 'y_2011_d_001.tif' ... 'y_2021_d_365.tif'
  
## combine datasets and save for models
  jobRunScript(paste0(script_path,"PrepLongData.R"),importEnv=T)
  # inputs: cleaned sensor data ('cleaned_sensordata.csv' in intermediate_path)
  #         elevation, tci, radiation, totrad, and strdist rasters in gis_path/gsmnp_ascii/
  #         interpolated PRISM vpd data in gis_path/PRISM/daily_vpd
  #         PRISM 4km precip data in gis_path/PRISM/Precip
  #         PRISM 4km DEM in gis_path/PRISM/elev
  #         microclimate temperature data in gis/Microclimate/MeanTs
  #         EVI seasonality and peak data in gis/Seasonality/
  # outputs: in intermediate_path
  #          demand_model_data.csv has daily data when soil moisture was declining for demand model
  #          demand_fracmodel_data.csv has daily data when soil moisture declined by at least 1% of previous day moisture
  #          model_data.csv has daily data for all days for moisture availability model


  
###########################################################
######### model variation across the landscape ############
###########################################################
## model overall spatial trends in annual dataset
  jobRunScript(paste0(script_path,"Model_topotrends.R"),importEnv=T)
  # inputs: 'model_data.csv' in intermediate_path
  # outputs: 'summer_drivers_lmer.RData' in model_path with model for making figs
  #          'scaled_summer_vmc_drivers.csv' with scaled model_data from summer in intermediat_path
  
  
  
  
###########################################################
######################### figures #########################
###########################################################
## sensor deployment description
  source(paste0(script_path,"Fig_sensordeployments.R"),local=T)
  
## effects of topography on moisture
  source(paste0(script_path,"fig_topomodels.R"),local=T)
  