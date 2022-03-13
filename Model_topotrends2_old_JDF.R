#### modelling broad topographic trends in sensor data ####
## edited for manuscript, Jordan Stark, Jan 2022

#### setup ####
## packages
library(lubridate)
library(lme4)
library(MuMIn)
library(lmerTest)
library(car)

## logit functions
logit <- function(x) log(x/(1-x))
ilogit <- function(x) exp(x)/(1+exp(x))

## import data
# full dataset
fulldf <- read.csv(paste0(intermediate_path,"model_data.csv"))
fulldf$SiteID <- as.factor(fulldf$SiteID)
fulldf$date <- ymd(fulldf$date)

# #### select sensor readings with matching deployment times ####
## select data from year with most sensors deployed
main_deployment <- fulldf[fulldf$date>"2020-07-10" &
                            fulldf$date<"2021-07-11",]
sensorlen <- aggregate(cbind(vmc_Surf,vmc_Deep) ~ SiteID + SensorID,
                       data=main_deployment,FUN=function(x)sum(!is.na(x)),
                       na.action=na.pass)
sensorlen$min_len <- apply(sensorlen[,c("vmc_Surf","vmc_Deep")],
                           1,FUN=function(x) min(x)/365)

good_id <- sensorlen[sensorlen$min_len>0.9,c("SensorID","SiteID")]
deep_id <- sensorlen[sensorlen$vmc_Deep/365 > 0.9,c("SensorID","SiteID") ]

SelectSensors <- function(data,IDs){
  data$fullid <- paste(data$SensorID,data$SiteID,sep="_")
  IDs$fullid <- paste(IDs$SensorID,IDs$SiteID,sep="_")
  
  data <- data[which(data$fullid %in% IDs$fullid),]
  data$fullid <- NULL
  
  return(data)
}

deep_all <- SelectSensors(main_deployment,deep_id)
deep_all[,c("vmc_Surf","delsurf","fracdem_surf")] <- list(NULL)
names(deep_all)[names(deep_all) == "vmc_Deep"] <- "vmc"

#write.csv(deep_all,paste0(intermediate_path,"deep_timeseries.csv"),row.names=F)

#### effect sizes of drivers on summer data#### JDF: remove validation sensors
scaledat <- fulldf[fulldf$doy>fulldf$mat &
                     fulldf$doy<fulldf$sen & fulldf$val==F,
                   c("vmc_Deep","SiteID","elev","slope","log_tci",
                     "tpi","prec","vpd","rad","meant","maxEVI","prec_dm1","doy")]
names(scaledat)[names(scaledat)=="vmc_Deep"] <- "vmc"

scale_values <-  data.frame(var=c("elev","slope","log_tci","tpi","prec",
                                   "prec_dm1","vpd","rad","meant","maxEVI"),
                             means=NA,
                             sds=NA)

for(i in 1:length(scale_values$var)){
  var <- scale_values$var[i]
  
  summermean <- mean(scaledat[,var],na.rm=T)
  summersd <- sd(scaledat[,var],na.rm=T)
  
  scale_values$means[i] <- summermean
  scale_values$sds[i] <- summersd
  
  scaledat[,var] <- (scaledat[,var] - summermean)/summersd
  
  
}


scaledat$logit_vmc <- logit(scaledat$vmc)

fullmod <- lmer(logit_vmc ~ elev + slope + tpi + prec_dm1 + prec*meant*rad + (1|SiteID), scaledat)

qqnorm(resid(fullmod))
qqline(resid(fullmod)) #overdispersion

summary(fullmod)
r.squaredGLMM(fullmod)
vif(fullmod)
anova(fullmod)

#BRT
library(dismo)
library(gbm)
gbm1 = gbm.step(data=na.omit(scaledat),gbm.x = c(3,4,6,7,9,10),gbm.y=15,family="gaussian",tree.complexity=2,learning.rate=.01,bag.fraction=.7)
summary(gbm1)


########## JDF: Predict validation sensors

# get mean and sd for transformations



# function to scale based on original dataset
ScaleVar <- function(value, varname, scaledat){
  return((value - scaledat$means[which(scaledat$var==varname)])/scaledat$sds[which(scaledat$var==varname)])
}

# get mean and sd for transformations

scaledatV <- fulldf[fulldf$doy>fulldf$mat &
                      fulldf$doy<fulldf$sen & 
                      fulldf$val==T,
                    c("vmc_Deep","SiteID","elev","slope","log_tci",
                      "tpi","prec","vpd","rad","meant","maxEVI","prec_dm1","doy")]
names(scaledatV)[names(scaledatV)=="vmc_Deep"] <- "vmc"

scaledatV$elev <- ScaleVar(scaledatV$elev,"elev",scale_values)
scaledatV$slope <-  ScaleVar(scaledatV$slope,"slope",scale_values)
scaledatV$log_tci <-  ScaleVar(scaledatV$log_tci,"log_tci",scale_values)
scaledatV$tpi <- ScaleVar(scaledatV$tpi,"tpi",scale_values)
scaledatV$prec <- ScaleVar(scaledatV$prec,"prec",scale_values)
scaledatV$prec_dm1 <- ScaleVar(scaledatV$prec_dm1,"prec_dm1",scale_values)
scaledatV$rad <- ScaleVar(scaledatV$rad,"rad",scale_values)
scaledatV$meant <- ScaleVar(scaledatV$meant,"meant",scale_values)

scaledatV$logit_vmc <- logit(scaledatV$vmc)

#predict with lmer model
scaledatV$logit_vmc.pred <- predict(fullmod,newdata=scaledatV,re.form=NA)

#back transform predicted values
scaledatV$vmc.pred = ilogit(scaledatV$logit_vmc.pred)

#examine predicted-observed  
plot(scaledatV$logit_vmc.pred,scaledatV$logit_vmc)
plot(scaledatV$vmc.pred,scaledatV$vmc,xlim=c(0,.4),ylim=c(0,.4))
abline(0,1)  

#model bias
accuracy = abs(scaledatV$vmc.pred-scaledatV$vmc)
bias = (scaledatV$vmc.pred-scaledatV$vmc)
plot(scaledatV$vmc,bias); abline(h=0)
median(accuracy,na.rm=T) #0.035

#validation stats with plot
plot(scaledatV$vmc,bias, 
       col=as.numeric(as.factor(scaledatV$SiteID)),
       xlim=c(0,.4),ylim=c(-.2,.2))
abline(h=0,col="gray",lty=2)
mtext("Observed VMC",side=1,line=3,cex=1.5)
mtext("Absolute Error",side=2,line=3,cex=1.5)




# save data and model for figures
save(fullmod,file=paste0(model_out_path,"summer_drivers_lmer_2.RData"))
write.csv(scale_values,paste0(intermediate_path,"scaling_values_nonval_2.csv"),row.names=F)



# save data and model for figures
save(fullmod,file=paste0(model_out_path,"summer_drivers_lmer_2.RData"))
write.csv(scaledat,paste0(intermediate_path,"scaled_summer_vmc_drivers_2.csv"))

#### summarize data annually ####
## variables of interest
topo_vars <- c("elev","log_tci","slope","tpi","maxEVI")
daily_vars <- c("prec","vpd","meant")

## annual medians
deep_ann <- aggregate(as.formula(paste0("vmc ~ SiteID +", paste(topo_vars,collapse=" + "))),
                      deep_all,median,na.rm=T,na.action=na.pass)
deep_ann$logit_vmc <- logit(deep_ann$vmc)

write.csv(deep_ann, paste0(intermediate_path,"annual_vmc_drivers.csv"),row.names=F)

## annual coefficient of variation
deep_var <- aggregate(as.formula(paste0("vmc ~ SiteID +", paste(topo_vars,collapse=" + "))),
                      deep_all,FUN=function(x) sd(x,na.rm=T)/mean(x,na.rm=T),na.action=na.pass)
names(deep_var)[names(deep_var)=="vmc"] <- "cv_vmc"

## scaled effects of drivers on moisture variability
summary(lm(cv_vmc ~ scale(elev) + scale(log_tci) + scale(slope) + scale(tpi) + scale(maxEVI),deep_var))

## effect of elevation on moisture
elev.mod <- lm(logit_vmc ~ elev, deep_ann)

# soil moisture increases by 0.00011 per masl

ilogit(predict(elev.mod,data.frame(elev=267)))
# low elevation mean vmc = 0.053

ilogit(predict(elev.mod,data.frame(elev=2025)))
# high elevation mean vmc = 0.25

0.25/0.053
# 4.7x wetter at high elevation than low elevation, compared to 1.5x more rainfall!


## full annual model
annmod <- lm(logit_vmc ~ scale(elev) + scale(log_tci) + 
               scale(slope) +
               scale(tpi) + scale(maxEVI),deep_ann,na.action=na.fail)
# slope has largest effect after elevation

slope.elev.mod <- lm(logit_vmc ~ elev+slope,deep_ann)
slope.elev.sc <- lm(logit_vmc ~ scale(elev)+scale(slope),deep_ann)

slope.elev.intxn <-  lm(logit_vmc ~ scale(elev)*scale(slope),deep_ann)

AIC(slope.elev.sc)
AIC(slope.elev.intxn)
# adding interaction makes model worse

save(slope.elev.mod,file=paste0(model_out_path,"med_slope_elev.RData"))


#### topo effect over time fig ####
moddat <- deep_all
moddat$logit_vmc <- logit(moddat$vmc)


moddat$elev <- scale(moddat$elev)
moddat$slope <- scale(moddat$slope)
moddat$log_tci <- scale(moddat$log_tci)

deepmod <- lmer(logit_vmc ~ (cos(doy*0.0172)+sin(doy*0.0172))*(elev + slope) + (1|SiteID),moddat)
r.squaredGLMM(deepmod)
summary(deepmod)
vif(deepmod)

save(deepmod,file=paste0(model_out_path,"topo_over_time.RData"))
