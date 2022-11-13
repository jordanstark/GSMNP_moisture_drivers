#### script to calculate API across stack of precip data 
## based on API_tests.R
## packages
library(raster)
library(rgdal)
library(lubridate)

# dates to predict -- starting 15 days prior to desired dates for vmc stack to allow API calc through entire stack
dates_2020 <- c(seq(ymd("2020-04-15"),ymd("2020-09-01"),by="1 day"))
dates_2021 <- c(seq(ymd("2021-04-15"),ymd("2021-07-10"),by="1 day")) # last day of microclim pred

prec_files <- c(list.files(paste0(gis_path,"PRISM/Precip/2019/"),
                           pattern="*.bil$",full.names=T),
                list.files(paste0(gis_path,"PRISM/Precip/2020/"),
                           pattern="*.bil$",full.names=T),
                list.files(paste0(gis_path,"PRISM/Precip/2021s/"),
                           pattern="*.bil$",full.names=T),
                list.files(paste0(gis_path,"PRISM/Precip/2021p/"),
                           pattern="*.bil$",full.names=T))
prec_dates <- ymd(simplify2array(strsplit(prec_files,"_"))[6,])

prec_2020 <- stack(prec_files[which(prec_dates %in% dates_2020)])
names(prec_2020) <- dates_2020

prec_2021 <- stack(prec_files[which(prec_dates %in% dates_2021)])
names(prec_2021) <- dates_2021


# calculate API
k <- 0.786

for(i in 14:nlayers(prec_2020)){
  
 out <- prec_2020[[i]] * k^(-1) +
    prec_2020[[(i-1)]] * k^(-2) +
    prec_2020[[i-2]] * k^(-3) +
    prec_2020[[i-3]] * k^(-4) +
    prec_2020[[i-4]] * k^(-5) +
    prec_2020[[i-5]] * k^(-6) +
    prec_2020[[i-6]] * k^(-7) +
    prec_2020[[i-7]] * k^(-8) +
    prec_2020[[i-8]] * k^(-9) +
    prec_2020[[i-9]] * k^(-10) +
    prec_2020[[i-10]] * k^(-11) +
    prec_2020[[i-11]] * k^(-12) +
    prec_2020[[i-12]] * k^(-13) +
    prec_2020[[i-13]] * k^(-14)
 
 writeRaster(out,paste0(gis_path,"/API/API_",dates_2020[i],".tiff"))
    
}


for(i in 14:nlayers(prec_2021)){
  
  out <- prec_2021[[i]] * k^(-1) +
    prec_2021[[(i-1)]] * k^(-2) +
    prec_2021[[i-2]] * k^(-3) +
    prec_2021[[i-3]] * k^(-4) +
    prec_2021[[i-4]] * k^(-5) +
    prec_2021[[i-5]] * k^(-6) +
    prec_2021[[i-6]] * k^(-7) +
    prec_2021[[i-7]] * k^(-8) +
    prec_2021[[i-8]] * k^(-9) +
    prec_2021[[i-9]] * k^(-10) +
    prec_2021[[i-10]] * k^(-11) +
    prec_2021[[i-11]] * k^(-12) +
    prec_2021[[i-12]] * k^(-13) +
    prec_2021[[i-13]] * k^(-14)
  
  writeRaster(out,paste0(gis_path,"API/API_",dates_2021[i],".tiff"))
  
}
