#Antecedent Precipitation Index
#JDF 9-9-22

#API = sum (Pt * k^(-t))
  #ie, summation of previous daily rainfall, but each day's rain total is downweighted according to a constant (k) taken to that power -- seems similar to lag-1 autocorrelation

#library(hydroBE) #hydro package with getApi function

#import precipitation data, long form, daily for each site from PRISM estimates
  #modified code by JDF in PrepLongData.R file on 9-9-22

prec <- read.csv(paste0(intermediate_path,"prec_long.csv"))

#one strategy is to create an array, with third dimension the antecedent 14 days of precip

Dates = sort(unique(prec$date))
Sites = sort(unique(prec$SiteID))
Pmat = matrix(0,nrow=length(Sites),ncol=length(Dates))
rownames(Pmat) = Sites; colnames(Pmat) = Dates
for(i in 1:length(Sites)) {
  for(j in 1:length(Dates)) {
    Pmat[i,j] = prec$prec[prec$SiteID==Sites[i]&prec$date==Dates[j]]
  }
}

#reduce dates, start on 10-1-19
Pmat = Pmat[,274:944]

#create array
Pa = array(0,dim=c(dim(Pmat),14),dimnames=list(rownames(Pmat),colnames(Pmat),1:14))
Pa[,,1] = Pmat
nullcol = rep(NA,length(Sites))
Pa[,,2] = cbind(nullcol,Pa[,-671,1])
for(i in 3:14) {
  Pa[,,i] = cbind(nullcol,Pa[,-671,(i-1)])
}
  
#checking...
prec[prec$SiteID=="BX1",][c(274:330),]
Pa["BX1",1:50,]
  #seems accurate

#ignore first two weeks (which include NAs)
Pa = Pa[,-c(1:14),]

#calculate API using apply function across third dimension (and include 'current' day)
k = .9
APIa1 = apply(Pa,c(1,2),function(x){ sum (x * (k^(1:14))) } )
k = .8
APIa2 = apply(Pa,c(1,2),function(x){ sum (x * (k^(1:14))) } )
plot(APIa1,APIa2)
plot(Pa[,,1],APIa1); abline(0,1)
  
#test using getApi function 
# sit = "BX1"
# x = Pa[sit,,1]
# apix =  getApi(x,k=.9,n=14)
# plot(APIa1[sit,],apix); abline(0,1)
  #close but no cigar: the getApi function uses a slightly different method based on recursion; ignore

#melt into dataframe
library(reshape2)
apdat = as.data.frame(APIa1)
apdat = cbind(rownames(APIa1),apdat)

APIlong <- melt(apdat,id.var=1)
names(APIlong) = c("SiteID","date","API")

#estimate best fit of k using correlation with VMC
dat <- read.csv(paste0(intermediate_path,"model_data.csv"))
dat2 = merge(dat,APIlong)

plot(dat2$prec,dat2$API)
plot(dat2$API,dat2$vmc_Deep,col=as.numeric(as.factor(dat2$SiteID)))

summary(lmer(vmc_Deep~scale(API)+(scale(API)|SiteID),dat2))

sites = unique(dat2$SiteID)
for(i in 1:length(sites)) {
  site = sites[i]
  plot(dat2$API[dat2$SiteID==site],dat2$vmc_Deep[dat2$SiteID==site])
  readline()
}

#find k value producing maximum correlation of VMC (deep) and API
kvec = seq(.5,1,length=50)
Nvec = c(5,10,14) #***14 days a much higher correlation than 5 or 10 days
Dcorvec = rep(0,length(kvec)) #deep
Scorvec = rep(0,length(kvec)) #shallow
for(i in 1:length(kvec)) {
    N=14
  print(i)
  
  #calculate API values for each site
  k = kvec[i]
  APIa1 = apply(Pa,c(1,2),function(x){ sum (x * (k^(1:N))) } )
  
  #melt and merge
  apdat = as.data.frame(APIa1)
  apdat = cbind(rownames(APIa1),apdat)
  APIlong <- melt(apdat,id.var=1)
  names(APIlong) = c("SiteID","date","API")
  dat2 = merge(dat,APIlong)
  
  #spearman rank correlation of API with VMC
  Dcorvec[i] <- cor.test(dat2$vmc_Deep, dat2$API, method = 'spearman')$estimate
  Scorvec[i] <- cor.test(dat2$vmc_Surf, dat2$API, method = 'spearman')$estimate
  
}
plot(kvec,Dcorvec,ylim=c(.07,.14),pch=19)
points(kvec,Scorvec,col="blue",pch=19)

kvec[which(Dcorvec==max(Dcorvec))] #best fit k = .786 (rho = .0963)
kvec[which(Scorvec==max(Scorvec))] #best fit k = .673 (rho = .1260)

#API datasets for use in analysis
API.deep = apply(Pa,c(1,2),function(x){ sum (x * (.786^(1:14))) } )
  apdat = as.data.frame(API.deep)
  apdat = cbind(rownames(API.deep),apdat)
  APIlong.deep <- melt(apdat,id.var=1)
  names(APIlong.deep) = c("SiteID","date","APIdeep")  
  
API.surf = apply(Pa,c(1,2),function(x){ sum (x * (.673^(1:14))) } )
  apdat = as.data.frame(API.surf)
  apdat = cbind(rownames(API.surf),apdat)
  APIlong.surf <- melt(apdat,id.var=1)
  names(APIlong.surf) = c("SiteID","date","APIsurf")   

APIout = merge(APIlong.deep,APIlong.surf)  
#write.csv(APIout,file="APIdataset.csv")


# 
# 
# #API patterns: example trend of rainfall
# 
# rain = c(1,3,0,0,0,1,6,3,0,0,10,8,3,1,0,0,1,7,0,0,4,4,6,2,1,0,0,0,0)
# 
# i = 14 #look back window
# 
# library(hydroBE)
# x = getApi(rain,k=.9,n=14)
# x2= rev(rain[1:14])
# sum (x2 * (.9^(1:14)))
# 
# y = getApi(rain,k=.8,n=5)
# z = getApi(rain,k=.5,n=5)
# 
# plot(rain,x)
# par(mfrow=c(1,1))
# plot(1:length(rain),rain,type="l")
# plot(1:length(rain),x,type="l")
# lines(1:length(rain),y,type="l",col="green")
# lines(1:length(rain),z,type="l",col="red")
# 
# #plot showing effect of k on API
# kvec = seq(0,1,by=.1)
# plot(1:length(rain),getApi(rain,k=kvec[11],n=5),type="l",ylim=c(0,25))
# for(i in 2:length(kvec)) {
#   lines(1:length(rain),getApi(rain,k=kvec[i],n=5),type="l")
# }
# 
