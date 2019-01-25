#set working directory
setwd("C:/Users/Tom/Documents/GitRepos/TomAnderson_WIdeer")

#load in raw data deer data
cty.list<-read.rds("Results/wlm_abun.rds")

#clean deer abundance data
deer.clean<-cleandat(cty.list$Abun,clev=5,times=1981:2016)

#load in raw PDO/MEI data
climindex<-read.rds("Results/climindex.rds")

#clean climate index data
climindex.clean<-lapply(climindex,function(x){x<-cleandat(x,times=1981:2016,clev=5)$cdat;x})
pdo.clean<-climindex.clean$WinterPDO
mei.clean<-climindex.clean$WinterMEI

#load in raw snow depth data
winter.clim<-read.rds("Results/winter.clim.rds")

#drop rows with NAs in snow data
snow.tmp<-winter.clim[['Snwd']][!is.na(rowMeans(winter.clim[['Snwd']])),]

#clean snow data
snow.clean<-cleandat(snow.tmp,clev=5,times=minyear:maxyear)$cdat

#drop rows from deer data set where snow depth had NAs so that spatial dimensions match
#note: would need to repeat this process on PDO and MEI to run wlm models so that they all match in spatial extent
deer.clean_snowsubset<-deer.clean[!is.na(rowMeans(winter.clim[['Snwd']])),]

#load in wavelet linear model of deer abundance (predictors snow depth, PDO and MEI)
wlm_abun<-read.rds("Results/wlm_abun.rds")
