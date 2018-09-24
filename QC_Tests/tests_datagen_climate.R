#test for equivalence in data sets using old and new aggregation methods
rm(list=ls())
#Load new data
setwd("C:/Users/Tom/Documents/GitRepos/mbsync/AndersDeerShort")
library(wsyn)
minyear<-1981
maxyear<-2016
dens.flag<-"dnr"
source("Code/Climate_DataGeneration.R")
winter.climNew<-winter.clim
winter.clim.usdaNew<-winter.clim.usda

#Load old data
library(Reumannplatz)
library(tidyr)
library(plyr)
library(dplyr)
detr.flag<-"none"
samplesize.flag<-"n"
minyear<-1981
maxyear<-2016
dens.flag<-"dnr"
scale.flag<-"both"
source("QC_Tests/wi_climate.R")
winter.climOld<-winter.clim
winter.clim.usdaOld<-winter.clim.usda

for(i in names(winter.climNew)){
  print(all.equal(winter.climNew[[i]],winter.climOld[[i]]))
  print(identical(winter.climNew[[i]],winter.climOld[[i]]))
}
for(i in names(winter.clim.usdaNew)){
  print(all.equal(winter.clim.usdaNew[[i]],winter.clim.usdaOld[[i]]))
  print(identical(winter.clim.usdaNew[[i]],winter.clim.usdaOld[[i]]))
}
