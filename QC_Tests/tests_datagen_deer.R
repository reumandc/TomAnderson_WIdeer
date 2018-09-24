#test for equivalence in data sets using old and new aggregation methods
rm(list=ls())
#Load new data
setwd("C:/Users/Tom/Documents/GitRepos/mbsync/AndersDeerShort")
library(wsyn)
library(tools)
library(reshape2)
minyear<-1981
maxyear<-2016
dens.flag<-"dnr"
source("Code/Deer_DataGeneration.R")
cty.listNew<-cty.list
usda.listNew<-usda.list

#Load old data
library(Reumannplatz)
library(tidyr)
library(plyr)
library(dplyr)
detr.flag<-"none"
minyear<-1981
maxyear<-2016
dens.flag<-"dnr"
source("Tests/wi.deer.cleandata.R")
cty.listOld<-cty.list
usda.listOld<-usda.list

for(i in names(cty.listNew)){
  print(all.equal(cty.listNew[[i]],cty.listOld[[i]]))
  print(identical(cty.listNew[[i]],cty.listOld[[i]]))
}

for(i in names(usda.listNew)){
  print(all.equal(usda.listNew[[i]],usda.listOld[[i]]))
  print(identical(usda.listNew[[i]],usda.listOld[[i]]))
}
