setwd("C:/Users/Tom/Documents/GitRepos/TomAnderson_WIdeer/")

#pull in the raw (un-transformed) deer data
d<-readRDS(file="Results/cty.list.rds")
deer<-d$Abun
totdeer<-apply(FUN=sum,X=deer,MARGIN=2)
deeryr<-1981:2016

#extract the dvc data from file above
dvcs<-d$Crashes
dvcs<-dvcs[,7:(dim(dvcs)[2])]
totdvcs<-apply(FUN=sum,X=dvcs,MARGIN=2)
dvcyr<-1987:2016

#pull in the Fourier surrogates developed for Fig 5
abunsurr<-read.csv("Data/abunsurrsum.csv")
abunsurr<-as.matrix(abunsurr)

#pull in the Fourier surrogates developed for Fig 5
dvcsurr<-read.csv("Data/dvcsurrsum.csv")
dvcsurr<-as.matrix(dvcsurr)

#pull in raw (un-transformed) winter weather data
ww<-readRDS(file="Results/winter.clim.rds")
snow<-ww$Snwd

#pull in raw (un-transformed) winter climate index data
wc<-readRDS(file="Results/climindex.rds")
mei<-wc$WinterMEI #is the matrix of repeating values, counties X years

#filter data so all variables are the same size, due to missing values in snow data
snow_noNA<-snow[!is.na(rowMeans(snow)),] #take out values of NA where data is missing
deer_noNA<-deer[!is.na(rowMeans(snow)),] #take out counties where snow is NA
mei_noNA<-mei[!is.na(rowMeans(snow)),]

#store data in list
all.dat<-list(deer=deer_noNA,snow=snow_noNA,mei=mei_noNA)

