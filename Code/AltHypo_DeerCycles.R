#set working directory
setwd("C:/Users/Tom/Documents/GitRepos/TomAnderson_WIdeer")
#setwd("/mnt/hgfs/C/Reuman/gitrepos/TomAnderson_WIdeer")

#load wsyn
library(wsyn)

#load in raw data deer data
cty.list<-readRDS("Results/cty.list.rds")

#***DAN: Tom, pls add some lines that load in lat and lon for counties. Please use
#the same ordering as the rows of cty.list above - I need to be able to correspond
#the lats/lons to the counties themselves, which correspond to the rows of cty.list
#***TOM: see below- it reads in the coordinate and sorts them. I then add a test
#to make sure the counties are in the same order for the gps and deer files.

#read in GPS coordinates
gps <- read.csv("Data/WI_Cty_CentroidDD.csv")

#sort alphabetically by county name
gps<-gps[gps$COUNTY_NAM!="Menominee",] # drop menominee county because there's no data
levels(gps$COUNTY_NAM)[levels(gps$COUNTY_NAM)=="Saint Croix"] <- "St. Croix" #change name to match deer.clean
levels(gps$COUNTY_NAM)[levels(gps$COUNTY_NAM)=="Fond du Lac"] <- "Fond Du Lac" #change name to match deer.clean
gps<-gps[order(gps$COUNTY_NAM),]

#clean deer abundance data
deer.clean<-cleandat(cty.list$Abun,clev=5,times=1981:2016)$cdat
class(deer.clean)
dim(deer.clean) #counties by years
deer.clean
rownames(deer.clean)
colnames(deer.clean)
colnames(deer.clean)<-1981:2016

#check to make sure that row.names in in same order between deer and gps
cbind(as.character(gps$COUNTY_NAM),rownames(deer.clean))
rownames(deer.clean)==gps$COUNTY_NAM

#load in raw PDO/MEI data
climindex<-readRDS("Results/climindex.rds")
class(climindex)
length(climindex)
names(climindex)
for (counter in 1:6){print(class(climindex[[counter]]))}
for (counter in 1:6){print(dim(climindex[[counter]]))} #so these have been copied across counties already

#clean climate index data
climindex.clean<-lapply(climindex,function(x){x<-cleandat(x,times=1981:2016,clev=5)$cdat;x})
pdo.clean<-climindex.clean$WinterPDO
mei.clean<-climindex.clean$WinterMEI

#load in raw snow depth data
winter.clim<-readRDS("Results/winter.clim.rds")
class(winter.clim)
names(winter.clim)
class(winter.clim$Snwd)
dim(winter.clim$Snwd) #71 by 36
rowSums(!is.finite(winter.clim$Snwd))

#drop/fill NAs in snow data
alt<-1
if (alt==1)
{
  #drop all counties that have any NAs
  snow.tmp<-winter.clim[['Snwd']][!is.na(rowMeans(winter.clim[['Snwd']])),]
  dim(snow.tmp)
}
if (alt==2)
{
  #drop counties with many NAs, fill remaining NAs with local median value
  snow.tmp<-winter.clim[['Snwd']][rowSums(!is.finite(winter.clim[['Snwd']]))<10,]  
  dim(snow.tmp)
  for (rc in 1:dim(snow.tmp)[1])
  {
    for (cc in 1:dim(snow.tmp)[2])
    {
      if (!is.finite(snow.tmp[rc,cc]))
      {
        snow.tmp[rc,cc]<-median(snow.tmp[rc,],na.rm=TRUE)
      }
    }
  }
}
dim(snow.tmp)
sum(!is.finite(snow.tmp))

#clean snow data
snow.clean<-cleandat(snow.tmp,clev=5,times=1981:2016)$cdat

#drop rows from deer data set where snow depth had NAs so that spatial dimensions match
if (alt==1)
{
  deer.clean_snowsubset<-deer.clean[!is.na(rowMeans(winter.clim[['Snwd']])),]
}
if (alt==2)
{
  deer.clean_snowsubset<-deer.clean[rowSums(!is.finite(winter.clim[['Snwd']]))<10,]
}
dim(deer.clean_snowsubset)
sum(!is.finite(deer.clean_snowsubset))

#spatially subset PDO and MEI to run wlm models so that all the variables match in spatial 
#extent
pdo.clean<-pdo.clean[1:dim(deer.clean_snowsubset)[1],]
mei.clean<-mei.clean[1:dim(deer.clean_snowsubset)[1],]

#check dims and no NAs
dim(deer.clean_snowsubset)
sum(!is.finite(deer.clean_snowsubset))
dim(snow.clean)
sum(!is.finite(snow.clean))
dim(pdo.clean)
sum(!is.finite(pdo.clean))
dim(mei.clean)
sum(!is.finite(mei.clean))

#get coherence between deer and mei in 3-7 range for each location
cohresmei37<-NA*numeric(dim(deer.clean_snowsubset)[1])
for (counter in 1:dim(deer.clean_snowsubset)[1])
{
  h<-coh(deer.clean_snowsubset[counter,],mei.clean[counter,],times=1981:2016,norm="powall")
  cohresmei37[counter]<-mean(Mod(h$coher[h$timescales>=3 & h$timescales<=7]))
}
range(cohresmei37)
sd(cohresmei37) #they don't actualy vary that much for te alt==1 option

#do the same for 4-7 range
cohresmei47<-NA*numeric(dim(deer.clean_snowsubset)[1])
for (counter in 1:dim(deer.clean_snowsubset)[1])
{
  h<-coh(deer.clean_snowsubset[counter,],mei.clean[counter,],times=1981:2016,norm="powall")
  cohresmei47[counter]<-mean(Mod(h$coher[h$timescales>=4 & h$timescales<=7]))
}
range(cohresmei47)
sd(cohresmei47) #they still don't actualy vary that much for alt==1, but maybe a bit more

#get coherence between deer and pdo in 3-7 range for each location
cohrespdo37<-NA*numeric(dim(deer.clean_snowsubset)[1])
for (counter in 1:dim(deer.clean_snowsubset)[1])
{
  h<-coh(deer.clean_snowsubset[counter,],pdo.clean[counter,],times=1981:2016,norm="powall")
  cohrespdo37[counter]<-mean(Mod(h$coher[h$timescales>=3 & h$timescales<=7]))
}
range(cohrespdo37)
sd(cohrespdo37) #they don't actualy vary that much for te alt==1 option

#do the same for 4-7 range
cohrespdo47<-NA*numeric(dim(deer.clean_snowsubset)[1])
for (counter in 1:dim(deer.clean_snowsubset)[1])
{
  h<-coh(deer.clean_snowsubset[counter,],pdo.clean[counter,],times=1981:2016,norm="powall")
  cohrespdo47[counter]<-mean(Mod(h$coher[h$timescales>=4 & h$timescales<=7]))
}
range(cohrespdo47)
sd(cohrespdo47) #they still don't actualy vary that much for alt==1

#wavelet power for deer in 3-7 range and in 4 to 7 range
wavpow37<-NA*numeric(dim(deer.clean_snowsubset)[1])
wavpow47<-NA*numeric(dim(deer.clean_snowsubset)[1])
for (counter in 1:dim(deer.clean_snowsubset)[1])
{
  h<-wt(deer.clean_snowsubset[counter,],1981:2016)
  h2<-(colMeans(Mod(h$values)^2,na.rm=TRUE))/h$timescales
  wavpow37[counter]<-mean(h2[h$timescales>=3 & h$timescales<=7])
  wavpow47[counter]<-mean(h2[h$timescales>=4 & h$timescales<=7])  
}

#get coherence between snow depth and deer in 3-4 range
cohressnow34<-NA*numeric(dim(deer.clean_snowsubset)[1])
for (counter in 1:dim(deer.clean_snowsubset)[1])
{
  h<-coh(deer.clean_snowsubset[counter,],snow.clean[counter,],times=1981:2016,norm="powall")
  cohressnow34[counter]<-mean(Mod(h$coher[h$timescales>=3 & h$timescales<=4]))
}
range(cohressnow34)
sd(cohressnow34)

#get coherence between snow depth and deer in 3-7 range
cohressnow37<-NA*numeric(dim(deer.clean_snowsubset)[1])
for (counter in 1:dim(deer.clean_snowsubset)[1])
{
  h<-coh(deer.clean_snowsubset[counter,],snow.clean[counter,],times=1981:2016,norm="powall")
  cohressnow37[counter]<-mean(Mod(h$coher[h$timescales>=3 & h$timescales<=7]))
}
range(cohressnow34)
sd(cohressnow34)

#get wavelet power of snow depth in 3-4 range
wavpowsnow34<-NA*numeric(dim(deer.clean_snowsubset)[1])
for (counter in 1:dim(deer.clean_snowsubset)[1])
{
  h<-wt(snow.clean[counter,],1981:2016)
  h2<-(colMeans(Mod(h$values)^2,na.rm=TRUE))/h$timescales
  wavpowsnow34[counter]<-mean(h2[h$timescales>=3 & h$timescales<=4])
}

#get wavelet power of snow depth in 3-7 range
wavpowsnow37<-NA*numeric(dim(deer.clean_snowsubset)[1])
for (counter in 1:dim(deer.clean_snowsubset)[1])
{
  h<-wt(snow.clean[counter,],1981:2016)
  h2<-(colMeans(Mod(h$values)^2,na.rm=TRUE))/h$timescales
  wavpowsnow37[counter]<-mean(h2[h$timescales>=3 & h$timescales<=7])
}


#see what is related to what

#--deer power v coherence with mei
plot(cohresmei37,wavpow37,type='p')
cor.test(cohresmei37,wavpow37) #p=0.0463
plot(cohresmei47,wavpow47,type='p')
cor.test(cohresmei47,wavpow47) #p=0.075
plot(cohresmei47,wavpow37,type='p')
cor.test(cohresmei47,wavpow37) #p=0.026 

#--deer power v coherence with pdo
plot(cohrespdo37,wavpow37,type='p')
cor.test(cohrespdo37,wavpow37) #p=0.073
plot(cohrespdo47,wavpow47,type='p')
cor.test(cohrespdo47,wavpow47) #p=0.01488 
plot(cohrespdo47,wavpow37,type='p')
cor.test(cohrespdo47,wavpow37) #p=0.009806 

#--build to a linear model
cor.test(wavpow37,cohresmei47) #deer power, 3-7 with deer coherence with mei on 4-7: p=0.026
cor.test(wavpow37,cohrespdo47) #deer power, 3-7 with deer coherence with pdo on 4-7; p=0.0098
cor.test(wavpow37,wavpowsnow34) #deer power, 3-7 with snow depth power 3-4; p=0.47
cor.test(wavpow37,cohressnow34) #deer power, 3-7 with deer coherence with snow depth 3-4; p=0.89
mod1<-lm(wavpow37~cohresmei47)
mod2<-lm(wavpow37~cohrespdo47)
summary(mod1)
summary(mod2)
mod3<-lm(wavpow37~cohresmei47+cohrespdo47)
anova(mod3,mod1)
anova(mod3,mod2) #so best model so far is mod2, try to add snow depth stuff
mod4<-lm(wavpow37~cohrespdo47+wavpowsnow34+cohressnow34)
anova(mod4,mod2) #so the snow depth stuff does not help

#--build a linear model, simple version
cor.test(wavpow37,cohresmei37) #deer power, 3-7 with deer coherence with mei on 3-7: p=0.04631
cor.test(wavpow37,cohrespdo37) #deer power, 3-7 with deer coherence with pdo on 3-7; p=0.07335
cor.test(wavpow37,wavpowsnow37) #deer power, 3-7 with snow depth power 3-7; p=0.63
cor.test(wavpow37,cohressnow37) #deer power, 3-7 with deer coherence with snow depth 3-4; p=0.6687

#So one approach would be to argue the most important thing in our wavelet model
#is mei/pdo over 4-7 timescales, so test for relationship between deer wavelet power
#on 3-7 timescales with coherence of deer with mei/pdo on 4-7 timescales. Get two 
#signifcant values (see above).

#--use all 71 locs for this same thing 
#wavelet power for deer in 3-7 range and in 4 to 7 range, all 71 locs
wavpowall37<-NA*numeric(dim(deer.clean)[1])
wavpowall47<-NA*numeric(dim(deer.clean)[1])
for (counter in 1:dim(deer.clean)[1])
{
  h<-wt(deer.clean[counter,],1981:2016)
  h2<-(colMeans(Mod(h$values)^2,na.rm=TRUE))/h$timescales
  wavpowall37[counter]<-mean(h2[h$timescales>=3 & h$timescales<=7])
  wavpowall47[counter]<-mean(h2[h$timescales>=4 & h$timescales<=7])  
}

#get coherence between deer and mei in 3-7 range for all 71 location
cohresmeiall37<-NA*numeric(dim(deer.clean)[1])
cohresmeiall47<-NA*numeric(dim(deer.clean)[1])
for (counter in 1:dim(deer.clean)[1])
{
  h<-coh(deer.clean[counter,],mei.clean[1,],times=1981:2016,norm="powall")
  cohresmeiall37[counter]<-mean(Mod(h$coher[h$timescales>=3 & h$timescales<=7]))
  cohresmeiall47[counter]<-mean(Mod(h$coher[h$timescales>=4 & h$timescales<=7]))
}

#get coherence between deer and pdo in 3-7 range for each location
cohrespdoall37<-NA*numeric(dim(deer.clean)[1])
cohrespdoall47<-NA*numeric(dim(deer.clean)[1])
for (counter in 1:dim(deer.clean)[1])
{
  h<-coh(deer.clean[counter,],pdo.clean[1,],times=1981:2016,norm="powall")
  cohrespdoall37[counter]<-mean(Mod(h$coher[h$timescales>=3 & h$timescales<=7]))
  cohrespdoall47[counter]<-mean(Mod(h$coher[h$timescales>=4 & h$timescales<=7]))
}

cor.test(wavpowall37,cohresmeiall47) #p=0.037
cor.test(wavpowall37,cohrespdoall47) #p=0.00715
#pretty similar to the case with the reduced number of time series

#now try what Lawrence said
dat<-list(deer=deer.clean_snowsubset,mei=mei.clean,pdo=pdo.clean,snow=snow.clean)
wlm_abun<-wlm(dat=dat,times=1981:2016,resp=1,pred=2:4,norm="powall")
wlmvals<-wlm_abun$modval
#names(wlm_abun)
#class(wlm_abun$modval)
#dim(wlm_abun$modval)
#length(wlm_abun$timescales)
#2016-1981
#for each location "counter", get coherence between wlm_abun[counter,,] and
#deer on 3-7-year timescales
cohreslws<-NA*numeric(dim(wlmvals)[1])
for (loccounter in 1:dim(wlmvals)[1])
{
  #get deer wt for this loc
  dwto<-wt(deer.clean_snowsubset[loccounter,],1981:2016)
  dwt<-dwto$values
  
  #get model value for this loc
  mv<-wlmvals[loccounter,,]
  
  #use the coherence formula
  normdenomdwt<-sqrt(colMeans((Mod(dwt))^2,na.rm=TRUE))
  normdenommv<-sqrt(colMeans((Mod(mv))^2,na.rm=TRUE))
  for (counter in 1:dim(mv)[2])
  {
    dwt[,counter]<-dwt[,counter]/normdenomdwt[counter]
    mv[,counter]<-mv[,counter]/normdenommv[counter]
  }
  h<-colMeans(dwt*Conj(mv),na.rm=TRUE)
  
  #average appropriately and store
  cohreslws[loccounter]<-mean(Mod(h[dwto$timescales>=3 & dwto$timescales<=7]))
}

cor.test(cohreslws,wavpow37) #p=0.0036, cor=0.36986

#So use the 60 locs and take a dual approahc: argue the most important thing in our 
#wavelet model is mei/pdo over 4-7 timescales, so test for relationship between deer 
#wavelet power on 3-7 timescales with coherence of deer with mei/pdo on 4-7 timescales. 
#Get two signifcant values (see above). Then do Lawrence's approach, which is also
#significant (see above)

