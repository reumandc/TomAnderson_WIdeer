#This script is to do some lagged regressions to help with making sure the ARX models we use are reasonable.

#***
#Do some lagged regressions of deer against lagged deer, snow and mei
#***

#***get the data

#pull in the raw (un-transformed) deer data
d<-readRDS(file="Results/cty.list.rds")
deer<-d$Abun
totdeer<-apply(FUN=sum,X=deer,MARGIN=2)
deeryr<-1981:2016

#pull in raw (un-transformed) winter weather data, same dimensions as deer, though with more NAs (see below)
ww<-readRDS(file="Results/winter.clim.rds")
snow<-ww$Snwd
snow<-snow[c(1:54,56:59,55,60:71),]

#pull in raw (un-transformed) winter climate index data
wc<-readRDS(file="Results/climindex.rds")
mei<-wc$WinterMEI #is the matrix of repeating values, counties X years, same dimensions as deer
mei<-mei[c(1:54,56:59,55,60:71),] #except for getting the rownames right, this makes no difference, because 
#data are the same in all locations 

#***clean the data

#Trim all three datasets to the same locations where all data are available
goodinds<-unname(which(!is.na(rowMeans(snow))))
snow<-snow[goodinds,]
deer<-deer[goodinds,]
mei<-mei[goodinds,]

#clean and transform
deer<-wsyn::cleandat(deer,clev=5,times=deeryr)$cdat
snow<-wsyn::cleandat(snow,clev=5,times=deeryr)$cdat
mei<-wsyn::cleandat(mei,clev=5,times=deeryr)$cdat

#***set up the regression formula and the right data frame containing all the transitions to use

#lags to use
ARo_deer<-3 #0 for no deer effect, larger values for lagged-year deer effects
ARo_snow<-3 #-1 for no snow effect, 0 for a current-year effect, larger values for additional lagged effects
ARo_mei<-3 #-1 for no mei effect, 0 for a current-year effect, larger values for additional lagged effects
maxlag<-max(ARo_deer,ARo_snow,ARo_mei)

regdat<-data.frame(deer_resp=as.vector(deer[,(maxlag+1):(ncol(deer))]))
regform<-"deer_resp~"
d2rd<-1
if (ARo_deer>0)
{
  #in this case we need to add deer predictors
  for (counter in 1:ARo_deer)
  {
    regdat<-cbind(regdat,as.vector(deer[,(maxlag+1-counter):(ncol(deer)-counter)]))
    d2rd<-d2rd+1
    colnames(regdat)[d2rd]<-paste0("deer_pred_",counter)
    regform<-paste0(regform,"deer_pred_",counter,"+")
  }
}  
if (ARo_snow>-1)
{
  #in this case we need to add snow predictors
  for (counter in 0:ARo_snow)
  {
    regdat<-cbind(regdat,as.vector(snow[,(maxlag+1-counter):(ncol(deer)-counter)]))
    d2rd<-d2rd+1
    colnames(regdat)[d2rd]<-paste0("snow_pred_",counter)
    regform<-paste0(regform,"snow_pred_",counter,"+")
  }
}
if (ARo_mei>-1)
{
  #in this case we need to add mei predictors
  for (counter in 0:ARo_mei)
  {
    regdat<-cbind(regdat,as.vector(mei[,(maxlag+1-counter):(ncol(deer)-counter)]))
    d2rd<-d2rd+1
    colnames(regdat)[d2rd]<-paste0("mei_pred_",counter)
    regform<-paste0(regform,"mei_pred_",counter,"+")
  }
}
regform<-paste0(regform,"(-1)") #no intercept for this regression
lmres_deer_vs_lags<-lm(formula=as.formula(regform),data=regdat)
#summary(lmres_deer_vs_lags)
coefs<-coef(lmres_deer_vs_lags)

if (any(Mod(polyroot(c(1,-coefs[1:ARo_deer])))<=1)) {stop("Error in Review2_LaggedRegressions.R: the filter associated with deer is not invertible, as assumed")}

#***
#Do some lagged regressions of DVCs against lagged DVCs and deer
#***

#I was going to do this, but then I realized we don't really have any reason to think DVCs depend on lagged deer - 
#they probably only depend on current deer, and we already did a cross correlation analysis in a separate file.
