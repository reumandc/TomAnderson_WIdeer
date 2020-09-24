#This script is to do some ML fitting of ARMA-type models to deer data, along with some prelims. This is part of a much
#larger effort to respond to one of the Ecol Lett reviews, second round of reviews. But this script is just focussed
#on checking that deer can reasonably be modelled as an ARX model driven by MEI and snow depth. Reuman.

#***
#***Load the data
#***

#pull in the raw (un-transformed) deer data
d<-readRDS(file="Results/cty.list.rds")
deer<-d$Abun
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

#***
#Clean data
#***

#Trim all three datasets to the same locations where all data are available
goodinds<-unname(which(!is.na(rowMeans(snow))))
snow<-snow[goodinds,]
deer<-deer[goodinds,]
mei<-mei[goodinds,]

#clean and transform
deer<-wsyn::cleandat(deer,clev=5,times=deeryr)$cdat
snow<-wsyn::cleandat(snow,clev=5,times=deeryr)$cdat
mei<-wsyn::cleandat(mei,clev=5,times=deeryr)$cdat

# #***
# #Make some plots of the cleaned, transformed data - for purposes of orienting myself
# #***
# 
# plot(deeryr,deer[1,],type="l",ylim=range(deer),main="Deer, cleaned and transformed",
#      xlab="Year",ylab="Transformed deer")
# for (counter in 2:(dim(deer)[1]))
# {
#   lines(deeryr,deer[counter,],type="l")
# }
# 
# plot(deeryr,snow[1,],type="l",ylim=range(snow),main="Snow depth, cleaned and transformed",
#      xlab="Year",ylab="Transformed snow depth")
# for (counter in 2:(dim(snow)[1]))
# {
#   lines(deeryr,snow[counter,],type="l")
# }
# 
# plot(deeryr,mei[1,],type="l",ylim=range(mei),main="MEI, cleaned and transformed",
#      xlab="Year",ylab="Transformed MEI")

#***
#make some partial autocorrelation diagrams for deer, to get a reasonable bound for AR lags to consider later
#***

#For making a partial autocorrelation function for this case, which just averages the 
#partial autocorrelation functions of the locations.
#
#d        A matrix with time series in the rows
#maxlag   The maximum lag to use
#
#Output - a vector of length maxlag
#
pautocorrfun<-function(d,maxlag)
{
  dd1<-dim(d)[1]

  res<-0*numeric(maxlag)
  for (counter in 1:dd1)
  {
    h<-stats::pacf(d[counter,],lag.max=maxlag,plot=FALSE)
    res<-res+h$acf
  }
  res<-res/dd1
  
  return(res)
}

#get a partial autocorrelation function for deer
maxlag<-15
deer_pauto<-pautocorrfun(deer,maxlag=15)

#test against reshuffled data
nsurr<-1000
deer_pauto_s<-matrix(NA,nsurr,length(deer_pauto))
set.seed(101)
for (counter in 1:nsurr)
{
  deer_s<-deer[,sample(dim(deer)[2],dim(deer)[2],replace=TRUE)]
  deer_pauto_s[counter,]<-pautocorrfun(deer_s,maxlag)
}
qs<-apply(FUN=quantile,X=deer_pauto_s,MARGIN=2,probs=c(.025,.975))

tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier_Step2_part2_pautocorr_deer.png",res=600,units="in",width = tot.wd,height = tot.ht)
plot(1:maxlag,deer_pauto,type="b",ylim=range(deer_pauto,qs))
lines(c(1,maxlag),c(0,0),type="l",lty="solid")
lines(1:maxlag,qs[1,],lty="dotted")
lines(1:maxlag,qs[2,],lty="dotted")
dev.off()
#suggests lags of up to 2 are all that are needed, maybe use up to 3 for thoroughness

#***
#make some cross correlation diagrams, to get reasonable bounds on lags to consider later
#***

#For making cross correlations functions for this case by averaging across locations
#
#d1         A matrix with time series in the rows
#d2         Another such, considered the "driver" so the only lags considered are ones where this
#             variable "goes first"
#maxlag     The maximum lag to use
#
#Output - a data frame with these entries
#lag        0 up to maxlag
#cor        The lagged correlation
#
crosscorrfun<-function(d1,d2,maxlag)
{
  if (any(dim(d1)!=dim(d2)))
  {
    stop("Error in crosscorrfun: nonconformable arguments d1 and d2")
  }
  
  dd2<-dim(d1)[2]
  
  res<-data.frame(lag=0:maxlag,cor=NA)
  for (lag in 0:maxlag)
  {
    res[res$lag==lag,2]<-mean(diag(cor(t(d1[,(lag+1):dd2]),t(d2[,1:(dd2-lag)]))))
  }
  
  return(res)
}

#get the cross correlation function for deer and snow
deer_snow_crosscor<-crosscorrfun(deer,snow,maxlag)

#get the cross correlation function for deer and mei
deer_mei_crosscor<-crosscorrfun(deer,mei,maxlag)

#test against surrogate data
deer_s<-wsyn::surrog(dat=deer,nsurrogs=nsurr,surrtype="fft",syncpres=TRUE)
mei_s<-wsyn::surrog(dat=mei,nsurrogs=nsurr,surrtype="fft",syncpres=TRUE)
snow_s<-wsyn::surrog(dat=snow,nsurrogs=nsurr,surrtype="fft",syncpres=TRUE)
deer_mei_crosscor_s<-matrix(NA,nsurr,dim(deer_mei_crosscor)[1])
deer_snow_crosscor_s<-matrix(NA,nsurr,dim(deer_snow_crosscor)[1])
for (counter in 1:nsurr)
{
  h<-crosscorrfun(deer_s[[counter]],snow_s[[counter]],maxlag)
  deer_snow_crosscor_s[counter,]<-h$cor
  h<-crosscorrfun(deer_s[[counter]],mei_s[[counter]],maxlag)
  deer_mei_crosscor_s[counter,]<-h$cor
}
qs_snow<-apply(FUN=quantile,MARGIN=2,X=deer_snow_crosscor_s,prob=c(.025,.975))
qs_mei<-apply(FUN=quantile,MARGIN=2,X=deer_mei_crosscor_s,prob=c(.025,.975))

#make plots
png("Results/Fourier_Step2_part2_crosscorr_deer_snow.png",res=600,units="in",width = tot.wd,height = tot.ht)
plot(deer_snow_crosscor$lag,deer_snow_crosscor$cor,type="b",ylim=range(deer_snow_crosscor$cor,qs_snow))
lines(range(deer_snow_crosscor$lag),c(0,0),type="l",lty="dashed")
lines(deer_snow_crosscor$lag,qs_snow[1,],lty="dotted")
lines(deer_snow_crosscor$lag,qs_snow[2,],lty="dotted")
dev.off()
#suggests lags up to 2 are enough, maybe use 3 for thoroughness

png("Results/Fourier_Step2_part2_crosscorr_deer_mei.png",res=600,units="in",width = tot.wd,height = tot.ht)
plot(deer_mei_crosscor$lag,deer_mei_crosscor$cor,type="b",ylim=range(deer_mei_crosscor$cor,qs_mei))
lines(range(deer_mei_crosscor$lag),c(0,0),type="l",lty="dashed")
lines(deer_mei_crosscor$lag,qs_mei[1,],lty="dotted")
lines(deer_mei_crosscor$lag,qs_mei[2,],lty="dotted")
dev.off()
#suggests lags up to 2 are enough, maybe use 3 for thoroughness

#***
#Now do some fitting
#***

#***A likelihood function and testing

#A likelihood function for one class of models, for which deer can be autoregressive, and past snow and 
#MEI can influence deer, and innovations are incorporated into the model without any moving averaging, 
#and have the same standard deviation everywhere and the same correlation between all pairs of locations. 
#
#Args
#parms        All the parameters stacked into a vector, in the order indicated by the below arguments
#ARo_deer     AR order for deer. Zero means past deer don't influence current deer.
#ARo_snow     AR order for snow - actually just the number of lags of snow that are allowed to influence
#               deer in the model. Zero means only the present year of snow influences deer. -1 means snow
#               does not influence deer at all.
#ARo_mei      AR order for MEI - similar to ARo_snow, but for MEI.
#sd_innov     Standard deviation parameter for the innovations. NOT AN ARGUMENT BUT CONTAINED IN parms.
#cor_innov    Correlation parameters for the innovations. NOT AN ARGUMENT BUT CONTAINED IN parms.
#fixedminyr   Used for making sure the same amount of data are used for different lags, for comparability of AICs
#deer         Deer data, pre-treated as above
#snow         Snow data, pre-treated as above
#mei          MEI data, pre-treated as above
#deeryr       Years for the data, as above
#innovs       TRUE to also return the innovations. Default FALSE.
#
#Output - the ln likelihood of the model
#
logLik_class1<-function(parms,ARo_deer,ARo_snow,ARo_mei,fixedminyr,deer,snow,mei,deeryr,innovs=FALSE)
{
  #***no systematic error checking is done, just a few things
  numparms<-ARo_deer+2
  if (ARo_snow>=0)
  {
    numparms<-numparms+ARo_snow+1
  }
  if (ARo_mei>=0)
  {
    numparms<-numparms+ARo_mei+1
  }
  if (length(parms)!=numparms)
  {
    stop("Error in logLik_class1: parms has wrong length")
  }
  
  #***it is convenient to put recent years on the left in matrices
  deeryr<-rev(deeryr)
  deer<-unname(deer[,(ncol(deer)):1])
  snow<-unname(snow[,(ncol(snow)):1])
  mei<-unname(mei[,(ncol(mei)):1])
  
  #***solve for the innovations, do it for all years for which it's possible
  
  #the actual deer numbers and corresponding years
  actual<-deer
  actual_yr<-deeryr
  
  #the numbers predicted from the AR terms for deer, and corresponding years
  if (ARo_deer>0)
  {
    ARparms_deer<-parms[1:ARo_deer]
    parms<-parms[-(1:ARo_deer)]
    pred_deer<-matrix(NA,nrow(deer),ncol(deer)-ARo_deer)
    pred_deer_yr<-deeryr[1:(length(deeryr)-ARo_deer)]
    for (counter in 1:(ncol(deer)-ARo_deer))
    {
      h1<-deer[,(counter+1):(counter+ARo_deer),drop=FALSE]
      h2<-matrix(ARparms_deer,nrow(deer),ARo_deer,byrow=TRUE)
      pred_deer[,counter]<-apply(FUN=sum,X=h1*h2,MARGIN=1)
    }
  }else
  {
    pred_deer<-matrix(0,nrow(deer),ncol(deer))
    pred_deer_yr<-deeryr
  }
  
  #the numbers predicted from the snow terms, and corresponding years
  if (ARo_snow>-1)
  {
    ARparms_snow<-parms[1:(ARo_snow+1)]
    parms<-parms[-(1:(ARo_snow+1))]
    pred_snow<-matrix(NA,nrow(deer),ncol(deer)-ARo_snow)
    pred_snow_yr<-deeryr[1:(length(deeryr)-ARo_snow)]
    for (counter in 1:(ncol(deer)-ARo_snow))
    {
      h1<-snow[,counter:(counter+ARo_snow),drop=FALSE]
      h2<-matrix(ARparms_snow,nrow(deer),ARo_snow+1,byrow=TRUE)
      pred_snow[,counter]<-apply(FUN=sum,X=h1*h2,MARGIN=1)
    }
  }else
  {
    pred_snow<-matrix(0,nrow(deer),ncol(deer))
    pred_snow_yr<-deeryr
  }
  
  #the numbers predicted from the MEI terms, and corresponding years
  if (ARo_mei>-1)
  {
    ARparms_mei<-parms[1:(ARo_mei+1)]
    parms<-parms[-(1:(ARo_mei+1))]
    pred_mei<-matrix(NA,nrow(deer),ncol(deer)-ARo_mei)
    pred_mei_yr<-deeryr[1:(length(deeryr)-ARo_mei)]
    for (counter in 1:(ncol(deer)-ARo_mei))
    {
      h1<-mei[,counter:(counter+ARo_mei),drop=FALSE]
      h2<-matrix(ARparms_mei,nrow(deer),ARo_mei+1,byrow=TRUE)
      pred_mei[,counter]<-apply(FUN=sum,X=h1*h2,MARGIN=1)
    }
  }else
  {
    pred_mei<-matrix(0,nrow(deer),ncol(deer))
    pred_mei_yr<-deeryr
  }
  
  #now add up all the components of the prediction and subtract from the actual to get the
  #innovations, but only for the years controlled by fixedminyr
  minyr<-max(min(actual_yr),min(pred_deer_yr),min(pred_snow_yr),min(pred_mei_yr)) #the min year for which it would be possible to get the innovation
  if (minyr>fixedminyr){stop("Error in logLik_class1: lags too large for fixedminyr")}
  minyr<-fixedminyr
  resids<-actual[,which(actual_yr>=minyr)]-pred_deer[,which(pred_deer_yr>=minyr)]-
    pred_snow[,which(actual_yr>=minyr)]-pred_mei[,which(actual_yr>=minyr)]
  
  #***get the ln likelihood using the innovations and return - this is where model assumptions 
  #about the innovations come in
  sd_innov<-parms[1]
  cor_innov<-parms[2]
  sig<-matrix(cor_innov*sd_innov^2,nrow(resids),nrow(resids))
  diag(sig)<-sd_innov^2
  myfun<-function(x){return(mvtnorm::dmvnorm(x,mean=rep(0,length(x)),sigma=sig,log=TRUE))}
  
  h<-sum(apply(FUN=myfun,MARGIN=2,X=resids))
  if (!innovs)
  {
    return(h)
  }else
  {
    return(list(innovs=resids[,ncol(resids):1],logLik=h)) #not forgetting to turn the innovs back around in time
  }
}

#for line-by-line testing, 1
#ARo_deer<-1
#ARo_snow<-1
#ARo_mei<-1
#parms<-c(.9,1,1,1,1,.5,.4)
#fixedminyr<-1986

#for line-by-line testing, 2
#ARo_deer<-2
#ARo_snow<-1
#ARo_mei<-1
#parms<-c(.9,.2,1,1,1,1,.5,.4)
#fixedminyr<-1982

#just try calling this new likelihood function a few times
#logLik_class1(parms=c(.9,.2,1,1,1,1,.5,.4),
#              ARo_deer=2,
#              ARo_snow=1,
#              ARo_mei=1,
#              fixedminyr=1986,
#              deer=deer,
#              snow=snow,
#              mei=mei,
#              deeryr=deeryr)
#logLik_class1(parms=c(.9,.2,.1,1,1,.5,1,1,.4,.5,.4),
#              ARo_deer=3,
#              ARo_snow=2,
#              ARo_mei=2,
#              fixedminyr=1986,
#              deer=deer,
#              snow=snow,
#              mei=mei,
#              deeryr=deeryr)

#***now write a function to automate model comparison across different choices of lags

getAICgivenlags<-function(ARo_deer,ARo_snow,ARo_mei,fixedminyr,numoptims)
{
  #get a startpar using regression
  maxlag<-max(ARo_deer,ARo_snow,ARo_mei)
  if (deeryr[maxlag+1]>fixedminyr){stop("Error in getAICgivenlags: lags too big for fixedminyr")}
  maxlag<-fixedminyr-min(deeryr) #this is to make different lags use the same data 
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
  regform<-paste0(regform,"(-1)")
  lmres<-lm(formula=as.formula(regform),data=regdat)
  coefs<-coef(lmres)
  startpar<-c(unname(coefs),.5,.2) #random choices for the last two parameters
  
  #now do an optimization, prep to do several more, and save the results from the first one
  print(paste0("Optimization ",1," of ",numoptims))
  optres1<-optim(par=startpar,
                 fn=logLik_class1,
                 method="Nelder-Mead",
                 control=list(trace=0,fnscale=-1,maxit=10000),
                 ARo_deer=ARo_deer,
                 ARo_snow=ARo_snow,
                 ARo_mei=ARo_mei,
                 fixedminyr=fixedminyr,
                 deer=deer,
                 snow=snow,
                 mei=mei,
                 deeryr=deeryr)
  alloptres<-matrix(NA,numoptims,2+length(optres1$par))
  colnames(alloptres)<-c("convergence",'value',paste0(rep("par",length(optres1$par)),1:length(optres1$par)))
  alloptres[1,1]<-optres1$convergence
  alloptres[1,2]<-optres1$value
  alloptres[1,3:(dim(alloptres)[2])]<-optres1$par
  
  #now try optimizing from other start locations
  for (counter in 2:numoptims)
  {
    print(paste0("Optimization ",counter," of ",numoptims))
    startpar_p<-startpar*runif(length(startpar),min=.5,max=1.5)
    optresn<-optim(par=startpar_p,
                   fn=logLik_class1,
                   method="Nelder-Mead",
                   control=list(trace=0,fnscale=-1,maxit=10000),
                   ARo_deer=ARo_deer,
                   ARo_snow=ARo_snow,
                   ARo_mei=ARo_mei,
                   fixedminyr=fixedminyr,
                   deer=deer,
                   snow=snow,
                   mei=mei,
                   deeryr=deeryr)
    alloptres[counter,1]<-optresn$convergence
    alloptres[counter,2]<-optresn$value
    alloptres[counter,3:(dim(alloptres)[2])]<-optresn$par
  }
  
  return(list(alloptres=alloptres,numparams=length(startpar)))
}

#***now iterate through various combinations of lags to find out which has the lowest AIC,
#in order to pick a model

#to make it possible to use mclapply in the parallel package
getAICgivenlags_wrapper<-function(x)
{
  return(c(ARo_deer=x$ARo_deer,
           ARo_snow=x$ARo_snow,
           ARo_mei=x$ARo_mei,
           fixedminyr=x$fixedminyr,
           numoptims=x$numoptims,
           getAICgivenlags(x$ARo_deer,x$ARo_snow,x$ARo_mei,x$fixedminyr,x$numoptims)))
}

#set up the argument list
numoptims<-5
ARo_deer_vals<-c(0,1,2,3,4,5)
ARo_snow_vals<-c(-1,0,1,2,3,4)
ARo_mei_vals<-c(-1,0,1,2,3,4)
fixedminyr<-min(deeryr)+max(ARo_deer_vals,ARo_snow_vals,ARo_mei_vals)
arglist<-list()
counter<-1
for (ARo_deer in ARo_deer_vals)
{
  for (ARo_snow in ARo_snow_vals)
  {
    for (ARo_mei in ARo_mei_vals)
    {
      arglist[[counter]]<-list(ARo_deer=ARo_deer,
                               ARo_snow=ARo_snow,
                               ARo_mei=ARo_mei,
                               fixedminyr=fixedminyr,numoptims=numoptims)
      counter<-counter+1
    }
  }
}

allres<-parallel::mclapply(arglist,getAICgivenlags_wrapper,mc.cores=11)

dars<-length(allres)
allres_summary<-data.frame(ARo_deer=NA*numeric(dars),ARo_snow=NA*numeric(dars),ARo_mei=NA*numeric(dars),
                   AIC=NA*numeric(dars))
datapts<-length(deeryr[deeryr>=fixedminyr])
for (counter in 1:dars)
{
  allres_summary$ARo_deer[counter]<-allres[[counter]]$ARo_deer
  allres_summary$ARo_snow[counter]<-allres[[counter]]$ARo_snow
  allres_summary$ARo_mei[counter]<-allres[[counter]]$ARo_mei
  allres_summary$AIC[counter]<-2*allres[[counter]]$numparams-2*max(allres[[counter]]$alloptres[,"value"])
  allres_summary$BIC[counter]<-allres[[counter]]$numparams*log(datapts)-2*max(allres[[counter]]$alloptres[,"value"])
}
allres_summary<-cbind(allres_summary,DeltaAIC=allres_summary$AIC-min(allres_summary$AIC))
allres_summary<-cbind(allres_summary,DeltaBIC=allres_summary$BIC-min(allres_summary$BIC))
allres_summary
allres_summary[allres_summary$DeltaAIC<=4,]
allres_summary[allres_summary$DeltaBIC<=4,]

#save.image(file="OvernightRun20200922v02.RData")
#load(file="OvernightRun20200922v02.RData") #AR0_deer went 0:5, AR0_snow and mei went -1:4
