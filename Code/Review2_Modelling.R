#This script is to do some ML fitting of ARMA-type models to deer data, along with some prelims. Reuman.

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

#***
#Make some plots of the cleaned, transformed data
#***

plot(deeryr,deer[1,],type="l",ylim=range(deer),main="Deer, cleaned and transformed",
     xlab="Year",ylab="Transformed deer")
for (counter in 2:(dim(deer)[1]))
{
  lines(deeryr,deer[counter,],type="l")
}

plot(deeryr,snow[1,],type="l",ylim=range(snow),main="Snow depth, cleaned and transformed",
     xlab="Year",ylab="Transformed snow depth")
for (counter in 2:(dim(snow)[1]))
{
  lines(deeryr,snow[counter,],type="l")
}

plot(deeryr,mei[1,],type="l",ylim=range(mei),main="MEI, cleaned and transformed",
     xlab="Year",ylab="Transformed MEI")

#***
#make some autocorrelation diagrams
#***

#For making the autocorrelation function that is appropriate in this case
#
#d        A matrix with time series in the rows
#maxlag   The maximum lag to use
#
#Output - a vector of length maxlag
#
autocorrfun<-function(d,maxlag)
{
  dd2<-dim(d)[2]
  
  res<-NA*numeric(maxlag)
  for (counter in 1:maxlag)
  {
    res[counter]<-mean(diag(cor(t(d[,1:(dd2-counter)]),t(d[,(1+counter):dd2]))))
  }
  
  return(res)
}

#get the autocorrelation function for deer
maxlag<-15
deer_autocor<-autocorrfun(deer,maxlag)

#test against reshuffled data, two types of reshuffling
nsurr<-1000
deer_autocor_s1<-matrix(NA,nsurr,length(deer_autocor))
deer_autocor_s2<-matrix(NA,nsurr,length(deer_autocor))
for (counter in 1:nsurr)
{
  deer_s1<-deer[,sample(dim(deer)[2],dim(deer)[2],replace=TRUE)]
  deer_autocor_s1[counter,]<-autocorrfun(deer_s1,maxlag)
  
  newinds<-matrix(sample(dim(deer)[2],prod(dim(deer)),replace=TRUE),dim(deer)[1],dim(deer)[2])
  deer_s2<-deer
  for (dr in 1:(dim(deer)[1]))
  {
    deer_s2[dr,]<-deer[dr,newinds[dr,]]
  }
  deer_autocor_s2[counter,]<-autocorrfun(deer_s2,maxlag)
}
qs1<-apply(FUN=quantile,X=deer_autocor_s1,MARGIN=2,probs=c(.025,.975))
qs2<-apply(FUN=quantile,X=deer_autocor_s2,MARGIN=2,probs=c(.025,.975))

plot(1:maxlag,deer_autocor,type="b",ylim=range(deer_autocor,qs))
lines(c(1,maxlag),c(0,0),type="l",lty="solid")
lines(1:maxlag,qs1[1,],lty="dotted")
lines(1:maxlag,qs1[2,],lty="dotted")
lines(1:maxlag,qs2[1,],lty="dashed")
lines(1:maxlag,qs2[2,],lty="dashed")

#Note: these are not really appropriate for judging what deer AR lags to use in a model,
#since deer autocorrelation can result from density dependence or from the influence of 
#external drivers. It's not entirely clear to me what this plot is good for, I just made
#it because when one with an ARX model one often starts with an autocorrelation plot. 

#***
#make some cross correlation diagrams
#***

#For making cross correlations functions that are appropriate in this case
#
#d1         A matrix with time series in the rows
#d2         Another such
#maxlag     The maximum lag to use
#
#Output - a data frame with these entries
#lag        -maxlag up to maxlag
#cor        The lagged correlation
#
crosscorrfun<-function(d1,d2,maxlag)
{
  if (any(dim(d1)!=dim(d2)))
  {
    stop("Error in crosscorrfun: nonconformable arguments d1 and d2")
  }
  
  dd2<-dim(d1)[2]
  
  res<-data.frame(lag=(-maxlag):maxlag,cor=NA)
  for (lag in (-maxlag):maxlag)
  {
    if (lag<0)
    {
      res[res$lag==lag,2]<-mean(diag(cor(t(d1[,1:(dd2+lag)]),t(d2[,(-lag+1):dd2]))))
    }
    if (lag>0)
    {
      res[res$lag==lag,2]<-mean(diag(cor(t(d1[,(lag+1):dd2]),t(d2[,1:(dd2-lag)]))))
    }
    if (lag==0)
    {
      res[res$lag==lag,2]<-mean(diag(cor(t(d1),t(d2))))
    }
  }
  
  return(res)
}

#get the cross correlation function for deer and snow
maxlag<-15
deer_snow_crosscor<-crosscorrfun(deer,snow,maxlag)

#get the cross correlation function for deer and mei
maxlag<-15
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
plot(deer_snow_crosscor$lag,deer_snow_crosscor$cor,type="b",ylim=range(deer_snow_crosscor$cor,qs_snow))
lines(range(deer_snow_crosscor$lag),c(0,0),type="l",lty="dashed")
lines(deer_snow_crosscor$lag,qs_snow[1,],lty="dotted")
lines(deer_snow_crosscor$lag,qs_snow[2,],lty="dotted")

plot(deer_mei_crosscor$lag,deer_mei_crosscor$cor,type="b",ylim=range(deer_mei_crosscor$cor,qs_mei))
lines(range(deer_mei_crosscor$lag),c(0,0),type="l",lty="dashed")
lines(deer_mei_crosscor$lag,qs_mei[1,],lty="dotted")
lines(deer_mei_crosscor$lag,qs_mei[2,],lty="dotted")

#***
#Now do some fitting, class 1
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
#fixedminyr   Used for making sure the same amount of data are used for different lags
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
  
  #***get the ln likelihood using the innovations and return
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
    return(list(innovs=resids,logLik=h))
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
  #get startpar using regression
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
    startpar_p<-startpar*runif(length(startpar),min=.07,max=1.3)
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
  #set.seed(101)
  return(c(ARo_deer=x$ARo_deer,
           ARo_snow=x$ARo_snow,
           ARo_mei=x$ARo_mei,
           fixedminyr=x$fixedminyr,
           numoptims=x$numoptims,
           getAICgivenlags(x$ARo_deer,x$ARo_snow,x$ARo_mei,x$fixedminyr,x$numoptims)))
}

#set up the argument list
fixedminyr<-1985
numoptims<-5
ARo_deer_vals<-c(0,1,2)
ARo_snow_vals<-c(-1,0,1,2,3,4)
ARo_mei_vals<-c(-1,0,1,2,3,4)
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

allres<-parallel::mclapply(arglist,getAICgivenlags_wrapper,mc.cores=5)

dars<-length(allres)
allres_summary<-data.frame(ARo_deer=NA*numeric(dars),ARo_snow=NA*numeric(dars),ARo_mei=NA*numeric(dars),
                   AIC=NA*numeric(dars))
for (counter in 1:dars)
{
  allres_summary$ARo_deer[counter]<-allres[[counter]]$ARo_deer
  allres_summary$ARo_snow[counter]<-allres[[counter]]$ARo_snow
  allres_summary$ARo_mei[counter]<-allres[[counter]]$ARo_mei
  allres_summary$AIC[counter]<-2*allres[[counter]]$numparams-2*max(allres[[counter]]$alloptres[,"value"])
}
allres_summary<-cbind(allres_summary,DeltaAIC=allres_summary$AIC-min(allres_summary$AIC))
allres_summary
allres_summary[allres_summary$DeltaAIC<=4,]
save.image(file="OvernightRun20200807v02.RData")




#***Sim the AIC-best model and see if the results look like the deer when plotted. Use the actual
#snow and MEI time series, and the actual initial conditions.



#***See if the AIC-best model has good properties of the residuals, inc. temporal independence,
#normal marginals of equal variance, etc.



#***See if synchrony of sims of the AIC-best model look like synchrony of deer, using my spectral and 
#maybe wavelet measures of synchrony. 
#Maybe also see if you can change something appropriate and make the synchrony of deer change. One
#thing you might want to change is to simply eliminate the influence of snow or mei (drop from the
#model), or else replace snow or mei with asynchronous versions of themselves using asynchronous
#Fourier surrogates. Not clear either of these is the right thing to do. You could drop both from the
#model, at once, but this seems like it will clearly result in some synchrony, but not timescale-specific.
#So maybe that makes sense to do.



#***Can analyze the new model semi-analytically (some terms of the theory can be computed directly
#from model coefficients, but you have to estimate spectral quantities of the noise) to illustrate
#how the synchrony passes from the noise to the deer, and hopefully this will parallel the direct
#empirical analysis (using the new theory) of the same thing. Actually can't do that because the model
#has two drivers.



#***Some things that maybe still need to be done, some currently vague, the most important ones 
#explained in greater depth above.
#1) see if the AIC-best model is a good fit using a bootstrapped ML approach
#2) see if the AIC-best model has good properties of the residuals, inc. temporal independence,
#normal marginals of equal variance, etc.
#3) ***see if sims of the AIC-best model look like sims of deer
#4) see if synchrony of sims of the AIC-best model look like sims of deer, and then see if you can
#change something appropriate and make the synchrony of deer change
#5) can analyze the new model semi-analytically (some terms of the theory can be computed directly
#from model coefficients, but you have to estimate spectral quantities of the noise) to illustrate
#how the synchrony passes from the noise to the deer, and hopefully this will parallel the direct
#empirical analysis (using the new theory) of the same thing.
#DONE 6) Examine parameters of the AIC-best model
#RUNNING 7) Do we need to go up to an MEI lag of 5, since a lag of 4 is best so far? Run that over night?