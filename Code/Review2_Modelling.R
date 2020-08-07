#***
#Thinking about models
#***

#The purpose here is to think about using ARMA models of deer driven by snow and MEI.
#What models? To what purpose?

#What models:
#1) In each location model deer with an TBD number of AR lags (though the above results may shed some light) and
#using noise drivers snow depth and MEI. You'll probably also need another error term, but not sure what kinds of 
#assumptions to make about it. 

#To what purpose:
#1) Broadly, to substantiate the argument that Moran effects from snow and MEI are causing the timescale-specific
#synchrony we see. 
#2) It would be nice if we could set up a model and argue it's a good fit, and then simulate the model and the sims
#have the same or similar patterns of timescale-specific synchrony as the real deer.
#3) The referee seems to want us to establish the drivers and then show they give the synchrony patterns in deer. He 
#even wants us to explain the lack of synchrony in deer at long timescales. 
#4) What if we cannot do anything about the long timescales, but we can get a model using snow and MEI and it generates 
#the right synchrony at 3-7-year timescales and then when you unsynchronize snow and simulate, you loose the deer 
#synchrony at 3-4 and when you unsynchronize MEI you loose the deer synchrony at 4-7?
#5) If you have an ARMA model driven by some noises, what's the coherence between noise and population and how does it
#depend on the ARMA coefficients?



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
#deer         Deer data, pre-treated as above
#snow         Snow data, pre-treated as above
#mei          MEI data, pre-treated as above
#deeryr       Years for the data, as above
#innovs       TRUE to also return the innovations. Default FALSE.
#
#Output - the ln likelihood of the model
#
logLik_class1<-function(parms,ARo_deer,ARo_snow,ARo_mei,deer,snow,mei,deeryr,innovs=FALSE)
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
  
  #***solve for the innovations
  
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
  #innovations
  minyr<-max(min(actual_yr),min(pred_deer_yr),min(pred_snow_yr),min(pred_mei_yr))
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

#for line-by-line testing, 2
#ARo_deer<-2
#ARo_snow<-1
#ARo_mei<-1
#parms<-c(.9,.2,1,1,1,1,.5,.4)

#just try calling this new likelihood function a few times
#logLik_class1(parms=c(.9,.2,1,1,1,1,.5,.4),
#              ARo_deer=2,
#              ARo_snow=1,
#              ARo_mei=1,
#              deer=deer,
#              snow=snow,
#              mei=mei,
#              deeryr=deeryr)
#logLik_class1(parms=c(.9,.2,.1,1,1,.5,1,1,.4,.5,.4),
#              ARo_deer=3,
#              ARo_snow=2,
#              ARo_mei=2,
#              deer=deer,
#              snow=snow,
#              mei=mei,
#              deeryr=deeryr)

#***now do some optimization, starting from coordinates determined by regression, using lags of
#1 for deer (AR(1) model for deer), 2 for snow (so this year's snow, last year's, and the 
#year before can influence deer), and 3 for MEI

#the lags
#ARo_deer<-1
#ARo_snow<-2
#ARo_mei<-3

#a random guess, just to get things working
#startpar<-c(.9,1,1,1,1,1,1,1,.5,.2) 
#
#optres1<-optim(par=startpar,
#               fn=logLik_class1,
#               method="Nelder-Mead",
#               control=list(trace=10000,fnscale=-1,maxit=1000),
#               ARo_deer=ARo_deer,
#               ARo_snow=ARo_snow,
#               ARo_mei=ARo_mei,
#               deer=deer,
#               snow=snow,
#               mei=mei,
#               deeryr=deeryr)
#optres1$value
#logLik_class1(parms=optres1$par,
#              ARo_deer=ARo_deer,
#              ARo_snow=ARo_snow,
#              ARo_mei=ARo_mei,
#              deer=deer,
#              snow=snow,
#              mei=mei,
#              deeryr=deeryr)

#now get startpar using regression
#regdat<-data.frame(deer_resp=as.vector(deer[,4:(ncol(deer))]),
#                   deer_pred_1=as.vector(deer[,3:(ncol(deer)-1)]),
#                   snow_pred_0=as.vector(snow[,4:(ncol(deer))]),
#                   snow_pred_1=as.vector(snow[,3:(ncol(deer)-1)]),
#                   snow_pred_2=as.vector(snow[,2:(ncol(deer)-2)]),
#                   mei_pred_0=as.vector(mei[,4:(ncol(deer))]),
#                   mei_pred_1=as.vector(mei[,3:(ncol(deer)-1)]),
#                   mei_pred_2=as.vector(mei[,2:(ncol(deer)-2)]),
#                   mei_pred_3=as.vector(mei[,1:(ncol(deer)-3)]))
#lmres<-lm(deer_resp~deer_pred_1+snow_pred_0+snow_pred_1+snow_pred_2+mei_pred_0+mei_pred_1+mei_pred_2+mei_pred_3-1,
#   data=regdat)
#coefs<-coef(lmres)
#startpar<-c(unname(coefs),.5,.2)

#now do the optimization
#optres1<-optim(par=startpar,
#      fn=logLik_class1,
#      method="Nelder-Mead",
#      control=list(trace=10000,fnscale=-1,maxit=10000),
#      ARo_deer=ARo_deer,
#      ARo_snow=ARo_snow,
#      ARo_mei=ARo_mei,
#      deer=deer,
#      snow=snow,
#      mei=mei,
#      deeryr=deeryr)
#optres1$convergence
#optres1$value
#optres1$par

#now try optimizing from other start locations
#set.seed(101)

#numoptims<-9
#alloptres1<-matrix(NA,numoptims+1,2+length(optres1$par))
#colnames(alloptres1)<-c("convergence",'value',paste0(rep("par",length(optres1$par)),1:length(optres1$par)))
#alloptres1[1,1]<-optres1$convergence
#alloptres1[1,2]<-optres1$value
#alloptres1[1,3:12]<-optres1$par
#for (counter in 1:numoptims)
#{
#  print(paste0("Optimization ",counter," of ",numoptims))
#  startpar_p<-startpar*runif(length(startpar),min=.07,max=1.3)
#  optresn<-optim(par=startpar_p,
#                 fn=logLik_class1,
#                 method="Nelder-Mead",
#                 control=list(trace=0,fnscale=-1,maxit=10000),
#                 ARo_deer=ARo_deer,
#                 ARo_snow=ARo_snow,
#                 ARo_mei=ARo_mei,
#                 deer=deer,
#                 snow=snow,
#                 mei=mei,
#                 deeryr=deeryr)
#  alloptres1[counter+1,1]<-optresn$convergence
#  alloptres1[counter+1,2]<-optresn$value
#  alloptres1[counter+1,3:12]<-optresn$par
#}

#get the AIC
#AIC1<-2*length(startpar)-2*max(alloptres1[,"value"])

#***now try lags of 2 for deer (AR(2) model for deer), 2 for snow (so this year's snow, last 
#year's, and the year before can influence deer), and 3 for MEI

#the lags
#ARo_deer<-2
#ARo_snow<-2
#ARo_mei<-3

#now get startpar using regression
#regdat<-data.frame(deer_resp=as.vector(deer[,4:(ncol(deer))]),
#                   deer_pred_1=as.vector(deer[,3:(ncol(deer)-1)]),
#                   deer_pred_2=as.vector(deer[,2:(ncol(deer)-2)]),
#                   snow_pred_0=as.vector(snow[,4:(ncol(deer))]),
#                   snow_pred_1=as.vector(snow[,3:(ncol(deer)-1)]),
#                   snow_pred_2=as.vector(snow[,2:(ncol(deer)-2)]),
#                   mei_pred_0=as.vector(mei[,4:(ncol(deer))]),
#                   mei_pred_1=as.vector(mei[,3:(ncol(deer)-1)]),
#                   mei_pred_2=as.vector(mei[,2:(ncol(deer)-2)]),
#                   mei_pred_3=as.vector(mei[,1:(ncol(deer)-3)]))
#lmres<-lm(deer_resp~deer_pred_1+deer_pred_2+snow_pred_0+snow_pred_1+snow_pred_2+mei_pred_0+mei_pred_1+mei_pred_2+mei_pred_3-1,
#          data=regdat)
#coefs<-coef(lmres)
#startpar<-c(unname(coefs),.5,.2)

#now do the optimization
#optres1<-optim(par=startpar,
#               fn=logLik_class1,
#               method="Nelder-Mead",
#               control=list(trace=10000,fnscale=-1,maxit=10000),
#               ARo_deer=ARo_deer,
#               ARo_snow=ARo_snow,
#               ARo_mei=ARo_mei,
#               deer=deer,
#               snow=snow,
#               mei=mei,
#               deeryr=deeryr)
#optres1$convergence
#optres1$value
#optres1$par

#now try optimizing from other start locations
#set.seed(101)

#numoptims<-9
#alloptres2<-matrix(NA,numoptims+1,2+length(optres1$par))
#colnames(alloptres2)<-c("convergence",'value',paste0(rep("par",length(optres1$par)),1:length(optres1$par)))
#alloptres2[1,1]<-optres1$convergence
#alloptres2[1,2]<-optres1$value
#alloptres2[1,3:13]<-optres1$par
#for (counter in 1:numoptims)
#{
#  print(paste0("Optimization ",counter," of ",numoptims))
#  startpar_p<-startpar*runif(length(startpar),min=.07,max=1.3)
#  optresn<-optim(par=startpar_p,
#                 fn=logLik_class1,
#                 method="Nelder-Mead",
#                 control=list(trace=0,fnscale=-1,maxit=10000),
#                 ARo_deer=ARo_deer,
#                 ARo_snow=ARo_snow,
#                 ARo_mei=ARo_mei,
#                 deer=deer,
#                 snow=snow,
#                 mei=mei,
#                 deeryr=deeryr)
#  alloptres2[counter+1,1]<-optresn$convergence
#  alloptres2[counter+1,2]<-optresn$value
#  alloptres2[counter+1,3:13]<-optresn$par
#}

#get the AIC
#AIC2<-2*length(startpar)-2*max(alloptres2[,"value"])

#***now write a function to automate things for different choices of lags

getAICgivenlags<-function(ARo_deer,ARo_snow,ARo_mei)
{
  #get startpar using regression
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
  regform<-paste0(regform,"(-1)")
  lmres<-lm(formula=as.formula(regform),data=regdat)
  coefs<-coef(lmres)
  startpar<-c(unname(coefs),.5,.2) #random choices for the last two parameters
  
  #now do an optimization, prep to do several, and save the results from the first one
  numoptims<-5
  print(paste0("Optimization ",1," of ",numoptims))
  optres1<-optim(par=startpar,
                 fn=logLik_class1,
                 method="Nelder-Mead",
                 control=list(trace=0,fnscale=-1,maxit=10000),
                 ARo_deer=ARo_deer,
                 ARo_snow=ARo_snow,
                 ARo_mei=ARo_mei,
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

#***now call the function on the two cases done manually above, for comparison with the above
#results, as a test

#set.seed(101)
#res1_f<-getAICgivenlags(ARo_deer=1,ARo_snow=2,ARo_mei=3)
#alloptres1_f<-res1_f$alloptres
#AIC1_f<-2*res1_f$numparams-2*max(res1_f$alloptres[,"value"])

#alloptres1
#alloptres1_f
#max(abs(alloptres1-alloptres1_f))
#AIC1
#AIC1_f
#AIC1-AIC1_f

#set.seed(101)
#res2_f<-getAICgivenlags(ARo_deer=2,ARo_snow=2,ARo_mei=3)
#alloptres2_f<-res2_f$alloptres
#AIC2_f<-2*res2_f$numparams-2*max(res2_f$alloptres[,"value"])

#alloptres2
#alloptres2_f
#max(abs(alloptres2-alloptres2_f))
#AIC2
#AIC2_f
#AIC2-AIC2_f

#***now iterate through various combinations of lags to find out which has the lowest AIC,
#in order to pick a model

allres<-list()
allres_summary<-data.frame(ARo_deer=NA*numeric(27),ARo_snow=NA*numeric(27),ARo_mei=NA*numeric(27),
                   AIC=NA*numeric(27))
counter<-1
for (ARo_deer in c(0,1,2))
{
  print(paste0("ARo_deer=",ARo_deer))
  for (ARo_snow in c(-1,0,1,2,3))
  {
    print(paste0("  ARo_snow=",ARo_snow))
    for (ARo_mei in c(-1,0,1,2,3,4,5))
    {
      print(paste0("    ARo_mei=",ARo_mei))
      allres[[counter]]<-getAICgivenlags(ARo_deer,ARo_snow,ARo_mei)
      allres_summary[counter,]<-c(ARo_deer,ARo_snow,ARo_mei,
                                  2*allres[[counter]]$numparams-2*max(allres[[counter]]$alloptres[,"value"]))
      counter<-counter+1
    }
  }
}
allres_summary<-cbind(allres_summary,delta_AIC=allres_summary$AIC-min(allres_summary$AIC))
allres_summary
save.image(file="OvernightRun20200806.RData")

#***Sim the AIC-best model and see if the results look like the deer when plotted. Use the actual
#snow and MEI time series, and the actual initial conditions.



#***See if the AIC-best model has good properties of the residuals, inc. temporal independence,
#normal marginals of equal variance, etc.



#***See if synchrony of sims of the AIC-best model look like synchrony of deer, using my spectral and 
#maybe wavelet measures of synchrony. 
#Maybe also see if you can change something appropriate and make the synchrony of deer change. One
#thing you might want to change is to simply eliminate the influence of snow or mei (drop from the
#model), or else replace snow or mei with asynchronous versions of themselves using asynchronous
#Fourier surrogates. Not clear either of these is the right thing to do.



#***Can analyze the new model semi-analytically (some terms of the theory can be computed directly
#from model coefficients, but you have to estimate spectral quantities of the noise) to illustrate
#how the synchrony passes from the noise to the deer, and hopefully this will parallel the direct
#empirical analysis (using the new theory) of the same thing.



#***Some things that maybe still need to be done, some currently vague, the most important ones 
#explained in greater depth below.
#1) see if the AIC-best model is a good fit using a bootstrapped ML approach
#2) see if the AIC-bets model has good properties of the residuals, inc. temporal independence,
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