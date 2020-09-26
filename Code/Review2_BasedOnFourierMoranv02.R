#This script does some calculations based on the Fourier Moran approach developed in August 2020, for deer. Reuman

theseed<-101

#***
#functions
#***

source("Code/SpectralTools.R")

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
#Apply the Fourier Moran theorem to deer and snow
#***

#This spectral matrix should have all the spectral and cross-spectral matrices which are needed
x<-rbind(deer,snow)
sres<-myspecmatbrill(x,detrend=FALSE,BiasVariance=0.75)
freq<-sres$freq
sres<-sres$specmat

#Now make surrogates and get spectral information about them 
nsurrog<-1000
set.seed(theseed)
snow_s<-wsyn::surrog(snow,nsurrog,"fft",TRUE)
sres_s<-list()
for (counter in 1:nsurrog)
{
  print(paste0("Making surrogate ",counter," of ",nsurrog))
  x<-rbind(deer,snow_s[[counter]])
  sres_s[[counter]]<-myspecmatbrill(x,detrend=TRUE,BiasVariance=0.75)$specmat  
}

#Numeric difficulties prevent doing the desired analysis for all locations. However, we can do it for 
#a subset. This function does that.
#
#Args
#locstouse          A vector of numbers between 1 and 60 specifying the locations to use.
#The function also uses some of the other variables created above.
#
#Output - a list with these named entries (see the math description for variable names, which correspond
#to the notation used there)
#rho_wOmegai        Total population synchrony for the locations in Omegai=locstouse. This is a vector of 
#                     length equal to that of freq (above), and is a function of frequency.
#rho_wOmegai_eps1   Portion of rho_wOmegai explained by epsilon1, the measured noise. This is a vector of 
#                     length equal to that of freq (above), and is a function of frequency.
#rho_wOmegai_delj1  Portion of rho_wOmegai explained by delta_j^(1), the jth surrogate of epsilon^(1). 
#                     This is an nsurrog by length(freq) matrix. 
#
do_analysis_some_locs<-function(locstouse,ploton=FALSE)
{
  #extract the needed spectral and cross-spectral matrices for the analysis of the real data
  numlocs<-length(locstouse)
  dd1<-dim(deer)[1]
  Sdeerdeer<-sres[locstouse,locstouse,]
  Ssnowsnow<-sres[locstouse+dd1,locstouse+dd1,]
  Sdeersnow<-sres[locstouse,locstouse+dd1,]
  
  #compute the actual deer synchrony for these sites, rho_wOmegai
  rho_wOmegai<-0*Sdeerdeer[1,1,] #receptacle that will store the result
  for (ci in 1:numlocs)
  {
    for (cj in 1:numlocs)
    {
      if (ci != cj)
      {
        rho_wOmegai<-rho_wOmegai+Sdeerdeer[ci,cj,]
      }
    }
  }
  rho_wOmegai<-rho_wOmegai/(numlocs*(numlocs-1))

  #compute the deer synchrony statistically explainable by snow for these sites, rho_wOmegai_eps1
  FMTmat<-NA*Sdeerdeer #receptacle for the matrix M in the Fourier Moran theorem
  for (counter in 1:length(freq))
  {
    if (kappa(Ssnowsnow[,,counter],exact=TRUE,method="qr")>1e9 || kappa(Ssnowsnow[,,counter],exact=TRUE,method="direct")>1e9)
    {
      stop("Error in do_analysis_some_locs: problem 1 finding the Fourier Moran matrix")
    }
    h<-solve(Ssnowsnow[,,counter],Conj(t(Sdeersnow[,,counter])))
    if (any(Mod(Ssnowsnow[,,counter] %*% h - Conj(t(Sdeersnow[,,counter])))>1e-8))
    {
      stop("Error in do_analysis_some_locs: problem 2 finding the Fourier Moran matrix")
    }
    FMTmat[,,counter]<-Sdeersnow[,,counter] %*% h
  }
  rho_wOmegai_eps1<-0*rho_wOmegai #receptacle for the result
  for (ci in 1:numlocs)
  {
    for (cj in 1:numlocs)
    {
      if (ci != cj)
      {
        rho_wOmegai_eps1<-rho_wOmegai_eps1+FMTmat[ci,cj,]
      }
    }
  }
  rho_wOmegai_eps1<-rho_wOmegai_eps1/(numlocs*(numlocs-1))
  
  #compute the deer synchrony explainable by surrogate snow data, rho_wOmegai_delj1
  rho_wOmegai_delj1<-matrix(complex(0,0),nsurrog,length(freq))
  for (scounter in 1:nsurrog)
  {
    Ssnowsnow_s<-sres_s[[scounter]][locstouse+dd1,locstouse+dd1,]
    Sdeersnow_s<-sres_s[[scounter]][locstouse,locstouse+dd1,]
    
    FMTmat<-NA*Sdeerdeer #receptacle for the matrix M in the Fourier Moran theorem
    for (counter in 1:length(freq))
    {
      if (kappa(Ssnowsnow_s[,,counter],exact=TRUE,method="qr")>1e9 || kappa(Ssnowsnow_s[,,counter],exact=TRUE,method="direct")>1e9)
      {
        stop("Error in do_analysis_some_locs: problem 1 finding the Fourier Moran matrix, surrogates")
      }
      h<-solve(Ssnowsnow_s[,,counter],Conj(t(Sdeersnow_s[,,counter])))
      if (any(Mod(Ssnowsnow_s[,,counter] %*% h - Conj(t(Sdeersnow_s[,,counter])))>1e-8))
      {
        stop("Error in do_analysis_some_locs: problem 2 finding the Fourier Moran matrix, surrogates")
      }
      FMTmat[,,counter]<-Sdeersnow_s[,,counter] %*% h
    }
    for (ci in 1:numlocs)
    {
      for (cj in 1:numlocs)
      {
        if (ci != cj)
        {
          rho_wOmegai_delj1[scounter,]<-rho_wOmegai_delj1[scounter,]+FMTmat[ci,cj,]
        }
      }
    }
  }
  rho_wOmegai_delj1<-rho_wOmegai_delj1/(numlocs*(numlocs-1))
  
  #check for non-zero imaginary parts and if you don't find them take the real part
  if (any(Im(rho_wOmegai)>1e-10))
  {
    stop("Error in do_analysis_some_locs: rho_wOmegai was not real")
  }
  if (any(Im(rho_wOmegai_eps1)>1e-10))
  {
    stop("Error in do_analysis_some_locs: rho_wOmegai_eps1 was not real")
  }
  if (any(Im(rho_wOmegai_delj1)>1e-10))
  {
    stop("Error in do_analysis_some_locs: rho_wOmegai_delj1 was not real")
  }
  rho_wOmegai<-Re(rho_wOmegai)
  rho_wOmegai_eps1<-Re(rho_wOmegai_eps1)
  rho_wOmegai_delj1<-Re(rho_wOmegai_delj1)
  
  #plot the quantities just computed on the same axes, if desired
  if (ploton==TRUE)
  {
    rho_wOmegai_delj1_quant<-apply(FUN=quantile,MARGIN=2,X=rho_wOmegai_delj1,probs=c(.025,.975))
    ylimits<-range(rho_wOmegai,rho_wOmegai_eps1,rho_wOmegai_delj1_quant)
    plot(freq,rho_wOmegai,type="l",ylim=ylimits,
         xlab="Frequency",ylab="Synchrony")
    lines(freq,rho_wOmegai_eps1,type="l",lty="dashed",col="red")
    lines(freq,rho_wOmegai_delj1_quant[1,],lty="dotted",col="red")
    lines(freq,rho_wOmegai_delj1_quant[2,],lty="dotted",col="red")
    lines(range(freq),c(0,0),type="l",lty="dotted")
    lines(rep(1/7,2),ylimits,type="l",lty="dotted")
    lines(rep(1/3,2),ylimits,type="l",lty="dotted")
    lines(rep(1/4,2),ylimits,type="l",lty="dotted")
  }
  
  return(list(rho_wOmegai=rho_wOmegai,rho_wOmegai_eps1=rho_wOmegai_eps1,rho_wOmegai_delj1=rho_wOmegai_delj1))
}

#now partition the sites into groups of five and do the analysis for each group
set.seed(theseed) #this is so that changing the number of surrogates and other changes to the upstream code
                  #does not also change the location grouping
locsrdmzd<-matrix(sample.int(dim(deer)[1],dim(deer)[1]),5,12)
allres<-list()
for (counter in 1:(dim(locsrdmzd)[2]))
{
  print(paste0("Analyzing location subset ",counter," of ",dim(locsrdmzd)[2]))
  allres[[counter]]<-do_analysis_some_locs(locsrdmzd[,counter],ploton=FALSE)
}

#now average everything across the 12 location groups
rho_wOmega<-allres[[1]]$rho_wOmegai
rho_wOmega_eps1<-allres[[1]]$rho_wOmegai_eps1
rho_wOmega_delj1<-allres[[1]]$rho_wOmegai_delj1
for (counter in 2:length(allres))
{
  rho_wOmega<-rho_wOmega+allres[[counter]]$rho_wOmegai
  rho_wOmega_eps1<-rho_wOmega_eps1+allres[[counter]]$rho_wOmegai_eps1
  rho_wOmega_delj1<-rho_wOmega_delj1+allres[[counter]]$rho_wOmegai_delj1
}
rho_wOmega<-rho_wOmega/length(allres)
rho_wOmega_eps1<-rho_wOmega_eps1/length(allres)
rho_wOmega_delj1<-rho_wOmega_delj1/length(allres)

#now make a plot
rho_wOmega_delj1_quant<-apply(FUN=quantile,MARGIN=2,X=rho_wOmega_delj1,probs=c(.025,.975))
ylimits<-range(rho_wOmega,rho_wOmega_eps1,rho_wOmega_delj1_quant)

tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier_Step2_deersnow_results.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
plot(freq,rho_wOmega,type="l",ylim=ylimits)
mtext("Frequency",side=1,line=2,cex=1)
mtext("Synchrony",side=2,line=2,cex=1)
lines(freq,rho_wOmega_eps1,type="l",lty="dashed",col="red")
lines(freq,rho_wOmega_delj1_quant[1,],lty="dotted",col="red")
lines(freq,rho_wOmega_delj1_quant[2,],lty="dotted",col="red")
lines(rep(1/7,2),ylimits,type="l",lty="dotted")
lines(rep(1/3,2),ylimits,type="l",lty="dotted")
lines(rep(1/4,2),ylimits,type="l",lty="dotted")
dev.off()

#now get fractions explained in each band
fracexplFMT_snow_3to7<-sum(rho_wOmega_eps1[freq>=1/7 & freq<=1/3])/sum(rho_wOmega[freq>=1/7 & freq<=1/3])
fracexplFMT_snow_3to4<-sum(rho_wOmega_eps1[freq>=1/4 & freq<=1/3])/sum(rho_wOmega[freq>=1/4 & freq<=1/3])
fracexplFMT_snow_4to7<-sum(rho_wOmega_eps1[freq>=1/7 & freq<=1/4])/sum(rho_wOmega[freq>=1/7 & freq<=1/4])

#same thing for each of the 12 sets of locations, just for reference
fracexplFMT_snow_3to7_i<-NA*numeric(length(allres))
fracexplFMT_snow_3to4_i<-NA*numeric(length(allres))
fracexplFMT_snow_4to7_i<-NA*numeric(length(allres))
for (counter in 1:length(allres))
{
  fracexplFMT_snow_3to7_i[counter]<-sum(allres[[counter]]$rho_wOmegai_eps1[freq>=1/7 & freq<=1/3])/sum(allres[[counter]]$rho_wOmegai[freq>=1/7 & freq<=1/3])
  fracexplFMT_snow_3to4_i[counter]<-sum(allres[[counter]]$rho_wOmegai_eps1[freq>=1/4 & freq<=1/3])/sum(allres[[counter]]$rho_wOmegai[freq>=1/4 & freq<=1/3])
  fracexplFMT_snow_4to7_i[counter]<-sum(allres[[counter]]$rho_wOmegai_eps1[freq>=1/7 & freq<=1/4])/sum(allres[[counter]]$rho_wOmegai[freq>=1/7 & freq<=1/4])
}

#now do the rank of ranks procedure
rho_wOmega_eps1_fl<-NA*rho_wOmega_eps1
rho_wOmega_delj1_fl<-NA*rho_wOmega_delj1
for (counter in 1:length(freq))
{
  print(paste0("Rank of ranks procedure, frequency ",counter," of ",length(freq)))
  rho_wOmega_eps1_fl[counter]<-sum(rho_wOmega_delj1[,counter]<rho_wOmega_eps1[counter])/nsurrog
  h<-rho_wOmega_delj1[,counter]
  if (length(h)!=length(unique(h))) {stop("Error at location X: ties where there should not be")}
  rh<-rank(h)
  rho_wOmega_delj1_fl[,counter]<-(rh-1)/(nsurrog-1)
}
mnrk3to7<-mean(rho_wOmega_eps1_fl[freq>=1/7 & freq<=1/3])
mnrk3to4<-mean(rho_wOmega_eps1_fl[freq>=1/4 & freq<=1/3])
mnrk4to7<-mean(rho_wOmega_eps1_fl[freq>=1/7 & freq<=1/4])
mnrk3to7_s<-apply(FUN=mean,X=rho_wOmega_delj1_fl[,freq>=1/7 & freq<=1/3],MARGIN=1)
mnrk3to4_s<-apply(FUN=mean,X=rho_wOmega_delj1_fl[,freq>=1/4 & freq<=1/3],MARGIN=1)
mnrk4to7_s<-apply(FUN=mean,X=rho_wOmega_delj1_fl[,freq>=1/7 & freq<=1/4],MARGIN=1)
pvalresFMT_snow_3to7<-1-sum(mnrk3to7>=mnrk3to7_s)/nsurrog
pvalresFMT_snow_3to4<-1-sum(mnrk3to4>=mnrk3to4_s)/nsurrog
pvalresFMT_snow_4to7<-1-sum(mnrk4to7>=mnrk4to7_s)/nsurrog

#An alternative significance testing procedure - tests significance of the actual fractions of synchrony explained in
#each band that we may end up focussing on 
rho_wOmega_eps1_3to7<-sum(rho_wOmega_eps1[freq>=1/7 & freq<=1/3])
rho_wOmega_delj1_3to7<-apply(FUN=sum,X=rho_wOmega_delj1[,freq>=1/7 & freq<=1/3],MARGIN=1)
pvalresFMTalt_snow_3to7<-1-sum(rho_wOmega_eps1_3to7>=rho_wOmega_delj1_3to7)/nsurrog

rho_wOmega_eps1_3to4<-sum(rho_wOmega_eps1[freq>=1/4 & freq<=1/3])
rho_wOmega_delj1_3to4<-apply(FUN=sum,X=rho_wOmega_delj1[,freq>=1/4 & freq<=1/3],MARGIN=1)
pvalresFMTalt_snow_3to4<-1-sum(rho_wOmega_eps1_3to4>=rho_wOmega_delj1_3to4)/nsurrog

rho_wOmega_eps1_4to7<-sum(rho_wOmega_eps1[freq>=1/7 & freq<=1/4])
rho_wOmega_delj1_4to7<-apply(FUN=sum,X=rho_wOmega_delj1[,freq>=1/7 & freq<=1/4],MARGIN=1)
pvalresFMTalt_snow_4to7<-1-sum(rho_wOmega_eps1_4to7>=rho_wOmega_delj1_4to7)/nsurrog

#***
#Apply the Fourier Moran theorem to deer and MEI
#***

#This spectral matrix should have all the spectral and cross-spectral matrices which are needed
mei<-mei[1,,drop=FALSE] #mei is the same in all locations, so it is now convenient to cut it down
x<-rbind(deer,mei)
sres<-myspecmatbrill(x,detrend=FALSE,BiasVariance=0.75)
freq<-sres$freq
sres<-sres$specmat

#Now make surrogates and get spectral information about them 
set.seed(theseed)
mei_s<-wsyn::surrog(mei,nsurrog,"fft",TRUE)
sres_s<-list()
for (counter in 1:nsurrog)
{
  print(paste0("Making surrogate ",counter," of ",nsurrog))
  x<-rbind(deer,mei_s[[counter]])
  sres_s[[counter]]<-myspecmatbrill(x,detrend=TRUE,BiasVariance=0.75)$specmat  
}

#compute the actual deer synchrony, rho_w
dd1<-dim(deer)[1]
Sdeerdeer<-sres[1:dd1,1:dd1,]
rho_w<-0*Sdeerdeer[1,1,] #receptacle that will store the result
for (ci in 1:dd1)
{
  for (cj in 1:dd1)
  {
    if (ci != cj)
    {
      rho_w<-rho_w+Sdeerdeer[ci,cj,]
    }
  }
}
rho_w<-rho_w/(dd1*(dd1-1))

#compute the portion explained by mei
Smeimei<-sres[dd1+1,dd1+1,]
Sdeermei<-sres[1:dd1,dd1+1,]
rho_w_eps1<-0*rho_w
for (ci in 1:dd1)
{
  for (cj in 1:dd1)
  {
    if (ci != cj)
    {
      rho_w_eps1<-rho_w_eps1+Sdeermei[ci,]*Conj(Sdeermei[cj,])
    }
  }
}
rho_w_eps1<-rho_w_eps1/(dd1*(dd1-1)*Smeimei)

#compute the portion explained by surrogates
rho_w_delj1<-matrix(complex(0,0),nsurrog,length(freq))
for (scounter in 1:nsurrog)
{
  Smeimei<-sres_s[[scounter]][dd1+1,dd1+1,]
  Sdeermei<-sres_s[[scounter]][1:dd1,dd1+1,]
  for (ci in 1:dd1)
  {
    for (cj in 1:dd1)
    {
      if (ci != cj)
      {
        rho_w_delj1[scounter,]<-rho_w_delj1[scounter,]+Sdeermei[ci,]*Conj(Sdeermei[cj,])
      }
    }
  }
  rho_w_delj1[scounter,]<-rho_w_delj1[scounter,]/(dd1*(dd1-1)*Smeimei)
}

#check for non-zero imaginary parts and if you don't find them take the real part
if (any(Im(rho_w)>1e-10))
{
  stop("Error in do_analysis_some_locs: rho_w was not real")
}
if (any(Im(rho_w_eps1)>1e-10))
{
  stop("Error in do_analysis_some_locs: rho_w_eps1 was not real")
}
if (any(Im(rho_w_delj1)>1e-10))
{
  stop("Error in do_analysis_some_locs: rho_w_delj1 was not real")
}
rho_w<-Re(rho_w)
rho_w_eps1<-Re(rho_w_eps1)
rho_w_delj1<-Re(rho_w_delj1)

#make a plot
png("Results/Fourier_Step2_deermei_results.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
rho_w_delj1_quant<-apply(FUN=quantile,MARGIN=2,X=rho_w_delj1,probs=c(.025,.975))
ylimits<-range(rho_w,rho_w_eps1,rho_w_delj1_quant)
plot(freq,rho_w,type="l",ylim=ylimits)
mtext("Frequency",side=1,line=2,cex=1)
mtext("Synchrony",side=2,line=2,cex=1)
lines(freq,rho_w_eps1,type="l",lty="dashed",col="red")
lines(freq,rho_w_delj1_quant[1,],lty="dotted",col="red")
lines(freq,rho_w_delj1_quant[2,],lty="dotted",col="red")
lines(rep(1/7,2),ylimits,type="l",lty="dotted")
lines(rep(1/3,2),ylimits,type="l",lty="dotted")
lines(rep(1/4,2),ylimits,type="l",lty="dotted")
dev.off()

#now get fractions explained in each band
fracexplFMT_mei_3to7<-sum(rho_w_eps1[freq>=1/7 & freq<=1/3])/sum(rho_w[freq>=1/7 & freq<=1/3])
fracexplFMT_mei_3to4<-sum(rho_w_eps1[freq>=1/4 & freq<=1/3])/sum(rho_w[freq>=1/4 & freq<=1/3])
fracexplFMT_mei_4to7<-sum(rho_w_eps1[freq>=1/7 & freq<=1/4])/sum(rho_w[freq>=1/7 & freq<=1/4])

#rank of ranks significance procedure
#now do the rank of ranks procedure
rho_w_eps1_fl<-NA*rho_w_eps1
rho_w_delj1_fl<-NA*rho_w_delj1
for (counter in 1:length(freq))
{
  print(paste0("Rank of ranks procedure, frequency ",counter," of ",length(freq)))
  rho_w_eps1_fl[counter]<-sum(rho_w_delj1[,counter]<rho_w_eps1[counter])/nsurrog
  h<-rho_w_delj1[,counter]
  if (length(h)!=length(unique(h))) {stop("Error at location Y: ties where there should not be")}
  rh<-rank(h)
  rho_w_delj1_fl[,counter]<-(rh-1)/(nsurrog-1)
}
mnrk3to7<-mean(rho_w_eps1_fl[freq>=1/7 & freq<=1/3])
mnrk3to4<-mean(rho_w_eps1_fl[freq>=1/4 & freq<=1/3])
mnrk4to7<-mean(rho_w_eps1_fl[freq>=1/7 & freq<=1/4])
mnrk3to7_s<-apply(FUN=mean,X=rho_w_delj1_fl[,freq>=1/7 & freq<=1/3],MARGIN=1)
mnrk3to4_s<-apply(FUN=mean,X=rho_w_delj1_fl[,freq>=1/4 & freq<=1/3],MARGIN=1)
mnrk4to7_s<-apply(FUN=mean,X=rho_w_delj1_fl[,freq>=1/7 & freq<=1/4],MARGIN=1)
pvalresFMT_mei_3to7<-1-sum(mnrk3to7>=mnrk3to7_s)/nsurrog
pvalresFMT_mei_3to4<-1-sum(mnrk3to4>=mnrk3to4_s)/nsurrog
pvalresFMT_mei_4to7<-1-sum(mnrk4to7>=mnrk4to7_s)/nsurrog

#An alternative significance testing procedure - tests significance of the actual fractions of synchrony explained in
#each band that we may end up focussing on 
rho_w_eps1_3to7<-sum(rho_w_eps1[freq>=1/7 & freq<=1/3])
rho_w_delj1_3to7<-apply(FUN=sum,X=rho_w_delj1[,freq>=1/7 & freq<=1/3],MARGIN=1)
pvalresFMTalt_mei_3to7<-1-sum(rho_w_eps1_3to7>=rho_w_delj1_3to7)/nsurrog

rho_w_eps1_3to4<-sum(rho_w_eps1[freq>=1/4 & freq<=1/3])
rho_w_delj1_3to4<-apply(FUN=sum,X=rho_w_delj1[,freq>=1/4 & freq<=1/3],MARGIN=1)
pvalresFMTalt_mei_3to4<-1-sum(rho_w_eps1_3to4>=rho_w_delj1_3to4)/nsurrog

rho_w_eps1_4to7<-sum(rho_w_eps1[freq>=1/7 & freq<=1/4])
rho_w_delj1_4to7<-apply(FUN=sum,X=rho_w_delj1[,freq>=1/7 & freq<=1/4],MARGIN=1)
pvalresFMTalt_mei_4to7<-1-sum(rho_w_eps1_4to7>=rho_w_delj1_4to7)/nsurrog

rm(sres_s) #to save memory for later processes
