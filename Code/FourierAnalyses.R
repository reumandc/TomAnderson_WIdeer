#***
#Some spectral tools
#***

#A raw, unsmoothed periodigram. Detrending is done by this function. 
#
#Args
#x      A time series as a vector
#
myspecraw<-function(x)
{
  tforx<-1:length(x)
  x<-residuals(lm(x~tforx))
  h<-fft(x)
  h<-Re(h*Conj(h))
  freq<-(0:length(x))/length(x)
  h<-h[-1]
  freq<-freq[-1]
  h<-h[freq<=0.5]
  freq<-freq[freq<=0.5]
  return(list(freq=freq,spec=h))
}

#The power spectrum, Brillinger's consistent estimator (5.6 of Brillinger's 2001 
#book). The only difference from what is described there is frequency is here in 
#units of cycles per sampling interval in the output here, and was in radians per 
#sampling interval in Brillinger. Detrending is done by this function. The function 
#actually returns the log scale power spectrum.
#
#Args
#x      A time series as a vector
#
myspecbrill<-function(x,detrend=TRUE)
{
  Tx<-length(x)
  
  #detrend, if desired
  if (detrend==TRUE)
  {
    tforx<-1:Tx
    x<-residuals(lm(x~tforx))
  }
  
  #get the raw periodogram
  fftx<-fft(x)
  I<-(Mod(fftx))^2/(2*pi*Tx)
  I[1]<-0 #Set zero frequency to 0. Should be zero anyway, yto within rounding error, because of the detrending
  freq<-(0:(Tx-1))/Tx
  freq<-2*pi*freq #to make frequencies be in units of radians per sampling interval
  
  #now do the smoothing that makes the estimator
  BTx<-0.5/sqrt(Tx) #adjust BT to adjust the bias-variance tradeoff
  TxBTx<-Tx*BTx
  xforW<-(0:floor(TxBTx))/TxBTx
  WT<-(15/(16*2*pi))*((xforW-1)^2)*((xforW+1)^2)
  intWsquare<-(15/(32*pi))^2*4*pi*(1/9-4/7+6/5-4/3+1)
  spec<-WT[1]*I
  lenI<-length(I)
  for (AbsShift in 1:(length(WT)-1))
  {
    temp<-WT[AbsShift+1]*I
    #spec<-spec + temp([(AbsShift+1):end 1:AbsShift],:) + temp([(end-AbsShift+1):end 1:(end-AbsShift)],:);
    spec<-spec+temp[c((AbsShift+1):lenI,1:AbsShift)]+temp[c((lenI-AbsShift+1):lenI,1:(lenI-AbsShift))]
  }
  
  #remove the 0 frequency, and change the units back to radians per sampling interval,
  #and cut the redundant part of the spectrum
  freq<-freq[-1]
  spec<-spec[-1]
  freq<-freq/(2*pi)
  spec<-spec[freq<=0.5]
  freq<-freq[freq<=0.5]
  
  #put in some normalization factors
  spec<-spec*2*pi/TxBTx
  
  #now get confidence intervals
  p<-0.95 #for 95% confidence intervals
  conf<-qnorm((1+p)/2,mean=0,sd=1)*0.4343*sqrt(2*pi*intWsquare/(TxBTx));
  
  return(list(freq=freq,log10spec=log10(spec),conf=conf))
}

#***
#Comparison of Fourier transforms of deer and the surrogates developed for fig. 5
#***

#pull in the deer data
d<-readRDS(file="Results/cty.list.rds")
deer<-d$Abun
totdeer<-apply(FUN=sum,X=deer,MARGIN=2)
deeryr<-1981:2016

#pull in the Fourier surrogates developed for Fig 5
abunsurr<-read.csv("Data/abunsurrsum.csv")
abunsurr<-as.matrix(abunsurr)

#get raw periodograms of data and surrogates
deerspecraw<-myspecraw(totdeer)
freqs<-deerspecraw$freq
deerspecraw<-deerspecraw$spec
surrspecsraw<-matrix(NA,dim(abunsurr)[1],length(freqs))
for (counter in 1:dim(abunsurr)[1])
{
  x<-unname(abunsurr[counter,])
  h<-myspecraw(x)
  surrspecsraw[counter,]<-h$spec
}

#sum over 3-7 year timescales, compare data and surrogate sums using a histogram
realsum<-sum(deerspecraw[freqs>=1/7 & freqs<=1/3])
surrsums<-apply(FUN=sum,X=surrspecsraw[,freqs>=1/7 & freqs<=1/3],MARGIN=1)
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier1.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
hist(surrsums,50,xlim=range(realsum,surrsums),main="")
mtext("Power in 3-7 yr band",side=1,line=2,cex=1)
mtext("Count",side=2,line=2,cex=1)
points(realsum,0,col="red",pch=20)
#sum(surrsums>realsum) 
dev.off()

#Now detrend in prep for the next step. This is done using the same detrending 
#line for the real and all the surrogate time series
y<-c(totdeer,t(abunsurr))
x<-rep(deeryr,times=1+dim(abunsurr)[1])
dtmod<-lm(y~x)
totdeerres<-totdeer-coef(dtmod)[2]*deeryr-coef(dtmod)[1]
abunsurrres<-abunsurr-coef(dtmod)[2]*matrix(rep(deeryr,each=dim(abunsurr)[1]),dim(abunsurr)[1],length(deeryr))-coef(dtmod)[1]

#Now get Brillinger consistent estimators of the spectra of data and surrogates
deerlog10specbrill<-myspecbrill(totdeerres,detrend=FALSE)
freqs<-deerlog10specbrill$freq
deerlog10specbrill<-deerlog10specbrill$log10spec
surrlog10specsbrill<-matrix(NA,dim(abunsurr)[1],length(freqs))
for (counter in 1:dim(abunsurr)[1])
{
  x<-unname(abunsurrres[counter,])
  h<-myspecbrill(x,detrend=FALSE)
  surrlog10specsbrill[counter,]<-h$log10spec
}

#Now make a plot
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier2.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(10^deerlog10specbrill,10^surrlog10specsbrill)
plot(freqs,10^deerlog10specbrill,type="l",lwd=2,ylim=ylimits)
for (counter in 1:(dim(surrlog10specsbrill)[1]))
{
  lines(freqs,10^surrlog10specsbrill[counter,],type="l",col="grey",lwd=.5)
}
lines(rep(1/3,2),ylimits,lty="dashed")
lines(rep(1/7,2),ylimits,lty="dashed")
mtext("Frequency",side=1,line=2,cex=1)
mtext("Power spectrum",side=2,line=2,cex=1)
dev.off()

#***
#Now compare to AR(1) surrogates
#***

#Fit an AR(1) model to the detrended data
totdeerres<-residuals(lm(totdeer~deeryr))
h<-ar.mle(totdeerres,order.max=1)
#so it gives a model with AR(1) coefficient 0.2516 and variance of the 
#inovations 9.098E9 
rho<-h$ar
sig<-sqrt(h$var.pred)

#get surrogates based on those
getarsurr<-function(rho,sig,len,burn,numrun)
{
  res<-matrix(NA,numrun,len)
  res[,1]<-0
  for(counter in 2:len)
  {
    res[,counter]<-res[,counter-1]*rho+rnorm(numrun,mean=0,sd=sig)
  }
  return(res[,(burn+1):len])
}
numrun<-10000
len<-1000
burn<-len-length(totdeerres)
abunsurr<-getarsurr(rho,sig,len,burn,numrun)

#get raw periodograms of data and surrogates
deerspecraw<-myspecraw(totdeer) #recall this detrends automatically
freqs<-deerspecraw$freq
deerspecraw<-deerspecraw$spec
surrspecsraw<-matrix(NA,dim(abunsurr)[1],length(freqs))
for (counter in 1:dim(abunsurr)[1])
{
  x<-unname(abunsurr[counter,])
  h<-myspecraw(x)
  surrspecsraw[counter,]<-h$spec
}

#sum over 3-7 year timescales, compare data and surrogate sums using a histogram
realsum<-sum(deerspecraw[freqs>=1/7 & freqs<=1/3])
surrsums<-apply(FUN=sum,X=surrspecsraw[,freqs>=1/7 & freqs<=1/3],MARGIN=1)
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier3.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
hist(surrsums,50,xlim=range(realsum,surrsums),main="")
mtext("Power in 3-7 yr band",side=1,line=2,cex=1)
mtext("Count",side=2,line=2,cex=1)
points(realsum,0,col="red",pch=20)
dev.off()
Fourier3PvalRes<-sum(surrsums>realsum)/numrun #The p-value of the test of the null
#hypothesis that power in the 3-7yr band is no greater for the real data than it is
#for surrogates.
saveRDS(Fourier3PvalRes,file="Results/Fourier3PvalResult.rds")
Fourier3NumEquivTests<-length(freqs)/sum(freqs>=1/7 & freqs<=1/3)
saveRDS(Fourier3NumEquivTests,file="Results/Fourier3NumberEquivalentTests.rds")
Fourier3BonferThresh<-0.05/ceiling(Fourier3NumEquivTests)
saveRDS(Fourier3BonferThresh,file="Results/Fourier3BonferroniThreshold.rds")

#Now get Brillinger consistent estimators of the spectra of data and surrogates
deerlog10specbrill<-myspecbrill(totdeer)
freqs<-deerlog10specbrill$freq
deerlog10specbrill<-deerlog10specbrill$log10spec
surrlog10specsbrill<-matrix(NA,dim(abunsurr)[1],length(freqs))
for (counter in 1:dim(abunsurr)[1])
{
  x<-unname(abunsurr[counter,])
  h<-myspecbrill(x)
  surrlog10specsbrill[counter,]<-h$log10spec
}

#Now make a plot
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier4.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(10^deerlog10specbrill,10^surrlog10specsbrill)
plot(freqs,10^deerlog10specbrill,type="l",lwd=2,ylim=ylimits,xlab="",ylab="")
for (counter in 1:(dim(surrlog10specsbrill)[1]))
{
  lines(freqs,10^surrlog10specsbrill[counter,],type="l",col="grey",lwd=.5)
}
lines(freqs,10^deerlog10specbrill,type="l",lwd=2,ylim=ylimits)
lines(rep(1/3,2),ylimits,lty="dashed")
lines(rep(1/7,2),ylimits,lty="dashed")
qsurrspec<-apply(FUN=quantile,X=surrlog10specsbrill,MARGIN=2,prob=c(0.95,1-.05/sum(freqs>=1/7 & freqs<=1/3)))
lines(freqs,10^qsurrspec[1,],type="l",lty="dashed")
lines(freqs,10^qsurrspec[2,],type="l",lty="dotted") #Bonferroni corrected threshold
mtext("Frequency",side=1,line=2,cex=1)
mtext("Power spectrum",side=2,line=2,cex=1)
dev.off()

#Now make a plot that includes the stricted multiple-testing correction,
#which is Bonferroni counting all frequencies as independent tests. 
#Probably this is just for use in the referee responses.
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier5.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(10^deerlog10specbrill,10^surrlog10specsbrill)
plot(freqs,10^deerlog10specbrill,type="l",lwd=2,ylim=ylimits,xlab="",ylab="")
for (counter in 1:(dim(surrlog10specsbrill)[1]))
{
  lines(freqs,10^surrlog10specsbrill[counter,],type="l",col="grey",lwd=.5)
}
lines(freqs,10^deerlog10specbrill,type="l",lwd=2,ylim=ylimits)
lines(rep(1/3,2),ylimits,lty="dashed")
lines(rep(1/7,2),ylimits,lty="dashed")
qsurrspec<-apply(FUN=quantile,X=surrlog10specsbrill,MARGIN=2,prob=c(0.95,1-.05/sum(freqs>=1/7 & freqs<=1/3),1-0.05/length(freqs)))
lines(freqs,10^qsurrspec[1,],type="l",lty="dashed")
lines(freqs,10^qsurrspec[2,],type="l",lty="dotted") #Bonferroni corrected threshold, freqs in 3-7yr timescale band
lines(freqs,10^qsurrspec[3,],type="l",lty="dotted") #Bonferroni corrected threshold, all freqs
mtext("Frequency",side=1,line=2,cex=1)
mtext("Power spectrum",side=2,line=2,cex=1)
dev.off()

#Now do the AR argument with AIC
h1<-ar.mle(totdeerres)
saveRDS(h1$aic,file="Results/AR_AICs.rds") #save the AICs of the best-fitted AR model of each order
#make the spectrum of the best one, which is order 2, and the second-best (because it has a similar 
#AIC), which is order 3
AR_AIC_Table<-data.frame(ARorder=0:12,AIC=h1$aic)
h2o2<-spec.ar(totdeerres,method="mle",aic=FALSE,order.max=2,plot=FALSE) 
h2o3<-spec.ar(totdeerres,method="mle",aic=FALSE,order.max=3,plot=FALSE) 
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier6.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
plot(h2o2$freq,h2o2$spec,type="l",lty="solid")
lines(h2o3$freq,h2o3$spec,type="l",lty="dashed",col="red")
mtext("Frequency",side=1,line=2,cex=1)
mtext("Power spectrum",side=2,line=2,cex=1)
dev.off()

#Now finally just show the time series
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier7.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
plot(deeryr,totdeerres,type="l",xlab="",ylab="")
mtext("Year",side=1,line=2,cex=1)
mtext("Detrended total WI deer pop.",side=2,line=2,cex=1)
dev.off()
