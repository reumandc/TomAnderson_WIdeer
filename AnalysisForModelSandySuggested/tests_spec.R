#***libraries and other pre work

rm(list=ls())
library(mvtnorm)
source("my.spec.pgram.R")

#***see if the raw spectrum, summed, equals the variance

#get some data
set.seed(101)
lents<-1024
x<-rnorm(lents)
#x<-x-mean(x)

#take spectra in two ways to compare - should be the same
pg<-my.spec.pgram(x,detrend=FALSE,plot=FALSE,taper=0)
length(pg$spec)
y<-((Mod(fft(x)))^2)/lents
y<-y[2:length(y)]
ycut<-y[1:(length(pg$spec))]
cbind(ycut,pg$spec)
max(abs(ycut-pg$spec))
#so, without smoothing or tapering, spec.pgram gives 
#((Mod(fft))^2)/lents (with the 0 frequency removed
#and the redundant half of the spectrum removed)

#compare to the variance
var(x)*(lents-1)/lents #use the biased estimator of variance
sum(y)/lents
sum(y)/lents-(var(x)*(lents-1)/lents)

#now use the unbiased estimator
var(x)
sum(y)/(lents-1)
var(x)-sum(y)/(lents-1)      

#***now look at covariance and the cospectrum

#get some data
set.seed(101)
lents<-1024
rho<-.8
Sig<-matrix(c(1,.8,.8,1),2,2)
x<-rmvnorm(lents,mean=c(0,0),sigma=Sig)

#take cross spectra in two ways to compare - should be the same
pg<-my.spec.pgram(x,detrend=FALSE,plot=FALSE,taper=0)
dim(pg$pgram)
y<-(fft(x[,1])*Conj(fft(x[,2])))/lents
y<-y[2:length(y)]
ycut<-y[1:(dim(pg$pgram)[1])]
cbind(ycut,pg$pgram[,1,2])
max(Mod(ycut-pg$pgram[,1,2]))
#so, without smoothing or tapering, spec.pgram gives 
#(fft(x[,1])*Conj(fft(x[,2])))/lents (with the 0 frequency removed
#and the redundant half of the cross spectrum removed)

#compare to the covariance, go straight to the unbiased estimator
cov(x[,1],x[,2])
sum(y)/(lents-1)
cov(x[,1],x[,2])-Re(sum(y)/(lents-1))

#now look at the cospectrum across a bunch of runs
numts<-50
x<-rmvnorm(lents*numts,mean=c(0,0),sigma=Sig)
dim(x)<-c(lents,numts,2)
x<-aperm(x,c(1,3,2))
cs1<-my.spec.pgram(x[,,1],detrend=FALSE,plot=FALSE,taper=0,spans=c(91,71))
allcosp<-matrix(NA,dim(cs1$pgram)[1],numts)
allcosp[,1]<-Re(cs1$pgram[,1,2])
for (counter in 2:numts)
{
  cs1<-my.spec.pgram(x[,,counter],detrend=FALSE,plot=FALSE,taper=0,spans=c(91,71))
  allcosp[,counter]<-Re(cs1$pgram[,1,2])
}
plot(cs1$freq,allcosp[,1],type='l',col=rgb(0,0,0,alpha=.3),ylim=c(-1,1))
for (counter in 2:numts)
{
  lines(cs1$freq,allcosp[,counter],type='l',col=rgb(0,0,0,alpha=.3))
}
lines(cs1$freq,rep(rho,length(cs1$freq)),col='red')
meds<-apply(FUN=median,X=allcosp,MARGIN=1)
lines(cs1$freq,meds,col='red',lty='dashed')
