#***libraries and other pre work

rm(list=ls())
library(mvtnorm)
source("my.spec.pgram.R")

#***settings

N<-10 #number of sampling locations
lensim<-1024 #length of simulations
lensimshort<-128

#***the model

#Simulates the model
#
#Args
#N        The number of sampling locations
#rho1     Covariance between locations for the inovations that go into producing 
#           the first environmental noise
#rho2     Covariance between locations used directly for the second environmental
#           noise
#c1
#c2       These two args specify the AR(2) dynamics for first environmental noise
#rhopop   Lag-1 autocorrelation coefficient for population dynamics
#lensim   Length of simulations to do
#numsims  Number of simulations to do
#
#Output - a list with these elements 
#inovations - a lensim by N by numsims array, used for the first environmental noise
#epsilon1 - a lensim by N by numsims array, the first environmental noise
#epsilon2 - a lensim by N by numsims array, the second environmental noise
#pops - a lensim by N by numsims array
simmodel2<-function(N,rho1,rho2,c1,c2,rhopop,lensim,numsims)
{
  #**generate the inovations for the generation of the first environmental noise
  
  Sig<-matrix(rho1,N,N)
  diag(Sig)<-1
  inovations<-rmvnorm(3*lensim*numsims,sigma=Sig)
  dim(inovations)<-c(3*lensim,numsims,N)
  inovations<-aperm(inovations,c(1,3,2))
  
  #**now use these to generate the first environmental noise
  
  #now make the noise
  epsilon1<-array(0,dim=c(3*lensim,N,numsims))
  for (counter in 3:(3*lensim))
  {
    epsilon1[counter,,]<-c1*epsilon1[counter-1,,]+
      c2*epsilon1[counter-2,,]+inovations[counter,,]
  }
  
  #throw away transients
  epsilon1<-epsilon1[(dim(epsilon1)[1]-2*lensim+1):(dim(epsilon1)[1]),,]
  inovations<-inovations[(dim(inovations)[1]-2*lensim+1):(dim(inovations)[1]),,]
  
  #**now get the second environmental noise
  
  Sig<-matrix(rho2,N,N)
  diag(Sig)<-1
  epsilon2<-rmvnorm(2*lensim*numsims,sigma=Sig)
  dim(epsilon2)<-c(2*lensim,numsims,N)
  epsilon2<-aperm(epsilon2,c(1,3,2))
  
  #**now add the two environmental noises and use them to simulate pops
  
  #add the noises
  epsilon<-epsilon1+epsilon2
  
  #do the sims  
  pops<-array(0,dim=c(2*lensim,N,numsims))
  for (counter in 2:(2*lensim))
  {
    pops[counter,,]<-rhopop*pops[counter-1,,]+epsilon[counter,,]
  }
  
  #throw away transients  
  pops<-pops[(dim(pops)[1]-lensim+1):(dim(pops)[1]),,]
  epsilon1<-epsilon1[(dim(epsilon1)[1]-lensim+1):(dim(epsilon1)[1]),,]
  epsilon2<-epsilon2[(dim(epsilon2)[1]-lensim+1):(dim(epsilon2)[1]),,]
  inovations<-inovations[(dim(inovations)[1]-lensim+1):(dim(inovations)[1]),,]

  return(list(inovations=inovations,
              epsilon1=epsilon1,
              epsilon2=epsilon2,
              pops=pops))
}

#***do the sims 

#do the sims, case 1
rho1_cs1<-0
rho2_cs1<-0.8
r<-1.5
theta<-pi/2
c1<-(1/r)*2*cos(theta)
c2<-(-1/r^2) 
rhopop<-0.15
numsims<-100
set.seed(101)
sims_case1<-simmodel2(N=N,rho1=rho1_cs1,rho2=rho2_cs1,c1=c1,c2=c2,rhopop=rhopop,
                      lensim=lensim,numsims=numsims)
  
#do the sims, case 2 - have to figure out rho1 for this case first
intgrd<-function(f,c1,c2)
{
  mu<-exp(-2*pi*complex(real=0,imaginary=1)*f)
  return((Mod(1/(1-c1*mu-c2*mu^2)))^2)
}
intres<-integrate(intgrd,-.5,.5,c1=c1,c2=c2)
if (intres$abs.error>1e-4 || intres$message!="OK")
{
  stop("Error location 1: problem with numeric integration")
}
rho1_cs2<-rho2_cs1/intres$value 
rho2_cs2<-0
sims_case2<-simmodel2(N=N,rho1=rho1_cs2,rho2=rho2_cs2,c1=c1,c2=c2,rhopop=rhopop,
                        lensim=lensim,numsims=numsims)

#***check that covariances between locations are the same in expectation 
#in the two cases

allcovs_case1<-NA*numeric(numsims)
allcovs_case2<-NA*numeric(numsims)
epsilon_cs1<-sims_case1$epsilon1+sims_case1$epsilon2
epsilon_cs2<-sims_case2$epsilon1+sims_case2$epsilon2
for (counter in 1:numsims)
{
  allcovs_case1[counter]<-cov(epsilon_cs1[,1,counter],epsilon_cs1[,2,counter])
  allcovs_case2[counter]<-cov(epsilon_cs2[,1,counter],epsilon_cs2[,2,counter])
}
pdf(file="Example2_case1_HistCovs.pdf")
hist(allcovs_case1)
dev.off()
pdf(file="Example2_case2_HistCovs.pdf")
hist(allcovs_case2)
dev.off()
mean(allcovs_case1)
mean(allcovs_case2)
#these are pretty similar so it seems to be working

#***show the spectra of noise in the two cases

#get ready
epsilon_cs1<-sims_case1$epsilon1+sims_case1$epsilon2
specnoise_cs1<-my.spec.pgram(epsilon_cs1[,1,],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-specnoise_cs1$freq
y_cs1<-specnoise_cs1$spec
epsilon_cs2<-sims_case2$epsilon1+sims_case2$epsilon2
specnoise_cs2<-my.spec.pgram(epsilon_cs2[,1,],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
y_cs2<-specnoise_cs2$spec

#now make the plot for case 1, numerics
pdf(file="Examp2_NoiseSpectra_cs1.pdf")
plot(x,y_cs1[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(y_cs1,y_cs2),
     xlab="Frequency",ylab="Power",main="Spectra, noise, loc 1, cs 1")
for (counter in 1:numsims)
{
  lines(x,y_cs1[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=y_cs1,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#add the theory
mu<-exp(-2*pi*complex(real=0,imaginary=1)*x)
y_cs1_an<-(Mod(1/(1-c1*mu-c2*mu^2)))^2+1
lines(x,y_cs1_an,col="red")
dev.off()

#make the plot for case 2, same axis extent
pdf(file="Examp2_NoiseSpectra_cs2.pdf")
plot(x,y_cs2[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(y_cs1,y_cs2),
     xlab="Frequency",ylab="Power",main="Spectra, noise, loc 1, cs 2")
for (counter in 1:numsims)
{
  lines(x,y_cs2[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=y_cs2,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#add the theory
y_cs2_an<-(Mod(1/(1-c1*mu-c2*mu^2)))^2+1
lines(x,y_cs2_an,col="red")
dev.off()

#**show the cospectra of noise in the two cases

#get ready, also getting ready for the quad spectra at the same time, which come next
crsp_cs1<-my.spec.pgram(epsilon_cs1[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
allcosp_cs1<-matrix(NA,dim(crsp_cs1$pgram)[1],numsims)
allquad_cs1<-matrix(NA,dim(crsp_cs1$pgram)[1],numsims)
allcosp_cs1[,1]<-Re(crsp_cs1$pgram[,1,2])
allquad_cs1[,1]<-Im(crsp_cs1$pgram[,1,2])
crsp_cs2<-my.spec.pgram(epsilon_cs2[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
allcosp_cs2<-matrix(NA,dim(crsp_cs2$pgram)[1],numsims)
allquad_cs2<-matrix(NA,dim(crsp_cs2$pgram)[1],numsims)
allcosp_cs2[,1]<-Re(crsp_cs2$pgram[,1,2])
allquad_cs2[,1]<-Im(crsp_cs2$pgram[,1,2])
for (counter in 2:numsims)
{
  crsp_cs1<-my.spec.pgram(epsilon_cs1[,1:2,counter],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  crsp_cs2<-my.spec.pgram(epsilon_cs2[,1:2,counter],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  allcosp_cs1[,counter]<-Re(crsp_cs1$pgram[,1,2])
  allquad_cs1[,counter]<-Im(crsp_cs1$pgram[,1,2])
  allcosp_cs2[,counter]<-Re(crsp_cs2$pgram[,1,2])
  allquad_cs2[,counter]<-Im(crsp_cs2$pgram[,1,2])
}

#show case 1, numerics
pdf(file="Examp2_NoiseCospectra_cs1.pdf")
plot(x,allcosp_cs1[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(allcosp_cs1,allcosp_cs2),
     xlab="Frequency",ylab="Cospectrum",main="Cospectra, noise, locs 1 and 2, cs 1")
for (counter in 2:numsims)
{
  lines(x,allcosp_cs1[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=allcosp_cs1,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#add the theory
lines(x,rep(rho2_cs1,length(x)),col="red")
dev.off()

#shows case 2, same y axis extent
pdf(file="Examp2_NoiseCospectra_cs2.pdf")
plot(x,allcosp_cs2[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(allcosp_cs1,allcosp_cs2),
     xlab="Frequency",ylab="Cospectrum",main="Cospectra, noise, locs 1 and 2, cs 2")
for (counter in 2:numsims)
{
  lines(x,allcosp_cs2[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=allcosp_cs2,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#add the theory
y_cs2_an<-rho1_cs2*(Mod(1/(1-c1*mu-c2*mu^2)))^2+rho2_cs2
lines(x,y_cs2_an,col="red")
dev.off()

#**show the quad spectra of noise in the two cases

#show case 1, numerics
pdf(file="Examp2_NoiseQuadspectra_cs1.pdf")
plot(x,allquad_cs1[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(allquad_cs1,allquad_cs2),
     xlab="Frequency",ylab="Quad. spect.",main="Quad. spect., noise, locs 1 and 2, cs 1")
for (counter in 2:numsims)
{
  lines(x,allquad_cs1[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=allquad_cs1,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#add the theory
lines(x,rep(0,length(x)),col="red")
dev.off()

#shows case 2, same y axis extent
pdf(file="Examp2_NoiseQuadspectra_cs2.pdf")
plot(x,allquad_cs2[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(allquad_cs1,allquad_cs2),
     xlab="Frequency",ylab="Quad. spect.",main="Quad. spect., noise, locs 1 and 2, cs 2")
for (counter in 2:numsims)
{
  lines(x,allquad_cs2[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=allquad_cs2,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#add the theory
lines(x,rep(0,length(x)),col="red")
dev.off()

#**show the spectra of pops in the two cases


#**show the cospectra of pops in the two cases


#**show the quad spectra of pops in the two cases


#**show the spectra of total populations in the two cases

#get ready
totpops_cs1<-apply(FUN=sum,X=sims_case1$pops,MARGIN=c(1,3))
totpops_cs2<-apply(FUN=sum,X=sims_case2$pops,MARGIN=c(1,3))
specs_cs1<-my.spec.pgram(totpops_cs1,detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
specs_cs2<-my.spec.pgram(totpops_cs2,detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-specs_cs1$freq
y_cs1<-specs_cs1$spec
y_cs2<-specs_cs2$spec

#make the plot for case 1, numerics
pdf(file="Examp2_TotpopsSpectra_cs1.pdf")
plot(x,y_cs1[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(y_cs1,y_cs2),
     xlab="Frequency",ylab="Power",main="Spectra, total pop, cs 1")
for (counter in 1:numsims)
{
  lines(x,y_cs1[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=y_cs1,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#add the theory
mu<-exp(-2*pi*complex(real=0,imaginary=1)*x)
y_cs1_an<-(Mod(1/(1-rhopop*mu)))^2*((N^2-N)*rho2_cs1+N+N*(Mod(1/(1-c1*mu-c2*mu^2)))^2)
lines(x,y_cs1_an,col="red")
dev.off()

#make the plot for case 2, y axis extent the same
pdf(file="Examp2_TotpopsSpectra_cs2.pdf")
plot(x,y_cs2[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(y_cs1,y_cs2),
     xlab="Frequency",ylab="Power",main="Spectra, total pop, cs 2")
for (counter in 1:numsims)
{
  lines(x,y_cs2[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=y_cs2,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#add the theory
y_cs2_an<-(Mod(1/(1-rhopop*mu)))^2*(((N^2-N)*rho1_cs2+N)*(Mod(1/(1-c1*mu-c2*mu^2)))^2+N)
lines(x,y_cs2_an,col="red")
dev.off()

#plot the total pop versus time for the two cases
pdf(file="Examp2_Totpops.pdf")
x<-101:151
xp<-x-min(x)+1
ycs1<-totpops_cs1[x,1]
ycs2<-totpops_cs2[x,1]
plot(xp,ycs1,type='l',ylim=range(ycs1,ycs2))
points(xp,ycs1,type='p',pch=20,cex=.5)
lines(xp,ycs2,type='l',col='green')
points(xp,ycs2,type='p',pch=20,cex=.5,col='green')
dev.off()

#do something similar, but subtract a long-term moving average
#to get rid of the red shift
pdf(file="Examp2_TotpopsMinusMA.pdf")
avgsz<-3
ycs1<-NA*numeric(length(x))
ycs2<-NA*numeric(length(x))
for (counter in 1:length(x))
{
  ycs1[counter]<-totpops_cs1[x[counter],1]-mean(totpops_cs1[(x[counter]-avgsz):(x[counter]+avgsz),1])
  ycs2[counter]<-totpops_cs2[x[counter],1]-mean(totpops_cs2[(x[counter]-avgsz):(x[counter]+avgsz),1])
}
plot(xp,ycs1,type='l',ylim=range(ycs1,ycs2))
points(xp,ycs1,type='p',pch=20,cex=.5)
lines(xp,ycs2,type='l',col='green')
points(xp,ycs2,type='p',pch=20,cex=.5,col='green')
dev.off()

#determine the indices of peaks in the ts and look at peak spacing
ycs1<-NA*numeric(lensim-2*avgsz)
ycs2<-NA*numeric(lensim-2*avgsz)
for (counter in 1:length(ycs1))
{
  ycs1[counter]<-totpops_cs1[counter+avgsz,1]-mean(totpops_cs1[counter:(counter+2*avgsz),1])
  ycs2[counter]<-totpops_cs2[counter+avgsz,1]-mean(totpops_cs2[counter:(counter+2*avgsz),1])
}
pycs1<-c()
pycs2<-c()
for (counter in 2:(length(ycs1)-1))
{
  if (ycs1[counter]>=ycs1[counter-1] && ycs1[counter+1]<=ycs1[counter])
  {
    pycs1<-c(pycs1,counter)
  }
  if (ycs2[counter]>=ycs2[counter-1] && ycs2[counter+1]<=ycs2[counter])
  {
    pycs2<-c(pycs2,counter)
  }
}
dpycs1<-diff(pycs1)
dpycs2<-diff(pycs2)
hist(dpycs1)
hist(dpycs2)

cc1<-c()
cc2<-c()
for (counter in 1:10)
{
  cc1[counter]<-sum(dpycs1==counter)
  cc2[counter]<-sum(dpycs2==counter)
}
plot(1:10,cc1,type='l')
lines(1:10,cc2,type='l',col='green')
