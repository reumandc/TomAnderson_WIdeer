#work on a first example with AR(1) dynamics, noise with weak periodicity 
#and possibly some synchrony (which is therefore timescale-specific synchrony,
#if there is synchrony)

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
#rho      Covariance between locations for the inovations that go into producing 
#           the environmental noise
#r
#theta    These two args specify the AR(2) dynamics for noise
#rhopop   Lag-1 autocorrelation coefficient for dynamics
#lensim   Length of simulations to do
#numsims  Number of simulations to do
#
#Output - a list with these elements, each 
#inovations - a lensim by N by numsims array
#epsilon - a lensim by N by numsims array
#pops - a lensim by N by numsims array
simmodel1<-function(N,rho,r,theta,rhopop,lensim,numsims)
{
  #**generate the inovations for the generation of the noise
  
  Sig<-matrix(rho,N,N)
  diag(Sig)<-1
  inovations<-rmvnorm(3*lensim*numsims,sigma=Sig)
  dim(inovations)<-c(3*lensim,numsims,N)
  inovations<-aperm(inovations,c(1,3,2))
  
  #**now use these to generate the noise
  
  #get c1 and c2
  c1<-(1/r)*2*cos(theta)
  c2<-(-1/r^2) 

  #now make the noise
  epsilon<-array(0,dim=c(3*lensim,N,numsims))
  for (counter in 3:(3*lensim))
  {
    epsilon[counter,,]<-c1*epsilon[counter-1,,]+
      c2*epsilon[counter-2,,]+inovations[counter,,]
  }
  
  #throw away transients
  epsilon<-epsilon[(dim(epsilon)[1]-2*lensim+1):(dim(epsilon)[1]),,]
  inovations<-inovations[(dim(inovations)[1]-2*lensim+1):(dim(inovations)[1]),,]

  #**now put the noise into an AR(1) to get pops

  #do the sims  
  pops<-array(0,dim=c(2*lensim,N,numsims))
  for (counter in 2:(2*lensim))
  {
    pops[counter,,]<-rhopop*pops[counter-1,,]+epsilon[counter,,]
  }

  #throw away transients  
  pops<-pops[(dim(pops)[1]-lensim+1):(dim(pops)[1]),,]
  epsilon<-epsilon[(dim(epsilon)[1]-lensim+1):(dim(epsilon)[1]),,]
  inovations<-inovations[(dim(inovations)[1]-lensim+1):(dim(inovations)[1]),,]


  return(list(inovations=inovations,
              epsilon=epsilon,
              pops=pops))
}

#***do the sims, show spectra, cospectra, etc.

#do the sims
rho<-0.8
r<-3
theta<-pi/2
rhopop<-0.15
numsims<-100
set.seed(101)
sims<-simmodel1(N=N,rho=rho,r=r,theta=theta,rhopop=rhopop,
                lensim=lensim,numsims=numsims)

#show the spectrum of the inovations, compare with theory
pdf(file="Examp1_InovationsSpectra.pdf")
spec_inovations_N1<-my.spec.pgram(sims$inovations[,1,],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-spec_inovations_N1$freq
plot(x,spec_inovations_N1$spec[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(spec_inovations_N1$spec),
     main="Inov spectra, loc 1") #plots for all sims
for (simc in 2:numsims)
{
  lines(x,spec_inovations_N1$spec[,simc],type="l",col=rgb(0,0,0,alpha=.1))
}
lines(x,rep(1,length(x)),col="red") #the analytic expectation
medn<-apply(FUN=median,X=spec_inovations_N1$spec,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the cospectrum of the innovations, compare with theory
pdf(file="Examp1_InovationsCospectra.pdf")
cs_inovations_N12<-my.spec.pgram(sims$inovations[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-cs_inovations_N12$freq
y<-Re(cs_inovations_N12$pgram[,1,2])
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(-1,1),main="Inov cospectrum, locs 1 and 2")
allvals<-matrix(NA,length(y),numsims)
allvals[,1]<-y
for (simc in 2:numsims)
{
  cs_inovations_N12<-my.spec.pgram(sims$inovations[,1:2,simc],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  y<-Re(cs_inovations_N12$pgram[,1,2])
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
  allvals[,simc]<-y
}
lines(x,rep(rho,length(x)),col="red") #the analytic expectation
medn<-apply(FUN=median,X=allvals,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the quad spectrum of the innovations, compare with theory
pdf(file="Examp1_InovationsQuadspectra.pdf")
cs_inovations_N12<-my.spec.pgram(sims$inovations[,1:2,1],spans=c(91,71),detrend=FALSE,taper=0,plot=FALSE)
x<-cs_inovations_N12$freq
y<-Im(cs_inovations_N12$pgram[,1,2])
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(-1,1),main="Inov quadspectrum, locs 1 and 2")
allvals<-matrix(NA,length(y),numsims)
allvals[,1]<-y
for (simc in 2:numsims)
{
  cs_inovations_N12<-my.spec.pgram(sims$inovations[,1:2,simc],spans=c(91,71),plot=FALSE)
  y<-Im(cs_inovations_N12$pgram[,1,2])
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
  allvals[,simc]<-y
}
lines(x,rep(0,length(x)),col="red") #the analytic expectation
medn<-apply(FUN=median,X=allvals,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the spectrum of the noise, compare with theory
pdf(file="Examp1_NoiseSpectra.pdf")
spec_noise_N1<-my.spec.pgram(sims$epsilon[,1,],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-spec_noise_N1$freq
plot(x,spec_noise_N1$spec[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(spec_noise_N1$spec),
     main="Noise spectra, loc 1") #plots for all sims
for (simc in 2:numsims)
{
  lines(x,spec_noise_N1$spec[,simc],type="l",col=rgb(0,0,0,alpha=.1))
}
mu<-exp(-2*pi*complex(real=0,imaginary=1)*x)
c1<-(1/r)*2*cos(theta)
c2<-(-1/r^2) 
y1<-(Mod(1/(1-c1*mu-c2*mu^2)))^2
lines(x,y1,col="red") #the analytic expectation
medn<-apply(FUN=median,X=spec_noise_N1$spec,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the cospectrum of the noise, compare with theory
pdf(file="Examp1_NoiseCospectra.pdf")
cs_noise_N12<-my.spec.pgram(sims$epsilon[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-cs_noise_N12$freq
y<-Re(cs_noise_N12$pgram[,1,2])
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(0,1.5),main="Noise cospectrum, locs 1 and 2")
allvals<-matrix(NA,length(y),numsims)
allvals[,1]<-y
for (simc in 2:numsims)
{
  cs_noise_N12<-my.spec.pgram(sims$epsilon[,1:2,simc],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  y<-Re(cs_noise_N12$pgram[,1,2])
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
  allvals[,simc]<-y
}
y1<-rho*(Mod(1/(1-c1*mu-c2*mu^2)))^2
lines(x,y1,col="red") #the analytic expectation
medn<-apply(FUN=median,X=allvals,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the quad spectrum of the noise, compare with theory
pdf(file="Examp1_NoiseQuadspectra.pdf")
cs_noise_N12<-my.spec.pgram(sims$epsilon[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-cs_noise_N12$freq
y<-Im(cs_noise_N12$pgram[,1,2])
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(-.5,.5),main="Noise quadspectrum, locs 1 and 2")
allvals<-matrix(NA,length(y),numsims)
allvals[,1]<-y
for (simc in 2:numsims)
{
  cs_noise_N12<-my.spec.pgram(sims$epsilon[,1:2,simc],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  y<-Im(cs_noise_N12$pgram[,1,2])
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
  allvals[,simc]<-y
}
lines(x,rep(0,length(x)),col="red") #the analytic expectation
medn<-apply(FUN=median,X=allvals,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the spectrum of the pops, compare with theory
pdf(file="Examp1_PopsSpectra.pdf")
spec_pops_N1<-my.spec.pgram(sims$pops[,1,],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-spec_pops_N1$freq
plot(x,spec_pops_N1$spec[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(spec_pops_N1$spec),
     main="Pops spectra, loc 1") #plots for all sims
for (simc in 2:numsims)
{
  lines(x,spec_pops_N1$spec[,simc],type="l",col=rgb(0,0,0,alpha=.1))
}
y1<-(Mod(1/(1-c1*mu-c2*mu^2)))^2
y2<-(Mod(1/(1-rhopop*mu)))^2
lines(x,y1*y2,col="red") #the analytic expectation
medn<-apply(FUN=median,X=spec_pops_N1$spec,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the cospectrum of the pops, compare with theory
pdf(file="Examp1_PopsCospectra.pdf")
cs_pops_N12<-my.spec.pgram(sims$pops[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-cs_pops_N12$freq
y<-Re(cs_pops_N12$pgram[,1,2])
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(0,1.5),main="Pops cospectrum, locs 1 and 2")
allvals<-matrix(NA,length(y),numsims)
allvals[,1]<-y
for (simc in 2:numsims)
{
  cs_pops_N12<-my.spec.pgram(sims$pops[,1:2,simc],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  y<-Re(cs_pops_N12$pgram[,1,2])
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
  allvals[,simc]<-y
}
y1<-rho*(Mod(1/(1-c1*mu-c2*mu^2)))^2
lines(x,y1*y2,col="red") #the analytic expectation
medn<-apply(FUN=median,X=allvals,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the quadspectrum of the pops, compare with theory
pdf(file="Examp1_PopsQuadspectra.pdf")
cs_pops_N12<-my.spec.pgram(sims$pops[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-cs_pops_N12$freq
y<-Im(cs_pops_N12$pgram[,1,2])
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(-.5,.5),main="Pops quadspectrum, locs 1 and 2")
allvals<-matrix(NA,length(y),numsims)
allvals[,1]<-y
for (simc in 2:numsims)
{
  cs_pops_N12<-my.spec.pgram(sims$pops[,1:2,simc],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  y<-Im(cs_pops_N12$pgram[,1,2])
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
  allvals[,simc]<-y
}
lines(x,rep(0,length(x)),col="red") #the analytic expectation
medn<-apply(FUN=median,X=allvals,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#now compute the total population for each sim, and look at spectra,
#and compare with theory
pdf(file="Examp12_TotpopsSpectra.pdf")
totpops<-apply(FUN=sum,X=sims$pops,MARGIN=c(1,3))
spec_totpops<-my.spec.pgram(totpops,detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-spec_totpops$freq
plot(x,spec_totpops$spec[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(0,max(spec_totpops$spec)),
     main="Totpops spectra") #plots for all sims
for (simc in 2:numsims)
{
  lines(x,spec_totpops$spec[,simc],type="l",col=rgb(0,0,0,alpha=.1))
}
y1<-(Mod(1/(1-c1*mu-c2*mu^2)))^2
y2<-(Mod(1/(1-rhopop*mu)))^2
prefact<-(N+(N^2-N)*rho)
lines(x,prefact*y1*y2,col="red") #the analytic expectation
medn<-apply(FUN=median,X=spec_totpops$spec,MARGIN=1)
lines(x,medn,col="red",lty="dashed")

#repeat the sims, but with rho=0
rho<-0
simsRho0<-simmodel1(N=N,rho=rho,r=r,theta=theta,rhopop=rhopop,
                lensim=lensim,numsims=numsims)

#now compute the total population for each sim fpor rho=0, and look at spectra,
#and compare with theory, and add all this to the previous plot
totpopsRho0<-apply(FUN=sum,X=simsRho0$pops,MARGIN=c(1,3))
spec_totpopsRho0<-my.spec.pgram(totpopsRho0,detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-spec_totpopsRho0$freq
for (simc in 1:numsims)
{
  lines(x,spec_totpopsRho0$spec[,simc],type="l",col=rgb(0,0,0,alpha=.1))
}
y1<-(Mod(1/(1-c1*mu-c2*mu^2)))^2
y2<-(Mod(1/(1-rhopop*mu)))^2
prefact<-(N+(N^2-N)*rho)
lines(x,prefact*y1*y2,col="red") #the analytic expectation
medn<-apply(FUN=median,X=spec_totpops$spec,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#plot the total pop versus time for rho=0.8 and rho=0
pdf(file="Examp1_Totpops.pdf")
x<-51:100
xp<-x-min(x)+1
yp8<-totpops[x]
y0<-totpopsRho0[x]
plot(xp,yp8,type='l',ylim=range(yp8,y0))
points(xp,yp8,type='p',pch=20,cex=.5)
lines(xp,y0,type='l',col='green')
points(xp,y0,type='p',pch=20,cex=.5,col='green')
dev.off()

#do something similar, but subtract a long-term moving average
#to get rid of the red shift
pdf(file="Examp1_TotpopsMinusMA.pdf")
avgsz<-3
yp8<-NA*numeric(length(x))
y0<-NA*numeric(length(x))
for (counter in 1:length(x))
{
  yp8[counter]<-totpops[counter]-mean(totpops[(x[counter]-avgsz):(x[counter]+avgsz)])
  y0[counter]<-totpopsRho0[counter]-mean(totpopsRho0[(x[counter]-avgsz):(x[counter]+avgsz)])
}
plot(xp,yp8,type='l',ylim=range(yp8,y0))
points(xp,yp8,type='p',pch=20,cex=.5)
lines(xp,y0,type='l',col='green')
points(xp,y0,type='p',pch=20,cex=.5,col='green')
dev.off()

#show the spectrum of the noise for rho=0, compare with theory
pdf(file="Examp1_NoiseSpectra_rho0.pdf")
spec_noise_N1<-my.spec.pgram(simsRho0$epsilon[,1,],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-spec_noise_N1$freq
plot(x,spec_noise_N1$spec[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(spec_noise_N1$spec),
     main="Noise spectra, loc 1") #plots for all sims
for (simc in 2:numsims)
{
  lines(x,spec_noise_N1$spec[,simc],type="l",col=rgb(0,0,0,alpha=.1))
}
mu<-exp(-2*pi*complex(real=0,imaginary=1)*x)
c1<-(1/r)*2*cos(theta)
c2<-(-1/r^2) 
y1<-(Mod(1/(1-c1*mu-c2*mu^2)))^2
lines(x,y1,col="red") #the analytic expectation
medn<-apply(FUN=median,X=spec_noise_N1$spec,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the cospectrum of the noise for rho=0, compare with theory
pdf(file="Examp1_NoiseCospectra_rho0.pdf")
cs_noise_N12<-my.spec.pgram(simsRho0$epsilon[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-cs_noise_N12$freq
y<-Re(cs_noise_N12$pgram[,1,2])
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(-.5,.5),main="Noise cospectrum, locs 1 and 2")
allvals<-matrix(NA,length(y),numsims)
allvals[,1]<-y
for (simc in 2:numsims)
{
  cs_noise_N12<-my.spec.pgram(simsRho0$epsilon[,1:2,simc],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  y<-Re(cs_noise_N12$pgram[,1,2])
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
  allvals[,simc]<-y
}
y1<-rho*(Mod(1/(1-c1*mu-c2*mu^2)))^2
lines(x,y1,col="red") #the analytic expectation
medn<-apply(FUN=median,X=allvals,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the quad spectrum of the noise for rho=0, compare with theory
pdf(file="Examp1_NoiseQuadspectra_rho0.pdf")
cs_noise_N12<-my.spec.pgram(simsRho0$epsilon[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-cs_noise_N12$freq
y<-Im(cs_noise_N12$pgram[,1,2])
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(-.5,.5),main="Noise quadspectrum, locs 1 and 2")
allvals<-matrix(NA,length(y),numsims)
allvals[,1]<-y
for (simc in 2:numsims)
{
  cs_noise_N12<-my.spec.pgram(simsRho0$epsilon[,1:2,simc],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  y<-Im(cs_noise_N12$pgram[,1,2])
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
  allvals[,simc]<-y
}
lines(x,rep(0,length(x)),col="red") #the analytic expectation
medn<-apply(FUN=median,X=allvals,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the spectrum of the pops for rho=0, compare with theory
pdf(file="Examp1_PopsSpectra_rho0.pdf")
spec_pops_N1<-my.spec.pgram(simsRho0$pops[,1,],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-spec_pops_N1$freq
plot(x,spec_pops_N1$spec[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(spec_pops_N1$spec),
     main="Pops spectra, loc 1") #plots for all sims
for (simc in 2:numsims)
{
  lines(x,spec_pops_N1$spec[,simc],type="l",col=rgb(0,0,0,alpha=.1))
}
y1<-(Mod(1/(1-c1*mu-c2*mu^2)))^2
y2<-(Mod(1/(1-rhopop*mu)))^2
lines(x,y1*y2,col="red") #the analytic expectation
medn<-apply(FUN=median,X=spec_pops_N1$spec,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the cospectrum of the pops for rho=0, compare with theory
pdf(file="Examp1_PopsCospectra_rho0.pdf")
cs_pops_N12<-my.spec.pgram(simsRho0$pops[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-cs_pops_N12$freq
y<-Re(cs_pops_N12$pgram[,1,2])
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(0,1.5),main="Pops cospectrum, locs 1 and 2")
allvals<-matrix(NA,length(y),numsims)
allvals[,1]<-y
for (simc in 2:numsims)
{
  cs_pops_N12<-my.spec.pgram(simsRho0$pops[,1:2,simc],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  y<-Re(cs_pops_N12$pgram[,1,2])
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
  allvals[,simc]<-y
}
y1<-rho*(Mod(1/(1-c1*mu-c2*mu^2)))^2
lines(x,y1*y2,col="red") #the analytic expectation
medn<-apply(FUN=median,X=allvals,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()

#show the quadspectrum of the pops for rho=0, compare with theory
pdf(file="Examp1_PopsQuadspectra_rho0.pdf")
cs_pops_N12<-my.spec.pgram(simsRho0$pops[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-cs_pops_N12$freq
y<-Im(cs_pops_N12$pgram[,1,2])
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=c(-.5,.5),main="Pops quadspectrum, locs 1 and 2")
allvals<-matrix(NA,length(y),numsims)
allvals[,1]<-y
for (simc in 2:numsims)
{
  cs_pops_N12<-my.spec.pgram(simsRho0$pops[,1:2,simc],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  y<-Im(cs_pops_N12$pgram[,1,2])
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
  allvals[,simc]<-y
}
lines(x,rep(0,length(x)),col="red") #the analytic expectation
medn<-apply(FUN=median,X=allvals,MARGIN=1)
lines(x,medn,col="red",lty="dashed")
dev.off()
