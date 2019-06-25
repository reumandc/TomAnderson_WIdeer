#***libraries and other pre work

rm(list=ls())
library(mvtnorm)
source("my.spec.pgram.R")

#***settings

N<-70 #number of sampling locations
lensim<-1024 #length of simulations
lensimshort<-100

#***work on a first example with AR(1) dynamics, noise with weak periodicity 
#and synchrony (which is therefore timescale-specific synchrony)

#generate the inovations for the generation of the noise
rho<-0.7
Sig<-matrix(rho,N,N)
diag(Sig)<-1
delta<-rmvnorm(3*lensim,sigma=Sig)

#now use these to generate the noise

#coefficients for the AR(2) model used to generate the environmental noise, option 1
r<-3 
theta<-pi/2

#coefficients for the AR(2) model used to generate the environmental noise, option 2
#r<-2
#theta<-pi/3

#coefficients for the AR(2) model used to generate the environmental noise, option 3
#r<-2
#theta<-2*pi/3

#get c1 and c2
c1<-(1/r)*2*cos(theta)
c2<-(-1/r^2) 

epsilon<-matrix(0,3*lensim,N)
for (counter in 3:(3*lensim))
{
  epsilon[counter,]<-c1*epsilon[counter-1,]+c2*epsilon[counter-2,]+delta[counter,]
}
epsilon<-epsilon[(dim(epsilon)[1]-2*lensim+1):(dim(epsilon)[1]),]

#now put the noise into an AR(1) to get pops
pops<-matrix(0,2*lensim,N)
popcoef<-0.15
for (counter in 2:(2*lensim))
{
  pops[counter,]<-popcoef*pops[counter-1,]+epsilon[counter,]
}
pops<-pops[(dim(pops)[1]-lensim+1):(dim(pops)[1]),]
epsilon<-epsilon[(dim(epsilon)[1]-lensim+1):(dim(epsilon)[1]),]

#make a plots that show the analytic and estimated spectra of noise 
spec_noise<-my.spec.pgram(epsilon,spans=c(91,71),plot=FALSE)
x<-spec_noise$freq
mu<-exp(-2*pi*complex(real=0,imaginary=1)*x)
y1<-(Mod(1/(1-c1*mu-c2*mu^2)))^2
plot(x,spec_noise$spec[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(spec_noise$spec,y1),
     main="Noise spectra") #plots for all locations
for (loc in 2:N)
{
  lines(x,spec_noise$spec[,loc],type="l",col=rgb(0,0,0,alpha=.1))
}
lines(x,y1,type="l",col="red")
x[which(y1==max(y1))] #where is the peak of the analytic result?

#make a few plots that show the cospectra for noise
plot(x,Re(spec_noise$pgram[,1,2]),col=rgb(0,0,0,alpha=0),type='l',
     ylim=range(Re(spec_noise$pgram)),
     main="Noise cospectra")
for (count1 in 1:(N-1))
{
  for (count2 in (count1+1):N)
  {
    lines(x,Re(spec_noise$pgram[,count1,count2]),col=rgb(0,0,0,alpha=0.01),type='l')
  }
}
lines(x,rho*y1,type="l",col="red") #***multiplying y1 by rho might be not-quite-right

#make a few plots that show the quadspectra for noise, should be close to 0
plot(x,Im(spec_noise$pgram[,1,2]),col=rgb(0,0,0,alpha=0),type='l',
     ylim=range(Im(spec_noise$pgram)),
     main="Noise quadspectra")
for (count1 in 1:(N-1))
{
  for (count2 in (count1+1):N)
  {
    lines(x,Im(spec_noise$pgram[,count1,count2]),col=rgb(0,0,0,alpha=0.01),type='l')
  }
}
lines(x,rep(0,length(x)),type='l',col='red')

#make a few plots that show the spectra of pops 
spec_pops<-my.spec.pgram(pops,spans=c(91,71),plot=FALSE)
y1p<-y1*(Mod(1/(1-popcoef*mu)))^2
plot(x,spec_pops$spec[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(spec_pops$spec,y1p),
     main="Pops spectra") #plots for all locations
for (loc in 2:N)
{
  lines(x,spec_pops$spec[,loc],type="l",col=rgb(0,0,0,alpha=.1))
}
lines(x,y1p,type="l",col="red")
x[which(y1p==max(y1p))] #where is the peak of the analytic result?

#cospectra of pops
plot(x,Re(spec_pops$pgram[,1,2]),col=rgb(0,0,0,alpha=0),type='l',
     ylim=range(Re(spec_pops$pgram)),
     main="Pops cospectra")
for (count1 in 1:(N-1))
{
  for (count2 in (count1+1):N)
  {
    lines(x,Re(spec_pops$pgram[,count1,count2]),col=rgb(0,0,0,alpha=0.01),type='l')
  }
}
lines(x,rho*y1p,type="l",col="red") #***multiplying y1 by rho might be not-quite-right

#quadspectra of pops
plot(x,Im(spec_pops$pgram[,1,2]),col=rgb(0,0,0,alpha=0),type='l',
     ylim=range(Im(spec_pops$pgram)),
     main="Pops quadspectra")
for (count1 in 1:(N-1))
{
  for (count2 in (count1+1):N)
  {
    lines(x,Im(spec_pops$pgram[,count1,count2]),col=rgb(0,0,0,alpha=0.01),type='l')
  }
}
lines(x,rep(0,length(x)),type='l',col='red')

#plot showing the spectrum of the total pop
totpop<-apply(FUN=sum,X=pops,MARGIN=1)
plot(totpop[1:100],type='l',cex=.5)
points(totpop[1:100],type='p',cex=.5,pch=20)
spec_totpop<-my.spec.pgram(totpop,spans=c(91,71),plot=FALSE)
plot(x,spec_totpop$spec,type="l",
     #ylim=range(spec_totpop$spec,y1p),
     main="totpop spectrum") 
