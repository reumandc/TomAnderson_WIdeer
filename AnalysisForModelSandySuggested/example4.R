#***libraries and other pre work

rm(list=ls())
library(mvtnorm)
source("my.spec.pgram.R")

#***settings

N<-71 #number of sampling locations
lensim<-1024 #length of simulations
lensimshort<-128

#***the models

#Simulates the model for case a
#
#Args
#N        The number of sampling locations
#w        The period of the oscillations
#rho      Written as rho_a in my notes on p. 14 and before that
#rhopop   Lag-1 autocorrelation coefficient for population dynamics
#lensim   Length of simulations to do
#numsims  Number of simulations to do
#
#Output - a list with these elements 
#epsilon1 - a lensim by N by numsims array, the first environmental noise (the periodic one)
#epsilon2 - a lensim by N by numsims array, the second environmental noise (the nonperiodic one)
#pops - a lensim by N by numsims array
simmodel4a<-function(N,w,rho,rhopop,lensim,numsims)
{
  #**generate the first environmental noise (the periodic one)
  tims<-seq(from=0,to=2*lensim-1,by=1)
  epsilon1<-sin((2*pi/w)*array(tims,dim=c(2*lensim,N,numsims))+
                  array(2*pi*rep(runif(N*numsims),each=2*lensim),dim=c(2*lensim,N,numsims)))

  #generate the second environmental noise (the non-periodic one)
  Sig<-matrix(rho,N,N)
  diag(Sig)<-1
  epsilon2<-rmvnorm(2*lensim*numsims,sigma=Sig)
  dim(epsilon2)<-c(2*lensim,numsims,N)
  epsilon2<-aperm(epsilon2,c(1,3,2))

  #**now add the two environmental noises and use them to simulate pops

  #do the sims  
  pops<-array(0,dim=c(2*lensim,N,numsims))
  for (counter in 2:(2*lensim))
  {
    pops[counter,,]<-rhopop*pops[counter-1,,]+epsilon1[counter,,]+epsilon2[counter,,]
  }
  
  #throw away transients  
  pops<-pops[(dim(pops)[1]-lensim+1):(dim(pops)[1]),,]
  epsilon1<-epsilon1[(dim(epsilon1)[1]-lensim+1):(dim(epsilon1)[1]),,]
  epsilon2<-epsilon2[(dim(epsilon2)[1]-lensim+1):(dim(epsilon2)[1]),,]

  return(list(epsilon1=epsilon1,
              epsilon2=epsilon2,
              pops=pops))
}

#Simulates the model for cases b, c
#
#Args
#N        The number of sampling locations
#w        The period of the oscillations
#rho      Written as rho_bc in my notes on p. 14 and before that
#rhopop   Lag-1 autocorrelation coefficient for population dynamics
#lensim   Length of simulations to do
#numsims  Number of simulations to do
#
#Output - a list with these elements 
#epsilon1 - a lensim by N by numsims array, the first environmental noise (the periodic one)
#epsilon2 - a lensim by N by numsims array, the second environmental noise (the nonperiodic one)
#pops - a lensim by N by numsims array
simmodel4bc<-function(N,w,rho,rhopop,lensim,numsims)
{
  #**generate the first environmental noise (the periodic one)
  tims<-seq(from=0,to=2*lensim-1,by=1)
  epsilon1<-sin((2*pi/w)*array(tims,dim=c(2*lensim,N,numsims))+
                  array(rep(2*pi*runif(numsims),each=2*lensim*N),dim=c(2*lensim,N,numsims)))
  
  #generate the second environmental noise (the non-periodic one)
  Sig<-matrix(rho,N,N)
  diag(Sig)<-1
  epsilon2<-rmvnorm(2*lensim*numsims,sigma=Sig)
  dim(epsilon2)<-c(2*lensim,numsims,N)
  epsilon2<-aperm(epsilon2,c(1,3,2))
  
  #**now add the two environmental noises and use them to simulate pops
  
  #do the sims  
  pops<-array(0,dim=c(2*lensim,N,numsims))
  for (counter in 2:(2*lensim))
  {
    pops[counter,,]<-rhopop*pops[counter-1,,]+epsilon1[counter,,]+epsilon2[counter,,]
  }
  
  #throw away transients  
  pops<-pops[(dim(pops)[1]-lensim+1):(dim(pops)[1]),,]
  epsilon1<-epsilon1[(dim(epsilon1)[1]-lensim+1):(dim(epsilon1)[1]),,]
  epsilon2<-epsilon2[(dim(epsilon2)[1]-lensim+1):(dim(epsilon2)[1]),,]
  
  return(list(epsilon1=epsilon1,
              epsilon2=epsilon2,
              pops=pops))
}

#***do the sims

#do the sims case a
wa<-3
rhoa<-0 #use rhoa=0.75 and rhobc (below) equal to 0.25 for an example with equal 
#covariances (in expectation) between locations, just use rhoa=rhobc if you just
#want to turn on and off the periodic component of synchrony while leaving all else
#the same
rhopop<-0.15
numsims<-100
set.seed(101)
sims_casea<-simmodel4a(N,wa,rhoa,rhopop,lensim,numsims)

#do the sims case b
wb<-3
rhobc<-0
rhopop<-0.15
numsims<-100
set.seed(101)
sims_caseb<-simmodel4bc(N,wb,rhobc,rhopop,lensim,numsims)

#do the sims case c
wc<-10
rhopop<-0.15
numsims<-100
set.seed(101)
sims_casec<-simmodel4bc(N,wc,rhobc,rhopop,lensim,numsims)

#***check that noise variances are the same in expectation in the three cases
allvars_casea<-NA*numeric(numsims)
allvars_caseb<-NA*numeric(numsims)
allvars_casec<-NA*numeric(numsims)
epsilon_csa<-sims_casea$epsilon1+sims_casea$epsilon2
epsilon_csb<-sims_caseb$epsilon1+sims_caseb$epsilon2
epsilon_csc<-sims_casec$epsilon1+sims_casec$epsilon2
for (counter in 1:numsims)
{
  allvars_casea[counter]<-var(epsilon_csa[,1,counter])
  allvars_caseb[counter]<-var(epsilon_csb[,1,counter])
  allvars_casec[counter]<-var(epsilon_csc[,1,counter])
}
hist(allvars_casea)
hist(allvars_caseb)
hist(allvars_casec)
mean(allvars_casea)
mean(allvars_caseb)
mean(allvars_casec)
#these are very similar, so it seems to be working

#***check that noise covariances between locations 1 and 2 are the same in expectation 
#in the three cases - this is only true if rhoa=0.5+rhobc

allcovs_casea<-NA*numeric(numsims)
allcovs_caseb<-NA*numeric(numsims)
allcovs_casec<-NA*numeric(numsims)
for (counter in 1:numsims)
{
  allcovs_casea[counter]<-cov(epsilon_csa[,1,counter],epsilon_csa[,2,counter])
  allcovs_caseb[counter]<-cov(epsilon_csb[,1,counter],epsilon_csb[,2,counter])
  allcovs_casec[counter]<-cov(epsilon_csc[,1,counter],epsilon_csc[,2,counter])
}
hist(allcovs_casea)
hist(allcovs_caseb)
hist(allcovs_casec)
mean(allcovs_casea)
mean(allcovs_caseb)
mean(allcovs_casec)

#**show the spectra of total populations in the three cases

#get ready
totpops_csa<-apply(FUN=sum,X=sims_casea$pops,MARGIN=c(1,3))
totpops_csb<-apply(FUN=sum,X=sims_caseb$pops,MARGIN=c(1,3))
totpops_csc<-apply(FUN=sum,X=sims_casec$pops,MARGIN=c(1,3))
specs_csa<-my.spec.pgram(totpops_csa,detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
specs_csb<-my.spec.pgram(totpops_csb,detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
specs_csc<-my.spec.pgram(totpops_csc,detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
x<-specs_csa$freq
y_csa<-specs_csa$spec
y_csb<-specs_csb$spec
y_csc<-specs_csc$spec

#make the plot for case a, numerics
pdf(file="Examp4_TotpopsSpectra_csa.pdf")
plot(x,y_csa[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(y_csa,y_csb,y_csc),
     xlab="Frequency",ylab="Power",main="Spectra, total pop, cs a")
for (counter in 1:numsims)
{
  lines(x,y_csa[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=y_csa,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#construct a replacement for the impulse function for plotting, for case a
basea<-0.15
triareaa<-N/4
heighta<-2*triareaa/basea
alphaa<-heighta/basea
impa<-NA*length(x)
impa[x<1/wa-basea/2]<-0
impa[1/wa-basea/2<=x & x<=1/wa]<-2*alphaa*(x[1/wa-basea/2<=x & x<=1/wa]-(1/wa-basea/2))
impa[1/wa<=x & x<=1/wa+basea/2]<-2*alphaa*(1/wa+basea/2-x[1/wa<=x & x<=1/wa+basea/2])
impa[x>1/wa+basea/2]<-0
#plot(x,impa,type="l")

#add the theory
mu<-exp(-2*pi*complex(real=0,imaginary=1)*x)
y_csa_an<-(Mod(1/(1-rhopop*mu)))^2*
  (N+impa+(N^2-N)*rhoa)
lines(x,y_csa_an,col="red")
#lines(c(1/3,1/3),c(0,1000))
dev.off()

#make the plot for case b, y axis extent the same
pdf(file="Examp4_TotpopsSpectra_csb.pdf")
plot(x,y_csb[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(y_csa,y_csb,y_csc),
     xlab="Frequency",ylab="Power",main="Spectra, total pop, cs b")
for (counter in 1:numsims)
{
  lines(x,y_csb[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=y_csb,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#construct a replacement for the impulse function for plotting, for case b
baseb<-basea
triareab<-N^2/4
heightb<-2*triareab/baseb
alphab<-heightb/baseb
impb<-NA*length(x)
impb[x<1/wb-baseb/2]<-0
impb[1/wb-baseb/2<=x & x<=1/wb]<-2*alphab*(x[1/wb-baseb/2<=x & x<=1/wb]-(1/wb-baseb/2))
impb[1/wb<=x & x<=1/wb+baseb/2]<-2*alphab*(1/wb+baseb/2-x[1/wb<=x & x<=1/wb+baseb/2])
impb[x>1/wb+baseb/2]<-0
#plot(x,impb,type="l")

#add the theory
y_csb_an<-(Mod(1/(1-rhopop*mu)))^2*
  (N+impb+(N^2-N)*rhobc)
lines(x,y_csb_an,col="red")
lines(c(1/3,1/3),c(0,1000))
dev.off()

#make the plot for case c, y axis extent the same
pdf(file="Examp4_TotpopsSpectra_csc.pdf")
plot(x,y_csc[,1],type="l",col=rgb(0,0,0,alpha=.1),
     ylim=range(y_csa,y_csb,y_csc),
     xlab="Frequency",ylab="Power",main="Spectra, total pop, cs c")
for (counter in 1:numsims)
{
  lines(x,y_csc[,counter],col=rgb(0,0,0,alpha=.1))
}
mdn<-apply(FUN=median,X=y_csc,MARGIN=1)
lines(x,mdn,col="red",lty="dashed")

#construct a replacement for the impulse function for plotting, for case c
basec<-basea
triareac<-N^2/4
heightc<-2*triareac/basec
alphac<-heightc/basec
impc<-NA*length(x)
impc[x<1/wc-basec/2]<-0
impc[1/wc-basec/2<=x & x<=1/wc]<-2*alphac*(x[1/wc-basec/2<=x & x<=1/wc]-(1/wc-basec/2))
impc[1/wc<=x & x<=1/wc+basec/2]<-2*alphac*(1/wc+basec/2-x[1/wc<=x & x<=1/wc+basec/2])
impc[x>1/wc+basec/2]<-0
#plot(x,impc,type="l")

#add the theory
y_csc_an<-(Mod(1/(1-rhopop*mu)))^2*
  (N+impc+(N^2-N)*rhobc)
lines(x,y_csc_an,col="red")
dev.off()

#plot the total pop versus time for the three cases
x<-101:151
xp<-x-min(x)+1
ycsa<-totpops_csa[x,2]
ycsb<-totpops_csb[x,2]
ycsc<-totpops_csc[x,2]
pdf(file="Examp4_Totpopsa.pdf")
plot(xp,ycsa,type='l',ylim=range(ycsa,ycsb,ycsc))
points(xp,ycsa,type='p',pch=20,cex=.5)
dev.off()
pdf(file="Examp4_Totpopsb.pdf")
plot(xp,ycsb,type='l',col='green',ylim=range(ycsa,ycsb,ycsc))
points(xp,ycsb,type='p',pch=20,cex=.5,col='green')
dev.off()
pdf(file="Examp4_Totpopsc.pdf")
plot(xp,ycsc,type='l',col='red',ylim=range(ycsa,ycsb,ycsc))
points(xp,ycsc,type='p',pch=20,cex=.5,col='red')
dev.off()

var(ycsa)
var(ycsb)
var(ycsc)


