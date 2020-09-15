
#****
#Some spectral tools
#****

source("./Code/SpectralTools.R")

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
deerlog10specbrill<-myspecbrill(totdeerres,detrend=FALSE,BiasVariance = 0.25)
freqs<-deerlog10specbrill$freq
deerlog10specbrill<-deerlog10specbrill$log10spec
surrlog10specsbrill<-matrix(NA,dim(abunsurr)[1],length(freqs))
for (counter in 1:dim(abunsurr)[1])
{
  x<-unname(abunsurrres[counter,])
  h<-myspecbrill(x,detrend=FALSE,BiasVariance = 0.25)
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
#Now compare to AR(1) surrogates (deer again)
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
deerlog10specbrill<-myspecbrill(totdeer,BiasVariance = 0.25)
freqs<-deerlog10specbrill$freq
deerlog10specbrill<-deerlog10specbrill$log10spec
surrlog10specsbrill<-matrix(NA,dim(abunsurr)[1],length(freqs))
for (counter in 1:dim(abunsurr)[1])
{
  x<-unname(abunsurr[counter,])
  h<-myspecbrill(x,BiasVariance = 0.25)
  surrlog10specsbrill[counter,]<-h$log10spec
}

#Now make a plot with no multiple-testing correction
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier3p5.png",res=600,units="in",width = tot.wd,height = tot.ht)
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
#lines(freqs,10^qsurrspec[2,],type="l",lty="dotted") #Bonferroni corrected threshold
mtext("Frequency",side=1,line=2,cex=1)
mtext("Power spectrum",side=2,line=2,cex=1)
dev.off()

#Now make a plot with one type of multiple-testing correction
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

#Now make a plot that includes the strictest multiple-testing correction,
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

#***
#Now do an AR argument with AIC
#***

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

#***
#Now just show the time series, which visually look like they have cycling
#***

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

#***
#Now compare spectra of deer for individual counties with the spectrum for the whole state deer time series
#***

#pull in the deer data
d<-readRDS(file="Results/cty.list.rds")
deer<-d$Abun
totdeer<-apply(FUN=sum,X=deer,MARGIN=2)
deeryr<-1981:2016

#get standardized time series (detrend, demean, standardized variance)
totdeerSV<-residuals(lm(totdeer~deeryr))
totdeerSV<-totdeerSV/sd(totdeerSV)
deerSV<-matrix(NA,dim(deer)[1],dim(deer)[2])
for (counter in 1:dim(deer)[1])
{
  countydeer<-deer[counter,]
  countydeerSV<-residuals(lm(countydeer~deeryr))
  countydeerSV<-countydeerSV/sd(countydeerSV)
  deerSV[counter,]<-countydeerSV
}

#now compute the spectral power in the 3-7 year band for each time series
totspecraw<-myspecraw(totdeerSV)
freqs<-totspecraw$freq
totspecraw37<-sum(totspecraw$spec[freqs>=1/7 & freqs<=1/3])
countyspecraw37<-NA*numeric(dim(deerSV)[1])
for (counter in 1:length(countyspecraw37))
{
  h<-myspecraw(deerSV[counter,])
  countyspecraw37[counter]<-sum(h$spec[freqs>=1/7 & freqs<=1/3])
}
#sum(totspecraw37>=countyspecraw37)
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier8.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
hist(countyspecraw37,xlim=range(countyspecraw37,totspecraw37),main="")
points(totspecraw37,0,col="red",pch=20)
mtext("Power in 3-7 yr band",side=1,line=2,cex=1)
mtext("Count",side=2,line=2,cex=1)
dev.off()
#so there is more power in the 3-7yr band for the standardized state-total
#time series than there is for any of the county-level time series, where
#the standardization (detrend and demean, standard variance) means that 
#what we are seeing here is an appropriate comparison of periodicity in
#that band.

#Now look directly at the spectra. Compute the spectra (Brillinger consistent estimator).
#Again using the standardized time series.
totspec<-myspecbrill(totdeerSV,detrend=FALSE,BiasVariance = 0.25)
countyspec<-list()
for (counter in 1:dim(deer)[1])
{
  countyspec[[counter]]<-myspecbrill(deerSV[counter,],detrend=FALSE,BiasVariance = 0.25)
}

#now isolate what you want to plot
freqs<-totspec$freq
yvals<-matrix(NA,length(countyspec)+1,length(freqs))
yvals[1,]<-10^totspec$log10spec
for (counter in 1:length(countyspec))
{
  yvals[counter+1,]<-10^countyspec[[counter]]$log10spec
}

#now plot it
png("Results/Fourier9.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(yvals)
plot(freqs,yvals[1,],type="l",lwd=2,ylim=ylimits,col="red",xlab="",ylab="")
for (counter in 2:(dim(yvals)[1]))
{
  lines(freqs,yvals[counter,])
}
lines(rep(1/3,2),ylimits,lty="dashed")
lines(rep(1/7,2),ylimits,lty="dashed")
mtext("Power spectrum",side=1,line=2,cex=1)
mtext("Frequency",side=2,line=2,cex=1)
dev.off()

#one of the counties has a peak at a similar location, similar height to the total,
#let's investigate that county
ind<-which(yvals[2:dim(yvals)[1],7]==max(yvals[2:dim(yvals)[1],7]))
plot(freqs,yvals[1,],type="l",col="red",ylim=range(yvals[1,],yvals[ind+1,]))
lines(freqs,yvals[ind+1,])
mean(totdeer)/mean(deer[38,]) #So the mean pop in this county is less than 1/28th 
#of the state total. This is more than an average county, but not enough to
#be driving things by itself.

#to make sure, however, remove that county and repeat one of the analyses above
deer<-deer[-38,]
totdeer<-apply(FUN=sum,X=deer,MARGIN=2)

#get standardized time series (detrend, demean, standardized variance)
totdeerSV<-residuals(lm(totdeer~deeryr))
totdeerSV<-totdeerSV/sd(totdeerSV)
deerSV<-matrix(NA,dim(deer)[1],dim(deer)[2])
for (counter in 1:dim(deer)[1])
{
  countydeer<-deer[counter,]
  countydeerSV<-residuals(lm(countydeer~deeryr))
  countydeerSV<-countydeerSV/sd(countydeerSV)
  deerSV[counter,]<-countydeerSV
}

#now compute the spectral power in the 3-7 year band for each time series
totspecraw<-myspecraw(totdeerSV)
freqs<-totspecraw$freq
totspecraw37<-sum(totspecraw$spec[freqs>=1/7 & freqs<=1/3])
countyspecraw37<-NA*numeric(dim(deerSV)[1])
for (counter in 1:length(countyspecraw37))
{
  h<-myspecraw(deerSV[counter,])
  countyspecraw37[counter]<-sum(h$spec[freqs>=1/7 & freqs<=1/3])
}
#sum(totspecraw37>=countyspecraw37)
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier10.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
hist(countyspecraw37,xlim=range(countyspecraw37,totspecraw37),main="")
points(totspecraw37,0,col="red",pch=20)
mtext("Power in 3-7 yr band",side=1,line=2,cex=1)
mtext("Count",side=2,line=2,cex=1)
dev.off()
#so there is STILL more power in the 3-7yr band for the standardized state-total
#time series than there is for any of the county-level time series

#Some text in an email, saved here for later in case needed
#Hi guys,
#
#The referee also asked us to check the spectral peak for the state-total deer time series in the 3-7 year band was 
#bigger than that for county time series.
#I first detrended and demeaned all county-level deer time series, and standardized their variance, and did the same 
#to the state-total time series. Then I computed total spectral power in the 3-7yr band for all these time series. 
#Fourier8.png, attached, shows the histogram, of county values, with the red dot being the state-total value. The 
#fact that I standardized all these time series means this is a valid comparison – there is more periodicity in the 
#3-7yr timescale band for the state total than for any of the counties.
#Fourier9.png, attached, shows the spectra for the counties (black) and for the state total, again using standardized 
#time series for a valid comparison. You can see the red peak in 3-7 is clearly above the black lines, except for one. 
#That one is Marquette county. The mean deer population in Marquette is less than 1/28th the state total mean deer 
#population, so Marquette could not reasonable have been driving things by itself. To be on the safe side, I removed 
#Marquette, recomputed the state total, re-standardized time series, and then recomputed total power in 3-7. The result 
#is Fourier10.png, attached, which shows the same pattern again.
#So there is definitely more periodicity in the state total than in the counties.
#I have not added this to any docs, I’ve just added the code that makes these results to the code suite. I imagine the best thing to do is probably just to add the result to the Sup Mat that there is more spectral power in 3-7 for the state total than for any of the counties, using standardized time series. And then add the other details to the response to referees. I’ll take a look at where Tom has got to in his edits and if it seems appropriate I’ll add this new material as indicated, probably later today. Before that, though, I am going to do some analyses of DVC data, similar to what I have already done for deer.
#Dan
