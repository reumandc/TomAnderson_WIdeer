#***
#Comparison of Fourier transforms of DVCs and the surrogates developed for fig. 5
#***

#pull in the dvc data
d<-readRDS(file="Results/cty.list.rds")
dvcs<-d$Crashes
dvcs<-dvcs[,7:(dim(dvcs)[2])]
totdvcs<-apply(FUN=sum,X=dvcs,MARGIN=2)
dvcyr<-1987:2016

#pull in the Fourier surrogates developed for Fig 5
dvcsurr<-read.csv("Data/dvcsurrsum.csv")
dvcsurr<-as.matrix(dvcsurr)

#get raw periodograms of data and surrogates
dvcspecraw<-myspecraw(totdvcs)
freqs<-dvcspecraw$freq
dvcspecraw<-dvcspecraw$spec
surrspecsraw<-matrix(NA,dim(dvcsurr)[1],length(freqs))
for (counter in 1:dim(dvcsurr)[1])
{
  x<-unname(dvcsurr[counter,])
  h<-myspecraw(x)
  surrspecsraw[counter,]<-h$spec
}

#sum over 3-7 year timescales, compare data and surrogate sums using a histogram
realsum<-sum(dvcspecraw[freqs>=1/7 & freqs<=1/3])
surrsums<-apply(FUN=sum,X=surrspecsraw[,freqs>=1/7 & freqs<=1/3],MARGIN=1)
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Fourier11.png",res=600,units="in",width = tot.wd,height = tot.ht)
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
y<-c(totdvcs,t(dvcsurr))
x<-rep(dvcyr,times=1+dim(dvcsurr)[1])
dtmod<-lm(y~x)
totdvcres<-totdvcs-coef(dtmod)[2]*dvcyr-coef(dtmod)[1]
dvcsurrres<-dvcsurr-coef(dtmod)[2]*matrix(rep(dvcyr,each=dim(dvcsurr)[1]),dim(dvcsurr)[1],length(dvcyr))-coef(dtmod)[1]

#Now get Brillinger consistent estimators of the spectra of data and surrogates
dvclog10specbrill<-myspecbrill(totdvcres,detrend=FALSE)
freqs<-dvclog10specbrill$freq
dvclog10specbrill<-dvclog10specbrill$log10spec
surrlog10specsbrill<-matrix(NA,dim(dvcsurr)[1],length(freqs))
for (counter in 1:dim(dvcsurr)[1])
{
  x<-unname(dvcsurrres[counter,])
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
png("Results/Fourier12.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(10^dvclog10specbrill,10^surrlog10specsbrill)
plot(freqs,10^dvclog10specbrill,type="l",lwd=2,ylim=ylimits)
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
#I decided to stop here for DVCs for two reasons: 1) the above is the "right" analysis,
#and 2) the other analyses I did for deer are even less appropriate for DVCs. For
#instance, the AR(1) surrogate analysis for deer does not make much sense for DVCs
#since there is low-frequency content in the state-total DVC time series. This either
#has to be removed (sacrifices statistical power and it would be ad hoc how to do - 
#maybe remove a cubic but that's pretty ad hoc), or else when you fit an AR(1) it's
#going to be heavily influenced by the low-frequency content, which means the AR(1)
#surrogates are going to be bad surrogates. The AIC-based AR argument I did for deer
#also does not make sense for DVCs for similar reasons. The comparison between state-
#and county-level analyses for deer needed detrending and variance standardization
#for deer, and again not clear how to do that for DVCs so it gets messy. Throw on 
#top of this the fact that all these analyses were just things we did to answer the 
#referee, and were predicated on the referee's faulty analysis. So it makes sense
#to just stick with the above for DVCs.
#
#I furthermore hesitate to do the county to state-level comparison because that comparison
#relies on detrending and standardizing the variance of both the state- and county-level
#time series. The variance standardization is going to mean we look at the proportion of
#variance in the 3-7-year timescale band, and this is also influenced by the amount
#of low-frequency variation (since it is a proportion). In fact since there is a lot of 
#low-frequency variation, it is probably influenced a lot by it. I think figuring out how 
#to improve on this analysis eventually leads right back to the analysis where one 
#compared the state-level spectral power in the 3-7-year band to what it would be if 
#there were no synchrony.
#***

