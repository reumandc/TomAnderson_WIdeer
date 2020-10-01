#This code is for implementing a suggestion of a reviewer, on second revision at Ecology Letters, to explicitly
#decompose the power spectrum of the state-total deer time series into components due to local spectra and 
#between site cospectra. An analogous analysis with DVCs is also done. 

#****
#Some spectral tools
#****

source("./Code/SpectralTools.R")

#***
#Now do the analysis for deer
#***

#pull in the deer data
d<-readRDS(file="Results/cty.list.rds")
deer<-d$Abun
deeryr<-1981:2016

#get the spectral matrix - recall the function myspecmatbrill detrends each time series
specmat_deer<-myspecmatbrill(deer,BiasVariance = 0.25)
freqs<-specmat_deer$freq
specmat_deer<-specmat_deer$specmat
#class(specmat_deer)
#dim(specmat_deer)

#compute the sum of the local spectra
sum_loc_spec<-specmat_deer[1,1,]
for (counter in 2:(dim(deer)[1]))
{
  sum_loc_spec<-sum_loc_spec+specmat_deer[counter,counter,]
}
sum_loc_spec<-Re(sum_loc_spec)

#compute the sum of all the cross spectra
sum_cross_spec<-complex(real=rep(0,dim(specmat_deer)[3]),imaginary=rep(0,dim(specmat_deer)[3]))
for (a in 1:(dim(deer)[1]))
{
  for (b in 1:(dim(deer)[1]))
  {
    if (a!=b)
    {
      sum_cross_spec<-sum_cross_spec+specmat_deer[a,b,]
    }
  }
}
sum_cross_spec<-Re(sum_cross_spec)

#compute the sum of all the spectra and cross spectra and compare to the spectrum of the 
#total time series, as an additional test
spectot1<-sum_loc_spec+sum_cross_spec
spectot2<-myspecmatbrill(matrix(colSums(deer),1,dim(deer)[2]),BiasVariance = 0.25)$specmat
spectot2<-Re(spectot2)
#class(spectot1)
#class(spectot2)
#spectot1
#spectot2
#max(abs(spectot1-spectot2))
#spectot1-spectot2

#do the resampling based on randomizing phases, to get confidence intervals on the sum 
#of the cross spectra
deer_detrend<-matrix(NA,dim(deer)[1],dim(deer)[2])
for (counter in 1:(dim(deer)[1]))
{
  y<-deer[counter,]
  deer_detrend[counter,]<-stats::residuals(stats::lm(y~deeryr))
}
nsurr<-100
deer_s<-wsyn::surrog(deer_detrend,nsurr,"fft",FALSE)
sum_cross_spec_s<-matrix(NA,nsurr,length(sum_cross_spec))
for (scounter in 1:nsurr)
{
  specmat_deer_s<-myspecmatbrill(deer_s[[scounter]],detrend=FALSE,BiasVariance = 0.25)$specmat
  
  h<-complex(real=rep(0,dim(specmat_deer)[3]),imaginary=rep(0,dim(specmat_deer)[3]))
  for (a in 1:(dim(deer)[1]))
  {
    for (b in 1:(dim(deer)[1]))
    {
      if (a!=b)
      {
        h<-h+specmat_deer_s[a,b,]
      }
    }
  }
  sum_cross_spec_s[scounter,]<-Re(h)
}
sum_cross_spec_q<-apply(FUN=quantile,X=sum_cross_spec_s,MARGIN=2,prob=c(0.025,0.975))

#make a plot
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Review2_Fourier1.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(sum_loc_spec,sum_cross_spec,spectot1,freqs,sum_cross_spec_q)
plot(freqs,sum_loc_spec,type="l",col="blue",ylim=ylimits)
lines(freqs,sum_cross_spec,type="l",col="red")
lines(freqs,spectot1,type="l",col="black")
lines(freqs,sum_cross_spec_q[1,],type="l",col="red",lty="dashed")
lines(freqs,sum_cross_spec_q[2,],type="l",col="red",lty="dashed")
lines(rep(1/3,2),ylimits,type="l",lty="dashed")
lines(rep(1/7,2),ylimits,type="l",lty="dashed")
mtext("Frequency",side=1,line=2,cex=1)
mtext("Power",side=2,line=2,cex=1)
text(.5,ylimits[2],"A)",cex=1.5,adj=c(1,1))
dev.off()

#now do an analogue where you only randomize phases in the 3-7-year band

#get surrogates for which only the 3-7-year timescale phases have been randomized
ft_deer_detrend<-myfft(deer_detrend)
tfreqs<-(0:(dim(deer_detrend)[2]-1))/(dim(deer_detrend)[2])
tfreqinds1<-which(tfreqs>=1/7 & tfreqs<=1/3)
tfreqinds2<-rev(dim(deer_detrend)[2]-tfreqinds1+2)  #for keeping the symmetry
deer_s<-list()
for (scounter in 1:nsurr)
{
  h<-ft_deer_detrend
  rphas<-exp(1i*2*pi*runif(length(tfreqinds1)*(dim(h)[1])))
  rphas<-matrix(rphas,dim(h)[1],length(tfreqinds1))
  h[,tfreqinds1]<-h[,tfreqinds1]*rphas
  h[,tfreqinds2]<-h[,tfreqinds2]*Conj(rphas[,seq(from=dim(rphas)[2],to=1,by=-1)])
  h<-imyfft(h)
  if (max(abs(Im(h)))>1e-10){stop("Error in randomizing phases in 3-7-year timescale band")}
  deer_s[[scounter]]<-Re(h)
}

#get the sum of cross spectra for each surrogate
sum_cross_spec_s<-matrix(NA,nsurr,length(sum_cross_spec))
for (scounter in 1:nsurr)
{
  specmat_deer_s<-myspecmatbrill(deer_s[[scounter]],detrend=FALSE,BiasVariance = 0.25)$specmat
  
  h<-complex(real=rep(0,dim(specmat_deer)[3]),imaginary=rep(0,dim(specmat_deer)[3]))
  for (a in 1:(dim(deer)[1]))
  {
    for (b in 1:(dim(deer)[1]))
    {
      if (a!=b)
      {
        h<-h+specmat_deer_s[a,b,]
      }
    }
  }
  sum_cross_spec_s[scounter,]<-Re(h)
}
sum_cross_spec_q<-apply(FUN=quantile,X=sum_cross_spec_s,MARGIN=2,prob=c(0.025,0.975))

#make a plot
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Review2_Fourier2.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(sum_loc_spec,sum_cross_spec,spectot1,freqs,sum_cross_spec_q)
plot(freqs,sum_loc_spec,type="l",col="blue",ylim=ylimits)
lines(freqs,sum_cross_spec,type="l",col="red")
lines(freqs,spectot1,type="l",col="black")
lines(freqs,sum_cross_spec_q[1,],type="l",col="red",lty="dashed")
lines(freqs,sum_cross_spec_q[2,],type="l",col="red",lty="dashed")
lines(rep(1/3,2),ylimits,type="l",lty="dashed")
lines(rep(1/7,2),ylimits,type="l",lty="dashed")
mtext("Frequency",side=1,line=2,cex=1)
mtext("Power",side=2,line=2,cex=1)
text(.5,ylimits[2],"B)",cex=1.5,adj=c(1,1))
dev.off()


#***
#Now do the analysis for dvcs
#***

#pull in the dvc data
d<-readRDS(file="Results/cty.list.rds")
dvcs<-d$Crashes
dvcs<-dvcs[,7:(dim(dvcs)[2])]
dvcyr<-1987:2016

#get the spectral matrix - recall the function myspecmatbrill detrends each time series
specmat_dvcs<-myspecmatbrill(dvcs,BiasVariance = 0.25)
freqs<-specmat_dvcs$freq
specmat_dvcs<-specmat_dvcs$specmat
#class(specmat_dvcs)
#dim(specmat_dvcs)

#compute the sum of the local spectra
sum_loc_spec<-specmat_dvcs[1,1,]
for (counter in 2:(dim(dvcs)[1]))
{
  sum_loc_spec<-sum_loc_spec+specmat_dvcs[counter,counter,]
}
sum_loc_spec<-Re(sum_loc_spec)

#compute the sum of all the cross spectra
sum_cross_spec<-complex(real=rep(0,dim(specmat_dvcs)[3]),imaginary=rep(0,dim(specmat_dvcs)[3]))
for (a in 1:(dim(dvcs)[1]))
{
  for (b in 1:(dim(dvcs)[1]))
  {
    if (a!=b)
    {
      sum_cross_spec<-sum_cross_spec+specmat_dvcs[a,b,]
    }
  }
}
sum_cross_spec<-Re(sum_cross_spec)

#compute the sum of all the spectra and cross spectra and compare to the spectrum of the 
#total time series, as an additional test
spectot1<-sum_loc_spec+sum_cross_spec
spectot2<-myspecmatbrill(matrix(colSums(dvcs),1,dim(dvcs)[2]),BiasVariance = 0.25)$specmat
spectot2<-Re(spectot2)
#class(spectot1)
#class(spectot2)
#spectot1
#spectot2
#max(abs(spectot1-spectot2))
#spectot1-spectot2

#do the resampling based on randomizing phases, to get confidence intervals on the sum 
#of the cross spectra
dvcs_detrend<-matrix(NA,dim(dvcs)[1],dim(dvcs)[2])
for (counter in 1:(dim(dvcs)[1]))
{
  y<-dvcs[counter,]
  dvcs_detrend[counter,]<-stats::residuals(stats::lm(y~dvcyr))
}
nsurr<-100
dvcs_s<-wsyn::surrog(dvcs_detrend,nsurr,"fft",FALSE)
sum_cross_spec_s<-matrix(NA,nsurr,length(sum_cross_spec))
for (scounter in 1:nsurr)
{
  specmat_dvcs_s<-myspecmatbrill(dvcs_s[[scounter]],detrend=FALSE,BiasVariance = 0.25)$specmat
  
  h<-complex(real=rep(0,dim(specmat_dvcs)[3]),imaginary=rep(0,dim(specmat_dvcs)[3]))
  for (a in 1:(dim(dvcs)[1]))
  {
    for (b in 1:(dim(dvcs)[1]))
    {
      if (a!=b)
      {
        h<-h+specmat_dvcs_s[a,b,]
      }
    }
  }
  sum_cross_spec_s[scounter,]<-Re(h)
}
sum_cross_spec_q<-apply(FUN=quantile,X=sum_cross_spec_s,MARGIN=2,prob=c(0.025,0.975))

#make a plot
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Review2_Fourier3.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(sum_loc_spec,sum_cross_spec,spectot1,freqs,sum_cross_spec_q)
plot(freqs,sum_loc_spec,type="l",col="blue",ylim=ylimits)
lines(freqs,sum_cross_spec,type="l",col="red")
lines(freqs,spectot1,type="l",col="black")
lines(freqs,sum_cross_spec_q[1,],type="l",col="red",lty="dashed")
lines(freqs,sum_cross_spec_q[2,],type="l",col="red",lty="dashed")
lines(rep(1/3,2),ylimits,type="l",lty="dashed")
lines(rep(1/7,2),ylimits,type="l",lty="dashed")
mtext("Frequency",side=1,line=2,cex=1)
mtext("Power",side=2,line=2,cex=1)
text(.5,ylimits[2],"A)",cex=1.5,adj=c(1,1))
dev.off()

#now do an analogue where you only randomize phases in the 3-7-year band

#get surrogates for which only the 3-7-year timescale phases have been randomized
ft_dvcs_detrend<-myfft(dvcs_detrend)
tfreqs<-(0:(dim(dvcs_detrend)[2]-1))/(dim(dvcs_detrend)[2])
tfreqinds1<-which(tfreqs>=1/7 & tfreqs<=1/3)
tfreqinds2<-rev(dim(dvcs_detrend)[2]-tfreqinds1+2)  #for keeping the symmetry
dvcs_s<-list()
for (scounter in 1:nsurr)
{
  h<-ft_dvcs_detrend
  rphas<-exp(1i*2*pi*runif(length(tfreqinds1)*(dim(h)[1])))
  rphas<-matrix(rphas,dim(h)[1],length(tfreqinds1))
  h[,tfreqinds1]<-h[,tfreqinds1]*rphas
  h[,tfreqinds2]<-h[,tfreqinds2]*Conj(rphas[,seq(from=dim(rphas)[2],to=1,by=-1)])
  h<-imyfft(h)
  if (max(abs(Im(h)))>1e-10){stop("Error in randomizing phases in 3-7-year timescale band, loc 2")}
  dvcs_s[[scounter]]<-Re(h)
}

#get the sum of cross spectra for each surrogate
sum_cross_spec_s<-matrix(NA,nsurr,length(sum_cross_spec))
for (scounter in 1:nsurr)
{
  specmat_dvcs_s<-myspecmatbrill(dvcs_s[[scounter]],detrend=FALSE,BiasVariance = 0.25)$specmat
  
  h<-complex(real=rep(0,dim(specmat_dvcs)[3]),imaginary=rep(0,dim(specmat_dvcs)[3]))
  for (a in 1:(dim(dvcs)[1]))
  {
    for (b in 1:(dim(dvcs)[1]))
    {
      if (a!=b)
      {
        h<-h+specmat_dvcs_s[a,b,]
      }
    }
  }
  sum_cross_spec_s[scounter,]<-Re(h)
}
sum_cross_spec_q<-apply(FUN=quantile,X=sum_cross_spec_s,MARGIN=2,prob=c(0.025,0.975))

#make a plot
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/Review2_Fourier4.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(sum_loc_spec,sum_cross_spec,spectot1,freqs,sum_cross_spec_q)
plot(freqs,sum_loc_spec,type="l",col="blue",ylim=ylimits)
lines(freqs,sum_cross_spec,type="l",col="red")
lines(freqs,spectot1,type="l",col="black")
lines(freqs,sum_cross_spec_q[1,],type="l",col="red",lty="dashed")
lines(freqs,sum_cross_spec_q[2,],type="l",col="red",lty="dashed")
lines(rep(1/3,2),ylimits,type="l",lty="dashed")
lines(rep(1/7,2),ylimits,type="l",lty="dashed")
mtext("Frequency",side=1,line=2,cex=1)
mtext("Power",side=2,line=2,cex=1)
text(.5,ylimits[2],"B)",cex=1.5,adj=c(1,1))
dev.off()

