#This script shows, using Fourier methods only, that patterns of synchrony in deer and DVCs are timescale specific.
#This is part of the comprehensive Fourier-only analysis compiled after second review at Ecology Letters.
#Reuman

#***
#functions
#***

source("./Code/SpectralTools.R")

#Computes the Fourier metric for timescale-specific synchrony
#
#Args
#dat        A matrix with time series along the rows. Assumed cleaned, including detrending.
#
#Output - a named list with these elements
#freq     frequencies
#sync     synchrony at those frequencies
#
FSsync<-function(dat)
{
  sres<-myspecmatbrill(dat,detrend=FALSE,BiasVariance = 0.25)
  freq<-sres$freq
  sres<-sres$specmat
  sync<-complex(rep(0,length(freq)),rep(0,length(freq)))
  dd<-dim(dat)[1]
  for (ci in 1:dd)
  {
    for (cj in 1:dd)
    {
      if (ci != cj)
      {
        sync<-sync+sres[ci,cj,]
      }
    }
  }
  sync<-sync/(dd*(dd-1))
  
  if (any(Im(sync)>1e-10))
  {
    stop("Error in FSsync: sync was not real")
  }
  sync<-Re(sync)
  
  return(list(freq=freq,sync=sync))
}

#***
#Load the right data and clean it
#***

#pull in the raw (un-transformed) deer data
d<-readRDS(file="Results/cty.list.rds")
deer<-d$Abun
deeryr<-1981:2016

#pull in the dvc data
dvcs<-d$Crashes
dvcs<-dvcs[,7:(dim(dvcs)[2])]
dvcyr<-1987:2016

#clean and transform
deer<-wsyn::cleandat(deer,clev=5,times=deeryr)$cdat
dvcs<-wsyn::cleandat(dvcs,clev=5,times=dvcyr)$cdat

#***
#Compute the Fourier metric for timescale specific synchrony
#***

#first do deer
deersync<-FSsync(deer)
freq_deer<-deersync$freq
deersync<-deersync$sync

#now do dvcs
dvcsync<-FSsync(dvcs)
freq_dvcs<-dvcsync$freq
dvcsync<-dvcsync$sync

#***
#Now do the same for surrogates - these are phase randomizations, independent for each sampling locations
#***

#make the surrogates
nsurrog<-10000
deer_s<-wsyn::surrog(deer,nsurrog,"fft",FALSE)
dvcs_s<-wsyn::surrog(dvcs,nsurrog,"fft",FALSE)

#get synchrony for each
deersync_s<-matrix(NA,nsurrog,length(deersync))
dvcsync_s<-matrix(NA,nsurrog,length(dvcsync))
for (counter in 1:nsurrog)
{
  #print(paste0("Surrogates 1, ",counter," of ",nsurrog))
  
  h<-deer_s[[counter]]
  h<-wsyn::cleandat(h,deeryr,clev=2)$cdat #detrend, so surrogates are treated the same way as data
  h<-FSsync(h)
  deersync_s[counter,]<-h$sync
  
  h<-dvcs_s[[counter]]
  h<-wsyn::cleandat(h,dvcyr,clev=2)$cdat
  h<-FSsync(h)
  dvcsync_s[counter,]<-h$sync
}

#***
#Now do the same for a second type of surrogates - these permute time 
#I ended up deciding these are not so useful
#***

#make the surrogates
deer_s2<-list()
dvcs_s2<-list()
for (counter in 1:nsurrog)
{
  deer_s2[[counter]]<-deer[,sample.int(length(deeryr),length(deeryr))]
  dvcs_s2[[counter]]<-dvcs[,sample.int(length(dvcyr),length(dvcyr))]
}

#get synchrony for each
deersync_s2<-matrix(NA,nsurrog,length(deersync))
dvcsync_s2<-matrix(NA,nsurrog,length(dvcsync))
for (counter in 1:nsurrog)
{
  print(paste0("Surrogates 2, ",counter," of ",nsurrog))

  h<-deer_s2[[counter]]
  h<-wsyn::cleandat(h,deeryr,clev=2)$cdat #detrend, so surrogates are treated the same way as data
  h<-FSsync(h)
  deersync_s2[counter,]<-h$sync

  h<-dvcs_s2[[counter]]
  h<-wsyn::cleandat(h,dvcyr,clev=2)$cdat
  h<-FSsync(h)
  dvcsync_s2[counter,]<-h$sync
}

#**
#Use the surrogates of type 2 to judge if "peakiness" is greater for deersync than it is
#for deersync_s2, and same for dvcs
#***

deer_peakiness<-deersync[7]-(deersync[4]+deersync[10])/2
deer_peakiness_s2<-deersync_s2[,7]-(deersync_s2[,4]+deersync_s2[,10])/2
#hist(deer_peakiness_s2)
pval_deer_peakiness<-sum(deer_peakiness_s2>deer_peakiness)/nsurrog

dvc_peakiness<-dvcsync[7]-(dvcsync[5]+dvcsync[10])/2
dvc_peakiness_s2<-dvcsync_s2[,7]-(dvcsync_s2[,5]+dvcsync_s2[,10])/2
#hist(dvc_peakiness_s2)
pval_dvcs_peakiness<-sum(dvc_peakiness_s2>dvc_peakiness)/nsurrog

#***
#Now make the plots
#***

#first do deer, surrogates of type 1
nsurrog_toplot<-min(1000,nsurrog)
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/ComprehensiveFourier_SynchronyIsTimescaleSpecific_01.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(deersync,deersync_s[1:nsurrog_toplot,])
plot(freq_deer,deersync,type="l",lwd=2,ylim=ylimits)
for (counter in 1:nsurrog_toplot)
{
  lines(freq_deer,deersync_s[counter,],type="l",col="grey",lwd=.5)
}
deersync_s_q<-apply(FUN=quantile,X=deersync_s,MARGIN=2,prob=.95)
lines(freq_deer,deersync_s_q,type="l",lwd=2,lty="dashed")
lines(freq_deer,deersync,type="l",lwd=2)
lines(rep(1/3,2),ylimits,lty="dashed")
lines(rep(1/7,2),ylimits,lty="dashed")
mtext("Frequency",side=1,line=2,cex=1)
mtext("Synchrony",side=2,line=2,cex=1)
text(0.5,ylimits[2],"A)",adj=c(1,1),cex=1.5)
dev.off()

#now do deer, surrogates of type 2
nsurrog_toplot<-min(1000,nsurrog)
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/ComprehensiveFourier_SynchronyIsTimescaleSpecific_02.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(deersync,deersync_s[1:nsurrog_toplot,],deersync_s2[1:nsurrog_toplot,])
plot(freq_deer,deersync,type="l",lwd=2,ylim=ylimits)
for (counter in 1:nsurrog_toplot)
{
  lines(freq_deer,deersync_s2[counter,],type="l",col="grey",lwd=.5)
}
deersync_s2_q<-apply(FUN=quantile,X=deersync_s2,MARGIN=2,prob=c(.025,.975))
lines(freq_deer,deersync_s2_q[1,],type="l",lwd=2,lty="dashed")
lines(freq_deer,deersync_s2_q[2,],type="l",lwd=2,lty="dashed")
lines(freq_deer,deersync,type="l",lwd=2)
lines(rep(1/3,2),ylimits,lty="dashed")
lines(rep(1/7,2),ylimits,lty="dashed")
mtext("Frequency",side=1,line=2,cex=1)
mtext("Synchrony",side=2,line=2,cex=1)
text(0.5,ylimits[2],"A)",adj=c(1,1),cex=1.5)
dev.off()

#now do dvcs, surrogates of type 1
nsurrog_toplot<-min(1000,nsurrog)
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/ComprehensiveFourier_SynchronyIsTimescaleSpecific_03.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(dvcsync,dvcsync_s[1:nsurrog_toplot,])
plot(freq_dvcs,dvcsync,type="l",lwd=2,ylim=ylimits)
for (counter in 1:nsurrog_toplot)
{
  lines(freq_dvcs,dvcsync_s[counter,],type="l",col="grey",lwd=.5)
}
dvcsync_s_q<-apply(FUN=quantile,X=dvcsync_s,MARGIN=2,prob=.95)
lines(freq_dvcs,dvcsync_s_q,type="l",lwd=2,lty="dashed")
lines(freq_dvcs,dvcsync,type="l",lwd=2)
lines(rep(1/3,2),ylimits,lty="dashed")
lines(rep(1/7,2),ylimits,lty="dashed")
mtext("Frequency",side=1,line=2,cex=1)
mtext("Synchrony",side=2,line=2,cex=1)
text(0.5,ylimits[2],"B)",adj=c(1,1),cex=1.5)
dev.off()

#now do dvcs, surrogates of type 2
nsurrog_toplot<-min(1000,nsurrog)
tot.wd<-4
tot.ht<-4
ywd<-.65
xht<-.65
gap<-0.1
pan.wd<-tot.wd-ywd-gap
pan.ht<-tot.ht-xht-gap
png("Results/ComprehensiveFourier_SynchronyIsTimescaleSpecific_04.png",res=600,units="in",width = tot.wd,height = tot.ht)
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd)/tot.wd,
          (xht)/tot.ht,
          (xht+pan.ht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.75,0))
ylimits<-range(dvcsync,dvcsync_s[1:nsurrog_toplot,],dvcsync_s2[1:nsurrog_toplot,])
plot(freq_dvcs,dvcsync,type="l",lwd=2,ylim=ylimits)
for (counter in 1:nsurrog_toplot)
{
  lines(freq_dvcs,dvcsync_s2[counter,],type="l",col="grey",lwd=.5)
}
dvcsync_s2_q<-apply(FUN=quantile,X=dvcsync_s2,MARGIN=2,prob=c(.025,.975))
lines(freq_dvcs,dvcsync_s2_q[1,],type="l",lwd=2,lty="dashed")
lines(freq_dvcs,dvcsync_s2_q[2,],type="l",lwd=2,lty="dashed")
lines(freq_dvcs,dvcsync,type="l",lwd=2)
lines(rep(1/3,2),ylimits,lty="dashed")
lines(rep(1/7,2),ylimits,lty="dashed")
mtext("Frequency",side=1,line=2,cex=1)
mtext("Synchrony",side=2,line=2,cex=1)
text(0.5,ylimits[2],"B)",adj=c(1,1),cex=1.5)
dev.off()
