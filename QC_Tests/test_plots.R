#Generate figures using wsyn base code
setwd("C:/Users/Tom/Documents/GitRepos/mbsync/AndersDeerShort")
library(reshape2)
library(tools)
library(fields)
library(shape)
library(devtools)
library(wsyn)

minyear<-1981
maxyear<-2016
dens.flag<-"dnr"
nsurrogs<-1000

source("Code/Deer_DataGeneration.R")
source("Code/Climate_DataGeneration.R")
source("Code/DeerCoherence_wsyn.R")

pdf("Results/test.plots.pdf")
#Figure 2
abun.dt<-wsyn::cleandat(cty.list$Abun,clev=5,times=1981:2016)$cdat
abun.wmf<-wsyn::wmf(abun.dt,times=1981:2016)
abun.wpmf<-wsyn::wpmf(abun.dt,times=1981:2016,sigmethod = "quick")
plotmag(abun.wmf)
title("Fig2A",adj=0)
plotmag(abun.wpmf,sigthresh=0.999,zlims = c(0,1))
title("Fig2B",adj=0)
plotmag(predsync(wlm_abun))
contour(x=abun.wmf$times,y=log2(abun.wmf$timescales),z=Mod(abun.wmf$values),add=T,
        drawlabels=T,lwd=2,frame=F)
title("Fig2C",adj=0)

#Figure 3
dvc.dt<-wsyn::cleandat(cty.list$Crashes[,!is.na(colSums(cty.list$Crashes))],clev=5,times=1987:2016)$cdat
dvc.wmf<-wsyn::wmf(dvc.dt,times=1987:2016)
dvc.wpmf<-wsyn::wpmf(dvc.dt,times=1987:2016,sigmethod = "quick")

plotmag(dvc.wmf)
title("Fig3A",adj=0)
plotmag(dvc.wpmf,sigthresh = 0.999,zlims = c(0,1))
title("Fig3B",adj=0)
plotmag(predsync(wlm_dvc))
contour(x=dvc.wmf$times,y=log2(dvc.wmf$timescales),z=Mod(dvc.wmf$values),add=T,
        drawlabels=T,lwd=2,frame=F)
title("Fig3C",adj=0)

#Figure S3
par(mfrow=c(4,2),mar=c(3,3,1.2,0.2),mgp=c(1.5,0.5,0))
plotphase(clim.res$WinterMEI.Abun)
mtext("FigS3A",adj=0,font=2)
plotphase(clim.res$SummerMEI.Abun)
mtext("FigS3B",adj=0,font=2)
plotphase(clim.res$WinterPDO.Abun)
mtext("FigS3C",adj=0,font=2)
plotphase(weath.res$Snwd.Abun)
mtext("FigS3D",adj=0,font=2)
plotphase(dvc.res)
mtext("FigS3E",adj=0,font=2)
plotphase(adjdvc.res)
mtext("FigS3F",adj=0,font=2)
plotphase(hunter.res)
mtext("FigS3G",adj=0,font=2)

#Figure S4
par(mfrow=c(4,2),mar=c(3,3,1.2,0.2),mgp=c(1.5,0.5,0))
plotrank(clim.res$WinterMEI.Abun)
mtext("FigS4A",adj=0,font=2)
plotrank(dvc.res)
mtext("FigS4B",adj=0,font=2)
plotrank(clim.res$SummerMEI.Abun)
mtext("FigS4A",adj=0,font=2)
plotrank(adjdvc.res)
mtext("FigS4B",adj=0,font=2)
plotrank(clim.res$WinterPDO.Abun)
mtext("FigS4A",adj=0,font=2)
plotrank(hunter.res)
mtext("FigS4B",adj=0,font=2)
plotrank(weath.res$Snwd.Abun)
mtext("FigS4A",adj=0,font=2)

#Figure S5-hunter plots
hunters.tmp<-cty.list$Hunters[!(row.names(cty.list$Hunters)%in%cwd),12:36]
hunters.tmp<-hunters.tmp[!is.na(rowSums(hunters.tmp)),]
abun.tmp<-cty.list$Abun[row.names(cty.list$Abun)%in%c(row.names(hunters.tmp)),12:36]
hunter.dt<-wsyn::cleandat(hunters.tmp,clev=5,times=1992:2016)$cdat
abun.tmp.dt<-wsyn::cleandat(abun.tmp,clev=5,times=1992:2016)$cdat

hunter.wmf<-wsyn::wmf(hunter.dt,times=1992:2016)
hunter.wpmf<-wsyn::wpmf(hunter.dt,times=1992:2016,sigmethod = "quick")
abun.tmp.wmf<-wsyn::wmf(abun.tmp.dt,times=1992:2016)
par(mfrow=c(3,1),mar=c(3,3,1.2,0.2),mgp=c(1.5,0.5,0))
plotmag(hunter.wmf)
mtext("FigS5A",adj=0)
plotmag(hunter.wpmf,zlims = c(0,1),sigthresh=0.999)
mtext("FigS5B",adj=0)
plotmag(predsync(wlm_hunters))
mtext("FigS5C",adj=0)
contour(x=abun.tmp.wmf$times,y=log2(abun.tmp.wmf$timescales),z=Mod(abun.tmp.wmf$values),add=T,
        drawlabels=T,lwd=2,frame=F)

#Figure S6
usda.abun.dt<-cleandat(usda.list$Abun,clev=5,times=minyear:maxyear)$cdat
usda.dvc.dt<-cleandat(usda.list$Crashes[,!is.na(colSums(usda.list$Crashes))],clev=5,times=1987:maxyear)$cdat

abun.wpmf.usda<-wsyn::wpmf(usda.abun.dt,times = minyear:maxyear,sigmethod = "quick")
abun.wmf.usda<-wsyn::wmf(usda.abun.dt,times = minyear:maxyear)
dvc.wpmf.usda<-wsyn::wpmf(usda.dvc.dt,times = 1987:maxyear,sigmethod = "quick")
dvc.wmf.usda<-wsyn::wmf(usda.dvc.dt,times = 1987:maxyear)

par(mfrow=c(3,2),mar=c(2.5,3,0.5,4),mgp=c(1.5,0.5,0))
plotmag(abun.wmf.usda)
mtext("FigS6A",adj=0.05,line=-1.2,side=3,font=2)
plotmag(dvc.wmf.usda)
mtext("FigS6B",adj=0.05,line=-1.2,side=3,font=2)
plotmag(abun.wpmf.usda,sigthresh = 0.999,zlims = c(0,1))
mtext("FigS6C",adj=0.05,line=-1.2,side=3,font=2)
plotmag(dvc.wpmf.usda,sigthresh = 0.999,zlims = c(0,1))
mtext("FigS6D",adj=0.05,line=-1.2,side=3,font=2)
plotmag(predsync(usda.wlm_abun))
mtext("FigS6E",adj=0.05,line=-1.2,side=3,font=2)
contour(x=abun.wmf.usda$times,y=log2(abun.wmf.usda$timescales),z=Mod(abun.wmf.usda$values),add=T,
        drawlabels=T,lwd=2,frame=F)
plotmag(predsync(usda.wlm_dvc))
mtext("FigS6F",adj=0.05,line=-1.2,side=3,font=2)
contour(x=dvc.wmf.usda$times,y=log2(dvc.wmf.usda$timescales),z=Mod(dvc.wmf.usda$values),add=T,
        drawlabels=T,lwd=2,frame=F)

#Figure S7
winter.clim.dt<-lapply(winter.clim,function(x){x<-cleandat(x[!(is.na(rowMeans(x))),],clev=5,times=minyear:maxyear)$cdat;x})
winter.clim.wmf<-lapply(winter.clim.dt,function(x){x<-wsyn::wmf(x,times=1981:2016);x})
winter.clim.wpmf<-lapply(winter.clim.dt,function(x){x<-wsyn::wpmf(x,times=1981:2016,sigmethod="quick");x})

par(mfrow=c(4,2),mar=c(2.5,3,0.5,4),mgp=c(1.5,0.5,0))
plotmag(winter.clim.wmf$Tmin,colorbar=T)
mtext("FigS7A)",adj=0.05,line=-1.2,side=3,font=2)
plotmag(winter.clim.wpmf$Tmin,colorbar=T,sigthresh = 0.999,zlims = c(0,1))
mtext("FigS7B)",adj=0.05,line=-1.2,side=3,font=2)
plotmag(winter.clim.wmf$Tmax,colorbar=T)
mtext("FigS7C)",adj=0.05,line=-1.2,side=3,font=2)
plotmag(winter.clim.wpmf$Tmax,colorbar=T,sigthresh = 0.999,zlims = c(0,1))
mtext("FigS7D)",adj=0.05,line=-1.2,side=3,font=2)
plotmag(winter.clim.wmf$Snwd,colorbar=T)
mtext("FigS7E)",adj=0.05,line=-1.2,side=3,font=2)
plotmag(winter.clim.wpmf$Snwd,colorbar=T,sigthresh = 0.999,zlims = c(0,1))
mtext("FigS7F)",adj=0.05,line=-1.2,side=3,font=2)
plotmag(winter.clim.wmf$Prcp,colorbar=T)
mtext("FigS7G)",adj=0.05,font=2,line=-1.2,side=3)
plotmag(winter.clim.wpmf$Prcp,colorbar=T,sigthresh = 0.999,zlims = c(0,1))
mtext("FigS7H)",adj=0.05,font=2,line=-1.2,side=3)

#Figure S8
climindex.dt<-lapply(climindex,function(x){x<-cleandat(x[1,],clev=5,times=minyear:maxyear)$cdat;x})
climindex.wt<-lapply(climindex.dt,function(x){x<-wsyn::wt(x,times=1981:2016);x})

par(mfrow=c(3,2),mar=c(2.5,3,0,4),mgp=c(1.5,0.5,0))
plotmag(climindex.wt$WinterPDO,colorbar=T)
mtext("FigS8A)",adj=0.05,line=-1.2,side=3,font=2)
plotmag(climindex.wt$SummerPDO,colorbar=T)
mtext("FigS8B)",adj=0.05,line=-1.2,side=3,font=2)
plotmag(climindex.wt$WinterNAO,colorbar=T)
mtext("FigS8C)",adj=0.05,line=-1.2,side=3,font=2)
plotmag(climindex.wt$SummerNAO,colorbar=T)
mtext("FigS8D)",adj=0.05,font=2,line=-1.2,side=3)
plotmag(climindex.wt$WinterMEI,colorbar=T)
mtext("FigS8E)",adj=0.05,line=-1.2,side=3,font=2)
plotmag(climindex.wt$SummerMEI,colorbar=T)
mtext("FigS8F)",adj=0.05,line=-1.2,side=3,font=2)
dev.off()
