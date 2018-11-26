#Make figures for the manuscript
source("Functions/Fn_syncexpplot.R")
source("Functions/deer_plotting_functions.R")
source("Functions/Fn_phaseplot.R")
source("Functions/Fn_rankplot.R")

#Make Fig. 1
source("Code/PedagogFig.R")

#clean data for Figs 2 and 3
abun.dt<-cleandat(cty.list$Abun,clev=5,times=minyear:maxyear)$cdat
dvc.dt<-cleandat(cty.list$Crashes[,!is.na(colSums(cty.list$Crashes))],clev=5,times=1987:maxyear)$cdat

#load model results needed
wlm_abun<-readRDS(file="Results/wlm_abun.rds")
wlm_dvc<-readRDS(file="Results/wlm_dvc.rds")
wlm_hunters<-readRDS(file="Results/wlm_hunters.rds")

# Set up dimensions of wavelet mean field and phasor mean fields for Figs 2 and 3
tot.wd<-4.75
xht<-0.75   #height of x axis label region
ywd<-0.5    #width of y axis label region
zwd<-0.5    #width of z-axes label region
gap<-.2   #small gap 
pan.wd<-(tot.wd-ywd-gap-zwd) #large panel width parameter
pan.ht<-0.75*pan.wd+zwd #big ones are square
tot.ht<-3*pan.ht+xht+3*gap

abun.wpmf<-wsyn::wpmf(abun.dt,times = minyear:maxyear,sigmethod = "quick")
abun.wmf<-wsyn::wmf(abun.dt,times = minyear:maxyear)
dvc.wpmf<-wsyn::wpmf(dvc.dt,times = 1987:maxyear,sigmethod = "quick")
dvc.wmf<-wsyn::wmf(dvc.dt,times = 1987:maxyear)

#Make Figure 2
png("Results/Fig2.png",res=600,height=tot.ht,width=tot.wd,unit="in")
par(mfrow=c(3,1),mgp=c(3.5,1.25,0),mai=c(1,0.75,0.2,0))
deer_wmfplot(abun.wmf,xlab="",ylab="Timescale (yrs)",cex.lab=3,cex.axis=2.25,las=1)
abline(h=c(log2(3),log2(7)),lty=2)
mtext(text = "A)",font=2,side = 3,adj =0.05,line=-1,cex=2)
deer_wpmfplot(abun.wpmf,sigthresh = 0.001,xlab="",ylab="Timescale (yrs)",cex.lab=3,cex.axis=2.25,las=1)
par(new=T)
q<-stats::quantile(abun.wpmf$signif[[2]],0.999)
contour(x=abun.wpmf$times,y=log2(abun.wpmf$timescales),z=Mod(abun.wpmf$values),levels=q,drawlabels=F,lwd=2,
        xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5),las = 1,frame=F)
abline(h=c(log2(3),log2(7)),lty=2)
mtext(text = "B)",font=2,side = 3,adj =0.05,line=-1,cex=2)
syncexpplot(resp.wmf=abun.wmf$values,exp.sync = predsync(wlm_abun)[[3]],1981:2016,
            wlm_abun$timescales,xlab="Year",smallplot=c(0.95,0.99,0.05,0.95),ylab="Timescale (yrs)",
            cex.lab=3,cex.axis=2.25)
mtext(text = "C)",font=2,side = 3,adj =0.05,line=-1,cex=2)
dev.off()

#Make Figure 3
png("Results/Fig3.png",res=600,height=tot.ht,width=tot.wd,unit="in")
par(mfrow=c(3,1),mgp=c(3.5,1.25,0),mai=c(1,0.75,0.2,0))
deer_wmfplot(dvc.wmf,xlab="",ylab="Timescale (yrs)",cex.lab=3,cex.axis=2.25,las=1)
mtext(text = "A)",font=2,side = 3,adj =0.05,line=-1,cex=2)
abline(h=c(log2(3),log2(7)),lty=2)
deer_wpmfplot(dvc.wpmf,sigthresh = 0.001,xlab="",ylab="Timescale (yrs)",cex.lab=3,cex.axis=2.25,las=1)
par(new=T)
q<-stats::quantile(dvc.wpmf$signif[[2]],0.999)
contour(x=dvc.wpmf$times,y=log2(dvc.wpmf$timescales),z=Mod(dvc.wpmf$values),levels=q,drawlabels=F,lwd=2,
        xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5),las = 1,frame=F)
abline(h=c(log2(3),log2(7)),lty=2)
mtext(text = "B)",font=2,side = 3,adj =0.05,line=-1,cex=2)
syncexpplot(resp.wmf=dvc.wmf$values,exp.sync = predsync(wlm_dvc)[[3]],1987:2016,
            wlm_dvc$timescales,xlab="Year",smallplot=c(0.95,0.99,0.05,0.95),ylab="Timescale (yrs)",
            cex.lab=3,cex.axis=2.25)
mtext(text = "C)",font=2,side = 3,adj =0.05,line=-1,cex=2)
dev.off()

#make plots of phase by timescale for all significant pairs of variables (Fig S4)
png("Results/FigS4.png",res=600,height=4800,width=3200)
TabS1_results<-readRDS(file = "Results/TableS1.rds")
TabS1_results<-cbind(TabS1_results[,1:3],apply(TabS1_results[,4:7],2,function(x){round(x,4)}))
wmeiabun_phase47<-TabS1_results$MeanPhase[TabS1_results$Predictor=="Winter MEI" &  TabS1_results$Timescale=="4-7"]
wmeiabun_phase37<-TabS1_results$MeanPhase[TabS1_results$Predictor=="Winter MEI" &  TabS1_results$Timescale=="3-7"]
wpdoabun_phase47<-TabS1_results$MeanPhase[TabS1_results$Predictor=="Winter PDO" &  TabS1_results$Timescale=="4-7"]
wpdoabun_phase37<-TabS1_results$MeanPhase[TabS1_results$Predictor=="Winter PDO" &  TabS1_results$Timescale=="3-7"]
snowabun_phase34<-TabS1_results$MeanPhase[TabS1_results$Predictor=="Snow Depth" & TabS1_results$Timescale=="3-4"]
snowabun_phase37<-TabS1_results$MeanPhase[TabS1_results$Predictor=="Snow Depth" & TabS1_results$Timescale=="3-7"]
hunterphase2_2.5<-TabS1_results$MeanPhase[TabS1_results$Predictor=="Hunters" & TabS1_results$Timescale=="2-2.5"]

par(mfrow=c(4,2),mar=c(2.5,4,2.5,0.5),mgp=c(1.5,0.5,0),cex.lab=1.5,las=1)
phaseplot(clim.coher$WinterMEI.Abun,abun.wmf$timescales,showphase = F,type="pi",tsrange=c(3,7),xlab="",ylab="Phase")
rect(xleft=4,ybottom = -pi,xright = 7,ytop = pi,col=rgb(1,0,0,0.5),density=25,lwd=2,angle=-45)
mtext(line=0,adj=1,cex=0.75,bquote(paste("3-7, ",bar(theta)==.(round(wmeiabun_phase37,3)))))
mtext(line=1,adj=1,cex=0.75,bquote(paste("4-7, ",bar(theta)==.(round(wmeiabun_phase47,3)))))
mtext("A)",font=2,adj=0.05)
phaseplot(clim.coher$SummerMEI.Abun,abun.wmf$timescales,type="pi",tsrange=c(3,7),xlab="",ylab="")
mtext("B)",font=2,adj=0.05)
phaseplot(clim.coher$WinterPDO.Abun,abun.wmf$timescales,showphase = F,type="pi",tsrange=c(3,7),xlab="",ylab="Phase")
rect(xleft=4,ybottom = -pi,xright = 7,ytop = pi,col=rgb(1,0,0,0.5),density=25,lwd=2,angle=-45)
mtext(line=0,adj=1,cex=0.75,bquote(paste("3-7, ",bar(theta)==.(round(wpdoabun_phase37,3)))))
mtext(line=1,adj=1,cex=0.75,bquote(paste("4-7, ",bar(theta)==.(round(wpdoabun_phase47,3)))))
mtext("C)",font=2,adj=0.05)
phaseplot(snow.coher,abun.wmf$timescales,showphase=F,type="pi",tsrange=c(3,7),ylab="",xlab="")
mtext("D)",font=2,adj=0.05)
rect(xleft=3,ybottom = -pi,xright = 4,ytop = pi,col=rgb(1,0,0,0.5),density=25,lwd=2,angle=-45)
mtext(line=0,adj=1,cex=0.75,bquote(paste("3-7, ",bar(theta)==.(round(snowabun_phase37,3)))))
mtext(line=1,adj=1,cex=0.75,bquote(paste("3-4, ",bar(theta)==.(round(snowabun_phase34,3)))))
phaseplot(dvc.coher,dvc.wmf$timescales,type="pi",tsrange=c(3,7),xlab="",ylab="Phase")
mtext("E)",font=2,adj=0.05)
phaseplot(adjdvc.coher,adjdvc.res$timescales,type="pi",tsrange=c(3,7),xlab="Timescale (yrs)",ylab="")
mtext("F)",font=2,adj=0.05)
phaseplot(get_coher(hunter.res),hunter.res$timescales,showphase=F,type="pi",tsrange=c(3,7),xlab="Timescale (yrs)",ylab="Phase",colrect=F)
mtext("G)",font=2,adj=0.05)
rect(xleft=2,ybottom = -pi,xright = 2.5,ytop = pi,col=rgb(1,0,0,0.5),density=25,lwd=2,angle=-45)
mtext(line=0,adj=1,cex=0.75,bquote(paste("2-2.5, ",bar(theta)==.(round(hunterphase2_2.5,3)))))
dev.off()

#Make plot of ranks by timescale for all significant pairs of variables (Fig S6)
png("Results/FigS6.png",res=600,height=3000,width=4800)
par(mfrow=c(1,2),mar=c(3.5,4,1.5,0),mgp=c(2.5,0.5,0))
plot(abun.wmf$timescales,clim.ranks$WinterMEI.Abun$coher,xaxt="n",las=1,type="l",lwd=2,ylab="Rank",xlab="Timescale (yrs)")
axis(1,at = seq(2,12,1),labels=c(2,"","4","",6,"",8,"",10,"",12))
axis(3,at = seq(2,12,1),labels=c(2,"","4","",6,"",8,"",10,"",12))
abline(h=0.95,col="red",lty=2)
lines(abun.wmf$timescales,clim.ranks$SummerMEI.Abun$coher,lwd=2,col="green")
lines(abun.wmf$timescales,clim.ranks$WinterPDO.Abunks$coher,lty=1,lwd=2,col="orange")
lines(abun.wmf$timescales,snow.ranks$coher,lty=1,lwd=2,col="blue")
legend("topright",lty=c(1,1),c("Winter MEI","Winter PDO","Summer MEI","Snow Depth"),cex=0.75,bty = "n",col=c("black","orange","green","blue"))
mtext("A)",font=2,adj=-0.1)

plot(dvc.res$timescales,dvc.res$ranks$coher,type="l",lwd=2,ylab="",xlab="Timescale (yrs)",xaxt="n",las=1)
lines(adjdvc.res$timescales,adjdvc.res$ranks$coher,lty=1,lwd=2,col="purple")
axis(1,at = seq(2,10,1),labels=c(2,"",4,"",6,"",8,"",10))
axis(3,at = seq(2,10,1),labels=c(2,"",4,"",6,"",8,"",10))
abline(h=0.95,col="red",lty=2)
lines(hunter.res$timescales,hunter.res$ranks$coher,lty=1,lwd=2,col="darkgray")
legend("topright",c("Abun-DVCs","Abun-Adj. DVCs","Hunters-Abun"),lty=c(1,1,1),
       col=c("black","purple","darkgray"),bty="n",cex=0.75)
mtext("B)",font=2,adj=-0.1)
dev.off()

#Make wavelet mean field, wavelet phasor meanfield and predicted synchrony plot for hunters (Fig S7)
png("Results/FigS7.png",res=600,height=tot.ht,width=tot.wd,unit="in")
hunters.tmp<-cty.list$Hunters[,12:dim(cty.list$Hunters)[2]]
hunters.tmp<-hunters.tmp[(!row.names(hunters.tmp)%in%cwd) & !is.na(rowMeans(hunters.tmp)),]
abun.tmp<-cty.list$Abun[row.names(cty.list$Abun)%in%row.names(hunters.tmp),12:36]
hunters.dt<-cleandat(hunters.tmp,clev=5,times=1992:2016)$cdat
abun.dt<-cleandat(abun.tmp,clev=5,times=1992:2016)$cdat

hunters.wpmf<-wsyn::wpmf(hunters.dt,times=1992:2016,sigmethod = "quick")
hunters.wmf<-wsyn::wmf(hunters.dt,times=1992:2016)
abun.wmf.short<-wsyn::wmf(abun.dt,times=1992:2016)

par(mfrow=c(3,1),mgp=c(3.5,1.25,0),mai=c(1,0.75,0.2,0))
deer_wmfplot(hunters.wmf,colorbar=T,ylab="",xlab="",las=1,cex.lab=3,cex.axis=2.25)
abline(h=c(log2(3),log2(7)),lty=2)
mtext("A)",adj=0.05,line=-1.2,side=3,font=2,cex=2)
deer_wpmfplot(hunters.wpmf,colorbar=T,zlims = c(0,1),xlab="",sigthresh=0.999,ylab="Timescale (yrs)",las=1,cex.lab=3,cex.axis=2.25)
mtext("B)",adj=0.05,line=-1.2,side=3,font=2,cex=2)
abline(h=c(log2(3),log2(7)),lty=2)
syncexpplot(resp.wmf=abun.wmf.short$values,exp.sync = predsync(wlm_hunters)[[3]],times=1992:2016,
            wlm_hunters$timescales,ylab = "",xlab="Year",cex.lab=3,cex.axis=2.25)
mtext("C)",adj=0.05,line=-1.2,side=3,font=2,cex=2)
dev.off()

#Make winter weather wavelet mean field and wavelet phasor mean field plots
png("Results/FigS9.png",res=600,height=9,width=7,unit="in")
winter.clim.dt<-lapply(winter.clim,function(x){x<-cleandat(x[!(is.na(rowMeans(x))),],clev=5,times=minyear:maxyear)$cdat;x})
par(mfrow=c(4,2),mar=c(2.5,3,0,4),mgp=c(1.5,0.5,0))
tmin.wmf<-wsyn::wmf(winter.clim.dt$Tmin,times=1981:2016)
deer_wmfplot(tmin.wmf,colorbar=T,xlab="",ylab="Timescale (yrs)",las=1)
mtext("A)",adj=0.05,line=-1.2,side=3,font=2)
tmin.wpmf<-wsyn::wpmf(winter.clim.dt$Tmin,times=1981:2016,sigmethod="quick")
deer_wpmfplot(tmin.wpmf,colorbar=T,ylab="Timescale (yrs)",xlab="",sigthresh = 0.999,las=1)
mtext("B)",adj=0.05,line=-1.2,side=3,font=2)
tmax.wmf<-wsyn::wmf(winter.clim.dt$Tmax,times=1981:2016)
deer_wmfplot(tmax.wmf,colorbar=T,xlab="",ylab="Timescale (yrs)",las=1)
mtext("C)",adj=0.05,line=-1.2,side=3,font=2)
tmax.wpmf<-wsyn::wpmf(winter.clim.dt$Tmax,times=1981:2016,sigmethod="quick")
deer_wpmfplot(tmax.wpmf,colorbar=T,ylab="Timescale (yrs)",xlab="",sigthresh = 0.999,las=1)
mtext("D)",adj=0.05,line=-1.2,side=3,font=2)
snwd.wmf<-wsyn::wmf(winter.clim.dt$Snwd,times=1981:2016)
deer_wmfplot(snwd.wmf,colorbar=T,xlab="",ylab="Timescale (yrs)",las=1)
mtext("E)",adj=0.05,line=-1.2,side=3,font=2)
snwd.wpmf<-wsyn::wpmf(winter.clim.dt$Snwd,times=1981:2016,sigmethod="quick")
deer_wpmfplot(snwd.wpmf,colorbar=T,ylab="Timescale (yrs)",xlab="",sigthresh = 0.999,las=1)
mtext("F)",adj=0.05,line=-1.2,side=3,font=2)
prcp.wmf<-wsyn::wmf(winter.clim.dt$Prcp,times=1981:2016)
deer_wmfplot(prcp.wmf,colorbar=T,xlab="Year",ylab="Timescale (yrs)",las=1)
mtext("G)",adj=0.05,font=2,line=-1.2,side=3)
prcp.wpmf<-wsyn::wpmf(winter.clim.dt$Prcp,times=1981:2016,sigmethod="quick")
deer_wpmfplot(prcp.wpmf,colorbar=T,ylab="Timescale (yrs)",xlab="Year",sigthresh = 0.999,las=1)
mtext("H)",adj=0.05,font=2,line=-1.2,side=3)
dev.off()

#Make Figure S8: plots of wavelet mean field, wavelet phasor mean field, and predicted synchrony of 
#deer abundance and DVCs using USDA district data
usda.abun.dt<-cleandat(usda.list$Abun,clev=5,times=minyear:maxyear)$cdat
usda.dvc.dt<-cleandat(usda.list$Crashes[,!is.na(colSums(usda.list$Crashes))],clev=5,times=1987:maxyear)$cdat

#generate mean fields
abun.wpmf.usda<-wsyn::wpmf(usda.abun.dt,times = minyear:maxyear,sigmethod = "quick")
abun.wmf.usda<-wsyn::wmf(usda.abun.dt,times = minyear:maxyear)
dvc.wpmf.usda<-wsyn::wpmf(usda.dvc.dt,times = 1987:maxyear,sigmethod = "quick")
dvc.wmf.usda<-wsyn::wmf(usda.dvc.dt,times = 1987:maxyear)

#make figure
png("Results/FigS8.png",res=600,height=9,width=7,unit="in")
par(mfrow=c(3,2),mar=c(2.5,3,0,4),mgp=c(1.5,0.5,0))
deer_wmfplot(abun.wmf.usda,xlab="",ylab="Timescale (yrs)",las=1)
mtext("A)",adj=0.05,line=-1.2,side=3,font=2)
deer_wmfplot(dvc.wmf.usda,xlab="",ylab="Timescale (yrs)",las=1)
mtext("B)",adj=0.05,line=-1.2,side=3,font=2)
deer_wpmfplot(abun.wpmf.usda,sigthresh = 0.999,xlab="",ylab="Timescale (yrs)",las=1)
par(new=T)
q<-stats::quantile(abun.wpmf.usda$signif[[2]],0.999)
contour(x=abun.wpmf.usda$times,y=log2(abun.wpmf.usda$timescales),z=Mod(abun.wpmf.usda$values),levels=q,drawlabels=F,lwd=2,
        xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5),las = 1,frame=F)
mtext("C)",adj=0.05,line=-1.2,side=3,font=2)
deer_wpmfplot(dvc.wpmf.usda,sigthresh = 0.999,xlab="",ylab="Timescale (yrs)",las=1)
par(new=T)
q<-stats::quantile(dvc.wpmf.usda$signif[[2]],0.999)
contour(x=dvc.wpmf.usda$times,y=log2(dvc.wpmf.usda$timescales),z=Mod(dvc.wpmf.usda$values),levels=q,drawlabels=F,lwd=2,
        xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5),las = 1,frame=F)
mtext("D)",adj=0.05,line=-1.2,side=3,font=2)
syncexpplot(resp.wmf=abun.wmf.usda$values,exp.sync = predsync(usda.wlm_abun)[[3]],1981:2016,
            usda.wlm_abun$timescales,xlab="Year",smallplot=c(0.95,0.99,0.05,0.95),ylab="Timescales (yrs)")
mtext("E)",adj=0.05,line=-1.2,side=3,font=2)
syncexpplot(resp.wmf=dvc.wmf.usda$values,exp.sync = predsync(usda.wlm_dvc)[[3]],1987:2016,
            usda.wlm_dvc$timescales,xlab="Year",smallplot=c(0.95,0.99,0.05,0.95),ylab="Timescales (yrs)")
mtext("F)",adj=0.05,line=-1.2,side=3,font=2)
dev.off()

#Make raw data plots (Fig S1)
png("Results/FigS1.png",res=600,height=3600,width=3600)
par(mfrow=c(3,2),mar=c(2.5,5,1.2,0),mgp=c(3,0.5,0))
plot(1981:2016,rep(NA,36),las=1,ylim=c(min(cty.list$Abun),max(cty.list$Abun)),xlab="",ylab="Abundance")
for(i in 1:nrow(cty.list$Abun)){
  lines(1981:2016,cty.list$Abun[i,],col=rainbow(71)[i])
}
mtext(text = "A)",side=3,adj=0.05,font=2)
plot(1987:2016,rep(NA,30),las=1,ylim=c(min(cty.list$Crashes,na.rm=T),max(cty.list$Crashes,na.rm=T)),xlim=c(1981,2016),xlab="Year",ylab="DVCs")
for(i in 1:nrow(cty.list$Crashes)){
  lines(1987:2016,na.omit(cty.list$Crashes[i,]),col=rainbow(71)[i])
}
mtext(text = "B)",side=3,adj=0.05,font=2)

winter.clim.tmp<-lapply(winter.clim,function(x){x[(!is.na(rowMeans(x))),]})
plot(1981:2016,rep(NA,36),las=1,ylim=c(min(winter.clim.tmp$Snwd),max(winter.clim.tmp$Snwd)),xlab="",ylab="Winter Snow Depth (cm)")
for(i in 1:nrow(winter.clim.tmp$Snwd)){
  lines(1981:2016,winter.clim.tmp$Snwd[i,],col=rainbow(nrow(winter.clim.tmp$Snwd))[i])
}
mtext(text = "C)",side=3,adj=0.05,font=2)

plot(mei$YEAR[mei$YEAR>1980],win.mei.mat[1,],las=1,col="blue",type="l",xlab="Year",ylab="Index Value",lwd=2)
lines(mean.pdo$Year[mean.pdo$Year>1980],win.pdo.mat[1,],col="green",lwd=2)
legend("topright",bty="n",c("Winter MEI","Winter PDO"),col=c("blue","green"),lty=c(1,1),cex=0.75)
mtext(text = "D)",side=3,adj=0.05,font=2)
mtext(text="Year",side=1,line=1.5,cex=0.75)

cty.list.tmp<-cty.list$Hunters[,12:dim(cty.list$Hunters)[2]]
cty.list.tmp<-cty.list.tmp[(!row.names(cty.list.tmp)%in%cwd) & !is.na(rowMeans(cty.list.tmp)),]
plot(1992:2016,rep(NA,25),ylim=c(min(cty.list.tmp,na.rm=T),max(cty.list.tmp,na.rm=T)),xlab="Year",ylab="Hunters",xlim=c(1981,2016),las=1)
for(i in 1:nrow(cty.list.tmp)){
  lines(1992:2016,na.omit(cty.list.tmp[i,]),col=rainbow(53)[i])
}
mtext(text = "E)",side=3,adj=0.05,font=2)
mtext(text="Year",side=1,line=1.5,cex=0.75)
dev.off()

#make Lawrence's plot showing how synchronous impacts statewide fluctuations
abunsurr<-read.csv("Data/abunsurrsum.csv")
dvcsurr<-read.csv("Data/dvcsurrsum.csv")

#Make plots showing statewide magnitude of deer and DVC fluctuatoins (Fig 4)
#***Set up plotting dimensions, units are inches
tot.wd<-3.5
xht<-.5   #height of x axis label region
ywd<-1    #width of y axis label region
gap<-.1   #small gap 
pan.wd.big<-(tot.wd-ywd-gap) #large panel width parameter
pan.ht.big<-pan.wd.big #big ones are square
pan.wd.small<-pan.wd.big #small panel width param
pan.ht.small<-0.33*pan.ht.big #small panel height param
tot.ht<-2*pan.ht.big+2*pan.ht.small+2*xht+4*gap

png("Results/Fig4.png",res=600,units="in",width = tot.wd,height = tot.ht)

#Deer- little panel
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd.small)/tot.wd,
          (2*pan.ht.big+pan.ht.small+2*xht+2*gap)/tot.ht,
          (2*pan.ht.big+2*pan.ht.small+2*xht+2*gap)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.5,0))
colors<-gray.colors(25)[sample(1:25, 999, replace=TRUE)]
plot(1981:2016,rep(NA,36),ylim=c(min(abunsurr),max(abunsurr)),ylab="",xlab="",las=1,axes=F)
for(i in 1:nrow(abunsurr)){
  lines(1981:2016,abunsurr[i,],col=colors[i])
}
lines(1981:2016,apply(cty.list$Abun,2,sum),lwd=2)
axis(2,labels= format(seq(600000,1300000,200000),scientific=T),at = seq(600000,1300000,200000),las=1,tck=-.05,cex.axis=0.75)
axis(1,labels=c(rep("",8)),at = seq(1980,2016,5),tck=-0.05,cex.axis=0.75)
mtext("Deer",side=2,line=3.5,cex=0.75)
mtext("A)",font=2,side=3,line=0,adj=0.05)
box()
#Deer- big panel
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd.big)/tot.wd,
          (2*xht+pan.ht.big+gap+pan.ht.small)/tot.ht,
          (2*xht+2*pan.ht.big+gap+pan.ht.small)/tot.ht),
    mai=c(0,0,0,0),mgp=c(3,0.5,0),new=T)
plot(1981:2016,rep(NA,36),ylab="",xlab="",ylim=c(-180000,200000),las=1)
for(i in 1:nrow(dvcsurr)){
  lines(1981:2016,abunsurr[i,]-apply(abunsurr,2,mean),col=colors[i])
}
lines(1981:2016,(apply(cty.list$Abun,2,sum)-apply(abunsurr,2,mean)),lwd=3,type="b",pch=19)
Arrows(x0 = 1999,y0 = 0,y1=140000,x1=1999,arr.type = "triangle",arr.adj=1,arr.length=0.2,lwd=2,col="red")
Arrows(x0 = 1997,y0 = 0,y1=-150000,x1=1997,arr.type = "triangle",arr.adj=1,arr.length=0.2,lwd=2,col="red")
text(x=2001,y=160000,labels="159054 deer",font=2,cex=0.75,adj=0.05,col="red")
text(x=1999,y=-174000,labels="174339 deer",font=2,cex=0.75,adj=0.05,col="red")
mtext("Departure from Surrogate Mean",side=2,line=3.5,cex=0.9)

#DVC little panel
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd.small)/tot.wd,
          (pan.ht.big+xht+gap)/tot.ht,
          (pan.ht.big+pan.ht.small+xht+gap)/tot.ht),
    mai=c(0,0,0,0),new=T)
plot(1987:2016,rep(NA,30),ylim=c(min(dvcsurr),max(dvcsurr)),ylab="DVCs",xlab="",las=1,axes=F,cex.axis=0.75)
for(i in 1:nrow(dvcsurr)){
  lines(1987:2016,dvcsurr[i,],col=colors[i])
}
lines(1987:2016,apply(cty.list$Crashes,2,sum,na.rm=T)[-c(1:6)],lwd=2)
axis(1,labels=c(rep("",7)),at = seq(1985,2016,5),tck=-0.05,cex.axis=0.75)
axis(2,labels= format(seq(16000,24000,2000),scientific=T),at = seq(16000,24000,2000),las=1,tck=-.05,cex.axis=0.75)
mtext("DVCs",side=2,line=3.5,cex=0.75)
mtext("B)",font=2,side=3,line=0,adj=0.05)
box()
#DVC big panel
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd.big)/tot.wd,
          xht/tot.ht,
          (pan.ht.big+xht)/tot.ht),
    mai=c(0,0,0,0),mgp=c(1,0.5,0),new=T)
plot(1987:2016,rep(NA,30),type="b",xlab="Year",ylab="",ylim=c(-2000,2000),las=1)
for(i in 1:nrow(dvcsurr)){
  lines(1987:2016,dvcsurr[i,]-apply(dvcsurr,2,mean),col=colors[i])
}
lines(1987:2016,(apply(cty.list$Crashes[,-c(1:6)],2,sum)-apply(dvcsurr,2,mean)),type="b",pch=19,lwd=3)
Arrows(x0 = 1999,y0 = 0,y1=1500,x1=1999,arr.type = "triangle",arr.adj=1,arr.length=0.2,lwd=2,col="red")
Arrows(x0 = 1997,y0 = 0,y1=-1300,x1=1997,arr.type = "triangle",arr.adj=1,arr.length=0.2,lwd=2,col="red")
text(x=1999.5,y=1800,labels="1597 DVCs",font=2,cex=0.75,adj=0.05,col="red")
text(x=1997.5,y=-1800,labels="1421 DVCs",font=2,cex=0.75,adj=0.05,col="red")
mtext("Departure from Surrogate Mean",side=2,line=3.5,cex=0.9)
mtext("Year",side=1,line=1.5)
dev.off()
