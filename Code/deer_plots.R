#Make figures for the manuscript
source("Functions/Fn_wsurfplot_modified.R")
source("Functions/Fn_wmfwt.R")

#Fig. 1
source("Code/PedagogFig.R")

#clean data for Figs 2 and 3
abun.dt<-Reumannplatz::CleanData(cty.list$Abun)$cleandat
dvc.dt<-Reumannplatz::CleanData(cty.list$Crashes[,!is.na(colSums(cty.list$Crashes))])$cleandat
#generate wavelet transforms
abun.wt<-warray(abun.dt,times=1981:2016)
dvc.wt<-warray(dvc.dt,times=1987:2016)

# Set up dimensions of wavelet mean field and phasor mean fields for Figs 2 and 3
tot.wd<-4.75
xht<-0.75   #height of x axis label region
ywd<-0.5    #width of y axis label region
zwd<-0.5    #width of z-axes label region
gap<-.2   #small gap 
pan.wd<-(tot.wd-ywd-gap-zwd) #large panel width parameter
pan.ht<-0.75*pan.wd+zwd #big ones are square
tot.ht<-3*pan.ht+xht+3*gap

png(filename="Results/Fig2.png",res=600,height=tot.ht,width=tot.wd,unit="in")
#tiff("Results/Fig2.tiff",res=600,compression=c("lzw"),height=tot.ht,width=tot.wd,unit="in")
#pdf("Results/Fig2.tiff",res=600,compression=c("lzw"),height=tot.ht,width=tot.wd,unit="in")
par(fig=c(ywd/tot.wd,(ywd+pan.wd)/tot.wd,(xht+2*pan.ht+2*gap)/tot.ht,(xht+3*pan.ht+2*gap)/tot.ht),mai=c(0,0,0,0),mgp=c(1,0.5,0))
wsurfplotTLA(abun.dt,times=1981:2016,colorbar=T,type="wmf",tsrange=c(3,7),
             xlab="Year",ylab="Timescales",xtcklab=rep("",9),smallplot=c(0.95,0.99,0.05,0.95))
par(fig=c(ywd/tot.wd,(ywd+pan.wd)/tot.wd,(xht+pan.ht+gap)/tot.ht,(xht+2*pan.ht+gap)/tot.ht),mai=c(0,0,0,0),new=T)
wsurfplotTLA(abun.dt,times=1981:2016,colorbar=T,type="wpmf",tsrange=c(3,7),siglevel = c(0.999),zlims = c(0,1),xlab="",xtcklab=rep("",9),ylab="Timescales",smallplot=c(0.95,0.99,0.05,0.95))
abline(h=c(3,7),lty=c(2,2))
par(fig=c(ywd/tot.wd,(ywd+pan.wd)/tot.wd,(xht+pan.ht+gap)/tot.ht,(xht+2*pan.ht+gap)/tot.ht),mai=c(0,0,0,0),new=T)
contour(Mod(wpmf(abun.dt,times=1981:2016)$wpmf),levels = 0.31,drawlabels=F,lwd=2,xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5),frame=F,font=2)
mtext(text = "Timescales",line = 1.3,side = 2,cex=1.25)
par(fig=c(ywd/tot.wd,(ywd+pan.wd)/tot.wd,xht/tot.ht,(xht+pan.ht)/tot.ht),mai=c(0,0,0,0),new=T)
syncexpplot(resp.wmf=wmfwt(abun.wt$wave.array),exp.sync = model.es$exp.sync,times=1981:2016,abun.wt$timescales,smallplot=c(0.95,0.99,0.05,0.95))
mtext(text = "Year",line = -2,side = 1,cex=1.25)
text(x=c(rep(1985.5,3)),y=c(log2(14),log2(7.5),log2(4)),labels=c("A)","B)","C)"),font=2)
dev.off()

png("Results/Fig3.png",res=600,height=tot.ht,width=tot.wd,unit="in")
#tiff("Results/Fig3.tiff",res=600,compression=c("lzw"),height=tot.ht,width=tot.wd,unit="in")
#pdf("Results/Fig3.tiff",res=600,compression=c("lzw"),height=tot.ht,width=tot.wd,unit="in")
par(fig=c(ywd/tot.wd,(ywd+pan.wd)/tot.wd,(xht+2*pan.ht+2*gap)/tot.ht,(xht+3*pan.ht+2*gap)/tot.ht),mai=c(0,0,0,0),mgp=c(1,0.5,0))
wsurfplotTLA(dvc.dt,times=1987:2016,colorbar=T,type="wmf",tsrange=c(3,7),ylab="TS",xlab="",xtcklab = rep("",8),smallplot=c(0.95,0.99,0.05,0.95))
par(fig=c(ywd/tot.wd,(ywd+pan.wd)/tot.wd,(xht+pan.ht+gap)/tot.ht,(xht+2*pan.ht+gap)/tot.ht),mai=c(0,0,0,0),new=T)
wsurfplotTLA(dvc.dt,times=1987:2016,colorbar=T,type="wpmf",tsrange=c(3,7),zlims = c(0,1),xlab="",siglevel = c(0.999),xtcklab = rep("",8),smallplot=c(0.95,0.99,0.05,0.95))
par(fig=c(ywd/tot.wd,(ywd+pan.wd)/tot.wd,(xht+pan.ht+gap)/tot.ht,(xht+2*pan.ht+gap)/tot.ht),mai=c(0,0,0,0),new=T)
contour(Mod(wpmf(dvc.dt,times=1987:2016)$wpmf),levels = 0.31,drawlabels=F,lwd=2,xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5),frame=F)
mtext(text = "Timescales",line = 1.3,side = 2,cex=1.25)
par(fig=c(ywd/tot.wd,(ywd+pan.wd)/tot.wd,xht/tot.ht,(xht+pan.ht)/tot.ht),mai=c(0,0,0,0),new=T)
syncexpplot(resp.wmf=wmfwt(dvc.wt$wave.array),exp.sync = dvcabun.es$pred.wmf,1987:2016,dvcabun.es$timescales,ylab = "",xlab="Year",smallplot=c(0.95,0.99,0.05,0.95))
mtext(text = "Year",line = -1.5,side = 1,cex=1.25)
text(x=c(rep(1990.5,3)),y=c(log2(11.3),log2(6.5),log2(3.7)),labels=c("A)","B)","C)"),font=2)
dev.off()

png("Results/FigS3.png",res=600,height=4800,width=3200)
#tiff("Results/FigS3.tiff",res=600,height=4800,width=3200,compression=c("lzw"))
source("Functions/Fn_phaseplot.R")
par(mfrow=c(4,2),mar=c(2.5,4,1.5,0.5),mgp=c(1.5,0.5,0),cex.lab=1.5,las=1)
phaseplot(climate.spcoh$WinterMEI.Abun$empirical,climate.spcoh$WinterMEI.Abun$timescales,type="pi",tsrange=c(3,7),xlab="",ylab="Phase")
mtext("A)",font=2,adj=0)
phaseplot(climate.spcoh$SummerMEI.Abun$empirical,climate.spcoh$SummerMEI.Abun$timescales,type="pi",tsrange=c(3,7),xlab="",ylab="")
mtext("B)",font=2,adj=0)
phaseplot(climate.spcoh$WinterPDO.Abun$empirical,climate.spcoh$WinterPDO.Abun$timescales,type="pi",tsrange=c(3,7),xlab="",ylab="Phase")
mtext("C)",font=2,adj=0)
phaseplot(winter.spcoh$Snwd.Abun$empirical,winter.spcoh$Snwd.Abun$timescales,showphase=F,type="pi",tsrange=c(3,7),ylab="",xlab="")
mtext("D)",font=2,adj=0)
rect(xleft=3,ybottom = -pi,xright = 5,ytop = pi,col=rgb(1,0,0,0.5),density=25,lwd=2,angle=-45)
mtext(line=0,adj=1,cex=0.75,bquote(paste("3-7, ",bar(theta)==.(0.854))))
mtext(line=1,adj=1,cex=0.75,bquote(paste("3-5, ",bar(theta)==.("0.870"))))
phaseplot(dvc.spcoh$empirical,dvc.spcoh$timescales,type="pi",tsrange=c(3,7),xlab="",ylab="Phase")
mtext("E)",font=2,adj=0)
phaseplot(adjdvc.spcoh$empirical,adjdvc.spcoh$timescales,type="pi",tsrange=c(3,7),xlab="Timescale",ylab="")
mtext("F)",font=2,adj=0)
phaseplot(hunter.spcoh$empirical,hunter.spcoh$timescales,showphase=F,type="pi",tsrange=c(3,7),xlab="Timescale",ylab="Phase")
mtext("G)",font=2,adj=0)
rect(xleft=2,ybottom = -pi,xright = 2.5,ytop = pi,col=rgb(1,0,0,0.5),density=25,lwd=2,angle=-45)
mtext(line=0,adj=1,cex=0.75,bquote(paste("3-7, ",bar(theta)==.(0.151))))
mtext(line=1,adj=1,cex=0.75,bquote(paste("2-2.5, ",bar(theta)==.(-0.836))))
dev.off()

#make plot of ranks by timescale
source("Functions/Fn_rankplot.R")
png("Results/FigS4.png",res=600,height=3000,width=4800)
#tiff("Results/FigS4.tiff",res=600,height=3000,width=4800,compression=c("lzw"))
par(mfrow=c(1,2),mar=c(3.5,4,1,0),mgp=c(2.5,0.5,0))
plot(climate.res$WinterMEI.Abun$timescale,climate.res$WinterMEI.Abun$emp.rank,las=1,ylim=c(0.5*nsurrogs,max(nsurrogs)),type="l",lwd=2,ylab="Rank",xlab="Timescale")
abline(h=0.95*nsurrogs,col="red",lty=2)
lines(climate.res$SummerMEI.Abun$timescales,climate.res$SummerMEI.Abun$emp.rank,lwd=2,col="green")
lines(climate.res$WinterPDO.Abun$timescales,climate.res$WinterPDO.Abun$emp.rank,lty=1,lwd=2,col="orange")
lines(winter.res$Snwd.Abun$timescales,winter.res$Snwd.Abun$emp.rank,lty=1,lwd=2,col="blue")
legend("topright",lty=c(1,1),c("Winter MEI","Winter PDO","Summer MEI","Snow Depth"),cex=0.75,bty = "n",col=c("black","orange","green","blue"))
mtext("A)",font=2,adj=0)

plot(dvc.res$timescales,dvc.res$emp.rank,type="l",lwd=2,ylab="",xlab="Timescale",las=1,ylim=c(0.5*nsurrogs,max(nsurrogs)))
lines(adjdvc.res$timescales,adjdvc.res$emp.rank,lty=1,lwd=2,col="purple")
abline(h=0.95*nsurrogs,col="red",lty=2)
lines(hunter.res3_7$timescales,hunter.res3_7$emp.rank,lty=1,lwd=2,col="darkgray")
legend("topright",c("Abun-DVCs","Abun-Adj. DVCs","Hunters-Abun"),lty=c(1,1,1),
       col=c("black","purple","darkgray"),bty="n",cex=0.75)
mtext("B)",font=2,adj=0)
dev.off()

# #plot weather and climate indices ranks to see whether relationship exists outside of the 3-7 year timescales
# plot(2:14,rep(NA,13),ylim=c(0.5*nsurrogs,max(nsurrogs)))
# for(i in 1:length(weath.climind.spcoh)){
#   lines(weath.climind.res$WinterMEI.Snwd$timescales,weath.climind.res[[i]]$emp.rank,col=rainbow(40)[i],main=i)
#   abline(h=0.95*nsurrogs,col="red",lty=2)
# }
# plot(weath.climind.res$WinterMEI.Snwd$timescales,weath.climind.res$WinterMEI.Snwd$emp.rank,type="l")
# lines(weath.climind.res$WinterPDO.Snwd$timescales,weath.climind.res$WinterPDO.Snwd$emp.rank,col="green")
# abline(h=0.95*nsurrogs,col="red",lty=2)

#Hunter Plots
png("Results/FigS5.png",res=600,height=4800,width=3000)
#tiff("Results/FigS5.tiff",res=600,compression=c("lzw"),height=4800,width=3000)
hunters.tmp<-cty.list$Hunters[,12:dim(cty.list$Hunters)[2]]
hunters.tmp<-hunters.tmp[(!row.names(hunters.tmp)%in%cwd) & !is.na(rowMeans(hunters.tmp)),]
hunters.dt<-Reumannplatz::CleanData(hunters.tmp)$cleandat
par(mfrow=c(3,1),mar=c(2.5,3,0,4),mgp=c(1.5,0.5,0))
wsurfplotTLA(hunters.dt,times=1992:2016,colorbar=T,type="wmf",ylab="",xlab="")
mtext("A)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(hunters.dt,times=1992:2016,colorbar=T,zlims = c(0,1),type="wpmf",xlab="",siglevel=0.999)
mtext("B)",adj=0.05,line=-1.2,side=3,font=2)
abun.dt1<-Reumannplatz::CleanData(cty.list$Abun[rownames(hunters.tmp),12:dim(cty.list$Hunters)[2]])
abun.wt1<-warray(abun.dt1$cleandat,times=1992:2016)
syncexpplot(resp.wmf=wmfwt(abun.wt1$wave.array),exp.sync = hunterabun.es2_2.5$pred.wmf,1992:2016,hunterabun.es2_2.5$timescales,ylab = "",xlab="Year")
mtext("C)",adj=0.05,line=-1.2,side=3,font=2)
dev.off()

#climate surface plots
png("Results/FigS8.png",res=600,height=9,width=7,unit="in")
#tiff("Results/FigS8.tiff",res=600,compression=c("lzw"),height=9,width=7,unit="in")
par(mfrow=c(3,2),mar=c(2.5,3,0,4),mgp=c(1.5,0.5,0))
wsurfplotTLA(Reumannplatz::CleanData(climindex$WinterPDO[1,])$cleandat,times=1981:2016,colorbar=T,type="power",xlab="")
mtext("A)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(Reumannplatz::CleanData(climindex$SummerPDO[1,])$cleandat,times=1981:2016,colorbar=T,type="power",xlab="",ylab="")
mtext("B)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(Reumannplatz::CleanData(climindex$WinterNAO[1,])$cleandat,times=1981:2016,colorbar=T,type="power",xlab="")
mtext("C)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(Reumannplatz::CleanData(climindex$SummerNAO[1,])$cleandat,times=1981:2016,colorbar=T,type="power",ylab="",xlab="")
mtext("D)",adj=0.05,font=2,line=-1.2,side=3)
wsurfplotTLA(Reumannplatz::CleanData(climindex$WinterMEI[1,])$cleandat,times=1981:2016,colorbar=T,type="power",xlab="Year",ylab="Timescale")
mtext("E)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(Reumannplatz::CleanData(climindex$SummerMEI[1,])$cleandat,times=1981:2016,colorbar=T,type="power",xlab="Year",ylab="Timescale")
mtext("F)",adj=0.05,line=-1.2,side=3,font=2)
dev.off()

#weather surface plots
png("Results/FigS7.png",res=600,height=9,width=7,unit="in")
#tiff("Results/FigS7.tiff",res=600,compression=c("lzw"),height=9,width=7,unit="in")
winter.clim.dt<-lapply(winter.clim,function(x){x<-Reumannplatz::CleanData(x[!(is.na(rowMeans(x))),])$cleandat;x})
par(mfrow=c(4,2),mar=c(2.5,3,0,4),mgp=c(1.5,0.5,0))
wsurfplotTLA(winter.clim.dt$Tmin,times=1981:2016,colorbar=T,type="wmf",xlab="")
mtext("A)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(winter.clim.dt$Tmin,times=1981:2016,colorbar=T,type="wpmf",xlab="",siglevel = 0.999)
mtext("B)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(winter.clim.dt$Tmax,times=1981:2016,colorbar=T,type="wmf",ylab="",xlab="")
mtext("C)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(winter.clim.dt$Tmax,times=1981:2016,colorbar=T,type="wpmf",ylab="",xlab="",siglevel = 0.999)
mtext("D)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(winter.clim.dt$Snwd,times=1981:2016,colorbar=T,type="wmf")
mtext("E)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(winter.clim.dt$Snwd,times=1981:2016,colorbar=T,type="wpmf",siglevel = 0.999)
mtext("F)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(winter.clim.dt$Prcp,times=1981:2016,colorbar=T,type="wmf",ylab="")
mtext("G)",adj=0.05,font=2,line=-1.2,side=3)
wsurfplotTLA(winter.clim.dt$Prcp,times=1981:2016,colorbar=T,type="wpmf",ylab="",siglevel = 0.999)
mtext("H)",adj=0.05,font=2,line=-1.2,side=3)
dev.off()

png("Results/FigS6.png",res=600,height=9,width=7,unit="in")
#tiff("Results/FigS6.tiff",res=600,compression=c("lzw"),height=9,width=7,unit="in")
source("Functions/Fn_wsurfplot_modified.R")
par(mfrow=c(3,2),mar=c(2.5,3,0,4),mgp=c(1.5,0.5,0))
wsurfplotTLA(Reumannplatz::CleanData(usda.list$Abun)$cleandat,times=1981:2016,colorbar=T,type="wmf",xlab="",ylab="")
mtext("A)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(Reumannplatz::CleanData(usda.list$Crashes[,!is.na(colSums(usda.list$Crashes))])$cleandat,times=1987:2016,colorbar=T,type="wmf",xlab="",ylab="")
mtext("B)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(Reumannplatz::CleanData(usda.list$Abun)$cleandat,times=1981:2016,colorbar=T,type="wpmf",siglevel = 0.999,xlab="")
mtext("C)",adj=0.05,line=-1.2,side=3,font=2)
wsurfplotTLA(Reumannplatz::CleanData(usda.list$Crashes[,!is.na(colSums(usda.list$Crashes))])$cleandat,times=1987:2016,colorbar=T,type="wpmf",siglevel = 0.999,xlab="",ylab="")
mtext("D)",adj=0.05,line=-1.2,side=3,font=2)
syncexpplot(resp.wmf=wmfwt(abun.wt$wave.array),exp.sync = usda.model.es$exp.sync,times=1981:2016,abun.wt$timescales,ylab="")
mtext("E)",adj=0.05,line=-1.2,side=3,font=2)
syncexpplot(resp.wmf=wmfwt(dvc.wt$wave.array),exp.sync = dvcabun.es$pred.wmf,times=1987:2016,dvc.wt$timescales,ylab="")
mtext("F)",adj=0.05,line=-1.2,side=3,font=2)
dev.off()

#raw data plots
png("Results/FigS1.png",res=600,height=3600,width=3600)
#tiff("Results/FigS1.tiff",res=600,compression=c("lzw"),height=3600,width=3600)
par(mfrow=c(3,2),mar=c(2.5,5,1.2,0),mgp=c(3,0.5,0))
plot(1981:2016,rep(NA,36),las=1,ylim=c(min(cty.list$Abun),max(cty.list$Abun)),xlab="",ylab="Abundance")
for(i in 1:nrow(cty.list$Abun)){
  lines(1981:2016,cty.list$Abun[i,],col=rainbow(71)[i])
}
mtext(text = "A)",side=3,adj=0,font=2)
plot(1987:2016,rep(NA,30),las=1,ylim=c(min(cty.list$Crashes,na.rm=T),max(cty.list$Crashes,na.rm=T)),xlim=c(1981,2016),xlab="Year",ylab="DVCs")
for(i in 1:nrow(cty.list$Crashes)){
  lines(1987:2016,na.omit(cty.list$Crashes[i,]),col=rainbow(71)[i])
}
mtext(text = "B)",side=3,adj=0,font=2)

winter.clim.tmp<-lapply(winter.clim,function(x){x[(!is.na(rowMeans(x))),]})
plot(1981:2016,rep(NA,36),las=1,ylim=c(min(winter.clim.tmp$Snwd),max(winter.clim.tmp$Snwd)),xlab="",ylab="Winter Snow Depth (cm)")
for(i in 1:nrow(winter.clim.tmp$Snwd)){
  lines(1981:2016,winter.clim.tmp$Snwd[i,],col=rainbow(nrow(winter.clim.tmp$Snwd))[i])
}
mtext(text = "C)",side=3,adj=0,font=2)

plot(mei$YEAR[mei$YEAR>1980],win.mei.mat[1,],las=1,col="blue",type="l",xlab="Year",ylab="Index Value",lwd=2)
lines(mean.pdo$Year[mean.pdo$Year>1980],win.pdo.mat[1,],col="green",lwd=2)
legend("topright",bty="n",c("Winter MEI","Winter PDO"),col=c("blue","green"),lty=c(1,1),cex=0.75)
mtext(text = "D)",side=3,adj=0,font=2)
mtext(text="Year",side=1,line=1.5,cex=0.75)

cty.list.tmp<-cty.list$Hunters[,12:dim(cty.list$Hunters)[2]]
cty.list.tmp<-cty.list.tmp[(!row.names(cty.list.tmp)%in%cwd) & !is.na(rowMeans(cty.list.tmp)),]
plot(1992:2016,rep(NA,25),ylim=c(min(cty.list.tmp,na.rm=T),max(cty.list.tmp,na.rm=T)),xlab="Year",ylab="Hunters",xlim=c(1981,2016),las=1)
for(i in 1:nrow(cty.list.tmp)){
  lines(1992:2016,na.omit(cty.list.tmp[i,]),col=rainbow(53)[i])
}
mtext(text = "E)",side=3,adj=0,font=2)
mtext(text="Year",side=1,line=1.5,cex=0.75)
dev.off()

#make Lawrence's plot showing how synchronous impacts statewide fluctuations
abunsurr<-read.csv("Data/abunsurrsum.csv")
dvcsurr<-read.csv("Data/dvcsurrsum.csv")

#Figure 4
#***plotting dimensions, units inches
tot.wd<-3.5
xht<-.5   #height of x axis label region
ywd<-1    #width of y axis label region
gap<-.1   #small gap 
pan.wd.big<-(tot.wd-ywd-gap) #large panel width parameter
pan.ht.big<-pan.wd.big #big ones are square
pan.wd.small<-pan.wd.big #small panel width param
pan.ht.small<-0.33*pan.ht.big #small panel height param
tot.ht<-2*pan.ht.big+2*pan.ht.small+2*xht+4*gap
#pdf("Results/Fig4.pdf",width = tot.wd,height = tot.ht)
#tiff("Results/Fig4.tiff",res=600,units="in",width = tot.wd,height = tot.ht,compression=c("lzw"))
png("Results/Fig4.png",res=600,units="in",width = tot.wd,height = tot.ht)

#Deer- little panel
par(fig=c(ywd/tot.wd,
          (ywd+pan.wd.small)/tot.wd,
          (2*pan.ht.big+pan.ht.small+2*xht+2*gap)/tot.ht,
          (2*pan.ht.big+2*pan.ht.small+2*xht+2*gap)/tot.ht),
    mai=c(0,0,0,0))
colors<-gray.colors(25)[sample(1:25, 999, replace=TRUE)]
plot(1981:2016,rep(NA,36),ylim=c(min(abunsurr),max(abunsurr)),ylab="",xlab="",las=1,axes=F)
for(i in 1:nrow(abunsurr)){
  lines(1981:2016,abunsurr[i,],col=colors[i])
}
lines(1981:2016,apply(cty.list$Abun,2,sum),lwd=2)
axis(2,labels= format(seq(600000,1300000,200000),scientific=T),at = seq(600000,1300000,200000),las=1,tck=-.05,cex.axis=0.75)
axis(1,labels=c(rep("",8)),at = seq(1980,2016,5),tck=-0.05,cex.axis=0.75)
mtext("Deer",side=2,line=3.5,cex=0.75)
mtext("A)",font=2,side=3,line=0,adj=0)
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
text(x=2001,y=160000,labels="159054 deer",font=2,cex=0.75,adj=0,col="red")
text(x=1999,y=-174000,labels="174339 deer",font=2,cex=0.75,adj=0,col="red")
mtext("Departure from Surrogate Mean",side=2,line=3.5,cex=0.9)

#DVC Plot
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
mtext("B)",font=2,side=3,line=0,adj=0)
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
text(x=1999.5,y=1800,labels="1597 DVCs",font=2,cex=0.75,adj=0,col="red")
text(x=1997.5,y=-1800,labels="1421 DVCs",font=2,cex=0.75,adj=0,col="red")
mtext("Departure from Surrogate Mean",side=2,line=3.5,cex=0.9)
mtext("Year",side=1,line=1.5)
dev.off()
