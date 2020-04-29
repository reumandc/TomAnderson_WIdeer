cty.list<-readRDS("Results/cty.list.rds")
winter.clim<-readRDS("Results/winter.clim.rds")
climindex<-readRDS("Results/climindex.rds")

#clean data for middle and right panels
abun.dt<-wsyn::cleandat(cty.list$Abun,clev=5,times=minyear:maxyear)$cdat
dvc.dt<-wsyn::cleandat(cty.list$Crashes[,!is.na(colSums(cty.list$Crashes))],clev=5,times=1987:maxyear)$cdat

#generate wmfs and wpmfs
abun.wmf<-wsyn::wmf(abun.dt,times = 1981:maxyear)
dvc.wmf<-wsyn::wmf(dvc.dt,times = 1987:maxyear)

#4 x 3 panel figure

#parameters for figure layout, units inches
panwd.b<-1.25
panht.b<-panwd.b
#panwd.s<-panwd.b
#panht.s<-panwd.s/2.25
xht<-0.5
ywd<-0.5
gap<-0.1
totwd<-3*ywd+3*panwd.b+3*gap
totht<-xht+4*panht.b+7*gap
adjmt<-0.5

pdf("Results/FigXXX.pdf",height=totht,width=totwd)
##Left Column
#Panel A
par(fig=c((ywd)/totwd, #left side
          (ywd+panwd.b)/totwd, #right side
          (xht+3*panht.b+3*gap)/totht, #bottom
          (xht+4*panht.b+3*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)

#Insert fig

#Panel D
par(fig=c((ywd)/totwd, #left side
          (ywd+panwd.b)/totwd, #right side
          (xht+2*panht.b+2*gap)/totht, #bottom
          (xht+3*panht.b+2*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#Insert fig

#Panel G
par(fig=c((ywd)/totwd, #left side
          (ywd+panwd.b)/totwd, #right side
          (xht+1*panht.b+1*gap)/totht, #bottom
          (xht+2*panht.b+1*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#Insert fig

#Panel J
par(fig=c((ywd)/totwd, #left side
          (ywd+panwd.b)/totwd, #right side
          (xht)/totht, #bottom
          (xht+panht.b)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#Insert fig

##Middle Column
#Panel B
par(fig=c((2*ywd+panwd.b+gap)/totwd, #left side
          (2*ywd+2*panwd.b+gap)/totwd, #right side
          (xht+3*panht.b+3*gap)/totht, #bottom
          (xht+4*panht.b+3*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#driver (snow) time series

winter.clim.tmp<-lapply(winter.clim,function(x){x[(!is.na(rowMeans(x))),]})
plot(1981:2016,rep(NA,36),las=1,ylim=c(min(winter.clim.tmp$Snwd),max(winter.clim.tmp$Snwd)),xlab="",ylab="Winter Snow Depth (cm)",axes=F,frame=T)
for(i in 1:nrow(winter.clim.tmp$Snwd)){
  lines(1981:2016,winter.clim.tmp$Snwd[i,],col=rainbow(nrow(winter.clim.tmp$Snwd))[i])
}
mtext(text = "B)",side=3,adj=0.05,font=2,line=-1)
mtext("Snow Depth (m)",side=2,line=1.2)
axis(1,labels=rep(NA,6),at = seq(1980,2016,5))
axis(2,labels=seq(0,1000,200)/100,at = seq(0,1000,200),las=1)

par(new=T)

plot(minyear:maxyear,climindex$WinterMEI[1,],las=1,axes=F,type="l",xlab="",ylab="",lwd=2)
mtext("Winter MEI",side=4,line=1) 
axis(4, ylim=c(min(climindex$WinterMEI[1,]),max(climindex$WinterMEI[1,])), las=1)

#Panel E
par(fig=c((2*ywd+panwd.b+gap)/totwd, #left side
          (2*ywd+2*panwd.b+gap)/totwd, #right side
          (xht+2*panht.b+2*gap)/totht, #bottom
          (xht+3*panht.b+2*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#Deer raw time series
plot(1981:2016,rep(NA,36),las=1,ylim=c(min(cty.list$Abun),max(cty.list$Abun)),axes=F,frame=T,xlab="",ylab="")
for(i in 1:nrow(cty.list$Abun)){
  lines(1981:2016,cty.list$Abun[i,],col=rainbow(71)[i])
}
mtext("Deer (K)",side=2,line=1.2)
mtext(text = "E)",side=3,adj=0.05,font=2,line=-1)
axis(1,labels=rep(NA,6),at = seq(1980,2016,5))
axis(2,labels= format(seq(0,48,10),scientific=F),at = seq(0,48000,10000),las=1)

#Panel H
par(fig=c((2*ywd+panwd.b+gap)/totwd, #left side
          (2*ywd+2*panwd.b+gap)/totwd, #right side
          (xht+1*panht.b+1*gap)/totht, #bottom
          (xht+2*panht.b+1*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#Deer wavelet mean field 
zlimits<-range(Mod(abun.wmf$values),na.rm=T)
xlimits<-range(1981:2016)
l2ts<-log2(abun.wmf$timescales)
jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
colorfill<-grDevices::colorRampPalette(jetcolors)

#plot wmf
image(x=abun.wmf$times,y=l2ts,z=Mod(abun.wmf$values),xlim=xlimits,
      zlim=zlimits,col=colorfill(100),yaxt='n',xaxs='r',xaxt="n",yaxs='r',ylab="",xlab="")
ylocs <- pretty(abun.wmf$timescales, n = 8)
xlocs <- pretty(abun.wmf$times, n = 8)
#axis(1, at = xlocs, labels = rep(NA,length=length(xlocs)))
axis(1,labels=rep(NA,6),at = seq(1980,2016,5))
axis(2, at = log2(ylocs), labels = ylocs,las=1)
#mtext("Year",side=1,line=1.2)
mtext("Timescale (yrs)",side=2,line=1.2)
text(xlimits[1],max(l2ts),'H)',adj=c(0,1),font=2)
abline(h=c(log2(3),log2(7)),lty=2)

#add color bar 
par(new=T,fig=c((2*ywd+1.93*panwd.b+1*gap)/totwd, #left side
                (2*ywd+1.98*panwd.b+1*gap)/totwd, #right side
                (xht+1*panht.b+6*gap)/totht, #bottom
                (xht+2*panht.b-0.2*gap)/totht), #top
    mai=c(0,0,0,0))
cut.pts <- seq(zlimits[1], zlimits[2], length = length(colorfill(100)) + 1)
z <- (cut.pts[1:length(colorfill(100))] + cut.pts[2:(length(colorfill(100)) + 1)])/2
image(x = 1, y = z, z = matrix(z, ncol = length(colorfill(100)), nrow= 1),
      col = colorfill(100), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1), 
     mgp = c(3, 0.2, 0), las = 1, cex.axis = 0.75, tcl = -0.1)

#Panel K
par(fig=c((2*ywd+panwd.b+gap)/totwd, #left side
          (2*ywd+2*panwd.b+gap)/totwd, #right side
          (xht)/totht, #bottom
          (xht+panht.b)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#Regional total
plot(minyear:maxyear,apply(cty.list$Abun,2,sum),lwd=2,xlab="Year",ylab="Total Deer",axes=F,type="l",frame=T)
axis(2,labels= format(seq(0.6,1.3,0.2),scientific=F),at = seq(600000,1300000,200000),las=1)
axis(1,labels=seq(1980,2016,5),at = seq(1980,2016,5))
mtext("Deer (M)",side=2,line=1.2)
mtext("K)",font=2,side=3,line=-1,adj=0.05)

##Right Column: DVCs
#Panel C
par(fig=c((3*ywd+2*panwd.b+2*gap)/totwd, #left side
          (3*ywd+3*panwd.b+2*gap)/totwd, #right side
          (xht+3*panht.b+3*gap)/totht, #bottom
          (xht+4*panht.b+3*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

# Driver time series (deer repeated)
plot(1987:2016,rep(NA,30),las=1,ylim=c(min(cty.list$Abun[,-(1:6)]),max(cty.list$Abun[,-(1:6)])),xlab="",ylab="Abundance",axes=F,frame=T)
for(i in 1:nrow(cty.list$Abun)){
  lines(1981:2016,cty.list$Abun[i,],col=rainbow(71)[i])
}
mtext(text = "C)",side=3,adj=0.05,font=2,line=-1)
axis(1,labels=rep(NA,6),at = seq(1990,2016,5))
axis(2,labels= format(seq(0,48,10),scientific=F),at = seq(0,48000,10000),las=1)
mtext("Deer (K)",side=2,line=1.2)

#Panel F
par(fig=c((3*ywd+2*panwd.b+2*gap)/totwd, #left side
          (3*ywd+3*panwd.b+2*gap)/totwd, #right side
          (xht+2*panht.b+2*gap)/totht, #bottom
          (xht+3*panht.b+2*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#DVC raw time series
plot(1987:2016,rep(NA,30),las=1,frame=T,ylim=c(min(cty.list$Crashes,na.rm=T),max(cty.list$Crashes,na.rm=T)),xlim=c(1987,2016),axes=F,xlab="Year",ylab="DVCs")
for(i in 1:nrow(cty.list$Crashes)){
  lines(1987:2016,na.omit(cty.list$Crashes[i,]),col=rainbow(71)[i])
}
axis(1,labels=rep(NA,6),at = seq(1990,2016,5))
axis(2,labels= format(seq(0,1.4,0.2),scientific=F),at = seq(0,1400,200),las=1)
mtext("DVCs (K)",side=2,line=1.2)
mtext(text = "F)",side=3,line=-1,adj=0.05,font=2)

#Panel I
par(fig=c((3*ywd+2*panwd.b+2*gap)/totwd, #left side
          (3*ywd+3*panwd.b+2*gap)/totwd, #right side
          (xht+1*panht.b+1*gap)/totht, #bottom
          (xht+2*panht.b+1*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#DVC wavelet mean field
zlimits<-range(Mod(dvc.wmf$values),na.rm=T)
xlimits<-range(1987:2016)
l2ts<-log2(dvc.wmf$timescales)

#plot wmf
image(x=dvc.wmf$times,y=l2ts,z=Mod(dvc.wmf$values),xlim=xlimits,
      zlim=zlimits,col=colorfill(100),yaxt='n',xaxs='r',xaxt="n",yaxs='r',ylab="",xlab="")
ylocs <- pretty(dvc.wmf$timescales, n = 8)
xlocs <- pretty(dvc.wmf$times, n = 8)
axis(1, at = xlocs, labels = rep(NA,length=length(xlocs)))
axis(2, at = log2(ylocs), labels = ylocs,las=1)
#mtext("Year",side=1,line=1.2)
mtext("Timescale (yrs)",side=2,line=1.2)
text(xlimits[1],max(l2ts),'I)',adj=c(0,1),font=2)
abline(h=c(log2(3),log2(7)),lty=2)

#add color bar (NEED TO UPDATE COLOR BAR PAR )
par(new=T,fig=c((3*ywd+2.93*panwd.b+2*gap)/totwd, #left side
                (3*ywd+2.98*panwd.b+2*gap)/totwd, #right side
                (xht+1*panht.b+6*gap)/totht, #bottom
                (xht+2*panht.b-0.2*gap)/totht), #top
    mai=c(0,0,0,0),new=T)
cut.pts <- seq(zlimits[1], zlimits[2], length = length(colorfill(100)) + 1)
z <- (cut.pts[1:length(colorfill(100))] + cut.pts[2:(length(colorfill(100)) + 1)])/2
image(x = 1, y = z, z = matrix(z, ncol = length(colorfill(100)), nrow= 1),
      col = colorfill(100), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1), 
     mgp = c(3, 0.2, 0), las = 1, cex.axis = 0.75, tcl = -0.1)

#Panel L
par(fig=c((3*ywd+2*panwd.b+2*gap)/totwd, #left side
          (3*ywd+3*panwd.b+2*gap)/totwd, #right side
          (xht)/totht, #bottom
          (xht+panht.b)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)

#Total DVCs
plot(1987:2016,apply(cty.list$Crashes,2,sum,na.rm=T)[-c(1:6)],frame=T,type="l",lwd=2,ylab="",xlab="",las=1,axes=F,cex.axis=0.75)
axis(1,labels=seq(1985,2016,5),at = seq(1985,2016,5))
axis(2,labels= format(seq(16,24,2),scientific=F),at = seq(16000,24000,2000),las=1)
mtext("DVCs (K)",side=2,line=1.2)
mtext("L)",font=2,side=3,line=-1,adj=0.05)
dev.off()