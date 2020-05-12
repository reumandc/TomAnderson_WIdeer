#Load results and data
cty.list<-readRDS("Results/cty.list.rds")
winter.clim<-readRDS("Results/winter.clim.rds")
climindex<-readRDS("Results/climindex.rds")
abunsurr<-read.csv("Data/abunsurrsum.csv")
dvcsurr<-read.csv("Data/dvcsurrsum.csv")

#clean data for middle and right panels
abun.dt<-wsyn::cleandat(cty.list$Abun,clev=5,times=minyear:maxyear)$cdat
dvc.dt<-wsyn::cleandat(cty.list$Crashes[,!is.na(colSums(cty.list$Crashes))],clev=5,times=1987:maxyear)$cdat

#generate wmfs and wpmfs for middle and right panels
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

#do the sims for the left column of panels
w_here<-5 #period of synchronous part of noise
eta<-3
rhobc<-0
rhopop<-0.4
numsims<-3
lensim<-1024
N<-71 
set.seed(101)
#The below function comes from PedagogFig.R, so assuming this will have been loaded
#and retained in a previous knitting step
sims_casec<-simmodel4bc(N,w_here,rhobc,eta,rhopop,lensim,numsims)

tslen<-35
numts<-15 #length and number of time series actually displayed

#prep the stuff to plot on panels A, D, G and J, Example c=3 - time series with synchrony on long timescales
simnumtouse<-1
tm<-0:(tslen-1)
d3ns<-matrix(NA,N,length(tm))
d3<-matrix(NA,N,length(tm))
d3mm<-matrix(NA,N,length(tm))
for (counter in 1:N)
{
  d3ns[counter,]<-sims_casec$epsilon1[1:tslen,counter,simnumtouse]+
    sims_casec$epsilon2[1:tslen,counter,simnumtouse]+7*counter-1
  d3[counter,]<-sims_casec$pops[1:tslen,counter,simnumtouse]+7*counter-1
  d3mm[counter,]<-d3[counter,]-mean(d3[counter,])
}
d3mn<-apply(FUN=mean,X=d3,MARGIN=2)
d3sm<-apply(FUN=sum,X=d3,MARGIN=2)
wmf3<-wsyn::wmf(dat=d3mm,times=tm)

#Panel A
par(fig=c((ywd)/totwd, #left side
          (ywd+panwd.b)/totwd, #right side
          (xht+3*panht.b+3*gap)/totht, #bottom
          (xht+4*panht.b+3*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,.25,0),tcl=-.25)
xlimits<-range(tm)
ylimits<-range(d3ns[1:numts,])
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
colmap<-rep(c("black","red"),times=numts)#rainbow(numts)
adjmt<-0.5
plot(tm,d3ns[1,],type='l',xlim=xlimits,ylim=ylimits,col=colmap[1],
     xaxt='n',yaxt="n")
for (counter in 2:numts)
{
  lines(tm,d3ns[counter,],type='l',col=colmap[counter])
}
mtext(text = "A)",side=3,adj=0.05,font=2,line=-1)
axis(side=1,labels=FALSE)
axis(2,labels=TRUE,cex.axis=0.75)
mtext("Env. noise",side=2,line=1.2,cex=0.75)

#Panel D
par(fig=c((ywd)/totwd, #left side
          (ywd+panwd.b)/totwd, #right side
          (xht+2*panht.b+2*gap)/totht, #bottom
          (xht+3*panht.b+2*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)
xlimits<-range(tm)
ylimits<-range(d3[1:numts,])
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(tm,d3[1,],type='l',xlim=xlimits,ylim=ylimits,col=colmap[1],
     xaxt='n',yaxt="n")
for (counter in 2:numts)
{
  lines(tm,d3[counter,],type='l',col=colmap[counter])
}
mtext(text = "D)",side=3,adj=0.05,font=2,line=-1)
axis(side=1,labels=FALSE)
axis(2,labels=TRUE,cex.axis=0.75)
mtext("Populations",side=2,line=1.2,cex=0.75)

#Panel G
par(fig=c((ywd)/totwd, #left side
          (ywd+panwd.b)/totwd, #right side
          (xht+1*panht.b+1*gap)/totht, #bottom
          (xht+2*panht.b+1*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)

zlimits<-range(Mod(wmf3$values),na.rm=T)
jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
colorfill <- colorRampPalette(jetcolors)
l2ts<-log2(wmf3$timescales)
image(x=tm,y=l2ts,z=Mod(wmf3$values),xlim=xlimits,
      zlim=zlimits,col=colorfill(100),yaxt='n',xaxt="n",xaxs='r',yaxs='r')
ylocs <- pretty(wmf3$timescales, n = 8)
axis(2, at = log2(ylocs), labels = ylocs, las=1, cex.axis=0.75)
mtext("Timescale (yrs)",side=2,line=1.2,cex=0.75)
mtext(text = "G)",side=3,adj=0.05,font=2,line=-1)
axis(side=1,labels=FALSE)

#axis(2, at = log2(ylocs), labels = ylocs,las=1,cex.axis=0.75)

#Panel J
par(fig=c((ywd)/totwd, #left side
          (ywd+panwd.b)/totwd, #right side
          (xht)/totht, #bottom
          (xht+panht.b)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)

xlimits<-range(tm)
ylimits<-range(d3sm[1:numts])
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(tm,d3sm,type='l',xlim=xlimits,ylim=ylimits,
     yaxt='n',xaxt="n")
mtext("J)",font=2,side=3,line=-1,adj=0.05)
axis(1,labels=TRUE,las=1,cex.axis=0.75)
mtext("Tot. pop.",side=2,line=1.2,cex=0.75)
axis(2,labels=FALSE,cex.axis=0.75)
mtext("Year",side=1,line=1.5)

##Middle Column
#Panel B
par(fig=c((2*ywd+panwd.b+gap)/totwd, #left side
          (2*ywd+2*panwd.b+gap)/totwd, #right side
          (xht+3*panht.b+3*gap)/totht, #bottom
          (xht+4*panht.b+3*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)

#driver (snow) time series
ylimits<-c(0,1000)
winter.clim.tmp<-lapply(winter.clim,function(x){x[(!is.na(rowMeans(x))),]})
plot(1981:2016,rep(NA,36),las=1,ylim=ylimits,xlab="",ylab="",axes=F,frame=T)
for(i in 1:nrow(winter.clim.tmp$Snwd)){
  lines(1981:2016,winter.clim.tmp$Snwd[i,],col=rainbow(nrow(winter.clim.tmp$Snwd))[i])
}
mtext(text = "B)",side=3,adj=0.05,font=2,line=-1)
mtext("Snow Depth (m)",side=2,line=1.2,cex=0.75)
axis(1,labels=rep(NA,6),at = seq(1980,2016,5),cex.axis=0.75)
axis(2,labels=seq(0,1000,200)/100,at = seq(0,1000,200),las=1,cex.axis=0.75,tcl=-0.15)

par(new=T)

ylimits<-range(climindex$WinterMEI[1,])
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(minyear:maxyear,climindex$WinterMEI[1,],ylim=ylimits,las=1,axes=F,type="l",xlab="",ylab="",lwd=2)
mtext("Winter MEI",side=4,line=0.35,cex=0.75) 
axis(4, ylim=c(min(climindex$WinterMEI[1,]),max(climindex$WinterMEI[1,])), las=1,cex=0.75,cex.axis=0.75)

#Panel E
par(fig=c((2*ywd+panwd.b+gap)/totwd, #left side
          (2*ywd+2*panwd.b+gap)/totwd, #right side
          (xht+2*panht.b+2*gap)/totht, #bottom
          (xht+3*panht.b+2*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)

#Deer raw time series
plot(1981:2016,rep(NA,36),las=1,ylim=c(min(cty.list$Abun),max(cty.list$Abun)),axes=F,frame=T,xlab="",ylab="")
for(i in 1:nrow(cty.list$Abun)){
  lines(1981:2016,cty.list$Abun[i,],col=rainbow(71)[i])
}
mtext("Deer (K)",side=2,line=1.2,cex=0.75)
mtext(text = "E)",side=3,adj=0.05,font=2,line=-1)
axis(1,labels=rep(NA,6),at = seq(1980,2016,5),cex.axis=0.75)
axis(2,labels= format(seq(0,48,10),scientific=F),at = seq(0,48000,10000),las=1,cex.axis=0.75)

#Panel H
par(fig=c((2*ywd+panwd.b+gap)/totwd, #left side
          (2*ywd+2*panwd.b+gap)/totwd, #right side
          (xht+1*panht.b+1*gap)/totht, #bottom
          (xht+2*panht.b+1*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)

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
axis(1,labels=rep(NA,6),at = seq(1980,2016,5),cex.axis=0.75)
axis(2, at = log2(ylocs), labels = ylocs,las=1,cex.axis=0.75)
#mtext("Year",side=1,line=1.2)
mtext("Timescale (yrs)",side=2,line=1.2,cex=0.75)
#text(xlimits[1],max(l2ts),'H)',adj=c(0,1),font=2)
mtext(text = "H)",side=3,adj=0.05,font=2,line=-1)
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
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)

#Regional total
plot(1981:2016,rep(NA,36),ylab="",xlab="Year",ylim=c(-180000,200000),las=1,yaxt="n",cex.axis=0.75,cex=0.75)
lines(1981:2016,(apply(cty.list$Abun,2,sum)-apply(abunsurr,2,mean)),lwd=1,cex=0.75,type="l")
axis(2,labels= format(seq(-200,200,50),scientific=F),at = seq(-200000,200000,50000),las=1,cex.axis=0.75)
mtext(expression(Delta~"from Surrog. Mean (K)"),side=2,line=1.2,cex=0.75)
mtext("K)",font=2,side=3,line=-1,adj=0.05)
mtext("Year",side=1,line=1.5)

##Right Column: DVCs
#Panel C
par(fig=c((3*ywd+2*panwd.b+2*gap)/totwd, #left side
          (3*ywd+3*panwd.b+2*gap)/totwd, #right side
          (xht+3*panht.b+3*gap)/totht, #bottom
          (xht+4*panht.b+3*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)

# Driver time series (deer repeated)
ylimits<-c(min(cty.list$Abun[,-(1:6)]),max(cty.list$Abun[,-(1:6)]))
ylimits[2]<-ylimits[2]+.1*diff(ylimits)           
plot(1987:2016,rep(NA,30),las=1,ylim=ylimits,xlab="",ylab="Abundance",axes=F,frame=T)
for(i in 1:nrow(cty.list$Abun)){
  lines(1981:2016,cty.list$Abun[i,],col=rainbow(71)[i])
}
mtext(text = "C)",side=3,adj=0.05,font=2,line=-1)
axis(1,labels=rep(NA,6),at = seq(1990,2016,5),cex.axis=0.75)
axis(2,labels= format(seq(0,58,10),scientific=F),at = seq(0,58000,10000),las=1,cex.axis=0.75)
mtext("Deer (K)",side=2,line=0.85,cex=0.75)

#Panel F
par(fig=c((3*ywd+2*panwd.b+2*gap)/totwd, #left side
          (3*ywd+3*panwd.b+2*gap)/totwd, #right side
          (xht+2*panht.b+2*gap)/totht, #bottom
          (xht+3*panht.b+2*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)

#DVC raw time series
plot(1987:2016,rep(NA,30),las=1,frame=T,ylim=c(min(cty.list$Crashes,na.rm=T),max(cty.list$Crashes,na.rm=T)),xlim=c(1987,2016),axes=F,xlab="Year",ylab="DVCs")
for(i in 1:nrow(cty.list$Crashes)){
  lines(1987:2016,na.omit(cty.list$Crashes[i,]),col=rainbow(71)[i])
}
axis(1,labels=rep(NA,6),at = seq(1990,2016,5),cex.axis=0.75)
axis(2,labels= format(seq(0,1.4,0.2),scientific=F),at = seq(0,1400,200),las=1,cex.axis=0.75)
mtext("DVCs (K)",side=2,line=1.2,cex=0.75)
mtext(text = "F)",side=3,line=-1,adj=0.05,font=2)

#Panel I
par(fig=c((3*ywd+2*panwd.b+2*gap)/totwd, #left side
          (3*ywd+3*panwd.b+2*gap)/totwd, #right side
          (xht+1*panht.b+1*gap)/totht, #bottom
          (xht+2*panht.b+1*gap)/totht), #top
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)

#DVC wavelet mean field
zlimits<-range(Mod(dvc.wmf$values),na.rm=T)
xlimits<-range(1987:2016)
l2ts<-log2(dvc.wmf$timescales)

#plot wmf
image(x=dvc.wmf$times,y=l2ts,z=Mod(dvc.wmf$values),xlim=xlimits,
      zlim=zlimits,col=colorfill(100),yaxt='n',xaxs='r',xaxt="n",yaxs='r',ylab="",xlab="")
ylocs <- pretty(dvc.wmf$timescales, n = 6)
xlocs <- pretty(dvc.wmf$times, n = 8)
axis(1, at = xlocs, labels = rep(NA,length=length(xlocs)),cex.axis=0.75)
axis(2, at = log2(ylocs), labels = ylocs,las=1,cex.axis=0.75)
#mtext("Year",side=1,line=1.2)
mtext("Timescale (yrs)",side=2,line=1.2,cex=0.75)
#text(xlimits[1],max(l2ts),'I)',adj=c(0,1),font=2)
mtext(text = "I)",side=3,adj=0.05,font=2,line=-1)
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
    mai=c(0,0,0,0),mgp=c(3,0.25,0),tcl=-.25,new=T)

#Total DVCs
ylimits<-c(-2000,2000)
ylimits[2]<-ylimits[2]+.15*diff(ylimits)
plot(1987:2016,(apply(cty.list$Crashes[,-c(1:6)],2,sum)-apply(dvcsurr,2,mean)),type="l",lwd=1,cex=0.75,
     xlab="",ylab="",ylim=ylimits,las=1,yaxt="n")
axis(2,labels= format(seq(-2000,2500,500)/1000,scientific=F),at = seq(-2000,2500,500),las=1,cex.axis=0.75)
mtext(expression(Delta~"from Surrog. Mean (K)"),side=2,line=1.2,cex=0.75)
mtext("Year",side=1,line=1.5)
mtext("L)",font=2,side=3,line=-1,adj=0.05)

dev.off()
