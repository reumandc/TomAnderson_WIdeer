#***
#For making a pedagogical figure in support of Tom Anderson's paper intended for Science about
#synchrony in deer and DVCs.
#***

#parameters for the example itself
tslen<-35
numts<-15
ts1<-3
ts2<-10
nsstr<-1
sinstr<-1
set.seed(101) 

#parameters for figure layout, units inches
panwd.b<-1.25
panht.b<-panwd.b
panwd.s<-panwd.b
panht.s<-panwd.s/2.25
xht<-0.5
ywd<-0.5
gap<-0.1
totwd<-ywd+3*panwd.b+3*gap
totht<-xht+2*panht.b+2*gap+panht.s+gap
colmap<-rep(c("black","red"),times=numts)#rainbow(numts)
adjmt<-0.5
#pdf(file=paste("Fig1.pdf"),width=totwd,height=totht)
#tiff(file=paste("Fig1.tiff"),width=totwd,height=totht,compression=c("lzw"),units="in",res=600)
png(file="Results/Fig1.png",width=totwd,height=totht,units="in",res=600)

#Example 1 - time series with no synchrony
tm<-0:(tslen-1)
d1<-matrix(NA,numts,length(tm))
d1mm<-matrix(NA,numts,length(tm))
for (counter in 1:numts)
{
  d1[counter,]<-sinstr*sin(2*pi*tm/ts1+runif(1,0,2*pi))+nsstr*rnorm(tslen)+3*counter-1  
  d1mm[counter,]<-d1[counter,]-mean(d1[counter,])
}
d1mn<-apply(FUN=mean,X=d1,MARGIN=2)
d1sm<-apply(FUN=sum,X=d1,MARGIN=2)
wmf1<-wmf(dat=d1mm,times=tm)
allcors1<-c()
for (i in 1:(numts-1))
{
  for (j in (i+1):numts)
  {
    allcors1<-c(allcors1,cor(d1[i,],d1[j,]))
  }
}

#Example 2- time series with synchrony on short timescales
d2<-matrix(NA,numts,length(tm))
d2mm<-matrix(NA,numts,length(tm))
for (counter in 1:numts)
{
  d2[counter,]<-sinstr*sin(2*pi*tm/ts1)+nsstr*rnorm(tslen)+3*counter-1  
  d2mm[counter,]<-d2[counter,]-mean(d2[counter,])
}
d2mn<-apply(FUN=mean,X=d2,MARGIN=2)
d2sm<-apply(FUN=sum,X=d2,MARGIN=2)
wmf2<-wmf(dat=d2mm,times=tm)
allcors2<-c()
for (i in 1:(numts-1))
{
  for (j in (i+1):numts)
  {
    allcors2<-c(allcors2,cor(d2[i,],d2[j,]))
  }
}

#Example 3 - time series with synchrony on long timescales
d3<-matrix(NA,numts,length(tm))
d3mm<-matrix(NA,numts,length(tm))
for (counter in 1:numts)
{
  d3[counter,]<-sinstr*sin(2*pi*tm/ts2)+nsstr*rnorm(tslen)+3*counter-1  
  d3mm[counter,]<-d3[counter,]-mean(d3[counter,])
}
d3mn<-apply(FUN=mean,X=d3,MARGIN=2)
d3sm<-apply(FUN=sum,X=d3,MARGIN=2)
wmf3<-wmf(dat=d3mm,times=tm)
allcors3<-c()
for (i in 1:(numts-1))
{
  for (j in (i+1):numts)
  {
    allcors3<-c(allcors3,cor(d3[i,],d3[j,]))
  }
}

#Plot example 1 - time series
xlimits<-range(tm)
ylimits.b<-range(d1,d2,d3)
#ss<-diff(ylimits.b)*panht.s/panht.b/4.5
#ylimits.s<-mean(ylimits.b)+c(-ss/2-adjmt,ss/2-adjmt)
#ylimits.s<-range(d1sm,d2sm,d3sm)
ylimits.b[2]<-ylimits.b[2]+.1*diff(ylimits.b)
#ylimits.s[2]<-ylimits.s[2]+.2*diff(ylimits.s)
ylimits.s<-c(300,420)
par(fig=c(ywd/totwd,
          (ywd+panwd.b)/totwd,
          (xht+panht.b+panht.s+2*gap)/totht,
          (xht+panht.s+2*gap+2*panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
plot(tm,d1[1,],type='l',xlim=xlimits,ylim=ylimits.b,col=colmap[1],
     xaxt='n')
mtext("Populations",side=2,line=1.2)
for (counter in 2:numts)
{
  lines(tm,d1[counter,],type='l',col=colmap[counter])
}
text(xlimits[1],ylimits.b[2],'A)',adj=c(0,1),font=2)

#plot example 1, mean or sum time series
par(fig=c(ywd/totwd,
          (ywd+panwd.b)/totwd,
          (xht+panht.b+gap)/totht,
          (xht+panht.b+gap+panht.s)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
plot(tm,d1sm,type='l',xlim=xlimits,ylim=ylimits.s,
     xaxt='n',yaxt='n')
axis(side=2,at=c(300,400),labels=c("300","400"),cex.axis=1)
text(xlimits[1],ylimits.s[2],'D)',adj=c(0,1),font=2)
mtext("Tot. pop.",side=2,line=1.2)

#plot example 1 - wavelet mean fields
par(fig=c(ywd/totwd,
          (ywd+panwd.b)/totwd,
          (xht)/totht,
          (xht+panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
zlimits<-range(Mod(wmf1$values),Mod(wmf2$values),Mod(wmf3$values),na.rm=T)
jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
colorfill <- colorRampPalette(jetcolors)
l2ts<-log2(wmf1$timescales)
image(x=tm,y=l2ts,z=Mod(wmf1$values),xlim=xlimits,
      zlim=zlimits,col=colorfill(100),yaxt='n',xaxs='r',yaxs='r')
ylocs <- pretty(wmf1$timescales, n = 8)
axis(2, at = log2(ylocs), labels = ylocs)
mtext("Timescale (yrs)",side=2,line=1.2)
text(xlimits[1],max(l2ts),'G)',adj=c(0,1),font=2)

#Plot example 2 - time series
par(fig=c((ywd+panwd.b+gap)/totwd,
          (ywd+2*panwd.b+gap)/totwd,
          (xht+panht.s+2*gap+panht.b)/totht,
          (xht+panht.s+2*gap+2*panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
plot(tm,d2[1,],type='l',xlim=xlimits,ylim=ylimits.b,col=colmap[1],
     yaxt='n',xaxt='n')
for (counter in 2:numts)
{
  lines(tm,d2[counter,],type='l',col=colmap[counter])
}
text(xlimits[1],ylimits.b[2],'B)',adj=c(0,1),font=2)

#plot example 2 - mean or sum time series
par(fig=c((ywd+panwd.b+gap)/totwd,
          (ywd+2*panwd.b+gap)/totwd,
          (xht+panht.b+gap)/totht,
          (xht+panht.b+gap+panht.s)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
plot(tm,d2sm,type='l',xlim=xlimits,ylim=ylimits.s,
     yaxt='n',xaxt='n')
text(xlimits[1],ylimits.s[2],'E)',adj=c(0,1),font=2)

#plot example 2 - wavelet mean fields
par(fig=c((ywd+panwd.b+gap)/totwd,
          (ywd+2*panwd.b+gap)/totwd,
          (xht)/totht,
          (xht+panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
image(x=tm,y=l2ts,z=Mod(wmf2$values),xlim=xlimits,
      zlim=zlimits,col=colorfill(100),yaxt='n',xaxs='r',yaxs='r')
text(xlimits[1],max(l2ts),'H)',adj=c(0,1),font=2)
mtext("Time step (yrs)",side=1,line=1.2)

#Plot example 3 - time series
par(fig=c((ywd+2*panwd.b+2*gap)/totwd,
          (ywd+3*panwd.b+2*gap)/totwd,
          (xht+panht.s+2*gap+panht.b)/totht,
          (xht+panht.s+2*gap+2*panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
plot(tm,d3[1,],type='l',xlim=xlimits,ylim=ylimits.b,col=colmap[1],
     xaxt='n',yaxt='n')
for (counter in 2:numts)
{
  lines(tm,d3[counter,],type='l',col=colmap[counter])
}
text(xlimits[1],ylimits.b[2],'C)',adj=c(0,1),font=2)

#plot example 3 - mean or total time series
par(fig=c((ywd+2*panwd.b+2*gap)/totwd,
          (ywd+3*panwd.b+2*gap)/totwd,
          (xht+panht.b+gap)/totht,
          (xht+panht.b+gap+panht.s)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
plot(tm,d3sm,type='l',xlim=xlimits,ylim=ylimits.s,
     yaxt='n',xaxt='n')
text(xlimits[1],ylimits.s[2],'F)',adj=c(0,1),font=2)

#plot example 3 - wavelet mean fields
par(fig=c((ywd+2*panwd.b+2*gap)/totwd,
          (ywd+3*panwd.b+2*gap)/totwd,
          (xht)/totht,
          (xht+panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
image(x=tm,y=l2ts,z=Mod(wmf3$values),xlim=xlimits,
      zlim=zlimits,col=colorfill(100),yaxt='n',xaxs='r',yaxs='r')
text(xlimits[1],max(l2ts),'I)',adj=c(0,1),font=2)

#tom's color bar
par(new=T,fig=c((ywd+2.93*panwd.b+2*gap)/totwd,
          (ywd+2.98*panwd.b+2*gap)/totwd,
          (xht+6*gap)/totht,
          (xht+panht.b-0.5*gap)/totht),
    mai=c(0,0,0,0))
cut.pts <- seq(zlimits[1], zlimits[2], length = length(colorfill(100)) + 1)
z <- (cut.pts[1:length(colorfill(100))] + cut.pts[2:(length(colorfill(100)) + 1)])/2
image(x = 1, y = z, z = matrix(z, ncol = length(colorfill(100)), nrow= 1),
      col = colorfill(100), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(2, at=round(seq(min(zlimits),max(zlimits),length=6),digits=1), 
     mgp = c(3, 0.2, 0), las = 1, cex.axis = 0.75, tcl = -0.1)

dev.off()



