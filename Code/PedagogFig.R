#***
#For making the pedagogical figure, Fig. 1
#***

#***the models

#Simulates the model for case a
#
#Args
#N        The number of sampling locations
#w        The period of the oscillations
#rho      Written as rho_a in my notes on p. 14 and before that. This
#         is the pairwise correlation between locations of the white
#         part of the noise
#eta      This is the variance of the white part of the noise
#rhopop   Lag-1 autocorrelation coefficient for population dynamics
#lensim   Length of simulations to do
#numsims  Number of simulations to do
#
#Output - a list with these elements 
#epsilon1 - a lensim by N by numsims array, the first environmental noise (the periodic one)
#epsilon2 - a lensim by N by numsims array, the second environmental noise (the nonperiodic one)
#pops - a lensim by N by numsims array
simmodel4a<-function(N,w,rho,eta,rhopop,lensim,numsims)
{
  #**generate the first environmental noise (the periodic one)
  tims<-seq(from=0,to=2*lensim-1,by=1)
  epsilon1<-sin((2*pi/w)*array(tims,dim=c(2*lensim,N,numsims))+
                  array(2*pi*rep(runif(N*numsims),each=2*lensim),dim=c(2*lensim,N,numsims)))
  
  #generate the second environmental noise (the non-periodic one)
  Sig<-matrix(eta*rho,N,N)
  diag(Sig)<-eta
  epsilon2<-rmvnorm(2*lensim*numsims,sigma=Sig)
  dim(epsilon2)<-c(2*lensim,numsims,N)
  epsilon2<-aperm(epsilon2,c(1,3,2))
  
  #**now add the two environmental noises and use them to simulate pops
  
  #do the sims  
  pops<-array(0,dim=c(2*lensim,N,numsims))
  for (counter in 2:(2*lensim))
  {
    pops[counter,,]<-rhopop*pops[counter-1,,]+epsilon1[counter,,]+epsilon2[counter,,]
  }
  
  #throw away transients  
  pops<-pops[(dim(pops)[1]-lensim+1):(dim(pops)[1]),,]
  epsilon1<-epsilon1[(dim(epsilon1)[1]-lensim+1):(dim(epsilon1)[1]),,]
  epsilon2<-epsilon2[(dim(epsilon2)[1]-lensim+1):(dim(epsilon2)[1]),,]
  
  return(list(epsilon1=epsilon1,
              epsilon2=epsilon2,
              pops=pops))
}

#Simulates the model for cases b, c
#
#Args
#N        The number of sampling locations
#w        The period of the oscillations
#rho      Written as rho_bc in my notes on p. 14 and before that. This
#         is the pairwise correlation between locations of the white
#         part of the noise
#eta      This is the variance of the white part of the noise
#rhopop   Lag-1 autocorrelation coefficient for population dynamics
#lensim   Length of simulations to do
#numsims  Number of simulations to do
#
#Output - a list with these elements 
#epsilon1 - a lensim by N by numsims array, the first environmental noise (the periodic one)
#epsilon2 - a lensim by N by numsims array, the second environmental noise (the nonperiodic one)
#pops - a lensim by N by numsims array
simmodel4bc<-function(N,w,rho,eta,rhopop,lensim,numsims)
{
  #**generate the first environmental noise (the periodic one)
  tims<-seq(from=0,to=2*lensim-1,by=1)
  epsilon1<-sin((2*pi/w)*array(tims,dim=c(2*lensim,N,numsims))+
                  array(rep(2*pi*runif(numsims),each=2*lensim*N),dim=c(2*lensim,N,numsims)))
  
  #generate the second environmental noise (the non-periodic one)
  Sig<-matrix(eta*rho,N,N)
  diag(Sig)<-eta
  epsilon2<-rmvnorm(2*lensim*numsims,sigma=Sig)
  dim(epsilon2)<-c(2*lensim,numsims,N)
  epsilon2<-aperm(epsilon2,c(1,3,2))
  
  #**now add the two environmental noises and use them to simulate pops
  
  #do the sims  
  pops<-array(0,dim=c(2*lensim,N,numsims))
  for (counter in 2:(2*lensim))
  {
    pops[counter,,]<-rhopop*pops[counter-1,,]+epsilon1[counter,,]+epsilon2[counter,,]
  }
  
  #throw away transients  
  pops<-pops[(dim(pops)[1]-lensim+1):(dim(pops)[1]),,]
  epsilon1<-epsilon1[(dim(epsilon1)[1]-lensim+1):(dim(epsilon1)[1]),,]
  epsilon2<-epsilon2[(dim(epsilon2)[1]-lensim+1):(dim(epsilon2)[1]),,]
  
  return(list(epsilon1=epsilon1,
              epsilon2=epsilon2,
              pops=pops))
}

#***do the sims

#parameters
lensim<-1024
N<-71 #length and number of time series simulated

#do the sims case a
wa<-3
rhoa<-0 #use rhoa=0.75 and rhobc (below) equal to 0.25 for an example with equal 
#covariances (in expectation) between locations, just use rhoa=rhobc if you just
#want to turn on and off the periodic component of synchrony while leaving all else
#the same
eta<-3
rhopop<-0.4
numsims<-3
set.seed(101)
sims_casea<-simmodel4a(N,wa,rhoa,eta,rhopop,lensim,numsims)

#do the sims case b
wb<-3
rhobc<-0
eta<-3
rhopop<-0.4
numsims<-3
set.seed(101)
sims_caseb<-simmodel4bc(N,wb,rhobc,eta,rhopop,lensim,numsims)

#do the sims case c
wc<-10
rhopop<-0.4
numsims<-3
set.seed(101)
sims_casec<-simmodel4bc(N,wc,rhobc,eta,rhopop,lensim,numsims)

#***make the figure for the main text, which only show the time
#series and some wavelet results

#parameters for what to display
tslen<-35
numts<-15 #length and number of time series actually displayed

#parameters for figure layout, units inches
panwd.b<-1.25
panht.b<-panwd.b
panwd.s<-panwd.b
panht.s<-panwd.s/2.25
xht<-0.5
ywd<-0.5
gap<-0.1
totwd<-ywd+3*panwd.b+3*gap
totht<-xht+3*panht.b+panht.s+7*gap
colmap<-rep(c("black","red"),times=numts)#rainbow(numts)
adjmt<-0.5

#prep the stuff to plot, Example a=1 - time series with no synchrony
simnumtouse<-1
tm<-0:(tslen-1)
d1ns<-matrix(NA,N,length(tm))
d1<-matrix(NA,N,length(tm))
d1mm<-matrix(NA,N,length(tm))
for (counter in 1:N)
{
  d1ns[counter,]<-sims_casea$epsilon1[1:tslen,counter,simnumtouse]+
    sims_casea$epsilon2[1:tslen,counter,simnumtouse]+7*counter-1
  d1[counter,]<-sims_casea$pops[1:tslen,counter,simnumtouse]+7*counter-1
  d1mm[counter,]<-d1[counter,]-mean(d1[counter,])
}
d1mn<-apply(FUN=mean,X=d1,MARGIN=2)
d1sm<-apply(FUN=sum,X=d1,MARGIN=2)
wmf1<-wsyn::wmf(dat=d1mm,times=tm)
allcors1<-c()
for (i in 1:(numts-1))
{
  for (j in (i+1):numts)
  {
    allcors1<-c(allcors1,cor(d1[i,],d1[j,]))
  }
}

#prep the stuff to plot, Example b=2- time series with synchrony on short timescales
simnumtouse<-1
d2ns<-matrix(NA,N,length(tm))
d2<-matrix(NA,N,length(tm))
d2mm<-matrix(NA,N,length(tm))
for (counter in 1:N)
{
  d2ns[counter,]<-sims_caseb$epsilon1[1:tslen,counter,simnumtouse]+
    sims_caseb$epsilon2[1:tslen,counter,simnumtouse]+7*counter-1
  d2[counter,]<-sims_caseb$pops[1:tslen,counter,simnumtouse]+7*counter-1
  d2mm[counter,]<-d2[counter,]-mean(d2[counter,])
}
d2mn<-apply(FUN=mean,X=d2,MARGIN=2)
d2sm<-apply(FUN=sum,X=d2,MARGIN=2)
wmf2<-wsyn::wmf(dat=d2mm,times=tm)
allcors2<-c()
for (i in 1:(numts-1))
{
  for (j in (i+1):numts)
  {
    allcors2<-c(allcors2,cor(d2[i,],d2[j,]))
  }
}

#prep the stuff to plot, Example c=3 - time series with synchrony on long timescales
simnumtouse<-1
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
allcors3<-c()
for (i in 1:(numts-1))
{
  for (j in (i+1):numts)
  {
    allcors3<-c(allcors3,cor(d3[i,],d3[j,]))
  }
}

#get the graphics device ready
#pdf(file=paste("Results/Fig1.pdf"),width=totwd,height=totht)
#tiff(file=paste("Results/Fig1.tif"),width=totwd,height=totht,compression=c("lzw"),units="in",res=600)
png(file="Results/Fig1.png",width=totwd,height=totht,units="in",res=600)

#plot example a=1 - noise time series
xlimits<-range(tm)
ylimits.bns<-range(d1ns[1:numts,],d2ns[1:numts,],d3ns[1:numts,])
ylimits.bns[2]<-ylimits.bns[2]+.1*diff(ylimits.bns)
par(fig=c((ywd)/totwd,
          (ywd+panwd.b)/totwd,
          (xht+2*panht.b+panht.s+3*gap)/totht,
          (xht+3*panht.b+panht.s+3*gap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
plot(tm,d1ns[1,],type='l',xlim=xlimits,ylim=ylimits.bns,col=colmap[1],
     xaxt='n')
mtext("Env. noise",side=2,line=1.2)
for (counter in 2:numts)
{
  lines(tm,d1ns[counter,],type='l',col=colmap[counter])
}
text(xlimits[1],ylimits.bns[2],'A)',adj=c(0,1),font=2)
mtext("Scenario 1:\n env. asynchrony",side=3,cex = 0.75)

#Plot example a=1 - population time series
xlimits<-range(tm)
ylimits.b<-range(d1[1:numts,],d2[1:numts,],d3[1:numts,])
ylimits.b[2]<-ylimits.b[2]+.1*diff(ylimits.b)
par(fig=c((ywd)/totwd,
          (ywd+panwd.b)/totwd,
          (xht+1*panht.b+panht.s+2*gap)/totht,
          (xht+2*panht.b+panht.s+2*gap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(tm,d1[1,],type='l',xlim=xlimits,ylim=ylimits.b,col=colmap[1],
     xaxt='n')
mtext("Populations",side=2,line=1.2)
for (counter in 2:numts)
{
  lines(tm,d1[counter,],type='l',col=colmap[counter])
}
text(xlimits[1],ylimits.b[2],'D)',adj=c(0,1),font=2)

#plot example 1, mean or sum time series
par(fig=c(ywd/totwd,
          (ywd+panwd.b)/totwd,
          (xht)/totht,
          (xht+panht.s)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
ylimits.s<-range(d1sm,d2sm,d3sm)
ylimits.s[2]<-ylimits.s[2]+.25*diff(ylimits.s)
plot(tm,d1sm,type='l',xlim=xlimits,ylim=ylimits.s,
     yaxt='n')
axis(side=2,at=c(300,400),labels=c("300","400"),cex.axis=1)
text(xlimits[1],ylimits.s[2],'J)',adj=c(0,1),font=2)
mtext("Tot. pop.",side=2,line=1.2)

#plot example 1 - wavelet mean fields
par(fig=c(ywd/totwd,
          (ywd+panwd.b)/totwd,
          (xht+panht.s+gap)/totht,
          (xht+panht.s+gap+panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
zlimits<-range(Mod(wmf1$values),Mod(wmf2$values),Mod(wmf3$values),na.rm=T)
jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
colorfill <- colorRampPalette(jetcolors)
l2ts<-log2(wmf1$timescales)
image(x=tm,y=l2ts,z=Mod(wmf1$values),xlim=xlimits,
      zlim=zlimits,col=colorfill(100),yaxt='n',xaxt="n",xaxs='r',yaxs='r')
ylocs <- pretty(wmf1$timescales, n = 8)
axis(2, at = log2(ylocs), labels = ylocs)
mtext("Timescale (yrs)",side=2,line=1.2)
text(xlimits[1],max(l2ts),'G)',adj=c(0,1),font=2)

#Plot example 2 - noise time series
par(fig=c((ywd+panwd.b+gap)/totwd,
          (ywd+2*panwd.b+gap)/totwd,
          (xht+panht.s+3*gap+2*panht.b)/totht,
          (xht+panht.s+3*gap+3*panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
plot(tm,d2ns[1,],type='l',xlim=xlimits,ylim=ylimits.bns,col=colmap[1],
     yaxt='n',xaxt='n')
for (counter in 2:numts)
{
  lines(tm,d2ns[counter,],type='l',col=colmap[counter])
}
text(xlimits[1],ylimits.bns[2],'B)',adj=c(0,1),font=2)
mtext("Scenario 2:\n 3-year synchrony",side=3,cex = 0.75)

#Plot example 2 - population time series
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
text(xlimits[1],ylimits.b[2],'E)',adj=c(0,1),font=2)

#plot example 2 - mean or sum time series
par(fig=c((ywd+panwd.b+gap)/totwd,
          (ywd+2*panwd.b+gap)/totwd,
          (xht)/totht,
          (xht+panht.s)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
plot(tm,d2sm,type='l',xlim=xlimits,ylim=ylimits.s,
     yaxt='n')
text(xlimits[1],ylimits.s[2],'K)',adj=c(0,1),font=2)
mtext("Time step (yrs)",side=1,line=1.2)

#plot example 2 - wavelet mean fields
par(fig=c((ywd+panwd.b+gap)/totwd,
          (ywd+2*panwd.b+gap)/totwd,
          (xht+panht.s+gap)/totht,
          (xht+panht.s+gap+panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
image(x=tm,y=l2ts,z=Mod(wmf2$values),xlim=xlimits,
      zlim=zlimits,col=colorfill(100),yaxt='n',xaxt="n",xaxs='r',yaxs='r')
text(xlimits[1],max(l2ts),'H)',adj=c(0,1),font=2)

#Plot example 3 - noise time series
par(fig=c((ywd+2*panwd.b+2*gap)/totwd,
          (ywd+3*panwd.b+2*gap)/totwd,
          (xht+panht.s+3*gap+2*panht.b)/totht,
          (xht+panht.s+3*gap+3*panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
plot(tm,d3ns[1,],type='l',xlim=xlimits,ylim=ylimits.bns,col=colmap[1],
     xaxt='n',yaxt='n')
for (counter in 2:numts)
{
  lines(tm,d3ns[counter,],type='l',col=colmap[counter])
}
text(xlimits[1],ylimits.bns[2],'C)',adj=c(0,1),font=2)
mtext("Scenario 3:\n 10-year synchrony",side=3,cex = 0.75)

#Plot example 3 - population time series
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
text(xlimits[1],ylimits.b[2],'F)',adj=c(0,1),font=2)

#plot example 3 - mean or total time series
par(fig=c((ywd+2*panwd.b+2*gap)/totwd,
          (ywd+3*panwd.b+2*gap)/totwd,
          (xht)/totht,
          (xht+panht.s)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
plot(tm,d3sm,type='l',xlim=xlimits,ylim=ylimits.s,
     yaxt='n')
text(xlimits[1],ylimits.s[2],'L)',adj=c(0,1),font=2)

#plot example 3 - wavelet mean fields
par(fig=c((ywd+2*panwd.b+2*gap)/totwd,
          (ywd+3*panwd.b+2*gap)/totwd,
          (xht+panht.s+gap)/totht,
          (xht+panht.s+gap+panht.b)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
image(x=tm,y=l2ts,z=Mod(wmf3$values),xlim=xlimits,
      zlim=zlimits,col=colorfill(100),yaxt='n',xaxt="n",xaxs='r',yaxs='r')
text(xlimits[1],max(l2ts),'I)',adj=c(0,1),font=2)

#tom's color bar
par(new=T,fig=c((ywd+2.93*panwd.b+2*gap)/totwd,
          (ywd+2.98*panwd.b+2*gap)/totwd,
          (xht+6*gap+panht.s+gap)/totht,
          (xht+panht.b-0.5*gap+panht.s+gap)/totht),
    mai=c(0,0,0,0))
cut.pts <- seq(zlimits[1], zlimits[2], length = length(colorfill(100)) + 1)
z <- (cut.pts[1:length(colorfill(100))] + cut.pts[2:(length(colorfill(100)) + 1)])/2
image(x = 1, y = z, z = matrix(z, ncol = length(colorfill(100)), nrow= 1),
      col = colorfill(100), xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(2, at=round(seq(min(zlimits),max(zlimits),length=6),digits=1), 
     mgp = c(3, 0.2, 0), las = 1, cex.axis = 0.75, tcl = -0.1)

dev.off()



