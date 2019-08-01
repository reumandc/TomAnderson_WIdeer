#Work on a first theoretical example with AR(1) population dynamics, noise with 
#weak periodicity and possibly some synchrony (which is therefore timescale-specific 
#synchrony, if there is synchrony) depending on whether you are in case 1 or 2

#***setup

source("Code/my.spec.pgram.R")

#***settings

N<-70 #number of sampling locations
lensim<-1024 #length of simulations

#***the model

#Simulates the model
#
#Args
#N        The number of sampling locations
#rho      Covariance between locations for the inovations that go into producing 
#           the environmental noise
#r
#theta    These two args specify the AR(2) dynamics for noise
#rhopop   Lag-1 autocorrelation coefficient for population dynamics
#lensim   Length of simulations to do
#numsims  Number of simulations to do
#
#Output - a list with these elements, each 
#inovations - a lensim by N by numsims array, inovations for env noise
#epsilon - a lensim by N by numsims array, env noise
#pops - a lensim by N by numsims array, populations
#
simmodel1<-function(N,rho,r,theta,rhopop,lensim,numsims)
{
  #**generate the inovations for the generation of the noise
  
  Sig<-matrix(rho,N,N)
  diag(Sig)<-1
  inovations<-rmvnorm(3*lensim*numsims,sigma=Sig)
  dim(inovations)<-c(3*lensim,numsims,N)
  inovations<-aperm(inovations,c(1,3,2))
  
  #**now use these to generate the noise using an AR(2) process
  
  #get c1 and c2, the AR(2) parameters
  c1<-(1/r)*2*cos(theta)
  c2<-(-1/r^2) 
  
  #now make the noise
  epsilon<-array(0,dim=c(3*lensim,N,numsims))
  for (counter in 3:(3*lensim))
  {
    epsilon[counter,,]<-c1*epsilon[counter-1,,]+
      c2*epsilon[counter-2,,]+inovations[counter,,]
  }
  
  #throw away transients
  epsilon<-epsilon[(dim(epsilon)[1]-2*lensim+1):(dim(epsilon)[1]),,]
  inovations<-inovations[(dim(inovations)[1]-2*lensim+1):(dim(inovations)[1]),,]
  
  #**now put the noise into an AR(1) to get pops
  
  #do the sims  
  pops<-array(0,dim=c(2*lensim,N,numsims))
  for (counter in 2:(2*lensim))
  {
    pops[counter,,]<-rhopop*pops[counter-1,,]+epsilon[counter,,]
  }
  
  #throw away transients  
  pops<-pops[(dim(pops)[1]-lensim+1):(dim(pops)[1]),,]
  epsilon<-epsilon[(dim(epsilon)[1]-lensim+1):(dim(epsilon)[1]),,]
  inovations<-inovations[(dim(inovations)[1]-lensim+1):(dim(inovations)[1]),,]
  
  return(list(inovations=inovations,
              epsilon=epsilon,
              pops=pops))
}

#***do the sims

#do the sims
rhoC1<-0.8
r<-3
theta<-pi/2 #produces c1=0, c2=-1/9
rhopop<-0.15
numsims<-100
set.seed(101)
simsC1<-simmodel1(N=N,rho=rhoC1,r=r,theta=theta,rhopop=rhopop,
                lensim=lensim,numsims=numsims)
rhoC2<-0
simsC2<-simmodel1(N=N,rho=rhoC2,r=r,theta=theta,rhopop=rhopop,
                    lensim=lensim,numsims=numsims)

#***compute things you need for the figure

#first compute spectra of the noise, loc 1
spec_noise_C1<-my.spec.pgram(simsC1$epsilon[,1,],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
spec_noise_C2<-my.spec.pgram(simsC2$epsilon[,1,],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)

#theory for spectra
mu<-exp(-2*pi*complex(real=0,imaginary=1)*x)
c1<-(1/r)*2*cos(theta)
c2<-(-1/r^2) 
ys_C1_th<-(Mod(1/(1-c1*mu-c2*mu^2)))^2
ys_C2_th<-(Mod(1/(1-c1*mu-c2*mu^2)))^2

#and cospectra of the noise, locs 1 and 2
hC1<-my.spec.pgram(simsC1$epsilon[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
hC2<-my.spec.pgram(simsC2$epsilon[,1:2,1],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
xcs_C1<-hC1$freq
xcs_C2<-hC2$freq
ycs_C1<-matrix(NA,length(hC1$pgram[,1,2]),numsims)
ycs_C2<-matrix(NA,length(hC2$pgram[,1,2]),numsims)
ycs_C1[,1]<-Re(hC1$pgram[,1,2])
ycs_C2[,1]<-Re(hC2$pgram[,1,2])
for (counter in 2:numsims)
{
  hC1<-my.spec.pgram(simsC1$epsilon[,1:2,counter],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  hC2<-my.spec.pgram(simsC2$epsilon[,1:2,counter],detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
  ycs_C1[,counter]<-Re(hC1$pgram[,1,2])
  ycs_C2[,counter]<-Re(hC2$pgram[,1,2])
}

#theory for cospectra
ycs_C1_th<-rhoC1*(Mod(1/(1-c1*mu-c2*mu^2)))^2
ycs_C2_th<-rhoC2*(Mod(1/(1-c1*mu-c2*mu^2)))^2
  
#spectra of total pop
totpopsC1<-apply(FUN=sum,X=simsC1$pops,MARGIN=c(1,3))
totpopsC2<-apply(FUN=sum,X=simsC2$pops,MARGIN=c(1,3))
spec_totpops_C1<-my.spec.pgram(totpopsC1,detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)
spec_totpops_C2<-my.spec.pgram(totpopsC2,detrend=FALSE,taper=0,spans=c(91,71),plot=FALSE)

#theory for the spectra of the tot pop
ytp_C1_th<-(N+(N^2-N)*rhoC1)*(Mod(1/(1-c1*mu-c2*mu^2)))^2*(Mod(1/(1-rhopop*mu)))^2
ytp_C2_th<-(N+(N^2-N)*rhoC2)*(Mod(1/(1-c1*mu-c2*mu^2)))^2*(Mod(1/(1-rhopop*mu)))^2


  

  
#***make the figure

#setup
xht<-0.5
ywd<-0.5
gap<-0.3
panwd<-1.25
panht<-panwd
totwd<-4*ywd+4*panwd+4*gap
totht<-xht+2*panht+2*gap

#get the graphics device ready
png(file="Results/TheoryExample1.png",width=totwd,height=totht,units="in",res=600)
#pdf(file="Results/TheoryExample1.pdf",width=totwd,height=totht)

#spectra of noise in loc 1 for the different sims, C1
xpos<-0
ypos<-1
par(fig=c((ywd+xpos*(panwd+gap+ywd))/totwd,
          (ywd+xpos*(panwd+gap+ywd)+panwd)/totwd,
          (xht+ypos*(panht+gap))/totht,
          (xht+ypos*(panht+gap)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
x<-spec_noise_C1$freq
y<-spec_noise_C1$spec[,1]
xlimits<-c(0,.5)
ylimits<-range(spec_noise_C1$spec,spec_noise_C2$spec)
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=ylimits,xlim=xlimits) 
for (simc in 2:numsims)
{
  y<-spec_noise_C1$spec[,simc]
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
}
#axis(side=1,at=seq(from=0,to=.5,by=.1),labels=FALSE)
mtext(text="Noise spectra, loc 1",side=2,line=1.2,adj=3.75)
text(xlimits[1],ylimits[2],'A)',adj=c(0,1),font=2)

#add theory
lines(x,ys_C1_th,col="red") 

#spectra of noise in loc 1 for the different sims, C2
xpos<-0
ypos<-0
par(fig=c((ywd+xpos*(panwd+gap+ywd))/totwd,
          (ywd+xpos*(panwd+gap+ywd)+panwd)/totwd,
          (xht+ypos*(panht+gap))/totht,
          (xht+ypos*(panht+gap)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
x<-spec_noise_C2$freq
y<-spec_noise_C2$spec[,1]
ylimits<-range(spec_noise_C1$spec,spec_noise_C2$spec)
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=ylimits,xlim=xlimits) 
for (simc in 2:numsims)
{
  y<-spec_noise_C2$spec[,simc]
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
}
text(xlimits[1],ylimits[2],'B)',adj=c(0,1),font=2)
mtext(text="Frequency",side=1,line=1.2)

#add theory
lines(x,ys_C2_th,col="red") 

#cospectra of noise between locs 1 and 2 for the different sims, C1
xpos<-1
ypos<-1
par(fig=c((ywd+xpos*(panwd+gap+ywd))/totwd,
          (ywd+xpos*(panwd+gap+ywd)+panwd)/totwd,
          (xht+ypos*(panht+gap))/totht,
          (xht+ypos*(panht+gap)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
x<-xcs_C1
y<-ycs_C1[,1]
xlimits<-c(0,.5)
ylimits<-range(ycs_C1,ycs_C2)
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=ylimits,xlim=xlimits)
for (simc in 2:numsims)
{
  y<-ycs_C1[,simc]
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
}
text(xlimits[1],ylimits[2],'C)',adj=c(0,1),font=2)
mtext(text="Noise cospectra, locs 1, 2",side=2,line=1.2,adj=1.75)

#add theory
lines(x,ycs_C1_th,col="red") 

#cospectra of noise between locs 1 and 2 for the different sims, C2
xpos<-1
ypos<-0
par(fig=c((ywd+xpos*(panwd+gap+ywd))/totwd,
          (ywd+xpos*(panwd+gap+ywd)+panwd)/totwd,
          (xht+ypos*(panht+gap))/totht,
          (xht+ypos*(panht+gap)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
x<-xcs_C2
y<-ycs_C2[,1]
xlimits<-c(0,.5)
ylimits<-range(ycs_C1,ycs_C2)
ylimits[2]<-ylimits[2]+.2*diff(ylimits)
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=ylimits,xlim=xlimits)
for (simc in 2:numsims)
{
  y<-ycs_C2[,simc]
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
}
text(xlimits[1],ylimits[2],'D)',adj=c(0,1),font=2)
mtext(text="Frequency",side=1,line=1.2)

#add theory
lines(x,ycs_C2_th,col="red") 

#spectra of tot pop for the different sims, C1
xpos<-2
ypos<-1
par(fig=c((ywd+xpos*(panwd+gap+ywd))/totwd,
          (ywd+xpos*(panwd+gap+ywd)+panwd)/totwd,
          (xht+ypos*(panht+gap))/totht,
          (xht+ypos*(panht+gap)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
x<-spec_totpops_C1$freq
y<-spec_totpops_C1$spec[,1]
ylimits<-range(spec_totpops_C1$spec,spec_totpops_C2$spec)
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=ylimits) 
for (simc in 2:numsims)
{
  y<-spec_totpops_C1$spec[,simc]
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
}
mtext(text="Tot. pop. spectra",side=2,line=1.2,adj=-70)
text(xlimits[2],ylimits[2],'E)',adj=c(1,1),font=2)

#add theory for spectra of totpops for C1
lines(x,ytp_C1_th,col="red") 

#spectra of tot pop for the different sims, C2
xpos<-2
ypos<-0
par(fig=c((ywd+xpos*(panwd+gap+ywd))/totwd,
          (ywd+xpos*(panwd+gap+ywd)+panwd)/totwd,
          (xht+ypos*(panht+gap))/totht,
          (xht+ypos*(panht+gap)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
x<-spec_totpops_C2$freq
y<-spec_totpops_C2$spec[,1]
ylimits<-range(spec_totpops_C1$spec,spec_totpops_C2$spec)
plot(x,y,type="l",col=rgb(0,0,0,alpha=.1),
     ylim=ylimits)
for (simc in 2:numsims)
{
  y<-spec_totpops_C2$spec[,simc]
  lines(x,y,type="l",col=rgb(0,0,0,alpha=.1))
}
mtext(text="Frequency",side=1,line=1.2)
text(xlimits[2],ylimits[2],'F)',adj=c(1,1),font=2)

#add theory for spectra of totpops for C2
lines(x,ytp_C2_th,col="red") #the analytic expectation





#total pop time series, sim 1, C1
xpos<-3
ypos<-1
par(fig=c((ywd+xpos*(panwd+gap+ywd))/totwd,
          (ywd+xpos*(panwd+gap+ywd)+panwd)/totwd,
          (xht+ypos*(panht+gap))/totht,
          (xht+ypos*(panht+gap)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(1:3,1:3)


#total pop time series, sim 1, C2
xpos<-3
ypos<-0
par(fig=c((ywd+xpos*(panwd+gap+ywd))/totwd,
          (ywd+xpos*(panwd+gap+ywd)+panwd)/totwd,
          (xht+ypos*(panht+gap))/totht,
          (xht+ypos*(panht+gap)+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(1:3,1:3)


dev.off()