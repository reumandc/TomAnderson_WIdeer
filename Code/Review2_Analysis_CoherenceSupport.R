
#***
#Functions
#***

#Takes the fft of each row of the matrix x
#
#Args
#x      A matrix
#
myfft<-function(x)
{
  return(t(mvfft(t(x))))
}

#Inverse of the above
#
#Args
#x      A matrix
#
imyfft<-function(x)
{
  return(t(mvfft(t(x),inverse=TRUE))/dim(x)[2])
}

#Removes periodicities outside the specified range by zero-ing out the appropriate
#components of the Fourier transform
#
#Args
#x      A matrix, the operation is performed on each row
#ts     Range of frequencies to keep, inclusive. Units cycles per sampling interval.
#
#Note: The argument ts must be a length-2 vector c(x,y) with 0<x<y<=1/2, but no error
#checking is done to enforce this. Time series means are also removed.
remper<-function(x,ts)
{
  #Transform
  fx<-myfft(x)
 
  #Remove components outside the specified range
  freqs<-((1:(dim(x)[2]))-1)/(dim(x)[2])
  freqinds1<-which(freqs>=ts[1] & freqs<=ts[2])
  freqinds2<-rev(dim(x)[2]-freqinds1+2)  #keeps the symmetry
  freqinds<-c(freqinds1,freqinds2)
  fxnew<-matrix(0,dim(fx)[1],dim(fx)[2])
  fxnew[,freqinds]<-fx[,freqinds]
  
  #Reverse transform and return the result
  h<-imyfft(fxnew)
  if (max(abs(Im(h)))>1e-12){stop("Error in remper")}
  return(Re(h)) 
}

#***
#Tests of functions
#***

#x<-matrix(rnorm(10*100),10,100)
#fx<-myfft(x)
#max(abs(fx[1,]-fft(x[1,])))
#max(abs(fx[2,]-fft(x[2,])))
#max(abs(fx[3,]-fft(x[3,])))

#ifx<-imyfft(fx)
#max(abs(x-ifx))

#x<-matrix(rnorm(10*100),10,100)
#xrem<-remper(x,c(1/4,1/3))
#fx<-myfft(x)
#fxrem<-myfft(xrem)
#freqs<-(0:99)/100
#max(abs(fx[,freqs>=1/4 & freqs<=1/3]-fxrem[,freqs>=1/4 & freqs<=1/3]))
#max(abs(fx[,freqs>=2/3 & freqs<=3/4]-fxrem[,freqs>=2/3 & freqs<=3/4]))
#max(abs(fxrem[,freqs<1/4]))
#max(abs(fxrem[,freqs>1/3 & freqs<2/3]))
#max(abs(fxrem[,freqs>3/4]))
#max(Im(xrem))

#x<-matrix(rnorm(10*36),10,36)
#xrem<-remper(x,c(1/4,1/3))
#fx<-myfft(x)
#fxrem<-myfft(xrem)
#freqs<-(0:35)/36
#max(abs(fx[,freqs>=1/4 & freqs<=1/3]-fxrem[,freqs>=1/4 & freqs<=1/3]))
#max(abs(fx[,freqs>=2/3 & freqs<=3/4]-fxrem[,freqs>=2/3 & freqs<=3/4]))
#max(abs(fxrem[,freqs<1/4]))
#max(abs(fxrem[,freqs>1/3 & freqs<2/3]))
#max(abs(fxrem[,freqs>3/4]))
#max(Im(xrem))

#***
#***Load the data
#***

#pull in the raw (un-transformed) deer data
d<-readRDS(file="Results/cty.list.rds")
deer<-d$Abun
totdeer<-apply(FUN=sum,X=deer,MARGIN=2)
deeryr<-1981:2016

#pull in the raw (un-transformed) DVC data
dvcs<-d$Crashes
dvcs<-dvcs[,7:(dim(dvcs)[2])]
totdvcs<-apply(FUN=sum,X=dvcs,MARGIN=2)
dvcyr<-1987:2016

#pull in raw (un-transformed) winter weather data, same dimensions as deer, though with more NAs (see below)
ww<-readRDS(file="Results/winter.clim.rds")
snow<-ww$Snwd
snow<-snow[c(1:54,56:59,55,60:71),]

#pull in raw (un-transformed) winter climate index data
wc<-readRDS(file="Results/climindex.rds")
mei<-wc$WinterMEI #is the matrix of repeating values, counties X years, same dimensions as deer
mei<-mei[c(1:54,56:59,55,60:71),] #except for getting the rownames right, this makes no difference, because 
                                  #data are the same in all locations 

#cbind(rownames(deer),rownames(dvcs),rownames(snow),rownames(mei))

#***
#***Deer/snow comparisons
#***

#***Data prep

#Trim snow and deer data so both variables are the same size, due to missing values in snow data.
#These datasets are fit for making comparisons between deer and snow depth.
snow_noNA<-snow[!is.na(rowMeans(snow)),] #take out locations missing snow data
deer_forsnow<-deer[!is.na(rowMeans(snow)),] #take out counties where snow is NA
#dim(snow_noNA)
#dim(deer_forsnow)
#sum(is.na(deer))
#sum(is.na(snow_noNA))

#clean and transform
snow_noNA_cl<-wsyn::cleandat(snow_noNA,clev=5,times=deeryr)$cdat
deer_forsnow_cl<-wsyn::cleandat(deer_forsnow,clev=5,times=deeryr)$cdat

#***Now do a cross correlation analysis - the goal here is to compute statistics having to 
#do with cross correlation, for both real data and surrogates, and compare

#compute some cross correlation stats for data
allcors<-c()
for (counter in 1:dim(deer_forsnow_cl)[1])
{
  allcors[counter]<-cor(deer_forsnow_cl[counter,],snow_noNA_cl[counter,])
}
stat1_for_dat<-mean(allcors) #This is one statistic that we can compare to the same 
#statistic computed for appropriate surrogates. Since the wavelet analysis revealed
#an almost out of phase relationship between deer and snow at timescales 3-4 years, 
#we hypothesize we'll get a negative correlation, and it'll be significant compared
#to surrogates. This is also the hypothesis one would make based purely on the biological
#consideration that deer do poorly when snow depth is large. Suggests a 1-tailed test.

#%%%DAN: Lawrence pointed out that, because of the B-C normalization, the correlation of the 
#two matrices #should be the same as the mean correlation of the time series (check). Might be 
#good to point out because it means no arbitrary choice has been made.

#Now we want another statistic that corresponds to a lag that corresponds to the 
#phase diff we got from the wavelet analysis. The above was just for exactly out-of-phase,
#but the wavelet analysis was not quite out of phase. It was 0.8426*pi radians phase 
#difference, with deer leading snow depth (on average over 3-4-year timescales). So for 
#3-year timescale, this corresponds to a lag of 3*0.8426/2=1.2639 and for a 4-year timescale,
#it corresponds to a lag of 4*0.8426/2=1.6852. Neither is an integer, so we take the two 
#closest integers, so we should see if deer at time t-1 or t-2 is correlated with snow depth.
#So we hypothesize these statistics will be positive, and significantly different from 
#their values on surrogates. If we take the wavelet analyses as prior knowledge, it suggests
#a 1-tailed test. And we sort of have to take them as prior knowledge since there's no 
#reasonable scenario I can think of where one would have thought to do these tests
#without having first done the wavelet analyses.
allcorslag1<-c()
allcorslag2<-c()
for (counter in 1:dim(deer_forsnow_cl)[1])
{
  x<-deer_forsnow_cl[counter,]
  y<-snow_noNA_cl[counter,]
  allcorslag1[counter]<-cor(x[1:(length(x)-1)],y[2:length(y)])
  allcorslag2[counter]<-cor(x[1:(length(x)-2)],y[3:length(y)])
}
stat2_for_dat<-mean(allcorslag1) 
stat3_for_dat<-mean(allcorslag2)

#get the same statistics for appropriate surrogates
nsurrogs<-10000
snow_noNA_cl_s<-wsyn::surrog(dat=snow_noNA_cl,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
deer_forsnow_cl_s<-wsyn::surrog(dat=deer_forsnow_cl,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
stat1_for_s<-NA*numeric(nsurrogs)
stat2_for_s<-NA*numeric(nsurrogs)
stat3_for_s<-NA*numeric(nsurrogs)
for (scounter in 1:nsurrogs)
{
  d<-deer_forsnow_cl_s[[scounter]]
  s<-snow_noNA_cl_s[[scounter]]
  allcors<-c()
  allcorslag1<-c()
  allcorslag2<-c()
  for (lcounter in 1:(dim(d)[1]))
  {
    x<-d[lcounter,]
    y<-s[lcounter,]
    allcors[lcounter]<-cor(x,y)
    allcorslag1[lcounter]<-cor(x[1:(length(x)-1)],y[2:length(y)])
    allcorslag2[lcounter]<-cor(x[1:(length(x)-2)],y[3:length(y)])
  }
  stat1_for_s[scounter]<-mean(allcors)
  stat2_for_s[scounter]<-mean(allcorslag1)
  stat3_for_s[scounter]<-mean(allcorslag2)
}

#compare the data and surrogate statistics
p_onetailed_stat1<-1-sum(stat1_for_dat<stat1_for_s)/nsurrogs 
p_onetailed_stat1 #significant
p_onetailed_stat2<-sum(stat2_for_dat<stat2_for_s)/nsurrogs
p_onetailed_stat2 #not quite significant
p_onetailed_stat3<-sum(stat3_for_dat<stat3_for_s)/nsurrogs
p_onetailed_stat3 #significant
stat4_for_dat<-(stat2_for_dat+stat3_for_dat)/2
stat4_for_s<-(stat2_for_s+stat3_for_s)/2 #This last statistic actually 
#makes the most sense of all of them. Probably we can dispense with 
#stats 2 and 3.
p_onetailed_stat4<-sum(stat4_for_dat<stat4_for_s)/nsurrogs 
p_onetailed_stat4 #significant

#make some histograms, just to see
hist(stat1_for_s,100,xlim=range(stat1_for_s,stat1_for_dat),
     main=paste0("p=",round(p_onetailed_stat1,4)))
points(stat1_for_dat,0,type="p",col="red",pch=20,cex=2)
hist(stat2_for_s,100,xlim=range(stat2_for_s,stat2_for_dat),
     main=paste0("p=",round(p_onetailed_stat2,4)))
points(stat2_for_dat,0,type="p",col="red",pch=20,cex=2)
hist(stat3_for_s,100,xlim=range(stat3_for_s,stat3_for_dat),
     main=paste0("p=",round(p_onetailed_stat3,4)))
points(stat3_for_dat,0,type="p",col="red",pch=20,cex=2)
hist(stat4_for_s,100,xlim=range(stat4_for_s,stat4_for_dat),
     main=paste0("p=",round(p_onetailed_stat4,4)))
points(stat4_for_dat,0,type="p",col="red",pch=20,cex=2)

#***now do another similar analysis after filtering

#filter to remove variation at timescales outside the 3-4-year range
snow_noNA_cl_filt<-remper(x=snow_noNA_cl,ts=c(1/4,1/3))
deer_forsnow_cl_filt<-remper(x=deer_forsnow_cl,ts=c(1/4,1/3))

#compute some cross correlation stats for data
allcors_filt<-c()
for (counter in 1:dim(deer_forsnow_cl_filt)[1])
{
  allcors_filt[counter]<-cor(deer_forsnow_cl_filt[counter,],snow_noNA_cl_filt[counter,])
}
stat1_for_dat_filt<-mean(allcors_filt) 

allcorslag1_filt<-c()
allcorslag2_filt<-c()
for (counter in 1:dim(deer_forsnow_cl_filt)[1])
{
  x<-deer_forsnow_cl_filt[counter,]
  y<-snow_noNA_cl_filt[counter,]
  allcorslag1_filt[counter]<-cor(x[1:(length(x)-1)],y[2:length(y)])
  allcorslag2_filt[counter]<-cor(x[1:(length(x)-2)],y[3:length(y)])
}
stat2_for_dat_filt<-mean(allcorslag1_filt) 
stat3_for_dat_filt<-mean(allcorslag2_filt)

#get the same statistics for appropriate surrogates
snow_noNA_cl_filt_s<-wsyn::surrog(dat=snow_noNA_cl_filt,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
deer_forsnow_cl_filt_s<-wsyn::surrog(dat=deer_forsnow_cl_filt,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
stat1_for_s_filt<-NA*numeric(nsurrogs)
stat2_for_s_filt<-NA*numeric(nsurrogs)
stat3_for_s_filt<-NA*numeric(nsurrogs)
for (scounter in 1:nsurrogs)
{
  d<-deer_forsnow_cl_filt_s[[scounter]]
  s<-snow_noNA_cl_filt_s[[scounter]]
  allcors_filt<-c()
  allcorslag1_filt<-c()
  allcorslag2_filt<-c()
  for (lcounter in 1:(dim(d)[1]))
  {
    x<-d[lcounter,]
    y<-s[lcounter,]
    allcors_filt[lcounter]<-cor(x,y)
    allcorslag1_filt[lcounter]<-cor(x[1:(length(x)-1)],y[2:length(y)])
    allcorslag2_filt[lcounter]<-cor(x[1:(length(x)-2)],y[3:length(y)])
  }
  stat1_for_s_filt[scounter]<-mean(allcors_filt)
  stat2_for_s_filt[scounter]<-mean(allcorslag1_filt)
  stat3_for_s_filt[scounter]<-mean(allcorslag2_filt)
}

#compare the data and surrogate statistics
p_onetailed_stat1_filt<-1-sum(stat1_for_dat_filt<stat1_for_s_filt)/nsurrogs 
p_onetailed_stat1_filt #significant
p_onetailed_stat2_filt<-sum(stat2_for_dat_filt<stat2_for_s_filt)/nsurrogs
p_onetailed_stat2_filt #not significant
p_onetailed_stat3_filt<-sum(stat3_for_dat_filt<stat3_for_s_filt)/nsurrogs
p_onetailed_stat3_filt #not significant
stat4_for_dat_filt<-(stat2_for_dat_filt+stat3_for_dat_filt)/2
stat4_for_s_filt<-(stat2_for_s_filt+stat3_for_s_filt)/2 
p_onetailed_stat4_filt<-sum(stat4_for_dat_filt<stat4_for_s_filt)/nsurrogs 
p_onetailed_stat4_filt #significant 

#make some histograms, just to see
hist(stat1_for_s_filt,100,xlim=range(stat1_for_s_filt,stat1_for_dat_filt),
     main=paste0("p=",round(p_onetailed_stat1_filt,4)))
points(stat1_for_dat_filt,0,type="p",col="red",pch=20,cex=2)
hist(stat2_for_s_filt,100,xlim=range(stat2_for_s_filt,stat2_for_dat_filt),
     main=paste0("p=",round(p_onetailed_stat2_filt,4)))
points(stat2_for_dat_filt,0,type="p",col="red",pch=20,cex=2)
hist(stat3_for_s_filt,100,xlim=range(stat3_for_s_filt,stat3_for_dat_filt),
     main=paste0("p=",round(p_onetailed_stat3_filt,4)))
points(stat3_for_dat_filt,0,type="p",col="red",pch=20,cex=2)
hist(stat4_for_s_filt,100,xlim=range(stat4_for_s_filt,stat4_for_dat_filt),
     main=paste0("p=",round(p_onetailed_stat4_filt,4)))
points(stat4_for_dat_filt,0,type="p",col="red",pch=20,cex=2)
#so we again get significance in the two stats we settled on

#store the results in a variable, for reporting in the paper
pvals_deer_snow<-c(p_onetailed_stat1=p_onetailed_stat1,
                   p_onetailed_stat2=p_onetailed_stat2,
                   p_onetailed_stat3=p_onetailed_stat3,
                   p_onetailed_stat4=p_onetailed_stat4,
                   p_onetailed_stat1_filt=p_onetailed_stat1_filt,
                   p_onetailed_stat2_filt=p_onetailed_stat2_filt,
                   p_onetailed_stat3_filt=p_onetailed_stat3_filt,
                   p_onetailed_stat4_filt=p_onetailed_stat4_filt)

#***
#***Deer/MEI comparisons
#***

#***Data prep

#throw out some locations because of NAs? no need
#sum(is.na(deer))
#sum(is.na(mei))

#clean and transform
deer_cl<-wsyn::cleandat(deer,clev=5,times=deeryr)$cdat
mei_cl<-wsyn::cleandat(mei,clev=5,times=deeryr)$cdat

#***Don't do a straight, unlagged cross correlation analysis, as was done for deer and snow, 
#because we already know from the wavelet analysis that the phase shift between deer and 
#winter mei is about a quarter phase on 4-7 or 3-7 year timescales, and that would give a 
#correlation of 0. So there is no reason to expect a non-zero unlagged correlation, so why 
#test for one? Instead move straight to a lagged analysis.

#On 4-7-year timescales, the average phase between deer and winter mei was -0.6827. So for a
#4-year timescale, this corresponds to a lag of 4*0.6827/2=1.3654, and for a 7-year timescale
#it corresponds to a lag of 7*0.6827/2=2.38945. So we will use lags of 1, 2 and 3. Whereas 
#in the deer/snow analysis, we looked to see if deer at time t-l (for lag l) was correlated
#with snow depth, because the sign of the phase shift is negative in this case, we look to
#see if deer at time t+l (for lag l) was correlated with mei. One-tailed test is appropriate,
#as before.
allcorslag1<-c()
allcorslag2<-c()
allcorslag3<-c()
for (counter in 1:dim(deer_cl)[1])
{
  x<-deer_cl[counter,]
  y<-mei_cl[counter,]
  allcorslag1[counter]<-cor(y[1:(length(x)-1)],x[2:length(y)])
  allcorslag2[counter]<-cor(y[1:(length(x)-2)],x[3:length(y)])
  allcorslag3[counter]<-cor(y[1:(length(x)-3)],x[4:length(y)])
}
stat1_for_dat<-mean(allcorslag1) 
stat2_for_dat<-mean(allcorslag2)
stat3_for_dat<-mean(allcorslag3)

#get the same statistics for appropriate surrogates
mei_cl_s<-wsyn::surrog(dat=mei_cl,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
deer_cl_s<-wsyn::surrog(dat=deer_cl,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
stat1_for_s<-NA*numeric(nsurrogs)
stat2_for_s<-NA*numeric(nsurrogs)
stat3_for_s<-NA*numeric(nsurrogs)
for (scounter in 1:nsurrogs)
{
  d<-deer_cl_s[[scounter]]
  s<-mei_cl_s[[scounter]]
  allcorslag1<-c()
  allcorslag2<-c()
  allcorslag3<-c()
  for (lcounter in 1:(dim(d)[1]))
  {
    x<-d[lcounter,]
    y<-s[lcounter,]
    allcorslag1[lcounter]<-cor(y[1:(length(x)-1)],x[2:length(y)])
    allcorslag2[lcounter]<-cor(y[1:(length(x)-2)],x[3:length(y)])
    allcorslag3[lcounter]<-cor(y[1:(length(x)-3)],x[4:length(y)])
  }
  stat1_for_s[scounter]<-mean(allcorslag1)
  stat2_for_s[scounter]<-mean(allcorslag2)
  stat3_for_s[scounter]<-mean(allcorslag3)
}

#compare the data and surrogate statistics
p_onetailed_stat1<-sum(stat1_for_dat<stat1_for_s)/nsurrogs 
p_onetailed_stat1 
p_onetailed_stat2<-sum(stat2_for_dat<stat2_for_s)/nsurrogs
p_onetailed_stat2 
p_onetailed_stat3<-sum(stat3_for_dat<stat3_for_s)/nsurrogs
p_onetailed_stat3 
stat4_for_dat<-(stat1_for_dat+stat2_for_dat+stat2_for_dat)/3
stat4_for_s<-(stat1_for_s+stat2_for_s+stat2_for_s)/3 
p_onetailed_stat4<-sum(stat4_for_dat<stat4_for_s)/nsurrogs 
p_onetailed_stat4 #this last one is the most appropriate one

#make some histograms, just to see
hist(stat1_for_s,100,xlim=range(stat1_for_s,stat1_for_dat),
     main=paste0("p=",round(p_onetailed_stat1,4)))
points(stat1_for_dat,0,type="p",col="red",pch=20,cex=2)
hist(stat2_for_s,100,xlim=range(stat2_for_s,stat2_for_dat),
     main=paste0("p=",round(p_onetailed_stat2,4)))
points(stat2_for_dat,0,type="p",col="red",pch=20,cex=2)
hist(stat3_for_s,100,xlim=range(stat3_for_s,stat3_for_dat),
     main=paste0("p=",round(p_onetailed_stat3,4)))
points(stat3_for_dat,0,type="p",col="red",pch=20,cex=2)
hist(stat4_for_s,100,xlim=range(stat4_for_s,stat4_for_dat),
     main=paste0("p=",round(p_onetailed_stat4,4)))
points(stat4_for_dat,0,type="p",col="red",pch=20,cex=2)
#this is all good - we get significance in the stat we care most about (stat4)

#***now do another similar analysis after filtering

#filter to remove variation at timescales outside the 4-7-year range
mei_cl_filt<-remper(x=mei_cl,ts=c(1/7,1/4))
deer_cl_filt<-remper(x=deer_cl,ts=c(1/7,1/4))

allcorslag1<-c()
allcorslag2<-c()
allcorslag3<-c()
for (counter in 1:dim(deer_cl_filt)[1])
{
  x<-deer_cl_filt[counter,]
  y<-mei_cl_filt[counter,]
  allcorslag1[counter]<-cor(y[1:(length(x)-1)],x[2:length(y)])
  allcorslag2[counter]<-cor(y[1:(length(x)-2)],x[3:length(y)])
  allcorslag3[counter]<-cor(y[1:(length(x)-3)],x[4:length(y)])
}
stat1_for_dat_filt<-mean(allcorslag1) 
stat2_for_dat_filt<-mean(allcorslag2)
stat3_for_dat_filt<-mean(allcorslag3)

#get the same statistics for appropriate surrogates
mei_cl_filt_s<-wsyn::surrog(dat=mei_cl_filt,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
deer_cl_filt_s<-wsyn::surrog(dat=deer_cl_filt,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
stat1_for_s_filt<-NA*numeric(nsurrogs)
stat2_for_s_filt<-NA*numeric(nsurrogs)
stat3_for_s_filt<-NA*numeric(nsurrogs)
for (scounter in 1:nsurrogs)
{
  d<-deer_cl_filt_s[[scounter]]
  s<-mei_cl_filt_s[[scounter]]
  allcorslag1<-c()
  allcorslag2<-c()
  allcorslag3<-c()
  for (lcounter in 1:(dim(d)[1]))
  {
    x<-d[lcounter,]
    y<-s[lcounter,]
    allcorslag1[lcounter]<-cor(y[1:(length(x)-1)],x[2:length(y)])
    allcorslag2[lcounter]<-cor(y[1:(length(x)-2)],x[3:length(y)])
    allcorslag3[lcounter]<-cor(y[1:(length(x)-3)],x[4:length(y)])
  }
  stat1_for_s_filt[scounter]<-mean(allcorslag1)
  stat2_for_s_filt[scounter]<-mean(allcorslag2)
  stat3_for_s_filt[scounter]<-mean(allcorslag3)
}

#compare the data and surrogate statistics
p_onetailed_stat1_filt<-sum(stat1_for_dat_filt<stat1_for_s_filt)/nsurrogs 
p_onetailed_stat1_filt 
p_onetailed_stat2_filt<-sum(stat2_for_dat_filt<stat2_for_s_filt)/nsurrogs
p_onetailed_stat2_filt 
p_onetailed_stat3_filt<-sum(stat3_for_dat_filt<stat3_for_s_filt)/nsurrogs
p_onetailed_stat3_filt 
stat4_for_dat_filt<-(stat1_for_dat_filt+stat2_for_dat_filt+stat2_for_dat_filt)/3
stat4_for_s_filt<-(stat1_for_s_filt+stat2_for_s_filt+stat2_for_s_filt)/3 
p_onetailed_stat4_filt<-sum(stat4_for_dat_filt<stat4_for_s_filt)/nsurrogs 
p_onetailed_stat4_filt

#make some histograms, just to see
hist(stat1_for_s_filt,100,xlim=range(stat1_for_s_filt,stat1_for_dat_filt),
     main=paste0("p=",round(p_onetailed_stat1_filt,4)))
points(stat1_for_dat_filt,0,type="p",col="red",pch=20,cex=2)
hist(stat2_for_s_filt,100,xlim=range(stat2_for_s_filt,stat2_for_dat_filt),
     main=paste0("p=",round(p_onetailed_stat2_filt,4)))
points(stat2_for_dat_filt,0,type="p",col="red",pch=20,cex=2)
hist(stat3_for_s_filt,100,xlim=range(stat3_for_s_filt,stat3_for_dat_filt),
     main=paste0("p=",round(p_onetailed_stat3_filt,4)))
points(stat3_for_dat_filt,0,type="p",col="red",pch=20,cex=2)
hist(stat4_for_s_filt,100,xlim=range(stat4_for_s_filt,stat4_for_dat_filt),
     main=paste0("p=",round(p_onetailed_stat4_filt,4)))
points(stat4_for_dat_filt,0,type="p",col="red",pch=20,cex=2)
#this is all good - we again get significance in the stat we care most about (stat4)

#store the results in a variable, for reporting in the paper
pvals_deer_mei<-c(p_onetailed_stat1=p_onetailed_stat1,
                   p_onetailed_stat2=p_onetailed_stat2,
                   p_onetailed_stat3=p_onetailed_stat3,
                   p_onetailed_stat4=p_onetailed_stat4,
                   p_onetailed_stat1_filt=p_onetailed_stat1_filt,
                   p_onetailed_stat2_filt=p_onetailed_stat2_filt,
                   p_onetailed_stat3_filt=p_onetailed_stat3_filt,
                   p_onetailed_stat4_filt=p_onetailed_stat4_filt)

#***
#***Deer/DVC comparisons
#***

#***Data prep

#Trim the deer data so it's the same size as DVCs, for comparisons between these variables 
deer_fordvcs<-deer[,7:(dim(deer)[2])]
#dim(deer_fordvcs)
#dim(dvcs)
#sum(is.na(dvcs))
#sum(is.na(deer_fordvcs))

#clean and transform
dvcs_cl<-wsyn::cleandat(dvcs,clev=5,times=dvcyr)$cdat
deer_fordvcs_cl<-wsyn::cleandat(deer_fordvcs,clev=5,times=dvcyr)$cdat

#***Now do a cross correlation analysis

#compute some cross correlation stats for data
allcors<-c()
for (counter in 1:dim(deer_fordvcs_cl)[1])
{
  allcors[counter]<-cor(deer_fordvcs_cl[counter,],dvcs_cl[counter,])
}
stat1_for_dat<-mean(allcors) 

#no need for the lagged statistics in this case, because the phase diff the wavelet
#analysis reveals between deer and dvcs is so close to 0 (in phase)

#get the same statistic for appropriate surrogates
dvcs_cl_s<-wsyn::surrog(dat=dvcs_cl,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
deer_fordvcs_cl_s<-wsyn::surrog(dat=deer_fordvcs_cl,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
stat1_for_s<-NA*numeric(nsurrogs)
for (scounter in 1:nsurrogs)
{
  d<-deer_fordvcs_cl_s[[scounter]]
  s<-dvcs_cl_s[[scounter]]
  allcors<-c()
  for (lcounter in 1:(dim(d)[1]))
  {
    x<-d[lcounter,]
    y<-s[lcounter,]
    allcors[lcounter]<-cor(x,y)
  }
  stat1_for_s[scounter]<-mean(allcors)
}

#compare the data and surrogate statistics
p_onetailed_stat1<-sum(stat1_for_dat<stat1_for_s)/nsurrogs 
p_onetailed_stat1

#make a histogram, just to see
hist(stat1_for_s,100,xlim=range(stat1_for_s,stat1_for_dat),
     main=paste0("p=",round(p_onetailed_stat1,4)))
points(stat1_for_dat,0,type="p",col="red",pch=20,cex=2)

#***now do another similar analysis after filtering

#filter to remove variation at timescales outside the 3-7-year range
dvcs_cl_filt<-remper(x=dvcs_cl,ts=c(1/7,1/3))
deer_fordvcs_cl_filt<-remper(x=deer_fordvcs_cl,ts=c(1/7,1/3))

#compute some cross correlation stats for data
allcors<-c()
for (counter in 1:dim(deer_fordvcs_cl_filt)[1])
{
  allcors[counter]<-cor(deer_fordvcs_cl_filt[counter,],dvcs_cl_filt[counter,])
}
stat1_for_dat_filt<-mean(allcors) 

#get the same statistics for appropriate surrogates
dvcs_cl_filt_s<-wsyn::surrog(dat=dvcs_cl_filt,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
deer_fordvcs_cl_filt_s<-wsyn::surrog(dat=deer_fordvcs_cl_filt,nsurrogs=nsurrogs,surrtype="fft",syncpres=TRUE)
stat1_for_s_filt<-NA*numeric(nsurrogs)
for (scounter in 1:nsurrogs)
{
  d<-deer_fordvcs_cl_filt_s[[scounter]]
  s<-dvcs_cl_filt_s[[scounter]]
  allcors_filt<-c()
  for (lcounter in 1:(dim(d)[1]))
  {
    x<-d[lcounter,]
    y<-s[lcounter,]
    allcors_filt[lcounter]<-cor(x,y)
  }
  stat1_for_s_filt[scounter]<-mean(allcors_filt)
}

#compare the data and surrogate statistics
p_onetailed_stat1_filt<-sum(stat1_for_dat_filt<stat1_for_s_filt)/nsurrogs 
p_onetailed_stat1_filt

#make a histogram, just to see
hist(stat1_for_s_filt,100,xlim=range(stat1_for_s_filt,stat1_for_dat_filt),
     main=paste0("p=",round(p_onetailed_stat1_filt,4)))
points(stat1_for_dat_filt,0,type="p",col="red",pch=20,cex=2)

#store the results in a variable, for reporting in the paper
pvals_deer_dvcs<-c(p_onetailed_stat1=p_onetailed_stat1,
                   p_onetailed_stat1_filt=p_onetailed_stat1_filt)
