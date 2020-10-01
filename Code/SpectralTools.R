#Some tools for the alternative Fourier analysis

#A raw, unsmoothed periodigram. Detrending is done by this function. 
#
#Args
#x      A time series as a vector
#
myspecraw<-function(x)
{
  tforx<-1:length(x)
  x<-residuals(lm(x~tforx))
  h<-fft(x)
  h<-Re(h*Conj(h))
  freq<-(0:length(x))/length(x)
  h<-h[-1]
  freq<-freq[-1]
  h<-h[freq<=0.5]
  freq<-freq[freq<=0.5]
  return(list(freq=freq,spec=h))
}

#The power spectrum, Brillinger's consistent estimator (5.6 of Brillinger's 2001 
#book). The only difference from what is described there is frequency is here in 
#units of cycles per sampling interval in the output here, and was in radians per 
#sampling interval in Brillinger. Detrending is optionally done by this function. 
#The function actually returns the log scale power spectrum.
#
#Args
#x      A time series as a vector
#BiasVariance     For adjusting the bias-variance tradeoff which comes from the 
#                   degree of smoothing selected
#
myspecbrill<-function(x,detrend=TRUE,BiasVariance=0.5)
{
  Tx<-length(x)
  
  #detrend, if desired
  if (detrend==TRUE)
  {
    tforx<-1:Tx
    x<-residuals(lm(x~tforx))
  }
  
  #get the raw periodogram
  fftx<-fft(x)
  I<-(Mod(fftx))^2/(2*pi*Tx)
  I[1]<-0 #Set zero frequency to 0. Should be zero anyway, to within rounding error, because of the detrending
  freq<-(0:(Tx-1))/Tx
  freq<-2*pi*freq #to make frequencies be in units of radians per sampling interval
  
  #now do the smoothing that makes the estimator
  BTx<-BiasVariance/sqrt(Tx) #adjust BT to adjust the bias-variance tradeoff
  TxBTx<-Tx*BTx
  xforW<-(0:floor(TxBTx))/TxBTx
  WT<-(15/(16*2*pi))*((xforW-1)^2)*((xforW+1)^2)
  intWsquare<-(15/(32*pi))^2*4*pi*(1/9-4/7+6/5-4/3+1)
  spec<-WT[1]*I
  lenI<-length(I)
  for (AbsShift in 1:(length(WT)-1))
  {
    temp<-WT[AbsShift+1]*I
    #spec<-spec + temp([(AbsShift+1):end 1:AbsShift],:) + temp([(end-AbsShift+1):end 1:(end-AbsShift)],:);
    spec<-spec+temp[c((AbsShift+1):lenI,1:AbsShift)]+temp[c((lenI-AbsShift+1):lenI,1:(lenI-AbsShift))]
  }
  
  #remove the 0 frequency, and change the units back to cycles per sampling interval,
  #and cut the redundant part of the spectrum
  freq<-freq[-1]
  spec<-spec[-1]
  freq<-freq/(2*pi)
  spec<-spec[freq<=0.5]
  freq<-freq[freq<=0.5]
  
  #put in some normalization factors
  spec<-spec*2*pi/TxBTx
  
  #now get confidence intervals
  p<-0.95 #for 95% confidence intervals
  conf<-qnorm((1+p)/2,mean=0,sd=1)*0.4343*sqrt(2*pi*intWsquare/(TxBTx));
  
  return(list(freq=freq,log10spec=log10(spec),conf=conf))
}

#Takes the fft of each row of the matrix x. Just a convenience function.
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

#The spectral matrix, Brillinger's consistent estimator (7.4 of Brillinger's 2001 
#book). The only difference from what is described there is frequency is here in 
#units of cycles per sampling interval in the output, and was in radians per sampling 
#interval in Brillinger. Linear detrending (and de-meaning) of the individual 
#component time series is optionally done.
#
#Args
#x      A vector-valued time series as a matrix, components by time (so time runs 
#         along the rows). 
#BiasVariance     For adjusting the bias-variance tradeoff which comes from the 
#                   degree of smoothing selected
#
myspecmatbrill<-function(x,detrend=TRUE,BiasVariance=0.5)
{
  Tx<-dim(x)[2] #length of time series
  N<-dim(x)[1] #number of components
  
  #detrend
  if (detrend)
  {
    tforx<-1:Tx
    for (counter in 1:N)
    {
      y<-x[counter,]
      x[counter,]<-residuals(lm(y~tforx))  
    }
  }
  
  #get the raw periodogram
  fftx<-myfft(x)
  I<-array(complex(real=0,imaginary=0),dim=c(N,N,Tx))
  for (a in 1:N)
  {
    for (b in 1:N)
    {
      I[a,b,]<-fftx[a,]*Conj(fftx[b,])
    }
  }
  I<-I/(2*pi*Tx)
  freq<-(0:(Tx-1))/Tx
  freq<-2*pi*freq #to make frequencies be in units of radians per sampling interval, for now
  
  #now do the smoothing that makes the estimator
  BTx<-BiasVariance/sqrt(Tx) #adjust BT to adjust the bias-variance tradeoff
  TxBTx<-Tx*BTx
  xforW<-(0:floor(TxBTx))/TxBTx
  WT<-(15/(16*2*pi))*((xforW-1)^2)*((xforW+1)^2)
  intWsquare<-(15/(32*pi))^2*4*pi*(1/9-4/7+6/5-4/3+1)
  specmat<-WT[1]*I
  lenI<-dim(I)[3]
  for (AbsShift in 1:(length(WT)-1))
  {
    temp<-WT[AbsShift+1]*I
    specmat<-specmat+temp[,,c((AbsShift+1):lenI,1:AbsShift)]+temp[,,c((lenI-AbsShift+1):lenI,1:(lenI-AbsShift))]
  }
  
  #remove the 0 frequency, and change the units back to radians per sampling interval,
  #and cut the redundant part 
  freq<-freq[-1]
  specmat<-specmat[,,-1,drop=FALSE]
  freq<-freq/(2*pi)
  specmat<-specmat[,,freq<=0.5]
  freq<-freq[freq<=0.5]
  
  #put in some normalization factors
  specmat<-specmat*2*pi/TxBTx
  
  return(list(freq=freq,specmat=specmat))
}
