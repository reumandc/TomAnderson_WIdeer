#Function for normalizing wavelet transforms of a time series
#
#Args
#wav  --wavelet transforms from a call to wt() function from the 'wsyn' package. Should be
#       in time X timescale format
#
#Output
##a matrix containing power normalized wavelet tranforms (i.e., wavelet transforms divided
##by the square root of the time-averaged squared wavelet transform magnitudes)
wavnorm<-function(wav){
  modwav <- Mod(wav)
  denom<-sqrt(apply(modwav^2,2,mean,na.rm=T))
  normwav<-sweep(modwav,2,denom,`/`)
}