rankplot<-function(rank,timescales,surrogates,sig=0.95,ylim=NULL){
  if(is.null(ylim)){
    ylim<-c(min(rank),surrogates)
  }
  plot(timescales,rank,ylim=ylim,xlab="Timescale",ylab="Rank",type="l")
  abline(h=sig*surrogates,col="red",lty=2)
}
