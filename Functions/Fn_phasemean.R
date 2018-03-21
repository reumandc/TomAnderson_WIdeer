
#Function for getting mean phases of significant spatial coherences
#
#Args
#spatcoh = a spatial coherence vector, made up of complex numbers
#timescales = vector of timescales reported from fast spatial coherence algorithm
#period = either long or short to return the average phase relationship over that timescale
#Output
##plots of all significant phase relationships for both short and long timescales

###CURRENTLY EMBEDDED WITHIN PLOTTING FUNCTION
phasemean<-function(spatcoh,timescales,tsrange,cutoff=NULL){
  if(is.null(cutoff)){
    tmp<-mean(spatcoh[timescales>=tsrange[1] & timescales <=tsrange[2]]/Mod(spatcoh[timescales>=tsrange[1] & timescales <=tsrange[2]]),na.rm=T)
    if(Mod(tmp)<0.2){
      warning("Average phase may be misleading")
    }
    phase<-round(Arg(tmp),digits=4)  
}
  if(!is.null(cutoff)){
    tmp1<-mean(spatcoh[timescales<=cutoff]/Mod(spatcoh[timescales<=cutoff]),na.rm=T)
    tmp2<-mean(spatcoh[timescales>=cutoff]/Mod(spatcoh[timescales>=cutoff]),na.rm=T)
    if(Mod(tmp1)<0.2||Mod(tmp2)<0.2){
      warning("Average phase may be misleading")
    }
    phase1<-round(Arg(tmp1),digits=4)
    phase2<-round(Arg(tmp2),digits=4)
    phase<-c(phase1,phase2)
  }
  return(phase)
}
