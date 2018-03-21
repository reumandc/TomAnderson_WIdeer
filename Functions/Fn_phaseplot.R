#Function for plotting phases of significant spatial coherences using fast algorithm
#
#Args
#spatcoh  --a spatial coherence vector, made up of complex numbers
#spatcoh.p --the results of the significance test for spatial coherence
#cutoff    --numerical value for the breakpoint for cutoffs in teh spatial coherence test
#pcut     --what alpha is in the spatial coherence significance test
#timescales --vector of timescales reported from fast spatial coherence algorithm
#
#Output
##plots of all significant phase relationships for both short and long timescales

phaseplot<-function(spatcoh,timescales,type,tsrange,xlab=xlab,ylab=ylab,showphase=T,colrect=T){
  if(type=="radian"){
    phasemean<-sprintf("%.3f",round(Arg(mean(spatcoh[timescales>=tsrange[1] & timescales <=tsrange[2]]/Mod(spatcoh[timescales>=tsrange[1] & timescales <=tsrange[2]]),na.rm=T)),digits=3)) 
    plot(Arg(spatcoh)~timescales,ylim=c(-pi,pi),xlab=xlab,ylab=ylab)
    abline(v=c(tsrange[1],tsrange[2]),lty=2)
    if(colrect){
      rect(xleft=tsrange[1],ybottom = -pi,xright = tsrange[2],ytop = pi,col=rgb(0,0,1,alpha=.5),density=25,lwd=2)
    }
    if(showphase){
      mtext(line=0,adj=1,cex=0.75,bquote(bar(theta)==.(phasemean)))
    }
  }
  if(type=="pi"){
    phasemean<-sprintf("%.3f",round(Arg(mean(spatcoh[timescales>=tsrange[1] & timescales <=tsrange[2]]/Mod(spatcoh[timescales>=tsrange[1] & timescales <=tsrange[2]]),na.rm=T))/3.14,digits=3)) 
    spatcohpi<-Arg(spatcoh)/3.14
    plot(spatcohpi~timescales,ylim=c(-1,1),xlab=xlab,ylab=ylab,pch=19)
    abline(v=c(tsrange[1],tsrange[2]),lty=2)
    if(colrect){
      rect(xleft=tsrange[1],ybottom = -pi,xright = tsrange[2],ytop = pi,col=rgb(0,0,1,alpha=.5),density=25,lwd=2)
    }
    if(showphase){
      mtext(line=0,adj=1,cex=0.75,bquote(bar(theta)==.(phasemean)))
    }
  }
}
