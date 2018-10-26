deer_wpmfplot<-function(object,zlims=NULL,neat=TRUE,xlab,ylab,colorfill=NULL,sigthresh=0.95,colorbar=TRUE,title=NULL,filename=NA,...)
{
  wav<-Mod(get_values(object))
  times<-get_times(object)
  timescales<-get_timescales(object)
  signif<-get_signif(object)

  if (any(sigthresh>=1 | sigthresh<=0))
  {
    stop("Error in plotmag.wpmf: inappropriate value for sigthresh")
  }
  if(is.null(zlims)){
    zlims<-range(wav,na.rm=T)
  }else
  {
    rg<-range(wav,na.rm=T)
    if (rg[1]<zlims[1] || rg[2]>zlims[2])
    {
      stop("Error in plotmag.wpmf: zlims must encompass the z axis range of what is being plotted")
    }
  }
  if(neat){
    inds<-which(!is.na(colMeans(wav,na.rm=T)))
    wav<-wav[,inds]
    timescales<-timescales[inds]
    if (!is.na(signif) && (signif[[1]] %in% c("fft","aaft")))
    {
      signif[[3]]<-signif[[3]][,inds]
    }
  }
  if(is.null(colorfill)){
    jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
    colorfill<-grDevices::colorRampPalette(jetcolors)
  }
  ylocs<-pretty(timescales,n=8)
  xlocs<-pretty(times,n=8)
  
  if (!is.na(filename))
  {
    grDevices::pdf(paste0(filename,".pdf"))
  }
  if (!colorbar)
  {
    graphics::image(x=times,y=log2(timescales),z=wav,xlab=xlab,zlim=zlims,
                    ylab=ylab,axes=F,col=colorfill(100),main=title,...)
    graphics::axis(1,at = xlocs,labels=xlocs)
    graphics::axis(2,at = log2(ylocs),labels = ylocs)
  }else
  {
    fields::image.plot(x=times,y=log2(timescales),z=wav,xlab=xlab,zlim=zlims,
                       ylab=ylab,axes=F,col=colorfill(100),main=title,...)
    graphics::axis(1,at = xlocs,labels=xlocs)
    graphics::axis(2,at = log2(ylocs),labels = ylocs)
  }
  if (!all(is.na(signif)))
  {
    graphics::par(new=T)
    if (signif[[1]]=="quick")
    {
      q<-stats::quantile(signif[[2]],sigthresh)
      graphics::contour(x=times,y=log2(timescales),z=wav,levels=q,drawlabels=F,lwd=2,
                        xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5),las = 1,frame=F)
    }
    if (signif[[1]] %in% c("fft","aaft"))
    {
      
      graphics::contour(x=times,y=log2(timescales),z=signif[[3]],levels=sigthresh,
                        drawlabels=F,lwd=2,xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5),
                        las = 1,frame=F)
    }
  }
  if (!is.na(filename))
  {
    grDevices::dev.off()
  }
  return(NULL) 
}

deer_wmfplot<-function(object,zlims=NULL,neat=TRUE,colorfill=NULL,power=F,xlab,ylab,colorbar=TRUE,title=NULL,filename=NA,...)
{
  times<-get_times(object)
  timescales<-get_timescales(object)
  if(power){
    wav<-wavnorm(object$values)
  }
  else{
    wav<-Mod(get_values(object))
  }
  if(is.null(zlims)){
    zlims<-range(wav,na.rm=T)
  }else
  {
    rg<-range(wav,na.rm=T)
    if (rg[1]<zlims[1] || rg[2]>zlims[2])
    {
      stop("Error in plotmag.tts: zlims must encompass the z axis range of what is being plotted")
    }
  }
  if(neat){
    inds<-which(!is.na(colMeans(wav,na.rm=T)))
    wav<-wav[,inds]
    timescales<-timescales[inds]
  }
  if(is.null(colorfill)){
    jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
    colorfill<-grDevices::colorRampPalette(jetcolors)
  }
  ylocs<-pretty(timescales,n=8)
  xlocs<-pretty(times,n=8)
  
  if (!is.na(filename))
  {
    grDevices::pdf(paste0(filename,".pdf"))
  }
  if (!colorbar)
  {
    graphics::image(x=times,y=log2(timescales),z=wav,xlab=xlab,zlim=zlims,
                    ylab=ylab,axes=F,col=colorfill(100),main=title,...)
    graphics::axis(1,at = xlocs,labels=xlocs)
    graphics::axis(2,at = log2(ylocs),labels = ylocs)
  }else
  {
    fields::image.plot(x=times,y=log2(timescales),z=wav,xlab=xlab,zlim=zlims,
                       ylab=ylab,axes=F,col=colorfill(100),main=title,...)
    graphics::axis(1,at = xlocs,labels=xlocs)
    graphics::axis(2,at = log2(ylocs),labels = ylocs)
  }
  if (!is.na(filename))
  {
    grDevices::dev.off()
  }
  return(NULL)
}

syncexpplot<-function(resp.wmf,exp.sync,times,timescales,xlab,ylab,smallplot=NULL){
  ylocs <- pretty(timescales, n = 8)
  xlocs <- pretty(times, n = 8)
  jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
  colorfill <- colorRampPalette(jetcolors)
  image(x = times, y = log2(timescales), z = Mod(exp.sync), xlab = xlab, 
        ylab = ylab, axes = F, col = colorfill(100),las=1)
  axis(1, at = xlocs, labels = xlocs)
  axis(2, at = log2(ylocs), labels = ylocs, las = 1)
  contour(x = times, y = log2(timescales),z=Mod(resp.wmf),add=T,frame=F,las=1,lwd=2)
  image.plot(x = times, y = timescales, z = Mod(exp.sync),
             legend.only = T, 100, col = colorfill(100),las=1)
}