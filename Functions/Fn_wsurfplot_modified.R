wsurfplotTLA<-function (dat, times, type, colorbar=F, zlims = NULL, colorfill = NULL,tsrange=NULL, 
          neat = T,siglevel=NULL,title = NULL, xlab=NULL,ylab=NULL,xtcklab=NULL,ytcklab=NULL,xlim=NULL,
          ...) 
{
  library(fields)
  if (is.null(type)) {
    if (is.matrix(dat)) {
      stop("Specify desired plot type: wmf or wpmf")
    }
    if (is.vector(dat)) {
      warning("Transforms of a single location detected; plotting wavelet power")
    }
  }
  if(is.null(xlab)){
    xlab<-"Time"
  }
  if(is.null(ylab)){
    ylab<-"Timescale"
  }
  if (type == "power") {
    dat.wav <- wt(dat, times = times)
    wav <- Mod(dat.wav$wave) #wavelet amplitude (could be power if we squared this)
    denom<-sqrt(apply(wav,2,mean,na.rm=T))
    wav<-sweep(wav,2,denom,`/`)
    timescales <- dat.wav$timescales
  }
  if (type == "wmf") {
    dat.wav <- wmf(dat, times = times)
    wav <- Mod(dat.wav$wmf)
    timescales <- dat.wav$timescales
  }
  if (type == "wpmf") {
    dat.wav <- wpmf(dat, times = times)
    wav <- Mod(dat.wav$wpmf)
    timescales <- dat.wav$timescales
  }
  if (is.null(zlims)) {
    zlims <- range(Mod(wav), na.rm = T)
  }
  if (neat) {
    wav <- wav[, !is.na(colMeans(wav, na.rm = T))]
    timescales <- timescales[1:ncol(wav)]
  }
  if (is.null(colorfill)) {
    jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
    colorfill <- colorRampPalette(jetcolors)
  }
  ylocs <- pretty(timescales, n = 8)
  xlocs <- pretty(times, n = 8)
  if(is.null(xtcklab)){
    xtcklab <- xlocs
  }
  if(is.null(ytcklab)){
    ytcklab <- ylocs
  }
  if(is.null(xlim)){
    xlim <- range(times,na.rm=T)
  }
  image(x = times, y = log2(timescales), z = wav, xlab = xlab, 
        zlim = zlims, ylab = ylab, axes = F, col = colorfill(100), 
        main = title, xlim=xlim,...)
  if(!is.null(tsrange)){
    abline(h=c(log2(tsrange[1]),log2(tsrange[2])),lty=c(2,2))
  }
  axis(1, at = xlocs, labels = xtcklab)
  axis(2, at = log2(ylocs), labels = ytcklab, las = 1)
  if(type=="wpmf"){
    surrogs=100000
    sites<-nrow(dat)
    surrog.mat<-matrix(runif(sites*surrogs,min = 0,max=2*pi),sites,surrogs)
    phase.mat<-matrix(NA,sites,surrogs)
    phase.mat<-matrix(complex(length.out = sites*surrogs,modulus=1,argument=surrog.mat),nrow =sites,ncol=surrogs)
    mean.phase<-colMeans(phase.mat)
    mod.mean.phase<-Mod(mean.phase)
    if(is.null(siglevel)){
      siglevel<-quantile(mod.mean.phase,probs=0.95)
    }
    else{
      siglevel<-quantile(mod.mean.phase,probs=siglevel)
    }
    for(i in 1:length(siglevel)){
      par(new=T)
      if(is.null(xlim)){
        contour(wav,levels=siglevel[i],drawlabels=F,lwd=2,lty=i,xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5), las = 1,frame=F)
      }
      else{
        contour(wav,levels=siglevel[i],drawlabels=F,lwd=2,lty=i,xaxs="i",xaxt="n",yaxt="n",xaxp=c(0,1,5), las = 1,frame=F,
                xlim=c(as.numeric(as.factor(times))[(times==xlim[1])]/length(times),as.numeric(as.factor(times))[(times==xlim[2])]/length(times)))
      }
    }
  }
  #Add color bar- must do last otherwise chaos ensues
  if(colorbar){
    if(is.null(xlim)){
      image.plot(x = times, y = timescales, z = wav, zlim = zlims, 
                 legend.only = T, 100, col = colorfill(100),...)
    }
    else{
      image.plot(x = times, y = timescales, z = wav, zlim = zlims, 
                 legend.only = T, 100, col = colorfill(100),
                 xlim=c(as.numeric(as.factor(times))[(times==xlim[1])]/length(times),as.numeric(as.factor(times))[(times==xlim[2])]/length(times)))
    }

  }

}