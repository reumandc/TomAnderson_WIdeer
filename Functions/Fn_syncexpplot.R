#Plot Expected synchrony based on output of wmrsig or modelsyncexp
#resp.array- wavelet transforms of the response variable in location x time x timescale format
#exp.sync- expected synchrony from either a model or modelsyncexp
#times = timestep vector
#time scales = timescale vector

syncexpplot<-function(resp.wmf,exp.sync,times,timescales,xlab,ylab,smallplot=NULL){
  ylocs <- pretty(timescales, n = 8)
  xlocs <- pretty(times, n = 8)
  jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
                 "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
  colorfill <- colorRampPalette(jetcolors)

  image(x = times, y = log2(timescales), z = Mod(exp.sync), xlab = "Year", 
        ylab = ylab, axes = F, col = colorfill(100),las=1)
  axis(1, at = xlocs, labels = xlocs)
  axis(2, at = log2(ylocs), labels = ylocs, las = 1)
  contour(x = times, y = log2(timescales),z=Mod(resp.wmf),add=T,frame=F,las=1,lwd=2)
  image.plot(x = times, y = timescales, z = Mod(exp.sync),
             legend.only = T, 100, col = colorfill(100),las=1,smallplot=smallplot)
}