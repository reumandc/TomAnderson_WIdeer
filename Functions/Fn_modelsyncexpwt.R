
#Expected synchrony based on a wmrsig model
#resp.array- wavelet transforms of the response variable in location x time x timescale format
#model- need to extract model1 from the wmrsig output
#tsrange- valus over which to calculate average synchrony
#plot the predicted and actual mean fields
#times = timestep vector
#time scales = timescale vector

syncexpwt<-function(resp.array,model,times,timescales,tsrange,plot=T){
  source("Functions/Fn_wmfwt.R")
  source("Functions/Fn_swcohwt.R")
  
  #convert model from wide 2d format to 3d array
  locs<-seq(1,dim(model)[2],dim(resp.array)[2])
  model.array<-array(NA,c(dim(resp.array)))
  for(i in 1:length(locs)){
    model.array[i,,]<-t(model[,c(locs[i]:(locs[i]+(dim(resp.array)[2]-1)))])
  }
  
  model.wmf<-wmfwt(model.array) #could replace with your own function to make a wmf out of transforms
  resp.wmf<-wmfwt(resp.array) #could replace with your own function to make a wmf out of transforms 
  model.spcoh<-swcohwt(bio.array = resp.array,model.array) #could replace with your own function to take spatial coherence of wavelet transforms
  spcoh.mat<-matrix(model.spcoh,nrow=dim(model.wmf)[1],ncol=dim(model.wmf)[2],byrow=T)
  
  exp.sync <- model.wmf * spcoh.mat#calculate expected synchrony
  ta.sync<-apply(Mod(resp.wmf)^2,2,mean,na.rm=T)#calculate time-averaged synchrony
  ta.predsync<-apply(Mod(exp.sync)^2,2,mean,na.rm=T)#calculate time-averaged predicted synchrony
  ta.diff<-apply(Mod(resp.wmf-exp.sync)^2,2,mean,na.rm=T)#calculate the difference in the expected synchrony from the response variable
  
  #calculate the average amount of synchrony explained over the timescales specified in range
  ave.exp<-mean(ta.predsync[timescales>tsrange[1] & timescales<tsrange[2]])/mean(ta.sync[timescales>tsrange[1] & timescales<tsrange[2]],na.rm=T) 
  ave.exp
  
  if(plot){
    ylocs <- pretty(timescales, n = 8)
    xlocs <- pretty(times, n = 8)
    jetcolors <- c("#00007F", "blue", "#007FFF", "cyan", 
                   "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
    colorfill <- colorRampPalette(jetcolors)
    par(mfrow=c(1,2))
    image(x = times, y = log2(timescales), z = Mod(exp.sync), xlab = "Time", 
          ylab = "Timescale", axes = F, col = colorfill(100), 
          main = "Predicted")
    axis(1, at = xlocs, labels = xlocs)
    axis(2, at = log2(ylocs), labels = ylocs)
    image.plot(x = times, y = timescales, z = Mod(exp.sync),
               legend.only = T, 100, col = colorfill(100))
    
    image(x = times, y = log2(timescales), z = Mod(resp.wmf), xlab = "Time", 
          ylab = "", axes = F, col = colorfill(100), 
          main = "Real")
    axis(1, at = xlocs, labels = xlocs)
    axis(2, at = log2(ylocs), labels = ylocs)
    image.plot(x = times, y = timescales, z = Mod(resp.wmf),
               legend.only = T, 100, col = colorfill(100))
    par(mfrow=c(1,1))
  }
  return(list(ave.syncexp=ave.exp,exp.sync=exp.sync))
}
# 
# abun.array<-warray(norm.dat[[1]],times=1981:2016)
# times<-abun.array$times
# timescales<-abun.array$timescales
# syncexp(resp.array=abun.array$wave.array,model = abunmod13a$predwt,times = times,timescales=timescales,tsrange=c(3,7))
