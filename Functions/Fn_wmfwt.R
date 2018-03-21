wmfwt<-function(wav.array){
  conj.array <- array(NA, dim = dim(wav.array))
  for (i in 1:dim(conj.array)[1]) {
    conj.array[i, , ] <- wav.array[i, , ] * Conj(wav.array[i, 
                                                           , ])
    conj.array <- Re(conj.array)
  }
  norm.denom <- sqrt(apply(conj.array, 3, mean, na.rm = T))
  norm.array = array(NA, dim = dim(wav.array))
  for (i in 1:dim(wav.array)[1]) {
    norm.array[i, , ] <- t(t(wav.array[i, , ])/norm.denom)
  }
  wavelet.mean.field <- apply(norm.array, c(2, 3), mean, na.rm = T)
  return(wmf = wavelet.mean.field)
}

