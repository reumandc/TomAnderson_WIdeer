#Do spatial coherence using wavelet transforms
#bio.array and env.arry are wavelet transforms in location x time x timescale format of a response and predictor variable

swcohwt <- function(bio.array, env.array) {
  if (sum(dim(bio.array) - dim(env.array)) != 0) {
    stop("Array dimensions do not match")
  }
  normalize <- function(in.array) {
    conj.array <- array(NA, dim = dim(in.array))
    for (i in 1:dim(conj.array)[1]) {
      conj.array[i, , ] <- in.array[i, , ] * Conj(in.array[i, 
                                                           , ])
      conj.array <- Re(conj.array)
    }
    norm.denom <- sqrt(apply(conj.array, 3, mean, na.rm = T))
    norm.array = array(NA, dim = dim(in.array))
    for (i in 1:dim(in.array)[1]) {
      norm.array[i, , ] <- t(t(in.array[i, , ])/norm.denom)
    }
    return(norm.array)
  }
  bio.norm <- normalize(bio.array)
  env.norm <- normalize(env.array)
  bXe <- array(NA, dim = dim(bio.array))
  for (i in 1:dim(bXe)[1]) {
    bXe[i, , ] <- bio.norm[i, , ] * Conj(env.norm[i, 
                                                  , ])
  }
  coh <- apply(bXe, 3, mean, na.rm = T)
  return(coh)
}