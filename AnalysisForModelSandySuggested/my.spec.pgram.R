#Exactly the same as the spec.pgram in the stats package in R, but results
#also include the full spectral matrix which is computed by spec.pgram but 
#not returned by the original function. The class is attribute is still
#set to "spec", in spite of this additional component of the output, so
#it is possible calling "spec" methods on the output of this function will
#work, and it is possible that will fail.

my.spec.pgram<-function (x, spans = NULL, kernel = NULL, taper = 0.1, pad = 0, 
          fast = TRUE, demean = FALSE, detrend = TRUE, plot = TRUE, 
          na.action = na.fail, ...) 
{
  series <- deparse(substitute(x))
  x <- na.action(as.ts(x))
  xfreq <- frequency(x)
  x <- as.matrix(x)
  N <- N0 <- nrow(x)
  nser <- ncol(x)
  if (!is.null(spans)) 
    kernel <- {
      if (is.tskernel(spans)) 
        spans
      else kernel("modified.daniell", spans%/%2)
    }
  if (!is.null(kernel) && !is.tskernel(kernel)) 
    stop("must specify 'spans' or a valid kernel")
  if (detrend) {
    t <- 1L:N - (N + 1)/2
    sumt2 <- N * (N^2 - 1)/12
    for (i in 1L:ncol(x)) x[, i] <- x[, i] - mean(x[, i]) - 
      sum(x[, i] * t) * t/sumt2
  }
  else if (demean) {
    x <- sweep(x, 2, colMeans(x), check.margin = FALSE)
  }
  x <- spec.taper(x, taper)
  u2 <- (1 - (5/8) * taper * 2)
  u4 <- (1 - (93/128) * taper * 2)
  if (pad > 0) {
    x <- rbind(x, matrix(0, nrow = N * pad, ncol = ncol(x)))
    N <- nrow(x)
  }
  NewN <- if (fast) 
    nextn(N)
  else N
  x <- rbind(x, matrix(0, nrow = (NewN - N), ncol = ncol(x)))
  N <- nrow(x)
  Nspec <- floor(N/2)
  freq <- seq.int(from = xfreq/N, by = xfreq/N, length.out = Nspec)
  xfft <- mvfft(x)
  pgram <- array(NA, dim = c(N, ncol(x), ncol(x)))
  for (i in 1L:ncol(x)) {
    for (j in 1L:ncol(x)) {
      pgram[, i, j] <- xfft[, i] * Conj(xfft[, j])/(N0 * 
                                                      xfreq)
      pgram[1, i, j] <- 0.5 * (pgram[2, i, j] + pgram[N, 
                                                      i, j])
    }
  }
  if (!is.null(kernel)) {
    for (i in 1L:ncol(x)) for (j in 1L:ncol(x)) pgram[, i, 
                                                      j] <- kernapply(pgram[, i, j], kernel, circular = TRUE)
    df <- df.kernel(kernel)
    bandwidth <- bandwidth.kernel(kernel)
  }
  else {
    df <- 2
    bandwidth <- sqrt(1/12)
  }
  df <- df/(u4/u2^2)
  df <- df * (N0/N)
  bandwidth <- bandwidth * xfreq/N
  pgram <- pgram[2:(Nspec + 1), , , drop = FALSE]
  spec <- matrix(NA, nrow = Nspec, ncol = nser)
  for (i in 1L:nser) spec[, i] <- Re(pgram[1L:Nspec, i, i])
  if (nser == 1) {
    coh <- phase <- NULL
  }
  else {
    coh <- phase <- matrix(NA, nrow = Nspec, ncol = nser * 
                             (nser - 1)/2)
    for (i in 1L:(nser - 1)) {
      for (j in (i + 1):nser) {
        coh[, i + (j - 1) * (j - 2)/2] <- Mod(pgram[, 
                                                    i, j])^2/(spec[, i] * spec[, j])
        phase[, i + (j - 1) * (j - 2)/2] <- Arg(pgram[, 
                                                      i, j])
      }
    }
  }
  for (i in 1L:nser) spec[, i] <- spec[, i]/u2
  spec <- drop(spec)
  spg.out <- list(freq = freq, spec = spec, pgram=pgram, coh = coh, phase = phase, 
                  kernel = kernel, df = df, bandwidth = bandwidth, n.used = N, 
                  orig.n = N0, series = series, snames = colnames(x), method = ifelse(!is.null(kernel), 
                                                                                      "Smoothed Periodogram", "Raw Periodogram"), taper = taper, 
                  pad = pad, detrend = detrend, demean = demean)
  class(spg.out) <- "spec"
  if (plot) {
    plot(spg.out, ...)
    return(invisible(spg.out))
  }
  else return(spg.out)
}