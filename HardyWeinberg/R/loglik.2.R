loglik.2 <- function(pa,z,f) {
  # we have to do some work to avoid negative arguments to the log function which can arise due to rounding if there is a zero homozygote count.
  paa <- pa^2 + pa*(1-pa)*f
  pab <- 2*pa*(1-pa)*(1-f)
  pbb <- (1-pa)^2 + pa*(1-pa)*f
  pvec <- c(paa,pab,pbb)
  ind <- !(z==0)
  logvec <- log(pvec[ind])
  loglik0 <- sum(z[ind]*logvec)
  nparam0 <- 2
  res <- c(loglik0,nparam0)
  return(res)
}
