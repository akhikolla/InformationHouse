loglik.M11 <-
function(pa,z) {
  paa <- pa*pa
  pab <- 2*pa*(1-pa)
  pbb <- (1-pa)*(1-pa)
  pvec <- c(paa,pab,pbb)
  ind <- !(z==0)
  logvec <- log(pvec[ind])
  loglik0 <- sum(z[ind]*logvec)
  nparam0 <- 1
  res <- c(loglik0,nparam0)
  return(res)
}
