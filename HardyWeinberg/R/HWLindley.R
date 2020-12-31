HWLindley <-
function(alphaseq=seq(-3,3,by=0.01),x) {
  N <- sum(x)
  nAA <- x[1]
  nAB <- x[2]
  nBB <- x[3]
  num <- lgamma(N) 
  den <- lgamma(nAA) + lgamma(nAB) + lgamma(nBB) + (nAA + nBB -1)*log(2)
  fac <- exp(num - den)
  fvalue <- numeric(length(alphaseq))
  for(i in 1:length(alphaseq)) {
    alpha <- alphaseq[i]
    z <- integrate(auxlindley,-Inf,Inf,alpha,nAA,nBB,N)
    value <- z$value
    fvalue[i] <-  fac*value
  }
  return(fvalue)
}
