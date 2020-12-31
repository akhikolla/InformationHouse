dgraffelmanweir.bi <- function(x,y) {
  nm <- sum(x)
  nf <- sum(y)
  n <- nm + nf
  mA <- x["A"]
  mB <- x["B"]
  nfAA <- y["AA"]
  nfAB <- y["AB"]
  nfBB <- y["BB"]
  fA <- 2*nfAA + nfAB
  fB <- 2*nfBB + nfAB
  nA <- mA + fA
  nB <- mB + fB
  nt <- nm + 2 * nf
  log2 <- log(2)
  K <- lgamma(nA + 1) + lgamma(nB + 1) + lgamma(nf + 1) + lgamma(nm + 
                                                                   1) - lgamma(mA + 1) - lgamma(mB + 1) - lgamma(nt + 1)
  quotient <- nfAB * log2 - lgamma(nfAA + 1) - lgamma(nfAB + 
                                                        1) - lgamma(nfBB + 1)
  prob <- exp(K + quotient)
  return(prob)
}
