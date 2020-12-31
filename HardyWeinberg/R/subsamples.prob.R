subsamples.prob <- function (nA, nB, nm, nf, nt, mA, mB, nAf, psamp, eps) 
{
  #
  # We now work with a relative epsilon
  #
  log2 <- log(2)
  nBf <- 2*nf - nAf
  if ((nAf%%2) == 0) nfAB <- seq(0, nAf, 2) else nfAB <- seq(1, nAf, 2)
  nfAA <- 0.5 * (nAf - nfAB)
  nfBB <- 0.5 * (nBf - nfAB)
  K <- lgamma(nA+1) + lgamma(nB+1) + lgamma(nf+1) + lgamma(nm+1) - lgamma(mA+1) - lgamma(mB+1) - lgamma(nt+1)
  quotient <- nfAB * log2 - lgamma(nfAA+1) - lgamma(nfAB+1) - lgamma(nfBB+1)
  prob <- exp(K + quotient)
  ii <- nearlyEqual(prob,rep(psamp,length(nfAB)),eps)
  pv <- sum(prob[ii]) # sum of all tied cases
  iii <- ((!ii) & (prob < psamp))
  pv <- pv + sum(prob[iii]) # sum of all cases with smaller prob
  return(pv)
}

