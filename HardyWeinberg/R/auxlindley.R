auxlindley <-
function(beta,alpha,nAA,nBB,N) {
  lnum <- nAA*(alpha + beta) + nBB*(alpha - beta)
  lden <- N*log((1+exp(alpha)*cosh(beta)))
  y <- exp(lnum-lden)
  return(y)
}
