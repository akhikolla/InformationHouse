dlevene <-
function(N) {
  n <- sum(N)                           # sample size
  nt <- 2*n
  ni <- rowSums(N) + colSums(N)         # allele counts
  gt <- N[lower.tri(N,diag=TRUE)]       # genotype counts  
  h  <- sum(N[lower.tri(N,diag=FALSE)]) # total of heterozygote counts
  num <- lgamma(n+1) + sum(lgamma(ni+1)) + h*log(2)
  den <- lgamma(nt+1) + sum(lgamma(gt+1))
  y <- exp(num-den)
  return(y)
}
