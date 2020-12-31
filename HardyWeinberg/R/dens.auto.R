dens.auto <- function(m,f) {
  #
  # density in autosomal case with multiple alleles that uses gender
  #
  
  #
  # f and m should be lower triangular. no need for sorting according to minor allele count
  #
  nm  <- sum(m)
  nf  <- sum(f)
  nc  <- rowSums(f)+colSums(f)+rowSums(m)+colSums(m) # allele counts
  nt  <- sum(nc)
  h   <- sum(f[lower.tri(f)]) + sum(m[lower.tri(m)])
  fem <- f[lower.tri(f,diag=TRUE)]
  mal <- m[lower.tri(m,diag=TRUE)]
  num <- lgamma(nm+1) + lgamma(nf+1) + sum(lgamma(nc+1)) + h*log(2)
  den <- lgamma(nt+1) + sum(lgamma(mal+1)) + sum(lgamma(fem+1))
  pro <- exp(num-den)
  return(pro)
}
