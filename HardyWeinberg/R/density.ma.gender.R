density.ma.gender <- function(m,f) {
#
# probability density for X chromosomal exact with multiple alleles
#
#
# f should be lower triangular. no need for sorting according to minor allele count
#
  nm  <- sum(m)
  nf  <- sum(f)
  nc  <- rowSums(f)+colSums(f)+m
  nt  <- sum(nc)
  h   <- sum(f[lower.tri(f)])
  fem <- f[lower.tri(f,diag=TRUE)]
  num <- lgamma(nm+1) + lgamma(nf+1) + sum(lgamma(nc+1)) + h*log(2)
  den <- lgamma(nt+1) + sum(lgamma(m+1)) + sum(lgamma(fem+1))
  pro <- exp(num-den)
  return(pro)
}
