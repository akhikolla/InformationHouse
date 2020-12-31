#' LD weights
#'
#' \code{ld.weights} computes LD weights for all markers, which is subsequently used to compute LD weighted GRM.
#'
#' \code{data} can either be the subject by marker numeric genotype matrix (with 0, 1 or 2 coding), or the matrix of marker genotypic correlations. The default option is to input genotype matrix.
#'
#' @param data numeric matrix.
#' @param input.genotype logical.
#' @importFrom kernlab ipop
#' @return A numeric vector of weights. Note that the sum of weights is not constrained to be 1. They should be scaled appropriately before computing the LD weighted GRM.
#' @examples
#' # simulate genotypes of 500 individuals at 100 markers
#' nsnp = 100 # number of SNPs
#' freq = runif(nsnp, 0.05, 0.95)
#' nhaplo = 1000 # number of founder haplotypes
#' haplo.mat = sim.haplotype(freq, nhaplo)
#' geno.mat = t(sapply(c(1:500), function(x) 4 - haplo.mat[2*x-1,] - haplo.mat[2*x,]))
#'
#' # compute unconstrained LD weights
#' ld.weights(geno.mat)
#' @export
ld.weights = function(data, input.genotype = TRUE){
  if(input.genotype){
    R = cor(as.matrix(data))^2
  }else{
    R = data^2
  }
  n = nrow(R)

  out = ipop(matrix(-1, n, 1), R, matrix(1, 1, n), matrix(0, 1, 1), matrix(0, n, 1), matrix(999, n, 1), matrix(999, 1, 1))
  w = out@primal
  return(w)
}
