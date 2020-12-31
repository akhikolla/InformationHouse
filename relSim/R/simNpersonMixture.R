#' Simulate and count unique alleles in N person mixtures
#' 
#' This function simulates N persons mixtures using the supplied frequencies and records
#' the number of times they share 1, 2, \ldots, 2N alleles locus by locus. 
#' 
#' @param freqs a set of allele frequencies. The format can be found in \code{\link{readFreqs}}
#' @param numContributors the number of contributors to each mixture. Must be >= 2.
#' @param numIterations the number of N person mixtures to simulate in total.
#' 
#' @return an object of class npmresult
#' @export
simNpersonMixture = function(freqs, numContributors, numIterations = 10000){
  if(numContributors < 2){
    stop("numContributors must be >= 2")
  }
  
  locusProbs = .simNpersonMixture(freqs$freq, numContributors, numIterations)
  rownames(locusProbs) = freqs$loci
  
  res = vector(length = numContributors - 1, mode = "list")
  names(res) = 1:(numContributors - 1)
  
  nc = 1
  while(nc <= (numContributors - 1)){
    i1 = 2 * nc
    byLoc = apply(locusProbs[,1:i1], 1, sum)
    prob = prod(byLoc)
    res[[nc]] = list(byLoc = byLoc, prob = prob)
    nc = nc + 1
  }

  result = list(raw = locusProbs, byNC = res, summary = sapply(res, function(x){x$prob}))
  class(result) = "npmresult"

  return(result)
}