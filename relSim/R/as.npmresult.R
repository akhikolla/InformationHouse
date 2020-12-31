as.npmresult = function(locusProbs, freqs){
  rownames(locusProbs) = freqs$loci
  numContributors = ncol(locusProbs) / 2
  
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