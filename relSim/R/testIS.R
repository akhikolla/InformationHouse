testIS = function(nc = c(3, 2), locus = 1, seed = 123456){
  set.seed(seed)
  np = 2 * nc[2]
  freqs = USCaucs$freqs
  for(i in 2:np){
    cmd = paste0("p",i," = IS(freqs[locus], numContributors = nc[1], maxPeaks =  ", i,", numIterations = 1e5)")
    cat(cmd, "\n")
    eval(parse(text = cmd))
    
    cmd = paste0("p",i,"$est = mean(exp(p",i,"$numerator - p",i,"$denominator))")
    eval(parse(text = cmd))
  }
  
  est = NULL
  cmd = paste0("est = c(", paste0("p", 1:np, "$est", collapse = ", "), ")")
  eval(parse(text = cmd))
  cat(paste(est, "\n"))
  ppoint = sum(est)
  
  set.seed(seed)
  ptail = IS(USCaucs$freqs[locus], numContributors = nc, numIterations = 1e5, bTail = TRUE)
  ptail = mean(exp(ptail$numerator - ptail$denominator))
  
  list(pointEstSum = ppoint, tailProb = ptail, ratio = ppoint / ptail)
}