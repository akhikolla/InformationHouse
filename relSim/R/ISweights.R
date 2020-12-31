#' @importFrom multicool initMC allPerm
IS = function(freqs, numContributors = 4, maxPeaks = NULL, numIterations = 100, bTail = FALSE){
  if(numContributors[1] < 2){
    stop("numContributors must be >= 2")
  }
  
  # m = initMC(1:4)
  # perms = list(allPerm(m))
  # r = .IS(USCaucs$freqs, numIterations, 4, 4, perms)
  
  if(bTail){
    if(length(numContributors) != 2){
      message = paste("If you want tail probabilities, then you must supply a vector of length 2 for ",
                      "the number of contributors (numContributors). These must be (in the following order):",
                      "the actual number of contributors (N), and the apparent number of contributors.\n")
      stop(message)
    }
    
    if(numContributors[2] >= numContributors[1]){
      stop("The number of apparent contributors must be less than the number of actual contributors\n")
    }
    
    if(numContributors[2] < 1){
      stop("The number of apparent contributors must be greater than or equal to one.\n")
    }
    
    maxPeaks = 2 * numContributors[2]
    perms = vector(length = maxPeaks, mode = "list")
    
    for(np in 1:maxPeaks){
      if(np != 1){
        m = initMC(1:np)
        perms[[np]] = allPerm(m)
      }else{
        perms[[np]] = matrix(1, 1, 1)
      }
    }
    
    r = .IS(freqs, numIterations, numContributors[1], maxPeaks, perms, TRUE)
    ## This correction implies the probabilities under the target density
    ## are too small by a factor of maxPeaks = 2 * apparent contribs
    ## but not sure why - this should not affect the target probabilities
    ## or the importance probabilities are too big by a factor of 2 * apparent contribs
    r$numerator = r$numerator + log(maxPeaks - 2) 
    
  }else{
    if(is.null(maxPeaks) || maxPeaks < 1 || maxPeaks > 2 * numContributors){
      stop("If you are not estimating tail probabilities, then you must specify maxPeaks, the number of peaks/alleles you are interested in. 
           This should be a number between 1 and 2 * numContributors.\n")
    }
    
    if(maxPeaks > 1){
      m = initMC(1:maxPeaks)
      perms = list(allPerm(m))
    }else{
      perms = list(matrix(1, 1, 1))
    }
    r = .IS(freqs, numIterations, numContributors, maxPeaks, perms)
  }
  # perms = vector(length = maxPeaks, mode = "list")
  # for(np in 1:maxPeaks)
  # if(np != 1){
  #   m = initMC(1:np)
  #   perms[[np]] = allPerm(m)
  # }else{
  #   perms[[np]] = matrix(1, 1, 1)
  # }
  # 
  # w = ISprob(X, perms) - log(maxPeaks)
  # p = mean(exp(r$probs - w))
  # 
  # return(list(Alleles = X, perms = perms, weights = w, probs = r$probs, est = p))
  #m = initMC(1:np)
  return(r)
}

ISprobR = function(X, perms){
  freqs = X$f
  numAlleles = X$n
  numPerms = nrow(perms)
  
  r = 0
  
  for(j in 1:numPerms){
    p = s = freqs[perms[j, 1]];
    
    for(k in 2:numAlleles){
      pk = freqs[perms[j, k]];
      p  = p *  pk / (1 - s);
      s  = s + pk;
    }
    
    r = r + p
  }
  
  r = log(r)
  
  freqs = freqs / sum(freqs)
  counts = X$c - 1
  p2 = sum(counts * log(freqs) - log(factorial(counts))) + log(factorial(sum(counts)))
  
  return(r + p2) 
}