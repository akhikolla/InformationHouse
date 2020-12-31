genPKcombo = function(kappaDistribution, kappaProbability, precisionDistribution, precisionProbability, baserate){
  # currKappa = sample(kappaDistribution, 1, replace = FALSE, prob = kappaProbability)
  # currPrecision = sample(precisionDistribution, 1, replace = FALSE, prob = precisionProbability)
  currKappa = kappaDistribution[sample.int(length(kappaDistribution), 1, replace = FALSE, prob = kappaProbability)]
  currPrecision = precisionDistribution[sample.int(length(precisionDistribution), 1, replace = FALSE, prob = precisionProbability)]

  if(!checkBRPKcombo(baserate, currPrecision, currKappa)){
    precisionMin = (2*baserate*currKappa - 2*baserate - currKappa)/(currKappa - 2)
    indicies = which(precisionDistribution > precisionMin)
    
    if (length(indicies) == 0) {
      return(genPKcombo(kappaDistribution, kappaProbability, precisionDistribution, precisionProbability, baserate))
    }
    
    precisionDistribution = precisionDistribution[indicies]
    precisionProbability = precisionProbability[indicies]
    # currPrecision = sample(precisionDistribution, 1, replace = FALSE, prob = precisionProbability)
    currPrecision = precisionDistribution[sample.int(length(precisionDistribution), 1, replace = FALSE, prob = precisionProbability)]
    return(list(precision = currPrecision, kappa = currKappa))
  } else {
    return(list(precision = currPrecision, kappa = currKappa))
  }
}
