###
# This function calculates rho given the parameters of the initial set and the kappa of the observed handset and is
# called by rhoK()
###
calcRho = function(
  x, 
  OcSBaserate, 
  testSetLength, 
  testSetBaserateInflation = 0,
  OcSLength = 10000, 
  replicates = 800, 
  ScSKappaThreshold = 0.9, 
  ScSKappaMin = 0.40, 
  ScSPrecisionMin = 0.6, 
  ScSPrecisionMax = 1, 
  KPs = NULL
) {
  
  if(is.null(KPs)) {
    KPs = generateKPs(replicates, OcSBaserate, ScSKappaMin, ScSKappaThreshold, ScSPrecisionMin, ScSPrecisionMax)
  }
  if(length(KPs) < replicates) {
    replicates = length(KPs)
  }
  
  #creates distribution from random kappa sets
  savedKappas = rep(0,replicates);
  
  # cat("Runs: ", runs, "\n")
  # cat("KPs length:, ", length(KPs), "\n")
  for (i in 1:replicates) {
    KP = KPs[i]
    kappa = KP[[1]]$kappa
    precision = KP[[1]]$precision
    recall = getR(kappa, OcSBaserate, precision)
    fullSet = prset(precision = precision, recall = recall, length = OcSLength, baserate = OcSBaserate);
    
    #computes kappa of handset to add to distribution
    savedKappas[i] = getHandSet(set = fullSet, handSetLength = testSetLength, handSetBaserate = testSetBaserateInflation);
  }

  # compare observedKappa to distribution
  return(getBootPvalue(savedKappas, x));
}
