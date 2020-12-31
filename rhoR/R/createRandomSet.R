###
# Create A Random Set
#
# This creates a set with given parameters
# @param setLength The length of the set to be created
# @param baserate The baserate of the set to be created
# @param kappa The kappa of the set to be created
# @param minPrecision The minimum precision accepted in this set. Defaults to 0
# @param maxPrecision The maximum precision accepted in this set. Defaults to 1
# @return This returns a set with the given parameters
###
createRandomSet = function(setLength, baserate, kappaMin, kappaMax, minPrecision = 0, maxPrecision = 1, type = "set"){
  #generates a valid precision with the given kappa and baserate
  KP = generateKPs(1, baserate, kappaMin, kappaMax, minPrecision, maxPrecision)[[1]]
  kappa = KP$kappa
  precision = KP$precision
  recall = getR(kappa, baserate, precision)
  
  if(type == "ct") {
    fullSet = pr.ct(precision = precision, recall = recall, length = setLength, baserate = baserate);
  } 
  else {
    fullSet = prset(precision = precision, recall = recall, length = setLength, baserate = baserate);
  }
  
  return(fullSet);
}
