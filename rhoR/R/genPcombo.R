###
# Generate Precision Combo
#
# This function generates a precision that fits within the range and will create a valid set
# @param kappa This is the kappa of the desired pair
# @param Pmin This is the minimum precision desired
# @param Pmax This is the maximum precision desired
# @param BR This is the baserate of the desired pair
# @return This returns a precision value that will generate a valid set with the other parameters
####

genPcombo = function(kappa, Pmin, Pmax, BR){
  #gets a random precision in the precision range
  currPrecision = stats::runif(1, Pmin, Pmax);
  #checks is the combination of precision, kappa, and baserate is valid.  If it is vaid, return the precision, if not, recall the function until a valid precision is found
  if(!checkBRPKcombo(BR, currPrecision, kappa)){
    genPcombo(kappa, Pmin, Pmax, BR);
  }else{
    return(currPrecision);
  }
}
