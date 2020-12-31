###
# Check Baserate, Precision, Kappa Combo
#
# This function checks to make sure the combination of BR, P, and K will create a valid set
# @param BR This is the supplied baserate
# @param P This is the supplied precision
# @param K This is the supplied kappa
# @return returns TRUE if the combination is valid, FALSE otherwise
###
checkBRPKcombo = function(BR, P, K) {
  #if right is less than P, then it is a valid combination. If not, then the pair of kappa, precision, and baserate does not correspond to a valid set
  right = (2*BR*K - 2*BR - K)/(K - 2);
  if(P > right){
    return(TRUE);
  }else{
    return(FALSE);
  }
}
