###
# Get Recall
#
# This gets the recall from a pair
# @param kappa This is the kappa of the pair
# @param BR This is the baserate of the pair
# @param P this is the precision of the pair
# @return returns the recall given the parameters
###
getR = function(kappa, BR, P){
  #gets recall that corresponds with the kappa, baserate, and precision
  top = kappa*P;
  R = top / (2*P-2*BR-kappa+2*BR*kappa);
  return (R);
}
