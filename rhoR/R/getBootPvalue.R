###
# Get Boot P Value
#
# This function gets the p value of the result in the distribution
# @param distribution The distribution that the result falls into
# @param result The value that you are evaluating where in the distribution it falls
# @return This returns the p value from the result into the distribution
###
getBootPvalue = function(distribution, result) {
  #returns the percentage of the time that the distribution was greater or equal to the observed kappa
  N = length(distribution);
  if(result < mean(distribution, na.rm = T)) {
    #if the result is less than the mean of the distribution, than the p value is 1
    return(1);
  } else {
    #return the number of times that the distribution is greater than the result as a percentage of the total number of items in the distribution
    return(length(which(distribution >= result)) / N);
  }
}
