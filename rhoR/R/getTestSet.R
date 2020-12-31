###
#' @title Get Test Set
#' 
#' @description
#' This function gets a \emph{testSet} from a larger \code{\link{codeSet}} given certain sampling parameters.
#' 
#' @details
#' A \emph{testSet} is a \code{\link{codeSet}} that is a subset of a larger \code{\link{codeSet}} with a given set of properties.  A \emph{testSet} is constructed by sampling (without replacement) P rows from rows in the larger \code{\link{codeSet}} where the first rater's code was 1, and then appending an additional sample (without replacement) of R rows taken at random from the larger \code{\link{codeSet}} excluding rows included in the first P rows sampled. P is computed as the minbaserate * length of the \emph{testset}. R is computed as testSetLength - P. The result of this sampling procedure is to create a sample with a minimum baserate regardless of the baserate of the larger \code{\link{codeSet}}.If \emph{testSetBaserateInflation} is set to zero, the function selects rows at random.
#' 
#' @param set The \code{\link{codeSet}} from which the \emph{testSet} is taken
#' @param testSetLength The length of the \emph{testSet} to be taken
#' @param testSetBaserateInflation The minimum guaranteed \code{\link{baserate}} of the \emph{testSet}. Default to 0
#' 
#' 
#' @export
#' @return A \code{\link{codeSet}} with the properties specified
###

getTestSet = function(set, testSetLength, testSetBaserateInflation = 0){
  if(length(unlist(set)) < testSetLength){stop("testSetLength must be less than the length of set")}
  if(testSetLength %% 1 | testSetLength < 1) {stop("testSetLength value must be a positive integer.")}
  if(testSetBaserateInflation < 0){stop("testSetBaserateInfation must be positive")}
  if(testSetBaserateInflation >= 1){stop("testSetBaserateInflation must be below 1")}
  if(length(which(set[,1] == 1 & set[,2] == 1)) == nrow(set)){
    return(set[1:testSetLength,])
  }else if(length(which(set[,1] == 0 & set[,2] == 0)) == nrow(set)){
    if(testSetBaserateInflation == 0){
      return(set[1:testSetLength,])
    }else{
      stop("Not enough positives in first rater to inflate to this level")
    }
  }
  return(getHandSet(set = set, handSetLength = testSetLength, handSetBaserate = testSetBaserateInflation, returnSet = T))
}