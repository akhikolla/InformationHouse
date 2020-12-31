#' Generate a Handset
#' 
#' @description Generate a vector representing indices of set, using the handSetBaserate
#' to determine the minimum number of indices that are positive
#' 
#' @param set matrix of two columns
#' @param handSetLength number of indices to find
#' @param handSetBaserate number between 0 and 1 to use as a minimum number of positive indices
#'
#' @return vector of indices from set
#' @export
getHandSetIndices = function(set, handSetLength = 20, handSetBaserate = 0.2) {
  positives = ceiling(handSetLength * handSetBaserate);
  posInd = which(set[, 1] == 1);
  
  if (positives > length(posInd)) {
    stop("Not enough positives in first rater to inflate to this level")
  }
  
  this.set = NULL
  if (positives > 0) {
    positiveIndices = posInd[sample.int(length(posInd), size=positives, replace=FALSE)]
    others = which(!(1:nrow(set) %in% positiveIndices))
    otherIndices = sample.int(length(others), size = (handSetLength - positives), replace = FALSE)
    this.set = c(positiveIndices, others[otherIndices])
  }
  else if (positives == 0) {
    theseIndices = sample.int(length(set),size=handSetLength,replace=FALSE)
    this.set = theseIndices
  }
  
  return(sample(this.set))
}
