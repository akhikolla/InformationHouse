###
# Calculate Kappa
#
# This function calculates kappa given a set
# @param set This is the set to calculate kappa from
# @param isSet TRUE if set FALSE if contigency table
# @param kappaThreshold if null then return kappa otherwise return list of kappa and above or not
# @return This returns the kappa of the set given
###
calcKappa <- function(set, isSet = TRUE, kappaThreshold = NULL) {
  if (!isSet) {
    set <- contingencyToSet(set[1, 1], set[2, 1], set[1, 2], set[2, 2])
  }

  if (
    length(which(set[, 1] == set[, 2] & set[, 1] == 1)) == nrow(set) |
    length(which(set[, 1] == set[, 2] & set[, 1] == 0)) == nrow(set)
  ) {
    return(1)
  }

  #calculates kappa by calculating the agreement and the baserates and then creating the adjacenty matrix
  agreement <- length(which(set[, 1] == set[, 2])) / nrow(set);
  baserate2 <- length(which(1 == set[, 2])) / nrow(set);
  baserate1 <- length(which(1 == set[, 1])) / nrow(set);
  randomAgreement <- baserate1 * baserate2 + (1 - baserate1) * (1 - baserate2);
  kappa <- (agreement - randomAgreement) / (1 - randomAgreement);

  if (is.null(kappaThreshold)) {
    return(kappa);
  } else {
    return(list(kappa = kappa, above = kappa > kappaThreshold))
  }
}
