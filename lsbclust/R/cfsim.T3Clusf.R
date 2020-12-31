#' Compare LSBCLUST Simulation Results
#' 
#' This function compares cluster membership and parameter estimates for the results of
#' \code{\link{T3Clusf}} on simulated data, using \code{\link{rlsbclust}}, 
#' to the true underlying values.
#' 
#' @param fitted An object of class \code{lsbclust} containing the fitted results.
#' @param actual An object of class \code{lsbclust_sim} containing the simulated data.
#' @inheritParams int.lsbclust
#' @method cfsim T3Clusf
#' @export
#' @importFrom clue solve_LSAP
#' @examples
#' ## Simulate LSBCLUST data, fit T3Clusf on double-centered data, and compare
#' set.seed(1)
#' dat <- rlsbclust(ndata = 1, nobs = 100, size = c(10, 8), nclust = c(5, 4, 6, 5))
#' dat[[1]]$data <- carray(dat[[1]]$data)
#' res <- T3Clusf(X = dat[[1]]$data, Q = 2, G = 5)
#' cfsim(res, dat[[1]])
cfsim.T3Clusf <- function(fitted, actual, method = c("diag", "cRand")) {
  
  ## Check
  if (!is(actual, "lsbclust_sim"))
    stop("Argument 'actual' must be of class 'lsbclust_sim'.")
  
  ## Reorder clusters to account for label switching
  len <- length(fitted$sizes)
  lst <- vector(mode = "list", length = len)
  
  ## Get confusion matrix with reordered columns (actual vs estimated)
  if (fitted$alpha == 1) {
    confmat <- table(actual$cluster[, "interactions"], fitted$cluster)    
  } else {
    confmat <- aggregate(x = fitted$G, by = list(actual$cluster[, "interactions"]), sum)
    confmat <- as.matrix(confmat[, -1])
  }
  if (nrow(confmat) <= ncol(confmat)) {
    ord <- solve_LSAP(confmat, maximum = TRUE)
    confmat <- confmat[, ord]
  } else {
    ord <- solve_LSAP(t(confmat), maximum = TRUE)
    confmat <- confmat[ord, ]
  } 
  
  ## Get clue states
  stats <- sapply(method, function(x) cl_agreement(fitted, actual$interactions, method = x))
  
  ## Get MSE for approximation of the data
  MSE_overall <- mean((fitted(fitted) - actual$data)^2)
  stats <- c(stats, loss = fitted$minloss, MSE = MSE_overall)

  out <- list(global = list(MSE_global = MSE_overall), 
              interactions = list(confusion = confmat, stats = stats))
  return(out)
}