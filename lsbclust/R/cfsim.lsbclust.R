#' Compare LSBCLUST Simulation Results
#' 
#' This function compares cluster membership and parameter estimates for the results of
#' \code{\link{lsbclust}} on simulated data to the true underlying values.
#' 
#' @param fitted An object of class \code{lsbclust} containing the fitted results.
#' @param actual An object of class \code{lsbclust_sim} containing the simulated data.
#' @param method The type of statistics to calculate, passed to \code{\link{cl_agreement}}
#' @inheritParams int.lsbclust
#' @method cfsim lsbclust
#' @export
#' @importFrom clue solve_LSAP
#' @examples 
#' ## Simulate LSBCLUST data, fit LSBCLUST, and compare
#' set.seed(1)
#' dat <- rlsbclust(ndata = 1, nobs = 100, size = c(10, 8), nclust = c(5, 4, 6, 5))
#' res <- lsbclust(data = dat[[1]]$data, nclust = c(5, 4, 6, 5))
#' cfsim(res, dat[[1]])
cfsim.lsbclust <- function(fitted, actual, method = c("diag", "cRand")) {
  
  ## Check
  if (!is(actual, "lsbclust_sim"))
    stop("Argument 'actual' must be of class 'lsbclust_sim'.")
  
  ## Check that C and D is of same length
  if (length(fitted$interactions$C) != length(actual$interactions$C))
    stop("Arguments 'fitted' and 'actual' contain differing numbers of C matrices.")
  if (length(fitted$interactions$D) != length(actual$interactions$D))
    stop("Arguments 'fitted' and 'actual' contain differing numbers of D matrices.")
  
  ## Reorder clusters to account for label switching
  len <- ncol(fitted$cluster)
  lst <- vector(mode = "list", length = len)
  
  ## Function to apply to each cluster type
  assess <- function(name) {
    
    ## Get confusion matrix with reordered columns (actual vs estimated)
    confmat <- table(actual$cluster[, name], fitted$cluster[, name])
    ord <- solve_LSAP(confmat, maximum = TRUE)
    confmat <- confmat[, ord]
    
    ## Get clue states
    stats <- sapply(method, function(x) cl_agreement(fitted[[name]], actual[[name]], method = x))
    
    ## Get MSE of parameter estimates (reorder before calculation)
    if (name == "interactions") {
      if (length(fitted$interactions$C) == 1L) {
        mse_c <- mean((unlist(fitted$interactions$C) - unlist(actual$interactions$C))^2)
        mse_c2 <- mean((unlist(fitted$interactions$C) + unlist(actual$interactions$C))^2)
      } else {
        mse_c <- mean((unlist(fitted$interactions$C[ord]) - unlist(actual$interactions$C))^2)
        mse_c2 <- mean((unlist(fitted$interactions$C[ord]) + unlist(actual$interactions$C))^2)
      }
      if (length(fitted$interactions$D) == 1L) {
        mse_d <- mean((unlist(fitted$interactions$D) - unlist(actual$interactions$D))^2)
        mse_d2 <- mean((unlist(fitted$interactions$D) + unlist(actual$interactions$D))^2)
      } else {
        mse_d <- mean((unlist(fitted$interactions$D[ord]) - unlist(actual$interactions$D))^2)
        mse_d2 <- mean((unlist(fitted$interactions$D[ord]) + unlist(actual$interactions$D))^2)
      }
      ## Determine which sign is the best
      if (mse_c2 < mse_c && mse_d2 < mse_d) 
        list(confusion = confmat, stats = c(stats, loss = fitted$interactions$minloss, 
                                            MSE_C = mse_c2, MSE_D = mse_d2))
      else 
        list(confusion = confmat, stats = c(stats, loss = fitted$interactions$minloss, 
                                            MSE_C = mse_c, MSE_D = mse_d))
    } else {
      if (name == "overall") {
        mse <- mean((fitted[[name]]$centers[ord] - drop(actual[[name]]$centers)^2))       
      } else {
        mse <- mean((fitted[[name]]$centers[ord, , drop = FALSE] - actual[[name]]$centers)^2)  
      }
      list(confusion = confmat, stats = c(stats, MSE = mse))
    }
  }

  ## Apply to all non-null components of fitted and actual
  ind_res <- !sapply(fitted[c("overall", "rows", "columns", "interactions")], is.null)
  ind_sim <- !sapply(actual[c("overall", "rows", "columns", "interactions")], is.null)
  if (any(ind_res != ind_sim))
    warning("Different components detected in 'fitted' and 'actual' -- using only those present in both.")
  out <- lapply(c("overall", "rows", "columns", "interactions")[ind_res & ind_sim], assess)
  names(out) <- c("overall", "rows", "columns", "interactions")[ind_res & ind_sim]
  
  ## Get global MSE (approximation to true data)
  mse_global <- mean((fitted(fitted) - actual$actual)^2)
  out <- c(global = list(c(MSE_global = mse_global)), out)
  
  return(out)
}