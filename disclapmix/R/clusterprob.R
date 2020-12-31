#' Cluster origin probabilities for haplotypes
#' 
#' \code{clusterprob} calculates the cluster origin probabilities for
#' haplotypes.
#' 
#' 
#' @param fit A \code{\link{disclapmixfit}} object.
#' @param newdata The haplotypes to predict the cluster origin probabilities
#' for.
#' @param ... Not used
#' @return A matrix where the rows correspond to the rows in \code{newdata} and
#' the sum of each row is 1.
#' @seealso \code{\link{disclapmix-package}} \code{\link{disclapmix}}
#' \code{\link{disclapmixfit}} \code{\link{clusterdist}}
#' \code{\link{predict.disclapmixfit}} \code{\link{print.disclapmixfit}}
#' \code{\link{summary.disclapmixfit}} \code{\link{simulate.disclapmixfit}}
#' %\code{\link{haplotype_diversity}}
#' \code{\link[disclap:disclap-package]{disclap}}
#' @keywords distance clusters
#' @export
clusterprob <- function(fit, newdata, ...) {
  if (!is(fit, "disclapmixfit")) stop("object must be a disclapmixfit")
  
  if (!is.matrix(newdata) || !is.integer(newdata)) {
    stop("newdata is not an integer matrix")
  }
  
  if (ncol(newdata) != ncol(fit$y)) {
    stop("newdata must have the same number of columns as the original model was fitted for")
  }
  
  wic <- rcpp_calculate_wic(newdata, fit$y, fit$disclap_parameters, fit$tau)      
  vic_matrix <- rcpp_calculate_vic(wic)

  return(vic_matrix)
}

