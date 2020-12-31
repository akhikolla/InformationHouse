#' Calculate distance between clusters
#' 
#' \code{clusterdist} calculates the distance between each pair of clusters.
#' The distance measure is based on a symmetric Kullback-Leibler divergence.
#' 
#' 
#' @param fit A \code{\link{disclapmixfit}} object.
#' @param ... Not used
#' @return A distance matrix
#' @seealso \code{\link{disclapmix-package}} \code{\link{disclapmix}}
#' \code{\link{disclapmixfit}} \code{\link{clusterprob}}
#' \code{\link{predict.disclapmixfit}} \code{\link{print.disclapmixfit}}
#' \code{\link{summary.disclapmixfit}} \code{\link{simulate.disclapmixfit}}
#' %\code{\link{haplotype_diversity}}
#' \code{\link[disclap:disclap-package]{disclap}}
#' @keywords distance clusters
#' @export
clusterdist <- function(fit, ...) {
  if (!is(fit, "disclapmixfit")) {
    stop("fit must be a disclapmixfit")
  }
  
  if (nrow(fit$y) < 2L) {
    stop("Fit must have at least two clusters") 
  }
  
  symDist <- function(p1, z1, p2, z2) {
      KL <- function(p1, p2, m) {
        a <- (2*p1*log(p1)) / (1 - p1^2)
        b <- log( ((1-p1)*(1+p2)) / ((1+p1)*(1-p2)) )
        c <- (2*p1^(m+1) - m*(p1^2 - 1) ) / (1 - p1^2)
        return(a + b - c*log(p2))
      }

      m <- abs(z1 - z2)
      return(KL(p1, p2, m) + KL(p2, p1, m))
  }
  
  clusters <- nrow(fit$disclap_parameters)
  distmat <- matrix(NA, nrow = clusters, ncol = clusters)
  
  for (j1 in 1:(clusters-1)) {
    for (j2 in (j1+1):clusters) {
      kldist <- 0
      
      for (k in 1:ncol(fit$disclap_parameters)) {
        pj1k <- fit$disclap_parameters[j1, k]
        pj2k <- fit$disclap_parameters[j2, k]
        kldist <- kldist + symDist(p1 = pj1k, z1 = fit$y[j1, k], p2 = pj2k, z2 = fit$y[j2, k])
      }
      
      distmat[j2, j1] <- kldist
    }
  }
  
  distmat <- as.dist(distmat)
  
  return(distmat)
}

