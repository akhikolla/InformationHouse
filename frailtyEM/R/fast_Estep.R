#' Fast fitting of the E step
#'
#' This function calculates the E step in a quicker way, without taking all the derivatives of the Laplace transform.
#' Such closed form solutions are only available for the gamma distribution (with or without left truncation) and
#' for the inverse Gaussian distribution (without left truncation).
#' For a data set with \code{K} clusters,
#' @param c Vector of length \code{K} of cumulative hazards, i.e. total accumulated hazards within a cluster
#' @param c_lt Vector of length \code{K} of cumulative hazard from 0 to the left truncation time
#' @param delta Vector of integers of length \code{K} of the number of events for each cluster
#' @param alpha,bbeta Parameters of the frailty distribution
#' @param pvfm Parameter for the PVF distribution, only matters in that case
#' @param dist One of 0 (for gamma), 1 (for stable) or 2 (for PVF)
#'
#' @return A \code{K x 4} matrix where the first column and the second column are the numerators
#' and the denominators of the frailty fraction, the third is the log-likelihood contribution, and the
#' last column is the expectation of the squared frailty (only used in calculating the information matrix)
#'
#' @keywords internal
fast_Estep <- function(c, c_lt = 0, delta, alpha, bbeta, pvfm, dist) {

  #if(!(dist == 0 | (dist == 2 & pvfm == -0.5))) stop("no fast option available here")

  res <- matrix(0, length(delta), 4)  # the 4th for Ex^2

  if(dist==0) {
    bbeta <- bbeta + c_lt
    res[,3] <- alpha * log(bbeta) - (alpha + delta)*log(bbeta + c) + lgamma(alpha + delta) - lgamma(alpha)
    res[,1] <- (alpha + delta)
    res[,2] <- (bbeta + c)
    res[,4] <- (alpha + delta) * (alpha + delta + 1) / (bbeta + c)^2
  }

  if(dist==2) {
   # warning("does not work with LT!")
    cc <- sqrt(2 * alpha * (c + alpha / 2))

    bk05 <-  besselK(cc, nu = delta - 0.5)

    res[,1] <- besselK(cc, nu = delta + 0.5)

    res[,2] <- sqrt(2 / alpha * c + 1) * bk05


    res[,3] <- (-delta / 2) * log(2/alpha * c + 1) +
      log(bk05)  - 0.5 * log(pi  / (2 * cc)) + cc +
      alpha * (1 - sqrt(1 + 2 / alpha * c))

    res[,4] <- besselK(cc, nu = delta + 1.5) /
      ((2 / alpha * c + 1) * bk05)

  }

  res
}
