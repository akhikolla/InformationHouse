#' Combine output of the swapping algorithm
#'
#' This (non-exported) function combines the output from the swapping algorithm (Puccetti,
#' Rüschendorf and Vanduffel, 2020).
#'
#' @param x a three dimensional array (rows = subsets, columns = par, slices
#'   = samples) containing posterior samples for all subsets
#'
#' @return A \code{wasp} object, which can be further analyzed using the
#'   associated function \code{\link{summary.wasp}}.
#'
#' @source Puccetti, G., Rüschendorf, L. & Vanduffel, S. (2020). On the
#'   computation of Wasserstein barycenters, Journal of Multivariate Analysis,
#'   176.
#'

combine <- function(x){

  out = apply(x, 2, colMeans)
  return(out)

}

#' Compute the mode
#'
#' @param x a numeric vector
#'
#' @return The mode of a numeric vector as computed by the methods from Venter
#'   (1967).
#'
#' @source Venter, J.H. (1967). On estimation of the mode, Annals of
#'   Mathematical Statistics, 38(5), 1446-1455.
#'
#' @examples
#' library(waspr)
#' mode_est(pois_logistic[1,1,])
#'
#' @export

mode_est <- function(x){


  if(!is.numeric(x)){stop("x is not numeric")}
  hmode(x, 0.1)

  }

#' Compute the 95 percent Highest Posterior Density interval
#'
#' @inheritParams mode_est
#'
#' @return A vector containing the lower and upper bound of the 96% Highest
#'   Posterior Density interval of a numeric vector as computed by the methods
#'   from Venter (1967).
#'
#' @source Venter, J.H. (1967). On estimation of the mode, Annals of
#'   Mathematical Statistics, 38(5), 1446-1455.
#'
#' @examples
#' library(waspr)
#' hpd_est(pois_logistic[1,1,])
#'
#' @export

hpd_est <- function(x){

  if(!is.numeric(x)){stop("x is not numeric")}
  hmodeci(x, 0.95)

  }
