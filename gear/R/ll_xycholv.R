#' Compute the log-likelihood of a model
#'
#' \code{ll_xycholv} computes the log-likelihood of
#' multivariate normal data using components typically found
#' in a \code{geolm}. For \code{ploglik_xycholv}, 
#' \code{cholv} is the cholesky decomposition of the 
#' covariance matrix for the observed data after dividing
#' the matrix by the (estimated) \code{psill}. See the
#' examples below.
#' Depending on parameter choices, the
#' function can return the log-likelihood, the restricted
#' log-likelihood, -2 times the log-likelihood or restricted
#' log-likelihood, or the estimated partial sill for both a
#' maximum likelihood and restricted maximum likelihood
#' setting. This is intended to be an internal function, so
#' minimal error checking is done.
#'
#' @inheritParams geolm_fit
#' @param cholv The cholesky decomposition of the covariance
#'   matrix of \code{y} (or of that matrix divided by \code{psill} for
#'   \code{loglik_xycholv}.
#' @param mu A single numeric value indicating the assumed
#' mean of the underlying process.
#' @param reml A logical value. Should the Restricted
#'   Maximum Likelihood be returned. The default is
#'   \code{FALSE}.
#' @param minus2 A logical value. Should -2 times the
#'   log-likelihood be returned. The default is \code{TRUE}.
#' @param return_ll A logical value. Should the
#'   log-liklihood be returned? Default is \code{TRUE}. If
#'   \code{FALSE}, the estimated partial sill is returned.
#'
#' @return A likelihood value, -2 times the likelihood
#'   value, or the estimated partial sill, depending on the
#'   user's argument choices.
#' @export
#' @references Statistical Methods for Spatial Data
#'   Analysis. Oliver Schabenberger and Carol A. Gotway
#'   (Chapman & Hall/CRC Press) 2005. pp. 259-263
#' @examples
#' y = rnorm(10)
#' x = matrix(rep(1, length(y)))
#' coords = matrix(runif(length(y) * 2), ncol = 2)
#' d = as.matrix(dist(coords))
#' pv = exp(-d/3) + 0.1 * diag(length(y))
#' est_psill = ploglik_xycholv(x, y, chol(pv), return_ll = FALSE)
#' v = pv * est_psill
#' # same result
#' ploglik_xycholv(x, y, chol(pv), minus2 = FALSE)
#' ll_xycholv(x, y, chol(v), minus2 = FALSE)
ll_xycholv = function(x, y, cholv, mu = NULL, reml = FALSE, minus2 = TRUE, return_ll = TRUE) {
  n = length(y)
  if (!is.null(mu) & reml) {
    stop("reml is only possible when mu = NULL")
  }
  if (!is.null(mu)) {
    resid = y - mu
    cholviresid = backsolve(cholv, resid, transpose = TRUE)
    ncolx = 0
  } else {
    cholvix = backsolve(cholv, x, transpose = TRUE)
    vix = forwardsolve(cholv, cholvix, upper.tri = TRUE)
    xtvix = crossprod(cholvix)
    cholxtvix = chol(xtvix)
    map2coeff = solve_chol(cholxtvix, t(vix))
    coeff = c(map2coeff %*% y)
    resid = y - x %*% coeff
    cholviresid = backsolve(cholv, resid, transpose = TRUE)
    ncolx = ncol(x)
  }
  rtvir = crossprod(cholviresid)[1,1]
  
  # both are -2 the loglikelihood
  if (!reml) {
    ll = n * log(2*pi) + rtvir + sum(2 * log(diag(cholv)))
  } else {
    # pp. 259-263 Schabenberger and Gotway
    ll = 
      # log(det(v))
      sum(2 * log(diag(cholv))) +  
      determinant(xtvix, logarithm = TRUE)$mod - 
      determinant(crossprod(x), logarithm = TRUE)$mod + 
      rtvir + (n - ncolx) * log(2*pi)
    attributes(ll) = NULL
  }
  if (!minus2) {
    ll = ll / (-2)
  }
  if (return_ll) {
    return(ll)
  } else {
    den = ifelse(reml, n - ncolx, n)
    return(rtvir/den)
  }
}
