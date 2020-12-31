#' Generate artificial data for ensemble verification using a signal-plus-noise model
#'
#' @param N number of forecasts and observations
#' @param mu.y expectation value of the observations
#' @param s.s standard deviation of the predictable signal
#' @param s.eps standard deviation of the unpredictable noise
#' @param mu.x expectation value of the ensemble
#' @param beta weighting parameter of the signal in the ensemble forecasting system
#' @param s.eta average spread of the ensemble
#' @param K number of members of the ensemble
#' @param mu.x.ref expectation value of the reference ensemble
#' @param beta.ref weighting parameter of the signal in the reference ensemble forecasting system
#' @param s.eta.ref average spread of the reference ensemble
#' @param K.ref number of members of the reference ensemble
#' @return A list with elements:
#' \describe{
#'   \item{obs}{N-vector of observations}
#'   \item{ens}{N*K matrix of ensemble members}
#'   \item{ens.ref}{N*K.ref matrix of reference ensemble members}
#' }
#' @details
#'  The function simulates data from the latent variable model:
#'
#'    y_t   = mu_y +     s_t    + eps_t
#'
#'    x_t,r = mu_x + beta * s_t + eta_t,r
#'
#' where y_t is the observation at time t, and x_t,r is the r-th ensemble member
#' at time t. The latent variable s_t is to be understood as the "predictable
#' signal" that generates correlation between observations and ensemble members.
#' If all arguments that end in ".ref" are specified, a reference ensemble is
#' returned to also test comparative verification.
#' @examples
#'  l <- GenerateToyData()
#'  with(l, EnsCrps(ens, obs))
#' @export

GenerateToyData <- function(N = 20, mu.y = 0, s.s = 7, s.eps = 6, mu.x = 0, beta = 0.2, s.eta = 8, K = 10, mu.x.ref = NA, beta.ref = NA, s.eta.ref = NA, K.ref = NA) {
  s <- rnorm(N, 0, s.s)
  y <- mu.y + s + rnorm(N, 0, s.eps)
  x1 <- mu.x + beta * s + matrix(rnorm(N*K, 0, s.eta), N, K)
  if (!any(is.na(c(mu.x.ref, beta.ref, s.eta.ref, K.ref)))) {
    x2 <- mu.x.ref + beta.ref * s + matrix(rnorm(N*K.ref, 0, s.eta.ref), N, K.ref)
    ret <- list(obs=y, ens=x1, ens.ref=x2)
  } else {
    ret <- list(obs=y, ens=x1)
  }
  return(ret)
}

