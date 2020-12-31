
# The true characteristic funtion, i.e. the one for the increment \Delta_{k, k} X
phi_true <- function(t, sigma, alpha, k, H) {
  hkrnorm = Norm_alpha(h_kr, alpha = alpha, k = k, H = H, r = 1, l = 0, rel.tol = 1.2e-3)$result
  exp(-abs(t * sigma * hkrnorm)^alpha)
}

# The integrand in the minimal contrast estimator, excluding the weight function
integrand <- function(t, path, sigma, alpha, k, H, r = 1) {
  ((phi(t = t, k = k, path = path, freq = "L") - phi_true(t = t, sigma = sigma, alpha = alpha, k = k, H = H)) * (t > 0))^2
}



#' Statistical estimator of sigma, alpha and H in low frequency setting based on minimal contrast estimation comparing the empirical characteristic function with the true one
#'
#' Estimates H using the \code{\link{H_hat}} function while sigma and alpha are obtained via
#' \deqn{\arg\min_{\sigma, \alpha} \int_{0}^{\infty} (\varphi_n(t) - \varphi_{\sigma, \alpha, H_{hat}}(t))^2 \exp(-t^2/2) d t},
#' where \eqn{\varphi_n} is the empirical characteristic function, see \code{\link{phi}}, and \eqn{\varphi_{\sigma, \alpha, H_{hat}}} is the characteristic function of the kth order increment wrt the parameters \eqn{\sigma, \alpha, H_{hat}}, see also \code{\link{increment}}.
#'
#' @param path low frequency sample path from which the parameters should be estimated.
#' @param k order of increments.
#' @param p any real number, the power used for \code{\link{H_hat}}.
#' @param order_GH number of weights in the Gauss-Hermite approximation of the integral, see the \code{gauss.hermite} function from the spatstat package.
#'
#' @details This algorithm approximates the above integral using Gauss-Hermite quadrature and uses the \code{L-BFGS-B} method from the \code{optim} function to minimize over the parameters sigma and alpha.
#' Due to numerical problems estimation of sigma below 0.01 and alpha or H below 0.05 is currently not possible.
#' @examples
#' m0 = 256
#' M0 = 600
#' alpha0 = 1.8
#' H0 = 0.8
#' sigma0 = 0.3
#' n = 100
#' X <- path(N = n, m = m0, M = M0, alpha = alpha0, H = H0, sigma = sigma0, freq = 'L')$lfsm
#' MinContrastEstim(path = X, k = 2, p = 0.4, order_GH = 8)
#'
#' @references \insertRef{LP2019}{rlfsm}
#' @export
#'
MinContrastEstim <- function(path, k, p, order_GH) {

  # Cutoffs, lower bounds for estimates on sigma and alpha
  sigma_low = 1e-3
  alpha_low = 5e-2

  # Estimate H
  H_est = H_hat(p = p, k = k, path = path)

  if (H_est < 0.05) {
    stop("H estimate too small")
  }

  # The integral \int_{0}{\infty} (phi_n(t) - phi(t))^2 exp(-t^2) d t, as a function of sigma and alpha
  func_aux <- function(z) {
    # Split argument into two, IIRC it's necessary for using the "optim" function
    a1 = z[1]
    a2 = z[2]
    # Evaluate the integrand function with the specific values
    integrand_eval <- function(t){
      integrand(t = t, path = path, sigma = a1, alpha = a2, k = k, H = H_est)
    }
    spatstat::gauss.hermite(integrand_eval, mu = 0, sd = 1, order = order_GH)
  }

  # Arg-minimize the integral over sigma and alpha, i.e., calculate the minimal contrast estimator
  theta_est = optim(c(2, 1), func_aux, method = 'L-BFGS-B', lower = c(sigma_low, alpha_low), upper = c(Inf, 2))

  # Return estimates, including H
  c(theta_est$par[1], theta_est$par[2], H_est)
  list(sigma=theta_est$par[1], alpha=theta_est$par[2], H=H_est)
}
