################################################################################
#' Compute autocovariances of an AR(p) process
#'
#' This functions returns the autocovariances \eqn{Cov(X_{t-k}, X_t)} of a
#' stationary time series \eqn{(Y_t)} that fulfills the following equation:
#' \deqn{Y_t = \sum_{j=1}^p a_j Y_{t-j} + \sigma \varepsilon_{t},}
#' where \eqn{\sigma > 0}, \eqn{\varepsilon_t} is white noise and
#' \eqn{a_1, \ldots, a_p} are real numbers satisfying that the roots
#' \eqn{z_0} of the polynomial \eqn{1 - \sum_{j=1}^p a_j z^j}
#' lie strictly outside the unit circle. 
#' 
#' @name acfARp
#' @aliases acfARp
#' @export
#' 
#' @param a     vector \eqn{(a_1, \ldots, a_p)} of coefficients; default NULL,
#'              corresponding to p = 0, white noise with variance \eqn{\sigma^2},
#' @param sigma standard deviation of \eqn{\varepsilon_t}; default 1,
#' @param k     lag for which to compute the autocovariances.
#' 
#' @return Returns autocovariance at lag k of the AR(p) process.
#'
#' @examples
#' ## Taken from Section 6 in Dahlhaus (1997, AoS)
#' a1 <- function(u) {1.8 * cos(1.5 - cos(4*pi*u))}
#' a2 <- function(u) {-0.81}
#' # local autocovariance for u === 1/2: lag 1
#' acfARp(a = c(a1(1/2), a2(1/2)), sigma = 1, k = 1)
#' # local autocovariance for u === 1/2: lag -2
#' acfARp(a = c(a1(1/2), a2(1/2)), sigma = 1, k = -1)
#' # local autocovariance for u === 1/2: the variance
#' acfARp(a = c(a1(1/2), a2(1/2)), sigma = 1, k = 0)
################################################################################


acfARp <- function(a = NULL, sigma, k) {
  
  if (!is.numeric(a) && !is.null(a)) {
    stop("a needs to be NULL or a vector of p real numbers.")
  }
  
  if (!is.null(a) && (max(abs(polyroot(c(1,-a)))) <= 1)) {
    stop("a(z)=1-a1*z-...-ap*z^p may not have roots z0 with |z0|>=1.")
  }
  
  if (!(is.numeric(sigma) && sigma > 0)) {
    stop("sigma needs to be a positive real number.")
  }
  
  if (!(is.numeric(k) && k%%1 == 0 )) {
    stop("k needs to be an integer.")
  }
  
  # note that acf is symmetric
  k <- abs(k)
  
  if (length(a) == 0) {
    res <- sigma^2 * (k == 0)
  } else if (length(a) == 1) {
    res <- sigma^2 * a^k / (1 - a^2)
  } else {
    maxLag <- k
    p <- length(a)
    Sigma <- diag(c(sigma, rep(0, p - 1)))
    a_vec <- matrix(a, ncol=1)
    A <- rbind(t(a_vec), cbind(diag(rep(1, p - 1)), rep(0, p - 1)) )  
    
    Gamma0 <- matrix(solve(diag(rep(1,p^2)) - A %x% A) %*% as.vector(Sigma), ncol=p)
    
    if (k == 0) {
      res <- Gamma0[1, 1]
    } else {
      Gamma <- array(0, dim = c(maxLag, p, p))
      
      Gamma[1,,] <- A %*% Gamma0
      if (maxLag > 1) {
        for (i in 2:maxLag) {
          Gamma[i,,] <- A %*% Gamma[i-1,,]
        }
      }
      res <- Gamma[k, 1, 1]
    }
  }
  
  return(res)
}

