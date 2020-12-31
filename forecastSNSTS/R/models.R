################################################################################
#' Simulation of an tvARMA(p,q) time series.
#'
#' Returns a simulated time series \eqn{Y_{1,T}, ..., Y_{T,T}} that fulfills
#' the following equation:
#' \deqn{Y_{t,T} = \sum_{j=1}^p a_j(t/T) Y_{t-j,T} + \sigma(t/T) \varepsilon_{t} + \sum_{k=1}^q \sigma((t-k)/T) b_k(t/T) \varepsilon_{t-k},}
#' where \eqn{a_1, \ldots, a_p, b_0, b_1, \ldots, b_q} are real-valued functions on \eqn{[0,1]},
#' \eqn{\sigma} is a positive function on \eqn{[0,1]} and \eqn{\varepsilon_t}
#' is white noise.
#'
#' @name ts-models-tvARMA
#' @aliases tvARMA
#' @export
#' 
#' @importFrom stats rnorm
#'
#' @param T length of the time series to be returned
#' @param a list of p real-valued functions defined on \eqn{[0,1]}
#' @param b list of q real-valued functions defined on \eqn{[0,1]} 
#' @param sigma function
#' @param innov a function with one argument \code{n} that simulates a vector of
#'              the \code{n} residuals \eqn{\varepsilon_t}.  
#' 
#' @return Returns a tvARMA(p,q) time series with specified parameters.
#'
#' @examples
#' ## Taken from Section 6 in Dahlhaus (1997, AoS)
#' a1 <- function(u) {1.8 * cos(1.5 - cos(4 * pi * u))}
#' a2 <- function(u) {-0.81}
#' plot(tvARMA(128, a = list(a1, a2), b = list()), type = "l")
################################################################################
tvARMA <- function(T = 128, a = list(), b = list(), sigma = function(u) {return(1)},
    innov = function(n) {rnorm(n, 0, 1)} ) {
  
  p <- length(a)
  q <- length(b)
  
  grid <- (0:T)/T
  A <- matrix(ncol=p+1, nrow=T+1)
  A[, 1] <- rep(1, T+1)
  if (p >= 1) {
    for (j in 1:p) { A[, j+1] <- Vectorize(a[[j]])(grid) }
  }
  
  B <- matrix(ncol=q+1, nrow=T+1)
  B[, 1] <- rep(1, T+1)
  if (q >= 1) {
    for (j in 1:q) { B[, j+1] <- Vectorize(b[[j]])(grid) }
  }
  Sigma <- Vectorize(sigma)(grid)
  
  z <- innov(T+q)
  x_init <- innov(max(p,q)) ## Warm up a stationary ARMA(p,q) with a(0) and b(0)?
  
  return(tvARMAcpp(z, x_init, A, B, Sigma))
}