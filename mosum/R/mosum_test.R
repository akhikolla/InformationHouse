#' Help function: asymptotic scaling.
#' @keywords internal
mosum.asymptoticA <- function(x) {
  return(sqrt(2*log(x)))
}

#' Help function: asymptotic shift
#' @keywords internal
mosum.asymptoticB <- function(x, K) {
  return(2*log(x) + 0.5*log(log(x)) + log((K^2+K+1)/(K+1)) - 0.5*log(pi))
}

#' MOSUM asymptotic critical value
#' 
#' Computes the asymptotic critical value for the MOSUM test.
#' @param n an integer value for the length of the input data
#' @param G.left,G.right integer values for the left and right moving sum bandwidth (G.left, G.right)
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}
#' @return a numeric value for the asymptotic critical value for the MOSUM test
#' @examples
#' x <- testData(lengths = rep(100, 3), means = c(0, 5, -2), sds = rep(1, 3), seed = 1234)$x
#' m <- mosum(x, G = 40)
#' par(mfrow = c(2, 1))
#' plot(m$stat, type = "l", xlab = "Time", ylab = "", main = "mosum")
#' abline(h = mosum.criticalValue(300, 40, 40, .1), col = 4)
#' abline(v = m$cpts, col = 2)
#' plot(m, display = "mosum") # identical plot is produced 
#' @export
mosum.criticalValue <- function(n, G.left, G.right, alpha) {
  G.min <- min(G.left, G.right)
  G.max <- max(G.left, G.right)
  K <- G.min / G.max
  return((mosum.asymptoticB(n/G.min,K) - log(log(1/sqrt(1-alpha))))/mosum.asymptoticA(n/G.min))
}

#' MOSUM asymptotic p-value
#' 
#' Computes the asymptotic p-value for the MOSUM test.
#' @param z a numeric value for the observation
#' @param n an integer value for the length of the input data
#' @param G.left,G.right integer values for the left moving sum bandwidth (G.left,G.right)
#' @return a numeric value for the asymptotic p-value for the asymmetric MOSUM test
#' @keywords internal
mosum.pValue <- function(z, n, G.left, G.right=G.left) {
  G.min <- min(G.left, G.right)
  G.max <- max(G.left, G.right)
  K <- G.min / G.max
  return(1-exp(-2*exp(mosum.asymptoticB(n/G.min,K) - mosum.asymptoticA(n/G.min)*z)))
}
