# TODO: Add comment
# 
# Author: kley
###############################################################################


tvARMA_naive <- function(T = 128, a = list(),
    b = list(),
    sigma = function(u) {return(1)},
    innov = function(n) {rnorm(n, 0, 1)} ) {
  p <- length(a)
  q <- length(b)
  z <- c()
  x <- c()
  z <- innov(T+q)
  x[1:max(p,q)] <- innov(max(p,q)) ## Warm up a stationary ARMA(p,q) with a(0) and b(0)?
  for(i in (max(p,q)+1):T) {
    x[i] <- sigma(i/T) * z[q+i]
    if (q > 0) {
      for (j in 1:q) {
        x[i] <- x[i] + b[[j]](i/T) * sigma((i-j)/T) * z[q+i-j]
      }
    }
    if (p > 0) {
      for (j in 1:p) {
        x[i] <- x[i] + a[[j]](i/T)*x[i-j]
      }
    }
  }
  return(x)
}

set.seed(2581)

a1 <- function(u) {1.8 * cos(1.5 - cos(4*pi*u))}
a2 <- function(u) {-0.81}
sigma <- function(u) {1 + 12*u}

X_old <- tvARMA_naive(128, a = list(a1, a2), b = list(a2), sigma = sigma)
plot(X_old, type="l")

set.seed(2581)
X_new <- tvARMA(128, a = list(a1,a2), b = list(a2), sigma = sigma)
lines(X_new, col="red")

sum( (X_new - X_old)^2 )


library(microbenchmark)

microbenchmark(
      {tvARMA_naive(5000, a = list(), b = list(a1,a2), sigma = sigma)},
      {tvARMA(5000, a = list(), b = list(a1,a2), sigma = sigma)}
)