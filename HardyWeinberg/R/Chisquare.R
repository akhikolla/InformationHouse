Chisquare <- function (X) 
{
  mono <- FALSE
  n <- sum(X)
  Xhom <- X[homozyg(X)]
  Xhet <- X[heterozyg(X)]
  X <- c(min(Xhom), Xhet, max(Xhom))
  names(X) <- c("AA", "AB", "BB")
  p <- (2 * X[1] + X[2])/(2 * n)
  names(p) <- NULL
  q <- 1 - p
  if ((p == 0) | (q == 0)) {
    mono <- TRUE
  }
  obs <- X
  expected <- c(n * p^2, n * 2 * p * q, n * q^2)
  names(expected) <- names(X)
  chi <- (abs(obs - expected))^2
  if (!mono) {
    chi2 <- chi/expected
    chisq <- sum(chi2)
  }
  else {
    chisq <- NA
  }
  return(chisq)
}
