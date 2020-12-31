lm.wfit.light <- function (x, y, w, tol = 1e-07) {
  x.asgn <- attr(x, "assign")
  zero.weights <- any(w == 0)
  save.r <- y
  save.f <- x

  if (zero.weights) {
    save.w <- w
    ok <- w != 0
    nok <- !ok
    w <- w[ok]
    x <- x[ok,  ,drop = FALSE]
    n <- nrow(x)
    y <-y[ok]
  }
  wts <- sqrt(w)
  out <- .lm.fit(x * as.numeric(wts), y * wts, tol)
  out$residuals <- save.r - save.f%*%out$coef
  return(out)
}
