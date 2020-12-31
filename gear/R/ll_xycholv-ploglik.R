#' @export
#' @rdname ll_xycholv
ploglik_xycholv = function(x, y, cholv, mu = NULL,
                           reml = FALSE, minus2 = TRUE,
                           return_ll = TRUE) {
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
  sigmasq = crossprod(cholviresid)[1,1]/n
  if (reml) {
    sigmasq = sigmasq * n / (n - ncolx)
  }
  
  # both are -2 the loglikelihood
  if (!reml) {
    ll = n * log(2*pi) + n * log(sigmasq) + sum(2 * log(diag(cholv))) + n
  } else {
    # pp. 259-263 Schabenberger and Gotway
    ll = 
      # log(det(v))
      sum(2 * log(diag(cholv))) +
      determinant(xtvix, logarithm = TRUE)$mod -
        determinant(crossprod(x), logarithm = TRUE)$mod +
      (n - ncolx) * log(sigmasq) + (n - ncolx) * log(2*pi) + n - ncolx
    attributes(ll) = NULL
  }
  if (!minus2) {
    ll = ll / (-2)
  }
  if (return_ll) {
    return(ll)
  } else {
    return(sigmasq)
  }
}
