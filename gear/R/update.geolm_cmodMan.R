#' @importFrom stats update
#' @rdname update.geolm
#' @export
update.geolm_cmodMan = function(object, mod, ...) {
  call = match.call()
  
  # create covariance matrix for observed data
  vediag = mod$evar / object$weights
  mod_evar0 = NULL
  v = mod$v
  cholv = chol(v)
  
  if (is.null(object$mu)) {
    ###compute matrix products for future use
    cholvix = backsolve(cholv, object$x, transpose = TRUE)
    vix = forwardsolve(cholv, cholvix, upper.tri = TRUE)
    xtvix = crossprod(cholvix)
    cholxtvix = chol(xtvix)
    vcov = chol2inv(cholxtvix)
    map2coeff = solve_chol(cholxtvix, t(vix))
    coeff = c(map2coeff %*% object$y)
    names(coeff) = names(object$coeff)
    resid = object$y - object$x %*% coeff
  } else {
    vix = NULL
    xtvix = NULL
    resid = object$resid
    coeff = NULL
    cholvix = NULL
    cholxtvix = NULL
    vcov = NULL
    map2coeff = NULL
  }
  
  cholviresid = backsolve(cholv, resid, transpose = TRUE)
  n = nrow(v)
  loglik = -n / 2 * log(2 * pi) -
    sum(2 * log(diag(cholv))) / 2 -
    crossprod(cholviresid)[1, 1] / 2
  
  return(structure(
    list(
      y = object$y,
      x = object$x,
      loglik = loglik,
      resid = resid,
      coeff = coeff,
      mu = object$mu,
      v = v,
      cholv = cholv,
      cholvix = cholvix,
      cholviresid = cholviresid,
      cholxtvix = cholxtvix,
      vcov = vcov,
      map2coeff = map2coeff,
      formula = object$formula,
      coordnames = object$coordnames,
      mod = mod,
      weights = object$weights,
      vediag = vediag,
      evar = mod$evar,
      coords = object$coords,
      # longlat = mod$longlat,
      call = call
    ),
    class = c("geolm_cmodMan",
              "geolm")
  ))
}
