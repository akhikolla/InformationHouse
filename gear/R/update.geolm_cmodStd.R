#' @rdname update.geolm
#' @export
update.geolm_cmodStd = function(object, mod, ...) {
  call = match.call()
  # create covariance matrix for observed data
  # mod_evar0 = mod
  # mod_evar0$evar = 0
  # diagonal of covariance matrix for errors
  vediag = mod$evar/object$weights

  # compute distance matrix
  if (mod$ratio < 1) {
    d = ganiso_d(object$coords, object$coords, invert = mod$invert)
  } else {
    d = geodist(as.matrix(object$coords), longlat = mod$longlat)
  }

  # create covariance matrix for observed data
  v = evaluate(mod, d = d, e = FALSE) + diag(vediag)

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
  loglik = -n/2*log(2*pi) -
    sum(2 * log(diag(cholv)))/2 -
    crossprod(cholviresid)[1,1]/2
  # loglik = ll_xycholv(x = object$x, y = object$y, cholv = cholv,
  #                     mu = object$mu, minus2 = FALSE)

  return(structure(list(y = object$y,
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
                        # cmod_evar0 = mod_evar0,
                        weights = object$weights,
                        vediag = vediag,
                        evar = mod$evar,
                        coords = object$coords,
                        # longlat = mod$longlat,
                        call = call),
                   class = c("geolm_cmodStd", "geolm")))
}
