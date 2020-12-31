#' @export
#' @noRd
geolm_fit.cmodMan = function(mod, x, y, coords, mu, weights,
                             formula, coordnames, n, call,
                             coeff_names){
  # diagonal of covariance matrix for errors
  vediag = mod$evar/weights

  # create covariance matrix for observed data
  v = mod$v

  # cholesky of v
  cholv = chol(v)

  # if doing simple kriging
  if (!is.null(mu)) {
    vix = NULL
    xtvix = NULL
    resid = y - mu
    coeff = NULL
    cholvix = NULL
    cholxtvix = NULL
    vcov = NULL
    map2coeff = NULL
  } else {# if doing universal kriging
    ###compute matrix products for future use
    cholvix = backsolve(cholv, x, transpose = TRUE)
    vix = forwardsolve(cholv, cholvix, upper.tri = TRUE)
    # range(vix - vix2)
    # xtvix2 <- crossprod(vix, x)
    xtvix = crossprod(cholvix)
    cholxtvix = chol(xtvix)
    vcov = chol2inv(cholxtvix)
    # range(xtvix - xtvix2)
    #compute gls estimates of regression coefficients
    map2coeff = solve_chol(cholxtvix, t(vix))
    coeff = c(map2coeff %*% y)
    names(coeff) = coeff_names
    # compute residuals
    resid = y - x %*% coeff
  }

  cholviresid = backsolve(cholv, resid, transpose = TRUE)
  loglik = -n/2*log(2*pi) -
    sum(2 * log(diag(cholv)))/2 -
    crossprod(cholviresid)[1,1]/2

  return(structure(list(y = y, x = x, loglik = loglik,
                        resid = resid, coeff = coeff, mu = mu,
                        v = v,
                        cholv = cholv,
                        cholvix = cholvix, # cholviy = cholviy,
                        # vix = vix,
                        cholviresid = cholviresid,
                        cholxtvix = cholxtvix,
                        vcov = vcov,
                        # xtvix = xtvix,
                        map2coeff = map2coeff,
                        formula = formula,
                        coordnames = coordnames,
                        mod = mod,
                        weights = weights,
                        vediag = vediag,
                        evar = mod$evar,
                        coords = coords,
                        # longlat = mod$longlat,
                        call = call),
                   class = c("geolm_cmodMan", "geolm")))
}


