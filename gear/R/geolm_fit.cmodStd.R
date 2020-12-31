## Linear model for geostatistical data when mod is of class modStd.
##
## \code{geolm} creates a geostatistical linear model object
## when \code{mod} is of class modStd.
## @param x The matrix of covariates.
## @param y The vector or matrix of observed data
## @param coords The coordinates of the observed data sets
## @param mu A single numeric value indicating the consant mean of the spatial process if simple kriging is desired.  Default is \code{NULL}, meaning that ordinary or universal kriging should be used.
## @param mod A covariance model object obtained of class modStd.
## @param weights An optional vector of weights for the errors to be used in the fitting process.  A vector that is proportional to the reciprocal variances of the errors, i.e., errors are assumed to be uncorrelated with variances \code{evar/weights}.  Default is \code{NULL}, meaning that the weights are uniformly 1.
## @param formula An object of class formula.  See Details.
## @param coordnames A vector of length 2 with the names of the columns in \code{data} containing the coordinates, e.g., \code{c("long", "lat")}.
## @param n The number of observations
geolm_fit.cmodStd <- function(x, y, coords, mu, mod, weights, formula, coordnames, n, call,
                              coeff_names){
  # diagonal of covariance matrix for errors
  vediag = mod$evar/weights

  # compute distance matrix
  if (mod$ratio < 1) {
    d = ganiso_d(coords, coords, radians = TRUE,
                 invert = mod$invert)
  } else {
    d = geodist(as.matrix(coords), longlat = mod$longlat)
  }

  # create covariance matrix for observed data
  v = evaluate(mod, d = d, e = FALSE) + diag(vediag)

  # cholesky of v
  cholv = chol(v)

  if (!is.null(mu)) { # if doing simple kriging
    vix = NULL
    xtvix = NULL
    resid = y - mu
    coeff = NULL
    cholvix = NULL
    cholxtvix = cholxtvix = NULL
    vcov = NULL
    map2coeff = map2coeff = NULL
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
    # map2coeff = forwardsolve(cholxtvix,
    #                          backsolve(cholxtvix, t(vix),
    #                                    transpose = TRUE),
    #                          upper.tri = TRUE)
    map2coeff = solve_chol(cholxtvix, t(vix))
    coeff = c(map2coeff %*% y)
    names(coeff) = coeff_names
    # coeff2 <- solve(xtvix, crossprod(vix, y))
    # range(coeff - coeff2)
    # compute residuals
    resid = y - x %*% coeff
  }

  cholviresid = backsolve(cholv, resid, transpose = TRUE)
  # viresid = solve(v, resid)
  # viresid2 = forwardsolve(cholv, cholviresid, upper.tri = TRUE)
  # range(viresid - viresid2)
  # compute log likelihood
  loglik = -n/2*log(2*pi) -
    #determinant(v, log = TRUE)$mod/2 -
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
                        formula = formula, coordnames = coordnames,
                        mod = mod,
                        weights = weights,
                        vediag = vediag,
                        evar = mod$evar,
                        coords = coords,
                        # longlat = mod$longlat,
                        call = call),
                   class = c("geolm_cmodStd", "geolm")))
}

## Compute the log likelihood for geostatistical data when mod is of class modStd.
## @param cholv The cholesky decomposition of the covariance matrix of the observed data
## @param cholviresid backsolve(cholv, resid, transpose = TRUE) where resid is y - x %*% coeff or y - mu for simple kriging
## @param n The number of observations
loglik.Std = function(cholv, cholviresid, n){
  -n/2*log(2*pi) - sum(2 * log(diag(cholv)))/2 - crossprod(cholviresid)[1,1]/2
}
