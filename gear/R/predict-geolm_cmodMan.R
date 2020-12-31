#' Predict method for geostatistical models
#'
#' \code{predict} calculates the predicted values at
#' specified locations.  The method can additionally provide
#' the mean square prediction error (mspe) and perform
#' conditional simulation.
#'
#' The \code{newdata} data frame must include the relevant
#' covariates for the prediction locations, where the
#' covariates are specified on the right side of the
#' \code{~} in \code{object$formula}.  \code{newdata} must
#' also include the coordinates of the prediction locations,
#' with these columns having the names provided in
#' \code{object$coordnames}.
#'
#' @param object An object produced by the \code{geolm}
#'   function.
#' @param newdata An optional data frame in which to look
#'   for the coordinates at which to predict. If omitted,
#'   the observed data locations are used.
#' @param nsim A non-negative integer indicating the number
#'   of realizations to sample at the specified coordinates
#'   using conditional simulation.
#' @param vop The cross-covariance matrix between the
#'   observed responses and the responses to predict.
#' @param vp The covariance matrix of the responses to
#'   predict.
#' @param return_type A character string indicating the type
#'   of object that should be returned. The default is
#'   \code{"\link[sp]{SpatialPointsDataFrame}"} for easy
#'   plotting of results (see Examples). Other options
#'   include \code{"\link[base]{data.frame}"},
#'   \code{"\link[gear]{geardf}"}, and \code{"\link[sf]{sf}"}.
#' @param dmethod The method used to decompose the
#'   covariance matrix for conditional simulation.  Valid
#'   options are \code{"chol"}, \code{"eigen"}, and
#'   \code{"svd"}.  The default is \code{"chol"}.
#' @param compute_mspe A logical value indicating whether
#'   the mean square prediction error should be calculated.
#'   Default is \code{TRUE}.
#' @param sp This argument will be deprecated in the future.
#'   Please use the \code{return_type} argument. A logical
#'   value indicating whether to object returned should be
#'   of class \code{\link[sp]{SpatialPointsDataFrame}} for
#'   easier plotting with the \code{sp} package.  Default is
#'   \code{NULL}.
#' @param ... Currently unimplemented.
#'
#' @return A \code{\link[base]{data.frame}},
#'   \code{\link[sp]{SpatialPointsDataFrame}},
#'   \code{\link[gear]{geardf}}, or \code{\link[sf]{sf}}
#'   object with the kriging predictions
#'   \code{pred}, kriging variance/mean-square prediction
#'   error (\code{mspe}), the root mean-square prediction
#'   error \code{mspe} (\code{rmspe}), and the conditional
#'   simulations \code{sim.1}, \code{sim.2}, etc.
#'   \code{sim.1}, \code{sim.2}, etc.
#'
#' @author Joshua French
#' @importFrom stats predict
#' @examples
#' # generate response
#' y = rnorm(10)
#' # generate coordinates
#' x1 = runif(10); x2 = runif(10)
#'
#' # data frame for observed data
#' data = data.frame(y, x1, x2)
#' coords = cbind(x1, x2)
#' d = as.matrix(dist(coords))
#' psill = 2 # partial sill
#' r = 4 # range parameter
#' evar = .1 # error variance
#' fvar = .1 # add finescale variance
#' # one can't generally distinguish between evar and fvar, but
#' # this is done for illustration purposes
#'
#' # manually specify an exponential covariance model
#' v = psill * exp(-d/r) + (evar + fvar) * diag(10)
#' mod_man = cmod_man(v = v, evar = evar)
#'
#' # coordinate names
#' cnames = c("x1", "x2")
#'
#' # geolm for universal kriging
#' gearmod_uk = geolm(y ~ x1 + x2, data = data, mod = mod_man,
#'                  coordnames = cnames)
#'
#' # newdata must have columns with prediction coordinates
#' # add 5 unsampled sites to sampled sites
#' newdata = data.frame(x1 = c(x1, runif(5)), x2 = c(x2, runif(5)))
#' newcoords = newdata[, cnames]
#' # create vop and vp using distances
#' dop = geodist(as.matrix(coords), as.matrix(newcoords))
#' dp = geodist(newcoords)
#'
#' # manually create cross-covariance and covariance for
#' # prediction locations
#' vop = psill * exp(-dop/r) + fvar * (dop == 0)
#' vp = psill * exp(-dp/r) + fvar * diag(nrow(newcoords))
#'
#' # prediction for universal kriging, with conditional simulation,
#' # using manual covariance matrices
#' pred_uk_man = predict(gearmod_uk, newdata, nsim = 2,
#'                       vop = vop, vp = vp, dmethod = "svd")
#'
#' # do the same thing, but using cmod_std
#'
#' # prediction for universal kriging, with conditional simulation
#' mod_std = cmod_std("exponential", psill = psill, r = r,
#'                    evar = evar, fvar = fvar)
#' gearmod_uk2 = geolm(y ~ x1 + x2, data = data, mod = mod_std,
#'                     coordnames = c("x1", "x2"))
#' pred_uk_std = predict(gearmod_uk2, newdata, nsim = 2, dmethod = "svd")
#'
#' # compare results
#' all.equal(pred_uk_man$pred, pred_uk_std$pred)
#' all.equal(pred_uk_man$mspe, pred_uk_std$mspe)
#' @export
predict.geolm_cmodMan =
  function(object, newdata, nsim = 0, vop, vp,
           return_type = "SpatialPointsDataFrame",
           dmethod = "chol", compute_mspe = TRUE, sp = NULL,
           ...) {
  if (!is.null(sp)) {
    arg_check_sp(sp)
    return_type = ifelse(sp, "SpatialPointsDataFrame", "data.frame")
  }
  arg_check_predict_geolm(coordnames = object$coordnames,
                          newdata = newdata, nsim = nsim,
                          return_type = return_type,
                          dmethod = dmethod,
                          compute_mspe = compute_mspe)
  if (return_type == "gearPredict") {
    warning("The 'gearPredict' return_type is being deprecated in favor of the 'geardf' return_type for easier plotting. See the Examples in predict.geolm_cmodStd or plot.geardf.")
    return_type = "geardf"
  }
  arg_check_predict_geolm_cmodMan(y = object$y, vop = vop,
                                  vp = vp)
  newcoords = as.matrix(newdata[,object$coordnames])

  # unneeded if simple kriging
  newx = NULL
  if (is.null(object$mu)) {
    newf = stats::delete.response(stats::terms(object$formula))
    newx = stats::model.matrix(newf, data = newdata)
    dimnames(newx) = NULL
  }

  # generate simulated data at observed and prediction locations
  if (nsim > 0) {
    n = nrow(object$coords)
    m = nrow(newcoords)

    # slightly different algorithm for simple kriging
    mu = 0
    if (!is.null(object$mu)) mu = object$mu

    newsim = decomp_cov(rbind(
      cbind(object$v - diag(object$vediag),vop),
      cbind(t(vop), vp)), method = dmethod) %*%
      matrix(stats::rnorm((n + m)*nsim), ncol = nsim) + mu

    # update various objects to include simulated data
    object$y = cbind(object$y, newsim[1:n,] +
                matrix(stats::rnorm(n*nsim,
                       sd = sqrt(object$vediag)),
                       nrow = n, ncol = nsim))

    if (!is.null(object$mu)) {
      object$cholviresid =
        backsolve(object$cholv, object$y - mu,
                  transpose = TRUE)
    } else {
      object$coeff = object$map2coeff %*% object$y
      object$cholviresid = backsolve(object$cholv,
                                     object$y - object$x %*% object$coeff,
                                     transpose = TRUE)
    }
  }

  cholvivop = forwardsolve(object$cholv, vop,
                           transpose = TRUE,
                           upper.tri = TRUE)

  # prediction
  pred = fitted(object, newx) + crossprod(cholvivop, object$cholviresid)

  # assume mspe is not computed
  mspe = NA
  if (compute_mspe) {
    #compute mspe
    # mspe1 = diag(vp)
    # mspe2 = colSums(cholvivop^2)
    mspe = diag(vp) - colSums(cholvivop^2)
    mspe3 = 0 # assume simple kriging
    if (is.null(object$mu)) {
      mspe3 = colSums(forwardsolve(object$cholxtvix,
                                    t(newx - crossprod(cholvivop,
                                                       object$cholvix)),
                                    transpose = TRUE, upper.tri = TRUE)^2)
    }
    # mspe = mspe1 - mspe2 + mspe3
    mspe = mspe + mspe3
    mspe[mspe < 0] = 0 # numerical imprecision correction
  }

  # data frame of kriging information
  kdtf = data.frame(pred = pred[,1], mspe = mspe, rmspe = sqrt(mspe))
  # update kdtf if nsim > 0
  if (nsim > 0) {
    kdtf$sim = newsim[-(1:n),] + (pred[,1] - pred[,-1])
  }

  # return kriging information in proper format
  if (return_type == "data.frame") {
    return(kdtf)
  } else {
    return(return_predict_geolm(kdtf, newcoords, return_type, object$coordnames))
  }
}

#' Check additional argument of predict.geolm_cmodMan
#'
#' @param y Vector of observed responses
#' @param vop Cross-covariance matrix
#' @param vp Covariance matrix of responses to be predicted
#' @noRd
arg_check_predict_geolm_cmodMan = function(y, vop, vp) {
  if (!is.matrix(vop)) {
    stop("vop must be a matrix")
  }
  if (length(dim(vop)) != 2) {
    stop("vop must be two-dimensional")
  }
  if (!is.matrix(vp)) {
    stop("vp must be a matrix")
  }
  if (length(dim(vp)) != 2) {
    stop("vp must be two-dimensional")
  }
  if (nrow(vp) != ncol(vp)) {
    stop("vp must be a square matrix")
  }
  if (nrow(vop) != length(y)) {
    stop("nrow(vop) != length(object$y)")
  }
  if (ncol(vop) != nrow(vp)) {
    stop("nrow(vop) != nrow(vp). They must match for compatibility.")
  }
}


