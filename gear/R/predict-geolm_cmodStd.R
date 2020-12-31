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
#' @inheritParams predict.geolm_cmodMan
#' @param ... Currently unimplemented.
#' @inherit predict.geolm_cmodMan return
#' @author Joshua French
#' @examples
#' # generate response
#' y = rnorm(10)
#' # generate coordinates
#' x1 = runif(10); x2 = runif(10)
#'
#' # data frame for observed data
#' data = data.frame(y, x1, x2)
#' # newdata must have columns with prediction coordinates
#' newdata = data.frame(x1 = runif(5), x2 = runif(5))
#'
#' # specify a standard covariance model
#' mod = cmod_std(model = "exponential", psill = 1, r = 1)
#'
#' # geolm for universal kriging
#' gearmod_uk = geolm(y ~ x1 + x2, data = data, mod = mod,
#'                    coordnames = c("x1", "x2"))
#' # prediction for universal kriging, with conditional simulation
#' pred_uk = predict(gearmod_uk, newdata, nsim = 2)
#'
#' # demonstrate plotting abilities if return_type == "geardf"
#'  pred_geardf = predict(gearmod_uk, newdata,
#'                return_type = "geardf")
#'  plot(pred_geardf, "pred")
#'  plot(pred_geardf, interp = TRUE)
#'
#' # demonstrate plotting abilities if sp package installed
#' if (requireNamespace("sp", quietly = TRUE)) {
#'  pred_spdf = predict(gearmod_uk, newdata,
#'              return_type = "SpatialPointsDataFrame")
#'  sp::spplot(pred_spdf, "pred")
#'  sp::spplot(pred_spdf)
#' }
#' # demonstrate plotting abilities if sf package installed
#' if (requireNamespace("sf", quietly = TRUE)) {
#'  pred_sfdf = predict(gearmod_uk, newdata,
#'              return_type = "sf")
#'  plot(pred_sfdf["pred"])
#'  plot(pred_sfdf)
#' }
#'
#' # geolm for ordinary kriging
#' gearmod_ok = geolm(y ~ 1, data = data, mod = mod,
#'                    coordnames = c("x1", "x2"))
#'# prediction for ordinary kriging
#' pred_ok = predict(gearmod_ok, newdata)
#'
#' # geolm for simple kriging
#' gearmod_sk = geolm(y ~ 1, data = data, mod = mod,
#'                  coordnames = c("x1", "x2"), mu = 1)
#'# prediction for simple kriging
#' pred_sk = predict(gearmod_sk, newdata)
#'@rdname predict.geolm_cmodStd
#'@export
predict.geolm_cmodStd = function(object, newdata, nsim = 0,
                            return_type = "SpatialPointsDataFrame",
                            dmethod = "chol",
                            compute_mspe = TRUE, sp = NULL, ...) {
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
  newcoords = as.matrix(newdata[, object$coordnames])
  if (object$mod$ratio < 1) {
    vop = evaluate(mod = object$mod, e = FALSE,
                   d = ganiso_d(object$coords, newcoords,
                                invert = object$mod$invert))
  } else {
    vop = evaluate(mod = object$mod, e = FALSE,
                   d = geodist(object$coords, newcoords,
                               longlat = object$mod$longlat))
  }

  # unneeded if simple kriging
  newx = NULL
  if (is.null(object$mu)) {
    newf = stats::delete.response(stats::terms(object$formula))
    newx = stats::model.matrix(newf, data = newdata)
    dimnames(newx) = NULL
  }

  # generate simulated data at observed and prediction locations
  if (nsim > 0) {
    if (object$mod$ratio < 1) {
      vp = evaluate(mod = object$mod, e = FALSE,
                    d = ganiso_d(newcoords,
                                 invert = object$mod$invert))
    } else {
      vp = evaluate(mod = object$mod, e = FALSE,
                    d = geodist(newcoords, longlat = object$mod$longlat))
    }

    n = nrow(object$coords); m = nrow(newcoords)
    # correct mean if simple kriging
    mu = 0
    if (!is.null(object$mu)) mu = object$mu

    newsim = decomp_cov(rbind(
      cbind(object$v - diag(object$vediag),vop),
      cbind(t(vop), vp)), method = dmethod) %*%
      matrix(stats::rnorm((n + m)*nsim), ncol = nsim) +
      mu

    # slightly different algorithm for simple kriging
    # update various object to include simulated data
    object$y = cbind(object$y, newsim[1:n,] +
                matrix(stats::rnorm(n * nsim,
                       sd = sqrt(object$vediag)),
                       nrow = n, ncol = nsim))

    if (!is.null(object$mu)) {
      object$cholviresid = backsolve(object$cholv, object$y - mu, transpose = TRUE)
    } else {
      object$coeff = object$map2coeff %*% object$y
      object$cholviresid = backsolve(object$cholv, object$y - object$x %*% object$coeff, transpose = TRUE)
    }
  }

  cholvivop = forwardsolve(object$cholv, vop, transpose = TRUE, upper.tri = TRUE)

  # predictions
  pred = fitted(object, newx) + crossprod(cholvivop, object$cholviresid)

  # assume mspe is not computed
  mspe = NA
  if (compute_mspe) {
    #compute mspe
    # mspe1 = object$mod$psill + object$mod$fvar
    # mspe2 = colSums(cholvivop^2)
    mspe = object$mod$psill + object$mod$fvar - colSums(cholvivop^2)
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

