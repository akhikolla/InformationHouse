#' Update cmodStd based as optimization
#'
#' @param optimx_out Optimization output
#' @param opt_fn The function used for optimization
#' @param object The geolm_cmodStd object
#' @param d The distance object
#' @param nugget The total nugget
#' @param scmod Scaled covariance model
#' @param reml Whether estimate is reml
#' @param noise_type Error or finescale variance
#' @return cmodStd object
#' @noRd
update_mod = function(optimx_out, opt_fn, object, d, nugget, scmod, reml, noise_type) {
  if (!is.null(optimx_out$lambda)) {
    opt_par = c(optimx_out$r, optimx_out$lambda)
    if (!is.null(optimx_out$angle)) {
      opt_par = c(opt_par, optimx_out$angle)
    }
    if (!is.null(optimx_out$ratio)) {
      opt_par = c(opt_par, optimx_out$ratio)
    }
    if (!is.null(optimx_out$par3)) {
      opt_par = c(opt_par, optimx_out$par3)
    }
    psill = do.call(opt_fn, 
                    list(par = opt_par, x = object$x,
                         y = object$y, d = d,
                         weights = object$weights, 
                         scmod = scmod, nugget = nugget,
                         mu = object$mu, reml = reml,
                         return_ll = FALSE))
    scmod$psill = psill
    if (noise_type == "e") {
      scmod$evar = psill * optimx_out$lambda
    } else {
      scmod$fvar = psill * optimx_out$lambda
    }
  } else {
    scmod$psill = optimx_out$psill
    if (noise_type == "e") {
      scmod$evar = nugget
    } else {
      scmod$fvar = nugget
    }
  }
  scmod$r = optimx_out$r
  if (!is.null(optimx_out$angle)) {
    scmod$angle = optimx_out$angle
  }
  if (!is.null(optimx_out$ratio)) {
    scmod$ratio = optimx_out$ratio
  }
  if (!is.null(optimx_out$par3)) {
    scmod$par3 = optimx_out$par3
  }
  return(scmod)
}
