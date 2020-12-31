#' Choose the upper argument for ploglik_* function for estimate.geolm_cmodStd
#'
#' @param est_nugget Logical values indicating which
#' parameters should be estimated
#' @param est_par3 
#' @param est_angle 
#' @param est_ratio 
#' @param upper bounds for relevant parameters
#' @noRd
assign_upper = function(est_nugget, est_par3, est_angle,
                        est_ratio, upper, radians, maxd, vary) {
  r = ifelse(is.null(upper$r), 5 * maxd, upper$r)
  arg_check_r(r)
  psill = ifelse(is.null(upper$psill), 5 * vary, upper$psill)
  arg_check_psill(psill)
  lambda = ifelse(is.null(upper$lambda), 10, upper$lambda)
  arg_check_lambda(lambda)
  angle = ifelse(is.null(upper$angle), ifelse(radians, pi - 0.001, (pi - 0.001)/pi * 180), upper$angle)
  arg_check_angle(angle, radians)
  ratio = ifelse(is.null(upper$ratio), 1, upper$ratio)
  arg_check_ratio(ratio)
  par3 = ifelse(is.null(upper$par3), 3, upper$par3)
  arg_check_par3(par3)
  name_angle = ifelse(radians, "angle (radians)", "angle (degrees)")
  
  if (!est_nugget & !est_par3 & !est_angle & !est_ratio) {
    upper = c(r, psill)
    names(upper) = c("r", "psill")
  } else if (!est_nugget & est_par3 & !est_angle & !est_ratio) {
    upper = c(r, psill, par3)
    names(upper) = c("r", "psill", "par3")
  } else if (!est_nugget & !est_par3 & est_angle & !est_ratio) {
    upper = c(r, psill, angle)
    names(upper) = c("r", "psill", name_angle)
  } else if (!est_nugget & !est_par3 & !est_angle & est_ratio) {
    upper = c(r, psill, ratio)
    names(upper) = c("r", "psill", "ratio (rminor/r)")
  } else if (!est_nugget & est_par3 & est_angle & !est_ratio) {
    upper = c(r, psill, angle, par3)
    names(upper) = c("r", "psill", name_angle, "par3")
  } else if (!est_nugget & est_par3 & !est_angle & est_ratio) {
    upper = c(r, psill, ratio, par3)
    names(upper) = c("r", "psill", "ratio (rminor/r)", "par3")
  } else if (!est_nugget & !est_par3 & est_angle & est_ratio) {
    upper = c(r, psill, angle, ratio)
    names(upper) = c("r", "psill", name_angle, "ratio (rminor/r)")
  } else if (!est_nugget & est_par3 & est_angle & est_ratio) {
    upper = c(r, psill, angle, ratio, par3)
    names(upper) = c("r", "psill", name_angle, "ratio (rminor/r)", "par3")
  } else if (est_nugget & !est_par3 & !est_angle & !est_ratio) {
    upper = c(r, lambda)
    names(upper) = c("r", "lambda")
  } else if (est_nugget & est_par3 & !est_angle & !est_ratio) {
    upper = c(r, lambda, par3)
    names(upper) = c("r", "lambda", "par3")
  } else if (est_nugget & !est_par3 & est_angle & !est_ratio) {
    upper = c(r, lambda, angle)
    names(upper) = c("r", "lambda", name_angle)
  } else if (est_nugget & !est_par3 & !est_angle & est_ratio) {
    upper = c(r, lambda, ratio)
    names(upper) = c("r", "lambda", "ratio (rminor/r)")
  } else if (est_nugget & est_par3 & est_angle & !est_ratio) {
    upper = c(r, lambda, angle, par3)
    names(upper) = c("r", "lambda", name_angle, "par3")
  } else if (est_nugget & est_par3 & !est_angle & est_ratio) {
    upper = c(r, lambda, ratio, par3)
    names(upper) = c("r", "lambda", "ratio (rminor/r)", "par3")
  } else if (est_nugget & !est_par3 & est_angle & est_ratio) {
    upper = c(r, lambda, angle, ratio)
    names(upper) = c("r", "lambda", name_angle, "ratio (rminor/r)")
  } else if (est_nugget & est_par3 & est_angle & est_ratio) {
    upper = c(r, lambda, angle, ratio, par3)
    names(upper) = c("r", "lambda", name_angle, "ratio (rminor/r)", "par3")
  } 
  return(upper)
}
