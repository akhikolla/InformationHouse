#' Choose the lower argument for ploglik_* function for estimate.geolm_cmodStd
#'
#' @param est_nugget Logical values indicating which
#' parameters should be estimated
#' @param est_par3 
#' @param est_angle 
#' @param est_ratio 
#' @param lr etc lower bounds for relevant parameters
#' @noRd
assign_lower = function(est_nugget, est_par3, est_angle,
                        est_ratio, lower, radians) {
  r = ifelse(is.null(lower$r), 0.001, lower$r)
  arg_check_r(r)
  psill = ifelse(is.null(lower$psill), 0.001, lower$psill)
  arg_check_psill(psill)
  lambda = ifelse(is.null(lower$lambda), 0, lower$lambda)
  arg_check_lambda(lambda)
  angle = ifelse(is.null(lower$angle), 0, lower$angle)
  arg_check_angle(angle, radians)
  ratio = ifelse(is.null(lower$ratio), 0.001, lower$ratio)
  arg_check_ratio(ratio)
  par3 = ifelse(is.null(lower$par3), 0.001, lower$par3)
  arg_check_par3(par3)
  name_angle = ifelse(radians, "angle (radians)", "angle (degrees)")
  
  if (!est_nugget & !est_par3 & !est_angle & !est_ratio) {
    lower = c(r, psill)
    names(lower) = c("r", "psill")
  } else if (!est_nugget & est_par3 & !est_angle & !est_ratio) {
    lower = c(r, psill, par3)
    names(lower) = c("r", "psill", "par3")
  } else if (!est_nugget & !est_par3 & est_angle & !est_ratio) {
    lower = c(r, psill, angle)
    names(lower) = c("r", "psill", name_angle)
  } else if (!est_nugget & !est_par3 & !est_angle & est_ratio) {
    lower = c(r, psill, ratio)
    names(lower) = c("r", "psill", "ratio (rminor/r)")
  } else if (!est_nugget & est_par3 & est_angle & !est_ratio) {
    lower = c(r, psill, angle, par3)
    names(lower) = c("r", "psill", name_angle, "par3")
  } else if (!est_nugget & est_par3 & !est_angle & est_ratio) {
    lower = c(r, psill, ratio, par3)
    names(lower) = c("r", "psill", "ratio (rminor/r)", "par3")
  } else if (!est_nugget & !est_par3 & est_angle & est_ratio) {
    lower = c(r, psill, angle, ratio)
    names(lower) = c("r", "psill", name_angle, "ratio (rminor/r)")
  } else if (!est_nugget & est_par3 & est_angle & est_ratio) {
    lower = c(r, psill, angle, ratio, par3)
    names(lower) = c("r", "psill", name_angle, "ratio (rminor/r)", "par3")
  } else if (est_nugget & !est_par3 & !est_angle & !est_ratio) {
    lower = c(r, lambda)
    names(lower) = c("r", "lambda")
  } else if (est_nugget & est_par3 & !est_angle & !est_ratio) {
    lower = c(r, lambda, par3)
    names(lower) = c("r", "lambda", "par3")
  } else if (est_nugget & !est_par3 & est_angle & !est_ratio) {
    lower = c(r, lambda, angle)
    names(lower) = c("r", "lambda", name_angle)
  } else if (est_nugget & !est_par3 & !est_angle & est_ratio) {
    lower = c(r, lambda, ratio)
    names(lower) = c("r", "lambda", "ratio (rminor/r)")
  } else if (est_nugget & est_par3 & est_angle & !est_ratio) {
    lower = c(r, lambda, angle, par3)
    names(lower) = c("r", "lambda", name_angle, "par3")
  } else if (est_nugget & est_par3 & !est_angle & est_ratio) {
    lower = c(r, lambda, ratio, par3)
    names(lower) = c("r", "lambda", "ratio (rminor/r)", "par3")
  } else if (est_nugget & !est_par3 & est_angle & est_ratio) {
    lower = c(r, lambda, angle, ratio)
    names(lower) = c("r", "lambda", name_angle, "ratio (rminor/r)")
  } else if (est_nugget & est_par3 & est_angle & est_ratio) {
    lower = c(r, lambda, angle, ratio, par3)
    names(lower) = c("r", "lambda", name_angle, "ratio (rminor/r)", "par3")
  } 
  return(lower)
}
