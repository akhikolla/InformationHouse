#' Choose the ploglik_* function for estimate.geolm_cmodStd
#'
#' @param est_nugget Logical values indicating which
#' parameters should be estimated
#' @param est_par3 
#' @param est_angle 
#' @param est_ratio 
#' @noRd
choose_ploglik = function(est_nugget, est_par3, est_angle,
                          est_ratio) {
  if (!est_nugget & !est_par3 & !est_angle & !est_ratio) {
    pll = ploglik_cmodStd_r_psill
  } else if (!est_nugget & est_par3 & !est_angle & !est_ratio) {
    pll = ploglik_cmodStd_r_psill_par3
  } else if (!est_nugget & !est_par3 & est_angle & !est_ratio) {
    pll = ploglik_cmodStd_r_psill_angle
  } else if (!est_nugget & !est_par3 & !est_angle & est_ratio) {
    pll = ploglik_cmodStd_r_psill_ratio
  } else if (!est_nugget & est_par3 & est_angle & !est_ratio) {
    pll = ploglik_cmodStd_r_psill_angle_par3
  } else if (!est_nugget & est_par3 & !est_angle & est_ratio) {
    pll = ploglik_cmodStd_r_psill_ratio_par3
  } else if (!est_nugget & !est_par3 & est_angle & est_ratio) {
    pll = ploglik_cmodStd_r_psill_angle_ratio
  } else if (!est_nugget & est_par3 & est_angle & est_ratio) {
    pll = ploglik_cmodStd_r_psill_angle_ratio_par3
  } else if (est_nugget & !est_par3 & !est_angle & !est_ratio) {
    pll = ploglik_cmodStd_r_lambda
  } else if (est_nugget & est_par3 & !est_angle & !est_ratio) {
    pll = ploglik_cmodStd_r_lambda_par3
  } else if (est_nugget & !est_par3 & est_angle & !est_ratio) {
    pll = ploglik_cmodStd_r_lambda_angle
  } else if (est_nugget & !est_par3 & !est_angle & est_ratio) {
    pll = ploglik_cmodStd_r_lambda_ratio
  } else if (est_nugget & est_par3 & est_angle & !est_ratio) {
    pll = ploglik_cmodStd_r_lambda_angle_par3
  } else if (est_nugget & est_par3 & !est_angle & est_ratio) {
    pll = ploglik_cmodStd_r_lambda_ratio_par3
  } else if (est_nugget & !est_par3 & est_angle & est_ratio) {
    pll = ploglik_cmodStd_r_lambda_angle_ratio
  } else if (est_nugget & est_par3 & est_angle & est_ratio) {
    pll = ploglik_cmodStd_r_lambda_angle_ratio_par3
  }
  return(pll)
}
