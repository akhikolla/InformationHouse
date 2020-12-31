#' Check arguments of predict_geolm*
#'
#' @param coordnames Vector of coordinate names
#' @param newdata A data.frame with newdata for predictions
#' @param nsim Integer of number of conditional simulation
#' @param return_type Indicate return type (sp, geardf, sf, data.frame)
#' @param dmethod Decomposition method
#' @param mspe Logical indicating whether mspe should be returned
#' @noRd
arg_check_predict_geolm = function(coordnames, newdata, nsim, return_type, dmethod, compute_mspe) {
  arg_check_newdata(newdata, coordnames)
  arg_check_nsim(nsim)
  arg_check_return_type(return_type)
  arg_check_dmethod(dmethod)
  arg_check_compute_mspe(compute_mspe)
}
