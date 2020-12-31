#' Argument check evar
#'
#' @param evar Variance of errors
#' @noRd
arg_check_evar = function(evar) {
  if (!is.numeric(evar)) {
    stop("evar must be numeric")
  }
  if (length(evar) != 1) {
    stop("evar must be a single value")
  }
  if (evar < 0) {
    stop("evar must be non-negative")
  }
  if (!is.finite(evar)) {
    stop("evar must be finite")
  }
}

#' Argument check split
#'
#' @param split Logical. Should plots be split?
#' @noRd
arg_check_split = function(split) {

  if (length(split) != 1) {
    stop("split must be a single value")
  }
  if (!is.logical(split)) {
    stop("split must be a logical value")
  }
}

#' Argument check add_legend
#'
#' @param add_legend Logical. Should we add_legend to plot?
#' @noRd
arg_check_add_legend = function(add_legend) {
  if (length(add_legend) != 1) {
    stop("add_legend must be a single value")
  }
  if (!is.logical(add_legend)) {
    stop("add_legend must be a logical value")
  }
}

#' Argument check partial sill parameter
#'
#' @param psill Partial sill
#' @noRd
arg_check_psill = function(psill) {
  if (length(psill) != 1) {
    stop("psill must be a single value")
  }
  if (!is.numeric(psill)) {
    stop("psill must be numeric")
  }
  if (!is.finite(psill)) {
    stop("psill must be a finite number")
  }
  if (psill <= 0) {
    stop("psill must be positive")
  }
}

#' Argument check range parameter
#'
#' @param r Range parameter
#' @noRd
arg_check_r = function(r) {
  if (length(r) != 1) {
    stop("r must be a single value")
  }
  if (!is.numeric(r)) {
    stop("r must be numeric")
  }
  if (!is.finite(r)) {
    stop("r must be a finite number")
  }
  if (r <= 0) {
    stop("r must be positive")
  }
}

#' Argument check fvar
#'
#' @param fvar Variance of errors
#' @noRd
arg_check_fvar = function(fvar) {
  if (!is.numeric(fvar)) {
    stop("fvar must be numeric")
  }
  if (length(fvar) != 1) {
    stop("fvar must be a single value")
  }
  if (fvar < 0) {
    stop("fvar must be non-negative")
  }
  if (!is.finite(fvar)) {
    stop("fvar must be finite")
  }
}

#' Argument check par3
#'
#' @param par3 3rd parameter
#' @noRd
arg_check_par3 = function(par3) {
  if (length(par3) != 1) {
    stop("par3 must be a single value")
  }
  if (!is.numeric(par3)) {
    stop("par3 must be numeric")
  }
  if (!is.finite(par3)) {
    stop("par3 must be a finite number")
  }
  if (par3 <= 0) {
    stop("par3 must be positive")
  }
}

#' Argument check longlat parameter
#'
#' @param longlat A single logical value
#' @noRd
arg_check_longlat = function(longlat) {
  if (length(longlat) != 1) {
    stop("longlat must be a single value")
  }
  if (!is.logical(longlat)) {
    stop("longlat must be a logical value")
  }
}

#' Argument check cmod_std
#'
#' @param model standard covariance model
#' @noRd
arg_check_cmod_std_model = function(model) {
  std_models = c("exponential",
                 "gaussian",
                 "matern",
                 "amatern",
                 "spherical",
                 "wendland1",
                 "wendland2",
                 "wu1",
                 "wu2",
                 "wu3")
  if (!is.element(model, std_models)) {
    stop("invalid model choice")
  }
}


#' Argument check cmod_std
#'
#' @param angle geometric anisotropy
#' @return NULL
#' @noRd
arg_check_angle = function(angle, radians) {
  if (length(angle) != 1) {
    stop("angle must be length 1")
  }
  if (!is.vector(angle)) {
    stop("angle must be a vector")
  }
  if (!is.numeric(angle)) {
    stop("angle must be numeric")
  }
  if (radians & (angle < 0 | angle >= pi)) {
    stop("angle must be in [0, pi) when radians = TRUE")
  }
  if (!radians & angle < 0 | angle >= 180) {
    stop("angle must be in [0, 180) when radians = FALSE")
  }
}

#' Argument check ratio
#'
#' @param ratio rmin/rmaj
#' @return NULL
#' @noRd
arg_check_ratio = function(ratio) {
  if (length(ratio) != 1) {
    stop("ratio must have length 1")
  }
  if (!is.vector(ratio)) {
    stop("ratio must be a vector")
  }
  if (!is.numeric(ratio)) {
    stop("ratio must be numeric")
  }
  if (ratio <= 0 | ratio > 1) {
    stop("ratio must be in (0, 1]")
  }
}

#' check distance object
#' @param d should be a numeric matrix
#' @noRd
arg_check_d = function(d) {
  if (!is.numeric(d)) {
    stop("d must be numeric")
  }
  if (!is.matrix(d)) {
    stop("d must be a matrix")
  }
}

#' check e argument of evaluate
#' @param e should be a single logical value
#' @noRd
arg_check_e = function(e) {
  if (length(e) != 1) {
    stop("e must be a single value")
  }
  if (!is.logical(e)) {
    stop("e must be a logical value")
  }
}

#' check f argument of evaluate
#' @param f should be a single logical value
#' @noRd
arg_check_f = function(f) {
  if (length(f) != 1) {
    stop("f must be a single value")
  }
  if (!is.logical(f)) {
    stop("f must be a logical value")
  }
}

#' Argument check formula
#'
#' Make sure the formula is valid
#'
#' @param formula A formula
#' @noRd
arg_check_formula = function(formula) {
  if (class(formula) != "formula") {
    stop("formula is not of class formula")
  }
  if (is.null(formula[[2]]) || length(formula[[2]]) > 2) {
    stop("formula should contain a single variable to the left of ~")
  }
}

#' Check verbose argument
#'
#' @param verbose A logical value indicating whether verbose
#' comments should be printed.
#' @noRd
arg_check_verbose = function(verbose) {
  if (length(verbose) != 1) {
    stop("verbose should be a single logical value")
  }
  if (!is.logical(verbose)) {
    stop("verbose should be a single logical value")
  }
}

#' Check mu argument
#'
#' @param mu A numeric value
#' @noRd
arg_check_mu = function(mu) {
  if (length(mu) != 1) {
    stop("mu must be a vector of length 1")
  }
  if (!is.numeric(mu)) {
    stop("mu must be a numeric value")
  }
  if (!is.finite(mu)) {
    stop("mu must be a finite value")
  }
}

#' Check formula argument for geolm
#'
#' @param formula A formula
#' @noRd
arg_check_geolm_formula = function(formula) {
  if (class(formula) != "formula" | is.null(formula))  {
    stop("formula must be a formula object or NULL")
  }
  if (length(formula) != 3) {
    stop("formula should have a single response to the left of the '~' and the covariates to the right")
  }
  if (length(formula[[2]]) != 1) {
    stop("There should only be a single response")
  }
}


#' Check weights argument
#'
#' @param weights A numeric vector with positive values
#' @param n The nrow(data)
#' @noRd
arg_check_weights = function(weights, n) {
  if (length(weights) != n) {
    stop("length(weights) != nrow(data)")
  }
  if (!is.numeric(weights)) {
    stop("weights must be numeric")
  }
  if (!is.vector(weights)) {
    stop("weights must be a vector")
  }
  if (min(weights) <= 0) {
    stop("all weights must be positive")
  }
  if (any(!is.finite(weights))) {
    stop("all weights must be finite")
  }
}

#' Check radians argument
#'
#' @param radians A logical value indicating whether the
#' angle is in radians
#' @noRd
arg_check_radians = function(radians) {
  if (length(radians) != 1) {
    stop("radians must be a single value")
  }
  if (!is.logical(radians)) {
    stop("radians must be logical")
  }
}

#' @param invert A logical value indicating whether the
#' angle is in invert
#' @noRd
arg_check_invert = function(invert) {
  if (length(invert) != 1) {
    stop("invert must be a single value")
  }
  if (!is.logical(invert)) {
    stop("invert must be logical")
  }
}

#' @param noise_type A character vector indicating whether
#' error or finescale variance should be estimated
#' @noRd
arg_check_noise_type = function(noise_type) {
  if (length(noise_type) != 1) {
    stop("noise_type must be a single value")
  }
  if (!is.character(noise_type)) {
    stop("noise_type must be type character")
  }
  if (!is.element(noise_type, c("e", "f"))) {
    stop("noise_type must be 'e' or 'f'")
  }
}

#' @param est_nugget A logical value indicating whether the
#' angle is in est_nugget
#' @noRd
arg_check_est_nugget = function(est_nugget) {
  if (length(est_nugget) != 1) {
    stop("est_nugget must be a single value")
  }
  if (!is.logical(est_nugget)) {
    stop("est_nugget must be logical")
  }
}

#' @param est_par3 A logical value indicating whether the
#' angle is in est_par3
#' @noRd
arg_check_est_par3 = function(est_par3) {
  if (length(est_par3) != 1) {
    stop("est_par3 must be a single value")
  }
  if (!is.logical(est_par3)) {
    stop("est_par3 must be logical")
  }
}

#' @param est_angle A logical value indicating whether the
#' angle is in est_angle
#' @noRd
arg_check_est_angle = function(est_angle) {
  if (length(est_angle) != 1) {
    stop("est_angle must be a single value")
  }
  if (!is.logical(est_angle)) {
    stop("est_angle must be logical")
  }
}

#' @param est_ratio A logical value indicating whether the
#' angle is in est_ratio
#' @noRd
arg_check_est_ratio = function(est_ratio) {
  if (length(est_ratio) != 1) {
    stop("est_ratio must be a single value")
  }
  if (!is.logical(est_ratio)) {
    stop("est_ratio must be logical")
  }
}

arg_check_lambda = function(lambda) {
  if (length(lambda) != 1) {
    stop("lambda must have length 1")
  }
  if (lambda < 0) {
    stop("lambda must be >= 0")
  }
}

arg_check_angle = function(angle, radians) {
  if (length(angle) != 1) {
    stop("angle must have length 1")
  }
  if (angle < 0) {
    stop("angle must be >= 0")
  }
  if (radians & angle > pi) {
    stop("angle (radians) must be < pi")
  }
  if (!radians & angle > 180) {
    stop("angle (degrees) must be < 180")
  }
}

arg_check_ratio = function(ratio) {
  if (length(ratio) != 1) {
    stop("ratio must have length 1")
  }
  if (ratio <= 0) {
    stop("ratio must be > 0")
  }
  if (ratio > 1) {
    stop("ratio must be <= 1")
  }
}

#' Argument check reml parameter
#'
#' @param reml A single logical value
#' @noRd
arg_check_reml = function(reml) {
  if (length(reml) != 1) {
    stop("reml must be a single value")
  }
  if (!is.logical(reml)) {
    stop("reml must be a logical value")
  }
}

#' Argument check sp parameter
#'
#' @param sp A single logical value
#' @noRd
arg_check_sp = function(sp) {
  if (length(sp) != 1) {
    stop("sp must be a single value")
  }
  if (!is.logical(sp)) {
    stop("sp must be a logical value")
  }
}


#' Argument check dmethod argument
#'
#' @param dmethod A single character
#' @noRd
arg_check_dmethod = function(dmethod) {
  if (length(dmethod) != 1) {
    stop("dmethod must be a single value")
  }
  if (!is.vector(dmethod)) {
    stop("dmethod must be a vector of length 1")
  }
  if (!is.character(dmethod)) {
    stop("dmethod must be of character type")
  }
  if (!is.element(dmethod, c("chol", "eigen", "svd"))) {
    stop("dmethod must be 'chol', 'eigen', or 'svd'")
  }
}

#' Argument check return_type argument
#'
#' @param return_type A single character
#' @noRd
arg_check_return_type = function(return_type) {
  if (length(return_type) != 1) {
    stop("return_type must be a single value")
  }
  if (!is.vector(return_type)) {
    stop("return_type must be a vector of length 1")
  }
  if (!is.character(return_type)) {
    stop("return_type must be of character type")
  }
  if (!is.element(return_type, c("data.frame", "gearPredict", "geardf", "SpatialPointsDataFrame", "sf"))) {
    stop("return_type must be 'data.frame', 'gearPredict', 'geardf', 'SpatialPointsDataFrame', or 'sf'")
  }
}

#' Check nsim argument
#'
#' @param nsim A non-negative integer
#' @noRd
arg_check_nsim = function(nsim) {
  if (length(nsim) != 1) {
    stop("nsim must have length 1")
  }
  if (!is.vector(nsim)) {
    stop("nsim must be a vector of length 1")
  }
  if (!is.numeric(nsim)) {
    stop("nsim must be numeric")
  }
  if (!is.finite(nsim)) {
    stop("nsim must be finite")
  }
  if (nsim < 0) {
    stop("nsim must be non-negative (and preferably an integer)")
  }
}

#' Check newdata argument
#'
#' @param nsim A two-d data.frame
#' @param coordnames coordinate names that must be in newdata
#' @noRd
arg_check_newdata = function(newdata, coordnames) {
  if (!is.data.frame(newdata)) {
    stop("newdata should be a data.frame")
  }
  df_names = names(newdata)
  if (min(is.element(coordnames, df_names)) == 0) {
    stop("The coordnames are not found in the column names of newdata.")
  }
}

#' Argument check compute_mspe argument
#'
#' @param compute_mspe A single logical value
#' @noRd
arg_check_compute_mspe = function(compute_mspe) {
  if (length(compute_mspe) != 1) {
    stop("compute_mspe must be a single value")
  }
  if (!is.logical(compute_mspe)) {
    stop("compute_mspe must be a logical value")
  }
}

#' Check coordnames argument and return vector of coordinate names
#'
#' @param coordnames A vector with the names of the coordinates
#' @param data_colnames A vector with the column names of data
#' @return A vector of coordinate names
#' @noRd
arg_check_coordnames = function(coordnames, data_colnames) {
  if (!(class(coordnames)[1] == "formula") &
      !is.character(coordnames) &
      !is.numeric(coordnames)) {
    error_message = paste("coordnames must indicate the coordinate columns of data through:\n",
                          "1. a formula (e.g., ~ lon + lat)\n",
                          "2. a character vector of column names (e.g., c('lon', 'lat'))\n",
                          "3. a numeric vector of column indices (e.g., 3:4)", sep = "")
    stop(error_message)
  }

  # get correct columns if coordnames is a formula
  if (class(coordnames) == "formula") {
    coordnames = labels(stats::terms(coordnames))
  }
  if (length(coordnames) != 2) {
    stop("coordnames must refer to only two columns of data")
  }
  if (is.numeric(coordnames)) {
    if (min(coordnames) < 1) {
      stop("The indices of numeric of coordnames must be >= 1")
    }
    coordnames = data_colnames[coordnames]
    if (max(is.na(coordnames)) == 1) {
      error_message = paste("The indices specified in coordnames appear to be invalid.\n",
                            "The command names(data)[coordnames] returns: ",
                            coordnames[1], " ", coordnames[2], sep = "")
      stop(error_message)
    }
  }
  if (is.character(coordnames)) {
    if (min(is.element(coordnames, data_colnames)) == 0) {
      stop("The names in coordnames are not found in the column names of data")
    }
  }
  return(coordnames)
}

