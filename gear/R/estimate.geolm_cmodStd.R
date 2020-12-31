#' Determine MLEs of model parameters for a geostatistical
#' model
#'
#' \code{estimate} estimates the parameters of a
#' geostatistical linear model of class \code{geolm_cmodStd}
#' using maximum likelihood estimation.
#'
#' The \code{\link[optimx]{optimx}} function is used to find
#' the MLEs.  The \code{control} argument of
#' \code{\link[optimx]{optimx}} has a parameter \code{kkt}
#' related to checking optimality conditions. This is
#' internally set to \code{FALSE}.  See
#' \code{\link[optimx]{optimx}} for Details.
#'
#' Only the sum of \code{evar} and \code{fvar} is
#' identifiable. Depending on the choice of
#' \code{noise_type}, the covariance model is internally
#' updated to estimate only one type of noise. e.g., if
#' \code{noise_type = "e"}, then internally we update
#' \code{evar} so that \code{evar = evar + fvar} and
#' \code{fvar = 0}. Estimation is then performed on
#' \code{evar} alone. Alternatively, the analagous estimated
#' would be made for \code{fvar} if \code{noise_type =
#' "fvar"}.
#'
#' When \code{est_nugget} is true, the likelihood is
#' profiled to simplify the optimization problem. In that
#' case a parameter \code{lambda = (evar + fvar)/psill} is
#' optimized. The optimal \code{psill} and noise variance
#' are then determined.
#'
#' The \code{lower} argument should be a named list with the
#' names of the parameters you wish to set lower bounds for.
#' If not specified, an attempt is made to specify
#' reasonable lower bounds. The current choices are \code{r
#' = 0.001}, \code{psill = 0.001}, \code{lambda = 0},
#' \code{angle = 0}, \code{ratio = 0.001}, \code{par3 =
#' 0.001}.
#'
#' The \code{upper} argument should be a named list with the
#' names of the parameters you wish to set upper bounds for.
#' If not specified, an attempt is made to specify
#' reasonable upper bounds. The current choices are \code{r
#' = 5 * } maximum intercentroid distance, \code{psill = 5 *
#' var(object$y)}, \code{lambda = 5}, \code{angle = 179.99},
#' \code{ratio = 1}, \code{par3 = 3}.
#'
#' @param object A geostatistical linear model object
#'   produced by the \code{geolm} function.
#' @param reml  A logical value indicating whether standard
#'   maximum likelihood estimation should be performed
#'   (\code{reml = FALSE}).  If \code{reml = TRUE}, then
#'   restricted maximum likelihood estimation is performed.  The
#'   default is \code{FALSE}.
#' @param noise_type A character vector indicating the type
#'   of noise (nugget) variance to estimate. The default is
#'   \code{est = "e"}, indicating the error variance should
#'   be estimated. Alternatively, the user can specify
#'   \code{est = "f"}, indicating that the finescale (microscale)
#'   variance should be estimated. The other type of noise
#'   variance is set to 0, otherwise the model is not
#'   identifiable. See Details.
#' @param lower A named list with the names of the
#' parameters you wish to set lower bounds for and the
#' associated value. See Details.
#' @param upper A named list with the names of the
#' parameters you wish to set upper bounds for and the
#' associated value.
#' @param method The optimization method.  The default is
#'   \code{"nlminb"}. \code{"L-BFGS-B"} is another
#'   acceptable choice.  See \code{\link[optimx]{optimx}}
#'   for further choices.
#' @param itnmax An integer indicating the maximum number of
#'   iterations to allow for the optimization procedure.
#' @param control A list of control parameters passed
#'   internally to \code{\link[optimx]{optimx}}.
#' @param est_nugget A logical value indicating whether the
#'   nugget variance (\code{evar} or \code{fvar}) should be
#'   estimated. The default is \code{TRUE}, indicating that
#'   the nugget should be estimated.
#' @param est_par3 A logical value indicating whether
#'   \code{par3} should be estimated (for an appropriate
#'   covariance model such as \code{"matern"} or
#'   \code{"amatern"}). The default is \code{TRUE},
#'   indicating that this parameter should be estimated.
#' @param est_angle A logical value indicating whether the
#'   geometric anisotropy angle should be estimated. The
#'   default is \code{FALSE}, indicating that this parameter
#'   should not be estimated.
#' @param est_ratio A logical value indicating whether the
#'   geometric anisotropy ratio of minor axis length to
#'   major axis length should be estimated. The default is
#'   \code{FALSE}, indicating that this parameter should not
#'   be estimated.
#' @param verbose A logical value indicating whether
#' potentially informative messages should be printed. The
#' default is \code{FALSE}.
#' @param ... Currently unimplemented
#'
#' @author Joshua French
#' @seealso \code{\link[gear]{cmod_std}},
#'   \code{\link[optimx]{optimx}}
#' @export
#' @rdname estimate.geolm_cmodStd
#' @examples
#' data(toydata, package = "gear")
#' # setup standard covariance model
#' mod_std = cmod_std("exponential", psill = 1, r = 1, evar = 0.1)
#' # setup  dataframe with data
#' # fit Std geolm
#' object = geolm(y ~ x1 + x2, data = toydata, mod = mod_std,
#'                   coordnames = c("x1", "x2"))
#' est_object = estimate(object, control = list(trace = 1),
#'                       verbose = TRUE,
#'                       lower = list(r = 0.05, lambda = 0.05))
estimate.geolm_cmodStd = function(object, reml = FALSE, noise_type = "e",
                        lower = NULL, upper = NULL,
                        method = "nlminb", itnmax = NULL,
                        control = list(), est_nugget = TRUE,
                        est_par3 = TRUE, est_angle = FALSE,
                        est_ratio = FALSE, verbose = FALSE,
                        ...) {
  arg_check_estimate_geolm_cmodStd(reml = reml,
                                   noise_type = noise_type,
                                   est_nugget = est_nugget,
                                   est_par3 = est_par3,
                                   est_angle = est_angle,
                                   est_ratio = est_ratio,
                                   control = control,
                                   verbose = verbose)
  # scaled covariance model for optimization
  scmod = object$mod
  nugget = scmod$evar + scmod$fvar
  scmod$evar = 0
  scmod$fvar = 0
  if (est_nugget) scmod$psill = 1
  control$kkt = FALSE

  if (est_angle || est_ratio || object$mod$ratio < 1) {
    d = ganiso_d(object$coords, object$coords,
                 radians = TRUE,
                 object$mod$invert)
    maxd = max(d$d)
  } else {
    d = geodist(as.matrix(object$coords), longlat = object$mod$longlat)
    maxd = max(d)
  }
  vary = stats::var(object$y)

  lower = assign_lower(est_nugget = est_nugget,
                       est_par3 = est_par3,
                       est_angle = est_angle,
                       est_ratio = est_ratio,
                       lower = lower,
                       radians = object$mod$radians)
  upper = assign_upper(est_nugget = est_nugget,
                       est_par3 = est_par3,
                       est_angle = est_angle,
                       est_ratio = est_ratio,
                       upper = upper,
                       radians = object$mod$radians,
                       maxd = maxd, vary = vary)
  parm = assign_parm(est_nugget = est_nugget,
                     est_par3 = est_par3,
                     est_angle = est_angle,
                     est_ratio = est_ratio,
                     mod = scmod, nugget = nugget)
  arg_check_lower_parm_upper(lower = lower, parm = parm,
                             upper = upper)
  if (verbose) {
    message("Assigned parameter values:")
    combine = rbind(lower, parm, upper)
    rownames(combine) = c("lower", "initial", "upper")
    colnames(combine) = names(lower)
    print(combine)
  }

  if (est_angle & !object$mod$radians) {
    if (verbose) message("Converting degrees to radians")
    lower[3] = lower[3] * pi / 180
    parm[3] = parm[3] * pi / 180
    upper[3] = upper[3] * pi / 180
    scmod$angle = scmod$angle * pi/180
    scmod$radians = TRUE
  }
  # choose optimization function
  opt_fn = choose_ploglik(est_nugget, est_par3, est_angle, est_ratio)

  out_call = list(parm = parm, fn = opt_fn,
                  lower = lower, upper = upper,
                  x = object$x, y = object$y, d = d,
                  weights = object$weights, scmod = scmod,
                  nugget = nugget, mu = object$mu,
                  reml = reml, return_ll = TRUE,
                  method = method, itnmax = itnmax,
                  control = control)

  # optimize parameters
  out = optimx::optimx(parm, fn = opt_fn,
                 lower = lower, upper = upper,
                 x = object$x, y = object$y, d = d,
                 weights = object$weights, scmod = scmod,
                 nugget = nugget, mu = object$mu,
                 reml = reml, return_ll = TRUE,
                 method = method, itnmax = itnmax,
                 control = control)

  if (est_angle & !object$mod$radians) {
    if (verbose) message("Converting radians to degrees")
    out$angle = out$angle * 180/pi
    scmod$angle = scmod$angle * 180/pi
    scmod$radians = FALSE
  }

  # update cmodStd
  new_mod = update_mod(out, opt_fn, object = object, d = d,
                       nugget = nugget, scmod = scmod,
                       reml = reml, noise_type = noise_type)

  if (est_angle & !scmod$radians) {
    out$angle = out$angle * pi/180 # convert back to radians for comparison
  }


  object = update(object, new_mod)
  object$optimx = out
  object$out_call = out_call
  # object$new_mod = new_mod
  return(object)
}

#' Argument check estimate.geolm_cmodStd
#' @param reml  A logical value indicating whether standard
#'   maximum likelihood estimation should be performed
#'   (\code{reml = FALSE}).  If \code{reml = TRUE}, then
#'   restricted maximum likelihood is performed.  Default is
#'   \code{FALSE}.
#' @param noise_type A character vector indicating the type
#'   of noise (nugget) variance to estimate. The default is
#'   (\code{est = "e"}), indicating the error variance
#'   should be estimated. Alternatively, the user can
#'   specify (\code{est = "f"}), indicating the finescale
#'   variance should be estimated. The other type of noise
#'   variance is set to 0, otherwise the model is not
#'   identifiable.
#' @param est_nugget A logical value indicating whether the
#'   nugget variance (\code{evar} or \code{fvar}) should be
#'   estimated. The default is \code{TRUE}.
#' @param est_par3 A logical value indicating whether
#'   \code{par3} should be estimated (for an appropriate
#'   covariance model such ash \code{"matern"} or
#'   \code{"amatern"}. The default is \code{TRUE}.
#' @param est_angle A logical value indicating whether the
#'   geometric anisotropy angle should be estimated. The
#'   default is \code{FALSE}. This argument is ignored of
#'   \code{mod$ganiso} is \code{NULL}.
#' @param est_ratio A logical value indicating whether the
#'   geometric anisotropy ratio of minor axis length to
#'   major axis length should be estimated. The default is
#'   \code{FALSE}. This argument is ignored of
#'   \code{mod$ganiso} is \code{NULL}.
#' @noRd
arg_check_estimate_geolm_cmodStd = function(reml, noise_type,
                                            est_nugget,
                                            est_par3,
                                            est_angle,
                                            est_ratio,
                                            control,
                                            verbose) {
  arg_check_reml(reml)
  arg_check_noise_type(noise_type)
  arg_check_est_nugget(est_nugget)
  arg_check_est_par3(est_par3)
  arg_check_est_angle(est_angle)
  arg_check_est_ratio(est_ratio)
  if (!is.list(control)) stop("control must be a list")
  arg_check_verbose(verbose)
}

#' Check lower, upper bounds for parm for estimate.geolm_cmodStd
#'
#' @param lower A vector of lower bounds for relevant parameters
#' @param parm A vector of initial starting values for relevant parameters
#' @param upper A vector of upper bounds for relevant parameters
#' @noRd
arg_check_lower_parm_upper = function(lower, parm, upper, mod) {
  if (length(lower) != length(parm)) {
    message("The length of lower must must match the number of parameters to be estimated.")
    message("lower:")
    print(lower)
    message("initial values:")
    names(parm) = names(lower)
    print(parm)
    stop("Please specify lower to have the same length as the number of parameters to be estimated.")
  }
  if (length(upper) != length(parm)) {
    message("The length of upper must must match the number of parameters to be estimated.")
    message("upper:")
    print(upper)
    message("initial values:")
    names(parm) = names(upper)
    print(parm)
    stop("Please specify upper to have the same length as the number of parameters to be estimated.")
  }
  for (i in seq_along(lower)) {
    if (lower[i] > parm[i] | parm[i] > upper[i]) {
      message("The initial parameter values are not within the specified bounds.")
      combine = rbind(lower, parm, upper)
      rownames(combine) = c("lower", "initial", "upper")
      colnames(combine) = names(lower)
      print(combine)
      stop("Please specify lower and upper so that the initial starting values are between them or change your initial values.")
    }
  }
}

