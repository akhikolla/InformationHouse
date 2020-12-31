#' Standard covariance models for geostatistical data.
#'
#' Creates a standard covariance model (\code{cmodStd})
#' object for geostatistical data.
#'
#' The general, isotropic form of the specified covariance function is
#' \code{psill} * \eqn{\rho}(\code{d}; \code{r}) +
#' (\code{evar} + \code{fvar}) * (\code{d == 0}), where
#' \eqn{\rho} is the correlation function of the parametric
#' models and \code{d} is the distance between the
#' relevant coordinates.
#'
#' For the \code{exponential} model, \eqn{\rho}(\code{d};
#' \code{r}) is exp(-\code{d}/\code{r}).
#'
#' For the \code{gaussian} model, \eqn{\rho}(\code{d};
#' \code{r}) is exp(-\code{d^2}/\code{r^2}).
#'
#' For the \code{matern} model, \eqn{\rho}(\code{d};
#' \code{r}) is
#' 2^(1-\code{par3})/\code{gamma}(\code{par3})*\code{sd}^\code{par3}*\code{besselK(sd,
#' nu = par3)}, where \code{sd = d/r}.
#'
#' For the \code{amatern} (alternative Matern) model,
#' \eqn{\rho}(\code{d}; \code{r}) is
#' \code{2^(1-par3)/gamma(par3)*sd^par3*besselK(sd, nu =
#' par3)}, where \code{sd = 2 * sqrt(par3) * d/r}.
#'
#' For the \code{spherical} model, \eqn{\rho}(\code{d};
#' \code{r}) is \code{1 - 1.5*sd + 0.5*(sd)^3} if \code{d <
#' r}, and 0 otherwise, with \code{sd = d/r}.
#'
#' For the \code{wendland1} model, \eqn{\rho}(\code{d};
#' \code{r}) is \code{(1 - sd)^4 * (4*sd + 1)} if \code{d <
#' r}, and 0 otherwise, with \code{sd = d/r}.
#'
#' For the \code{wendland2} model, \eqn{\rho}(\code{d};
#' \code{r}) is \code{(1 - sd)^6 * (35*sd^2 + 18*sd + 3))/3}
#' if \code{d < r}, and 0 otherwise, with \code{sd = d/r}.
#'
#' For the \code{wu1} model, \eqn{\rho}(\code{d}; \code{r})
#' is \code{(1 - sd)^3 * (1 + 3*sd + sd^2)} if \code{d < r},
#' and 0 otherwise, with \code{sd = d/r}.
#'
#' For the \code{wu2} model, \eqn{\rho}(\code{d}; \code{r})
#' is \code{(1 - sd)^4*(4 + 16*sd + 12*sd^2 + 3*sd^3))/4} if
#' \code{d < r}, and 0 otherwise, with \code{sd = d/r}.
#'
#' For the \code{wu3} model, \eqn{\rho}(\code{d}; \code{r})
#' is \code{(1 - sd)^6 * (1 + 6*sd + 41/3*sd^2 + 12*sd^3 +
#' 5*sd^4 + 5/6*sd^5)} if \code{d < r}, and 0 otherwise,
#' with \code{sd = d/r}.
#'
#' @param model A covariance model (e.g.,
#'   \code{"exponential"}). See Details for the complete
#'   list of choices.
#' @param psill The partial sill of the model.  Must be a
#'   positive number.
#' @param r The range parameter \code{r}.  Must be a
#'   positive number.
#' @param evar The variance of the errors.  Must be
#'   non-negative number.  The default is 0.
#' @param fvar The finescale variance (microscale error).
#'   Must be a non-negative number.  The default is 0.
#' @param par3 The value of the third parameter for 3
#'   parameter models.  Must be a positive number.  The
#'   default is 0.5.
#' @param longlat A logical value indicating whether great
#'   circle distance should be used. The default is
#'   \code{FALSE}.
#' @param angle The major axis of geometric anisotropy (the
#'   direction of strongest spatial dependence). Must be
#'   between [0, 180) if \code{radians = FALSE}, otherwise
#'   it must be between [0, \eqn{\pi}).
#' @param ratio The ratio of the minor axis range over the
#'   major axis range. The value must be between (0, 1].
#' @inheritParams angle2d
#' @return Returns a \code{cmodStd} object.
#'
#' @author Joshua French
#' @export
#' @references Waller, L. A., & Gotway, C. A. (2004).
#'   Applied Spatial Statistics for Public Health Data. John
#'   Wiley & Sons.
#' @examples
#' cmod_std(model = "exponential", psill = 1, r = 1)
cmod_std = function(model, psill, r, evar = 0,
                    fvar = 0, par3 = 0.5, longlat = FALSE,
                    angle = 0, ratio = 1,
                    radians = FALSE,
                    invert = TRUE) {
  arg_check_cmod_std(model, psill, r, evar, fvar, par3,
                     longlat, angle, ratio, radians, invert)
  structure(list(model = model, psill = psill, r = r,
                 evar = evar, fvar = fvar, par3 = par3,
                 longlat = longlat,
                 angle = angle,
                 ratio = ratio,
                 radians = radians,
                 invert = invert),
            class = "cmodStd")
}

#' Argument check cmod_std
#' @param model A standard covariance model type.
#' @param psill The partial sill of the model.  Must be a
#'   positive number.
#' @param r The range parameter r.  Must be a positive
#'   number.
#' @param evar The variance of the errors.  Must be a
#'   non-negative number.  The default is 0.
#' @param fvar The finescale variance (microscale error).
#'   Must be a non-negative number.  The default is 0.
#' @param par3 The value of the third parameter for 3
#'   parameter models.  Must be a positive number.  The
#'   default is 0.5.
#' @param longlat A logical value indicating whether great
#'   circle distance should be used. The default is
#'   \code{FALSE}.
#' @noRd
arg_check_cmod_std = function(model, psill, r, evar,
                              fvar, par3, longlat, angle,
                              ratio, radians, invert) {
  arg_check_cmod_std_model(model)
  arg_check_psill(psill)
  arg_check_r(r)
  arg_check_evar(evar)
  arg_check_fvar(fvar)
  arg_check_par3(par3)
  arg_check_longlat(longlat)
  arg_check_radians(radians)
  arg_check_angle(angle, radians)
  arg_check_ratio(ratio)
  arg_check_invert(invert)
}

