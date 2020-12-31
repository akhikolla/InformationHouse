#' Standard covariance models for geostatistical data.
#'
#' Creates a standard covariance model (\code{cmodStd})
#' object for geostatistical data. This function will be
#' deprecated in the future. Please update your code to use
#' the \code{\link[gear]{cmod_std}} function, which also
#' allows the user to specify geometric anisotropy.
#'
#' The general form of the specified covariance function is
#' \code{psill} * \eqn{\rho}(\code{d}; \code{r}) +
#' (\code{evar} + \code{fvar})*(\code{d==0}), where
#' \eqn{\rho} is the covariance function of the parametric
#' models.
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
#' @inheritParams cmod_std
#' @return Returns a \code{cmodStd} object.
#' @author Joshua French
#' @export
#' @references Waller, L. A., & Gotway, C. A. (2004).
#'   Applied Spatial Statistics for Public Health Data. John
#'   Wiley & Sons.
#' @examples
#' cmod.std(model = "exponential", psill = 1, r = 1)
cmod.std = function(model, psill, r, evar = 0, fvar = 0,
                    par3 = 0.5) {
  warning("This function will be deprecated in the future. Please update your code to use the cmod_std function, which also allows the user to specify geometric anisotropy.")
  cmod_std(model = model, psill = psill, r = r, evar = evar,
           fvar = fvar, par3 = par3)
}


