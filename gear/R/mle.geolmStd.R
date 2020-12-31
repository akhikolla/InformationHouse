#' Finds maximum likelihood estimates of model parameters
#' for a geostatistical model
#'
#' \code{mle} estimates the parameters of a geostatistical
#' linear model. The \code{mle} function will be deprecated
#' in the future. Please update your code to use the
#' estimate function.
#'
#' In the case of a \code{geolmStd} \code{object}, the
#' likelihood has been concentrated so that only the range
#' parameter \code{r} and a scale parameter \code{lambda =
#' nugget/psill} need to be estimated.
#'
#' If \code{object} is a \code{geolmStd}, then \code{lower}
#' is of length 2 if the covariance model of \code{cmod} is
#' not \code{matern} or \code{amatern}.  Otherwise, it
#' should be of length 3.  The first parameter is related to
#' the range parameter \code{r}, the second to the scale
#' parameter \code{lambda}, and the third to \code{par3}, if
#' applicable.  If \code{lower = NULL}, then the lower
#' bounds are 0.001, 0, and 0.1, respectively.  A similar
#' pattern holds for \code{upper}, with the default being
#' \code{3 * max(d)}, where \code{d} is the matrix of
#' distances between coordinates, \code{5}, and \code{2.5}.
#'
#' The \code{kkt} argument in the \code{control} list is set
#' to be \code{FALSE}.
#'
#'
#' @param object A geostatistical linear model object
#'   producted by the \code{geolm} function.
#' @param reml  A logical value indicating whether standard
#'   maximum likelihood estimation should be performed
#'   (\code{reml = FALSE}).  If \code{reml = TRUE}, then
#'   restricted maximum likelihood is performed.  Defaul is
#'   \code{FALSE}.
#' @param est A character vector indicator whether the error
#'   variance (\code{est="e"}) or finescale variance
#'   (\code{est = "f"}) should be estimated.  The other
#'   component of the nugget variance is held constant, and
#'   in the case of a \code{geolmStd} object, is set to 0.
#' @param lower A vector of 2 or 3 specifying the lowerbound
#'   of parameter values.  See Details.
#' @param upper lower A vector of 2 or 3 specifying the
#'   lowerbound of parameter values.  See Details.
#' @param method The optimization method.  Default is
#'   \code{"nlminb"}, with \code{"L-BFGS-B"} being another
#'   acceptable choice.  See \code{\link[optimx]{optimx}}
#'   for details.
#' @param itnmax An integer indicating the maximum number of
#'   iterations to allow for the optimization prodedure.
#' @param control A list of control parameters passed
#'   internally to \code{\link[optimx]{optimx}}.  See
#'   \code{\link[optimx]{optimx}} for details.
#' @param ... Currently unimplemented.
#'
#' @author Joshua French
#' @export
#' @examples
#' set.seed(10)
#' n = 100
#' @rdname mle
#' @export
mle = function(object, reml = FALSE, est = "e", lower = NULL,
               upper = NULL, method = "nlminb", itnmax = NULL,
               control = list(), ...) {
  warning("The mle function will be deprecated in the future. Please update your code to use the estimate function. The est argument is replaced by the noise_type argument in the estimate function.")
  estimate(object, reml = reml, noise_type = est, lower = lower, upper = upper,
           method = method, itnmax = itnmax, control = control, ...)
}
