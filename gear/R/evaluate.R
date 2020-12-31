#' Evaluate spatial dependence model
#'
#' \code{evaluate} evaluates the spatial dependence model
#' based on the provided arguments.
#'
#' If \code{mod} is of class \code{cmodStd} (from the
#' \code{cmod_std} function), then the function returns an
#' \eqn{n \times m} matrix with the evaluated standard
#' covariance function.
#'
#' @param mod A covariance or semivariogram model.
#' @param d An \eqn{n \times m} matrix of distances.
#' If \code{mod$ratio != 1}, i.e., if geometric anisotropy
#' has been specified, then \code{d} must be produced by the
#' \code{\link[gear]{ganiso_d}} function.
#' @param e A single logical value indicating whether the
#'   error variance should be added to the returned
#'   covariance matrix.  Default is \code{TRUE}.
#' @param f A single logical value indicating whether the
#'   finescale/microscale variance should be added to the
#'   returned covariance matrix.  Default is \code{TRUE}.
#'
#' @return Returns the evaluated model with necessary
#'   components needed for \code{\link[gear]{estimate}} and
#'   \code{\link{predict}}.
#' @author Joshua French
#' @export
#' @examples
#' n = 10
#' coords = matrix(runif(2*n), nrow = n, ncol = 2)
#' d = as.matrix(dist(coords))
#' cmod = cmod_std(model = "exponential", psill = 1, r = 1)
#' evaluate(cmod, d)
evaluate = function(mod, d, e = TRUE, f = TRUE) {
  UseMethod("evaluate")
}
