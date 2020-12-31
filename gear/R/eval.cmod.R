#' Evaluate covariance or semivariance model.
#'
#' \code{eval.cmod} evaluates the covariance 
#' of a model based on the provided arguments.  This
#' function will be deprecated in the future. Please use
#' the \code{\link[gear]{evaluate}} function.
#'
#' @param mod A covariance or semivariance model.
#' @param d An \eqn{n \times m} matrix of distances.
#' @param coords Not used. 
#'
#' @return Returns the evaluated model with necessary
#'   components needed for \code{\link[gear]{estimate}} and
#'   \code{\link{predict}}.
#' @author Joshua French
#' @examples
#' n = 10
#' coords = matrix(runif(2*n), nrow = n, ncol = 2)
#' d = as.matrix(dist(coords))
#' cmod = cmod_std(model = "exponential", psill = 1, r = 1)
#' eval.cmod(cmod, d)
#' @export
eval.cmod = function(mod, d, coords = NULL) {
  warning("This function will be deprecated in the future. Please update your code to use the evaluate function.")
  evaluate(mod, d = d)
}

