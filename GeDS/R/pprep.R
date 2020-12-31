#' Piecewise Polynomial Spline Representation
#'
#' @description The function converts a GeDS fit which has a  B-spline representation
#' to a piecewise polynomial form.
#'
#' @param object  the \code{\link{GeDS-class}} where the GeDS fit to be converted is found.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the GeDS fit which should be converted to a piecewise polynomial form.
#' By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @details
#' This function converts a selected GeDS fit from a \code{\link{GeDS-class}} object
#' represented in terms of B-splines into an object where the fit is represented in terms
#' of piecewise polynomials.
#'
#' The function  wraps \code{\link[splines]{polySpline}} in order to let it
#' accept \code{\link{GeDS-class}} objects as input.
#' Hence the function provides a useful link between the package \pkg{GeDS} and the package \pkg{splines}, allowing
#' the user to take advantage of the functions provided in the \pkg{splines} package.
#'
#'
#'
#' @note Let us note that the first \eqn{k+1} rows of the matrix contain the \code{n} coefficients
#'  of the \eqn{k+1} consecutive pieces of the piecewise polynomial representation. The last \eqn{(k+2)}-th
#'  row is extraneous and it appears as a result of the use of the function \code{\link[splines]{polySpline}}.
#'
#' @return An object that inherits from classes  \code{"spline"} and \code{"polySpline"}. It is a list
#' whose arguments are:
#'  \item{knots}{ a vector of size  \eqn{k + 2}
#'  containing the  complete set of  knots (internal knots plus the limits of the interval)
#'  of the GeDS fit.}
#'  \item{coefficients}{ a \eqn{(k + 2) \times n} matrix containing the coefficients of the  polynomials
#'  in the required piecewise polynomial representation. }
#'
#' @examples
#' # Generate a data sample for the response variable
#' # Y and the single covariate X
#' set.seed(123)
#' N <- 500
#' f_1 <- function(x) (10*x/(1+100*x^2))*4+4
#' X <- sort(runif(N, min = -2, max = 2))
#' # Specify a model for the mean of Y to include only
#' # a component non-linear in X, defined by the function f_1
#' means <- f_1(X)
#' # Add (Normal) noise to the mean of Y
#' Y <- rnorm(N, means, sd = 0.1)
#'
#' # Fit a Normal GeDS regression using NGeDS
#' Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995, Xextr = c(-2,2))
#'
#' # construct the PP representation of the cubic GeDS fit
#' # and apply some functions of the package splines
#' Polymod <- PPolyRep(Gmod, 4)
#' require(splines)
#' class(Polymod)
#' splineKnots(Polymod)
#' knots(Gmod, n = 4)
#' plot(Polymod)
#'
#'
#' # Generate a plot showing the PP representation
#' # based on the same example
#' knt <- splineKnots(Polymod)
#' coeffs <- coef(Polymod)
#' plot(Gmod, n = 4, legend = FALSE, main = "Cubic Curves")
#' cols <- sample(heat.colors(length(knt)), length(knt))
#' for(i in 1:(length(knt))){
#'   curve(coeffs[i,1] + coeffs[i,2]*(x - knt[i])+
#'           coeffs[i,3]*(x - knt[i])^2+
#'         coeffs[i,4]*(x - knt[i])^3,
#'         add = TRUE, col = cols[i])
#'   abline(v = knt[i])
#' }
#'
#' @export
PPolyRep <- function(object, n = 3){
  if(class(object)!="GeDS") stop("This function works only with GeDS class objects")
  n <- as.integer(n)
  kn <- knots(object, n = n, options="all")
  cf <- coef(object, n = n)
  newlist <- list(knots = kn, coefficients = cf, order = n)
  class(newlist) <- c("nbSpline", "bSpline", "spline")
  xname <- attr(object$terms,"specials")$f-1
  xname <- attr(object$terms,"term.labels")[xname]
  xname <- substr(xname,3,(nchar(xname)-1))
  yname <- rownames(attr(object$terms,"factors"))[1]
  fortmp <- paste0(yname," ~ ", xname)
  attr(newlist,"formula") <- as.formula(fortmp)
  out <- polySpline(newlist)
  return(out)
}
