



#' Defined integral of GeDS objects
#' @export
#' @description This function computes defined integrals of a fitted GeDS regression model.
#'
#'
#' @references De Boor, C. (2001). \emph{A Practical Guide to Splines (Revised Edition)}.
#' Springer, New York.
#'
#' @param object the \code{\link{GeDS-class}} object containing the  GeDS fit which should be integrated.
#' It should be the result of fitting a
#' univariate GeDS regression via \code{\link{NGeDS}} or \code{\link{GGeDS}}.
#' @param to numeric vector containing the upper limit(s) of integration.
#' @param from optional numeric vector containing the lower limit(s) of integration. It should be either
#' of size one or of the same size as the argument \code{to}. If left unspecified,
#' by default it is set to the left-most limit
#' of the interval embedding the observations of the independent variable.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the GeDS fit
#' to be integrated.
#' By default equal to \code{3L}.
#' Non-integer values will be passed to the function \code{\link{as.integer}}.
#' @details The function is based on the well known property (c.f. De Boor, 2001, Chapter X, formula (33))
#'  that the integral of a linear combination of appropriately normalized B-splines is equal to
#' the sum of its corresponding coefficients, noting that the GeDS regression is
#' in fact such a linear combination.
#'
#'
#' Since the function is based on this property, it is designed to work only on the predictor scale in
#' the GNM (GLM) framework.
#'
#' If the argument \code{from} is a single value, then it is taken as the lower limit of integration for all
#' the defined integrals required, whereas the upper limits of integration are the values contained
#' in the argument \code{to}.
#' If the arguments \code{from} and \code{to} are of similar size, the integrals (as many as the size)
#' are computed by sequentially taking the pairs of values in the \code{from} and \code{to} vectors as limits
#' of integration.
#'
#' @examples
#' # Generate a data sample for the response variable
#' # Y and the single covariate X
#' # see Dimitrova et al. (2017), section 4.1
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
#' # Fit GeDS regression using NGeDS
#' Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = .995, Xextr = c(-2,2))
#'
#' # Compute defined integrals (in TeX style) $\int_{1}^{-1} f(x)dx$
#' # and $\int_{1}^{1} f(x)dx$
#' # $f$ being the quadratic fit
#' Integrate(Gmod, to = c(-1,1), from = 1, n = 3)
#'
#' # Compute defined integrals (in TeX style) $\int_{1}^{-1} f(x)dx$
#' # and $\int_{-1}^{1} f(x)dx$
#' # $f$ being the quadratic fit
#' Integrate(Gmod, to = c(-1,1), from = c(1,-1), n = 3)
#'
#' \dontrun{
#' ## This gives an error
#' Integrate(Gmod, to = 1, from = c(1,-1), n = 3)
#' }
Integrate <- function(object, to, from, n = 3L){
  if(class(object)!= "GeDS") stop("incorrect object class")
  to <- as.numeric(to)
  l <- length(to)
  n <- as.numeric(n)
  if(missing(from)){
    from = min(knots(object, n=n, options="all"))
  } else {
    from <- as.numeric(from)
  }
  if(!length(from) %in% c(1,length(to))) stop("length of argument 'from' must be either 1 or length(to)")

  kn <- knots(object , options = "all", n = n)
  lastkn <- length(kn)
  thetas <- coef(object, n = n)

  newtheta <- numeric(length(thetas))
  knnew <- (kn[-(1:n)]-kn[-((lastkn-n+1):lastkn)])/n
  newtheta[1] <- thetas[1]*knnew[1]
  for(i in 2:(length(thetas))){
    newtheta[i] <- newtheta[i-1]+thetas[i]*knnew[i]
  }
  resTo <- sapply(to, gedsint, knts = kn, coefs = newtheta, n=n)
  resFrom <- sapply(from, gedsint, knts = kn, coefs = newtheta, n=n)

  res <- rowSums(cbind(resTo,-resFrom))
  return(res)
}



#' Derivative of GeDS objects
#' @export
#'
#' @description This function computes derivatives of a fitted GeDS regression model.
#'
#' @details The function is based on \code{\link[splines]{splineDesign}} and
#' it computes the exact derivative of the fitted GeDS regression.
#'
#' The function uses the property that the \eqn{m}-th derivative of a spline, \eqn{m= 1,2,...},
#'  expressed in terms of B-splines can be computed by differencing the corresponding B-spline coefficients
#'  (see e.g. De Boor, 2001, Chapter X, formula (15)).
#' Since the GeDS fit is a B-spline representation of the predictor,
#' it cannot work on the response scale in the GNM (GLM) framework.
#'
#' @param object the \code{\link{GeDS-Class}} object containing the GeDS fit which should be differentiated.
#' It should be the result of fitting a univariate GeDS regression via \code{\link{NGeDS}} or
#' \code{\link{GGeDS}}.
#' @param order integer value indicating the order of differentiation required (e.g. first, second or higher derivatives).
#' Note that \code{order} should be lower than \code{n} and that non-integer values will be passed
#' to the function \code{\link{as.integer}}.
#' @param x numeric vector containing values of the independent
#' variable at which the derivatives of order \code{order} should be computed.
#' @param n integer value (2, 3 or 4) specifying the order (\eqn{=} degree \eqn{+ 1}) of the GeDS fit
#' to be differentiated.
#' By default equal to \code{3L}.
#' @references De Boor, C. (2001). \emph{A Practical Guide to Splines (Revised Edition)}.
#' Springer, New York.
#'
#' @examples
#' # Generate a data sample for the response variable
#' # Y and the covariate X
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
#' # Fit GeDS regression with only a spline component in the predictor model
#' Gmod <- NGeDS(Y ~ f(X), beta = 0.6, phi = 0.995, Xextr = c(-2,2))
#'
#' # Compute the second derivative of the cubic GeDS fit
#' # at the points 0, -1 and 1
#' Derive(Gmod, x = c(0, -1, 1), order = 2, n = 4)
Derive <- function(object, order = 1L, x, n = 3L){
  if(class(object)!= "GeDS") stop("incorrect object class")
  if(!(object$Type %in% c("LM - Univ","GLM - Univ"))) stop("Implemented only for the univariate case")
  x <- as.numeric(x)
  l <- length(x)
  if(order >= n) stop("'order' must be less than 'n'")
  if(length(n)!=1) stop("'n' must have length 1")
  n <- as.integer(n)
  if(length(order)!=1) stop("'order' must have length 1")
  order <- as.integer(order)
  kn <- knots(object , options = "all", n = n)
  thetas <- coef(object, n = n)
  basis <- splines::splineDesign(knots=kn, x=x, ord = n, derivs = rep(order,l), outer.ok = TRUE)
  der <- as.numeric(basis%*%thetas)
  der
}

