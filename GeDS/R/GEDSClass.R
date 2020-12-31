#' GeDS Class
#' @name GeDS-class
#' @rdname GeDS-class
#' @aliases GeDS-Class GeDS-class
#'
#' @description A fitted GeDS object returned by functions \code{\link{NGeDS}} or \code{\link{GGeDS}}
#' inheriting the methods from class \code{"GeDS"}.
#'  Methods for functions \code{coef}, \code{knots}, \code{print}, \code{predict}, \code{plot},
#' and  \code{lines} are available.
#'
#' @slot Type Character string indicating the type of the regression performed.
#' One of \code{"LM - Univ"}, \code{"LM - Biv"} or \code{"GLM - Univ"} corresponding to
#' the Normal univariate GeDS, the Normal bivariate GeDS performed by \code{\link{NGeDS}} and the
#' generalized (GNM-GLM) univariate GeDS performed by \code{\link{GGeDS}}.
#' @slot Linear.Knots  Vector containing the
#'  locations of the knots of the second order GeDS spline fit generated at stage A.
#' @slot Quadratic.Knots  Vector containing the locations
#' of the knots of the third order GeDS spline fitted in stage B.
#' @slot Cubic.knots  Vector containing the locations of
#' the knots of the fourth order GeDS spline fitted in stage B.
#' @slot Dev.Linear Deviance of the second order GeD spline fit of stage A.
#' @slot Dev.Quadratic Deviance of the third order GeD spline fit of stage B.
#' @slot Dev.Cubic Deviance of the fourth order GeD spline fit of stage B.
#' @slot Linear List containing the results from running a \code{\link{SplineReg}}
#' function used to fit the second order spline of stage A.
#' @slot Quadratic List containing the results from running \code{\link{SplineReg}}
#' function used to fit the third order spline in stage B.
#' @slot Cubic List containing the results from a \code{\link{SplineReg}}
#' function used to fit the fourth order spline in stage B.
#' @slot Stored Matrix containing the knot locations estimated at each step of stage A.
#' @slot Args List containing the input arguments passed on the \code{\link{Fitters}} functions.
#' @slot Call \code{call} to the \code{\link{Fitters}} functions.
#' @slot Nintknots The final number of internal knots of the second order GeD spline fit of stage A.
#' @slot iters Number of iterations performed in stage A  of the GeDS fitting procedure.
#' @slot Guesses Initial values for the coefficients used at each
#' iteration of stage A in order to estimate the spline coefficients.
#' Since the initial values are used only in the IRLS procedure,
#' this slot is empty if the object is not created by \code{\link{GGeDS}} or \code{\link{GenUnivariateFitter}}
#'  functions.
#' @slot Coefficients Matrix containing the fitted coefficients of the GeD spline regression  component and the
#' parametric component at each iteration of stage A.
#' @slot deviance Vector containing the deviances of the second order spline fits computed at each IRLS
#' iteration in stage A.  Since the IRLS procedure is used only in \code{\link{GGeDS}} or \code{\link{GenUnivariateFitter}},
#' this slot is empty if the object is not created by one of these functions.
#' @slot iterIrls Vector containing the numbers of IRLS iterations for all iterations of stage A cumulatively.
#' Since the IRLS procedure is used only in \code{\link{GGeDS}} or \code{\link{GenUnivariateFitter}},
#' this slot is empty if the object is not created by one of these functions.
#' @slot stopinfo List of values providing information related to
#' the stopping rule of stage A of GeDS. The sub-slots of \code{stopinfo} are \code{phis},  \code{phis_star},
#' \code{oldintc} and \code{oldslp}. The sub-slot \code{phis} is a vector containing the values
#' of the ratios of deviances (or the difference of deviances if the \code{LR} stopping
#' rule was chosen). The sub-slots \code{phis_star}, \code{oldintc} and \code{oldslp} are non-empty slots
#' if the \code{SR} stopping rule was chosen. They contain respectively \eqn{\hat{\phi}_{\kappa}}, \eqn{\hat{\gamma}_0} and
#' \eqn{\hat{\gamma}_1} computed at each iteration of stage A, see Dimitrova et al. (2017).
#' @slot Formula The model \code{\link[=formula.GeDS]{formula}}.
#' @slot extcall \code{call} to the \code{\link{NGeDS}} or \code{\link{GGeDS}} functions.
#' @slot terms \code{terms} object containing information on the model frame.
#'
#' @references
#' Dimitrova, D.S., Kaishev, V.K., Lattuada A. and Verrall, R.J. (2017).
#' Geometrically designed, variable knot splines in Generalized (Non-)Linear Models.
#' Available at \href{http://openaccess.city.ac.uk/18460/}{openaccess.city.ac.uk}
#'
#'
#'
setClass(
  "GeDS",
  representation(
    Type = "character", Linear.Knots = "numeric", Quadratic.Knots = "numeric",
    Cubic.Knots = "numeric", RMS.Linear  = "numeric", RMS.Quadratic = "numeric",
    RMS.Cubic  = "numeric", Knots  = "numeric", RSS = "numeric",
    Linear = "list", Quadratic = "list", Cubic = "list", Stored = "matrix",
    Args = "list", Call = "call", Nintknots = "numeric", iters = "numeric", Guesses = "matrix",
    Coefficients = "matrix", deviance = "numeric", iter = "numeric", stopinfo = "list", Formula = "formula", extcall = "call"
    )
  )


