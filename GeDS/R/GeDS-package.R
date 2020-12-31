# #' Implemetation of the Geometrically Designed Spline (GeDS) Regression.
# #'  first on B-splines that automatically selects number and position of the
# #' knots and the order of the splines to be used.
# #' The package includes both
# #' the original version of the algorithm, developed in order to work with Gaussian error
# #' and its generalization to the broader GLM (GNM) context and to the bivariate framework. The
# #' latter, in version 0.1, is available only for the Normal case.
# #' The package includes two main functions: \code{\link{NGeDS}} and
# #' \code{\link{GGeDS}}. The first provides a user friendly implementation of the version
# #' of the algorithm that was orignally studied in Kaishev et al. (2016) in orer to work
# #' with Gaussian error, with the possibility to set two independent variables to be modelled via
# #' tensor product splines. The latter is a generalization that allows the distribution of the
# #' response variable to be one of the Exponential Family (see Dimitrova et al., 2017).


#' GeDS
#'
#' Geometrically Designed Splines (GeDS) regression is a non-parametric geometrically motivated method
#' for fitting variable knots spline predictor models which are generalized (non)-linear, (i.e. GNM (GLM))
#' models. The GeDS regression is fitted  based on a sample of \eqn{N} observations of a response variable \eqn{y}, dependent
#' on a set of (currently up to two) covariates, assuming \eqn{y} has a distribution from the exponential family.
#'
#'
#' @details
#' The GeDS estimation method is based on: first constructing a piecewise linear
#' fit (spline fit of order 2) at stage A  which captures the shape of the data and;
#' second approximating this fit with shape preserving
#' (variation diminishing) spline fits of higher orders 3, 4,\eqn{\dots} (i.e. degrees 2, 3,\eqn{\dots}) at stage B.
#' As a result of this,  GeDS estimates the number and location of the knots
#' and the order of the spline fit in a fast and efficient way.
#'
#' The GeDS method was originally developed by Kaishev et al. (2016) assuming the response \eqn{y} is normally
#' distributed and a corresponding \emph{Mathematica} code was provided.
#'
#' The GeDS method was recently extended by Dimitrova et al. (2017)  to cover any distribution from the exponential
#' family. The \pkg{GeDS} \R package presented here includes an enhanced \R implementation of the original
#' Normal GeDS \emph{Mathematica} code due to Kaishev et al. (2016), implemented as the \code{\link{NGeDS}} function and a generalization of it in the function
#' \code{\link{GGeDS}} which covers the case of any distribution from the exponential family.
#'
#' The \pkg{GeDS} package allows also to fit two dimensional response surfaces
#' currently implemented only in the Normal case via the function \code{\link{NGeDS}}. It also allows
#'  to construct multivariate (predictor) models with a GeD spline
#' component and a parametric component (see the functions \code{\link{f}}, \code{\link[=formula.GeDS]{formula}},
#' \code{\link{NGeDS}} and \code{\link{GGeDS}} for details).
#'
#' The outputs of both \code{\link{NGeDS}} and \code{\link{GGeDS}} functions are \code{\link{GeDS-class}} objects.
#' As described in Kaishev et al. (2016) and Dimitrova et al. (2017)
#' the final GeDS fit is the one whose order is chosen according to a strategy described
#' in stage B of the algorithm. However, \code{\link{GeDS-class}} objects contain second, third and fourth
#' order spline fits and the user has the possibility to choose among them.
#'
#' This package also includes some datasets where GeDS regression proves to be very efficient
#' and some user friendly functions that are designed to easily extract required
#' information.  Several methods are also provided to handle GeDS output results (see \code{\link{GeDS-class}}).
#'
#' %This package comes distributed with a \emph{vignette} that may help unexperienced
#' %users to install \R and this package and to explore its main features.
#'
#' Throughout this document, we use the terms GeDS predictor model, GeDS regression and GeDS fit
#' interchangeably.
#'
#' Please report any issue arising or bug in the code to andrea.lattuada@unicatt.it.
#'
#' \tabular{rl}{
#' Package: \tab GeDS\cr
#' Version: \tab 0.1.3 \cr
#' Date: \tab 2017-12-19\cr
#' License: \tab GPL-3 \cr
#' }
#'
#' @references
#' Kaishev, V.K., Dimitrova, D.S., Haberman, S., & Verrall, R.J. (2016).
#' Geometrically designed, variable knot regression splines.
#' \emph{Computational Statistics}, \strong{31}, 1079--1105. \cr
#' DOI: \href{https://doi.org/10.1007/s00180-015-0621-7}{doi.org/10.1007/s00180-015-0621-7}
#'
#' Dimitrova, D.S., Kaishev, V.K., Lattuada A. and Verrall, R.J. (2017).
#' Geometrically designed, variable knot splines in Generalized (Non-)Linear Models.
#' Available at \href{http://openaccess.city.ac.uk/18460/}{openaccess.city.ac.uk}
#'
#' @docType package
#' @name GeDS-package
#'
#' @author  Dimitrina S. Dimitrova <D.Dimitrova@city.ac.uk>,
#' Vladimir K. Kaishev <V.Kaishev@city.ac.uk>,
#' Andrea Lattuada <Andrea.Lattuada@unicatt.it> and
#' Richard J. Verrall <R.J.Verrall@city.ac.uk>
#' @keywords package
#'
#' @aliases GeDS
#'
#'
NULL

