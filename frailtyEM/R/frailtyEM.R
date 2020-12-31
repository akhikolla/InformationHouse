#' @details
#' This is an R package for fitting semiparametric shared frailty models with the EM algorithm.
#' You can check the "issues" section Github (link below)
#' For the gamma frailty model, the results are identical with those from the survival pacakage,
#' although frailtyEM provides a more readable output, including confidence intervals for the frailty variance.
#' Other supported distributions include the PVF, compound Poisson, inverse Gaussian, positive stable.
#' Univariate and multivariate data with left truncation are supported,
#' including recurrent events data in Andersen-Gill formulation.
#'
#' For background on the methods and basic usage of the package, see the paper in the references or the package vignette.
#' The main fitting function is \code{\link{emfrail}}. For prediction, see \code{\link{predict.emfrail}} and for
#' plotting, \code{\link{autoplot.emfrail}} (recommended, uses ggplot2) or \code{\link{plot.emfrail}}.
#'
#' @keywords internal
#' @author Theodor Balan \email{hello@@tbalan.com}
#' @references Balan TA, Putter H (2019) "frailtyEM: An R Package for Estimating Semiparametric Shared Frailty Models", \emph{Journal of Statistical Software} \strong{90}(7) 1-29. doi:10.18637/jss.v090.i07
"_PACKAGE"
