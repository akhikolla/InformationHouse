#' @include forecastSNSTS-package.R
NULL

################################################################################
#' Forecasting of Stationary and Non-Stationary Time Series
#'
#' Methods to compute linear \eqn{h}-step ahead prediction coefficients based
#' on localised and iterated Yule-Walker estimates and empirical mean squared
#' and absolute prediction errors for the resulting predictors. Also, functions
#' to compute autocovariances for AR(p) processes, to simulate tvARMA(p,q) time
#' series, and to verify an assumption from Kley et al. (2019).
#'
#' @details
#'  \tabular{ll}{
#'    \cr Package: \tab forecastSNSTS
#'    \cr Type:    \tab Package
#'    \cr Version: \tab 1.3-0
#'    \cr Date:    \tab 2019-09-02
#'    \cr License: \tab GPL (>= 2)
#'  }
#'
#' @section Contents:
#' The core functionality of this R package is accessable via the function
#' \code{\link{predCoef}}, which is used to compute the linear prediction
#' coefficients, and the functions \code{\link{MSPE}} and \code{\link{MAPE}},
#' which are used to compute the empirical mean squared or absolute prediction
#' errors. Further, the function \code{\link{f}} can be used to verify
#' condition (10) of Theorem 3.1 in Kley et al. (2019) for any given tvAR(p) model.
#' The function \code{\link{tvARMA}} can be used to simulate time-varying
#' ARMA(p,q) time series.
#' The function \code{\link{acfARp}} computes the autocovariances of a AR(p)
#' process from the coefficients and innovations standard deviation.
#' 
#'  
#'
#' @name forecastSNSTS-package
#' @aliases forecastSNSTS
#' @docType package
#' @author Tobias Kley
#'
#' @useDynLib forecastSNSTS
#' @importFrom Rcpp sourceCpp

#'
#' @references
#' Kley, T., Preuss, P. & Fryzlewicz, P. (2019).
#' Predictive, finite-sample model choice for time series under stationarity
#' and non-stationarity. Electronic Journal of Statistics, forthcoming.
#' [cf. \url{https://arxiv.org/abs/1611.04460}]
#'
NULL

# Taken from quantreg-package and adapted.
".onAttach" <- function(lib, pkg) {
  if(interactive() || getOption("verbose"))
    packageStartupMessage("Package forecastSNSTS loaded.\n     To cite, see citation(\"forecastSNSTS\").\n     For demos, see demo(package = \"forecastSNSTS\").")
}
