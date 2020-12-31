#' @title Ensemble Forecast Verification for Large Data Sets
#'
#' @name easyVerification
#' 
#' @description Set of tools to simplify application of atomic forecast
#' verification metrics for (comparative) verification of ensemble forecasts
#' to large data sets. The forecast metrics are imported from the
#' 'SpecsVerification' package, and additional forecast metrics are provided
#' with this package. Alternatively, new user-defined forecast scores can be
#' implemented using the example scores provided and applied using the
#' functionality of this package.
#' 
#' @docType package
#' @import SpecsVerification pbapply stats utils
#' @importFrom Rcpp sourceCpp
#' @useDynLib easyVerification, .registration=TRUE
NULL
