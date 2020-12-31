#' @useDynLib ffstream, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats qnorm
#' @importFrom stats rnorm
#' @importFrom methods new
#' @importFrom graphics abline
#' @importFrom graphics plot
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom grDevices dev.new
#' @importFrom methods new
#' @importFrom Rcpp cpp_object_initializer
#library(methods)
#library(Rcpp)
Rcpp::loadModule("fffmodule", TRUE)
Rcpp::loadModule("affmodule", TRUE)
Rcpp::loadModule("fffcdmodule", TRUE)
Rcpp::loadModule("affcdmodule", TRUE)
Rcpp::loadModule("cusumcdmodule", TRUE)
Rcpp::loadModule("ewmacdmodule", TRUE)
