#' Goodness of fit checks for binomial N-mixture models
#'
#' The package contains methods to compute overdispersion metrics, randomized quantile residuals,
#' and graphical diagnostics of model fit for binomial N-mixture models fitted using the \link[unmarked]{unmarked} package.
#' Details about the checks are given in Knape et al. (2018) and at \url{https://www.biorxiv.org/content/early/2017/09/27/194340}.
#'
#'
#' @docType package
#' @name nmixgof
#' @references Knape et al. 2018. Sensitivity of binomial N-mixture models to overdispersion: 
#' the importance of assessing model fit. Methods in Ecology and Evolution, in press.
#' @import unmarked
#' @useDynLib nmixgof
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics abline
#' @importFrom stats dnbinom dpois pbinom plogis pnbinom ppois qnorm qqnorm runif
NULL


