#' RZigZag
#' 
#' Implements various piecewise deterministic Monte Carlo methods, including the Zig-Zag Sampler (Bierkens, Fearnhead, Roberts, 2019, \url{https://arxiv.org/abs/1607.03188}) and the Bouncy Particle Sampler (BPS, Bouchard-Côté et al., 2017, \url{https://arxiv.org/abs/1510.02451}).
#' Typical usage consists of first creating a "skeleton" consisting of "events", which can be used directly for plotting trajectories.
#' The skeleton may be post-processed to extract information, such as as moment and covariance estimates, discrete samples at fixed time intervals along the trajectory, effective sample size and asymptotic variance.
#' @details 
#' This package currently consists of the following functions for generating skeletons:
#' \code{\link{ZigZagLogistic}} for logistic regression, \code{\link{ZigZagGaussian}} for multivariate Gaussian, \code{\link{ZigZagIIDGaussian}} for a IID Gaussian target using Zig-Zag, \code{\link{ZigZagStudentT}} for spherically symmetric or factorized Student-t distribution, \code{\link{BPSGaussian}} for multivariate Gaussian using BPS, \code{\link{BPSIIDGaussian}} for a IID Gaussian target using BPS, \code{\link{BPSStudentT}} for BPS applied to a spherically symmetric or factorized Student-t distribution.
#' Furthermore the package contains the following functions for post-processing:
#' \code{\link{EstimateESS}} (to estimate asymptotic variance and effective sample size for individual coordinates), \code{\link{EstimateMoment}}, \code{\link{EstimateCovarianceMatrix}} and \code{\link{DiscreteSamples}}.
#' 
#' @docType package
#' @author Joris Bierkens
#' @author With thanks to Matt Moores, \url{https://mattstats.wordpress.com/}, for his help in getting from C++ code to a CRAN-ready Rcpp based package.
#' @import Rcpp 
#' @importFrom Rcpp evalCpp
#' @useDynLib RZigZag
#' @name RZigZag
NULL  
