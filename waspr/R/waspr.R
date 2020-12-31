#' waspr: an R package for computing Wasserstein barycenters of subset
#' posteriors
#'
#' This package contains functions to compute Wasserstein barycenters of subset
#' posteriors using the swapping algorithm developed by Puccetti, Rüschendorf
#' and Vanduffel (2020). The Wasserstein barycenter is a geometric approach for
#' combining subset posteriors. It allows for parallel and distributed
#' computation of the posterior in case of complex models and/or big datasets,
#' thereby increasing computational speed tremendously.
#'
#' @section Functions:
#'
#'   The main function of the package is:
#'
#'   \code{\link{wasp}}, which runs the swapping algorithm developed by
#'   Puccetti, Rüschendorf and Vanduffel (2020), combines the output from the
#'   swapping algorithm and computes the Wasserstein barycenter. It returns an
#'   S3 object of type \code{wasp}.
#'
#' @source Puccetti, G., Rüschendorf, L. & Vanduffel, S. (2020). On the
#'   computation of Wasserstein barycenters, Journal of Multivariate Analysis,
#'   176.
#'
#' @useDynLib waspr, .registration = TRUE
#'
#' @importFrom Rcpp evalCpp
#'
#' @docType package
#' @name waspr

NULL
