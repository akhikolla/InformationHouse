#' PS-Integrated Methods for Incorporating RWE in Clinical Studies
#'
#' @docType   package
#' @name      psrwe-package
#' @aliases   psrwe
#' @useDynLib psrwe, .registration = TRUE
#'
#' @import Rcpp
#' @import methods
#' @import ggplot2
#'
#' @importFrom stats approxfun as.formula binomial cov density ecdf glm
#'     integrate optim predict quantile sd var
#'
#' @importFrom rstan sampling extract stanc rstan_options traceplot stan_rhat
#' @importFrom randomForest randomForest
#'
#' @importFrom grDevices colors
#' @importFrom graphics  axis box legend lines par plot points text
#'
#' @importFrom parallel detectCores
#' @importFrom cowplot plot_grid
#' @importFrom dplyr %>% group_by_ group_by summarize mutate count_ mutate_if
#'     rename_ filter select arrange
#'
#'
#' @description
#'
#' This package provide R functions for conducting clinical studies with
#' real-world evidence incorporated in the study design and analysis.
#'
#' @section PS-integrated power prior:
#'
#' We extend the Bayesian power prior approach for a single-arm study (the
#' current study) to leverage external real-world data. We use propensity score
#' methodology to pre-select a subset of real-world data containing patients
#' that are similar to those in the current study in terms of covariates, and to
#' stratify the selected patients together with those in the current study into
#' more homogeneous strata. The power prior approach is then applied in each
#' stratum to obtain stratum-specific posterior distributions, which are
#' combined to complete the Bayesian inference for the parameters of interest.
#'
#' @section PS-integrated composite likelihood:
#'
#' A propensity score-integrated composite likelihood (PSCL) approach is
#' developed for cases in which the control arm of a two-arm randomized
#' controlled trial (RCT) (treated vs control) is augmented with patients from
#' real-world data (RWD) containing both clinical outcomes and covariates at the
#' patient-level. The PSCL approach first estimates the propensity score for
#' every patient as the probability of the patient being in the RCT rather than
#' the RWD, and then stratifies all patients into strata based on the estimated
#' propensity scores. Within each propensity score stratum, a composite
#' likelihood function is specified and utilized to down-weight the information
#' contributed by the RWD source. Estimates of the stratum-specific parameters
#' are obtained by maximizing the composite likelihood function. These
#' stratum-specific estimates are then combined to obtain an overall
#' population-level estimate of the parameter of interest.
#'
#'
#' @references
#'
#' Chen, W.C., Wang, C., Li, H., Lu, N., Tiwari, R., Xu, Y. and Yue, L.Q., 2020.
#' Propensity score-integrated composite likelihood approach for augmenting the
#' control arm of a randomized controlled trial by incorporating real-world
#' data. Journal of Biopharmaceutical Statistics, 30(3), pp.508-520.
#'
#' Wang, C., Lu, N., Chen, W. C., Li, H., Tiwari, R., Xu, Y., & Yue, L. Q.
#' (2020). Propensity score-integrated composite likelihood approach for
#' incorporating real-world evidence in single-arm clinical studies. Journal of
#' biopharmaceutical statistics, 30(3), 495-507.
#'
#' Wang, C., Li, H., Chen, W. C., Lu, N., Tiwari, R., Xu, Y., & Yue, L. Q.
#' (2019). Propensity score-integrated power prior approach for incorporating
#' real-world evidence in single-arm clinical studies. Journal of
#' biopharmaceutical statistics, 29(5), 731-748.
#'
NULL


#' Example dataset
#'
#' Example dataset of a single arm study
#'
#' @name ex_dta
#'
#' @format A dataframe with the following variables:
#' \itemize{
#'   \item{Group}{current, rwd}
#'   \item{Y}{Binary outcome}
#'   \item{V1-V7}{Covariates}
#' }
NULL


#' Example dataset
#'
#' Example dataset of a randomized study
#'
#' @name ex_dta_rct
#'
#' @format A dataframe with the following variables:
#' \itemize{
#'   \item{Group}{0,1}
#'   \item{Arm}{0, 1}
#'   \item{Y}{Continuous outcome}
#'   \item{V1-V7}{Covariates}
#' }
NULL
