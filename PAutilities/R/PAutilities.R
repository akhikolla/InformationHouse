#' PAutilities: Streamline Physical Activity Research
#'
#' A collection of utilities that are useful for a broad range of tasks that are
#' common in physical activity research. The main features (with associated
#' functions in parentheses) are:
#'
#' * Bland-Altman plots (\code{\link{ba_plot}})
#' * Bout analysis for moderate-to-vigorous physical activity (\code{\link{bout_mvpa}})
#' * Formatted descriptive statistics \code{\link{descriptives}}
#' * Demographic calculations (\code{\link{get_age}} and \code{\link{get_BMI_percentile}})
#' * Metabolic calculations (\code{\link{get_bmr}}, \code{\link{weir_equation}}, and \code{\link{get_kcal_vo2_conversion}})
#' * Analysis of bout detection algorithm performance (\code{\link{get_transition_info}} and associated methods, e.g. \code{summary} and \code{plot})
#'
#' @docType package
#' @name PAutilities
NULL

#' @import ggplot2 dplyr methods
NULL

#' @importFrom magrittr %>% %T>%
NULL

#' @importFrom stats sd median quantile
NULL

#' @importFrom rlang :=
NULL

#' @importFrom graphics par plot title
NULL

#' @useDynLib PAutilities, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

utils::globalVariables(c(
  "x_label", "CI_low", "CI_high",
  "CI_sig", "low", "high", "vartype",
  "."
))
