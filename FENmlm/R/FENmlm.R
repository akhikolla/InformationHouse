#' Fixed Effects Nonlinear Maximum Likelihood Models
#'
#' Efficient estimation of multiple fixed-effects maximum likelihood models with, possibly, non-linear in parameters right hand sides. Standard-errors can easily be clustered. It also includes tools to seamlessly export (to Latex) the results of various estimations.
#'
#' @details
#' This package efficiently estimates maximum likelihood models with multiple fixed-effect (i.e. large factor variables).
#'
#' The core function is \code{\link[FENmlm]{femlm}} which estimates maximum likelihood models with, possibly, non-linear in parameters right hand sides. The ML families available are: poisson, negative binomial, logit and Gaussian.
#'
#'Several features are also included such as the possibility to easily compute different types of standard-errors (including multi-way clustering).
#'
#'It is possible to compare the results of severeal estimations by using the function \code{\link[FENmlm]{res2table}}, and to export them to Latex using \code{\link[FENmlm]{res2tex}}.
#'
#' @references
#' Berg\\'e, Laurent, 2018, "Efficient estimation of maximum likelihood models with multiple fixed-effects: the R package FENmlm." CREA Discussion Papers, 13 (\url{https://wwwen.uni.lu/content/download/110162/1299525/file/2018_13}).
#'
#'
"_PACKAGE"
