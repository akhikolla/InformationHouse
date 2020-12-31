#' landscapeR: A landscape simulator for R.
#'
#' This package is aimed at simulating categorical landscapes on actual geographical realms,
#' starting from either empty landscapes, or landscapes provided by the user (e.g. land use maps).
#' landscapeR allows to tweak or create landscapes while retaining a high degree of control on its features,
#' without the hassle of specifying each location attribute. In this it differs from other tools
#' which generate null or neutral landscape in a theoretical space. The basic algorithm currently
#' implemented uses a simple agent style/cellular automata growth model, with no rules
#' (apart from areas of exclusion). Outputs are raster dataset exportable to any common GIS format.
#'
#' @useDynLib landscapeR
#' @importFrom Rcpp sourceCpp
#'
#' @section landscapeR functions:
#' \itemize{
#' \item{\code{\link{makePatch}} creates a single patch in the landscape.}
#' \item{\code{\link{makeClass}} creates a group of patches belonging to the same class.}
#' \item{\code{\link{expandClass}} expands an existing class of patches.}
#' \item{\code{\link{makeLine}} creates a linear patch.}
#' \item{\code{\link{rmSingle}} removes single tones from patches and background.}
#' }
#' @details Check out the vignette illustrating the use of landscapeR.\cr
#' Also: \url{https://github.com/dariomasante/landscapeR}
#' @docType package
#' @name landscapeR-package
#' @aliases landscapeR
#' @author Dario Masante
NULL
