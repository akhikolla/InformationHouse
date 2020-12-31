#' Tools for climate- and ecological-niche factor analysis
#'
#' \code{CENFA} provides tools for performing ecological-niche factor analysis
#' (ENFA) and climate-niche factor analysis (CNFA).
#'
#' @details
#' This package was created with three goals in mind:
#'
#' - To update the ENFA method for use with large datasets and modern data formats.
#'
#' - To expand the application of ENFA in the context of climate change in order
#'   to quantify different aspects of species vulnerability to climate change,
#'   and to facilitate quantitative comparisons of vulnerability between species.
#'
#' - To correct a minor error in the ENFA method itself, that has persisted in the
#'   literature since Hirzel et al. first introduced ENFA in 2002.
#'
#'
#' \code{CENFA} takes advantage of the \code{raster} and \code{sp} packages,
#' allowing the user to conduct analyses directly with raster, shapefile, and
#' point data, and to handle large datasets efficiently via partial data loading
#' and parallelization.
#'
#' In addition, \code{CENFA} also contains a few functions that speed up some
#' basic `raster` functions considerably by parallelizing on a layer-by-layer
#' basis rather than a cell-by-cell basis.
#'
#' @author D. Scott Rinnan
#'
#' @references
#' Basille, Mathieu, et al. Assessing habitat selection using multivariate
#' statistics: Some refinements of the ecological-niche factor analysis. Ecological
#' Modelling 211.1 (2008): 233-240.
#'
#' Hirzel, Alexandre H., et al. Ecological-niche factor analysis: how to compute
#' habitat-suitability maps without absence data?. Ecology 83.7 (2002): 2027-2036.
#'
#' @seealso \code{\link[sp]{sp}}, \code{\link[raster]{raster-package}}
#'
#' @docType package
#' @name CENFA-package
#' @useDynLib CENFA
#' @importFrom Rcpp sourceCpp
# @importClassesFrom raster RasterBrick RasterStack RasterLayer
#' @import raster
#' @importClassesFrom sp SpatialPolygons SpatialPolygonsDataFrame SpatialPoints SpatialPointsDataFrame
# @importMethodsFrom raster raster
#' @importFrom magrittr %>%
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\nSupport women in science: 500womenscientists.org")
}

.onUnload <- function (libpath) {
  library.dynam.unload("CENFA", libpath)
  closeAllConnections()
}
