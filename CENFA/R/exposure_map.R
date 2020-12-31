#' Create an exposure map
#'
#' Creates a map of exposure to climate change in a species' habitat from a
#' \code{departure} object.
#'
#' @param dep Object of class \code{departure}
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Number of cores to use for calculation (optional)
#' @param filename character. Output filename (optional)
#' @param ... Additional arguments for file writing as for \code{\link[raster]{writeRaster}}
#'
#' @details
#' The values of the exposure raster are calculated by projecting onto the
#' departure factor \strong{d}, given by the formula
#'
#'   \eqn{\epsilon} = \bold{Fd}.
#'
#' @examples
#' dep <- departure(x = climdat.hist, y = climdat.fut, s.dat = ABPR)
#' exp.map <- exposure_map(dep)
#'
#' @return A RasterLayer of exposure values
#'
#' @seealso \code{\link{departure}}, \code{\link{sensitivity_map}}, \code{\link{vulnerability_map}}
#'
# @importFrom raster beginCluster clusterR endCluster
# @importMethodsFrom raster raster
#' @export

exposure_map <- function(dep, parallel = FALSE, n, filename = "", ...){

  ras <- dep@ras
  d <- dep@df
  filename <- trim(filename)
  if (!canProcessInMemory(ras) && filename == '') {
    filename <- rasterTmpFile()
  }

  f1 <- function(x) (x %*% d) #/ length(d)

  if(parallel) {
    beginCluster(n)
    exp.ras <- clusterR(ras, fun = .calc, args = list(fun = f1, forceapply = T, names = "Exposure"), filename = filename, ...)
    endCluster()
  } else {
    exp.ras <- .calc(ras, fun = f1, forceapply = T, filename = filename, names = "Exposure", ...)
  }

  return(exp.ras)
}
