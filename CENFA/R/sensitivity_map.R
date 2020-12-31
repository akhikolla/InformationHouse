#' Create a sensitivity map
#'
#' Creates a sensitivity map of species habitat from a \code{cnfa} object.
#'
#' @param cnfa Object of class \code{cnfa}
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Number of cores to use for calculation
#' @param filename character. Output filename (optional)
#' @param ... Additional arguments for file writing as for \code{\link[raster]{writeRaster}}
#'
#' @details The values of the sensitivity raster are calculated by centering the
#' habitat's climate data around the marginality factor \strong{m}
#' and projecting onto the sensitivity factor \strong{s}, given by the formula
#'
#'  \eqn{\sigma} = |\bold{S} - \bold{m}|\bold{s}.
#'
#' @examples
#' mod1 <- cnfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#' sens.map <- sensitivity_map(mod1)
#'
#' @return A RasterLayer of sensitivity values
#'
#' @seealso \code{\link{cnfa}}, \code{\link{exposure_map}}, \code{\link{vulnerability_map}}
#'
#' @export

sensitivity_map <- function(cnfa, parallel = FALSE, n = 1, filename = "", ...){

  ras <- cnfa@ras
  m <- cnfa@mf
  s <- cnfa@sf
  U.inv <- solve(cnfa@co)
  filename <- trim(filename)
  if (!canProcessInMemory(ras) && filename == '') {
    filename <- rasterTmpFile()
  }

  f1 <- function(x) x %*% U.inv
  f2 <- function(x) (abs(x - m) %*%  s) / length(s)

  if(parallel) {
    beginCluster(n)
    S <- clusterR(ras, fun = calc, args = list(fun = f1, forceapply = T))
    sens.ras <- clusterR(S, fun = .calc, args = list(fun = f2, forceapply = T, names = "Sensitivity"), filename = filename, ...)
    endCluster()
  } else {
    S <- calc(ras, fun = f1, forceapply = T)
    sens.ras <- .calc(S, fun = f2, forceapply = T, filename = filename, names = "Sensitivity", ...)
  }

  return(sens.ras)
}
