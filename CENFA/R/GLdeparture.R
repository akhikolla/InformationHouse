#' Climatic departure of reference study area
#'
#' This function is used to facilitate comparisons between species in the same
#' study area. It speeds up the computation of multiple departures by calculating
#' the global covariance matrix as a first step, which
#' can then be fed into the \code{\link{departure}} function as a first argument.
#' This saves the user from having to calculate the global covariance matrix for
#' each species, which can take quite a bit of time.
#'
#' @aliases print.GLdeparture show.GLdeparture
#'
#' @param x Raster* object of p historical climate layers
#' @param y Raster* object of p future climate layers, with the same names as
#'   \code{x}
#' @param center logical or numeric. If \code{TRUE}, centering is done by
#'   subtracting the layer means (omitting NAs), and if \code{FALSE}, no centering
#'   is done. If \code{center} is a numeric vector with length equal to the
#'   \code{nlayers(x)}, then each layer of \code{x} has the corresponding value
#'   from center subtracted from it
#' @param scale logical or numeric. If \code{TRUE}, scaling is done by dividing
#'   the (centered) layers of \code{x} by their standard deviations if center is
#'   \code{TRUE}, and the root mean square otherwise. If scale is \code{FALSE},
#'   no scaling is done. If scale is a numeric vector with length equal to
#'   \code{nlayers(x)}, each layer of \code{x} is divided by the corresponding
#'   value. Scaling is done after centering
#' @param filename character. Optional filename to save the RasterBrick output
#'   to file. If this is not provided, a temporary file will be created for large
#'   \code{x}
#' @param progress logical. If \code{TRUE}, messages and progress bar will be
#'   printed
#' @param parallel logical. If \code{TRUE} then multiple cores are utilized
#' @param n numeric. Number of CPU cores to utilize for parallel processing
#' @param cl optional cluster object
#' @param keep.open logical. If \code{TRUE} and \code{parallel = TRUE}, the
#'   cluster object will not be closed after the function has finished
#' @param ... Additional arguments for \code{\link[raster]{writeRaster}}
#'
#' @examples
#' gld <- GLdeparture(x = climdat.hist, y = climdat.fut)
#'
#' @return Returns an S4 object of class \code{GLcenfa} with the following components:
#' \describe{
#'   \item{global_difras}{Raster* \code{x} of p layers, possibly centered and scaled}
#'   \item{cov}{Global p x p covariance matrix}
#'   }
#'
#' @details
#' If there is too much correlation between the layers of \code{x}, the covariance
#' matrix will be singular, which will lead to later problems in computing the overall
#' departures of species. In this case, a warning will be issued, suggesting the
#' removal of correlated variables or a transformation of the data.
#'
#' @seealso \code{\link{departure}}
#' @export
# @importFrom raster compareRaster cellStats

setGeneric("GLdeparture", function(x, y, center = TRUE, scale = TRUE, filename = '', progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE, ...) {
  standardGeneric("GLdeparture")
})

#' @rdname GLdeparture
setMethod("GLdeparture",
          signature(x = "Raster", y = "Raster"),
          function(x, y, center = TRUE, scale = TRUE, filename = '', progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE, ...){

            if (!all.equal(names(x), names(y))) stop("historical and future raster layers do not match")
            if (!compareRaster(x, y)) stop("historical and future rasters resolutions or extent do not match")

            if (center || scale) {
              if (progress) cat("Scaling historical raster data...\n")
              means <- cellStats(x, mean)
              sds <- cellStats(x, sd)
              x <- parScale(x, center = center, scale = scale, filename = '', progress = progress, parallel = parallel, n = n, cl = cl, keep.open = keep.open)
              if (progress) cat("Scaling future raster data...\n")
              y <- parScale(y, center = means, scale = sds, filename = '', progress = progress, parallel = parallel, n = n, cl = cl, keep.open = keep.open)
            }

            filename <- trim(filename)
            if (parallel || !canProcessInMemory(x)) {
              if (filename == '') filename <- rasterTmpFile()
            }

            #x.dif <- abs(y - x)
            x.dif <- .overlay(x, y, fun = function(x, y) { return(abs(x - y)) }, forcefun = TRUE, filename = filename, names = names(x), ...)
            names(x.dif) <- names(x)



            if (progress) cat("Calculating global covariance matrix...\n")
            cov.mat <- parCov(x = x, parallel = parallel, n = n, progress = progress, cl = cl, keep.open = keep.open)
            tryCatch(solve(cov.mat, tol = 1e-10),
                     error = function(e) message("Warning: covariance matrix is not invertible. Consider removing correlated variables or transforming data and trying again."))

            GLdif <- methods::new("GLdeparture", global_difras = x.dif, cov = cov.mat)
            return(GLdif)
          }
)

#' @rdname GLdeparture
setMethod("GLdeparture",
          signature(x = "Raster", y = "missing"),
          function(x, filename = '', progress = FALSE, parallel = FALSE, n = 1, cl = NULL, keep.open = FALSE, ...){

              if (progress) cat("Calculating global covariance matrix...\n")
              cov.mat <- parCov(x = x, parallel = parallel, n = n, progress = progress, cl = cl, keep.open = keep.open)
              GLdif <- methods::new("GLdeparture", global_difras = x, cov = cov.mat)
              return(GLdif)

          }
)
