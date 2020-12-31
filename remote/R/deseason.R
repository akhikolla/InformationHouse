if ( !isGeneric("deseason") ) {
  setGeneric("deseason", function(x, ...)
    standardGeneric("deseason"))
}
#' Create seasonal anomalies
#' 
#' @description
#' The function calculates anomalies of a RasterStack by supplying a 
#' suitable seasonal window. E. g. to create monthly anomalies of a 
#' raster stack of 12 layers per year, use \code{cycle.window = 12}.
#' 
#' @param x An object of class 'RasterStack' (or 'RasterBrick') or, 
#' alternatively, a 'numeric' time series.
#' @param cycle.window Integer. The window for the creation of the anomalies.
#' @param use.cpp Logical. Determines whether or not to use \strong{Rcpp} 
#' functionality, defaults to \code{TRUE}. Only applies if \code{x} is a 
#' 'RasterStack' (or 'RasterBrick') object.
#' @param filename \code{character}. Output filename (optional).
#' @param ... Additional arguments passed on to \code{\link{writeRaster}}, only 
#' considered if \code{filename} is specified.
#' 
#' @return If \code{x} is a 'RasterStack' (or 'RasterBrick') object, a 
#' deseasoned 'RasterStack'; else a deseasoned 'numeric' vector.
#' 
#' @seealso
#' \code{\link{anomalize}}, \code{\link{denoise}}
#' 
#' @export deseason
#' @name deseason
#' 
#' @examples 
#' data("australiaGPCP")
#' 
#' aus_dsn <- deseason(australiaGPCP, 12)
#' 
#' opar <- par(mfrow = c(1,2))
#' plot(australiaGPCP[[1]], main = "original")
#' plot(aus_dsn[[1]], main = "deseasoned")
#' par(opar)


################################################################################
### function using 'RasterStack' or 'RasterBrick' ##############################
#' @aliases deseason,RasterStackBrick-method
#' @rdname deseason
setMethod("deseason",
          signature(x = "RasterStackBrick"),
          function(x, 
                   cycle.window = 12L,
                   use.cpp = FALSE,
                   filename = "", 
                   ...) {
            
            if (use.cpp) {
              ## raster to matrix
              mat <- as.matrix(x)
              
              ## deseasoning
              mat_mv <- monthlyMeansC(mat, cycle.window)
              x_mv <- x[[1:cycle.window]]
              x_mv <- setValues(x_mv, values = mat_mv)
            } else {
              # Calculate layer averages based on supplied seasonal window
              x_mv <- raster::stack(rep(lapply(1:cycle.window, function(i) {
                raster::calc(x[[seq(i, raster::nlayers(x), cycle.window)]], fun = mean)
              }), raster::nlayers(x) / cycle.window))
            }
            
            # Subtract monthly averages from actually measured values
            x_dsn <- x - x_mv
            
            # Write to file (optional)
            if (filename != "")
              x_dsn <- raster::writeRaster(x_dsn, filename = filename, ...)
            
            # Return output
            return(x_dsn)
          })


################################################################################
### function using 'numeric' ###################################################
#' @aliases deseason,numeric-method
#' @rdname deseason
setMethod("deseason",
          signature(x = "numeric"),
          function(x, 
                   cycle.window = 12L) {
            
            ## calculate long-term mean values
            x_mv <- sapply(1:cycle.window, function(i) {
              val <- x[seq(i, length(x), cycle.window)]
              mean(val, na.rm = TRUE)
            })
            
            ## create anomalies
            x_dsn <- x - x_mv
            
            # Return output
            return(x_dsn)
          })
