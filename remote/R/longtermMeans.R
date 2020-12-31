#' Calculate long-term means from a 'RasterStack'
#' 
#' @description 
#' Calculate long-term means from an input 'RasterStack' (or 'RasterBrick') 
#' object. Ideally, the number of input layers should be divisable by the 
#' supplied \code{cycle.window}. For instance, if \code{x} consists of monthly 
#' layers, \code{cycle.window} should be a multiple of 12.
#' 
#' @param x A 'RasterStack' (or 'RasterBrick') object.
#' @param cycle.window 'integer'. See \code{\link{deseason}}.
#' 
#' @return
#' If \code{cycle.window} equals \code{nlayers(x)} (which obviously doesn't make 
#' much sense), a 'RasterLayer' object; else a 'RasterStack' object.
#' 
#' @author 
#' Florian Detsch
#' 
#' @seealso 
#' \code{\link{deseason}}. 
#' 
#' @examples 
#' data("australiaGPCP")
#' 
#' longtermMeans(australiaGPCP)
#' 
#' @export longtermMeans
#' @name longtermMeans
longtermMeans <- function(x, cycle.window = 12L) {

  ## raster to matrix  
  mat <- raster::as.matrix(x)
  
  ## long-term means
  mat_ltm <- monthlyMeansC(mat, 12)
  
  ## insert values
  rst_ltm <- x[[1:(raster::nlayers(x) / cycle.window)]]
  rst_ltm <- raster::setValues(rst_ltm, NA)
  
  rst_ltm <- raster::setValues(rst_ltm, mat_ltm)
  return(rst_ltm)
}