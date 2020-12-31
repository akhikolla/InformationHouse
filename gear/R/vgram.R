#' Empirical variogram
#'
#' \code{vgram} calculates an empirical variogram.  Note that, by convention,
#' the empirical variogram actually estimates the semivariogram, not the
#' theoretical variogram (which is twice the semivariogram).
#'
#' Note that the directions may be different from other packages
#' (e.g., \code{gstat} or \code{geoR} packages) because those
#' packages calculate angles clockwise from the y-axis, which is
#' a convention frequently seen in geostatistics
#' (e.g., the GSLIB software library).
#'
#' Additionally, note that calculating the empirical variogram for
#' the residuals of lm(response ~ 1) will produce identical results
#' to simply computing the sample semivariogram from the original response.
#' However, if a trend is specified (the righthand side of ~ has non-trival
#' covariates), then the empirical variogram of the residuals will differ
#' from that of the original response.  A trend should be specified when
#' the mean is non-stationary over the spatial domain.
#'
#' @param formula A formula describing the relationship between the response and any covariates of interest, e.g., response ~ 1.  The variogram is computed for the residuals of the linear model \code{lm(formula, data)}.
#' @param data A \code{data.frame}, \code{SpatialPointsDataFrame}, \code{SpatialPixelsDataFrame}, or \code{SpatialGridDataFrame} object.
#' @param coordnames The columns of \code{data} containing the spatial coordinates,
#' provided as a formula (e.g., ~ x + y),
#' column numbers (e.g., c(1, 2)), or column names (e.g., c("x", "y")).
#' The default is NULL, but this must be specified if \code{data} is of class \code{data.frame}.
#' @param nbins The number of bins (tolerance regions) to use when estimating the sample semivariogram.
#' @param maxd The maximum distance used when calculating the semivariogram.  Default is NULL, in which case half the maximum distance between coordinates is used.
#' @param angle A single value (in degrees) indicating the starting direction for a directional variogram.  The default is 0.
#' @param ndir The number of directions for which to calculate a sample semivariogram.  The default is 1, meaning calculate an omnidirection semivariogram.
#' @param type The name of the estimator to use in the estimation process.  The default is "standard", the typical method-of-moments estimator.  Other options include "cressie" for the robust Cressie-Hawkins estimator, and "cloud" for a semivariogram cloud based on the standard estimator.  If "cloud" is specified, the \code{nbins} argument is ignored.
#' @param npmin The minimum number of pairs of points to use in the semivariogram estimator.  For any bins with fewer points, the estimate for that bin is dropped.
#' @param longlat A logical indicating whether Euclidean (\code{FALSE}) or Great Circle distance (WGS84 ellipsoid) (\code{longlat = TRUE}) should be used.  Default is \code{FALSE}.
#' @param verbose Print computation information.  Default is \code{TRUE}.
#' @param coords A deprecated argument.
#'
#' @return Returns an \code{evgram} object with components:
#' @author Joshua French
#' @export
#' @examples
#' data(co)
#' v = vgram(Al ~ 1, co, ~ easting + northing)
#' plot(v)
#' v2 = vgram(Al ~ 1, co, c("easting", "northing"), angle = 22.5, ndir = 4)
#' plot(v2)
vgram = function(formula, data, coordnames = NULL, nbins = 10,
                 maxd = NULL, angle = 0, ndir = 1,
                 type = "standard", npmin = 2,
                 longlat = FALSE,
                 verbose = TRUE,
                 coords = NULL) {
  warning("vgram has been deprecated in favor of evgram.  Please update your code to use evgram function.")
  if (!is.null(coords)) {
    warning("coords has been deprecated in favor of coordnames. Please update your code.")
    coordnames = coords
  }
  evgram(formula = formula, data = data, coordnames = coordnames, nbins = nbins,
         maxd = maxd, angle = angle, ndir = ndir, type = type,
         npmin = npmin, longlat = longlat, verbose = verbose)
}

