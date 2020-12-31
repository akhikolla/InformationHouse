#' Compute distance for geographic coordinates
#'
#' \code{geodist} computes the distance between the
#' coordinates in \code{coords1} and \code{coords2}. If
#' \code{coords2} isn't supplied, then the distances are
#' computed between the coordinates in \code{coords1} alone.
#' Otherwise, the pairwise distances between then points in
#' \code{coords1} and \code{coords2} is computed. If
#' \code{longlat = TRUE}, then the great circle distance is
#' computed.
#'
#' The algorithm used when \code{longlat = TRUE} is a C++
#' port of the C code written by Roger Bivand for
#' \code{\link[sp]{spDists}}, which appears to be based on a
#' special case of the Vincenty formula with a slight
#' correction based on the WGS84 flattening constant. See
#' \url{https://en.wikipedia.org/wiki/Great-circle_distance}.
#'
#' @param coords1 A two-dimensional matrix of coordinates.
#' @param coords2 A two-dimensional matrix of coordinates.
#' @param longlat A logical value indicating whether
#'   Euclidean distance (\code{longlat = FALSE}) or great
#'   circle distance (\code{longlat = FALSE}) should be
#'   computed. The default is \code{longlat = FALSE}.
#'
#' @return A matrix of distances
#' @export
#' @examples
#' coords = matrix(runif(10), ncol = 2)
#' d = geodist(coords)
#' all.equal(d, as.matrix(dist(coords)), check.attributes = FALSE)
geodist = function(coords1, coords2, longlat = FALSE) {
  if (missing(coords2) & !longlat) {
    eucdist1(x = coords1[,1], y = coords1[,2],
            eps = .Machine$double.eps)
  } else if (missing(coords2) & longlat) {
    gcdist1(lon = coords1[,1], lat = coords1[,2], eps = .Machine$double.eps)
  } else if (!missing(coords2) & !longlat) {
    eucdist2(x1 = coords1[,1], y1 = coords1[,2],
             x2 = coords2[,1], y2 = coords2[,2],
             eps = .Machine$double.eps)
  } else if (!missing(coords2) & longlat) {
    gcdist2(lon1 = coords1[,1], lat1 = coords1[,2],
            lon2 = coords2[,1], lat2 = coords2[,2],
            eps = .Machine$double.eps)
  }
}
