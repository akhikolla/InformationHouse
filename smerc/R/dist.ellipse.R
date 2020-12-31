#' Compute minor axis distance of ellipse
#'
#' \code{dist.ellipse} computes the length of the minor axis
#' needed for an ellipse of a certain \code{shape} and
#' \code{angle} to intersect each of the other coordinates
#' from a starting coordinate.
#'
#' @param coords An \eqn{N \times 2} matrix of coordinates
#' @param shape The ratio of the major axis to the minor
#'   axis of the ellipse
#' @param angle The angle of the ellipse in the range [0,
#'   180).
#'
#' @return A matrix of distances between each coordinate and
#'   all other coordinates (and itself).  Each row contains
#'   the distances for a coordinate.
#' @export
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("x", "y")])
#' d = dist.ellipse(coords, 4, 15)
dist.ellipse = function(coords, shape, angle) {
  arg_check_dist_ellipse(coords, shape, angle)
  rad = angle * pi / 180
  N = nrow(coords)

  sinr = rep(sin(rad), N)
  cosr = rep(cos(rad), N)
  t(sapply(seq_len(N), function(i) {
    coordsi = matrix(coords[i, ], byrow = TRUE,
                     ncol = 2, nrow = N)
    cdiff = (coords - coordsi)
    dx = cdiff[, 1]
    dy = cdiff[, 2]
    # formula based on https://math.stackexchange.com/
    # questions/426150/what-is-the-general-equation-of-the-
    # ellipse-that-is-not-in-the-origin-and-rotate
    m1 = (dx * cosr + dy * sinr) / shape
    m2 = (dx * sinr - dy * cosr)
    sqrt(m1 ^ 2 + m2 ^ 2)
  }))
}

#' Check argments if dist.ellipse
#' @return NULL
#' @export
#' @keywords internal
arg_check_dist_ellipse = function(coords, shape, angle) {
  if (!is.matrix(coords)) {
    stop("coords must be a matrix")
  }
  if (length(shape) != 1) {
    stop("shape should be a single value")
  }
  if (shape < 1) {
    stop("shape should be at least 1")
  }
  if (length(angle) != 1) {
    stop("angle should be a single value")
  }
  if (angle < 00 | angle >= 360) {
    stop("angle should be in [0, 360) degrees")
  }
}

#' Return all shapes and distances for each zone
#' @return NULL
#' @export
#' @keywords internal
all_shape_dists = function(s, na, coords) {
  angle = seq(90, 270, len = na + 1)[seq_len(na)]
  d = lapply(seq_along(angle), function(j) {
    dist.ellipse(coords, shape = s, angle = angle[j])
  })
  do.call(rbind, d)
}
