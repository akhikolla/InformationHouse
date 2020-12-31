#' Determine nearest neighbors based on maximum distance
#'
#' \code{nndist} determines the nearest
#' neighbors for a set of observations within a certain
#' radius.
#'
#' This function determines the nearest neighbors of each
#' centroid based on the intercentroid distance.  The number
#' of nearest neighbors is limited by the furthest distance
#' between the starting centroid and the farthest neighbor.
#'
#' @param d An \eqn{n\times n} square distance matrix
#'   containing the intercentroid distance between the
#'   \eqn{n} region centroids.
#' @inheritParams dmst.test
#'
#' @return Returns the indices of the nearest neighbors as a
#'   list.
#' @author Joshua French
#' @export
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' d = as.matrix(dist(coords))
#' nn = nndist(d, ubd = 0.01)
nndist = function(d, ubd) {
  if (is.null(dim(d))) stop("d must be matrix-like")
  if (nrow(d) != ncol(d)) stop("d must be square")
  if (length(ubd) != 1 | !is.numeric(ubd) | ubd <= 0) {
    stop("ubd must be a single positive value")
  }
  max_dist = ubd * max(d)
  # find all neighbors from each starting zone within distance upperbound
  lapply(seq_len(nrow(d)), function(i) unname(which(d[i, ] <= max_dist)))
}
