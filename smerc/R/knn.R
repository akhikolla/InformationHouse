#' K nearest neighbors
#'
#' \code{knn} returns the k nearest neighbors of the
#' n coordinates in \code{coords}.  The nearest neighbors
#' are constructed to be self-inclusive, i.e., an
#' observations is its closest neighbor.
#'
#' @inheritParams flex.test
#' @param d An n by n distance matrix.  If provided,
#' this is used instead of computing \code{d} based on
#' \code{coords} and \code{longlat}.
#'
#' @return An \eqn{n \times k} matrix of nearest neighbors.
#' @export
#'
#' @examples
#' data(nydf)
#' coords = nydf[,c("longitude", "latitude")]
#' knn(coords, longlat = TRUE, k = 4)
knn = function(coords, longlat = FALSE, k = 1, d = NULL) {
  arg_check_knn(coords, longlat, k, d)
  if (is.null(d)) {
    d = sp::spDists(as.matrix(coords), longlat = longlat)
  }
  # return list of sorted nn indices for each row of d
  lapply(seq_len(nrow(d)), function(i) order(d[i, ])[seq_len(k)])
}

arg_check_knn = function(coords, longlat, k, d) {
  if (is.null(dim(coords))) {
    stop("coords must be matrix-like")
  }
  if (length(longlat) != 1 | !is.logical(longlat)) {
    stop("longlat must be a logical value")
  }
  n = nrow(coords)
  if (k < 1 | k > n) {
    stop("k must be between 1 and n")
  }
  if (!is.null(d)) {
    if (nrow(d) != n | ncol(d) != n) {
      stop("nrow(coords) must match nrow(d) and ncol(d)")
    }
  }
}
