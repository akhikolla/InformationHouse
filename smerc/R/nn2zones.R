#' Convert nearest neighbors list to zones
#'
#' \code{nn2zones} converts a list of nearest neighbors to
#' a list of zones.  The list of nearest neighbors will
#' come from functions such as \code{\link{nnpop}} or
#' \code{knn}.
#'
#' @param nn A list of nearest neighbors
#'
#' @return A list of zones
#' @export
#'
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(x, y))
#' nn = knn(coords, k = 2)
#' nn2zones(nn)
nn2zones = function(nn) {
  unlist(lapply(nn, function(x) {
    sapply(seq_along(x), function(i) {
      x[seq_len(i)]
    })
  }), recursive = FALSE)
}
