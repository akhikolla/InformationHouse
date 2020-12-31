#' Determine duplicates in nearest neighbor list
#'
#' \code{nndup} determines the indices of duplicated
#' elements for a nearest neighbors list created by a
#' function such as \code{\link{nnpop}} or
#' \code{\link{knn}}.  The indices are related to the list
#' returned by \code{\link{nn2zones}}.
#' @param nn A list of nearest neighbors.
#' @param N The largest value in \code{nn}.
#'
#' @return A logical vector of indicating duplicate indices.
#' @export
#'
#' @examples
#' nn = list(1:3, c(2:1, 4))
#' nndup(nn, 4)
nndup = function(nn, N = max(unlist(nn))) {
  # get prime numbers
  pri = randtoolbox::get.primes(N)
  # determine the el
  duplicated(unlist(lapply(nn, function(x) {
    cumsum(log(pri[x]))
  })))
}
