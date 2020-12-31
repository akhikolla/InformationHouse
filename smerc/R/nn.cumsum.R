#' Cumulative sum over nearest neighbors
#'
#' \code{nn.cumsum} computes the cumulative sum of \code{y}
#' for the sequences of indices in each element of the list
#' contained in \code{nn}.
#'
#' @param nn A list of nearest neighbors in the format
#'   produced by \code{\link{nnpop}}.
#' @param y A numeric vector of values to be summed over.
#' @param simplify A logical value indicating whether the
#'   results should be simplified to a numeric vector.  The
#'   default is \code{TRUE}.
#'
#' @return A vector or list, depending on the value of
#'   \code{simplify}.
#' @export
#'
#' @examples
#' # show nn.cumsum example for a circular scan setting
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' cases = floor(nydf$cases)
#' d = sp::spDists(coords, longlat = TRUE)
#' # compute circular nearest neigbhors
#' nn = nnpop(d, pop = nydf$pop, ubpop = 0.1)
#' # compute cumulative sums over all nn
#' cnn = nn.cumsum(nn, cases)
#' # compute cumulative sums over just the first set of nn
#' cnn1 = cumsum(cases[nn[[1]]])
#' # check equality
#' all.equal(cnn1, cnn[seq_along(cnn1)])
nn.cumsum = function(nn, y, simplify = TRUE) {
  if (!is.list(nn)) stop("nn must be a list")
  if (!is.numeric(y)) stop("y must be a numeric vector")
  if (!is.logical(simplify)) stop("simplify must be a logical value")
  if (simplify) {
    unlist(lapply(nn, function(x) cumsum(y[x])), use.names = FALSE)
  } else {
    lapply(nn, function(x) cumsum(y[x]))
  }
}
