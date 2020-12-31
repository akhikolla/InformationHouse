#' Compute region weights for \code{cepp.test}
#'
#' @inheritParams cepp.test
#' @param nn A list of nearest neighbors produced by
#' \code{\link{casewin}}.
#'
#' @return A list with elements related to the weight
#' each nearest neighbor region will have in the
#' corresponding weighted sum used to compute the test
#' statistic
#' @export
#'
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(x, y))
#' pop = nydf$pop
#' # intercentroid distances
#' d = sp::spDists(coords)
#' # find smallest windows with cumulative population of
#' # at least n* = 1000
#' nn = casewin(d, pop, 1000)
#' # compute weights
#' w = cepp.weights(nn, pop, 1000)
cepp.weights = function(nn, pop, nstar) {
  if (length(nn) != length(pop)) {
    stop("length(nn) != length(pop)")
  }
  if (!is.numeric(pop)) {
    stop("pop must be numeric")
  }
  # determine cases and population in each window
  nn_cpop = nn.cumsum(nn, pop, simplify = FALSE)
  # prev cpop
  # determine next to last cpop for each element of nn_cpop
  # (or last of there is only one value in the element)
  n2l_cpop = sapply(nn_cpop, function(x) utils::tail(x, n = 2)[1],
                    USE.NAMES = FALSE)
  # determine pop for last region for each element of nn
  lpop = sapply(nn, function(x) pop[utils::tail(x, n = 1)],
                USE.NAMES = FALSE)
  # number of neighbors (inclusive) for each window
  nnn = sapply(nn, length, USE.NAMES = FALSE)

  w = lapply(seq_along(nnn), function(i) {
    if (nnn[i] == 1) {
      return(nstar / lpop[i])
    } else {
      out = rep(1, nnn[i])
      out[nnn[i]] = (nstar - n2l_cpop[i]) / lpop[i]
      return(out)
    }
  })
  return(w)
}
