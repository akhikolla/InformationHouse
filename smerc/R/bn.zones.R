#' Determine case windows (circles)
#'
#' \code{bn.zones} determines the case windows (circles) for
#' the Besag-Newell method.
#'
#' Using the distances provided in \code{d}, for each
#' observation, the nearest neighbors are included in
#' increasingly larger windows until at least \code{cstar}
#' cases are included in the window.  Each row of \code{d}
#' is matched with the same position in \code{cases}.
#' @inheritParams nnpop
#' @param cases A vector of length \eqn{n} containing the
#'   observed number of cases for the \eqn{n} region
#'   centroids.
#' @param cstar A non-negative integer indicating the
#'   minimum number of cases to include in each window.
#' @return Returns the indices of the regions in each case
#'   window as a list.  For each element of the list, the
#'   indices are ordered from nearest to farthest from each
#'   centroid (and include the starting region).
#' @references Besag, J. and Newell, J.  (1991). The
#'   detection of clusters in rare diseases, Journal of the
#'   Royal Statistical Society, Series A, 154, 327-333.
#' @export
#' @author Joshua French
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' d = sp::spDists(coords, longlat = FALSE)
#' cwins = bn.zones(d, cases = nydf$cases, cstar = 6)
bn.zones = function(d, cases, cstar) {
  # order distances for each region
  # results has each column showing order of indices from
  # smallest to largest distance
  od = apply(d, 2, order)

  # for each row of ordered distance matrix sum the
  # cumulative number of cases for the expanding collection
  # of regions return the smallest collection of regions for
  # which the cumulative population is less than the desired
  # proportion of the total popuation
  wins = apply(od, 2, FUN = function(x) cumsum(cases[x]))
  size_cwin = apply(wins >= cstar, 2, which.max)
  return(sapply(seq_along(size_cwin),
                function(i) od[seq_len(size_cwin[i]), i]
  ))
}

#' @rdname bn.zones
#' @export
casewin = bn.zones
