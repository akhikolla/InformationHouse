#' Determine sequence of fast subset scan zones
#'
#' \code{fast.zones} determines the unique zones obtained by
#' implementing the fast subset scan method of Neill (2012).
#'
#' The \code{simple} argument determines the formatting of
#' the returned zones.  If \code{simple = TRUE}, then a
#' vector containing the sequential indices of the regions
#' in each successive zones is returned.  If \code{simple =
#' FALSE}, then the complete list of all zones is returned
#' (which is the standard format of most of the other
#' \code{*.zones} functions.
#'
#' The zones returned must have a total population less than
#' \code{ubpop * sum(pop)} of all regions in the study area.
#'
#' @inheritParams uls.test
#' @param simple A logical value indicating whether a simple
#'   version of the fast zones should be returned.  See
#'   Details.
#' @return Returns a vector of regions to sequentially and
#'   cumulatively consider for clustering.
#' @author Joshua French
#' @export
#' @references Neill, D. B. (2012), Fast subset scan for
#'   spatial pattern detection. Journal of the Royal
#'   Statistical Society: Series B (Statistical
#'   Methodology), 74: 337-360.
#'   <doi:10.1111/j.1467-9868.2011.01014.x>
#' @examples
#' data(nydf)
#' cases = nydf$cases
#' pop = nydf$pop
#' # compare output format
#' fast.zones(cases, pop, ubpop = 0.05)
#' fast.zones(cases, pop, ubpop = 0.05, simple = FALSE)
fast.zones = function(cases, pop, ubpop = 0.5, simple = TRUE) {
  if (length(cases) != length(pop)) {
    stop("length(cases) != length(pop)")
  }
  if (length(ubpop) != 1 | !is.numeric(ubpop)) {
    stop("ubpop should be a single number")
  }
  if (ubpop <= 0 | ubpop > 1) {
    stop("ubpop not in (0, 1]")
  }
  if (length(simple) != 1 | !is.logical(simple)) {
    stop("simple must be a logical value")
  }

  # order rates from largest to smallest
  or = order(cases / pop, decreasing = TRUE)
  max_pop = sum(pop) * ubpop

  if (simple) {
    return(or[cumsum(pop[or]) <= max_pop])
  } else {
    or = or[cumsum(pop[or]) <= max_pop]
    return(lapply(seq_along(or), function(i) or[seq_len(i)]))
  }
}
