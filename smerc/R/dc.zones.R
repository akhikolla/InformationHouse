#' Determine zones for the Double Connected scan test
#'
#' \code{dc.zones} determines the zones for the Double
#' Connected scan test (\code{\link{dc.test}}).  The
#' function returns the zones, as well as the associated
#' test statistic, cases in each zone, the expected number
#' of cases in each zone, and the population in each zone.
#'
#' Every zone considered must have a total population less
#' than \code{ubpop * sum(pop)}.  Additionally, the maximum
#' intercentroid distance for the regions within a zone must
#' be no more than \code{ubd * the maximum intercentroid
#' distance across all regions}.
#' @inheritParams dmst.test
#' @inheritParams mst.all
#' @inheritParams flex.zones
#' @return Returns a list with elements: \item{zones}{A list
#'   contained the location ids of each potential cluster.}
#'   \item{loglikrat}{The loglikelihood ratio for each zone
#'   (i.e., the log of the test statistic).}
#'   \item{cases}{The observed number of cases in each
#'   zone.} \item{expected}{The expected number of cases
#'   each zone.} \item{pop}{The total population in each
#'   zone.}
#' @author Joshua French
#' @references Costa, M.A. and Assuncao, R.M. and Kulldorff,
#'   M. (2012) Constrained spanning tree algorithms for
#'   irregularly-shaped spatial clustering, Computational
#'   Statistics & Data Analysis, 56(6), 1771-1783.
#'   <doi:10.1016/j.csda.2011.11.001>
#' @export
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' # find zone with max statistic starting from each individual region
#' all_zones = dc.zones(coords, cases = floor(nydf$cases),
#'                      nydf$pop, w = nyw, ubpop = 0.25,
#'                      ubd = .25, longlat = TRUE)
dc.zones = function(coords, cases, pop, w,
                    ex = sum(cases) / sum(pop) * pop,
                    ubpop = 0.5, ubd = 1, longlat = FALSE,
                    cl = NULL, progress = TRUE) {
  # sanity checking
  arg_check_dmst_zones(coords = coords, cases = cases,
                       pop = pop, w = w, ex = ex,
                       ubpop = ubpop, ubd = ubd,
                       longlat = longlat, type = "all",
                       progress = progress)
  # setup various arguments and such
  ty = sum(cases)   # total number of cases
  # intercentroid distances
  d = sp::spDists(as.matrix(coords), longlat = longlat)
  # upperbound for population in zone
  max_pop = ubpop * sum(pop)
  # find all neighbors from each starting zone within distance upperbound
  nn = nndist(d, ubd)

  out = mst.all(neighbors = nn, cases = cases, pop = pop,
                w = w, ex = ex, ty = ty, max_pop = max_pop,
                type = "all", nlinks = "two", early = TRUE,
                cl = cl, progress = progress)
  # return results in a list
  prep.mst(out)
}
