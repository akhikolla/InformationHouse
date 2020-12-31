#' Minimum spanning tree for all regions
#'
#' \code{mst.all} finds the set of connected regions that
#' maximize the spatial scan statistic (the likelihood ratio
#' test statistic) from each starting region, subject to
#' relevant constraints.  The function can be used to
#' construct candidate zones for the dynamic minimum
#' spanning tree (dmst), early stopping dynamic minimum
#' spanning tree (edmst), double connected spatial scan test
#' (dc), and maximum linkage (mlink) spatial scan test.
#'
#' This function is not intended to be used by users
#' directly. Consequently, it prioritizes efficiency over
#' user friendliness.
#'
#' \code{type} is a character vector indicating what should
#' be returned by the function.  If \code{type = "maxonly"},
#' then the maximum test statistic from each starting region
#' is returned .  If \code{type = "pruned"}, the function
#' returns a list that includes the location ids, test
#' statistic, total cases, expected cases, and total
#' population for the zone with the maximum test statistic
#' for each starting region.  If \code{type = "all"}, the
#' function returns a list of lists that includes the
#' location ids, test statistic, total cases, expected
#' cases, and total population for the sequence of candidate
#' zones associated with each starting region.
#'
#' If \code{nlinks = "one"}, then a region only needs to be
#' connected to one other region in the current zone to be
#' considered for inclusion in the next zone.  If
#' \code{nlinks = "two"}, then the region must be connected
#' to at least two other regions in the current zone.  If
#' \code{nlinks = "max"}, then only regions with the maximum
#' number of connections to the current zone are considered
#' for inclusion in the next zone.
#'
#' @param neighbors A list containing the vector of
#'   neighbors for each region (in ascending order of
#'   distance from the region).  The starting region itself
#'   is included among the neighbors.
#' @inheritParams scan.test
#' @inheritParams scan.stat
#' @inheritParams uls.test
#' @param max_pop The population upperbound (in total
#'   population) for a candidate zone.
#' @param type One of \code{"maxonly"}, \code{"pruned"}, or
#'   \code{"all"}.  See Details.
#' @param nlinks A character vector.  The options are
#'   \code{"one"}, \code{"two"}, or \code{"max"}.  See
#'   Details.
#' @param early A logical value indicating whether the
#'   "early" stopping criterion should be used.  If
#'   \code{TRUE}, each sequence is stopped when the next
#'   potential zone doesn't produce a test statistic larger
#'   than the current zone. The default is \code{FALSE}.
#' @param progress A logical value indicating whether a
#'   progress bar should be displayed.  The default is
#'   \code{TRUE}.
#' @return Returns a list of relevant information.  See
#'   Details.
#' @author Joshua French
#' @export
#' @references Assuncao, R.M., Costa, M.A., Tavares, A. and
#'   Neto, S.J.F. (2006). Fast detection of arbitrarily
#'   shaped disease clusters, Statistics in Medicine, 25,
#'   723-742.  <doi:10.1002/sim.2411>
#'
#'   Costa, M.A. and Assuncao, R.M. and Kulldorff, M. (2012)
#'   Constrained spanning tree algorithms for
#'   irregularly-shaped spatial clustering, Computational
#'   Statistics & Data Analysis, 56(6), 1771-1783.
#'   <doi:10.1016/j.csda.2011.11.001>
#' @examples
#' # load data
#' data(nydf)
#' data(nyw)
#'
#' # create relevant data
#' coords = nydf[,c("longitude", "latitude")]
#' cases = floor(nydf$cases)
#' pop = nydf$population
#' w = nyw
#' ex = sum(cases)/sum(pop)*pop
#' ubpop = 0.5
#' ubd = 0.5
#' ty = sum(cases)   # total number of cases
#' # intercentroid distances
#' d = sp::spDists(as.matrix(coords), longlat = TRUE)
#' # upperbound for population in zone
#' max_pop = ubpop * sum(pop)
#' # upperbound for distance between centroids in zone
#' max_dist = ubd * max(d)
#' # create list of neighbors for each region
#' # (inclusive of region itself)
#' all_neighbors = nndist(d, ubd)
#' # find the dmst max zone
#' \dontrun{
#' out = mst.all(all_neighbors, cases, pop, w, ex, ty, max_pop,
#'               type = "maxonly")
#' head(out)
#'
#' out = mst.all(all_neighbors, cases, pop, w, ex, ty, max_pop,
#'               type = "pruned")
#' head(out)
#' }
mst.all = function(neighbors, cases, pop, w, ex, ty, max_pop,
                   type = "maxonly", nlinks = "one",
                   early = FALSE, cl = NULL, progress = FALSE) {
  if (!progress) {
    max_zones = lapply(seq_along(cases), function(i) {
      mst.seq(i, neighbors[[i]], cases,
              pop, w, ex, ty, max_pop, type,
              nlinks, early)
    })
  } else {
    max_zones = pbapply::pblapply(seq_along(cases), function(i) {
      mst.seq(i, neighbors[[i]], cases,
              pop, w, ex, ty, max_pop, type,
              nlinks, early)
    }, cl = cl)
  }
  if (type == "maxonly") {
    return(unlist(max_zones, use.names = FALSE))
  } else {
    max_zones
  }
}
