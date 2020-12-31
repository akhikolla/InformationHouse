#' Dynamic Minimum Spanning Tree spatial scan test
#'
#' \code{dmst.test} implements the dynamic Minimum Spanning
#' Tree scan test of Assuncao et al. (2006). Starting with a
#' single region as a current zone, new candidate zones are
#' constructed by combining the current zone with the
#' connected region that maximizes the resulting likelihood
#' ratio test statistic.  This procedure is repeated until
#' the population or distance upper bounds are reached.  The
#' same procedure is repeated for each region.  The clusters
#' returned are non-overlapping, ordered from most
#' significant to least significant. The first cluster is
#' the most likely to be a cluster. If no significant
#' clusters are found, then the most likely cluster is
#' returned (along with a warning).
#'
#' The maximum intercentroid distance can be found by
#' executing the command:
#' \code{sp::spDists(as.matrix(coords), longlat = longlat)},
#' based on the specified values of \code{coords} and
#' \code{longlat}.
#'
#' @inheritParams flex.test
#' @inheritParams scan.test
#' @param ubd A proportion in (0, 1].  The distance of
#'   potential clusters must be no more than \code{ubd * m},
#'   where \code{m} is the maximum intercentroid distance
#'   between all coordinates.
#' @return Returns a \code{smerc_cluster} object.
#' @author Joshua French
#' @export
#' @seealso \code{\link{print.smerc_cluster}},
#' \code{\link{summary.smerc_cluster}},
#' \code{\link{plot.smerc_cluster}},
#' \code{\link{scan.stat}}, \code{\link{scan.test}}
#' @references Assuncao, R.M., Costa, M.A., Tavares, A. and
#'   Neto, S.J.F. (2006). Fast detection of arbitrarily
#'   shaped disease clusters, Statistics in Medicine, 25,
#'   723-742.  <doi:10.1002/sim.2411>
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = dmst.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, w = nyw,
#'                 alpha = 0.12, longlat = TRUE,
#'                 nsim = 2, ubpop = 0.05, ubd = 0.1)
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
dmst.test = function(coords, cases, pop, w,
                   ex = sum(cases) / sum(pop) * pop,
                   nsim = 499, alpha = 0.1, ubpop = 0.5,
                   ubd = 1, longlat = FALSE, cl = NULL) {
  # sanity checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha,
                      nsim + 1, ubpop, longlat, FALSE,
                      k = 1, w = w)

  coords = as.matrix(coords) # ensure proper format
  N = nrow(coords) # number of regions
  ty = sum(cases) # sum of all cases

  # intercentroid distances
  d = sp::spDists(as.matrix(coords), longlat = TRUE)
  # upperbound for distance between centroids in zone
  max_pop = ubpop * sum(pop)
  # find all neighbors from each starting zone within distance upperbound
  nn = nndist(d, ubd)

  all_zones = mst.all(nn, cases = cases, pop = pop, w = w,
                      ex = ex, ty = ty, max_pop = max_pop,
                      type = "all", nlinks = "one", early = FALSE,
                      cl = cl, progress = FALSE)
  # extract relevant information
  nn2 = lapply(all_zones, getElement, name = "locids")
  tobs = unlist(lapply(all_zones, getElement, name = "loglikrat"),
                use.names = FALSE)
  # determine distinct zones
  wdup = nndup(nn2, N)

  # remove zones with a test statistic of 0 or don't have
  # min number of cases or are duplicted
  w0 = which(tobs == 0 | wdup)

  # determine zones
  zones = nn2zones(nn2)

  # remove zones with a test statistic of 0
  zones = zones[-w0]
  tobs = tobs[-w0]

  # compute test statistics for simulated data
  if (nsim > 0) {
    message("computing statistics for simulated data:")
    tsim = dmst.sim(nsim = nsim, nn = nn, ty = ty,
                    ex = ex, w = w, pop = pop,
                    max_pop = max_pop, cl = cl)
    pvalue = mc.pvalue(tobs, tsim)
  } else {
    pvalue = rep(1, length(tobs))
  }

  # determine which potential clusters are significant
  sigc = which(pvalue <= alpha, useNames = FALSE)

  # if there are no significant clusters, return most likely cluster
  if (length(sigc) == 0) {
    sigc = which.max(tobs)
    warning("No significant clusters.  Returning most likely cluster.")
  }

  # significant, ordered, non-overlapping clusters and
  # information
  pruned = sig_noc(tobs = tobs, zones = zones,
                   pvalue = pvalue, alpha = alpha,
                   order_by = "tobs")

  smerc_cluster(tobs = pruned$tobs, zones = pruned$zones,
                pvalue = pruned$pvalue, coords = coords,
                cases = cases, pop = pop, ex = ex,
                longlat = longlat, method = "dynamic minimum spanning tree",
                rel_param = list(type = "poisson",
                                 simdist = "multinomial",
                                 nsim = nsim,
                                 ubpop = ubpop,
                                 ubd = ubd),
                alpha = alpha,
                w = w, d = NULL)
}
