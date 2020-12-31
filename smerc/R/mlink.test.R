#' Maximum Linkage spatial scan test
#'
#' \code{mlink.test} implements the Maximum Linkage spatial
#' scan test of Costa et al. (2012). Starting with a single
#' region as a current zone, new candidate zones are
#' constructed by combining the current zone with the
#' connected region that maximizes the resulting likelihood
#' ratio test statistic, with the added constraint that the
#' region has the maximum connections (i.e., shares a border
#' with) with the regions in the current zone.  This
#' procedure is repeated until the population or distance
#' upper bounds constraints are reached.  The same procedure
#' is repeated for each region.  The clusters returned are
#' non-overlapping, ordered from most significant to least
#' significant. The first cluster is the most likely to be a
#' cluster. If no significant clusters are found, then the
#' most likely cluster is returned (along with a warning).
#'
#' The maximum intercentroid distance can be found by
#' executing the command:
#' \code{sp::spDists(as.matrix(coords), longlat = longlat)},
#' based on the specified values of \code{coords} and
#' \code{longlat}.
#'
#' @inheritParams dmst.test
#' @return Returns a \code{smerc_cluster} object.
#' @author Joshua French
#' @export
#' @seealso \code{\link{print.smerc_cluster}},
#' \code{\link{summary.smerc_cluster}},
#' \code{\link{plot.smerc_cluster}},
#' \code{\link{scan.stat}}, \code{\link{scan.test}}
#' @references Costa, M.A. and Assuncao, R.M. and Kulldorff,
#'   M. (2012) Constrained spanning tree algorithms for
#'   irregularly-shaped spatial clustering, Computational
#'   Statistics & Data Analysis, 56(6), 1771-1783.
#'   <doi:10.1016/j.csda.2011.11.001>
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = mlink.test(coords = coords, cases = floor(nydf$cases),
#'                  pop = nydf$pop, w = nyw,
#'                  alpha = 0.12, longlat = TRUE,
#'                  nsim = 2, ubpop = 0.05, ubd = 0.1)
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
mlink.test = function(coords, cases, pop, w,
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
  # upperbound for zone population
  max_pop = ubpop * sum(pop)
  # find all neighbors from each starting zone within distance upperbound
  nn = nndist(d, ubd)

  all_zones = mst.all(nn, cases = cases, pop = pop, w = w,
                      ex = ex, ty = ty, max_pop = max_pop,
                      type = "all", nlinks = "max", early = FALSE,
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
    tsim = mlink.sim(nsim = nsim, nn = nn, ty = ty,
                  ex = ex, w = w, pop = pop,
                  max_pop = max_pop, cl = cl)
    pvalue = mc.pvalue(tobs, tsim)
  } else {
    pvalue = rep(1, length(tobs))
  }

  # significant, ordered, non-overlapping clusters and
  # information
  pruned = sig_noc(tobs = tobs, zones = zones,
                   pvalue = pvalue, alpha = alpha,
                   order_by = "tobs")

  smerc_cluster(tobs = pruned$tobs, zones = pruned$zones,
                pvalue = pruned$pvalue, coords = coords,
                cases = cases, pop = pop, ex = ex,
                longlat = longlat, method = "maximum linkage",
                rel_param = list(type = "poisson",
                                 simdist = "multinomial",
                                 nsim = nsim,
                                 ubpop = ubpop,
                                 ubd = ubd),
                alpha = alpha,
                w = w, d = NULL)
}
