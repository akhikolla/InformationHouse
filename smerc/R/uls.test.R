#' Upper Level Set Spatial Scan Test
#'
#' \code{uls.test} performs the Upper Level Set (ULS)
#' spatial scan test of Patil and Taillie (2004).  The test
#' is performed using the spatial scan test based on a fixed
#' number of cases.  The windows are based on the Upper
#' Level Sets proposed by Patil and Taillie (2004).  The
#' clusters returned are non-overlapping, ordered from most
#' significant to least significant.  The first cluster is
#' the most likely to be a cluster.  If no significant
#' clusters are found, then the most likely cluster is
#' returned (along with a warning).
#'
#' The ULS method has a special (and time consuming)
#' construction when the observed rates aren't unique.  This
#' is unlikely to arise for real data, except with observed
#' rates of 0, which are of little interest.  The method can
#' take substantially if this is considered.
#'
#' @param w A binary spatial adjacency matrix for the
#'   regions.
#' @param check.unique A logical value indicating whether a
#'   check for unique values should be determined.  The
#'   default is \code{FALSE}.  This is unlikely to make a
#'   practical different for most real data sets.
#' @inheritParams scan.test
#'
#' @return Returns a list of length two of class scan. The
#'   first element (clusters) is a list containing the
#'   significant, non-ovlappering clusters, and has the the
#'   following components: \item{locids}{The location ids of
#'   regions in a significant cluster.} \item{pop}{The total
#'   population in the cluser window.} \item{cases}{The
#'   observed number of cases in the cluster window.}
#'   \item{expected}{The expected number of cases in the
#'   cluster window.} \item{smr}{Standarized mortaility
#'   ratio (observed/expected) in the cluster window.}
#'   \item{rr}{Relative risk in the cluster window.}
#'   \item{loglikrat}{The loglikelihood ratio for the
#'   cluster window (i.e., the log of the test statistic).}
#'   \item{pvalue}{The pvalue of the test statistic
#'   associated with the cluster window.} The second element
#'   of the list is the centroid coordinates.  This is
#'   needed for plotting purposes.
#' @author Joshua French
#' @seealso \code{\link{print.smerc_cluster}},
#' \code{\link{summary.smerc_cluster}},
#' \code{\link{plot.smerc_cluster}},
#' \code{\link{scan.stat}}, \code{\link{scan.test}}
#' @export
#' @references Patil, G.P. & Taillie, C. Upper level set
#'   scan statistic for detecting arbitrarily shaped
#'   hotspots. Environmental and Ecological Statistics
#'   (2004) 11(2):183-197.
#'   <doi:10.1023/B:EEST.0000027208.48919.7e>
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = uls.test(coords = coords, cases = floor(nydf$cases),
#'                pop = nydf$pop, w = nyw,
#'                alpha = 0.05, longlat = TRUE,
#'                nsim = 9, ubpop = 0.5)
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
uls.test = function(coords, cases, pop, w,
                    ex = sum(cases) / sum(pop) * pop,
                    nsim = 499, alpha = 0.1,
                    ubpop = 0.5, longlat = FALSE,
                    cl = NULL, type = "poisson",
                    check.unique = FALSE) {
  # sanity checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha,
                      nsim + 1, ubpop, longlat, TRUE,
                      k = 1, w = w, type = type)

  coords = as.matrix(coords)
  zones = uls.zones(cases, pop, w, ubpop)

  # compute needed information
  ty = sum(cases)
  yin = zones.sum(zones, cases)

  # compute test statistics for observed data
  if (type == "poisson") {
    ein = zones.sum(zones, ex)
    tobs = stat.poisson(yin, ty - yin, ein, ty - ein)
  } else if (type == "binomial") {
    tpop = sum(pop)
    popin = zones.sum(zones, pop)
    tobs = stat.binom(yin, ty - yin, ty, popin, tpop - popin, tpop)
  }

  # compute test statistics for simulated data
  if (nsim > 0) {
    message("computing statistics for simulated data:")
    tsim = uls.sim(nsim = nsim, ty = ty, ex = ex, w = w,
                   pop = pop, ubpop = ubpop, cl = cl)
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
                longlat = longlat, method = "upper level set",
                rel_param = list(type = type,
                                 simdist = "multinomial",
                                 nsim = nsim,
                                 ubpop = ubpop),
                alpha = alpha,
                w = w, d = NULL)
}
