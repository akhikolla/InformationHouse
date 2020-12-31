#' Restricted Flexibly-shaped Spatial Scan Test
#'
#' \code{rflex.test} performs the restricted flexibly shaped
#' spatial scan test of Tango and Takahashi (2012).
#'
#' The test is performed using the spatial scan test based
#' on the Poisson test statistic and a fixed number of
#' cases.  The first cluster is the most likely to be a
#' cluster.  If no significant clusters are found, then the
#' most likely cluster is returned (along with a warning).
#'
#' @inheritParams scan.test
#' @param w A binary spatial adjacency matrix for the
#'   regions.
#' @param k An integer indicating the maximum number of
#'   regions to inclue in a potential cluster.  Default is
#'   10
#' @param alpha1 The middle p-value threshold.
#'
#' @return Returns a list of length two of class scan. The
#'   first element (clusters) is a list containing the
#'   significant, non-ovlappering clusters, and has the the
#'   following components: \item{coords}{The centroid of the
#'   significant clusters.} \item{r}{The radius of the
#'   window of the clusters.} \item{pop}{The total
#'   population in the cluster window.} \item{cases}{The
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
#' @export
#' @seealso \code{\link{print.smerc_cluster}},
#' \code{\link{summary.smerc_cluster}},
#' \code{\link{plot.smerc_cluster}},
#' \code{\link{scan.stat}}, \code{\link{scan.test}}
#' @references Tango, T. and Takahashi, K. (2012), A
#'   flexible spatial scan statistic with a restricted
#'   likelihood ratio for detecting disease clusters.
#'   Statist. Med., 31: 4207-4218. <doi:10.1002/sim.5478>
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = rflex.test(coords = coords, cases = floor(nydf$cases),
#'                  w = nyw, k = 10,
#'                  pop = nydf$pop, nsim = 49,
#'                  alpha = 0.05, longlat = TRUE)
#'
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
rflex.test = function(coords, cases, pop, w, k = 50,
                     ex = sum(cases) / sum(pop) * pop,
                     type = "poisson",
                     nsim = 499,
                     alpha = 0.1,
                     longlat = FALSE,
                     alpha1 = 0.2,
                     cl = NULL) {
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha,
                      nsim + 1, 0.5, longlat, FALSE, k = k,
                      w = w, type = type)
  coords = as.matrix(coords)

  # compute k nearest neighbors
  nn = knn(coords = coords, longlat = longlat, k = k)

  # determine zones for observed data
  zones = rflex.zones(nn = nn, w = w,
                      cases = cases, ex = ex, alpha1 = alpha1,
                      cl = cl, verbose = FALSE)
  # compute needed information
  ty = sum(cases)
  yin = zones.sum(zones, cases)
  if (type == "binomial") tpop = sum(pop)

  # compute test statistics for observed data
  if (type == "poisson") {
    ein = zones.sum(zones, ex)
    tobs = stat.poisson(yin, ty - yin, ein, ty - ein)
  } else if (type == "binomial") {
    popin = zones.sum(zones, pop)
    tobs = stat.binom(yin, ty - yin, ty,
                      popin, tpop - popin, tpop)
  }

  # compute test statistics for simulated data
  if (nsim > 1) {
    tsim = rflex.sim(nsim = nsim, nn = nn, w = w, ex = ex,
                     alpha1 = alpha1, type = type,
                     pop = pop, cl = cl)
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
                longlat = longlat, method = "restricted flexible",
                rel_param = list(type = type,
                                 simdist = "multinomial",
                                 nsim = nsim,
                                 alpha1 = alpha1,
                                 k = k),
                alpha = alpha,
                w = w, d = NULL)
}
