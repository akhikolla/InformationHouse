#' Fast Subset Scan Test
#'
#' \code{fast.test} performs the fast subset scan test of
#' Neill (2012).
#'
#' The test is performed using the spatial scan test based
#' on the Poisson test statistic and a fixed number of
#' cases.  The windows are based on the Upper Level Sets
#' proposed by Patil and Taillie (2004).  The clusters
#' returned are non-overlapping, ordered from most
#' significant to least significant.  The first cluster is
#' the most likely to be a cluster.  If no significant
#' clusters are found, then the most likely cluster is
#' returned (along with a warning).
#'
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
#' @references Neill, D. B. (2012), Fast subset scan for
#'   spatial pattern detection. Journal of the Royal
#'   Statistical Society: Series B (Statistical
#'   Methodology), 74: 337-360.
#'   <doi:10.1111/j.1467-9868.2011.01014.x>
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = fast.test(coords = coords, cases = floor(nydf$cases),
#'                pop = nydf$pop,
#'                alpha = 0.05, longlat = TRUE,
#'                nsim = 49, ubpop = 0.5)
fast.test = function(coords, cases, pop,
                    ex = sum(cases) / sum(pop) * pop,
                    nsim = 499, alpha = 0.1,
                    ubpop = 0.5, longlat = FALSE,
                    cl = NULL, type = "poisson") {
  # sanity checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha,
                      nsim + 1, ubpop, longlat, TRUE,
                      k = 1, w = diag(nrow(coords)),
                      type = type)

  coords = as.matrix(coords)
  zones = fast.zones(cases, pop, ubpop)

  # compute needed information
  ty = sum(cases)
  yin = cumsum(cases[zones])

  # compute test statistics for observed data
  if (type == "poisson") {
    ein = cumsum(ex[zones])
    tobs = stat.poisson(yin, ty - yin, ein, ty - ein)
  } else if (type == "binomial") {
    tpop = sum(pop)
    popin = cumsum(pop[zones])
    tobs = stat.binom(yin, ty - yin, ty, popin, tpop - popin, tpop)
  }

  # compute test statistics for simulated data
  if (nsim > 0) {
    message("computing statistics for simulated data:")
    tsim = fast.sim(nsim = nsim, ty = ty, ex = ex,
                   pop = pop, ubpop = ubpop, cl = cl)
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

  # reformat zones
  zones = lapply(seq_along(zones), function(i) zones[seq_len(i)])

  # significant, ordered, non-overlapping clusters and
  # information
  pruned = sig_noc(tobs = tobs,
                   zones = zones,
                   pvalue = pvalue, alpha = alpha,
                   order_by = "tobs")

  smerc_cluster(tobs = pruned$tobs, zones = pruned$zones,
                pvalue = pruned$pvalue, coords = coords,
                cases = cases, pop = pop, ex = ex,
                longlat = longlat, method = "fast",
                rel_param = list(type = type,
                                 simdist = "multinomial",
                                 nsim = nsim,
                                 ubpop = ubpop),
                alpha = alpha,
                w = NULL, d = NULL)
}
