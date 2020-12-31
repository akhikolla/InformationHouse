#' Maxima Likelihood First Scan Test
#'
#' \code{mlf.test} implements the Maxima Likelihood First
#' scan test of Yao et al. (2011), which is actually a
#' special case of the Dynamic Minimum Spanning Tree of
#' Assuncao et al. (2006).  Find the single region that
#' maximizes the likelihood ratio test statistic.  Starting
#' with this single region as a current zone, new candidate
#' zones are constructed by combining the current zone with
#' the connected region that maximizes the likelihood ratio
#' test statisic.  This procedure is repeated until the
#' population and/or distance upper bound is reached.
#'
#' Only a single candidate zone is ever returned because the
#' algorithm only constructs a single sequence of starting
#' zones, and overlapping zones are not returned.  Only the
#' zone that maximizes the likelihood ratio test statistic
#' is returned.
#'
#' @inheritParams dmst.test
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
#'   associated with the cluster window.} \item{w}{The
#'   adjacency matrix of the cluster.} \item{r}{The maximum
#'   radius of the cluster (in terms of intercentroid
#'   distance from the starting region).} The second element
#'   of the list is the centroid coordinates.  This is
#'   needed for plotting purposes.
#' @author Joshua French
#' @seealso \code{\link{print.smerc_cluster}},
#' \code{\link{summary.smerc_cluster}},
#' \code{\link{plot.smerc_cluster}},
#' \code{\link{scan.stat}}, \code{\link{scan.test}}
#' @export
#' @references Yao, Z., Tang, J., & Zhan, F. B. (2011).
#'   Detection of arbitrarily-shaped clusters using a
#'   neighbor-expanding approach: A case study on murine
#'   typhus in South Texas. International journal of health
#'   geographics, 10(1), 1.
#'
#'   Assuncao, R.M., Costa, M.A., Tavares, A. and Neto,
#'   S.J.F. (2006). Fast detection of arbitrarily shaped
#'   disease clusters, Statistics in Medicine, 25, 723-742.
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = mlf.test(coords = coords, cases = floor(nydf$cases),
#'                   pop = nydf$pop, w = nyw,
#'                   alpha = 0.12, longlat = TRUE,
#'                   nsim = 10, ubpop = 0.1, ubd = 0.5)
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
mlf.test = function(coords, cases, pop, w,
                    ex = sum(cases) / sum(pop) * pop,
                    nsim = 499, alpha = 0.1,
                    ubpop = 0.5, ubd = 0.5,
                    longlat = FALSE, cl = NULL) {
  # sanity checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha,
                      nsim + 1, ubpop, longlat, FALSE,
                      k = 1, w = w)

  coords = as.matrix(coords)
  ty = sum(cases) # sum of all cases

  # calculate test statistics for each individual region
  # yin, ein, eout, ty
  eout = ty - ex
  tobs = scan.stat(cases, ex, eout, ty)

  # determine starting region for maxima likelihood first
  # algorithm
  start = which.max(tobs)

  # intercentroid distances
  d = sp::spDists(coords, longlat = longlat)

  # upperbound for population in zone
  max_pop = ubpop * sum(pop)
  # upperbound for distance between centroids in zone
  max_dist = ubd * max(d)

  # find neighbors for starting region
  all_neighbors = lapply(seq_along(cases),
                         function(i) which(d[i, ] <= max_dist))

  # return sequence of candidate zones (or a subset depending on type)
  max_zone = mst.seq(start, all_neighbors[[start]], cases,
                     pop, w, ex, ty, max_pop, "pruned")

  # determine which call for simulations
  fcall = pbapply::pblapply
  # setup list for call
  fcall_list = list(X = as.list(seq_len(nsim)), FUN = function(i) {
    # simulate new data set
    ysim = stats::rmultinom(1, size = ty, prob = ex)

    sim_tstat = scan.stat(ysim, ex, eout, ty)
    # determine starting region for maxima likelihood first algorithm
    sim_start = which.max(sim_tstat)

    # find max statistic for best candidate zone
    mst.seq(sim_start, all_neighbors[[sim_start]], cases,
            pop, w, ex, ty, max_pop, type = "maxonly")
  }, cl = cl)

  # get max statistics for simulated data sets
  tsim = unlist(do.call(fcall, fcall_list), use.names = FALSE)

  # p-values associated with these max statistics for each centroid
  pvalue = (sum(tsim >= max_zone$loglikrat) + 1) / (nsim + 1)

  # significant, ordered, non-overlapping clusters and
  # information
  pruned = sig_noc(tobs = max_zone$loglikrat,
                   zones = list(max_zone$locids),
                   pvalue = pvalue, alpha = alpha,
                   order_by = "tobs")

  smerc_cluster(tobs = pruned$tobs, zones = pruned$zones,
                pvalue = pruned$pvalue, coords = coords,
                cases = cases, pop = pop, ex = ex,
                longlat = longlat, method = "maxima likelihood first",
                rel_param = list(type = "poisson",
                                 simdist = "multinomial",
                                 nsim = nsim,
                                 ubpop = ubpop,
                                 ubd = ubd),
                alpha = alpha,
                w = w, d = NULL)
}
