#' Flexibly-shaped Spatial Scan Test
#'
#' \code{flex_test} performs the flexibly-shaped scan test
#' of Tango and Takahashi (2005).
#'
#' The test is performed using the spatial scan test based
#' on the Poisson test statistic and a fixed number of
#' cases.  The first cluster is the most likely to be a
#' cluster.  If no significant clusters are found, then the
#' most likely cluster is returned (along with a warning).
#'
#' @inheritParams rflex.test
#' @param lonlat Deprecated in favor of \code{longlat}.
#' @param ... Not used.
#'
#' @return Returns a list of length two of class scan. The
#'   first element (clusters) is a list containing the
#'   significant, non-ovlappering clusters, and has the the
#'   following components:
#' @author Joshua French
#' @export
#' @seealso \code{\link{print.smerc_cluster}},
#' \code{\link{summary.smerc_cluster}},
#' \code{\link{plot.smerc_cluster}},
#' \code{\link{scan.stat}}, \code{\link{scan.test}}
#' @references Tango, T., & Takahashi, K. (2005). A flexibly
#'   shaped spatial scan statistic for detecting clusters.
#'   International journal of health geographics, 4(1), 11.
#'   Kulldorff, M. (1997) A spatial scan statistic.
#'   Communications in Statistics -- Theory and Methods 26,
#'   1481-1496.
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = flex_test(coords = coords, cases = floor(nydf$cases),
#'                 w = nyw, k = 3,
#'                 pop = nydf$pop, nsim = 49,
#'                 alpha = 0.12, longlat = TRUE)
#'
#' data(nypoly)
#' library(sp)
#' # plot(nypoly, col = color.clusters(out))
flex_test = function(coords, cases, pop, w, k = 10,
                     ex = sum(cases) / sum(pop) * pop,
                     type = "poisson", nsim = 499,
                     alpha = 0.1, longlat = FALSE,
                     cl = NULL,
                     lonlat = longlat, ...) {
  if (!identical(lonlat, longlat)) {
    longlat = lonlat
    warning("lonlat is deprecated. Please use longlat.")
  }
  # arg_check_scan_test(coords, cases, pop, ex, nsim, alpha,
  #                     nsim + 1, 0.5, longlat, FALSE, k = k,
  #                     w = w, type = type)

  coords = as.matrix(coords)
  # get list of nearest neighbors
  nn = knn(coords, longlat = longlat, k = k)
  # get vector of log primes
  lprimes = log(randtoolbox::get.primes(nrow(w)))
  # convert type to integer
  itype = ifelse(type == "poisson", 0, 1)
  # run flex_test via cpp
  pruned = flex_test_cpp(nn = nn, cases = cases, pop = pop, w = w, k = k,
                         ex = ex, type = itype, nsim = nsim, alpha = alpha,
                         lprimes = lprimes, verbose = TRUE)
  pruned = sig_noc(tobs = pruned$tobs, zones = pruned$zones,
                   pvalue = pruned$pvalues, alpha = alpha,
                   order_by = "tobs")
  # pruned2 = NULL
  # pruned2$zones = logical2zones(pruned$zones, nn)
  # # convert tobs and pvalues from nested lists of vectors to vectors
  # pruned2$tobs = unlist(pruned$tobs)
  # pruned2$pvalue = unlist(pruned$pvalue)
  # pruned3 = sig_noc(tobs = pruned2$tobs, zones = pruned2$zones,
  #                   pvalue = pruned2$pvalue, alpha = 1,
  #                   order_by = "tobs")
  # return(pruned3)
  # convert z from nested list of logicals to list of zones
  # pruned$zones = logical2zones(pruned$zones, nn)
  # # convert tobs and pvalues from nested lists of vectors to vectors
  # pruned$tobs = unlist(pruned$tobs)
  # pruned$pvalue = unlist(pruned$pvalue)
  # pruned = sig_noc(tobs = pruned$tobs, zones = pruned$zones,
  #                  pvalue = pruned$pvalue, alpha = alpha,
  #                  order_by = "tobs")
  smerc_cluster(tobs = pruned$tobs, zones = pruned$zones,
                pvalue = pruned$pvalue, coords = coords,
                cases = cases, pop = pop, ex = ex,
                longlat = longlat, method = "flexible",
                rel_param = list(type = type,
                                 simdist = "multinomial",
                                 nsim = nsim,
                                 k = k,
                                 nn = nn),
                alpha = alpha,
                w = w, d = NULL)
}
