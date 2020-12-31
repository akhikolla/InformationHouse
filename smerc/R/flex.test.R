#' Flexibly-shaped Spatial Scan Test
#'
#' \code{flex.test} performs the flexibly-shaped scan test
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
#' out = flex.test(coords = coords, cases = floor(nydf$cases),
#'                 w = nyw, k = 3,
#'                 pop = nydf$pop, nsim = 49,
#'                 alpha = 0.12, longlat = TRUE)
#'
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
flex.test = function(coords, cases, pop, w, k = 10,
                     ex = sum(cases) / sum(pop) * pop,
                     type = "poisson", nsim = 499,
                     alpha = 0.1, longlat = FALSE,
                     cl = NULL,
                     lonlat = longlat, ...) {
  if (!identical(lonlat, longlat)) {
    longlat = lonlat
    warning("lonlat is deprecated. Please use longlat.")
  }
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha,
                      nsim + 1, 0.5, longlat, FALSE, k = k,
                      w = w, type = type)

  coords = as.matrix(coords)

  zones = flex.zones(coords, w, k, longlat)

  # compute needed information
  ty = sum(cases)
  yin = zones.sum(zones, cases)

  # compute test statistics for observed data
  if (type == "poisson") {
    ein = zones.sum(zones, ex)
    eout = ty - ein
    popin = NULL
    popout = NULL
    tpop = NULL
    tobs = stat.poisson(yin, ty - yin, ein, eout)
  } else if (type == "binomial") {
    ein = NULL
    eout = NULL
    tpop = sum(pop)
    popin = zones.sum(zones, pop)
    popout = tpop - popin
    tobs = stat.binom(yin, ty - yin, ty, popin, popout, tpop)
  }

  # compute test statistics for simulated data
  if (nsim > 1) {
    message("computing statistics for simulated data:")
    tsim = flex.sim(nsim = nsim, zones = zones, ty = ty,
                    ex = ex,
                    type = type, ein = ein, eout = eout,
                    popin = popin, popout = popout, tpop = tpop,
                    cl = cl)
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
                longlat = longlat, method = "flexible",
                rel_param = list(type = type,
                                 simdist = "multinomial",
                                 nsim = nsim,
                                 k = k),
                alpha = alpha,
                w = w, d = NULL)
}
