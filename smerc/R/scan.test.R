#' Spatial Scan Test
#'
#' \code{scan.test} performs the original spatial scan test
#' of Kulldorf (1997) based on a fixed number of cases.
#' Candidate zones are circular and extend from the observed
#' region centroids.  The clusters returned are
#' non-overlapping, ordered from most significant to least
#' significant.  The first cluster is the most  likely to be
#' a cluster.  If no significant clusters are found, then
#' the most likely cluster is returned (along with a
#' warning).
#'
#' @param coords An \eqn{n \times 2} matrix of centroid
#'   coordinates for the regions.
#' @param cases The number of cases observed in each region.
#' @param pop The population size associated with each
#'   region.
#' @param ex The expected number of cases for each region.
#'   The default is calculated under the constant risk
#'   hypothesis.
#' @param nsim The number of simulations from which to
#'   compute the p-value.
#' @param ubpop The upperbound of the proportion of the
#'   total population to consider for a cluster.
#' @param alpha The significance level to determine whether
#'   a cluster is signficant.  Default is 0.10.
#' @param longlat The default is \code{FALSE}, which
#'   specifies that Euclidean distance should be used. If
#'   \code{longlat} is \code{TRUE}, then the great circle
#'   distance is used to calculate the intercentroid
#'   distance.
#' @param type The type of scan statistic to compute. The
#'   default is \code{"poisson"}. The other choice
#'   is \code{"binomial"}.
#' @param min.cases The minimum number of cases required for
#'   a cluster.  The default is 2.
#' @param simdist Character string indicating the simulation
#' distribution. The default is \code{"multinomial"}, which
#' conditions on the total number of cases observed. The
#' other options are \code{"poisson"} and \code{"binomial"}
#' @inheritParams pbapply::pblapply
#'
#' @return Returns a \code{smerc_cluster} object.
#' @seealso \code{\link{print.smerc_cluster}},
#' \code{\link{summary.smerc_cluster}},
#' \code{\link{plot.smerc_cluster}},
#' \code{\link{scan.stat}}
#' @author Joshua French
#' @export
#' @references Kulldorff, M. (1997) A spatial scan
#'   statistic. Communications in Statistics - Theory and
#'   Methods, 26(6): 1481-1496,
#'   <doi:10.1080/03610929708831995>
#'
#' Waller, L.A. and Gotway, C.A. (2005). Applied Spatial
#' Statistics for Public Health Data. Hoboken, NJ: Wiley.
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 0,
#'                 alpha = 1, longlat = TRUE)
#' ## plot output for new york state
#' # specify desired argument values
#' mapargs = list(database = "county", region = "new york",
#' xlim = range(out$coords[,1]), ylim = range(out$coords[,2]))
#' # needed for "state" database (unless you execute library(maps))
#' data(countyMapEnv, package = "maps")
#' plot(out, usemap = TRUE, mapargs = mapargs)
#'
#' # a second example to match the results of Waller and Gotway (2005)
#' # in chapter 7 of their book (pp. 220-221).
#' # Note that the 'longitude' and 'latitude' used by them has
#' # been switched.  When giving their input to SatScan, the coords
#' # were given in the order 'longitude' and 'latitude'.
#' # However, the SatScan program takes coordinates in the order
#' # 'latitude' and 'longitude', so the results are slightly different
#' # from the example above.
#' coords = with(nydf, cbind(y, x))
#' out2 = scan.test(coords = coords, cases = floor(nydf$cases),
#'                   pop = nydf$pop, nsim = 0,
#'                   alpha = 1, longlat = TRUE)
#' # the cases observed for the clusters in Waller and Gotway: 117, 47, 44
#' # the second set of results match
#' sget(out2$clusters, name = "cases")[1:3]
scan.test = function(coords, cases, pop,
                     ex = sum(cases) / sum(pop) * pop,
                     nsim = 499, alpha = 0.1,
                     ubpop = 0.5, longlat = FALSE, cl = NULL,
                     type = "poisson",
                     min.cases = 2,
                     simdist = "multinomial") {
  # argument checking
  type = match.arg(type, c("poisson", "binomial"))
  simdist = match.arg(simdist, c("multinomial", "poisson", "binomial"))
  arg_check_scan_test(coords = coords, cases = cases,
                      pop = pop, ex = ex, nsim = nsim,
                      alpha = alpha, ubpop = ubpop,
                      longlat = longlat,
                      k = 1, w = diag(nrow(coords)),
                      type = type, simdist = simdist,
                      min.cases = min.cases)

  # convert to proper format
  coords = as.matrix(coords)
  N = nrow(coords)
  # compute inter-centroid distances
  d = sp::spDists(coords, longlat = longlat)

  # for each region, determine sorted nearest neighbors
  # subject to population constraint
  nn = scan.nn(d, pop, ubpop)

  # determine total number of cases in each successive
  # window, total number of cases
  yin = nn.cumsum(nn, cases)
  ty = sum(cases) # sum of all cases

  # compute test statistics for observed data
  if (type == "poisson") {
    ein = nn.cumsum(nn, ex)
    eout = sum(ex) - ein
    popin = NULL
    popout = NULL
    tpop = NULL
    tobs = stat.poisson(yin, ty - yin, ein, eout)
  } else if (type == "binomial") {
    ein = NULL
    eout = NULL
    tpop = sum(pop)
    popin = nn.cumsum(nn, pop)
    popout = tpop - popin
    tobs = stat.binom(yin, ty - yin, ty, popin, popout, tpop)
  }

  # determine distinct zones
  wdup = nndup(nn, N)

  # remove zones with a test statistic of 0 or don't have
  # min number of cases or are duplicted
  w0 = which(tobs == 0 | yin < min.cases | wdup)

  # determine zones
  zones = nn2zones(nn)

  # remove zones with a test statistic of 0
  zones = zones[-w0]
  tobs = tobs[-w0]

  # compute test statistics for simulated data
  if (nsim > 0) {
    message("computing statistics for simulated data:")
    tsim = scan.sim(nsim = nsim, nn = nn, ty = ty,
                    ex = ex, type = type, ein = ein,
                    eout = eout, popin = popin,
                    popout = popout, tpop = tpop, cl = cl,
                    simdist = simdist, pop = pop)
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
                longlat = longlat, method = "circular scan",
                rel_param = list(type = type,
                                 simdist = simdist,
                                 nsim = nsim,
                                 ubpop = ubpop,
                                 min.cases = min.cases),
                alpha = alpha,
                w = NULL, d = d)
}

#' Argument checking for scan tests
#'
#' @param coords A matrix of coordinates
#' @param cases A vector of numeric cases
#' @param pop A vector of population values
#' @param ex A vector of expected counts
#' @param nsim A non-negative integer
#' @param alpha A value greater than 0
#' @param nreport Not used
#' @param ubpop A value between 0 and 1
#' @param longlat A logical. TRUE is great circle distance.
#' @param parallel Not used.
#' @param k Number of nearest neighbors. Not always needed.
#' @param w A spatial proximity matrix
#' @param type Statistic type
#' @param simdist Distribution of simulation
#' @param min.cases Minimum number of cases. Only for scan.test.
#' @return NULL
#' @noRd
arg_check_scan_test =
  function(coords, cases, pop, ex, nsim, alpha,
           nreport = NULL,
           ubpop, longlat, parallel = NULL, k, w, type = NULL,
           simdist = NULL, min.cases = NULL) {
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_ex(ex, N)
  arg_check_nsim(nsim)
  arg_check_alpha(alpha)
  # nreport no check, deprecated
  arg_check_ubpop(ubpop)
  arg_check_longlat(longlat)
  # parallel no check, deprecated
  arg_check_k(k, N)
  arg_check_w(w, N)
  if (!is.null(type)) {
    arg_check_type(type)
  }
  if (!is.null(simdist)) {
    arg_check_simdist(simdist)
  }
  arg_check_simdist(simdist)
  if (!is.null(min.cases)) {
    arg_check_min_cases(min.cases)
  }
}

