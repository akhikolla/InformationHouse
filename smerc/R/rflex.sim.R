#' Perform \code{rflex.test} on simualated data
#'
#' \code{rflex.sim} efficiently performs
#' \code{\link{rflex.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{rflex.test}} function, but is informative for
#' better understanding the implementation of the test.
#'
#' @param nsim A positive integer indicating the number of
#'   simulations to perform.
#' @param nn A matrix of the k nearest neighbors for the
#'   regions described by \code{w}.
#' @inheritParams rflex.test
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' data(nyw)
#' # determine knn
#' coords = with(nydf, cbind(longitude, latitude))
#' nn = knn(coords, longlat = TRUE, k = 50)
#' # determine expected number of cases in each region
#' cases = floor(nydf$cases)
#' pop = nydf$pop
#' ex = pop * sum(cases)/sum(pop)
#' tsim = rflex.sim(nsim = 5, nn = nn, w = nyw, ex = ex)
rflex.sim = function(nsim = 1, nn, w, ex, alpha1 = 0.2,
                     type = "poisson", pop = NULL, cl = NULL) {
  if (length(nsim) != 1 | !is.numeric(nsim) | nsim < 1) {
    stop("nsim must be a positive integer")
  }
  arg_check_rflex_zones(nn = nn, w = w, cases = ex, ex = ex,
                        alpha1 = alpha1, type = type,
                        pop = pop, loop = FALSE,
                        verbose = FALSE, 1)

  # determine total cases and population (if binomial)
  ty = sum(ex)
  if (type == "binomial") tpop = sum(pop)

  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    ysim = stats::rmultinom(1, size = ty, prob = ex)
    # determine rflex zones
    zones = rflex.zones(nn = nn, w = w,
                        cases = ysim, ex = ex,
                        alpha1 = alpha1, pop = pop,
                        cl = NULL, verbose = FALSE)
    # compute test statistics for each zone
    yin = zones.sum(zones, ysim)
    if (type == "poisson") {
      ein = zones.sum(zones, ex)
      tall = stat.poisson(yin, ty - yin, ein, ty - ein)
    } else if (type == "binomial") {
      popin = zones.sum(zones, pop)
      tall = stat.binom(yin, ty - yin, ty,
                        popin, tpop - popin, tpop)
    }
    max(tall)
  })
  unlist(tsim, use.names = FALSE)
}
