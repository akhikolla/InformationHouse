#' Perform \code{uls.test} on simulated data
#'
#' \code{uls.sim} efficiently performs
#' \code{\link{uls.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{uls.test}} function, but is informative for
#' better understanding the implementation of the test.
#'
#' @inheritParams scan.sim
#' @inheritParams uls.test
#' @inherit scan.sim return
#' @export
#'
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' cases = floor(nydf$cases)
#' pop = nydf$pop
#' ty = sum(cases)
#' ex = ty/sum(pop) * pop
#' tsim = uls.sim(1, ty, ex, nyw, pop = pop, ubpop = 0.5)
uls.sim = function(nsim = 1, ty, ex, w, pop, ubpop,
                   type = "poisson", check.unique = FALSE,
                   cl = NULL) {
  tpop = sum(pop)
  arg_check_sim(nsim = nsim, ty = ty, ex = ex, type = type,
                tpop = tpop, w = w, ubpop = ubpop,
                static = FALSE)

  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    ysim = c(stats::rmultinom(1, size = ty, prob = ex))
    zones = uls.zones(cases = ysim, pop = pop, w = w,
                      ubpop = ubpop,
                      check.unique = check.unique)
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
