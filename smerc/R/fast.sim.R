#' Perform \code{fast.test} on simulated data
#'
#' \code{fast.sim} efficiently performs
#' \code{\link{fast.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{fast.test}} function, but is informative for
#' better understanding the implementation of the test.
#'
#' @inheritParams scan.sim
#' @inheritParams fast.test
#' @inherit scan.sim return
#' @export
#'
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' cases = floor(nydf$cases)
#' pop = nydf$pop
#' ty = sum(cases)
#' ex = ty/sum(pop) * pop
#' tsim = fast.sim(1, ty, ex, pop = pop, ubpop = 0.5)
fast.sim = function(nsim = 1, ty, ex, pop, ubpop,
                    type = "poisson", cl = NULL) {
  tpop = sum(pop)
  arg_check_sim(nsim = nsim, ty = ty, ex = ex, type = type,
                tpop = tpop, ubpop = ubpop,
                w = diag(length(ex)))

  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    ysim = stats::rmultinom(1, size = ty, prob = ex)
    zones = fast.zones(ysim, pop = pop, ubpop = ubpop)
    # compute test statistics for each zone
    yin = cumsum(ysim[zones])
    if (type == "poisson") {
      ein = cumsum(ex[zones])
      tall = stat.poisson(yin, ty - yin, ein, ty - ein)
    } else if (type == "binomial") {
      popin = cumsum(pop[zones])
      tall = stat.binom(yin, ty - yin, ty,
                        popin, tpop - popin, tpop)
    }
    max(tall)
  })
  unlist(tsim, use.names = FALSE)
}
