#' Perform \code{dc.test} on simulated data
#'
#' \code{dc.sim} efficiently performs \code{\link{dc.test}}
#' on a simulated data set.  The function is meant to be
#' used internally by the \code{\link{dc.test}} function,
#' but is informative for better understanding the
#' implementation of the test.
#' @param nn A list of distance-based nearest neighbors,
#'   preferably from the \code{\link{nndist}} function.
#' @inheritParams scan.sim
#' @inheritParams dc.test
#' @inheritParams mst.all
#' @inherit scan.sim return
#' @export
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' cases = floor(nydf$cases)
#' pop = nydf$pop
#' ty = sum(cases)
#' ex = ty/sum(pop) * pop
#' d = sp::spDists(coords, longlat = TRUE)
#' nn = nndist(d, ubd = 0.05)
#' max_pop = sum(pop) * 0.25
#' tsim = dc.sim(1, nn, ty, ex, nyw, pop = pop,
#'               max_pop = max_pop)
dc.sim = function(nsim = 1, nn, ty, ex, w, pop, max_pop,
                  cl = NULL) {
  arg_check_sim(nsim = nsim, ty = ty, ex = ex,
                type = "poisson", tpop = NULL, w = w,
                ubpop = 0.1, static = FALSE)

  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    ysim = stats::rmultinom(1, size = ty, prob = ex)
    tall = mst.all(nn, cases = ysim, pop = pop, w = w,
            ex = ex, ty = ty,
            max_pop = max_pop, type = "maxonly",
            early = TRUE, nlinks = "two", progress = FALSE)
    max(tall)
  })
  unlist(tsim, use.names = FALSE)
}
