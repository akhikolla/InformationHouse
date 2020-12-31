#' Perform \code{cepp.test} on simulated data
#'
#' \code{cepp.sim} efficiently performs
#' \code{\link{cepp.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{cepp.test}} function, but is informative for
#' better understanding the implementation of the test.
#'
#' @inheritParams flex.sim
#' @inheritParams cepp.test
#' @param nn A list of nearest neighbors produced by \code{\link{casewin}}.
#' @param wts A list that has the weights associated with each
#' region of each element of \code{nn}.
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' d = sp::spDists(as.matrix(coords), longlat = TRUE)
#' nn = casewin(d, cases = nydf$pop, cstar = 15000)
#' cases = floor(nydf$cases)
#' ty = sum(cases)
#' ex = ty/sum(nydf$pop) * nydf$pop
#' # find smallest windows with at least n* pop
#' nstar = 1000
#' nn = casewin(d, cases = nydf$pop, cstar = nstar)
#' # determine ts
#' wts = cepp.weights(nn, nydf$pop, nstar)
#' tsim = cepp.sim(1, nn = nn, ty = ty, ex = ex, wts = wts)
cepp.sim = function(nsim = 1, nn, ty, ex, wts, simdist = "multinomial") {
  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(idx) {
    # simulate new data
    if (simdist == "multinomial") {
      ysim = stats::rmultinom(1, size = ty, prob = ex)
    } else {
      ysim = stats::rpois(length(ex), lambda = ex)
    }
    cstar = sapply(seq_along(nn), function(i) {
      sum(ysim[nn[[i]]] * wts[[i]])
    })
    max(cstar)
  })
  unlist(tsim, use.names = FALSE)
}
