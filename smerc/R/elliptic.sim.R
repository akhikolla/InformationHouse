#' Perform \code{elliptic.test} on simulated data
#'
#' \code{elliptic.sim} efficiently performs
#' \code{\link{elliptic.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{elliptic.test}} function, but is informative
#' for better understanding the implementation of the test.
#'
#' @inheritParams scan.sim
#' @param nn A list of nearest neighbors produced by
#'   \code{\link{elliptic.nn}}.
#' @inheritParams elliptic.test
#' @param shape_all A vector of the shapes associated with
#'   all of the possible zones constructed from \code{nn}.
#'   This can be obtained from \code{\link{elliptic.nn}}.
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' pop = nydf$pop
#' enn = elliptic.nn(coords, pop, ubpop = 0.1,
#'                   shape = c(1, 1.5), nangle = c(1, 4))
#' cases = floor(nydf$cases)
#' ty = sum(cases)
#' ex = ty/sum(pop) * pop
#' yin = nn.cumsum(enn$nn, cases)
#' ein = nn.cumsum(enn$nn, ex)
#' tsim = elliptic.sim(nsim = 2, nn = enn$nn, ty = ty, ex = ex,
#'                     a = 0.5, shape_all = enn$shape_all,
#'                     ein = ein, eout = ty - ein)
elliptic.sim = function(nsim = 1, nn, ty, ex, a, shape_all,
                        ein, eout, cl = NULL) {
  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    ysim = stats::rmultinom(1, size = ty, prob = ex)
    # compute test statistics for each zone
    yin = nn.cumsum(nn, ysim)
    max(stat.poisson(yin, ty - yin, ein, eout, a, shape_all))
  })
  unlist(tsim, use.names = FALSE)
}

arg_check_elliptic_sim = function(nsim, nn, ty, ex, ein, eout, shape_all) {
  if (length(nsim) != 1 | !is.numeric(nsim) | nsim < 1) {
    stop("nsim must be a positive integer")
  }
  if (!is.list(nn)) stop("nn must be a list")
  if (length(ty) != 1) stop("ty must be a single number")
  if (!is.numeric(ex)) stop("ex must be a numeric vector")
  nstat = sum(sapply(nn, length))
  if (length(shape_all) != nstat) {
    stop("length(shapenn) is not compatible with the dimensionality of nn")
  }
  if (length(ein) != nstat) {
    stop("the length of ein is not compatible with the dimensionality of nn")
  }
  if (length(eout) != nstat) {
    stop("the length of eout is not compatible with the dimensionality of nn")
  }
}
