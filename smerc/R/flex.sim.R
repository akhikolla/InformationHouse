#' Perform \code{flex.test} on simualated data
#'
#' \code{flex.sim} efficiently performs
#' \code{\link{flex.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{flex.test}} function, but is informative for
#' better understanding the implementation of the test.
#'
#' @param nsim A positive integer indicating the number of
#'   simulations to perform.
#' @param zones A list of zones to compute the test statistic
#' over for each simulated data set.
#' @inheritParams flex.test
#' @inheritParams scan.stat
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' zones = flex.zones(coords, w = nyw, k = 3, longlat = TRUE)
#' cases = floor(nydf$cases)
#' ty = sum(cases)
#' ex = ty/sum(nydf$pop) * nydf$pop
#' ein = zones.sum(zones, ex)
#' tsim = flex.sim(nsim = 2, zones, ty, ex, ein = ein, eout = ty - ein)
flex.sim = function(nsim = 1, zones, ty, ex, type = "poisson",
                    ein = NULL, eout = NULL,
                    tpop = NULL, popin = NULL, popout = NULL,
                    cl = NULL) {

  arg_check_flex_sim(nsim, zones, ty, ex, ein, eout,
                     tpop, popin, popout, type)

  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    ysim = stats::rmultinom(1, size = ty, prob = ex)
    # compute test statistics for each zone
    yin = zones.sum(zones, ysim)
    if (type == "poisson") {
      tall = stat.poisson(yin, ty - yin, ein, eout)
    } else if (type == "binomial") {
      tall = stat.binom(yin, ty - yin, ty,
                        popin, popout, tpop)
    }
    max(tall)
  })
  unlist(tsim, use.names = FALSE)
}

arg_check_flex_sim = function(nsim, zones, ty, ex, ein, eout,
                              tpop, popin, popout, type) {
  if (length(nsim) != 1 | !is.numeric(nsim) | nsim < 1) {
    stop("nsim must be a positive integer")
  }
  if (!is.list(zones)) stop("zones must be a list")
  if (length(ty) != 1) stop("ty must be a single number")
  if (!is.numeric(ex)) stop("ex must be a numeric vector")
  if (length(type) != 1 | !is.element(type, c("poisson", "binomial"))) {
    stop("type must be 'poisson' or 'binomial'")
  }
  if (type == "poisson") {
    if (is.null(ein) | is.null(eout)) {
      stop("ein and eout must be provided when type is 'poisson'")
    }
    if (length(zones) != length(ein)) {
      stop("length(zones) != length(ein)")
    }
    if (length(zones) != length(eout)) {
      stop("length(zones) != length(eout)")
    }
  }
  if (type == "binomial") {
    if (is.null(popin) | is.null(popout) | is.null(tpop)) {
      stop("popin, popout, and tpop must be provided when type is 'binomial'")
    }
    if (length(zones) != length(popin)) {
      stop("length(zones) != length(popin)")
    }
    if (length(zones) != length(popout)) {
      stop("length(zones) != length(popout)")
    }
    if (length(tpop) != 1) {
      stop("tpop must be a single number")
    }
  }
}
