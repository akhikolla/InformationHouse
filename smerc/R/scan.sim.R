#' Perform \code{scan.test} on simulated data
#'
#' \code{scan.sim} efficiently performs
#' \code{\link{scan.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{scan.test}} function, but is informative for
#' better understanding the implementation of the test.
#'
#' @inheritParams flex.sim
#' @inheritParams scan.test
#' @param nn A list of nearest neighbors produced by \code{\link{nnpop}}.
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' d = sp::spDists(as.matrix(coords), longlat = TRUE)
#' nn = scan.nn(d, pop = nydf$pop, ubpop = 0.1)
#' cases = floor(nydf$cases)
#' ty = sum(cases)
#' ex = ty/sum(nydf$pop) * nydf$pop
#' yin = nn.cumsum(nn, cases)
#' ein = nn.cumsum(nn, ex)
#' tsim = scan.sim(nsim = 1, nn, ty, ex, ein = ein, eout = sum(ex) - ein)
scan.sim = function(nsim = 1, nn, ty, ex, type = "poisson",
                    ein = NULL, eout = NULL,
                    tpop = NULL, popin = NULL, popout = NULL,
                    cl = NULL,
                    simdist = "multinomial",
                    pop = NULL) {
  # match simdist with options
  simdist = match.arg(simdist, c("multinomial", "poisson", "binomial"))
  arg_check_sim(nsim = nsim, ty = ty, ex = ex, type = type,
                nn = nn, ein = ein, eout = eout, tpop = tpop,
                popin = popin, popout = popout, static = TRUE,
                simdist = simdist, pop = pop,
                w = diag(length(ex)))

  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    if (simdist == "multinomial") {
      ysim = stats::rmultinom(1, size = ty, prob = ex)
    } else if (simdist == "poisson") {
      ysim = stats::rpois(length(ex), lambda = ex)
      ty = sum(ysim)
      mult = ty / sum(ex)
      ein = ein * mult
      eout = eout * mult
    } else if (simdist == "binomial") {
      ysim = stats::rbinom(n = length(ex), size = pop,
                           prob = ex / pop)
      ty = sum(ysim)
    }
    # compute test statistics for each zone
    yin = nn.cumsum(nn, ysim)
    if (type == "poisson") {
      tall = stat.poisson(yin, ty - yin, ein, eout)
    } else if (type == "binomial") {
      tall = stat.binom(yin, ty - yin, ty, popin, popout, tpop)
    }
    max(tall)
  })
  unlist(tsim, use.names = FALSE)
}

#' Argument checking for *.sim functions
#'
#' Check the arguments of the \code{*.sim} functions.
#'
#' @param nsim Number of simulations
#' @param ty Total number of cases
#' @param ex Expected counts
#' @param type Type of statistic
#' @param nn List of nn (e.g., nnpop function)
#' @param zones List of zones (e.g., scan.zones)
#' @param ein List of expected in each zone
#' @param eout List of expected out of each zone
#' @param tpop Total population
#' @param popin Population in each zone
#' @param popout Population outside of each zone
#' @param w Spatial adjacency matrix
#' @param pop Vector of populations
#' @param ubpop Population upperbound
#' @param static Static zones. Logical. TRUE for scan.test.
#' FALSE for uls.test.
#' @param simdist Simulation distribution.
#' @return NULL
#' @noRd
arg_check_sim = function(nsim, ty, ex, type,
                         nn = NULL, zones = NULL,
                         ein = NULL, eout = NULL,
                         tpop = NULL, popin = NULL,
                         popout = NULL, w = NULL,
                         pop = NULL, ubpop = NULL,
                         static = FALSE,
                         simdist = "multinomial") {
  arg_check_nsim(nsim)
  arg_check_ty(ty)
  N = length(ex)
  arg_check_ex(ex, N)
  arg_check_type(type)
  if (!is.null(nn)) {
    if (!is.list(nn)) stop("nn must be a list")
    nz = sum(sapply(nn, length))
  }
  if (!is.null(zones)) {
    if (!is.list(zones)) stop("zones must be a list")
  }
  arg_check_w(w, N)
  if (!is.null(ubpop)) {
    arg_check_ubpop
  }
  if (!is.null(pop)) {
    arg_check_pop(pop, N)
  }
  arg_check_simdist(simdist)
  if (simdist == "binomial" & is.null(pop)) {
    stop("pop must be specified when simdist == 'binomial'")
  }
}

