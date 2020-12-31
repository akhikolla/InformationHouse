#' Tango's clustering detection test
#'
#' \code{tango.test} performs a test for clustering proposed
#' by Tango (1995).  The test uses Tango's chi-square
#' approximation for significance testing by default, but
#' also uses Monto Carlo simulation when \code{nsim > 0}.
#'
#' The \code{\link{dweights}} function can be used to
#' construct a weights matrix \code{w} using the method of
#' Tango (1995), Rogerson (1999), or a basic style.
#'
#' @inheritParams tango.stat
#' @param nsim The number of simulations for which to
#'   perform a Monto Carlo test of significance.  Counts are
#'   simulated according to a multinomial distribution with
#'   \code{sum(cases)} total cases and class probabilities
#'   \code{pop/sum(pop)}. \code{sum(cases)} .
#'
#' @seealso \code{\link{dweights}}
#' @return Returns a list of class \code{tango} with
#'   elements: \item{tstat}{Tango's index}
#'   \item{tstat.chisq}{The approximately chi-squared
#'   statistic proposed by Tango that is derived from
#'   \code{tstat}} \item{dfc}{The degrees of freedom of
#'   \code{tstat.chisq}} \item{pvalue.chisq}{The p-value
#'   associated with \code{tstat.chisq}}
#'   \item{tstat.sim}{The vector of test statistics from the
#'   simulated data if \code{nsim > 0}}
#'   \item{pvalue.sim}{The p-value associated with the Monte
#'   Carlo test of significance when \code{nsim > 0}}
#'   Additionally, the goodness-of-fit \code{gof} and
#'   spatial autocorrelation \code{sa} components of the
#'   Tango's index are provided (and for the simulated data
#'   sets also, if appropriate).
#' @references Tango, T.  (1995) A class of tests for
#' detecting "general" and "focused" clustering of rare
#' diseases.  Statistics in Medicine.  14, 2323-2334.
#'
#' Rogerson, P. (1999) The Detection of Clusters Using A
#' Spatial Version of the Chi-Square Goodness-of-fit Test.
#' Geographical Analysis. 31, 130-147
#'
#' Tango, T.  (2010) Statistical Methods for Disease
#' Clustering. Springer.
#'
#' Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial
#' Statistics for Public Health Data.  Hoboken, NJ: Wiley.
#' @author Joshua French
#' @export
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("x", "y")])
#' w = dweights(coords, kappa = 1)
#' results = tango.test(nydf$cases, nydf$pop, w, nsim = 49)
tango.test = function(cases, pop, w, nsim = 0) {
  arg_check_tango_test(cases, pop, w, nsim)
  N = length(cases)
  yplus = sum(cases)
  r = cases / yplus
  p = pop / sum(pop)
  ee = r - p
  gof = sum(ee ^ 2)
  sa = (crossprod(ee, w - diag(N)) %*% ee)[1, 1]
  tstat = gof + sa

  # compute standardized tstat
  vp = diag(p) - tcrossprod(p)

  wvp = w %*% vp
  wvp2 = wvp %*% wvp
  wvp3 = wvp2 %*% wvp
  ec = sum(diag(wvp)) / yplus
  vc = sum(diag(wvp2)) * 2 / yplus ^ 2
  skc = 2 * sqrt(2) * sum(diag(wvp3)) / (sum(diag(wvp2)) ^ (1.5))
  dfc = 8 / skc ^ 2

  tstat.std = (tstat - ec) / sqrt(vc)
  tstat.chisq = dfc + tstat.std * sqrt(2 * dfc)
  pvalue.chisq = 1 - stats::pchisq(tstat.chisq, dfc)
  out = list(tstat = tstat, gof = gof, sa = sa,
             tstat.chisq = tstat.chisq, pvalue.chisq = pvalue.chisq,
             dfc = dfc)
  if (nsim > 0) {
    # simulate new data
    ysim = stats::rmultinom(n = nsim, size = sum(cases),
                             prob = pop / sum(pop))
    # compute r and ee for simulated data
    rsim = ysim / yplus
    eesim = rsim - p

    # compute gof and sa components of statistic for
    # each simulated data
    gof.sim = colSums(eesim ^ 2)
    sa.sim = rowSums(crossprod(eesim, w - diag(N)) * t(eesim))
    tstat.sim = gof.sim + sa.sim
    out$gof.sim = gof.sim
    out$sa.sim = sa.sim
    out$tstat.sim = tstat.sim
    pvalue.sim = (1 + sum(tstat.sim >= tstat)) / (1 + nsim)
    out$pvalue.sim = pvalue.sim
  }
  class(out) = "tango"
  return(out)
}

#' Argument checking for tango.test
#'
#' Check the arguments of the tango.test function
#' @return NULL
#' @noRd
arg_check_tango_test = function(cases, pop, w, nsim) {
  N = length(cases)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_tango_w(w, N)
  arg_check_nsim(nsim)
}

