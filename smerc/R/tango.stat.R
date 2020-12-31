#' Tango's statistic
#'
#' \code{tango.stat} computes Tango's index (Tango, 1995),
#' including both the goodness-of-fit and spatial
#' autocorrelation components.  See Waller and Gotway
#' (2005).
#' @inheritParams scan.test
#' @param w An \eqn{n \times n} weights matrix.
#'
#' @return Returns a list with the test statistic
#'   (\code{tstat}), the goodness-of-fit component
#'   (\code{gof}), and the spatial autocorrelation component
#'   (\code{sa}).
#' @references Tango, T.  (1995) A class of tests for
#'   detecting "general" and "focused" clustering of rare
#'   diseases.  Statistics in Medicine.  14:2323-2334.
#'
#'   Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial
#'   Statistics for Public Health Data.  Hoboken, NJ: Wiley.
#' @author Joshua French
#' @export
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' w = dweights(coords, kappa = 1, type = "tango", longlat = TRUE)
#' tango.stat(nydf$cases, nydf$pop, w)
tango.stat = function(cases, pop, w) {
  arg_check_tango_stat(cases, pop, w)
  N = length(cases)
  pcases = cases / sum(cases)
  p = pop / sum(pop)
  pd = pcases - p
  gof = sum(pd ^ 2)
  sa = (crossprod(pd, w - diag(N)) %*% pd)[1, 1]
  tstat = gof + sa
  return(list(tstat = tstat, gof = gof, sa = sa))
}

arg_check_tango_stat = function(cases, pop, w) {
  N = length(cases)
  if (!is.numeric(cases)) {
    stop("cases should be a numeric vector")
  }
  if (length(pop) != N) stop("length(pop) != length(cases)")
  if (!is.numeric(pop)) {
    stop("pop should be a numeric vector")
  }
  if (!is.matrix(w)) stop("w should be a matrix")
  if (ncol(w) != N || nrow(w) != N) {
    stop("nrow(w) and ncol(w) != length(cases)")
  }
}
