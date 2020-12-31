#' Sum over zones
#'
#' \code{zones.sum} computes the sum of \code{y} for the
#' indices in each element of the list contained in \code{zones}.
#'
#' @param zones A list of nearest neighbors in the format
#' produced by \code{\link{scan.zones}}.
#' @inheritParams nn.cumsum
#'
#' @return A numeric vector.
#' @export
#'
#' @examples
#' # show nn.cumsum example for a circular scan setting
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' cases = floor(nydf$cases)
#' zones = scan.zones(coords, pop = nydf$pop, ubpop = 0.1)
#' # compute cumulative sums over all nn
#' szones = zones.sum(zones, cases)
#' # compute cumulative sums over just the first set of nn
#' szones2 = sapply(zones, function(x) sum(cases[x]))
#' # check equality
#' all.equal(szones, szones2)
zones.sum = function(zones, y) {
  if (!is.list(zones)) stop("zones must be a list")
  if (!is.numeric(y)) stop("y must be a numeric vector")
  vapply(zones, function(x) sum(y[x]), FUN.VALUE = numeric(1),
         USE.NAMES = FALSE)
}
