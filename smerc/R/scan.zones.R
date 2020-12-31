#' Determine zones for the spatial scan test
#'
#' \code{scan.zones} determines the unique candidate
#' zones to consider for the circular spatial scan test of
#' Kulldorff (1997).
#'
#' @inheritParams scan.test
#' @return Returns a list of zones to consider for
#'   clustering.  Each element of the list contains a vector
#'   with the location ids of the regions in that zone.
#' @author Joshua French
#' @export
#' @references Kulldorff, M. (1997) A spatial
#' scan statistic. Communications in Statistics - Theory and
#' Methods, 26(6): 1481-1496,
#' <doi:10.1080/03610929708831995>
#' @examples
#' data(nydf)
#' coords = cbind(nydf$longitude, nydf$latitude)
#' zones = scan.zones(coords = coords, pop = nydf$pop,
#'                    ubpop = 0.1, longlat = TRUE)
scan.zones = function(coords, pop, ubpop = 0.5, longlat = FALSE) {
  # argument checking
  arg_check_scan_zones(coords, pop, ubpop, longlat)

  # compute intercentroid distance
  d = sp::spDists(as.matrix(coords), longlat = longlat)

  # for each region, determine sorted nearest neighbors
  # subject to population constraint
  nn = nnpop(d, pop, ubpop)

  # use nearest neighbors to construct all zones
  zones = nn2zones(nn)

  # return only unique zones
  return(zones[distinct(zones)])
}

#' Argument checking for scan.zones
#'
#' Check the arguments of the scan.zones function
#' @param coords A matrix of coordinates
#' @param pop A vector of population for each centroid
#' @param ubpop The population upperbound
#' @param longlat TRUE means great circle distance
#' @return NULL
#' @noRd
# argument checking for all scan zones
arg_check_scan_zones = function(coords, pop,
                                ubpop, longlat) {
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_pop(pop, N)
  arg_check_ubpop(ubpop)
  arg_check_longlat(longlat)
}
