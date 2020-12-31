#' Combine distinct zones
#'
#' \code{combine.zones} combines the elements of \code{z1}
#' and \code{z2} into a single list, returning only the
#' unique zones.
#'
#' @param z1 A list of zones
#' @param z2 A list of zones
#'
#' @return A list of distinct zones
#' @export
#'
#' @examples
#' z1 = list(1:2, 1:3)
#' z2 = list(2:1, 1:4)
#' combine.zones(z1, z2)
combine.zones = function(z1, z2) {
  if (!is.list(z1) | !is.list(z2)) {
    stop("z1 and z2 should be lists of zones")
  }
  z1 = c(z1, z2)
  z1[distinct(z1)]
}
