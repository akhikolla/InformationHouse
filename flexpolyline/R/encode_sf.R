#' Wrapper function for encoding simple features
#'
#' A wrapper function for \code{\link{encode}} that converts simple feature geometries
#' of the sf package to flexible polyline encoded strings.
#'
#' @param geom simple feature, \code{sf}, \code{sfc} or \code{sfg} object with geometry type \code{"POINT"}, \code{"LINESTRING"} or \code{"POLYGON"}.
#' @param precision integer, precision to use in encoding (between 0 and 15, \code{default=5}).
#' @param third_dim integer, type of the third dimension (0: ABSENT, 1: LEVEL, 2: ALTITUDE, 3: ELEVATION, 4, 6: CUSTOM1, 7: CUSTOM2, \code{default=NULL}).
#' @param third_dim_precision integer, precision to use in encoding for the third dimension (between 1 and 15, \code{default=precision}).
#'
#' @return
#' The line as string in the flexible polyline encoding format.
#'
#' @export
#'
#' @examples
#' # 3D point
#' point3d <- sf::st_point(
#'   matrix(c(8.69821, 50.10228, 10), ncol = 3, byrow = TRUE), dim = "XYZ")
#' encode_sf(point3d)
#'
#' # 2D linestring
#' line2d <- sf::st_linestring(
#'   matrix(c(8.69821, 50.10228,
#'            8.69567, 50.10201,
#'            8.68752, 50.09878), ncol = 2, byrow = TRUE))
#' encode_sf(line2d)
#'
#' # 3D polygon
#' poly3d <- sf::st_polygon(list(
#'   matrix(c(8.69821, 50.10228, 10,
#'            8.69567, 50.10201, 20,
#'            8.69150, 50.10063, 30,
#'            8.69821, 50.10228, 10), ncol = 3, byrow = TRUE)), dim = "XYM")
#' encode_sf(poly3d)
encode_sf <- function(geom, precision = 5, third_dim = NULL,
                      third_dim_precision = precision) {
  UseMethod("encode_sf", geom)
}

#' @export
encode_sf.sfg <- function(geom, precision = 5, third_dim = NULL,
                       third_dim_precision = precision) {
  if(!sf::st_geometry_type(geom) %in% c("POINT", "LINESTRING", "POLYGON")){
    stop(
      sprintf(
        "Invalid geometry type '%s' of input, only 'POINT', 'LINESTRING' and 'POLYGON' is supported.",
        sf::st_geometry_type(geom)
      )
    )
  }
  if (class(geom)[1] == "XY") {
    third_dim <- 0
    encoded <- encode(
      sf::st_coordinates(geom)[, c(1:2), drop = FALSE],
      precision, third_dim, third_dim_precision
    )
  } else {
    if (is.null(third_dim)) {
      if (class(geom)[1] == "XYZ") {
        third_dim <- 3
      } else {
        third_dim <- 6
      }
    }
    encoded <- encode(
      sf::st_coordinates(geom)[, c(1:3), drop = FALSE],
      precision, third_dim, third_dim_precision
    )
  }
  return(encoded)
}

#' @export
encode_sf.sfc <- function(geom, precision = 5, third_dim = NULL,
                          third_dim_precision = precision) {
  return(
    sapply(geom, function(x) {
      encode_sf.sfg(x, precision, third_dim, third_dim_precision)
    })
  )
}

#' @export
encode_sf.sf <- function(geom, precision = 5, third_dim = NULL,
                         third_dim_precision = precision) {
  return(
    encode_sf.sfc(
      sf::st_geometry(geom),
      precision, third_dim, third_dim_precision
    )
  )
}
