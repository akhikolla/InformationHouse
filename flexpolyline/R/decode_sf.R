#' Wrapper function for decoding to simple features
#'
#' A wrapper function for \code{\link{decode}} that converts the input polylines,
#' encoded in the flexible polyline enoding, to simple feature geometries of the sf package.
#'
#' @note
#' The function returns an sf object, therefore the input set of encoded polylines
#' must be of consistent dimension (e.g \code{"XY"}, \code{"XYM"} or \code{"XYZ"})
#' to meet the requirements of the constructor of sf objects. For mixed dimensions
#' use the \code{\link{decode}} function directly.
#'
#' @param encoded character, encoded flexible polyline string.
#' @param crs integer or character, coordinate reference system to assign to the sf object (\code{default = sf::NA_crs_}).
#'
#' @return
#' An \code{sf} object, containing the geometries of the decoded lines (Geometry type: \code{"LINESTRING"}).
#'
#' @export
#'
#' @examples
#' decode_sf("B1Voz5xJ67i1Bgkh9B")
#' decode_sf("BFoz5xJ67i1B1B7PlU9yB")
#' decode_sf("BlXoz5xJ67i1Bgkh9B1B7Pgkh9BzIhagkh9BqK-pB_ni6D")
decode_sf <- function(encoded, crs = sf::NA_crs_) {
  UseMethod("decode_sf", encoded)
}

#' @export
decode_sf.character <- function(encoded, crs = sf::NA_crs_) {
  dim3 <- character(length(encoded))
  ind3 <- 2
  sfdi <- "XY"
  geom <- sf::st_sfc(
    lapply(1:length(encoded), function(x) {
      m <- decode(encoded[[x]])
      d3 <- colnames(m)[3]
      if (is.na(d3)) {
        dim3[x] <<- "ABSENT"
        ind3 <<- 2
        sfdi <<- "XY"
      } else {
        dim3[x] <<- d3
        ind3 <<- 3
        if (d3 %in% c("LEVEL", "ALTITUDE", "ELEVATION")) {
          sfdi <<- "XYZ"
        } else {
          sfdi <<- "XYM"
        }
      }
      if (nrow(m) <= 1) {
        sf::st_point(m, dim = sfdi)
      } else {
        if (all(m[1, ] == m[nrow(m), ])) {
          sf::st_polygon(list(m), dim = sfdi)
        } else {
          sf::st_linestring(m, dim = sfdi)
        }
      }
    }),
    crs = crs
  )
  dim3[is.na(dim3)] <- "ABSENT"
  return(
    sf::st_as_sf(
      data.frame(
        id = seq(1, length(encoded)),
        dim3 = dim3,
        geometry = geom
      )
    )
  )
}

#' @export
decode_sf.factor <- function(encoded, crs = sf::NA_crs_) {
  return(decode_sf.character(as.character(encoded)))
}
