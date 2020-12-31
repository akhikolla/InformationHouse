#' @title Match coordinates to nearest coordinates
#' @description When geocoding coordinates to known addresses, an efficient way to
#' match the given coordinates with the known is necessary. This function provides this
#' efficiency by using \code{C++} and allowing approximate matching.
#' @param lat,lon Coordinates to be geocoded. Numeric vectors of equal length.
#' @param addresses_lat,addresses_lon Coordinates of known locations. Numeric vectors of equal length
#' (likely to be a different length than the length of \code{lat}, except when \code{excl_self = TRUE}).
#' @param Index A vector the same length as \code{lat} to encode the match between \code{lat,lon}
#' and \code{addresses_lat,addresses_lon}. The default is to use the integer position
#' of the nearest match to
#' \code{addresses_lat,addresses_lon}.
#' @param cartesian_R The maximum radius of any address from the points to be geocoded. Used
#' to accelerate the detection of minimum distances. Note, as the argument name suggests,
#' the distance is in cartesian coordinates, so a small number is likely.
#'
#' @param close_enough The distance, in metres, below which a match will be considered to have occurred.
#' (The distance that is considered "close enough" to be a match.)
#'
#' For example, \code{close_enough = 10} means the first location within ten metres will be matched,
#' even if a closer match occurs later.
#'
#' May be provided as a string to emphasize the units, e.g. \code{close_enough = "0.25km"}.
#' Only \code{km} and \code{m} are permitted.
#'
#' @param excl_self (bool, default: \code{FALSE}) For each \eqn{x_i} of the first coordinates,
#' exclude the \eqn{y_i}-th point when determining closest match. Useful to determine the
#' nearest neighbour within a set of coordinates, \emph{viz.}
#' \code{match_nrst_haversine(x, y, x, y, excl_self = TRUE)}.
#'
#' @param as.data.table Return result as a \code{data.table}?
#' If \code{FALSE}, a list is returned. \code{TRUE} by default to
#' avoid dumping a huge list to the console.
#'
#' @param .verify_box Check the initial guess against other points within the
#' box of radius \eqn{\ell^\infty}.
#'
#'
#' @return A list (or \code{data.table} if \code{as.data.table = TRUE}) with two elements,
#' both the same length as \code{lat}, giving for point \code{lat,lon}:
#' \describe{
#' \item{\code{pos}}{the position (or corresponding value in \code{Table})
#'   in \code{addresses_lat,addresses_lon} nearest to \code{lat, lon}.}
#' \item{\code{dist}}{the distance, in kilometres, between the two points.}
#' }
#'
#' @examples
#' lat2 <- runif(5, -38, -37.8)
#' lon2 <- rep(145, 5)
#'
#' lat1 <- c(-37.875, -37.91)
#' lon1 <- c(144.96, 144.978)
#'
#' match_nrst_haversine(lat1, lon1, lat2, lon2, 0L)
#' match_nrst_haversine(lat1, lon1, lat1, lon1, 11:12, excl_self = TRUE)
#'
#' @export match_nrst_haversine
#'
#'

match_nrst_haversine <- function(lat,
                                 lon,
                                 addresses_lat,
                                 addresses_lon,
                                 Index = seq_along(addresses_lat),
                                 cartesian_R = NULL,
                                 close_enough = 10,
                                 excl_self = FALSE,
                                 as.data.table = TRUE,
                                 .verify_box = TRUE) {

  stopifnot(is.numeric(lat),
            is.numeric(lon),
            length(lat) == length(lon),
            length(lat) >= 1L,
            is.numeric(addresses_lat),
            is.numeric(addresses_lon),
            length(addresses_lat) == length(addresses_lon),
            length(addresses_lat) >= 1L)
  check_TF(excl_self)
  check_TF(as.data.table)

  if (ce_err_msg <- isnt_number(close_enough)) {
    if (length(close_enough) != 1L || !grepl("^[0-9]+(\\.[0-9]+)?\\s*k?m$", close_enough)) {
      stop(attr(ce_err_msg, "ErrorMessage"), "\n",
           "`close_enough` may be a string of numbers ",
           "ending in 'km' or 'm' to designate units.")
    }
    dist0_km <- units2km(close_enough)
  } else {
    dist0_km <- close_enough / 1000
  }


  min_lat <- min(lat)
  max_lat <- max(lat)
  min_lon <- min(lon)
  max_lon <- max(lon)
  verify_cartR <- FALSE
  if (is.null(cartesian_R)) {
    if (min_lat > -63 && max_lat < 63) {
      # Within 63 degrees of latitude, the cartesian distance
      # is at least 50 times the distance in km.
      cartesian_R <- 0.02  # within a km
      verify_cartR <- TRUE
    }
  }





  if (length(Index) != length(addresses_lat)) {
    warning("`Table` was provided, but was not the same ",
            "length as `addresses_lat`, ",
            "so returning positions.")
    Index <- seq_along(addresses_lat)
  }
  Index.int <-
    if (!is.integer(Index)) {
      seq_along(Index)
    } else {
      Index
    }

  out <- match_min_Haversine(lat,
                             lon,
                             addresses_lat,
                             addresses_lon,
                             Index.int,
                             cartR = cartesian_R,
                             verify_cartR = verify_cartR,
                             do_verify_box = .verify_box,
                             dist0_km = dist0_km,
                             ncores = 1L)
  if (!is.integer(Index)) {
    out[[1L]] <- Index[out[[1L]]]
  }

  if (as.data.table) {
    return(data.table::as.data.table(out))
  } else {
    return(out)
  }
}


