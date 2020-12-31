#' Find a binary pole of inaccessibility
#' @name poleInaccessibility
#' @param x,y Coordinates.
#' @param DT A \code{data.table} containing \code{LONGITUDE} and \code{LATITUDE} to define
#' the \code{x} and \code{y} coordinates.
#' @param x_range,y_range Numeric vectors of length-2; the range of \code{x} and \code{y}.
#' Use this rather than the default when the 'vicinity' of \code{x,y} is different from
#' the minimum closed rectangle covering the points.
#' @param test_both (logical, default: \code{TRUE}) For \code{3}, test both stretching vertically
#' then horizontally and horizontally then vertically.
#' @param copy_DT (logical, default: \code{TRUE}) Run \code{\link[data.table]{copy}} on \code{DT}
#' before proceeding. If \code{FALSE}, \code{DT} have additional columns updated by reference.
#' @return
#' \describe{
#' \item{\code{poleInaccessibility2}}{
#' A named vector containing the
#' \code{xmin}, \code{xmax} and
#' \code{ymin}, \code{ymax} coordinates of
#' the largest rectangle of width an integer power of two that is empty.}
#' \item{\code{poleInaccessibility3}}{
#' Starting with the rectangle formed by \code{poleInaccessibility2},
#' the rectangle formed by stretching it out vertically and horizontally until
#' the edges intersect the points \code{x,y}
#' }
#' }
#'
#' @examples
#' library(data.table)
#' library(hutils)
#' # A square with a 10 by 10 square of the northeast corner removed
#' x <- runif(1e4, 0, 100)
#' y <- runif(1e4, 0, 100)
#' DT <- data.table(x, y)
#' # remove the NE corner
#' DT_NE <- DT[implies(x > 90, y < 89)]
#' DT_NE[, poleInaccessibility2(x, y)]
#' DT_NE[, poleInaccessibility3(x, y)]
#'
#' @export poleInaccessibility2 poleInaccessibility3


#' @rdname poleInaccessibility
poleInaccessibility2 <- function(x = NULL,
                                 y = NULL,
                                 DT = NULL,
                                 x_range = NULL,
                                 y_range = NULL,
                                 copy_DT = TRUE) {
  if (is.null(DT) && is.null(x) && is.null(y)) {
    stop("`DT`, `x`, and `y` were all NULL. ",
         "`DT` or both `x` and `y` must be provided.")
  }
  if (is.null(DT)) {
    DT <- data.table(x, y)
    setnames(DT,
             c("x", "y"),
             c("LONGITUDE", "LATITUDE"))
  } else {
    if (!hasName(DT, "LATITUDE") || !hasName(DT, "LONGITUDE")) {
      stop("`DT` was supplied but did not have both a column 'LATITUDE' and a column 'LONGITUDE'.")
    }

    if (!is.null(x) || !is.null(y)) {
      warning("`DT` was provided, but `x` and `y` are not both NULL and will be ignored.")
    }
    if (copy_DT) {
      DT <- copy(DT)
    }
  }
  # CRAN NOTE avoidance
  LATITUDE <- LONGITUDE <- NULL

  # Get the ranges so we can remap to the original scale
  # after we identify the 'binary' position of the empty ball
  if (is.null(x_range)) {
    x_range <- DT[, range_rcpp(LONGITUDE)]
    global_minx <- x_range[1]
    global_maxx <- x_range[2]
  } else {
    global_minx <- x_range[1]
    global_maxx <- x_range[2]
  }
  if (is.null(y_range)) {
    y_range <- DT[, range_rcpp(LATITUDE)]
    global_miny <- y_range[1]
    global_maxy <- y_range[2]
  } else {
    global_miny <- y_range[1]
    global_maxy <- y_range[2]
  }
  npoints <- DT[, .N]
  d <- 0L
  while (npoints > 0) {
    d <- d + 1L
    cut_DT(DT, d)
    n_by_quarter <-
      DT[, .N, keyby = c(paste0("xbreaks", d), paste0("ybreaks", d))]
    if (n_by_quarter[, .N] < 4^d) {
      npoints <- 0L
      DT_expand <- CJ(xbreaks = seq_len(2^d),
                      ybreaks = seq_len(2^d))
      setnames(DT_expand,
               c("xbreaks", "ybreaks"),
               c(paste0("xbreaks", d), paste0("ybreaks", d)))
      setkeyv(DT_expand, key(n_by_quarter))

      zero_by_quarter <-
        # n = 1L because we only want/need the first?
        head(DT_expand[!n_by_quarter], 1L) # anti-join to identify

    }
  }


  # This returns the 'scaled' versions of each rectangle's vertex (0, 0) - (1, 1)
  z_by_q <- zero_by_quarter[, lapply(.SD, function(x) c({x - 1} / 2^d, x / 2^d))]

  unscale <- function(p, v0, v1) {
    v0 + (v1 - v0) * p
  }
  # Need to then unscale
  c("xmin" = unscale(z_by_q[1][[1]], global_minx, global_maxx),
    "xmax" = unscale(z_by_q[2][[1]], global_minx, global_maxx),
    "ymin" = unscale(z_by_q[1][[2]], global_miny, global_maxy),
    "ymax" = unscale(z_by_q[2][[2]], global_miny, global_maxy))
}


cut_DT <- function(DT, depth = 1L, x_range = NULL, y_range = NULL) {
  if (anyNA(match(c("LATITUDE", "LONGITUDE"), names(DT), nomatch = NA_integer_))) {
    stop("`DT` lacked columns 'LATITUDE' and 'LONGITUDE'.")
  }
  LONGITUDE <- LATITUDE <- NULL
  if (is.null(x_range)) {
    x_range <- range_rcpp(.subset2(DT, "LONGITUDE"))
  }
  if (is.null(y_range)) {
    y_range <- range_rcpp(.subset2(DT, "LATITUDE"))
  }
  DT[, "xbreaks" := .bincode(LONGITUDE,
                             include.lowest = TRUE,
                             breaks = seq.int(from = x_range[1],
                                              to = x_range[2],
                                              length.out = 2^depth + 1))]
  setnames(DT, "xbreaks", paste0("xbreaks", depth))
  DT[, "ybreaks" := .bincode(LATITUDE,
                             include.lowest = TRUE,
                             breaks = seq.int(from = y_range[1],
                                              to = y_range[2],
                                              length.out = 2^depth + 1))]
  setnames(DT, "ybreaks", paste0("ybreaks", depth))
  DT[]
}

#' @rdname poleInaccessibility
poleInaccessibility3 <- function(x = NULL,
                                 y = NULL,
                                 DT = NULL,
                                 x_range = NULL,
                                 y_range = NULL,
                                 copy_DT = TRUE,
                                 test_both = TRUE) {
  if (is.null(DT) && is.null(x) && is.null(y)) {
    stop("`DT`, `x`, and `y` were all NULL. ",
         "`DT` or both `x` and `y` must be provided.")
  }
  if (is.null(DT)) {
    DT <- data.table(x, y)
    setnames(DT,
             c("x", "y"),
             c("LONGITUDE", "LATITUDE"))
  } else {
    if (anyNA(match(c("LATITUDE", "LONGITUDE"), names(DT), nomatch = NA_integer_))) {
      stop("`DT` was supplied but did not have both a column 'LATITUDE' and a column 'LONGITUDE'.")
    }

    if (!is.null(x) || !is.null(y)) {
      warning("`DT` was provided, but `x` and `y` are not both NULL and will be ignored.")
    }
    if (copy_DT) {
      DT <- copy(DT)
    }
  }
  # CRAN NOTE avoidance
  LATITUDE <- LONGITUDE <- NULL
  if (is.null(x_range)) {
    x_range <- DT[, range_rcpp(LONGITUDE)]
    global_minx <- x_range[1]
    global_maxx <- x_range[2]
  } else {
    global_minx <- x_range[1]
    global_maxx <- x_range[2]
  }
  if (is.null(y_range)) {
    y_range <- DT[, range_rcpp(LATITUDE)]
    global_miny <- y_range[1]
    global_maxy <- y_range[2]
  } else {
    global_miny <- y_range[1]
    global_maxy <- y_range[2]
  }

  p2 <- poleInaccessibility2(DT = DT,
                             x_range = x_range,
                             y_range = y_range,
                             # Already copied
                             copy_DT = FALSE)
  stopifnot(hasName(DT, "LONGITUDE"))
  stopifnot(hasName(DT, "LATITUDE"))
  LONGITUDE <- LATITUDE <- NULL
  # range_dbl internals has four, which
  # will conflict with %between%
  y_box <- c(p2["ymin"], p2["ymax"])
  stopifnot(length(y_box) == 2L)

  # suppress warnings about no non-missing arguments to min
  # we don't mind +/-Inf as they will be excluded in the if
  # statements ante
  suppressWarnings({
    # Look outward, west and east, stop once we hit any point
    xmin_new <-  DT[LATITUDE %between% c(y_box) & LONGITUDE < p2["xmin"], max(LONGITUDE)]
    xmax_new <-  DT[LATITUDE %between% c(y_box) & LONGITUDE > p2["xmax"], min(LONGITUDE)]

    if (xmin_new < global_minx) {
      xmin_new <- global_minx
    }
    if (xmax_new > global_maxx) {
      xmax_new <- global_maxx
    }
    x_box <- c(xmin_new, xmax_new)
    # We need to take into account the new xmins immediately to ensure we take the
    # *intersection* of empty boxes
    ymin_new <- DT[LONGITUDE %between% c(x_box) & LATITUDE < p2["ymin"], max(LATITUDE)]
    ymax_new <- DT[LONGITUDE %between% c(x_box) & LATITUDE > p2["ymax"], min(LATITUDE)]

  })
  if (ymin_new < global_miny) {
    ymin_new <- global_miny
  }
  if (ymax_new > global_maxy) {
    ymax_new <- global_maxy
  }
  out <- c("xmin" = xmin_new,
           "xmax" = xmax_new,
           "ymin" = ymin_new,
           "ymax" = ymax_new)
  A <- (out[2] - out[1]) * (out[4] - out[3])

  if (test_both) {
    x_box <- c(p2["xmin"], p2["xmax"])
    suppressWarnings({
      ymin_new <- DT[LONGITUDE %between% c(x_box) & LATITUDE < p2["ymin"], max(LATITUDE)]
      ymax_new <- DT[LONGITUDE %between% c(x_box) & LATITUDE > p2["ymax"], min(LATITUDE)]

      if (ymin_new < global_miny) {
        ymin_new <- global_miny
      }
      if (ymax_new > global_maxy) {
        ymax_new <- global_maxy
      }
      y_box <- c(ymin_new, ymax_new)

      xmin_new <-  DT[LATITUDE %between% c(y_box) & LONGITUDE < p2["xmin"], max(LONGITUDE)]
      xmax_new <-  DT[LATITUDE %between% c(y_box) & LONGITUDE > p2["xmax"], min(LONGITUDE)]
    })
    if (xmin_new < global_minx) {
      xmin_new <- global_minx
    }
    if (xmax_new > global_maxx) {
      xmax_new <- global_maxx
    }
    out2 <- c("xmin" = xmin_new,
              "xmax" = xmax_new,
              "ymin" = ymin_new,
              "ymax" = ymax_new)
    A2 <- (out2[2] - out2[1]) * (out2[4] - out2[3])
    if (A2 > A) {
      return(out2)
    }
  }
  out
}




