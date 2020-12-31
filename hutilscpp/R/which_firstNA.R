#' First/last position of missing values
#' @description Introuduced in \code{v 1.6.0}
#'
#' @param x An atomic vector.
#' @return The position of the first/last missing value in \code{x}.
#'
#' @examples
#' N <- 1e8
#' N <- 1e6  # for CRAN etc
#' x <- c(1:1e5, NA, integer(N))
#' bench_system_time(which.max(is.na(x))) # 123ms
#' bench_system_time(Position(is.na, x))  #  22ms
#' bench_system_time(which_firstNA(x))    #  <1ms

#' @export
which_firstNA <- function(x) {
  stopifnot(is.atomic(x))
  return(R_xlen_t(do_which_firstNA(x)))
}

#' @rdname which_firstNA
#' @export
which_lastNA <- function(x) {
  stopifnot(is.atomic(x))
  return(R_xlen_t(do_which_lastNA(x)))
}

