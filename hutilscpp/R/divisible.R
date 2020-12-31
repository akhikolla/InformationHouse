#' Divisible
#'
#' @param x An integer vector
#' @param d \code{integer(1)}. The divisor.
#' @param nThread The number of threads to use.
#'
#' @return Logical vector: \code{TRUE} where \code{x} is divisible by \code{d}.
#'
#' \code{divisible16} is short for (and quicker than) \code{divisble(x, 16)}.
#'
#' @export

divisible <- function(x, d, nThread = getOption("hutilscpp.nThread", 1L)) {
  do_divisible(x, d, nThread)
}

#' @rdname divisible
#' @export
divisible16 <- function(x, nThread = getOption("hutilscpp.nThread", 1L)) {
  do_divisible16(x, nThread)
}
