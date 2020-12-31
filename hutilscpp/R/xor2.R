#' Exclusive or
#' @param x,y Logical vectors.
#' @param anyNAx,anyNAy Could \code{x} and \code{y} possibly contain \code{NA} values?
#' Only set to \code{FALSE} if known to be free of \code{NA}.
#' @export xor2

xor2 <- function(x, y, anyNAx = TRUE, anyNAy = TRUE) {
  stopifnot(length(x) == length(y))
  do_xor2(x, y, anyNAx = !isFALSE(anyNAx), anyNAy = !isFALSE(anyNAy))
}
