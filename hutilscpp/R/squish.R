#' @title Squish into a range
#' @name squish
#' @param x A numeric vector.
#' @param a,b Upper and lower bounds
#' @param in_place (logical, default: \code{FALSE}) Should the function operate on \code{x} in place?
#'
#' @return A numeric/integer vector with the values of \code{x} "squished" between \code{a}
#' and \code{b}; values above \code{b}
#' replaced with \code{b} and values below \code{a} replaced with \code{a}.
#'
#' @examples
#' squish(-5:5,-1L, 1L)
#'
#' @export

squish <- function(x, a, b, in_place = FALSE) {
  if (!length(x)) {
    return(x)
  }
  if (length(a) != 1L) {
    stop("`length(a) = ", length(a), "`, but must be length-one.")
  }
  if (length(b) != 1L) {
    stop("`length(b) = ", length(b), "`, but must be length-one.")
  }
  check_TF(in_place)

  if (is.integer(x)) {
    if (is.integer(a) && is.integer(b)) {
      squishi(x, a, b, in_place)
    } else {
      stop("`x` was type integer but `a` was type '", class(a), "' and ",
           "`b` was type '", class(b), "'. ",
           "Coerce `x`, `a`, `b` to the same type.")
    }
  } else if (is.double(x)) {
    squishn(x, as.double(a), as.double(b), in_place)
  } else {
    stop("`x` was type ", typeof(x), ", but must be numeric.")
  }
}
