#' Are any values outside the interval specified?
#' @param x A numeric vector.
#' @param a,b Single numeric values designating the interval.
#' @param nas_absent Are \code{NA}s \emph{known} to be absent from \code{x}?
#' If \code{nas_absent = NA}, the default, \code{x} will be searched for \code{NA}s;
#' if \code{nas_absent = TRUE}, \code{x} will not be checked;
#' if \code{nas_absent = FALSE}, the answer is \code{NA_integer_} if \code{na.rm = FALSE}
#' otherwise only non-NA values outside \code{[a, b]}.
#'
#' If \code{nas_absent = TRUE} but \code{x} has missing values then the result is unreliable.
#'
#' @param na_is_outside (logical, default: \code{NA}) How should \code{NA}s in \code{x} be treated?
#' \describe{
#' \item{If \code{NA}}{ the default, then the first value in \code{x} that is either
#' outside \code{[a, b]} or \code{NA} is detected: if it is \code{NA}, then
#' \code{NA_integer_} is returned; otherwise the position of that value is returned.}#'
#' \item{If \code{FALSE}}{then \code{NA} values are effectively skipped; the position of the first
#' \emph{known} value outside \code{[a, b]} is returned.}
#' \item{If \code{TRUE}}{the position of the first value that is either
#' outside \code{[a, b]} or \code{NA} is returned.}
#' }
#' @return \code{0L} if no values in \code{x} are outside \code{[a, b]}. Otherwise, the position
#' of the first value of \code{x} outside \code{[a, b]}.
#'
#' @examples
#' anyOutside(1:10, 1L, 10L)
#' anyOutside(1:10, 1L, 7L)
#'
#' # na_is_outside = NA
#' anyOutside(c(1:10, NA), 1L, 7L)     # Already outside before the NA
#' anyOutside(c(NA, 1:10, NA), 1L, 7L) # NA since it occurred first
#'
#' anyOutside(c(1:7, NA), 1L, 7L, na_is_outside = FALSE)
#' anyOutside(c(1:7, NA), 1L, 7L, na_is_outside = TRUE)
#'
#' ##
#' # N <- 500e6
#' N <- 500e3
#' x <- rep_len(hutils::samp(-5:6, size = 23), N)
#' bench_system_time(anyOutside(x, -5L, 6L))
#' #    process      real
#' #  453.125ms 459.758ms
#'
#' @export


anyOutside <- function(x, a, b, nas_absent = NA, na_is_outside = NA) {
  if (!length(x)) {
    return(0L)
  }
  if (amsg <- isnt_number(a)) {
    stop(attr(amsg, "ErrorMessage"))
  }
  if (bmsg <- isnt_number(b)) {
    stop(attr(bmsg, "ErrorMessage"))
  }

  if (!is.logical(nas_absent)) {
    stop("`nas_absent` was type ", class(nas_absent), ", but must be a length-one logical.")
  }
  if (length(nas_absent) != 1L) {
    stop("`nas_absent` had length ", length(nas_absent), ", but must be length-one.")
  }
  nas_present <-
    if (anyNA(nas_absent)) {
      anyNA(x)
    } else {
      !nas_absent
    }

  if (nas_present && isFALSE(na_is_outside)) {
    if (min(x, na.rm = TRUE) >= a && max(x, na.rm = TRUE) <= b) {
      return(0L)
    }
  }

  if (is.integer(x) && is.integer(a) && is.integer(b)) {
    if (anyNA(na_is_outside)) {
      o <- anyOutside_int(x, a, b, nas_present = nas_present, na_is_outside = TRUE)
      if (o && is.na(x[o])) {
        return(NA_integer_)
      }
    } else {
      o <- anyOutside_int(x, a, b, nas_present = nas_present, na_is_outside = na_is_outside)
    }
    return(o)
  }
  if (is.double(x) && is.double(a) && is.double(b)) {
    if (anyNA(na_is_outside)) {
      o <- anyOutside_dbl(x, a, b, nas_present = nas_present, na_is_outside = TRUE)
      if (o && is.na(x[o])) {
        return(NA_integer_)
      }
    } else {
      o <- anyOutside_dbl(x, a, b, nas_present = nas_present, na_is_outside = na_is_outside)
    }
    return(o)
  }
  stop("`x`, `a`, and `b` must be numeric vectors of the same type.")
}





