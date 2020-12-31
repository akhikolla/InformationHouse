#' Are elements of a vector even?
#' @param x An integer vector. Double vectors may also be used, but will
#' be truncated, with a warning if any element are not integers.
#' Long vectors are not supported unless \code{x} is integer and \code{keep_nas = FALSE}.
#' @param check_integerish (logical, default: \code{TRUE}) Should the
#' values in \code{x} be checked for non-integer values
#' if \code{x} is a double vector. If \code{TRUE}
#' and values are found to be non-integer a warning is emitted.
#'
#' @param keep_nas (logical, default: \code{TRUE}) Should \code{NA}s in \code{x}
#' return \code{NA} in the result? If \code{FALSE}, will return \code{TRUE} since
#' the internal reprsentation of \code{x} is even. Only applies if \code{is.integer(x)}.
#'
#' @param nThread Number of threads to use.
#'
#'
#' @return For \code{are_even}, a logical vector the same length as \code{x},
#' \code{TRUE} whenever \code{x} is even.
#'
#' For \code{which_are_even} the integer positions of even values in \code{x}.
#'
#'
#'
#'
#'
#' @export are_even which_are_even

are_even <- function(x,
                     check_integerish = TRUE,
                     keep_nas = TRUE,
                     nThread = getOption("hutilscpp.nThread", 1L)) {
  check_TF(keep_nas)
  check_omp(nThread)
  if (!keep_nas && is.integer(x)) {
    return(do_divisible2(x, nThread = nThread))
  }

  wb <- 0L

  if (is.integer(x)) {
    return(do_are_even(x, double(0), wb, nThread))
  }
  if (is.double(x)) {
    if (AND(check_integerish,
            wb <- which_isnt_integerish(x))) {
      warning("`x` was type double, but element ",
              wb, " = ", x[wb], " was not an integer value. ",
              "Will be coerced to integer.")
    }
    do_are_even(integer(0), x, wb, nThread)
  } else {
    stop("`x` was not an integer or double.")
  }
}

#' @rdname are_even
which_are_even <- function(x, check_integerish = TRUE) {
  if (length(x) >= .Machine$integer.max) {
    stop("`which_are_even() not implemented for long `x`.") # nocov
  }
  wb <- 0L
  if (is.integer(x)) {
    return(do_which_even(x, double(0), wb))
  }
  if (is.double(x)) {
    if (AND(check_integerish,
            wb <- which_isnt_integerish(x))) {
      warning("`x` was type double, but element ",
              wb, " = ", x[wb], " was not an integer value. ",
              "Will be coerced to integer: ", as.integer(x[wb]), ".")
    }
    return(do_which_even(integer(0), x, wb))
  } else {
    stop("`x` was not an integer or double.")
  }
}

