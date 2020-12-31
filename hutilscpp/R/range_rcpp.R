#' Range C++
#' @description Range of a vector using Rcpp.
#' @param x A vector for which the range is desired. Vectors with missing values
#' are not supported and have no definite behaviour.
#' @param anyNAx (logical, default: \code{anyNA(x)} lazily). Set to \code{TRUE}
#' only if \code{x} is known to contain no missing values (including \code{NaN}).
#'
#' @param warn_empty (logical, default: \code{TRUE}) If \code{x} is
#' empty (i.e. has no length), should a warning be emitted (like \code{\link[base]{range}})?
#'
#' @param integer0_range_is_integer (logical, default: \code{FALSE})
#' If \code{x} is a length-zero integer, should the result also be an integer?
#' Set to \code{FALSE} by default in order to be compatible with \code{\link[base]{range}}, but can
#' be set to \code{TRUE} if an integer result is desired, in which case
#' \code{range_rcpp(integer())} is \code{(INT_MAX, -INT_MAX)}.
#'
#' @return A length-4 vector, the first two positions give the range and
#' the next two give the positions in \code{x} where the max and min occurred.
#'
#' This is almost equivalent to \code{c(range(x), which.min(x), which.max(x))}.
#' Note that the type is not strictly preserved, but no loss should occur. In particular,
#' logical \code{x} results in an integer result, and a double \code{x} will
#' have double values for \code{which.min(x)} and \code{which.max(x)}.
#'
#' A completely empty, logical \code{x} returns \code{c(NA, NA, NA, NA)} as an integer vector.
#'
#' @examples
#' x <- rnorm(1e3) # Not noticeable at this scale
#' bench_system_time(range_rcpp(x))
#' bench_system_time(range(x))
#'
#'
#'
#' @export range_rcpp

range_rcpp <- function(x,
                       anyNAx = anyNA(x),
                       warn_empty = TRUE,
                       integer0_range_is_integer = FALSE) {
  if (!length(x)) {
    if (is.integer(x) && integer0_range_is_integer) {
      if (warn_empty) {
        warning("no non-missing arguments to range_rcpp; returning c(INT_MAX, -INT_MAX).")
      }
      out <- c(.Machine$integer.max, -.Machine$integer.max)
    } else {
      if (warn_empty) {
        warning("no non-missing arguments to range_rcpp; returning c(Inf, -Inf).")
      }
      out <- c(Inf, -Inf)
    }
    return(out)
  }
  if (is.integer(x)) {
    return(R_xlen_t(do_range_int(x)))
  }
  if (is.double(x)) {
    if (anyNAx) {
      return(do_range_dbl(x))
    } else {
      return(do_range_dbl_simple(x))
    }
  }
  if (is.logical(x)) {
    if (any(x, na.rm = TRUE)) {
      if (all(x, na.rm = TRUE)) {
        o <- c(TRUE, TRUE, rep(do_which_first(x), 2L))
      } else {
        o <- c(FALSE, TRUE, do_which_first_false(x), do_which_first(x))
      }
    } else {
      ff <- do_which_first_false(x)
      if (ff) {
        o <- c(FALSE, FALSE, rep(ff, 2L))
      } else {
        o <- as.integer(c(NA, NA, NA, NA))
      }
    }
    return(o)
  }
  range(x)
}


