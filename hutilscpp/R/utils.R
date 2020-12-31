check_TF <- function(x) {
  if (is.logical(x) && length(x) == 1L) {
    if (anyNA(x)) {
      xc <- deparse(substitute(x))
      stop("`", xc, " = NA` but must be TRUE or FALSE. ",
           "Change `", xc, "` to be TRUE or FALSE.")
    } else {
      return(NULL)
    }
  } else {
    xc <- deparse(substitute(x))
    if (length(x) != 1L) {
      stop("`", xc, "` had length ", length(x), " but must be length-one. ",
           "Change `", xc, "` to be TRUE or FALSE.")
    } else {
      stop("`", xc, "` was type ", typeof(x), " but must be logical. ",
           "Change `", xc, "` to be TRUE or FALSE.")
    }
  }
}

isnt_number <- function(a, na.bad = TRUE, infinite.bad = TRUE, int.only = FALSE) {
  if (!is.numeric(a)) {
    o <- TRUE
    ac <- deparse(substitute(a))
    attr(o, "ErrorMessage") <- paste0("`", ac, "` was a ", class(a), ", but must be numeric.")
    return(o)
  }
  if (length(a) != 1L) {
    o <- TRUE
    ac <- deparse(substitute(a))
    attr(o, "ErrorMessage") <- paste0("`", ac, "` had length ", length(a), ", but must be length-one.")
    return(o)
  }
  if (na.bad && is.na(a)) {
    o <- TRUE
    ac <- deparse(substitute(a))
    attr(o, "ErrorMessage") <- paste0("`", ac, "= NA`, but this is not permitted.")
    return(o)
  }
  if (infinite.bad && is.infinite(a)) {
    o <- TRUE
    ac <- deparse(substitute(a))
    attr(o, "ErrorMessage") <- paste0("`", ac, "` was not finite, but this is not permitted.")
    return(o)
  }
  if (int.only && !is.integer(a)) {
    if (is.nan(a)) {
      o <- TRUE
      ac <- deparse(substitute(a))
      attr(o, "ErrorMessage") <- paste0("`", ac, "` was not safely coercible to integer (NaN).")
      return(o)
    }
    if (is.na(a)) {
      return(FALSE)
    }
    if ((a > 2147483647) || (a < -2147483647)) {
      o <- TRUE
      ac <- deparse(substitute(a))
      attr(o, "ErrorMessage") <- paste0("`", ac, " = ", a, "` was not safely coercible to integer (out of range).")
      return(o)
    }
    if (abs(as.integer(a) - a) > sqrt(.Machine$double.eps)) {
      o <- TRUE
      ac <- deparse(substitute(a))
      attr(o, "ErrorMessage") <- paste0("`", ac, " = ", a, "` was not safely coercible to integer (not a whole number).")
      return(o)
    }
  }
  FALSE
}

AND <- `&&`
OR <- `||`

isFALSE <- function(x) {
  is.logical(x) && length(x) == 1L && !anyNA(x) && !x
}


firstNonNegativeRadix <- function(x, ...) {
  if (is.double(x)) {
    do_firstNonNegativeRadix_dbl(x, ...)
  } else {
    do_firstNonNegativeRadix_int(x, ...)
  }
}

g <- glue::glue


is_wholer <- function(dbl) {
  length(dbl) == 1L &&
  !is.na(dbl) &&
  dbl >= -2147483647 &&
  dbl <= 2147483647 &&
  dbl == as.integer(dbl)
}

# quiet double to int -- when passed to a C++ function that
# accepts int but only conditionally uses
qd2i <- function(x) {
  if (is_safe2int(x)) {
    as.integer(x)
  } else {
    NA_integer_
  }
}

# nocov start
is64bit <- function() .Machine$sizeof.pointer == 8L

is_covr <- function() {
  requireNamespace("covr", quietly = TRUE) &&
    requireNamespace("testthat", quietly = TRUE) &&
    covr::in_covr()
}

skip_if_covr <- function() {
  if (requireNamespace("testthat", quietly = TRUE)) {
    testthat::skip_if(is_covr())
  }
}

# nocov end


