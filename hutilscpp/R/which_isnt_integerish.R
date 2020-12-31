



which_isnt_integerish <- function(x) {
  if (is.integer(x) || length(x) == 0L) {
    return(0L)
  }
  if (!is.double(x)) {
    return(1L)
  }
  if (is_altrep(x)) {
    if (anyNA(x)) {
      return(1L) # nocov
    }
    if (first(x) > 2147483647 || first(x) < -2147483647) {
      return(1L)
    }
    if (last(x) > 2147483647) {
      o <- .Machine$integer.max - first(x) + 2L
    }
    if (last(x) < -2147483647) {
      o <- .Machine$integer.max + first(x) + 2L
    }
  } else {
    o <- do_which_isnt_integerish(x)
  }
  if (o <= 2147483647) {
    o <- as.integer(o)
  }
  o
}
