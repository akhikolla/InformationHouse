#' Cumulative sum unless reset
#' @param x A logical vector indicating when the sum should \emph{continue}.
#' @param y Optional: a numeric vector the same length as \code{x} to cumulatively sum.
#' @return If \code{y} is a double vector, a double vector of cumulative sums,
#' resetting whenever \code{x} is \code{FALSE}; otherwise an integer vector.
#'
#' If \code{length(x) == 0}, \code{y} is returned (i.e. \code{integer(0)} or \code{double(0)}.
#'
#' @examples
#' cumsum_reset(c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE))
#' cumsum_reset(c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
#'              c(1000, 1000, 10000,   10,   20,   33,     0))
#'
#' @export

cumsum_reset <- function(x, y = as.integer(x)) {
  if (!is.logical(x)) {
    stop("`x` must be a logical vector.")
  }
  if (anyNA(x)) {
    stop("`x` has missing values. This is not permitted.")
  }
  if (missing(y)) {
    if (length(x)) {
      return(do_cumsum_reset_logical(x))
    } else {
      return(integer(0L))
    }
  }
  if (length(y) != length(x)) {
    stop("`y` had length ", length(y), ", yet `length(x) = ", length(x), ". ",
         "`x` and `y` must have the same lengths. ",
         "Provide vectors of equal length.")
  }
  if (!length(y)) {
    return(y)
  }
  if (is.integer(y)) {
    return(do_cumsum_reset_integer(x, y))
  }
  if (is.double(y)) {
    return(do_cumsum_reset_double(x, y))
  }
  stop("`y` was type ", class(y), ", but must be an integer or double.")
}


