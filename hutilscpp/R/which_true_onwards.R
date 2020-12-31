#' At which point are all values true onwards
#' @param x A logical vector. \code{NA} values are not permitted.
#' @return The position of the first \code{TRUE} value in \code{x} at which all
#' the following values are \code{TRUE}.
#'
#' @examples
#' which_true_onwards(c(TRUE, FALSE, TRUE, TRUE, TRUE))
#'
#' @export

which_true_onwards <- function(x) {
  if (!is.logical(x)) {
    stop("`x` was type ", class(x), ". `x` must be a logical vector.")
  }
  if (anyNA(x)) {
    stop("`x` had missing values. This is not permitted.")
  }

  o <- do_which_true_onwards(x)
  if (o < .Machine$integer.max) {
    o <- as.integer(o)
  }
  o
}
