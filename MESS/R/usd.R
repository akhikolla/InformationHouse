#' Unbiased standard deviation
#'
#' This function computes the unbiased standard deviation of the values in x. If \code{na.rm} is TRUE then missing values are removed before computation proceeds.
#'
#' Like var this uses denominator n - 1.
#' The standard deviation of a length-one or zero-length vector is NA.
#' 
#' @param x a numeric vector or an R object but not a factor coercible to numeric by as.double(x)
#' @param na.rm logical. Should missing values be removed?
#'
#' @return A scalar
#'
#' @examples
#' sd(1:5)
#' usd(1:5)
#'
#' @export
usd <- function(x, na.rm = FALSE) {

  # Get the number of observations used for the computation
  xx <- if (is.vector(x) || is.factor(x)) x else as.double(x)
  n <- ifelse(na.rm, sum(!is.na(xx)), length(xx))

  sd(x, na.rm=na.rm) / (1 - 1/(4*n) - 7/(32*n^2) - 19/(128*n^3)) 

}
