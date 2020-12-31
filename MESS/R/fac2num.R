#' Convert factor to numeric vector
#'
#' Converts the factor labels to numeric values and returns the factor as a numeric vector
#'
#' Returns a vector of numeric values. Elements in the input factor that cannot be converted to numeric will produce NA.
#'
#' @param x A factor
#' @return Returns a numeric vector of the same length as x
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords manip
#' @examples
#'
#' f <- factor(c(1,2,1,3,2,1,2,3,1))
#' fac2num(f)
#'
#' @export
fac2num <- function(x) {
    if (!is.factor(x))
        stop("input needs to be a factor")
    as.numeric(levels(x))[x]
}

