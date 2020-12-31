#' Form row means conditional on number of non-missing
#'
#' Form row means for multiple vectors, numeric arrays (or data frames) conditional on the number of non-missing observations.
#' NA is returned unless a minimum number of observations is observed.
#'
#' @param ... a series of numeric vectors, arrays, or data frames that have can be combined with cbind
#' @param minobs an integer stating the minimum number of non-NA observations necessary to compute the row mean. Defaults to 1.
#'
#' @return A numeric vector containing the row sums or NA if not enough non-NA observations are present
#'
#' @examples
#' conditional_rowMeans(1:5, c(1:4, NA), c(1:3, NA, NA))
#' conditional_rowMeans(1:5, c(1:4, NA), c(1:3, NA, NA), minobs=0)
#' conditional_rowMeans(1:5, c(1:4, NA), c(1:3, NA, NA), minobs=2)
#'
#' @export
conditional_rowMeans <- function(..., minobs=1L) {
    m <- cbind(...)
    res <- rowMeans(m, na.rm=TRUE)
    res[rowSums(is.na(m))>ncol(m)-minobs] <- NA
    res
}
