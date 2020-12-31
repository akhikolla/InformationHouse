#' Simulate values from a perfect normal distribution
#' 
#' Random generation for a perfect normal distribution with mean equal to mean and standard deviation equal to sd.
#' 
#' The function will return the same set of quantiles for fixed n. In that sense there is not much randomness going on, and the function is mostly useful for illustrative purposes.
#' 
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param mean number of mean.
#' @param sd number of standard deviation.
#' @return Returns a vector of values from a perfect normal distribution
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords hplot
#' @examples
#'
#' rnorm_perfect(30, mean=10, sd=2)
#' 
#' @export 
rnorm_perfect <- function(n, mean=0, sd=1) {

    stopifnot(is.numeric(n))
    stopifnot(length(mean)==1)
    stopifnot(length(sd)==1)


    if (length(n)>1) {
        n <- length(n)
    }

    stopifnot(n>0)

    (qnorm(seq(0, 1, length=n+2)[-c(1, n+2)]) + mean)*sd
}
