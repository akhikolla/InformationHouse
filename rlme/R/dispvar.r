#' Rank-based dispersion estimate.
#' 
#' This is an unbiased estimator with a correction factor for standard
#' deviation when normal errors.
#' 
#' @param x vector
#' @param score score type - 1 or 2
#' @references T. P. Hettmansperger and J. W. McKean. Robust Nonparametric
#' Statistical Methods. Chapman Hall, 2012.
#' 
#' @export
dispvar <- function(x, score = 1) {
    n = length(x)
    rx = rank(x, ties.method = c("random"))
    if (score == 1) {
        sc = sqrt(12) * ((rx/(n + 1)) - 0.5)
        dispvar = sqrt(pi/3) * sum(x * sc)/n
    }
    if (score == 2) {
        sc = sign((rx - (n + 1)/2))
        dispvar = sqrt(pi/2) * sum(x * sc)/n
    }
    dispvar
}
