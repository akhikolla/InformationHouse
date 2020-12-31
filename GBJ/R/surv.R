#' surv.R
#'
#' Survival (1 minus the CDF) function of standard normal random variable.
#' 
#' @param x Vector of quantiles
#' 
#' @return Probability that a standard normal random variable is greater than x.
#'
#' @export
#' @examples 
#' surv(0)		# Should return 0.5

surv <- function(x) {
	1-pnorm(x, lower.tail=TRUE)
}