#' qnorm_mu.R
#'
#' Internal function to calculate Pr(-t < Z < t) - kkk/d for Z~N(mu,1).
#' Take the root of this function to find the mu such that P(|Z|>=t_k) = k/d.
#' 
#' @param mu The mean of the normal random variable.
#' @param t The threshold (boundaries) we are interested in.
#' @param kkk We decided to make you input the fraction in two parts.
#' @param d We decided to make you input the fraction in two parts.
#' 
#' @return Pr(-t < Z < t) - kkk/d for Z~N(mu,1).
#'
#' @keywords internal
#' @export
#'
#' @examples 
#' qnorm_mu(mu=0, t=1.96, kkk=1, d=5)		# Should return 0

##################################################################
# Use this to find the mu such that P(|Z|>=t_k)=k/d
# Calculates the normal cdf of t_k for a given mu
qnorm_mu <- function(mu, t, kkk, d) {
	1 - ( pnorm(t, mean=mu, lower.tail=TRUE) - pnorm(-t, mean=mu, lower.tail=TRUE) ) - kkk/d
}