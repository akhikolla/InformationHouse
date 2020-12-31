#' ebb_loglik.R
#'
#' An internal function only (doesn't do error checking or take care of boundary cases).
#' The log-likelihood (log-PMF) calculation for the Extended Beta-Binomial proposed by Prentice (1986).
#' Takes in a vectorized argument because we apply() it in run_BB_GBJ().
#' 
#' @param x A vector of length 3 with (1) value of outcome (2) mu parameter (3) gamma parameter
#' @param d The number of test stsatistics in the set
#' 
#' @keywords internal
#' @export
#'
#' @return log( Pr(V=x[1]) ) where V~EBB(mu,gamma; d)
#'
#' @examples 
#' ebb_loglik(x=c(1, 0.5, 0.1), d=10)	

ebb_loglik <- function(x, d)
{
	vec1 <- 0:(x[1]-1)				# x should never be zero in our application, don't have to check that case
	vec2 <- 0:(d-x[1]-1)			# x should also never be d since we are going from 1:(d/2)
	vec3 <- 0:(d-1)				 
	lchoose(d,x[1]) + sum( log(x[2]+x[3]*vec1) ) + sum( log(1-x[2]+x[3]*vec2) ) - sum( log(1+x[3]*vec3) )
}
