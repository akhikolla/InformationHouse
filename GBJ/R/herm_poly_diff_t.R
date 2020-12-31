#' herm_poly_diff_t.R
#'
#' Internal function to add up the infinite sum term in calculating the variance of S(t) when
#' the test statistics Z_1,...,Z_d have a non-zero mean.  See GBJ paper for more details.
#' 
#' @param t1 t1 and t2 are interchangeable.  One should be t-mu and one should be -t-mu.
#' @param t2 t1 and t2 are interchangeable.  One should be t-mu and one should be -t-mu.
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test 
#' statistics, where d is total number of test statistics in the set.
#' 
#' @return The result of the infinite sum.
#'
#' @keywords internal
#' @export
#'
#' @examples 
#' herm_poly_diff_t(t1=0, t2=0, pairwise_cors=rep(0.3, 10))

# The infinite sum terms when we have to multiply two different Hermite poly terms
# The 't' here refers to t-mu or -t-mu or whatever is inside the H
herm_poly_diff_t <- function(t1, t2, pairwise_cors) {
	# The r_bar terms
	rho_bar1 <- mean(pairwise_cors)
	rho_bar2 <- mean(pairwise_cors^2)
	rho_bar3 <- mean(pairwise_cors^3)
	rho_bar4 <- mean(pairwise_cors^4)
	rho_bar5 <- mean(pairwise_cors^5)
	rho_bar6 <- mean(pairwise_cors^6)
	rho_bar7 <- mean(pairwise_cors^7)
	rho_bar8 <- mean(pairwise_cors^8)
	rho_bar9 <- mean(pairwise_cors^9)
	rho_bar10 <- mean(pairwise_cors^10)
	
	He1 <- (t1) * (t2)
	He3 <- (t1^3 - 3*t1) * (t2^3 - 3*t2)
	He5 <- (t1^5 - 10*t1^3 + 15*t1) * (t2^5 - 10*t2^3 + 15*t2)	
	He7 <- (t1^7 - 21*t1^5 + 105*t1^3 - 105*t1) * (t2^7 - 21*t2^5 + 105*t2^3 - 105*t2)
	He9 <- (t1^9 - 36*t1^7 + 378*t1^5 - 1260*t1^3 + 945*t1) * (t2^9 - 36*t2^7 + 378*t2^5 - 1260*t2^3 + 945*t2)
	He0 <- 1
	He2 <- (t1^2 - 1) * (t2^2 - 1)
	He4 <- (t1^4 - 6*t1^2 + 3) * (t2^4 - 6*t2^2 + 3)
	He6 <- (t1^6 - 15*t1^4 + 45*t1^2 - 15) * (t2^6 - 15*t2^4 + 45*t2^2 - 15)
	He8 <- (t1^8 - 28*t1^6 + 210*t1^4 - 420*t1^2 + 105) * (t2^8 - 28*t2^6 + 210*t2^4 - 420*t2^2 + 105)
	odds <- ( He1*rho_bar2/2 + He3*rho_bar4/24 + He5*rho_bar6/720 + He7*rho_bar8/40320 + He9*rho_bar10/3628800 )
	evens <- ( He0*rho_bar1/1 + He2*rho_bar3/6 + He4*rho_bar5/120 + He6*rho_bar7/5040 + He8*rho_bar9/362880 )
	sum_term <- odds + evens
	return(sum_term)
}
