#' calc_var_mu_nonzero.R
#'
#' Internal function to calculate variance of S(t) when the Z_1,...,Z_d have nonzero mean.  
#' See GBJ paper for more details.  Vectorized in t.
#' 
#' @param d The number of test statistics in the set.
#' @param t The threshold which determines the properties of S(t).
#' @param mu The common mean of the test statistics Z_1,...,Z_d
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test 
#' statistics, where d is total number of test statistics in the set.
#' 
#' @keywords internal
#' @export
#'
#' @return The variance of S(t) when the Z_1,...,Z_d have nonzero mean mu.
#'
#' @examples 
#' calc_var_nonzero_mu(d=5, t=1, mu=1, pairwise_cors=rep(0.3, 10))

calc_var_nonzero_mu <- function(d, t, mu, pairwise_cors) {
	
	# The variance of S(t) under independence
	# Remember that P(|Z_i|>=t) is no longer 2*surv(t-mu)
	prob_greater <- 1 - ( pnorm(t, mean=mu, lower.tail=TRUE) - pnorm(-t, mean=mu, lower.tail=TRUE) )
	ind_term <- d*prob_greater - d*prob_greater^2
	
	# The 4 terms that make up the 'covariance term'
	# P(|Z_j|,|Z_k|>=t) = P(Z_j,Z_k>=t) + P(Z_j,Z_k<=-t) + 2*P(Z_j>=t,Z_k<=-t)
	cov_term_gg <- d*(d-1)* ( surv(t-mu)^2 + dnorm(t-mu)^2*herm_poly_diff_t(t1=(t-mu), t2=(t-mu), 
		pairwise_cors=pairwise_cors) ) 
	cov_term_ll <- d*(d-1)* ( 1-2*surv(-t-mu)+surv(-t-mu)^2 + dnorm(-t-mu)^2*herm_poly_diff_t(t1=(-t-mu), 
		t2=(-t-mu), pairwise_cors=pairwise_cors) )
	cov_term_diff <- d*(d-1)* ( surv(t-mu) - surv(t-mu)*surv(-t-mu) - dnorm(t-mu)*dnorm(-t-mu)*herm_poly_diff_t(t1=(t-mu),t2=(-t-mu),pairwise_cors=pairwise_cors) )
	
	# Remember to subtract the E[X1]E[X2] terms in the covariance term calculation!
	cov_term <- cov_term_gg + cov_term_ll + 2*cov_term_diff - d*(d-1)*prob_greater^2
	
	return(cov_term + ind_term)
}


