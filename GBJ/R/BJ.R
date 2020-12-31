#' BJ.R
#'
#' Calculate the Berk-Jones test statistic and p-value.
#'
#' @param test_stats Vector of test statistics for each factor in the set (i.e. marginal
#' test statistic for each SNP in a gene).
#' @param cor_mat d*d matrix of the correlations between all the test statistics in
#' the set, where d is the total number of test statistics in the set.
#' You only need to specify EITHER cor_mat OR pairwise_cors.
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test
#' statistics. You only need to specify EITHER cor_mat OR pairwise_cors.
#'
#' @return A list with the elements:
#' \item{BJ}{The observed Berk-Jones test statistic.}
#' \item{BJ_pvalue}{The p-value of this observed value, given the size of the set and
#' correlation structure.}
#'
#' @export
#' @examples
#' # Should return statistic = 1.243353 and p_value = 0.256618
#' set.seed(100)
#' Z_vec <- rnorm(5) + rep(1,5)
#' cor_Z <- matrix(data=0.2, nrow=5, ncol=5)
#' diag(cor_Z) <- 1
#' BJ(test_stats=Z_vec, cor_mat=cor_Z)

BJ <- function(test_stats, cor_mat=NULL, pairwise_cors=NULL) {

	# Parse inputs, do some error checking.
  param_list <- parse_input(test_stats=test_stats, cor_mat=cor_mat,
                            pairwise_cors=pairwise_cors)
  t_vec <- param_list$t_vec
  pairwise_cors <- param_list$pairwise_cors
  d <- length(t_vec)

	# Check for qualifying p-values under the null (the indicator part of the GBJ statistic)
	# and also that we are only considering 'first half' p-values
	p_values <- 1-pchisq(t_vec^2, df=1)
	BJ_indicator <- which( p_values < (1:d)/d )
	first_half <- 1:(ceiling(d/2))
	non_zero <- intersect(BJ_indicator, first_half)

	# If no indicies qualified, stop
	if (length(non_zero) == 0) {
		return ( list(BJ=0, BJ_pvalue=1) )
	}

	#################################
	# Some indicies passed
	i_vec <- 1:d
	BJ_stats <- rep(0, d)
	BJ_stats[non_zero] <- i_vec[non_zero] * log(i_vec[non_zero]/(d*p_values[non_zero])) +
								(d-i_vec[non_zero]) * log((1-i_vec[non_zero]/d)/(1-p_values[non_zero]))
	BJ_stats[d] <- 0

	# Observed BJ statistic
	b <- max(BJ_stats[1:(ceiling(d/2))])

	# Calculate p-value
	if (b<=0) {
		return ( list(BJ=0, BJ_pvalue=1) )
	}

	# BJ bounds
	BJ_p_bounds <- rep(NA, d)

	# Use uniroot to find the pvalue bounds.
	for ( jjj in 1:(ceiling(d/2)) ) {
		BJ_p_bounds[jjj] <- uniroot(f=function(x, k, d, b) {k*log(k/(d*x)) +
										(d-k)*log((1-k/d)/(1-x)) - b}, k=jjj, d=d,
										b=b, lower=0, upper=jjj/d, tol=(10^(-12)))$root
	}

	# The last half of the order statistic bounds
	BJ_p_bounds[(ceiling(d/2)+1):d] <- BJ_p_bounds[ceiling(d/2)]

	# Now put the bounds in terms of the Z statistics
	BJ_z_bounds <- qnorm(1 - BJ_p_bounds/2)
	BJ_z_bounds <- sort(BJ_z_bounds, decreasing=F)

	# qnorm can't handle more precision than 10^-16
	# Also crossprob_cor can only handle Z up to 8.2
	BJ_z_bounds[which(BJ_z_bounds > 8.2)]= 8.2

	# Send it to the C++.
	if (sum(abs(pairwise_cors)) == 0) {
		# For the independence flag in the c++, just have to send a number < -1.
		BJ_corp <- ebb_crossprob_cor_R(d=d, bounds=BJ_z_bounds, correlations=rep(-999,2))
	} else {
		BJ_corp <- ebb_crossprob_cor_R(d=d, bounds=BJ_z_bounds, correlations=pairwise_cors)
	}


	return ( list(BJ=b, BJ_pvalue=BJ_corp) )
}
