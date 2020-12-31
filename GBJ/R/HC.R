#' HC.R
#'
#' Calculate the Higher Criticism test statistic and p-value.
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
#' \item{HC}{The observed Higher Criticism test statistic.}
#' \item{HC_pvalue}{The p-value of this observed value, given the size of the set and
#' correlation structure.}
#'
#' @export
#' @examples
#' # Should return statistic = 2.067475 and p_value = 0.2755146
#' set.seed(100)
#' Z_vec <- rnorm(5) + rep(1,5)
#' cor_Z <- matrix(data=0.2, nrow=5, ncol=5)
#' diag(cor_Z) <- 1
#' HC(test_stats=Z_vec, cor_mat=cor_Z)

HC <- function(test_stats, cor_mat=NULL, pairwise_cors=NULL) {

  # Parse inputs, do some error checking.
  param_list <- parse_input(test_stats=test_stats, cor_mat=cor_mat,
                            pairwise_cors=pairwise_cors)
  t_vec <- param_list$t_vec
  pairwise_cors <- param_list$pairwise_cors
  d <- length(t_vec)

	# Calculate HC objectives
	p_values <- 1-pchisq(t_vec^2, df=1)
	i_vec <- 1:d
	HC_stats <- sqrt(d) * (i_vec/d - p_values) / sqrt(p_values*(1-p_values))

	# Observed HC statistic
	h <- max(HC_stats, na.rm=TRUE)

	# Calculate p-value
	if (h<=0) {
		return ( list(HC=0, HC_pvalue=1) )
	}

	# BJ bounds
	HC_p_bounds <- rep(NA, d)

	# Explicit inverse of HC to find the p-value bounds
	HC_p_bounds <- ((2*i_vec+h^2)/d - sqrt((2*i_vec/d+h^2/d)^2 - 4*i_vec^2/d^2 - 4*i_vec^2*h^2/d^3))/(2*(1+h^2/d))
	HC_z_bounds <- qnorm(1-HC_p_bounds/2)
	HC_z_bounds <- sort(HC_z_bounds, decreasing=F)

	# qnorm can't handle more precision than 10^-16
	# Also crossprob_cor can only handle Z up to 8.2
	HC_z_bounds[which(HC_z_bounds > 8.2)]= 8.2

	# Send it to the C++.
	if (sum(abs(pairwise_cors)) == 0) {
		# For the independence flag in the c++, just have to send a number < -1.
		HC_corp <- ebb_crossprob_cor_R(d=d, bounds=HC_z_bounds, correlations=rep(-999,2))
	} else {
		HC_corp <- ebb_crossprob_cor_R(d=d, bounds=HC_z_bounds, correlations=pairwise_cors)
	}


	return ( list(HC=h, HC_pvalue=HC_corp) )
}
