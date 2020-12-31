#' minP.R
#'
#' Given a vector of individual test statistics and their pairwise correlations, calculate
#' the MinimumP (see Conneely and Boehnke, 2007) second-level test statistic and it's p-value.
#'
#' @param test_stats Vector of test statistics for each factor in the set (i.e. marginal
#' test statistic for each SNP in a gene)
#' @param cor_mat d*d matrix of the correlations between all the test statistics in
#' the set, where d is the total number of test statistics in the set.
#' You only need to specify EITHER cor_mat OR pairwise_cors.
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test
#' statistics. You only need to specify EITHER cor_mat OR pairwise_cors.
#'
#' @return A list with the elements:
#' \item{minP}{The observed MinimumP test statistic.}
#' \item{minP_pvalue}{The p-value of this observed value, given the size of the set and
#' correlation structure.}
#'
#' @export
#' @examples
#' # Should return statistic = 0.05918928 and p_value = 0.2525972.
#' set.seed(100)
#' Z_vec <- rnorm(5) + rep(1,5)
#' cor_Z <- matrix(data=0.2, nrow=5, ncol=5)
#' diag(cor_Z) <- 1
#' minP(test_stats=Z_vec, cor_mat=cor_Z)


minP <- function(test_stats, cor_mat=NULL, pairwise_cors=NULL)
{
  # Parse inputs, do some error checking.
  param_list <- parse_input(test_stats=test_stats, cor_mat=cor_mat,
                            pairwise_cors=pairwise_cors)
  t_vec <- param_list$t_vec
  pairwise_cors <- param_list$pairwise_cors
  d <- length(t_vec)

	# minP bounds
  minP_stat = 1-pchisq(t_vec[1]^2, df=1)
	minP_p_bounds <- rep(minP_stat, d)
	minP_z_bounds <- qnorm(1-minP_p_bounds/2)
	minP_z_bounds <- sort(minP_z_bounds, decreasing=F)

	minP_z_bounds[which(minP_z_bounds > 8.2)]= 8.2

	# Send it to the C++.
	if (sum(abs(pairwise_cors)) == 0) {
		# For the independence flag in the c++, just have to send a number < -1.
		minP_corp <- ebb_crossprob_cor_R(d=d, bounds=minP_z_bounds, correlations=rep(-999,2))
	} else {
		minP_corp <- ebb_crossprob_cor_R(d=d, bounds=minP_z_bounds, correlations=pairwise_cors)
	}

	return ( list(minP=minP_stat, minP_pvalue=minP_corp) )
}

