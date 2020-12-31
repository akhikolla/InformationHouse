#' GBJ.R
#'
#' Calculate the Generalized Berk-Jones test statistic and p-value.
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
#' \item{GBJ}{The observed Generalized Higher Criticism test statistic.}
#' \item{GBJ_pvalue}{The p-value of this observed value, given the size of the set and
#' correlation structure.}
#' \item{err_code}{Sometimes if your p-value is very small (<10^(-12) usually), R/C++ do not
#' have enough precision in their standard routines to calculate the number accurately. In
#' these cases (and very rarely others) we switch to standard Berk-Jones instead (more stable
#' numerically) and let you know with a message here.}
#'
#'
#' @import stats BH
#' @importFrom Rcpp evalCpp
#' @useDynLib GBJ
#'
#' @export
#' @examples
#' # Should return statistic = 0.9248399 and p_value = 0.2670707
#' set.seed(100)
#' Z_vec <- rnorm(5) + rep(1,5)
#' cor_Z <- matrix(data=0.2, nrow=5, ncol=5)
#' diag(cor_Z) <- 1
#' GBJ(test_stats=Z_vec, cor_mat=cor_Z)


GBJ <- function(test_stats, cor_mat=NULL, pairwise_cors=NULL)
{
  # Parse inputs, do some error checking.
  param_list <- parse_input(test_stats=test_stats, cor_mat=cor_mat,
                            pairwise_cors=pairwise_cors)
  t_vec <- param_list$t_vec
  pairwise_cors <- param_list$pairwise_cors
  d <- length(t_vec)

  # Move to BJ if no correlation at all
  if (sum(abs(pairwise_cors)) == 0) {
    BJ_output <- BJ(test_stats=t_vec, pairwise_cors=pairwise_cors)
    return ( list(GBJ=BJ_output$BJ, GBJ_pvalue=BJ_output$BJ_pvalue, err_code=0) )
  }

	# Calculate the observed GBJ statistic
	GBJ_stats <- GBJ_objective(t_vec=t_vec, d=d, pairwise_cors=pairwise_cors)
	gbj <- max(GBJ_stats)

	# Calculate p_value
	GBJ_p_list <- GBJ_pvalue(observed_gbj=gbj, d=d, pairwise_cors=pairwise_cors)
	GBJ_corp=GBJ_p_list$GBJ_corp

	# If we get NA for the p-value, don't want to return NA to the user
	if (is.na(GBJ_corp)) {

	  # Give BJ as replacement
	  BJ_output <- BJ(test_stats=t_vec, pairwise_cors=pairwise_cors)

	  # If gbj >= 20, give them the BJ p-value with a disclaimer
	  if (gbj >= 20) {
	    GBJ_err_code <- '1: Pvalue likely less than 10^(-12), R/C++ not enough precision. Returning standard Berk-Jones test instead.'
	    return ( list(GBJ=BJ_output$BJ, GBJ_pvalue=BJ_output$BJ_pvalue, err_code=GBJ_err_code) )
	  }

	  # If evidence of underdispersion, again give them BJ p-value
	  else if (sum(pairwise_cors) < 0) {
	    GBJ_err_code <- '2: Error in numerical routines. Many apologies, please report to developer! Returning standard Berk-Jones test instead.'
	    return ( list(GBJ=BJ_output$BJ, GBJ_pvalue=BJ_output$BJ_pvalue, err_code=GBJ_err_code) )
	  }

	  # Any other errors, give them BJ p-value
	  else {
	    GBJ_err_code <- '3: Unknown error. Many apologies, please report to developer! Returning standard Berk-Jones test instead.'
	    return ( list(GBJ=BJ_output$BJ, GBJ_pvalue=BJ_output$BJ_pvalue, err_code=GBJ_err_code) )
	  }
	}

	# If we have a negative p-value, don't wnat to return that to the user either
	if (GBJ_corp < 0) {

	  # Give BJ as replacement
	  BJ_output <- BJ(test_stats=t_vec, pairwise_cors=pairwise_cors)

	  # If gbj >= 20, give them the BJ p-value with a disclaimer
	  if (gbj >= 20) {
	    GBJ_err_code <- '1: Pvalue likely less than 10^(-12), R/C++ not enough precision. Returning standard Berk-Jones test instead.'
	    return ( list(GBJ=BJ_output$BJ, GBJ_pvalue=BJ_output$BJ_pvalue, err_code=GBJ_err_code) )
	  } else {
	    GBJ_err_code <- '3: Unknown error. Many apologies, please report to developer! Returning standard Berk-Jones test instead.'
	    return ( list(GBJ=BJ_output$BJ, GBJ_pvalue=BJ_output$BJ_pvalue, err_code=GBJ_err_code) )
	  }
	}

	# If no NA in the p-value, then everything success and return GBJ output.
	return ( list(GBJ=gbj, GBJ_pvalue=GBJ_corp, err_code=0) )
}

