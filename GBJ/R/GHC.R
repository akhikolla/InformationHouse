#' GHC.R
#'
#' Calculate the Generalized Higher Criticism test statistic and p-value.
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
#' \item{GHC}{The observed Generalized Higher Criticism test statistic.}
#' \item{GHC_pvalue}{The p-value of this observed value, given the size of the set and
#' correlation structure.}
#'  \item{err_code}{Sometimes if your p-value is very small (<10^(-12) usually), R/C++ do not
#' have enough precision in their standard routines to calculate the number accurately. In
#' these cases (and very rarely others) we switch to standard Higher Criticism instead (more stable
#' numerically) and let you know with a message here.}
#'
#' @export
#' @examples
#' set.seed(100)
#' Z_vec <- rnorm(5)
#' cor_Z <- matrix(data=0.2, nrow=5, ncol=5)
#' diag(cor_Z) <- 1
#' GHC(test_stats=Z_vec, cor_mat=cor_Z)


GHC <- function(test_stats, cor_mat=NULL, pairwise_cors=NULL)
{
  # Parse inputs, do some error checking.
  param_list <- parse_input(test_stats=test_stats, cor_mat=cor_mat,
                            pairwise_cors=pairwise_cors)
  t_vec <- param_list$t_vec
  pairwise_cors <- param_list$pairwise_cors
  d <- length(t_vec)

  # Move to hC if no correlation at all
  if (sum(abs(pairwise_cors)) == 0) {
    HC_output <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
    return ( list(GHC=HC_output$HC, GHC_pvalue=HC_output$HC_pvalue, err_code=0) )
  }

	# Calculate GHC objectives
	i_vec <- 1:d
	p_values <- 1-pchisq(t_vec^2, df=1)
	GHC_stats <- (i_vec - d*p_values) / sqrt(calc_var_nonzero_mu(d=d, t=t_vec, mu=0,
				 pairwise_cors=pairwise_cors))

	# Observed GHC statistic - sometimes a Z-statistic is 0 and so we get NA for variance
	ghc <- max(GHC_stats, na.rm=TRUE)

	# Calculate p-value
	if (ghc <= 0) {
		return (list(GHC=ghc, GHC_pvalue=1, err_code=0))
	}

	# GHC bounds
	GHC_p_bounds <- rep(NA, d)

	# increase tolerance of uniroot for large ghc
	if(ghc>10) {
		my_tol <- (-100)
	} else {my_tol <- (-12)}

	# Use uniroot to find the pvalue bounds.
	GHC_lowerbound <- 10^(-20)
	for(kkk in 1:d) {
		# Sometimes run into precision errors, need to think about how
		# we can fix this so don't need the tryCatch
		temp_ghc <- tryCatch(uniroot(GHC_objective, k=kkk, d=d, offset=ghc,
				pairwise_cors=pairwise_cors, lower=GHC_lowerbound, upper=(1-10^(-12)),
				tol=(10^(my_tol))), error=function(e) e, warning=function(w) w)

		# If it doesn't work, just run the HC for them
		if(length(class(temp_ghc))>1) {
		  if (ghc >= 200000) {
		    HC_output <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
		    GHC_err_code <- '1: Pvalue likely less than 10^(-12), R/C++ not enough precision. Returning standard Higher Criticism test instead.'
		    return ( list(GHC=HC_output$HC, GHC_pvalue=HC_output$HC_pvalue, err_code=GHC_err_code) )
		  }

		  # If evidence of underdispersion, again give them BJ p-value
		  else if (sum(pairwise_cors) < 0) {
		    HC_output <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
		    GHC_err_code <- '2: Error in numerical routines. Many apologies, please report to developer! Returning standard Higher Criticism test instead.'
		    return ( list(GHC=HC_output$HC, GHC_pvalue=HC_output$HC_pvalue, err_code=GHC_err_code) )
		  }

		  # Any other errors, give them BJ p-value
		  else {
		    HC_output <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
		    GHC_err_code <- '3: Unknown error. Many apologies, please report to developer! Returning standard Higher Criticism test instead.'
		    return ( list(GHC=HC_output$HC, GHC_pvalue=HC_output$HC_pvalue, err_code=GHC_err_code) )
		  }
		}

		# It worked, keep going
		GHC_p_bounds[kkk] <- temp_ghc$root

		# small security measure to ensure that GHC bounds are increasing
		GHC_lowerbound <- GHC_p_bounds[kkk]
	}

	# now put the bounds in terms of the Z statistics
	GHC_z_bounds <- qnorm(1-GHC_p_bounds/2)
	GHC_z_bounds <- sort(GHC_z_bounds, decreasing=F)

	# qnorm can't handle more precision than 10^-16
	# Also crossprob_cor can only handle Z up to 8.2
	GHC_z_bounds[which(GHC_z_bounds > 8.2)]= 8.2

	# Send it to the C++.
	if (sum(abs(pairwise_cors)) == 0) {
		# For the independence flag in the c++, just have to send a number < -1.
		GHC_corp <- ebb_crossprob_cor_R(d=d, bounds=GHC_z_bounds, correlations=rep(-999,2))
	} else {
		GHC_corp <- ebb_crossprob_cor_R(d=d, bounds=GHC_z_bounds, correlations=pairwise_cors)
	}

	# If we get NA for the p-value, don't want to return NA to the user
	if (is.na(GHC_corp)) {
	  # If gbj >= 250,000, give them the HC p-value with a disclaimer
	  if (ghc >= 20) {
	    HC_output <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
	    GHC_err_code <- '1: Pvalue likely less than 10^(-12), R/C++ not enough precision. Returning standard Higher Criticism test instead.'
	    return ( list(GHC=HC_output$HC, GHC_pvalue=HC_output$HC_pvalue, err_code=GHC_err_code) )
	  }

	  # If evidence of underdispersion, again give them HC p-value
	  else if (sum(pairwise_cors) < 0) {
	    HC_output <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
	    GHC_err_code <- '2: Error in numerical routines. Many apologies, please report to developer! Returning standard Higher Criticism test instead.'
	    return ( list(GHC=HC_output$HC, GHC_pvalue=HC_output$HC_pvalue, err_code=GHC_err_code) )
	  }

	  # Any other errors, give them HC p-value
	  else {
	    HC_output <- HC(test_stats=t_vec, pairwise_cors=pairwise_cors)
	    GHC_err_code <- '3: Unknown error. Many apologies, please report to developer! Returning standard Higher Criticism test instead.'
	    return ( list(GHC=HC_output$HC, GHC_pvalue=HC_output$HC_pvalue, err_code=GHC_err_code) )
	  }
	}


	return ( list(GHC=ghc, GHC_pvalue=GHC_corp, err_code=0) )
}

