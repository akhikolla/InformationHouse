#' GHC_objective.R
#'
#' Internal function to calculate the value of the kth GHC objective (possibly minus an offset)
#' when the kth p-value order statistic is x, for a set of size d.
#' 
#' @param x The p-value of the kth p-value order statistic.
#' @param k Which objective to use.
#' @param d The size of the set.
#' @param offset Used to zero the correct value when we put this into uniroot
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test 
#' statistics, where d is total number of test statistics in the set.
#'
#' @return The value of the kth HC objective 
#'
#' @keywords internal
#' @export
#'
#' @examples 
#' GHC_objective(x=0.1, k=2, d=5, offset=0, pairwise_cors=rep(0.2,10))


GHC_objective <- function(x, k, d, offset, pairwise_cors) 
{
	# We get errors in qnorm once x is too small
	if(x < 10^(-16)) {
		temp_Z <- qnorm(1-10^(-15)/2)
	} else {
		temp_Z <- qnorm(1-x/2)
	}
	(k - d*x) / sqrt(calc_var_nonzero_mu(d=d, t=temp_Z, mu=0, pairwise_cors=pairwise_cors)) - offset
}

