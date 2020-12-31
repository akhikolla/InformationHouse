#' GBJ_pvalue.R
#'
#' Calculate the p-value for the Generalized Berk-Jones (GBJ) statistic.
#'
#' @param observed_gbj The observed value of the GBJ statistic.
#' @param d The number of test statistics in the set.
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test
#' statistics, where d is total number of test statistics in the set.
#' @param times_to_try Sometimes the numerical root-finder is finnicky, so we have to give
#' it extra chances to try and calculate the p-value if first time is failure.  Recommend
#' setting this parameter to 5.
#'
#' @import stats
#' @return The p-value of the GBJ test.
#'
#' @export
#' @examples
#' GBJ_pvalue(observed_gbj=2, d=5, pairwise_cors=rep(0.2,10))


GBJ_pvalue <- function(observed_gbj, d, pairwise_cors, times_to_try=5)
{
	if (observed_gbj<=0)
	{
		return ( list(GBJ_corp=1, eFlag=0) )
	} else
	{
		# Try to rerun multiple times if fails at first.
		times_tried <- 0
		eFlag <- 1
		while (times_tried<times_to_try & eFlag==1)
		{
			# Error flag, if we don't manage to throw this, then break out of while loop
			eFlag <- 0

			# Make gBJ_BB bounds
			GBJ_z_bounds <- rep(NA,d)

			# Sometimes starting too high gives an error
			if (times_tried == 2) {
			  prev_bound <- 5
			} else if (times_tried == 3) {
			  prev_bound <- 3
			} else if (times_tried == 4) {
			  prev_bound <- 1
			} else {
			  prev_bound <- 8.2
			}

			for ( kkk in 1:(ceiling(d/2)) )
			{
				temp_gbj <- tryCatch(uniroot(GBJ_objective, interval=c(0, prev_bound), d=d, k_vec=kkk, pairwise_cors=pairwise_cors, offset=observed_gbj), error=function(e) e, warning=function(w) w)

				# Sometimes, we can't go high enough in t, because pnorm,etc will just round to 0, and thus
				# the signs on both sides of the interval will be the same.  In this case, we will try again
				# a few times and then give up.
				if(length(class(temp_gbj))>1) {
					eFlag <- 1
					times_tried <- times_tried + 1
					break
				} else {
					GBJ_z_bounds[kkk] <- temp_gbj$root
				}
				prev_bound <- GBJ_z_bounds[kkk]
			}
		}

		# If eFlag still 1, then we tried multiple times times with no success.
		if(eFlag==1)
		{
			return ( list(GBJ_corp=NA, eFlag=eFlag) )
		}

		# Only need the first half
		# Make sure to sort in increasing order for crossprob_cor
		GBJ_z_bounds[(ceiling(d/2)+1):d] <- GBJ_z_bounds[ceiling(d/2)]
		GBJ_z_bounds <- sort(GBJ_z_bounds, decreasing=FALSE)

		# qnorm can't handle more precision than 10^-16
		# Also crossprob_cor can only handle Z up to 8.2
		GBJ_z_bounds[which(GBJ_z_bounds > 8.2)]= 8.2

		# Send it to the C++.
		if (sum(abs(pairwise_cors)) == 0) {
			# For the independence flag in the c++, just have to send a number < -1.
			GBJ_corp <- ebb_crossprob_cor_R(d=d, bounds=GBJ_z_bounds, correlations=rep(-999,2))
		} else {
			GBJ_corp <- ebb_crossprob_cor_R(d=d, bounds=GBJ_z_bounds, correlations=pairwise_cors)
		}

		return ( list(GBJ_corp=GBJ_corp, eFlag=eFlag) )
	}
}

