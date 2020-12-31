#' GBJ_objective.R
#'
#' Calculates the GBJ objective function at given threshold points.  Used both to calculate the
#' observed GBJ statistic and also to find the boundary points for p-value calculation (through uniroot).
#' Mostly for internal use.
#'
#' @param t_vec A scalar or vector of threshold points (magnitudes of test statistics)
#' to calculate the objective at.
#' @param d The total number of test statistics in the set.
#' @param k_vec If t_vec is not a vector all of the test statistics, let us know which ordered
#' objective functions we are calculating (can be a vector with the same length as t).
#' @param pairwise_cors A vector of all d(d-1)/2 pairwise correlations between the test
#' statistics, where d is total number of test statistics in the set.
#' @param offset Used if we want to subtract the observed gbj value for uniroot.
#'
#' @return The GBJ objective at t (for given d, kkk, pairwise_cors) minus the offset
#'
#' @keywords internal
#' @export
#'
#' @examples
#' GBJ_objective(t_vec=1:5, d=5, k_vec=NULL, pairwise_cors=rep(0.2,10), offset=0)

GBJ_objective <- function(t_vec, d, k_vec=NULL, pairwise_cors, offset=0) {

	# Ensure that the thresholds are sorted in descending order, largest first.
	t_vec <- sort(abs(t_vec), decreasing=TRUE)

	# If didn't pass in a kkk vector, create it now (assume passed in all test statistics)
	if (is.null(k_vec)) {
		k_vec <- 1:d
	}

	# Check for qualifying p-values under the null (the indicator part of the GBJ statistic)
	# and also that we are only considering 'first half' p-values
	p_values <- 1-pchisq(t_vec^2, df=1)
	GBJ_indicator <- which( p_values < k_vec/d )
	first_half <- which(k_vec <= ceiling(d/2))
	non_zero <- intersect(GBJ_indicator, first_half)

	# If no indicies qualified, stop
	if (length(non_zero) == 0) {
		return (rep(0,length(t_vec)) - rep(offset, length(t_vec)) )
	}

	#################################
	# Some indices passed
	# Calculate mean and variance of S(t) under the null for non-zero t
	mean_null <- rep(NA, length(t_vec))
	sigsq_null <- rep(NA, length(t_vec))
	mean_null[non_zero] <- 2*d*surv(t_vec[non_zero])
	sigsq_null <- calc_var_nonzero_mu(d=d, t=t_vec, mu=rep(0,length(t_vec)), pairwise_cors=pairwise_cors)

	# Now match moments for null (directly)
	mu_null <- rep(NA, length(t_vec))
	rho_null <- rep(NA, length(t_vec))
	gamma_null <- rep(NA, length(t_vec))
	mu_null[non_zero] <- mean_null[non_zero] / d
	rho_null[non_zero] <- (sigsq_null[non_zero] - d*mu_null[non_zero]*(1-mu_null[non_zero])) /
			(d*(d-1)*mu_null[non_zero]*(1-mu_null[non_zero]))
	gamma_null[non_zero] = rho_null[non_zero] / (1-rho_null[non_zero])


	################################
	# Calculate variance of S(t) under the alternative for non-zero t (the mean is given by k_vec)
	# First, for each t, need to find the fitted common mean of the underlying test statistics
	Z_common_means <- rep(NA, length(t_vec))
	for (iii in 1:length(t_vec)) {
		if(iii %in% non_zero)
		{
			# First get us the mean of the Z stats under the alternative
			Z_common_means[iii] <- uniroot(qnorm_mu, lower=0, upper=100, t=t_vec[iii], kkk=k_vec[iii], d=d)$root
		}
	}

	# Now we can calculate the variance
	alt_mean_vec <- k_vec
	alt_var_vec <- rep(NA, length(t_vec))
	alt_var_vec[non_zero] <- calc_var_nonzero_mu(d=d, t=t_vec[non_zero], mu=Z_common_means[non_zero],
		pairwise_cors=pairwise_cors)

	# Now match moments for alternative
	mu_alt <- rep(NA, length(t_vec))
	rho_alt <- rep(NA, length(t_vec))
	gamma_alt <- rep(NA, length(t_vec))
	mu_alt[non_zero] <- alt_mean_vec[non_zero] / d
	rho_alt[non_zero] <- (alt_var_vec[non_zero] - d*mu_alt[non_zero]*(1-mu_alt[non_zero])) /
		 (d*(d-1)*mu_alt[non_zero]*(1-mu_alt[non_zero]))
	gamma_alt[non_zero] = rho_alt[non_zero] / (1-rho_alt[non_zero])

	# Numerical issues may give us mu and gamma that don't respect the bounds.
	# We need gamma >= max{ -mu/(d-1), -(1-mu)/(d-1) }
	pq_mat_null <- matrix(data=NA, nrow=length(t_vec), ncol=2)
	pq_mat_alt <- matrix(data=NA, nrow=length(t_vec), ncol=2)
	pq_mat_null[non_zero, 1] <- -mu_null[non_zero] / (d-1)
	pq_mat_null[non_zero, 2] <- -(1-mu_null[non_zero]) / (d-1)
	gamma_check_null <- apply(pq_mat_null, 1, max)
	pq_mat_alt[non_zero, 1] <- -mu_alt[non_zero] / (d-1)
	pq_mat_alt[non_zero, 2] <- -(1-mu_alt[non_zero]) / (d-1)
	gamma_check_alt <- apply(pq_mat_alt, 1, max)

	# Recalibrate non_zero
	non_zero <- which(gamma_null >= gamma_check_null &
					gamma_alt >= gamma_check_alt)
	if (length(non_zero)==0 ) {
		return ( rep(0,length(t_vec)) - rep(offset, length(t_vec)) )
	}

	# Log-liklihood for null
	null_loglik <- rep(0, length(t_vec))
	null_param_mat <- cbind( k_vec[non_zero], mu_null[non_zero], gamma_null[non_zero])
	null_loglik[non_zero] <- apply(null_param_mat, 1, ebb_loglik, d=d)

	# Log-likelihood for alternative
	alt_loglik <- rep(0, length(t_vec))
	alt_param_mat <- cbind( k_vec[non_zero], mu_alt[non_zero], gamma_alt[non_zero])
	alt_loglik[non_zero] <- apply(alt_param_mat, 1, ebb_loglik, d=d)

	BB_GBJ_stats <- alt_loglik - null_loglik

	return(BB_GBJ_stats - rep(offset, length(t_vec)))
}
