#' calc_score_stats.R
#'
#' Starting with individual-level data on p factors, generate score test statistics for each
#' factor for input into GBJ/GHC/HC/BJ/minP.  Also get the correlations between these test statistics.
#' Designed to be used with linear or logistic or log-linear regression null models.
#'
#' @param null_model An R regression model fitted using glm().  Do not use lm(), even for linear regression!
#' @param factor_matrix An n*p matrix with each factor as one column.  There should be no missing data.
#' @param link_function Either "linear" or "logit" or "log"
#' @param P_mat The projection matrix used in calculation may be passed in to speed up the calculation.
#' See paper for details. Default is null.
#'
#' @return A list with the elements:
#' \item{test_stats}{The p score test statistics.}
#' \item{cor_mat}{The p*p matrix giving the pairwise correlation of every two test statistics.}
#'
#' @export
#' @examples
#' set.seed(0)
#' Y <- rbinom(n=100, size=1, prob=0.5)
#' null_mod <- glm(Y~1, family=binomial(link="logit"))
#' factor_mat <- matrix(data=rnorm(n=100*5), nrow=100)
#' calc_score_stats(null_mod, factor_mat, "logit")

calc_score_stats <- function(null_model, factor_matrix, link_function, P_mat=NULL) {

	X_mat <- model.matrix(null_model)
	d <- ncol(factor_matrix)
	fitted_Y <- null_model$fitted.values
	actual_Y <- null_model$y

	# Only difference between linear and logistic procedure
	if (link_function == 'logit') {
		W_vec <- fitted_Y * (1-fitted_Y)
	} else if (link_function == 'linear') {
		W_vec <- rep(summary(null_model)$dispersion, nrow(X_mat))
	} else if (link_function == 'log') {
		W_vec <- fitted_Y
	} else {
		stop("Invalid model type")
	}

	########################
	# EZ Mode if linear regression, no additional covariates except for intercept
	if (link_function == 'linear' & ncol(model.matrix(null_model)) == 1) {
		num_sub <- nrow(X_mat)
		sig_sq_hat <- sum( (actual_Y - fitted_Y)^2 ) / (num_sub-1)
		test_stats <- rep(NA, d)
		denominators <- rep(NA, d)
		for(kkk in 1:d)
		{
			tempF<- factor_matrix[,kkk]
			score_num <- t(tempF) %*% (actual_Y-fitted_Y)
			score_denom <- tryCatch(sqrt(sig_sq_hat * (sum(tempF^2) - mean(tempF)^2*num_sub)),
			                        warning=function(w) w, error=function(e) e)

			# We've been getting negative denominators with, for example, very rare SNPs
			if (!(class(score_denom)[1] %in% c("matrix", "array", "numeric"))) {
			  err_msg <- paste('Error in calculating test statistic for factor ', kkk,
			                   ' possibly it is constant?  Try removing and rerunning.', sep='')
			  stop(err_msg)
			}

			denominators[kkk] <- score_denom
			test_stats[kkk] <- score_num / score_denom
		}
		est_cor <- cor(factor_matrix)

		# Return from here
		return ( list(test_stats=test_stats, cor_mat=est_cor) )
	}

	########################
	# Regular mode
	if (is.null(P_mat)) {
	  W_mat <- diag(W_vec)
	  P_mat <- W_mat - W_mat%*%X_mat %*% solve(t(X_mat)%*%W_mat%*%X_mat) %*% t(X_mat)%*%W_mat
	} else {
	  # If they provided a P_mat, make sure it's the correct dimensions
	  if (nrow(P_mat) != ncol(P_mat) | nrow(P_mat) != nrow(factor_matrix)) {
	    stop('Your P_mat does not have the correct dimensions (n*n).')
	  }
	}

	# Now our score test
	test_stats <- rep(NA, d)
	denominators <- rep(NA, d)
	for(kkk in 1:d)
	{
		# Pick out next SNP, conduct score test (no additional covariates).
		tempF <- factor_matrix[,kkk]
		score_num <- t(tempF) %*% (actual_Y-fitted_Y)
		score_denom <- tryCatch(sqrt(tempF %*% P_mat %*% tempF), warning=function(w) w,
		                        error=function(e) e)

		# We've been getting negative denominators with, for example, very rare SNPs
		if (!(class(score_denom)[1] %in% c("matrix", "array", "numeric"))) {
		    err_msg <- paste('Error in calculating test statistic for factor ', kkk,
		                     ' - possibly it is constant?  Try removing and rerunning.', sep='')
		    stop(err_msg)
		}

		test_stats[kkk] <- score_num / score_denom
		denominators[kkk] <- score_denom
	}

	# Estimate the correlation matrix for the test statistics.
	# Same as the correlation matrix of the SNPs if linear regression and no additional covariates.
	est_cor <- matrix(data=NA, nrow=d, ncol=d)
	for (temp_row in 2:d)
	{
		for (temp_col in 1:(temp_row-1))
		{
			est_cor[temp_row, temp_col] <- t(factor_matrix[,temp_row]) %*% P_mat %*% factor_matrix[,temp_col] / 										(denominators[temp_row] * denominators[temp_col])
			est_cor[temp_col, temp_row] <- est_cor[temp_row, temp_col]
		}
	}

	return ( list(test_stats=test_stats, cor_mat=est_cor) )
}
