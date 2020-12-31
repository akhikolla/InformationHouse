#' score_stats_only.R
#'
#' Starting with individual-level data on p factors, generate score test statistics for each
#' factor for input into GBJ/GHC/HC/BJ/minP.  DOES NOT get the correlations (assumed known).
#'
#' @param null_model An R regression model fitted using glm().  Do not use lm(), even for linear regression!
#' @param factor_matrix An n*d matrix with each factor as one column.  There should be no missing data.
#' @param link_function Either "linear" or "logit" or "log".
#' @param P_mat The projection matrix used in calculation may be passed in to speed up the calculation.
#' See paper for details. Default is null.
#'
#' @return The d score test statistics.
#'
#' @export
#' @examples
#' Y <- rbinom(n=100, size=1, prob=0.5)
#' null_mod <- glm(Y~1, family=binomial(link="logit"))
#' factor_matrix <- matrix(data=rnorm(n=100*5), nrow=100)
#' score_stats_only(null_mod, factor_matrix, "logit")

score_stats_only <- function(null_model, factor_matrix, link_function, P_mat=NULL) {

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

    # Return from here
    return ( test_stats)
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

  return (test_stats)
}
