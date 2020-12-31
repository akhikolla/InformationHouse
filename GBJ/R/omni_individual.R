#' omni_individual.R
#'
#' Computes the omnibus test statistic combining GBJ, GHC, minP, and SKAT.
#' This version of the function assumes you have the individual factor data (i.e. genotypes)
#' for each subject.  If you only have summary statistics, use omni_ss().
#' You WILL NOT be able to use this function unless you have also loaded
#' the SKAT package (install.packages("SKAT"); library(SKAT)).
#'
#' @param null_model An R regression model fitted using glm().  Do not use lm(), even for linear regression!
#' @param factor_matrix An n*d matrix with each factor (i.e. each SNP) as one column.  There should be no missing data.
#' @param link_function Either "linear" or "logit" or "log".
#' @param num_boots Number of bootstrap repetitions to find correlation matrix of set-based statistics.
#'
#' @return A list with the elements:
#' \item{OMNI}{The observed omnibus test statistic.}
#' \item{OMNI_pvalue}{The p-value of the OMNI test}
#' \item{err_code}{Sometimes if your p-value is very small (< 1*10^(-10)), R may run into numerical
#' issues. This message will alert you if such a situation occurs.}
#'
#' @import mvtnorm
#' @import SKAT
#'
#' @export
#' @examples
#' factor_matrix <- matrix(data=rbinom(n=1000, size=2, prob=0.3), ncol=5)
#' Y <- rnorm(n=200)
#' null_mod <- glm(Y ~ 1)
#' OMNI_individual(null_model=null_mod, factor_matrix=factor_matrix,
#' link_function='linear', num_boots=5)

OMNI_individual <- function(null_model, factor_matrix, link_function, num_boots=100) {

  # Check link function
  if ( !(link_function %in% c('logit', 'linear', 'log')) ) {
    stop('Incorrect link function')
  }

  # Require SKAT
  if (!requireNamespace("SKAT", quietly = TRUE)) {
    stop("SKAT needed for this function to work. Please install it.",
         call. = FALSE)
  }

  # Calculate the SKAT p-value
  skat_Y <- null_model$y
  skat_dmat <- model.matrix(null_model)
  if (link_function == 'logit') {
    skat_obj <- SKAT::SKAT_Null_Model(skat_Y ~ skat_dmat - 1, out_type='D', Adjustment=FALSE)
  } else {
    skat_obj <- SKAT::SKAT_Null_Model(skat_Y ~ skat_dmat - 1, out_type='C', Adjustment=FALSE)
  }
  skat_pvalue <- SKAT::SKAT(factor_matrix, skat_obj)$p.value

  # Get all the marginal test statistics
  score_stats_output <- calc_score_stats(null_model=null_model, factor_matrix=factor_matrix,
                                  link_function=link_function)
  marginal_stats <- score_stats_output$test_stats
  cor_mat <- score_stats_output$cor_mat

  # Do GBJ, GHC, minP
  GHC_output <- GHC(test_stats=marginal_stats, cor_mat=cor_mat)
  GHC_pvalue <- GHC_output$GHC_pvalue
  if (!(is.numeric(GHC_pvalue))) {
    return(list(OMNI=NA, OMNI_pvalue=NA, err_code='Error with GHC'))
  }

  GBJ_output <- GBJ(test_stats=marginal_stats, cor_mat=cor_mat)
  GBJ_pvalue <- GBJ_output$GBJ_pvalue
  if (!(is.numeric(GBJ_pvalue))) {
    return(list(OMNI=NA, OMNI_pvalue=NA, err_code='Error with GBJ'))
  }

  minP_output <- minP(test_stats=marginal_stats, cor_mat=cor_mat)
  minP_pvalue <- minP_output$minP_pvalue
  if (!(is.numeric(minP_pvalue))) {
    return(list(OMNI=NA, OMNI_pvalue=NA, err_code='Error with minP'))
  }

  # Omnibus test statistic
  omni_stat <- min(skat_pvalue, GHC_pvalue, GBJ_pvalue, minP_pvalue)

  ############################################################################
  ############################################################################
  # Simulate correlation matrix under the null
  fitted_values <- null_model$fitted.values
  num_sub <- length(fitted_values)
  boot_pvalues <- matrix(data=NA, nrow=num_boots, ncol=4)
  ell <- rep(NA, num_boots)
  for (i in 1:num_boots) {

    # Simulate new outcome, fit new null models
    if (link_function == 'logit') {
      boot_Y <- rbinom(n=num_sub, size=1, prob=fitted_values)
      boot_skat_obj <- SKAT::SKAT_Null_Model(boot_Y ~ skat_dmat - 1, out_type = 'D', Adjustment=FALSE)
      boot_null_mod <- glm(boot_Y ~ skat_dmat - 1, family=binomial)
    } else if (link_function == 'linear') {
      boot_Y <- fitted_values + rnorm(n=num_sub)
      boot_skat_obj <- SKAT::SKAT_Null_Model(boot_Y ~ skat_dmat - 1, out_type = 'C', Adjustment=FALSE)
      boot_null_mod <- glm(boot_Y ~ skat_dmat - 1)
    } else if (link_function == 'log') {
      boot_Y <- rpois(n=num_sub, lambda=fitted_values)
      boot_skat_obj <- SKAT::SKAT_Null_Model(boot_Y ~ skat_dmat - 1, out_type = 'C', Adjustment=FALSE)
      boot_null_mod <- glm(boot_Y ~ skat_dmat - 1, family=poisson)
    }

    # Bootstrapped SKAT
    boot_skat_p <- SKAT::SKAT(Z=factor_matrix, obj=boot_skat_obj)$p.value

    # Get all the marginal test statistics
    boot_stats <- score_stats_only(null_model=boot_null_mod, factor_matrix=factor_matrix,
                                    link_function=link_function)

    # Do GBJ, GHC, minP
    boot_GHC_p <- GHC(test_stats=boot_stats, cor_mat=cor_mat)$GHC_pvalue
    boot_GBJ_p <- GBJ(test_stats=boot_stats, cor_mat=cor_mat)$GBJ_pvalue
    boot_minP_p <- minP(test_stats=boot_stats, cor_mat=cor_mat)$minP_pvalue

    # Errors?
    boot_vec <- c(boot_minP_p, boot_GHC_p, boot_GBJ_p, boot_skat_p)
    if ( length(which(is.na(boot_vec))) > 0) {
      next
    }

    # Can't transform p-value of 1
    if (length(which(boot_vec == 1))) {
      boot_vec[which(boot_vec == 1)] <- 0.99
    }

    # Transform and record
    boot_pvalues[i, ] <- qnorm(1-boot_vec)
    ell[i] <- min(boot_vec)
  }

  # Remove NA rows, find correlation matrix
  bad_boots <- which(is.na(boot_pvalues[,1]))
  if (length(bad_boots) > 0) {
    boot_pvalues <- boot_pvalues[-bad_boots, ]
  }
  setbased_cor <- cor(boot_pvalues)

  # Get pvalue of omnibus test
  omni_p <- 1 - mvtnorm::pmvnorm(lower=-Inf, upper=rep(qnorm(1-omni_stat), 4), sigma=setbased_cor)[1]

  return ( list(OMNI=omni_stat, OMNI_pvalue=omni_p, err_code=0) )
}
