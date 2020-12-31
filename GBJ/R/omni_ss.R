#' omni_ss.R
#'
#' Computes the omnibus test statistic combining GBJ, GHC, minP, and SKAT.
#' This version of the function assumes you are using GWAS summary statistics.
#' If you individual-level genotype data, use omni_individual().
#'
#' @param test_stats Vector of test statistics for each factor in the set (i.e. marginal
#' test statistic for each SNP in a gene)
#' @param cor_mat d*d matrix of the correlations between all the test statistics in
#' the set, where d is the total number of test statistics in the set.
#' @param num_boots Number of bootstrap repetitions to find correlation matrix of set-based statistics.
#'
#' @return A list with the elements:
#' \item{OMNI}{The observed omnibus test statistic.}
#' \item{OMNI_pvalue}{The p-value of the OMNI test}
#' \item{err_code}{Sometimes if your p-value is very small (< 1*10^(-10)), R may run into numerical
#' issues. This message will alert you if such a situation occurs.}
#'
#' @import mvtnorm
#'
#' @export
#' @examples
#' cor_mat <- matrix(data=0.3, nrow=5, ncol=5)
#' diag(cor_mat) <- 1
#' test_stats <- as.numeric(mvtnorm::rmvnorm(n=1, sigma=cor_mat))
#' OMNI_ss(test_stats=test_stats, cor_mat=cor_mat, num_boots=5)

OMNI_ss <- function(test_stats, cor_mat, num_boots=100) {

  # Make sure our correlation matrix actually has 1 on the diagonal
  # because we use it in rmvnorm here.
  diag(cor_mat) <- 1

  # Do GBJ, GHC, minP
  GHC_output <- GHC(test_stats=test_stats, cor_mat=cor_mat)
  GHC_pvalue <- GHC_output$GHC_pvalue
  if (!(is.numeric(GHC_pvalue))) {
    return(list(OMNI=NA, OMNI_pvalue=NA, err_code='Error with GHC'))
  }

  GBJ_output <- GBJ(test_stats=test_stats, cor_mat=cor_mat)
  GBJ_pvalue <- GBJ_output$GBJ_pvalue
  if (!(is.numeric(GBJ_pvalue))) {
    return(list(OMNI=NA, OMNI_pvalue=NA, err_code='Error with GBJ'))
  }

  minP_output <- minP(test_stats=test_stats, cor_mat=cor_mat)
  minP_pvalue <- minP_output$minP_pvalue
  if (!(is.numeric(minP_pvalue))) {
    return(list(OMNI=NA, OMNI_pvalue=NA, err_code='Error with minP'))
  }

  # Omnibus test statistic
  omni_stat <- min(GHC_pvalue, GBJ_pvalue, minP_pvalue)

  # Simulate correlation matrix under the null
  boot_pvalues <- matrix(data=NA, nrow=num_boots, ncol=3)
  for (i in 1:num_boots) {

    boot_stats <- as.numeric(mvtnorm::rmvnorm(n=1, mean=rep(0, length(test_stats)), sigma=cor_mat))

    # Do GBJ, GHC, minP
    boot_GHC_p <- GHC(test_stats=boot_stats, cor_mat=cor_mat)$GHC_pvalue
    boot_GBJ_p <- GBJ(test_stats=boot_stats, cor_mat=cor_mat)$GBJ_pvalue
    boot_minP_p <- minP(test_stats=boot_stats, cor_mat=cor_mat)$minP_pvalue

    # Errors?
    boot_vec <- c(boot_minP_p, boot_GHC_p, boot_GBJ_p)
    if ( length(which(is.na(boot_vec))) > 0) {
      next
    }

    # Can't transform p-value of 1
    if (length(which(boot_vec == 1))) {
      boot_vec[which(boot_vec == 1)] <- 0.99
    }

    # Transform and record
    boot_pvalues[i, ] <- qnorm(1-boot_vec)
  }

  # Remove NA rows, find correlation matrix
  bad_boots <- which(is.na(boot_pvalues[,1]))
  if (length(bad_boots) > 0) {
    boot_pvalues <- boot_pvalues[-bad_boots, ]
  }
  setbased_cor <- cor(boot_pvalues)

  # Get pvalue of omnibus test
  omni_p <- 1 - mvtnorm::pmvnorm(lower=-Inf, upper=rep(qnorm(1-omni_stat), 3), sigma=setbased_cor)[1]

  return ( list(OMNI=omni_stat, OMNI_pvalue=omni_p, err_code=0) )
}
