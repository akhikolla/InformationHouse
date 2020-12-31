#' estimate_ss_cor.R
#'
#' Estimate the correlations between GWAS summary statistics using reference panel eigenvectors
#' and reference panel genotypes.
#'
#'
#' @param ref_pcs An n*m matrix containing PCs calculated from the reference panel. Here n is the
#' number of subjects in the reference panel and m is roughly the number of PCs used in the original
#' analysis which produced the summary statistics.
#' @param ref_genotypes An n*d matrix holding the genotypes from the reference panel, where the d columns
#' correspond to the d SNPs for which we have summary statistics. No missing data allowed.
#' @param link_function Either "linear" or "logit" or "log".
#'
#' @return A list with the elements:
#' \item{cor_mat}{The d*d matrix giving the pairwise correlation of every two test statistics.}
#'
#' @export
#' @examples
#' ref_pcs <- matrix(data=runif(n=1000, min=-0.2, max=0.2), ncol=5)
#' ref_genotypes <- matrix(data=rbinom(n=2000, size=2, prob=0.3), ncol=10)
#' estimate_ss_cor(ref_pcs=ref_pcs, ref_genotypes=ref_genotypes, link_function="linear")

estimate_ss_cor <- function(ref_pcs, ref_genotypes, link_function) {

  # For the summary statistic correlation estimation, we take \mu_0_i to be
  # the same constant for each subject.  This subject gets cancelled out in the
  # estimation, so W is just the identity matrix
  X_mat <- as.matrix(cbind(1, ref_pcs))
  W_mat <- diag(x=1, nrow=nrow(ref_pcs), ncol=nrow(ref_pcs))
  P_mat <- W_mat - X_mat %*% solve(t(X_mat) %*% X_mat) %*% t(X_mat)

  # Estimate the correlation matrix for the test statistics.
  # Same as the correlation matrix of the SNPs if linear regression and no additional covariates.
  est_cor <- matrix(data=NA, nrow=ncol(ref_genotypes), ncol=ncol(ref_genotypes))
  denominators <- rep(NA, ncol(ref_genotypes))

  # First get all the denominator terms
  for (i in 1:ncol(ref_genotypes))
  {
    temp_G <- ref_genotypes[,i]
    denominators[i] <- sqrt(t(temp_G) %*% P_mat %*% temp_G)
  }

  for (temp_row in 2:ncol(ref_genotypes))
  {
    for (temp_col in 1:(temp_row-1))
    {
      est_cor[temp_row, temp_col] <- t(ref_genotypes[,temp_row]) %*% P_mat %*% ref_genotypes[,temp_col] /
          (denominators[temp_row] * denominators[temp_col])
      est_cor[temp_col, temp_row] <- est_cor[temp_row, temp_col]
    }
  }

  return (est_cor)
}
